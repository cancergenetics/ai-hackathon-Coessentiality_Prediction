"""Gated association runner and lockbox reader (DESIGN v3.1 §0 F4/F5, §2.5).

* ``run_association(analysis_id, Y, contexts, ...)`` refuses an analysis whose registry entry
  (docs/registry/<id>.yaml) is missing, invalid, untracked, modified since the last commit, or closed; refuses inputs
  that are not listed with a matching SHA-256 in a committed freeze manifest; refuses any model outside A. It fits HC3
  OLS per (context, unit) with Satterthwaite df, applies the registered FDR family (Benjamini–Bogomolov: families =
  contexts, Simes selection, within-family BH at q R / M), optionally re-tests screening hits by PatientID-blocked
  within-lineage Freedman–Lane permutation, and appends one summary line (no hits) to results/part2/trials.jsonl.
* ``load_A_effects()`` is the A-only effect matrix used for all exploration.
* ``read_lockbox(analysis_id, targets_file)`` is the ONLY code path that returns outcome data (CRISPR / DEMETER2) of
  lockbox models. It requires a committed registry entry and a committed target list naming the analysis, allows one
  read per entry, enforces the read order (all confirmatory directions before benchmark / bridge), and records the
  read in the git-tracked ledger results/part2/lockbox_ledger.jsonl, committed BEFORE any outcome is loaded (if the
  commit fails the ledger is restored from HEAD and nothing is read).
tests/test_part2_lockbox_firewall.py checks by AST/grep that no other Part 2 code reaches lockbox outcomes.
"""
from __future__ import annotations

import hashlib
import json
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd

from part2 import paths, registry, stats
from part2.registry import RegistryError, git


class FreezeError(Exception):
    pass


class FirewallError(Exception):
    pass


def _sha_file(p) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


def _sha_text(s: str) -> str:
    return hashlib.sha256(s.encode()).hexdigest()


def _default_outcome_loader(model_ids, genes, source: str) -> pd.DataFrame:
    """Outcome rows for the given models (CRISPR gene effect or DEMETER2), restricted to `genes` (DepMap column keys
    'SYMBOL (ENTREZ)' for CRISPR; 'SYMBOL (ENTREZ)' row keys for DEMETER2) when given."""
    from part2 import data
    if source == "crispr":
        ge = data.load_gene_effect()
        ids = [m for m in model_ids if m in set(ge.models)]
        X = ge.subset(ids)
        df = pd.DataFrame(X, index=ge.models[ge.rows(ids)], columns=ge.columns)
    elif source == "demeter2":
        d2 = pd.read_csv(paths.DEMETER2_CSV, index_col=0)
        m = data.load_models()
        name2id = pd.Series(m.index, index=m.CCLEName).dropna()
        name2id = name2id[~name2id.index.duplicated()]
        d2 = d2.T
        d2.index = [name2id.get(c) for c in d2.index]
        d2 = d2[pd.Index(d2.index).notna()]
        d2 = d2[~d2.index.duplicated()]
        df = d2.reindex([i for i in model_ids if i in d2.index])
    else:
        raise ValueError(source)
    if genes is not None:
        df = df.reindex(columns=list(genes))
    return df


@dataclass
class Project:
    """Paths + outcome loader; defaults are the real project. Tests use a temporary git repository."""
    root: Path = paths.PROJECT
    outcome_loader: Callable = field(default=_default_outcome_loader)

    @property
    def registry_dir(self):
        return self.root / "docs" / "registry"

    @property
    def lockbox_dir(self):
        return self.root / "results" / "part2" / "lockbox"

    @property
    def trials(self):
        return self.root / "results" / "part2" / "trials.jsonl"

    @property
    def ledger(self):
        return self.root / "results" / "part2" / "lockbox_ledger.jsonl"

    def rel(self, p) -> str:
        import os
        p = Path(p)
        p = p if p.is_absolute() else self.root / p
        return os.path.relpath(os.path.abspath(p), os.path.abspath(self.root))   # no symlink resolution


# ------------------------------------------------------------------------------------------------ freeze checks
def _manifest_files(proj: Project, manifest_rel: str) -> dict:
    """Committed, clean manifest -> {relpath: sha256}. Lockbox manifests list files relative to their folder."""
    mp = proj.root / manifest_rel
    if not mp.exists():
        raise FreezeError(f"freeze manifest {manifest_rel} does not exist")
    try:
        registry.committed_blob(proj.root, mp)
    except RegistryError as e:
        raise FreezeError(f"freeze manifest not committed: {e}") from e
    m = json.loads(mp.read_text())
    files = {}
    for k, v in m.get("files", {}).items():
        kp = Path(k)
        if len(kp.parts) == 1:                                   # relative to the manifest's folder
            kp = Path(manifest_rel).parent / kp
        files[str(kp)] = v
    return files


def check_frozen(proj: Project, inputs, manifests) -> dict:
    """Every input path must be listed with its current SHA-256 in one of the committed manifests. Returns
    {relpath: sha256}."""
    listed = {}
    for m in manifests:
        listed.update(_manifest_files(proj, proj.rel(proj.root / m) if not Path(m).is_absolute() else proj.rel(m)))
    out = {}
    for p in inputs:
        rel = proj.rel(p)
        if rel not in listed:
            raise FreezeError(f"input {rel} is not listed in any committed freeze manifest {list(manifests)}")
        h = _sha_file(proj.root / rel)
        if h != listed[rel]:
            raise FreezeError(f"input {rel} changed since freeze (sha256 {h[:12]} != {listed[rel][:12]})")
        out[rel] = h
    return out


def _model_list(proj: Project, name: str) -> list[str]:
    """A model list file (A or L) verified against the committed lockbox manifest."""
    man_rel = proj.rel(proj.lockbox_dir / "manifest.json")
    files = _manifest_files(proj, man_rel)
    rel = proj.rel(proj.lockbox_dir / name)
    txt = (proj.root / rel).read_text()
    if files.get(rel) != _sha_text(txt):
        raise FreezeError(f"{rel} does not match the committed lockbox manifest")
    return [l for l in txt.splitlines() if l]


def a_models(proj: Project | None = None) -> list[str]:
    return _model_list(proj or Project(), "A_models.txt")


def load_A_effects(genes=None, source: str = "crispr", proj: Project | None = None) -> pd.DataFrame:
    """A-only outcome matrix (exploration / discovery). Lockbox rows are never loaded here."""
    proj = proj or Project()
    A = a_models(proj)
    df = proj.outcome_loader(A, genes, source)
    if not set(df.index) <= set(A):
        raise FirewallError("outcome loader returned non-A models")
    return df


# ------------------------------------------------------------------------------------------------ git helpers
def _head(proj: Project) -> str | None:
    r = git(proj.root, "rev-parse", "HEAD", check=False)
    return r.stdout.strip() if r.returncode == 0 else None


def _append_jsonl(path: Path, rec: dict):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "a") as fh:
        fh.write(json.dumps(rec, sort_keys=True, default=str) + "\n")


# ------------------------------------------------------------------------------------------------ association
def _design(ctx: pd.Series, cov: pd.DataFrame | None) -> tuple[np.ndarray, list[str], np.ndarray]:
    """[1, context, covariates (constant columns dropped)] on complete cases. Returns X, names, row mask."""
    parts = [pd.Series(1.0, index=ctx.index, name="intercept"), ctx.astype(float).rename("context")]
    if cov is not None and cov.shape[1]:
        parts.append(cov.astype(float))
    D = pd.concat(parts, axis=1)
    ok = D.notna().all(axis=1).values
    D = D[ok]
    keep = [c for c in D.columns if c in ("intercept", "context") or D[c].nunique() > 1]
    return D[keep].values, keep, ok


def run_association(analysis_id: str, Y: pd.DataFrame, contexts: pd.DataFrame, covariates: pd.DataFrame | None = None,
                    *, inputs=(), proj: Project | None = None, blocks: pd.Series | None = None,
                    strata: pd.Series | None = None, n_perm: int = 0, perm_screen_q: float = 0.1,
                    perm_seed: int = 4000, min_context_n: int = 1) -> pd.DataFrame:
    """Y: A models x tested units (complex scores); contexts: A models x registered contexts (0/1, NaN = excluded);
    covariates: A models x numeric covariates (already expanded, e.g. lineage dummies). Returns the per-test table
    to the caller; the trials ledger receives counts only."""
    proj = proj or Project()
    t0 = time.time()
    entry, blob = registry.load_entry(proj.registry_dir, analysis_id, root=proj.root)
    if entry["status"] == "closed":
        raise RegistryError(f"{analysis_id} is closed")
    if entry["fdr_family"]["method"] not in ("benjamini_bogomolov", "bh", "none"):   # e.g. stagewise_bb (H2)
        raise RegistryError(f"fdr method {entry['fdr_family']['method']!r} is registered but not implemented here")
    # frozen inputs: the lockbox manifest is always required (A must be frozen)
    manifests = list(dict.fromkeys(["results/part2/lockbox/manifest.json"] + list(entry["frozen_inputs"])))
    in_hashes = check_frozen(proj, inputs, manifests)
    A = set(a_models(proj))
    for name, df in (("Y", Y), ("contexts", contexts)) + ((("covariates", covariates),) if covariates is not None
                                                          else ()):
        extra = set(df.index) - A
        if extra:
            raise FirewallError(f"{name} contains {len(extra)} models outside A (e.g. {sorted(extra)[:3]})")
    reg_ctx = registry.contexts_of(entry)
    unknown = set(contexts.columns) - set(reg_ctx)
    if unknown:
        raise RegistryError(f"contexts not registered in {analysis_id}: {sorted(unknown)[:5]}")
    if registry.set_hash(Y.columns) != entry["complex_set_hash"]:
        raise RegistryError("tested units do not match the registered complex_set_hash")
    idx = Y.index
    C = contexts.reindex(idx)
    Z = covariates.reindex(idx) if covariates is not None else None
    rows = []
    for c in C.columns:
        sign = reg_ctx[c]["sign"]
        alt = {"+": "less", "-": "greater", "two-sided": "two-sided"}[sign]   # '+' = more dependent = lower effect
        X, names, ok = _design(C[c], Z)
        if np.nansum(X[:, 1] == 1) < min_context_n:
            continue
        Yc = Y.values[ok]
        nanpat = np.isnan(Yc)
        groups = pd.Series(range(Yc.shape[1])).groupby([nanpat[:, k].tobytes() for k in range(Yc.shape[1])])
        for _, cols in groups:
            cols = cols.values
            rmask = ~nanpat[:, cols[0]]
            d = stats.HC3Design(X[rmask], 1)
            b, se, t = d.fit(Yc[rmask][:, cols])
            p = stats.pvalue_t(t, d.df, alt)
            for k, j in enumerate(cols):
                rows.append(dict(context=c, unit=Y.columns[j], n=int(rmask.sum()),
                                 n_context=int(X[rmask, 1].sum()), beta=b[k], se=se[k], t=t[k], df=d.df, p=p[k],
                                 alternative=alt))
    res = pd.DataFrame(rows)
    ff = entry["fdr_family"]
    q = float(ff["q"])
    if len(res):
        if ff["method"] == "benjamini_bogomolov":
            fam = res["context"] if ff["families"] == "contexts" else res["unit"].astype(str)
            bb = stats.benjamini_bogomolov(res.p.values, fam.values, q)
            res["q_adj"], res["reject"] = bb["adj"], bb["reject"]
            n_sel = bb["R"]
        elif ff["method"] == "bh":
            rj, adj = stats.bh(res.p.values, q)
            res["q_adj"], res["reject"], n_sel = adj, rj, None
        else:
            res["q_adj"], res["reject"], n_sel = np.nan, False, None
        res["p_perm"] = np.nan
        if n_perm > 0:
            hits = res.index[res.q_adj.fillna(1) < perm_screen_q]
            rng = np.random.default_rng(perm_seed)
            for i in hits:
                r = res.loc[i]
                X, names, ok = _design(C[r.context], Z)
                y = Y[r.unit].values[ok]
                m = ~np.isnan(y)
                ids = idx[ok][m]
                bl = None
                if blocks is not None:
                    bl = blocks.reindex(ids)
                    bl = bl.where(bl.notna(), pd.Series("__" + ids, index=ids)).values
                perms = stats.blocked_permutations(
                    n_perm, rng, blocks=bl,
                    strata=None if strata is None else strata.reindex(ids).fillna("NA").values, n=len(ids))
                res.loc[i, "p_perm"] = stats.freedman_lane(y[m], X[m], 1, perms, r.alternative)["p_perm"]
            if len(hits) and ff["method"] == "benjamini_bogomolov":
                pp = res.p.where(res.p_perm.isna(), res.p_perm).values      # permutation re-entry, BB re-run
                bb = stats.benjamini_bogomolov(pp, fam.values, q)
                res["q_adj_perm"], res["reject_perm"] = bb["adj"], bb["reject"]
    else:
        n_sel = 0
    rec = dict(time=time.strftime("%Y-%m-%dT%H:%M:%S%z"), analysis_id=analysis_id, registry_blob=blob,
               head=_head(proj), direction=entry["direction"], inputs=in_hashes,
               models_hash=registry.set_hash(idx), n_models=int(len(idx)), n_contexts=int(C.shape[1]),
               n_units=int(Y.shape[1]), n_tests=int(len(res)),
               n_p_lt_005=int((res.p < 0.05).sum()) if len(res) else 0,
               n_families_selected=n_sel, n_reject=int(res.reject.sum()) if len(res) else 0,
               n_perm=n_perm, n_perm_retested=int(res.p_perm.notna().sum()) if len(res) else 0,
               runtime_s=round(time.time() - t0, 2))
    _append_jsonl(proj.trials, rec)
    return res


# ------------------------------------------------------------------------------------------------ lockbox
def _ledger(proj: Project) -> list[dict]:
    if not proj.ledger.exists():
        return []
    return [json.loads(l) for l in proj.ledger.read_text().splitlines() if l.strip()]


def read_lockbox(analysis_id: str, targets_file, proj: Project | None = None, source: str = "crispr") -> pd.DataFrame:
    """Return outcome data of lockbox models for the genes of a committed target list (JSON:
    {"analysis_id": ..., "genes": [...], ...}). One read per registry entry; ledger committed before the read."""
    proj = proj or Project()
    entry, blob = registry.load_entry(proj.registry_dir, analysis_id, root=proj.root)
    tp = Path(targets_file)
    tp = tp if tp.is_absolute() else proj.root / tp
    try:
        t_blob = registry.committed_blob(proj.root, tp)
    except RegistryError as e:
        raise FreezeError(f"target list must be committed before a lockbox read: {e}") from e
    targets = json.loads(tp.read_text())
    if targets.get("analysis_id") != analysis_id:
        raise FreezeError("target list names a different analysis_id")
    genes = targets.get("genes")
    if not genes:
        raise FreezeError("target list has no genes")
    if proj.ledger.exists():
        try:
            registry.committed_blob(proj.root, proj.ledger)
        except RegistryError as e:
            raise FreezeError(f"lockbox ledger has uncommitted changes: {e}") from e
    led = _ledger(proj)
    if any(r["analysis_id"] == analysis_id and r.get("source", "crispr") == source for r in led):
        raise FirewallError(f"{analysis_id} already read the lockbox ({source}); one read per registry entry")
    if entry["direction"] in ("benchmark", "bridge"):
        confirm = [k for k, e in registry.load_all(proj.registry_dir).items() if e["direction"] in registry.CONFIRMATORY]
        done = {r["analysis_id"] for r in led}
        missing = [k for k in confirm if k not in done]
        if missing or not confirm:
            raise FirewallError(f"{entry['direction']} lockbox read before all confirmatory reads: missing {missing}")
    rec = dict(time=time.strftime("%Y-%m-%dT%H:%M:%S%z"), analysis_id=analysis_id, registry_blob=blob,
               targets_file=proj.rel(tp), targets_blob=t_blob, head=_head(proj), source=source,
               n_genes=len(genes), direction=entry["direction"])
    old = proj.ledger.read_text() if proj.ledger.exists() else None
    _append_jsonl(proj.ledger, rec)
    rel = proj.rel(proj.ledger)
    r1 = git(proj.root, "add", "--", rel, check=False)
    r2 = git(proj.root, "commit", "-m", f"lockbox read: {analysis_id} ({source})", "--", rel, check=False) \
        if r1.returncode == 0 else r1
    if r2.returncode != 0:
        if old is None:
            git(proj.root, "rm", "--cached", "-q", "--", rel, check=False)
            proj.ledger.unlink()
        else:
            proj.ledger.write_text(old)
            git(proj.root, "reset", "-q", "--", rel, check=False)
        raise FreezeError(f"could not commit the ledger entry; nothing read: {r2.stderr.strip()}")
    Lm = _model_list(proj, "lockbox_models.txt")
    df = proj.outcome_loader(Lm, genes, source)
    if not set(df.index) <= set(Lm):
        raise FirewallError("outcome loader returned non-lockbox models")
    return df
