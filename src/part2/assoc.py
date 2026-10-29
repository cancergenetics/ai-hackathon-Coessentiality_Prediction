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
* ``resume_lockbox_read(analysis_id)`` is the only other path: the crash-recovery copy that ``read_lockbox(...,
  persist=True)`` stores under data_cache/part2/lockbox_reads/ (gitignored; SHA-256 committed at once), for an entry
  that already read and is still open; it never calls the loader. Retention rule (PI, R1): the copies stay only in the
  local gitignored data_cache until the manuscript's final reproduction check; they are never copied, committed or read
  by any other code, and deleting them needs a PI note in the research tree.
tests/test_part2_lockbox_firewall.py checks by AST/grep that no other Part 2 code reaches lockbox outcomes, and that only
assoc.py, confirm.py and the lockbox driver name the persisted reads.
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
from part2.cscore import is_lineage_context
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

    @property
    def cache_dir(self):
        return self.root / "data_cache" / "part2"

    @property
    def results_dir(self):
        return self.root / "results" / "part2"

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


FE_PREFIXES = ("lineage_FE:", "ModelType:", "GrowthPattern:", "library:")


def _design(ctx: pd.Series, cov: pd.DataFrame | None, return_dropped: bool = False):
    """[1, context, covariates] on complete cases. Constant covariates are dropped; one-hot / multi-hot indicator
    columns (FE_PREFIXES) with < 2 models on either side are merged into the reference level; covariates that are
    collinear with the columns before them are dropped (greedy, in order; e.g. ModelType identical to
    ScreenType3DO). Returns X, names, row mask (, dropped names)."""
    parts = [pd.Series(1.0, index=ctx.index, name="intercept"), ctx.astype(float).rename("context")]
    if cov is not None and cov.shape[1]:
        parts.append(cov.astype(float))
    D = pd.concat(parts, axis=1)
    ok = D.notna().all(axis=1).values
    D = D[ok]
    n = len(D)
    keep, dropped = ["intercept", "context"], []
    for c in D.columns[2:]:
        v = D[c]
        if v.nunique() <= 1:
            dropped.append(c)
            continue
        if str(c).startswith(FE_PREFIXES) and (min(v.sum(), n - v.sum()) < 2):
            dropped.append(c)
            continue
        keep.append(c)
    X = D[keep].values
    if X.shape[1] > 2 and n > 0:
        sel = [0, 1]
        for k in range(2, X.shape[1]):
            r_new = np.linalg.matrix_rank(X[:, sel + [k]])
            if r_new == len(sel) + 1:
                sel.append(k)
            else:
                dropped.append(keep[k])
        X = X[:, sel]
        keep = [keep[k] for k in sel]
    return (X, keep, ok, dropped) if return_dropped else (X, keep, ok)


def _fdr(res: pd.DataFrame, pcol: str, ff: dict, parent_selection: dict | None):
    """Adjusted p (NaN where not tested / not eligible), rejections, number of selected families."""
    q = float(ff["q"])
    adj = np.full(len(res), np.nan)
    rej = np.zeros(len(res), bool)
    ev = res.evaluable.values.astype(bool) & np.isfinite(res[pcol].values)
    method = ff["method"]
    n_sel = None
    if method == "benjamini_bogomolov":
        if ev.any():
            fam = (res.context if ff.get("families", "contexts") == "contexts" else res.unit).astype(str).values
            bb = stats.benjamini_bogomolov(res[pcol].values[ev], fam[ev], q)
            adj[ev], rej[ev], n_sel = bb["adj"], bb["reject"], bb["R"]
        else:
            n_sel = 0
    elif method == "bh":
        if ev.any():
            r_, a_ = stats.bh(res[pcol].values[ev], q)
            adj[ev], rej[ev] = a_, r_
    elif method == "stagewise_bb":
        el = ev & res.eligible.values.astype(bool)
        fam = (res.context.astype(str) + "|" + res.parent.astype(str)).values
        lev = res.context.map(parent_selection["level"]).astype(float).values
        n_sel = 0
        for f in sorted(set(fam[el])):
            i = np.flatnonzero(el & (fam == f))
            level = float(lev[i[0]])
            if not level > 0:
                continue
            r_, a_ = stats.bh(res[pcol].values[i], level)
            adj[i] = np.minimum(a_ * q / level, 1.0)        # on the q scale: reject iff adj <= q
            rej[i] = r_
            n_sel += 1
    return adj, rej, n_sel


def _committed_text(proj: "Project", path: Path) -> str:
    r = git(proj.root, "show", f"HEAD:{proj.rel(path)}", check=False)
    return r.stdout if r.returncode == 0 else ""


def check_run_after(proj: "Project", parent_directions) -> None:
    """stagewise_bb: every parent direction must have a discovery trial line in the COMMITTED trials ledger."""
    txt = _committed_text(proj, proj.trials)
    dirs = {json.loads(l).get("direction") for l in txt.splitlines() if l.strip()}
    missing = [d for d in parent_directions if d not in dirs]
    if missing:
        raise RegistryError(f"parent directions {missing} have no committed discovery run in {proj.rel(proj.trials)}")


def parent_selection_from(results: dict, parent_of: dict, q: float, min_children: int = 2) -> dict:
    """Stage-wise selection for H2 from the parent directions' discovery tables.
    results: {direction: run_association table}; parent_of: {sub-module: parent}. A (context c, parent P) is selected
    when P was discovered for c (column `discovery`) and P has >= min_children sub-modules. The per-family level is
    q R_x / m_x of the parent direction x that tested c (R_x discoveries, m_x evaluable tests)."""
    nchild = pd.Series(parent_of).value_counts()
    selected, level, direction_of = set(), {}, {}
    for x, r in results.items():
        ev = r[r.evaluable]
        R, m = int(ev.discovery.sum()), int(len(ev))
        for c in sorted(set(ev.context)):
            level[c] = q * R / m if m else 0.0
            direction_of[c] = x
        for c, u in ev.loc[ev.discovery, ["context", "unit"]].itertuples(index=False):
            if nchild.get(u, 0) >= min_children:
                selected.add((c, u))
    return dict(selected=selected, level=level, parent_of=dict(parent_of), direction_of=direction_of)


def run_association(analysis_id: str, Y, contexts: pd.DataFrame, covariates=None, *, inputs=(),
                    proj: Project | None = None, blocks: pd.Series | None = None, strata: pd.Series | None = None,
                    n_perm: int = 0, perm_screen_q: float = 0.1, perm_screen_p: float | None = None,
                    perm_seed: int = 4000, min_context_n: int = 1,
                    min_ref_n: int = 1, pairs=None, lineage_contexts=None, parent_selection: dict | None = None,
                    tag: str = "discovery", save: bool = True, release_hits: bool = False) -> pd.DataFrame:
    """Gated association run (DESIGN §2.5, F5).

    Y: A models x tested units (complex scores) -- or {context: DataFrame} when scores depend on the context
       (member removals, H1 alignment); every frame has the same unit columns, whose set hash must equal the entry's
       complex_set_hash. A unit column that is all-NaN for a context is 'not evaluable' there.
    contexts: A models x registered contexts (1 / 0 / NaN = excluded).
    covariates: A models x numeric covariates (expanded, see covariates.assemble) -- or {context: DataFrame}.
    pairs: optional iterable of (context, unit) to test (H1: each loss column only with its own event components);
       default = every context x every unit.
    Evaluable test: >= min_context_n context and >= min_ref_n reference models after complete-case (registry: 12/12)
       and a non-empty score; non-evaluable tests are listed (evaluable=False) and excluded from m.
    FDR (entry fdr_family): benjamini_bogomolov, bh, stagewise_bb (needs parent_selection from
       parent_selection_from() and committed parent-direction trial lines), none.
    Permutation (n_perm > 0): tests with adjusted p < perm_screen_q are re-tested by exact Freedman-Lane permutation
       of the HC3 t, whole PatientID blocks exchanged within lineage strata (`strata`), unrestricted (strata ignored)
       for lineage contexts; the permutation p replaces the analytic p and the FDR procedure is re-run
       (q_adj_perm / reject_perm). `discovery` = reject_perm if permutations ran, else reject.
    Outputs: the full table is returned and written to data_cache/part2/assoc/<id>/<tag>.parquet (gitignored);
       results/part2/assoc/<id>/<tag>_summary.json holds counts (the discovery list only if release_hits); the
       trials ledger gets counts only (no units, no contexts)."""
    proj = proj or Project()
    t0 = time.time()
    entry, blob = registry.load_entry(proj.registry_dir, analysis_id, root=proj.root)
    if entry["status"] == "closed":
        raise RegistryError(f"{analysis_id} is closed")
    ff = entry["fdr_family"]
    if ff["method"] not in ("benjamini_bogomolov", "bh", "none", "stagewise_bb"):
        raise RegistryError(f"fdr method {ff['method']!r} is not implemented")
    if ff["method"] == "stagewise_bb":
        if parent_selection is None:
            raise RegistryError("stagewise_bb needs parent_selection (parent_selection_from the parent directions)")
        check_run_after(proj, ff["parent_directions"])
    manifests = list(dict.fromkeys(["results/part2/lockbox/manifest.json"] + list(entry["frozen_inputs"])))
    in_hashes = check_frozen(proj, inputs, manifests)
    A = set(a_models(proj))
    Yd = Y if isinstance(Y, dict) else None
    Zd = covariates if isinstance(covariates, dict) else None
    frames = [("contexts", contexts)]
    frames += [(f"Y[{k}]", v) for k, v in Yd.items()] if Yd is not None else [("Y", Y)]
    if covariates is not None:
        frames += [(f"covariates[{k}]", v) for k, v in Zd.items()] if Zd is not None else [("covariates", covariates)]
    for name, df in frames:
        extra = set(df.index) - A
        if extra:
            raise FirewallError(f"{name} contains {len(extra)} models outside A (e.g. {sorted(extra)[:3]})")
    reg_items = {c["id"]: c for c in entry["context_family"]["contexts"]}
    reg_ctx = registry.contexts_of(entry)
    unknown = set(contexts.columns) - set(reg_ctx)
    if unknown:
        raise RegistryError(f"contexts not registered in {analysis_id}: {sorted(unknown)[:5]}")
    if Yd is not None:
        missing = set(contexts.columns) - set(Yd)
        if missing:
            raise ValueError(f"no score frame for contexts {sorted(missing)[:5]}")
        colsets = {tuple(sorted(map(str, v.columns))) for v in Yd.values()}
        if len(colsets) != 1:
            raise ValueError("per-context score frames must share one unit column set")
        unit_cols = list(next(iter(Yd.values())).columns)
    else:
        unit_cols = list(Y.columns)
    if registry.set_hash(unit_cols) != entry["complex_set_hash"]:
        raise RegistryError("tested units do not match the registered complex_set_hash")
    allowed = None
    if pairs is not None:
        allowed = {}
        for c, u in pairs:
            if u not in unit_cols:
                raise ValueError(f"pair unit {u!r} is not a score column")
            allowed.setdefault(c, []).append(u)
    if lineage_contexts is None:
        lineage_contexts = {c for c in contexts.columns if is_lineage_context(c, reg_items.get(c))}
    parent_of = (parent_selection or {}).get("parent_of", {})
    selected = (parent_selection or {}).get("selected", set())
    if ff["method"] == "stagewise_bb":            # only the stage-wise eligible subset is tested (no open screen)
        allowed = {c: [u for u in (allowed.get(c, []) if allowed is not None else unit_cols)
                       if (c, parent_of.get(u)) in selected] for c in contexts.columns}
    rows = []
    designs = {}
    for c in contexts.columns:
        Yc = (Yd[c] if Yd is not None else Y)
        idx = Yc.index
        Zc = (Zd.get(c) if Zd is not None else covariates)
        Zc = Zc.reindex(idx) if Zc is not None else None
        sign = reg_ctx[c]["sign"]
        alt = {"+": "less", "-": "greater", "two-sided": "two-sided"}[sign]   # '+' = more dependent = lower effect
        X, names, ok, dropped = _design(contexts[c].reindex(idx), Zc, return_dropped=True)
        designs[c] = (idx, Zc)
        units = allowed.get(c, []) if allowed is not None else unit_cols
        if not units:
            continue
        Yv = Yc[units].values[ok]
        nanpat = np.isnan(Yv)
        base = dict(context=c, alternative=alt, covariates_dropped=";".join(map(str, dropped)))
        groups = pd.Series(range(Yv.shape[1])).groupby([nanpat[:, k].tobytes() for k in range(Yv.shape[1])])
        for _, cols in groups:
            cols = cols.values
            rmask = ~nanpat[:, cols[0]]
            n_ctx = int(X[rmask, 1].sum()) if rmask.any() else 0
            n_ref = int(rmask.sum()) - n_ctx
            reason = ""
            if not rmask.any():
                reason = "no score"
            elif n_ctx < min_context_n or n_ref < min_ref_n:
                reason = f"n_context {n_ctx} / n_ref {n_ref} below {min_context_n}/{min_ref_n}"
            d = None
            if not reason:
                try:
                    d = stats.HC3Design(X[rmask], 1)
                except (np.linalg.LinAlgError, ValueError) as e:
                    reason = f"design: {e}"
            if d is None:
                for j in cols:
                    rows.append(dict(base, unit=units[j], n=int(rmask.sum()), n_context=n_ctx, beta=np.nan,
                                     se=np.nan, t=np.nan, df=np.nan, p=np.nan, evaluable=False, reason=reason))
                continue
            b, se, t = d.fit(Yv[rmask][:, cols])
            p = stats.pvalue_t(t, d.df, alt)
            for k, j in enumerate(cols):
                rows.append(dict(base, unit=units[j], n=int(rmask.sum()), n_context=n_ctx, beta=b[k], se=se[k],
                                 t=t[k], df=d.df, p=p[k], evaluable=True, reason=""))
    cols_out = ["context", "unit", "n", "n_context", "beta", "se", "t", "df", "p", "alternative", "evaluable",
                "reason", "covariates_dropped"]
    res = pd.DataFrame(rows, columns=cols_out) if rows else pd.DataFrame(columns=cols_out)
    res["evaluable"] = res.evaluable.astype(bool)
    res["parent"] = [parent_of.get(u, "") for u in res.unit]
    if ff["method"] == "stagewise_bb":
        res["eligible"] = [(c, p) in selected for c, p in zip(res.context, res.parent)]
    res["q_adj"], res["reject"], n_sel = _fdr(res, "p", ff, parent_selection)
    res["p_perm"] = np.nan
    n_sel_perm = None
    if n_perm > 0 and len(res) and ff["method"] != "none":
        scr = res.q_adj.fillna(1) < perm_screen_q
        if perm_screen_p is not None:                 # pre-declared calibration response: also raw p < 0.01
            scr |= res.evaluable & (res.p.fillna(1) < perm_screen_p)
        hits = res.index[scr]
        for i in hits:
            r = res.loc[i]
            idx, Zc = designs[r.context]
            Yc = (Yd[r.context] if Yd is not None else Y)
            X, _, ok = _design(contexts[r.context].reindex(idx), Zc)
            y = Yc[r.unit].values[ok]
            m = ~np.isnan(y)
            ids = idx[ok][m]
            st = None if (strata is None or r.context in lineage_contexts) else \
                strata.reindex(ids).fillna("NA").astype(str).values
            res.loc[i, "p_perm"] = permutation_p(y[m], X[m], ids, r.alternative, n_perm,
                                                 np.random.default_rng(perm_seed + int(i)), blocks, st)
        pp = res.p.where(res.p_perm.isna(), res.p_perm)
        res["p_final"] = pp
        res["q_adj_perm"], res["reject_perm"], n_sel_perm = _fdr(res, "p_final", ff, parent_selection)
        res["discovery"] = res.reject_perm
    else:
        res["p_final"] = res.p
        res["discovery"] = res.reject
    res["discovery"] = res.discovery.astype(bool)
    rec = dict(time=time.strftime("%Y-%m-%dT%H:%M:%S%z"), analysis_id=analysis_id, registry_blob=blob,
               head=_head(proj), direction=entry["direction"], tag=tag, inputs=in_hashes,
               models_hash=registry.set_hash(set().union(*[set(df.index) for _, df in frames[1:2]])),
               n_models=int(len(contexts)), n_contexts=int(contexts.shape[1]), n_units=int(len(unit_cols)),
               n_tests=int(len(res)), n_evaluable=int(res.evaluable.sum()),
               n_p_lt_005=int((res.p < 0.05).sum()) if len(res) else 0,
               n_families_selected=n_sel, n_reject=int(res.reject.sum()) if len(res) else 0,
               n_perm=n_perm, n_perm_retested=int(res.p_perm.notna().sum()) if len(res) else 0,
               n_families_selected_perm=n_sel_perm, n_discovery=int(res.discovery.sum()) if len(res) else 0,
               fdr_method=ff["method"], runtime_s=round(time.time() - t0, 2))
    _append_jsonl(proj.trials, rec)
    if save:
        save_outputs(proj, analysis_id, res, tag, rec, release_hits=release_hits)
    return res


def permutation_p(y, X, ids, alternative: str, n_perm: int, rng, blocks: pd.Series | None = None, strata=None) -> float:
    """Exact Freedman-Lane permutation p of the HC3 t of column 1 (PatientID blocks within strata; GPD tail)."""
    bl = None
    if blocks is not None:
        b = blocks.reindex(ids)
        bl = b.where(b.notna(), pd.Series(["__" + str(i) for i in ids], index=ids)).astype(str).values
    perms = stats.blocked_permutations(n_perm, rng, blocks=bl, strata=strata, n=len(ids))
    return stats.freedman_lane(y, X, 1, perms, alternative)["p_perm"]


def save_outputs(proj: Project, analysis_id: str, res: pd.DataFrame, tag: str, rec: dict,
                 release_hits: bool = False) -> dict:
    """Full table -> data_cache (gitignored); counts (+ discovery list only if release_hits) -> results/ summary."""
    cdir = proj.cache_dir / "assoc" / analysis_id
    cdir.mkdir(parents=True, exist_ok=True)
    out = res.copy()
    for c in out.columns:
        if out[c].dtype == object:
            out[c] = out[c].astype(str)
    out.to_parquet(cdir / f"{tag}.parquet", index=False)
    ev = res[res.evaluable] if len(res) else res
    per_ctx = ev.groupby("context").agg(n_tests=("unit", "size"), n_p_lt_005=("p", lambda v: int((v < 0.05).sum())),
                                        n_discovery=("discovery", "sum")) if len(ev) else pd.DataFrame()
    summ = dict(analysis_id=analysis_id, tag=tag, registry_blob=rec["registry_blob"], head=rec["head"],
                fdr_method=rec["fdr_method"], n_tests=rec["n_tests"], n_evaluable=rec["n_evaluable"],
                n_not_evaluable=int(rec["n_tests"] - rec["n_evaluable"]), n_discovery=rec["n_discovery"],
                n_perm_retested=rec["n_perm_retested"], full_table=proj.rel(cdir / f"{tag}.parquet"),
                full_table_sha256=_sha_file(cdir / f"{tag}.parquet"),
                hits_released=bool(release_hits))
    if release_hits:
        summ["per_context"] = {k: {kk: int(vv) for kk, vv in v.items()} for k, v in per_ctx.to_dict("index").items()}
        hits = res[res.discovery].sort_values("p_final")
        summ["discoveries"] = hits[["context", "unit", "beta", "t", "p", "p_final"]].to_dict("records")
    sdir = proj.results_dir / "assoc" / analysis_id
    sdir.mkdir(parents=True, exist_ok=True)
    (sdir / f"{tag}_summary.json").write_text(json.dumps(summ, indent=1, sort_keys=True, default=float) + "\n")
    return summ


# ------------------------------------------------------------------------------------------------ lockbox
def _ledger(proj: Project) -> list[dict]:
    if not proj.ledger.exists():
        return []
    return [json.loads(l) for l in proj.ledger.read_text().splitlines() if l.strip()]


def read_lockbox(analysis_id: str, targets_file, proj: Project | None = None, source: str = "crispr",
                 persist: bool = False) -> pd.DataFrame:
    """Return outcome data of lockbox models for the genes of a committed target list (JSON:
    {"analysis_id": ..., "genes": [...], ...}). One read per registry entry (any source; review S1); ledger committed
    before the read. persist=True: the matrix is stored at once under data_cache/part2/lockbox_reads/ (gitignored) and
    its SHA-256 is committed to results/part2/lockbox_reads.jsonl (crash recovery: resume_lockbox_read)."""
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
    if genes is None or not isinstance(genes, list):
        raise FreezeError("target list has no 'genes' list")
    if not genes and not targets.get("zero_targets"):
        raise FreezeError("target list has no genes (a zero-target close needs \"zero_targets\": true)")
    if proj.ledger.exists():
        try:
            registry.committed_blob(proj.root, proj.ledger)
        except RegistryError as e:
            raise FreezeError(f"lockbox ledger has uncommitted changes: {e}") from e
    led = _ledger(proj)
    if any(r["analysis_id"] == analysis_id for r in led):
        raise FirewallError(f"{analysis_id} already read the lockbox; one read per registry entry (any source)")
    if entry["direction"] in ("benchmark", "bridge"):
        allent = registry.load_all(proj.registry_dir)
        confirm = [k for k, e in allent.items() if e["direction"] in registry.CONFIRMATORY]
        done = {r["analysis_id"] for r in led}
        missing = [k for k in confirm if k not in done]
        if missing or not confirm:
            raise FirewallError(f"{entry['direction']} lockbox read before all confirmatory reads: missing {missing}")
        not_closed = [k for k in confirm if allent[k]["status"] != "closed"]
        for k in confirm:
            try:
                registry.committed_blob(proj.root, proj.registry_dir / f"{k}.yaml")
            except RegistryError:
                not_closed.append(k)
        if not_closed:
            raise FirewallError(f"{entry['direction']} lockbox read before every confirmation is closed (committed "
                                f"status: closed): {sorted(set(not_closed))}")
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
    if not genes:                                               # zero-target close: ledger line only, nothing read
        return pd.DataFrame()
    Lm = _model_list(proj, "lockbox_models.txt")
    df = proj.outcome_loader(Lm, genes, source)
    if not set(df.index) <= set(Lm):
        raise FirewallError("outcome loader returned non-lockbox models")
    if persist:
        _persist_read(proj, analysis_id, source, df)
    return df


def _reads_dir(proj: Project) -> Path:
    return proj.cache_dir / "lockbox_reads"


def _reads_sidecar(proj: Project) -> Path:
    return proj.results_dir / "lockbox_reads.jsonl"


def _persist_read(proj: Project, analysis_id: str, source: str, df: pd.DataFrame) -> None:
    """Store a read matrix (gitignored) and commit its SHA-256 at once (PI crash-recovery decision, review MUST-2)."""
    d = _reads_dir(proj)
    d.mkdir(parents=True, exist_ok=True)
    p = d / f"{analysis_id}__{source}.parquet"
    if p.exists():
        raise FirewallError(f"{proj.rel(p)} already exists; a read is persisted once")
    df.to_parquet(p)
    rec = dict(time=time.strftime("%Y-%m-%dT%H:%M:%S%z"), analysis_id=analysis_id, source=source,
               file=proj.rel(p), sha256=_sha_file(p), n_models=int(df.shape[0]), n_genes=int(df.shape[1]))
    sc = _reads_sidecar(proj)
    _append_jsonl(sc, rec)
    rel = proj.rel(sc)
    r1 = git(proj.root, "add", "--", rel, check=False)
    r2 = git(proj.root, "commit", "-q", "-m", f"lockbox read persisted: {analysis_id} ({source}) sha256 "
             f"{rec['sha256'][:16]}", "--", rel, check=False) if r1.returncode == 0 else r1
    if r2.returncode != 0:
        raise FreezeError(f"read persisted at {proj.rel(p)} but its hash could not be committed: {r2.stderr.strip()}")


RELOOK_README = "results/part2/lockbox_confirm/README.md"
RELOOK_INTENT = "results/part2/lockbox_confirm/CORRECTED_REFIT_INTENT.md"


def resume_lockbox_read(analysis_id: str, proj: Project | None = None, source: str = "crispr",
                        relook: str | None = None) -> pd.DataFrame:
    """The persisted matrix of an entry that already read the lockbox (committed ledger line and committed hash);
    the outcome loader is never called. Refused if the hash differs or the entry is already closed -- unless `relook`
    names a PI ruling (e.g. 'R6') whose note is committed verbatim in RELOOK_README and whose committed RELOOK_INTENT
    names it: a disclosed, non-confirmatory relook (no ledger line; it is not a read)."""
    proj = proj or Project()
    led = [json.loads(l) for l in _committed_text(proj, proj.ledger).splitlines() if l.strip()] \
        if proj.ledger.exists() else []
    if not any(r["analysis_id"] == analysis_id for r in led):
        raise FirewallError(f"{analysis_id}: no committed ledger line; nothing to resume")
    entry, _ = registry.load_entry(proj.registry_dir, analysis_id, root=proj.root)
    if entry["status"] == "closed":
        if relook is None:
            raise FirewallError(f"{analysis_id} is closed; nothing to resume")
        note, intent = proj.root / RELOOK_README, proj.root / RELOOK_INTENT
        for p in (note, intent):
            try:
                registry.committed_blob(proj.root, p)
            except RegistryError as e:
                raise FirewallError(f"relook {relook}: {proj.rel(p)} must be committed ({e})") from e
        if f'"{relook} (PI' not in note.read_text() or relook not in intent.read_text():
            raise FirewallError(f"relook {relook}: no committed PI note / intent naming it")
    sc = _reads_sidecar(proj)
    recs = [json.loads(l) for l in _committed_text(proj, sc).splitlines() if l.strip()] if sc.exists() else []
    recs = [r for r in recs if r["analysis_id"] == analysis_id and r["source"] == source]
    if len(recs) != 1:
        raise FirewallError(f"{analysis_id}: {len(recs)} committed persisted-read records (need exactly 1)")
    p = proj.root / recs[0]["file"]
    if not p.exists() or _sha_file(p) != recs[0]["sha256"]:
        raise FirewallError(f"{analysis_id}: persisted read missing or its SHA-256 differs from the committed hash")
    return pd.read_parquet(p)
