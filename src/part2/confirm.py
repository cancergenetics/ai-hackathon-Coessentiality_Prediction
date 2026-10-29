"""Lockbox confirmation step (DESIGN v3.1 §0 F4/F5, §2.5; registry `lockbox` blocks).

Order for one registry entry (enforced here and in assoc.read_lockbox):
  1. ``write_target_list``: the target list (genes = members of every discovered unit, plus the discovered hits with
     their A sign) is written and COMMITTED before any lockbox outcome is touched. No discovery -> a zero-target list.
  2. ``confirm_entry``: assoc.read_lockbox (ledger line committed first, one read per entry) -> per hit, the same
     score / removals / covariates in L -> one-sided test in the A direction (HC3, Satterthwaite df; exact
     PatientID-blocked Freedman-Lane permutation, within lineage for non-lineage contexts, when < 20 context models
     are in L) -> confirmed iff same sign and one-sided p < 0.05 and >= 5 context and >= 5 reference L models;
     otherwise 'not lockbox-testable' (rare-context route computed in A by ``rare_route``). L effects with 95% CIs
     are reported for every hit. Results -> results/part2/assoc/<id>/confirmation.json.
  3. ``close_entry``: status open -> closed, committed (only after the entry's ledger line exists and, for a read
     with targets, after its confirmation output is committed -- an entry whose read crashed cannot be closed by hand).
BENCH / BRIDGE-L reads are refused by read_lockbox until every confirmatory entry is closed.
Review MUST-1/2 (results/part2/verify/lockbox_code_review.md): the target list records each hit's eligibility and
route; confirm_entry refuses hits that differ from the committed target list; persist=True stores the read matrix
(assoc.read_lockbox) and resume=True re-uses it (assoc.resume_lockbox_read, never the loader).
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sst

from part2 import assoc, registry, stats
from part2.registry import RegistryError, git
from part2.robust import fit_one, leave_one_model_out, loso


def commit_paths(proj, paths_, msg: str) -> None:
    rels = [proj.rel(p) for p in paths_]
    r1 = git(proj.root, "add", "--", *rels, check=False)
    r2 = git(proj.root, "commit", "-q", "-m", msg, "--", *rels, check=False) if r1.returncode == 0 else r1
    if r2.returncode != 0:
        raise RegistryError(f"commit failed: {r2.stderr.strip() or r2.stdout.strip()}")


def target_path(proj, analysis_id: str) -> Path:
    return proj.results_dir / "lockbox_targets" / f"{analysis_id}.json"


HIT_FIELDS = ("eligible", "route", "label_final", "headline_eligible", "event", "axis")


def _jsonable(v):
    if isinstance(v, (bool, np.bool_)):
        return bool(v)
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (float, np.floating)):
        return None if not np.isfinite(v) else float(v)
    return None if v is None else str(v)


def hit_keys(hits) -> list:
    """(context, unit, beta_A, alternative) of each hit, in order: the hypothesis identity checked at the read."""
    if hits is None or not len(hits):
        return []
    return [(str(r.context), str(r.unit), float(r.beta), str(r.alternative)) for r in hits.itertuples(index=False)]


def write_target_list(proj, analysis_id: str, genes, hits: pd.DataFrame | None = None, units: dict | None = None,
                      contexts=None, commit: bool = True, msg_suffix: str = "") -> Path:
    """Write (and commit) the target list. Refuses to overwrite an existing list (one list per entry). Each hit
    records (context, unit, beta_A, alternative) and, when the hits table carries them, its eligibility, route,
    E2 label and headline flag (HIT_FIELDS)."""
    p = target_path(proj, analysis_id)
    if p.exists():
        raise RegistryError(f"target list {proj.rel(p)} already exists (one target list per entry)")
    genes = sorted(set(map(str, genes)))
    rec = dict(analysis_id=analysis_id, genes=genes, units={k: list(v) for k, v in (units or {}).items()},
               contexts=sorted(set(contexts or [])))
    if hits is not None and len(hits):
        extra = [f for f in HIT_FIELDS if f in hits.columns]
        rec["hits"] = [dict(context=str(r["context"]), unit=str(r["unit"]), beta_A=float(r["beta"]),
                            alternative=str(r["alternative"]), **{f: _jsonable(r[f]) for f in extra})
                       for r in hits.to_dict("records")]
    if not genes:
        rec["zero_targets"] = True
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(rec, indent=1, sort_keys=True) + "\n")
    if commit:
        commit_paths(proj, [p], f"lockbox targets: {analysis_id}{msg_suffix}")
    return p


def lockbox_test(y: pd.Series, ctx: pd.Series, cov: pd.DataFrame | None, beta_A: float, *, blocks=None, strata=None,
                 n_perm: int = 10000, seed: int = 0, perm_below: int = 20, min_n: int = 5) -> dict:
    """One-sided test in the A direction on lockbox rows."""
    from part2.assoc import _design, permutation_p
    alt = "less" if beta_A < 0 else "greater"
    f = fit_one(y, ctx, cov, alt)
    out = dict(beta_L=f["beta"], se_L=f["se"], df_L=f["df"], n_L=f["n"], n_context_L=f["n_context"],
               p_L_hc3=f["p"], p_L=f["p"], p_L_perm=np.nan)
    out["n_ref_L"] = f["n"] - f["n_context"]
    out["testable"] = bool(f["n_context"] >= min_n and out["n_ref_L"] >= min_n)
    if f["ok"]:
        h = sst.t.ppf(0.975, f["df"]) * f["se"]
        out["ci95_L"] = [f["beta"] - h, f["beta"] + h]
        if f["n_context"] < perm_below and n_perm > 0:
            idx = y.index
            X, _, ok = _design(ctx.reindex(idx), None if cov is None else cov.reindex(idx))
            yy = y.values[ok]
            m = ~np.isnan(yy)
            ids = idx[ok][m]
            st = None if strata is None else strata.reindex(ids).fillna("NA").astype(str).values
            out["p_L_perm"] = permutation_p(yy[m], X[m], ids, alt, n_perm, np.random.default_rng(seed), blocks, st)
            out["p_L"] = out["p_L_perm"]
    else:
        out["ci95_L"] = [np.nan, np.nan]
    out["same_sign"] = bool(f["ok"] and np.sign(f["beta"]) == np.sign(beta_A))
    out["confirmed"] = bool(out["testable"] and out["same_sign"] and out["p_L"] < 0.05)
    return out


def rare_route(y_A, ctx_A, cov_A, groups_A, beta_A: float, y_D2_A: pd.Series | None = None) -> dict:
    """Rare-context route in A: LOSO same sign, DEMETER2 (A models) same sign, leave-one-model-out same sign."""
    sign = np.sign(beta_A)
    lo = loso(y_A, ctx_A, cov_A, groups_A, "two-sided")
    # LOSO same sign only (registry rare route): re-count without the p requirement
    g = groups_A.reindex(y_A.index)
    same = []
    for lev in sorted(g.dropna().unique()):
        f = fit_one(y_A[(g != lev).values], ctx_A, cov_A)
        if f["ok"]:
            same.append(np.sign(f["beta"]) == sign)
    lomo = leave_one_model_out(y_A, ctx_A, cov_A)
    d2 = fit_one(y_D2_A, ctx_A, cov_A) if y_D2_A is not None else None
    d2_same = None if d2 is None or not d2["ok"] else bool(np.sign(d2["beta"]) == sign)
    ok = bool(same and all(same)) and lomo["lomo_same_sign"] and bool(d2_same)
    return dict(route="not lockbox-testable", loso_same_sign=bool(same and all(same)), loso_folds=len(same),
                lomo_same_sign=lomo["lomo_same_sign"], demeter2_same_sign=d2_same, rare_route_pass=ok,
                loso_full=lo)


def confirm_entry(analysis_id: str, hits: pd.DataFrame, score_L, contexts_L: pd.DataFrame, covariates_L,
                  *, proj=None, blocks=None, strata=None, lineage_contexts=(), n_perm: int = 10000,
                  seed: int = 0, source: str = "crispr", write: bool = True, persist: bool = False,
                  resume: bool = False) -> pd.DataFrame:
    """Read the lockbox for the committed target list and test every hit.
    hits: columns context, unit, beta (A). score_L(E_L, context, unit) -> Series over L models (same score, removals,
    A-fitted parameters). contexts_L / covariates_L: L labels (DataFrame or {context: DataFrame})."""
    proj = proj or assoc.Project()
    tp = target_path(proj, analysis_id)
    try:
        registry.committed_blob(proj.root, tp)
    except RegistryError as e:
        raise assoc.FreezeError(f"target list must be committed before a lockbox read: {e}") from e
    want = [(h["context"], h["unit"], float(h["beta_A"]), h["alternative"])
            for h in json.loads(tp.read_text()).get("hits", [])]
    if hit_keys(hits) != want:
        raise RegistryError(f"{analysis_id}: the hits differ from the committed target list (nothing read)")
    if resume:
        E_L = assoc.resume_lockbox_read(analysis_id, proj=proj, source=source)
    else:
        E_L = assoc.read_lockbox(analysis_id, tp, proj=proj, source=source, persist=persist)
    rows = []
    for i, r in enumerate(hits.itertuples(index=False)):
        if E_L.empty:
            break
        y = score_L(E_L, r.context, r.unit)
        cov = covariates_L.get(r.context) if isinstance(covariates_L, dict) else covariates_L
        st = None if r.context in set(lineage_contexts) else strata
        t = lockbox_test(y, contexts_L[r.context], cov, float(r.beta), blocks=blocks, strata=st, n_perm=n_perm,
                         seed=seed + i)
        rows.append(dict(context=r.context, unit=r.unit, beta_A=float(r.beta), **t))
    res = pd.DataFrame(rows)
    if write:
        d = proj.results_dir / "assoc" / analysis_id
        d.mkdir(parents=True, exist_ok=True)
        recs = json.loads(res.to_json(orient="records")) if len(res) else []
        summ = dict(analysis_id=analysis_id, source=source, n_hits=int(len(hits)), n_tested=int(len(res)),
                    n_testable=int(res.testable.sum()) if len(res) else 0,
                    n_confirmed=int(res.confirmed.sum()) if len(res) else 0, hits=recs)
        (d / f"confirmation_{source}.json").write_text(json.dumps(summ, indent=1, sort_keys=True) + "\n")
    return res


def close_entry(proj, analysis_id: str, outputs=None) -> None:
    """status: open -> closed in a commit that changes nothing else; needs the entry's committed ledger line and, when
    that read had targets, its committed confirmation output (default results/part2/assoc/<id>/confirmation_crispr.json;
    the driver names its BENCH / BRIDGE-L outputs). An entry whose read crashed before its output was committed can
    therefore not be closed by hand (review S6 / MUST-2): it is finished with the resume mode or by a PI ruling."""
    led = [json.loads(l) for l in assoc._committed_text(proj, proj.ledger).splitlines() if l.strip()]
    mine = [r for r in led if r["analysis_id"] == analysis_id]
    if not mine:
        raise RegistryError(f"{analysis_id} has no committed lockbox ledger line; read (or zero-target close) first")
    if any(int(r.get("n_genes", 0)) > 0 for r in mine):
        outs = outputs or [proj.results_dir / "assoc" / analysis_id / "confirmation_crispr.json"]
        for o in outs:
            try:
                registry.committed_blob(proj.root, o)
            except RegistryError as e:
                raise RegistryError(f"{analysis_id}: confirmation output {proj.rel(o)} is not committed; an entry "
                                    f"read without committed output cannot be closed (resume or PI ruling): {e}") from e
    p = proj.registry_dir / f"{analysis_id}.yaml"
    registry.committed_blob(proj.root, p)
    old = p.read_text()
    new, k = re.subn(r"(?m)^status: open$", "status: closed", old)
    if k != 1:
        raise RegistryError(f"{analysis_id}: expected exactly one 'status: open' line")
    e_old, e_new = registry.load_yaml(old), registry.load_yaml(new)
    if registry.validate(e_new, analysis_id):
        raise RegistryError(f"closed entry does not validate: {registry.validate(e_new, analysis_id)}")
    diff = [k for k in set(e_old) | set(e_new) if e_old.get(k) != e_new.get(k)]
    if diff != ["status"]:
        raise RegistryError(f"close changed more than status: {diff}")
    p.write_text(new)
    commit_paths(proj, [p], f"close registry entry {analysis_id}")


def is_closed(proj, analysis_id: str) -> bool:
    """True iff the committed entry is closed (used by gates and tests)."""
    e, _ = registry.load_entry(proj.registry_dir, analysis_id, root=proj.root, allow_pending=True)
    return e["status"] == "closed"

