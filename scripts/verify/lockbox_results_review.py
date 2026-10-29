#!/usr/bin/env python
"""Independent review of the lockbox read results (results/part2/verify/lockbox_results_review.md).

Independent of src/part2/confirm.py, robust.py, stats.py and assoc._design: the L design (covariate one-hot
assembly, FE-merge rule, span), HC3 with Satterthwaite df, the PatientID-blocked Freedman-Lane permutation and the
confirmation rule are re-implemented here (numpy; statsmodels HC3 as a cross-check). Re-used from the project: the
registry parsing of covariate names (covariates.covariates_for), the frozen score constructions (Spec.score_L,
cscore / h1 / controls), the label tables, and the persisted L reads.

L ACCESS (PI-authorised): only assoc.resume_lockbox_read(aid, relook="R6"), the driver's R6/R7 resume path. It checks
the committed ledger line, the committed SHA-256 and the committed R6 note + intent; it never calls the loader and
writes no ledger line. The project's outcome loader is replaced by a guard that refuses any non-A CRISPR request, and
the A matrix is materialised from a read-only memmap of A rows only. No L value is written to disk except the
per-hit L statistics of the sampled hits (the same fields the committed confirmation.tsv already carries).

Run (from the worktree root):
  export PART1_DIR=".../coessentiality_protein complexes"; export PYTHONHASHSEED=0
  .venv/bin/python scripts/verify/lockbox_results_review.py --out results/part2/verify/lockbox_results_review
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sst

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from part2 import assoc, covariates as covmod, data, paths  # noqa: E402

OUTC = ROOT / "results" / "part2" / "lockbox_confirm"
SEED = 20261029
N_PERM = 10000
LOADER_CALLS: list = []


# ------------------------------------------------------------------------------------------------ guarded A loading
def _mmap_gene_effect(use_cache: bool = True):
    csv = paths.GENE_EFFECT_CSV
    h = data.sha256_file(csv)
    npy = paths.CACHE / f"CRISPRGeneEffect.{h[:16]}.f32.npy"
    m = json.loads(npy.with_suffix(".json").read_text())
    return data.GeneEffect(np.load(npy, mmap_mode="r"), pd.Index(m["models"]), pd.Index(m["columns"]), h)


def guarded_loader(A: set):
    base = assoc._default_outcome_loader

    def f(model_ids, genes, source):
        ids = list(model_ids)
        LOADER_CALLS.append(dict(source=source, n=len(ids), all_A=set(ids) <= A))
        if not set(ids) <= A:
            raise assoc.FirewallError("review guard: the outcome loader was asked for non-A models")
        return base(ids, genes, source)
    return f


def load_driver():
    spec = importlib.util.spec_from_file_location("drv51", ROOT / "scripts" / "51_lockbox_confirm.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules["drv51"] = m
    spec.loader.exec_module(m)
    return m


# ------------------------------------------------------------------------------------------------ my design
FE = {"lineage_FE": "OncotreeLineage", "ModelType": "ModelType"}
PREFIX = {"library": "Library", "GrowthPattern": "GrowthPattern"}
ALIAS = {"NNMD": "ScreenNNMD", "doubling_time": "ScreenDoublingTime"}


def my_covariates(names, src: pd.DataFrame, assembly_ids, rows) -> pd.DataFrame:
    """Numeric covariate columns on `rows`. FE one-hot: reference = first sorted level among `assembly_ids`
    (assembly_ids = all models reproduces the executed (defect) design; = L rows read is the corrected design)."""
    out = {}
    for n in names:
        if n in FE:
            col = FE[n]
            levels = sorted(src.loc[list(assembly_ids), col].dropna().unique())
            v = src.loc[list(rows), col]
            for lev in levels[1:]:
                out[f"FE|{n}:{lev}"] = np.where(v.isna(), np.nan, (v == lev).astype(float))
        elif n in PREFIX:
            cols = sorted(c for c in src.columns if c.startswith(PREFIX[n]))
            if n == "GrowthPattern" and len(cols) > 1:
                cols = cols[1:]
            for c in cols:
                out[f"FE|{n}:{c}"] = src.loc[list(rows), c].astype(float).values
        else:
            out[n] = src.loc[list(rows), ALIAS.get(n, n)].astype(float).values
    return pd.DataFrame(out, index=list(rows))


def build_design(ctx: pd.Series, cov: pd.DataFrame):
    """Registered conventions: complete cases over context + covariates; constant columns dropped; indicator columns
    with < 2 models on either side merged into the reference (dropped). Returns (x, Z, keep mask, names)."""
    D = pd.concat([ctx.astype(float).rename("__ctx"), cov], axis=1)
    ok = D.notna().all(axis=1).values
    D = D[ok]
    n = len(D)
    cols = []
    for c in D.columns[1:]:
        v = D[c].values
        if np.unique(v).size <= 1:
            continue
        if c.startswith("FE|") and min(v.sum(), n - v.sum()) < 2:
            continue
        cols.append(c)
    Z = np.column_stack([np.ones(n)] + [D[c].values for c in cols]) if n else np.zeros((0, 1))
    return D["__ctx"].values, Z, ok, cols


def _basis(Z, tol=1e-10):
    if Z.shape[0] == 0:
        return Z
    U, s, _ = np.linalg.svd(Z, full_matrices=False)
    r = int((s > tol * max(1.0, s[0])).sum())
    return U[:, :r]


def my_fit(y: np.ndarray, x: np.ndarray, Z: np.ndarray, alt: str):
    """OLS of y on [Z, x] by Frisch-Waugh; HC3 SE of x; Satterthwaite df of the HC3 variance (exact quadratic-form
    moments under homoskedastic normal errors)."""
    n = len(y)
    nc = int(x.sum())
    out = dict(n=n, n_context=nc, n_ref=n - nc, beta=np.nan, se=np.nan, t=np.nan, df=np.nan, p=np.nan, ok=False,
               reason="")
    if nc < 2 or n - nc < 2:
        out["reason"] = "< 2 context or reference"
        return out, None
    Qz = _basis(Z)
    xt = x - Qz @ (Qz.T @ x)
    sxx = float(xt @ xt)
    if sxx < 1e-10 * max(1.0, float(x @ x)):
        out["reason"] = "context in covariate span"
        return out, None
    h = (Qz ** 2).sum(1) + xt ** 2 / sxx
    if np.any(h > 1 - 1e-10):
        out["reason"] = "leverage 1"
        return out, None
    c = xt / sxx
    beta = float(c @ y)
    M = np.eye(n) - Qz @ Qz.T - np.outer(xt, xt) / sxx
    e = M @ y
    w2 = (c / (1 - h)) ** 2
    se = float(np.sqrt((w2 * e ** 2).sum()))
    trA = float((w2 * np.diag(M)).sum())
    trA2 = float(w2 @ (M ** 2) @ w2)
    df = trA ** 2 / trA2
    t = beta / se
    p = float(sst.t.sf(t, df) if alt == "greater" else sst.t.cdf(t, df))
    out.update(beta=beta, se=se, t=t, df=df, p=p, ok=True)
    return out, dict(Qz=Qz, M=M, w2=w2, c=c)


def sm_check(y, x, Z):
    import statsmodels.api as sm
    Qz = _basis(Z)
    r = sm.OLS(y, np.column_stack([x, Qz])).fit(cov_type="HC3")
    return float(r.params[0]), float(r.bse[0])


def my_perm_p(y, x, piece, t_obs, alt, blocks, strata, rng, B=N_PERM):
    """Freedman-Lane: residuals of the reduced (covariates-only) model permuted by exchanging whole PatientID blocks of
    equal size within a stratum; HC3 t recomputed; one-sided p = (1 + #{s t* >= s t_obs}) / (B + 1)."""
    Qz, M, w2, c = piece["Qz"], piece["M"], piece["w2"], piece["c"]
    f = Qz @ (Qz.T @ y)
    r = y - f
    n = len(y)
    bl = np.asarray(blocks, dtype=object)
    bl = np.array([f"__{i}" if (b is None or (isinstance(b, float) and np.isnan(b))) else str(b)
                   for i, b in enumerate(bl)])
    st = np.zeros(n, dtype=object) if strata is None else np.asarray(strata, dtype=object)
    st = np.array(["NA" if (s is None or (isinstance(s, float) and np.isnan(s))) else str(s) for s in st])
    pos = {}
    for i, b in enumerate(bl):
        pos.setdefault(b, []).append(i)
    groups = {}
    for b, ix in pos.items():
        groups.setdefault((st[ix[0]], len(ix)), []).append(np.array(ix))
    groups = [np.stack(v) for v in groups.values() if len(v) > 1]
    P = np.tile(np.arange(n), (B, 1))
    for G in groups:
        for k in range(B):
            P[k, G.ravel()] = G[rng.permutation(G.shape[0])].ravel()
    Ys = f[:, None] + r[P].T
    bs = c @ Ys
    E = M @ Ys
    ts = bs / np.sqrt((w2[:, None] * E ** 2).sum(0))
    s = 1.0 if alt == "greater" else -1.0
    return float((1 + np.sum(s * ts >= s * t_obs - 1e-12)) / (B + 1))


def test_hit(y: pd.Series, ctx: pd.Series, names, src, assembly_ids, beta_A, blocks, strata, rng, perm=True):
    L = y.index
    cov = my_covariates(names, src, assembly_ids, L)
    x, Z, ok, cols = build_design(ctx.reindex(L), cov)
    yy = y.values.astype(float)[ok]
    m = ~np.isnan(yy)
    alt = "less" if beta_A < 0 else "greater"
    fit, piece = my_fit(yy[m], x[m], Z[m], alt)
    fit["p_perm"] = np.nan
    fit["p_final"] = fit["p"]
    fit["n_y_nan_in_cc"] = int((~m).sum())
    if fit["ok"]:
        b2, se2 = sm_check(yy[m], x[m], Z[m])
        fit["sm_beta"], fit["sm_se"] = b2, se2
        if fit["n_context"] < 20 and perm:
            ids = L[ok][m]
            fit["p_perm"] = my_perm_p(yy[m], x[m], piece, fit["t"], alt, blocks.reindex(ids).values,
                                      None if strata is None else strata.reindex(ids).values, rng)
            fit["p_final"] = fit["p_perm"]
        h = sst.t.ppf(0.975, fit["df"]) * fit["se"]
        fit["ci95"] = [fit["beta"] - h, fit["beta"] + h]
    testable = fit["n_context"] >= 5 and fit["n_ref"] >= 5
    same = bool(fit["ok"] and np.sign(fit["beta"]) == np.sign(beta_A))
    fit["testable"] = bool(testable)
    fit["confirmed"] = bool(fit["ok"] and testable and same and fit["p_final"] < 0.05)
    fit["call"] = ("NOT EVALUABLE" if not fit["ok"] else ("confirmed" if fit["confirmed"] else "not confirmed"))
    fit["n_cov_cols"] = len(cols)
    fit["rank_Z"] = int(_basis(Z[m]).shape[1]) if fit["ok"] or len(Z[m]) else 0
    return fit, (x[m], Z[m])


def span_relation(Zd, Zc, x):
    A, B = np.column_stack([x, Zd]), np.column_stack([x, Zc])
    ra, rb = np.linalg.matrix_rank(A), np.linalg.matrix_rank(B)
    rab = np.linalg.matrix_rank(np.column_stack([A, B]))
    return dict(rank_defect=int(ra), rank_corrected=int(rb), same_span=bool(ra == rb == rab))


# ------------------------------------------------------------------------------------------------ main
def committed_rows(aid):
    t = pd.read_csv(OUTC / aid / "confirmation.tsv", sep="\t")
    for c in ("eligible", "lockbox_testable", "confirmed_final"):
        t[c] = t[c].astype(str).str.lower() == "true"
    return t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--n-perm", type=int, default=N_PERM)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    data.load_gene_effect = _mmap_gene_effect            # A rows only are materialised (memmap fancy indexing)
    drv = load_driver()
    Aset = set(assoc.a_models(assoc.Project()))
    # the project's default Project() is built inside real_env; guard every loader call it makes
    assoc.Project.__dataclass_fields__["outcome_loader"].default = guarded_loader(Aset)
    orig_init = assoc.Project.__init__

    def init(self, root=paths.PROJECT, outcome_loader=None):
        orig_init(self, root=root, outcome_loader=outcome_loader or guarded_loader(Aset))
    assoc.Project.__init__ = init
    world = drv.real_world()
    world.proj.outcome_loader = guarded_loader(Aset)
    print(f"[{time.time()-t0:.0f}s] world built; loader calls so far: {LOADER_CALLS}", flush=True)
    src = world.cov
    allids = src.index
    corr = pd.read_csv(OUTC / "corrected_refit_nonconfirmatory.tsv", sep="\t")
    ckey = corr.set_index(["entry", "hit_index"])
    rng_s = np.random.default_rng(SEED)
    rows, spans = [], []
    access = []

    def run(aid, idx, spec, E_L, tag):
        e, _ = drv.entry(world, aid)
        com = committed_rows(aid)
        for i in idx:
            r = spec.hits.iloc[i]
            cr = com.iloc[i]
            assert (cr.context, cr.unit) == (r.context, r.unit)
            names = covmod.covariates_for(e, r.context)[0]
            y = spec.score_L(E_L, r.context, r.unit)
            st = None if r.context in spec.lineage_contexts else world.cov.OncotreeLineage
            rng = np.random.default_rng([SEED, len(rows)])
            fd, (xd, Zd) = test_hit(y, spec.ctx_all[r.context], names, src, allids, float(r.beta),
                                    world.cov.PatientID, st, rng)
            fc, (xc, Zc) = test_hit(y, spec.ctx_all[r.context], names, src, E_L.index, float(r.beta),
                                    world.cov.PatientID, st, np.random.default_rng([SEED, len(rows), 1]))
            sp = span_relation(Zd, Zc, xd) if len(xd) == len(xc) and len(xd) else {}
            k = (aid, i)
            c = ckey.loc[k] if k in ckey.index else None
            row = dict(tag=tag, entry=aid, hit_index=i, context=r.context, unit=r.unit, beta_A=float(r.beta),
                       eligible=bool(cr.eligible), route=cr.route, names=";".join(names),
                       lineage_context=r.context in spec.lineage_contexts,
                       prod_beta=cr.beta_L, prod_se=cr.se_L, prod_df=cr.df_L, prod_p=cr.p_L, prod_p_perm=cr.p_L_perm,
                       prod_nctx=cr.n_context_L, prod_confirmed=bool(cr.confirmed_final),
                       prod_status=cr.status,
                       corr_cause="" if c is None else c.affected_cause,
                       corr_beta=cr.beta_L if c is None else c.beta_L_relook,
                       corr_se=cr.se_L if c is None else c.se_L_relook,
                       corr_p=cr.p_L if c is None else c.p_L_relook,
                       corr_confirmed=bool(cr.confirmed_final) if c is None else
                       str(c.confirmed_relook).lower() == "true",
                       **{f"d_{k2}": fd.get(k2) for k2 in ("ok", "reason", "n", "n_context", "beta", "se", "df", "p",
                                                           "p_perm", "p_final", "sm_se", "testable", "confirmed",
                                                           "call", "n_cov_cols", "n_y_nan_in_cc")},
                       **{f"c_{k2}": fc.get(k2) for k2 in ("ok", "reason", "n", "n_context", "beta", "se", "df", "p",
                                                           "p_perm", "p_final", "sm_se", "testable", "confirmed",
                                                           "call", "n_cov_cols")},
                       **sp)
            rows.append(row)

    def resume(aid):
        n0 = len(LOADER_CALLS)
        E = assoc.resume_lockbox_read(aid, proj=world.proj, relook="R6")
        access.append(dict(entry=aid, call="assoc.resume_lockbox_read(relook='R6')", n_models=int(E.shape[0]),
                           n_genes=int(E.shape[1]), loader_calls_during=len(LOADER_CALLS) - n0,
                           all_L=not (set(E.index) & Aset)))
        return E

    # ---- H3a / H3b samples
    for aid, k in (("H3a_lineage_atlas", 20), ("H3b_driver_family", 15)):
        spec = drv.build_spec(world, aid)
        com = committed_rows(aid)
        pool = np.flatnonzero((com.eligible & com.lockbox_testable).values)
        idx = sorted(rng_s.choice(pool, size=k, replace=False).tolist())
        if aid == "H3b_driver_family":
            # defect diagnosis: one failed fit (WGD Z0026) and one shifted fit (TP53_altered Z0096) added explicitly
            for c_, u_ in (("WGD", "Z0026"), ("TP53_altered", "Z0096")):
                j = int(np.flatnonzero(((com.context == c_) & (com.unit == u_)).values)[0])
                if j not in idx:
                    idx.append(j)
        E_L = resume(aid)
        run(aid, idx, spec, E_L, "sample" if aid.startswith("H3a") else "sample")
        del E_L
        print(f"[{time.time()-t0:.0f}s] {aid} done", flush=True)
    # ---- H1 eligible (5)
    aid = "H1_within_complex_GI"
    spec = drv.build_spec(world, aid)
    com = committed_rows(aid)
    idx = np.flatnonzero(com.eligible.values).tolist()
    E_L = resume(aid)
    run(aid, idx, spec, E_L, "H1")
    # EP300 rare route on A (LOSO, LOMO, DEMETER2 same sign), my fits
    e, _ = drv.entry(world, aid)
    rr = []
    for i in idx:
        r = spec.hits.iloc[i]
        if com.iloc[i].route == "lockbox":
            continue
        names = covmod.covariates_for(e, r.context)[0]
        A = world.A
        yA = spec.score_A(r.context, r.unit).reindex(A)
        ctxA = spec.ctx_all[r.context].reindex(A)
        lin = world.cov.OncotreeLineage.reindex(A)
        sign = np.sign(float(r.beta))

        def bsign(keep):
            fit, _ = test_hit(yA[keep], ctxA, names, src, A, float(r.beta), world.cov.PatientID, None, None,
                              perm=False)
            return fit
        full = bsign(np.ones(len(A), bool))
        loso = [bsign((lin != lev).values) for lev in sorted(lin.dropna().unique())]
        loso_ok = [np.sign(f["beta"]) == sign for f in loso if f["ok"]]
        pos = np.flatnonzero((ctxA == 1).values & yA.notna().values)
        lomo = []
        for j in range(len(A)):
            keep = np.ones(len(A), bool)
            keep[j] = False
            if ctxA.iloc[j] != 1:            # dropping a reference model cannot flip much; still fitted below
                continue
            f = bsign(keep)
            lomo.append(f["ok"] and np.sign(f["beta"]) == sign)
        d2 = spec.score_D2(r.context, r.unit)
        # DEMETER2 rows = the A models DEMETER2 screens (not reindexed to A: the FE merge counts x-complete rows)
        fD, _ = test_hit(d2, ctxA, names, src, A, float(r.beta), world.cov.PatientID, None, None,
                         perm=False) if d2 is not None else (dict(ok=False), None)
        rr.append(dict(unit=r.unit, full_beta_A=full["beta"], full_ok=full["ok"], loso_folds=len(loso_ok),
                       loso_same_sign=bool(loso_ok and all(loso_ok)), lomo_context_models=len(lomo),
                       lomo_same_sign=bool(lomo and all(lomo)), demeter2_beta=fD.get("beta"),
                       demeter2_same_sign=bool(fD.get("ok") and np.sign(fD["beta"]) == sign),
                       prod_loso=com.iloc[i].loso_same_sign, prod_lomo=com.iloc[i].lomo_same_sign,
                       prod_d2=com.iloc[i].demeter2_same_sign, prod_pass=com.iloc[i].rare_route_pass,
                       n_context_A=int(len(pos))))
    pd.DataFrame(rr).to_csv(a.out / "h1_rare_route.tsv", sep="\t", index=False)
    del E_L
    print(f"[{time.time()-t0:.0f}s] H1 done", flush=True)
    # ---- controls (5 entries, every unit)
    for aid in ("CTRL_PC1", "CTRL_PC5", "CTRL_PC9", "CTRL_PC17", "CTRL_PC21"):
        spec = drv.build_spec(world, aid)
        com = committed_rows(aid)
        idx = np.flatnonzero((com.eligible & com.lockbox_testable).values).tolist()
        if not idx:
            continue
        E_L = resume(aid)
        run(aid, idx, spec, E_L, "controls")
        del E_L
    print(f"[{time.time()-t0:.0f}s] controls done", flush=True)
    # ---- span audit over every fit of every entry (labels + the L row index only): which fits does the defect touch?
    span_rows = []
    order = drv.confirm_order(world.proj) + [drv.BRIDGE_ID]
    for aid in order:
        e, _ = drv.entry(world, aid)
        com = committed_rows(aid)
        Lidx = resume(aid).index
        h1 = aid in drv.H1_IDS
        ctab = world.loss if h1 else world.ctx
        cache = {}
        for i, r in enumerate(com.itertuples(index=False)):
            if r.context not in cache:
                names = covmod.covariates_for(e, r.context)[0]
                ctx = ctab[r.context].reindex(Lidx)
                out = {}
                for tag, asm in (("d", allids), ("c", Lidx)):
                    x, Z, ok, cols = build_design(ctx, my_covariates(names, src, asm, Lidx))
                    Qz = _basis(Z)
                    xt = x - Qz @ (Qz.T @ x)
                    sxx = float(xt @ xt)
                    lev = (Qz ** 2).sum(1) + (xt ** 2 / sxx if sxx > 1e-10 else 0)
                    out[tag] = (x, Z, bool(np.any(lev > 1 - 1e-10)))
                sp = span_relation(out["d"][1], out["c"][1], out["d"][0])
                cache[r.context] = dict(same_span=sp["same_span"], lev1_defect=out["d"][2], lev1_corrected=out["c"][2])
            span_rows.append(dict(entry=aid, hit_index=i, context=r.context, unit=r.unit, eligible=bool(r.eligible),
                                  route=r.route, prod_fitted=bool(np.isfinite(r.beta_L)), **cache[r.context]))
    resume(drv.BENCH_ID)                   # committed-hash check of the 30th persisted read (values unused)
    pd.DataFrame(span_rows).to_csv(a.out / "span_audit.tsv", sep="\t", index=False)
    print(f"[{time.time()-t0:.0f}s] span audit done", flush=True)
    res = pd.DataFrame(rows)
    res.to_csv(a.out / "recompute_calls.tsv", sep="\t", index=False, float_format="%.8g")
    pd.DataFrame(access).to_csv(a.out / "l_access_log.tsv", sep="\t", index=False)
    json.dump(dict(loader_calls=LOADER_CALLS, runtime_s=round(time.time() - t0, 1)),
              open(a.out / "run_meta.json", "w"), indent=1, default=str)
    print(f"done {time.time()-t0:.0f}s; loader calls {LOADER_CALLS}")


if __name__ == "__main__":
    main()
