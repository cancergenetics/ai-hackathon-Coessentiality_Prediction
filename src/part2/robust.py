"""Per-hit robustness (DESIGN v3.1 §2.4-2.5; reported, never selected) and negative-control pseudo-contexts.

* LOSO: leave one lineage out (one OncotreePrimaryDisease out for lineage contexts); robust iff every evaluable fold
  has the same sign as the full fit and p < 0.05 (the test's alternative). Folds leaving < 2 context or < 2 reference
  models, or a singular design, are skipped and counted.
* Leave-top-2-out: drop the 2 context models with the most extreme score in the direction of the effect; robust iff
  same sign and p < 0.05.
* Breadth: fraction of the unit's members (after removals) whose own gene-level test has one-sided p < 0.05 in the
  direction of the unit's effect.
* Leave-one-model-out (rare-context route): every context model left out in turn; same sign throughout.
* Pseudo-contexts: lineage-stratified, PatientID-blocked random relabelling of the A models with the real context's
  size per stratum (lineage contexts: random model sets of the same size).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from part2 import stats


def fit_one(y: pd.Series, ctx: pd.Series, cov: pd.DataFrame | None, alternative: str = "two-sided") -> dict:
    from part2.assoc import _design
    idx = y.index
    X, names, ok = _design(ctx.reindex(idx), None if cov is None else cov.reindex(idx))
    yy = y.values[ok]
    m = ~np.isnan(yy)
    n_ctx = int(X[m, 1].sum())
    out = dict(n=int(m.sum()), n_context=n_ctx, beta=np.nan, se=np.nan, t=np.nan, df=np.nan, p=np.nan, ok=False)
    if n_ctx < 2 or m.sum() - n_ctx < 2:
        return out
    try:
        d = stats.HC3Design(X[m], 1)
    except (np.linalg.LinAlgError, ValueError):
        return out
    b, se, t = d.fit(yy[m])
    out.update(beta=float(b), se=float(se), t=float(t), df=d.df, p=float(stats.pvalue_t(t, d.df, alternative)),
               ok=True)
    return out


def _passes(f: dict, sign: float) -> bool:
    return bool(f["ok"] and np.sign(f["beta"]) == sign and f["p"] < 0.05)


def loso(y, ctx, cov, groups: pd.Series, alternative: str = "two-sided") -> dict:
    full = fit_one(y, ctx, cov, alternative)
    sign = np.sign(full["beta"])
    g = groups.reindex(y.index)
    n_f = n_pass = n_skip = 0
    for lev in sorted(g.dropna().unique()):
        keep = (g != lev).values
        f = fit_one(y[keep], ctx, cov, alternative)
        if not f["ok"]:
            n_skip += 1
            continue
        n_f += 1
        n_pass += _passes(f, sign)
    return dict(loso_folds=n_f, loso_pass=n_pass, loso_skipped=n_skip, loso_robust=bool(n_f > 0 and n_pass == n_f))


def leave_top2_out(y, ctx, cov, alternative: str = "two-sided", k: int = 2) -> dict:
    full = fit_one(y, ctx, cov, alternative)
    sign = np.sign(full["beta"])
    c = ctx.reindex(y.index)
    pos = y[(c == 1) & y.notna()]
    drop = (pos.sort_values() if sign < 0 else pos.sort_values(ascending=False)).index[:k]
    f = fit_one(y.drop(drop), ctx, cov, alternative)
    return dict(top2_beta=f["beta"], top2_p=f["p"], top2_robust=_passes(f, sign))


def leave_one_model_out(y, ctx, cov, alternative: str = "two-sided") -> dict:
    full = fit_one(y, ctx, cov, alternative)
    sign = np.sign(full["beta"])
    c = ctx.reindex(y.index)
    ids = y.index[(c == 1).values & y.notna().values]
    same = [np.sign(fit_one(y.drop(i), ctx, cov, alternative)["beta"]) == sign for i in ids]
    return dict(lomo_n=len(same), lomo_same_sign=bool(len(same) and all(same)))


def breadth(E: pd.DataFrame, members, ctx, cov, sign: float) -> dict:
    alt = "less" if sign < 0 else "greater"
    ps = [fit_one(E[g].astype(float).clip(-3, 3), ctx, cov, alt)["p"] for g in members if g in E.columns]
    ps = np.array([p for p in ps if np.isfinite(p)])
    return dict(breadth=float(np.mean(ps < 0.05)) if len(ps) else np.nan, breadth_n=int(len(ps)))


def robustness(hits: pd.DataFrame, Y, contexts: pd.DataFrame, covariates, groups: dict, E: pd.DataFrame | None = None,
               members: dict | None = None) -> pd.DataFrame:
    """One row per hit (columns context, unit, alternative): LOSO, leave-top-2-out, breadth.
    Y / covariates: DataFrame or {context: DataFrame} as in run_association; groups: {context: Series} of the LOSO
    grouping (lineage; primary disease for lineage contexts)."""
    rows = []
    for r in hits.itertuples(index=False):
        Yc = Y[r.context] if isinstance(Y, dict) else Y
        Zc = covariates.get(r.context) if isinstance(covariates, dict) else covariates
        y, c = Yc[r.unit], contexts[r.context]
        full = fit_one(y, c, Zc, r.alternative)
        row = dict(context=r.context, unit=r.unit, beta=full["beta"], p=full["p"])
        row.update(loso(y, c, Zc, groups[r.context], r.alternative))
        row.update(leave_top2_out(y, c, Zc, r.alternative))
        if E is not None and members is not None and r.unit in members:
            row.update(breadth(E, members[r.unit], c, Zc, np.sign(full["beta"])))
        rows.append(row)
    return pd.DataFrame(rows)


def pseudo_contexts(ctx: pd.Series, n: int, seed: int, strata: pd.Series | None = None,
                    blocks: pd.Series | None = None) -> pd.DataFrame:
    """n pseudo-contexts matched to `ctx` in size per stratum (lineage composition), relabelling whole PatientID
    blocks where possible (a block is drawn whole; the per-stratum count is matched to within the block size)."""
    rng = np.random.default_rng(seed)
    c = ctx.dropna()
    st = (strata.reindex(c.index).fillna("NA") if strata is not None else pd.Series("all", index=c.index))
    bl = blocks.reindex(c.index) if blocks is not None else pd.Series(c.index, index=c.index)
    bl = bl.where(bl.notna(), pd.Series(c.index, index=c.index)).astype(str)
    out = {}
    for k in range(n):
        lab = pd.Series(0.0, index=c.index)
        for s in sorted(st.unique()):
            ids = c.index[st == s]
            need = int(c[ids].sum())
            if need == 0:
                continue
            blks = list(dict.fromkeys(bl[ids].values))
            got = 0
            for b in rng.permutation(len(blks)):
                if got >= need:
                    break
                mem = ids[bl[ids].values == blks[b]]
                lab[mem] = 1.0
                got += len(mem)
        out[f"pseudo{k}"] = lab.reindex(ctx.index)
    return pd.DataFrame(out)
