"""Benchmark and bridge-L estimands (DESIGN v3.1 §2.3, §2.5; BENCH_complex_vs_gene, BRIDGE-L_confirmation_rate).
Both read the lockbox only after every confirmatory entry is closed (enforced by assoc.read_lockbox).

BENCH: Delta(C, c) = z_L(S_C) - z_L(g*), z_L = HC3 t of the context coefficient in L signed so that positive = the
A direction; g* = member of C (after removals) with the largest A t in the A direction (ties: larger |beta|), chosen
in A only. Breadth terciles fixed on the primary pairs' A breadth. CI by PatientID-blocked bootstrap over L models
(B = 2000, seed 5100; pairs and g* fixed, t recomputed per resample).

BRIDGE-L: per map, confirmation rate (confirmed / discovered, lockbox-testable contexts) within size bins
{3-4, 5-8, 9-16, 17-60}, averaged with the pooled (both maps) discovered-cluster size distribution; difference
zavg - GLS with a bootstrap over discovered hits within map and bin (B = 2000, seed 5200).

BENCH secondary (al_correlation, effect_rsd, distinct_gstar): A-L correlation of context effects in residual-SD units
for complex scores vs matched random gene sets; the correlation method, R, pair set and summary are registry rulings
(no defaults). BRIDGE-L multiverse: bridge_unweighted (confirmed / discovered per map without size-bin weights).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from part2.robust import fit_one

SIZE_BINS = ((3, 4), (5, 8), (9, 16), (17, 60))


def choose_gstar(E_A: pd.DataFrame, members, ctx, cov, beta_A: float) -> str:
    s = np.sign(beta_A)
    best, key = None, None
    for g in members:
        f = fit_one(E_A[g].astype(float).clip(-3, 3), ctx, cov)
        if not f["ok"]:
            continue
        k = (s * f["t"], abs(f["beta"]))
        if key is None or k > key:
            best, key = g, k
    return best


def _signed_t(y, ctx, cov, sign) -> float:
    f = fit_one(y, ctx, cov)
    return sign * f["t"] if f["ok"] else np.nan


def bench_delta(pairs: pd.DataFrame, score_L, E_L: pd.DataFrame, contexts_L: pd.DataFrame, covariates_L,
                blocks_L: pd.Series | None = None, B: int = 2000, seed: int = 5100, tercile=None) -> dict:
    """pairs: columns context, unit, beta_A, gstar, breadth (A). Returns per-pair table, mean Delta and per-tercile
    means with 95% bootstrap CIs. tercile: per-pair labels fixed elsewhere (multiverse: the primary pairs' breadth
    terciles, review N8); default = terciles of this pair set's breadth (ties broken by pair order, N11)."""
    ys = {(r.context, r.unit): (score_L(E_L, r.context, r.unit), E_L[r.gstar].astype(float).clip(-3, 3))
          for r in pairs.itertuples(index=False)}
    cov_of = (lambda c: covariates_L.get(c)) if isinstance(covariates_L, dict) else (lambda c: covariates_L)
    if tercile is not None:
        terc = pd.Series(list(tercile), index=pairs.index)
        if len(terc) != len(pairs):
            raise ValueError("tercile labels do not match the pairs")
    else:
        terc = pd.qcut(pairs.breadth.rank(method="first"), 3, labels=["low", "mid", "high"]) if len(pairs) >= 3 \
            else pd.Series(["high"] * len(pairs), index=pairs.index)

    def deltas(rows_idx):
        out = []
        for r in pairs.itertuples(index=False):
            yc, yg = ys[(r.context, r.unit)]
            s = np.sign(r.beta_A)
            ix = E_L.index[rows_idx] if rows_idx is not None else E_L.index
            yc2, yg2 = yc.loc[ix], yg.loc[ix]
            yc2.index = yg2.index = pd.RangeIndex(len(ix))
            c2 = contexts_L[r.context].loc[ix].reset_index(drop=True)
            cv = cov_of(r.context)
            cv2 = None if cv is None else cv.loc[ix].reset_index(drop=True)
            out.append(_signed_t(yc2, c2, cv2, s) - _signed_t(yg2, c2, cv2, s))
        return np.array(out)

    d0 = deltas(None)
    rng = np.random.default_rng(seed)
    bl = (blocks_L.reindex(E_L.index) if blocks_L is not None else pd.Series(E_L.index, index=E_L.index))
    bl = bl.where(bl.notna(), pd.Series(E_L.index, index=E_L.index)).astype(str)
    groups = {b: np.flatnonzero(bl.values == b) for b in pd.unique(bl.values)}
    names = list(groups)
    boots = []
    for _ in range(B):
        pick = rng.integers(len(names), size=len(names))
        rows = np.concatenate([groups[names[i]] for i in pick])
        boots.append(deltas(rows))
    boots = np.array(boots)
    tab = pairs.assign(delta=d0, tercile=np.asarray(terc))
    summ = dict(mean=float(np.nanmean(d0)), ci95=np.nanquantile(np.nanmean(boots, 1), [0.025, 0.975]).tolist())
    for t in ("low", "mid", "high"):
        m = np.asarray(terc) == t
        if m.any():
            summ[f"mean_{t}"] = float(np.nanmean(d0[m]))
            summ[f"ci95_{t}"] = np.nanquantile(np.nanmean(boots[:, m], 1), [0.025, 0.975]).tolist()
    summ["top_tercile_wins"] = bool("ci95_high" in summ and summ["ci95_high"][0] > 0)
    return dict(table=tab, summary=summ)


def size_bin(n: int):
    for lo, hi in SIZE_BINS:
        if lo <= n <= hi:
            return f"{lo}-{hi}"
    return None


def bridge_rate(hits: pd.DataFrame, B: int = 2000, seed: int = 5200) -> dict:
    """hits: columns map ('zavg' | 'gls'), size (cluster size), confirmed (bool), lockbox-testable hits only."""
    h = hits.assign(bin=[size_bin(int(s)) for s in hits["size"]]).dropna(subset=["bin"])
    w = h.bin.value_counts(normalize=True)

    def rate(df):
        r = df.groupby("bin").confirmed.mean()
        ww = w.reindex(r.index)
        return float((r * ww).sum() / ww.sum()) if len(r) else np.nan

    obs = {m: rate(h[h["map"] == m]) for m in ("zavg", "gls")}
    rng = np.random.default_rng(seed)
    diffs = []
    cells = {k: g for k, g in h.groupby(["map", "bin"])}
    for _ in range(B):
        res = []
        for (m, b), g in cells.items():
            res.append(g.iloc[rng.integers(len(g), size=len(g))])
        bh = pd.concat(res)
        diffs.append(rate(bh[bh["map"] == "zavg"]) - rate(bh[bh["map"] == "gls"]))
    ci = np.nanquantile(diffs, [0.025, 0.975]).tolist()
    return dict(rate_zavg=obs["zavg"], rate_gls=obs["gls"], diff=obs["zavg"] - obs["gls"], ci95=ci,
                supports_bridge=bool(ci[0] > 0), bin_weights=w.to_dict())


# ------------------------------------------------------------------------------------------------ BENCH secondary
def effect_rsd(y: pd.Series, ctx: pd.Series, cov: pd.DataFrame | None) -> float:
    """Context effect in residual-SD units (BENCH secondary): OLS coefficient of the context column divided by the
    residual SD sqrt(RSS / (n - p)) of the same design (complete cases, assoc._design: same covariates and FE rules as
    the HC3 fits). NaN when < 2 context or < 2 reference models remain (robust.fit_one's rule)."""
    from part2.assoc import _design
    X, _, ok = _design(ctx.reindex(y.index), None if cov is None else cov.reindex(y.index))
    yy = y.values[ok].astype(float)
    m = ~np.isnan(yy)
    X, yy = X[m], yy[m]
    n_ctx = int(X[:, 1].sum()) if len(yy) else 0
    if n_ctx < 2 or len(yy) - n_ctx < 2 or len(yy) <= X.shape[1]:
        return np.nan
    b, *_ = np.linalg.lstsq(X, yy, rcond=None)
    r = yy - X @ b
    s = np.sqrt(float(r @ r) / (len(yy) - X.shape[1]))
    return float(b[1] / s) if s > 0 else np.nan


def _corr(a, b, method: str) -> float:
    a, b = np.asarray(a, float), np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3:
        return np.nan
    if method == "pearson":
        return float(np.corrcoef(a[ok], b[ok])[0, 1])
    if method == "spearman":
        from scipy.stats import spearmanr
        return float(spearmanr(a[ok], b[ok]).statistic)
    raise ValueError(f"correlation {method!r}: registry ruling must be 'pearson' or 'spearman'")


def al_correlation(pairs: pd.DataFrame, score_A, score_L, ctx_A: pd.DataFrame, ctx_L: pd.DataFrame, cov_A: dict,
                   cov_L: dict, sets_A=None, sets_L=None, R: int = 0, method: str | None = None) -> dict:
    """BENCH secondary: A-L correlation of context effects (residual-SD units) across (C, c) pairs for the complex
    scores, and for each matched random-set replicate r (sets_A(c, u, r) / sets_L(c, u, r) -> score Series or None
    when unit u has no r-th set); the correlation method is a registry ruling (no default). Returns the per-pair
    effects, r_complex, the R replicate correlations and the count of replicates with r_random >= r_complex (the
    summary statistic built from them is a registry ruling)."""
    if method is None:
        raise ValueError("BENCH secondary: the correlation method is not registered (ruling needed)")
    rows = []
    for p in pairs.itertuples(index=False):
        c, u = p.context, p.unit
        rows.append(dict(context=c, unit=u, eff_A=effect_rsd(score_A(c, u), ctx_A[c], cov_A.get(c)),
                         eff_L=effect_rsd(score_L(c, u), ctx_L[c], cov_L.get(c))))
    tab = pd.DataFrame(rows, columns=["context", "unit", "eff_A", "eff_L"])
    r_c = _corr(tab.eff_A, tab.eff_L, method)
    r_rand, n_used = [], []
    for r in range(R):
        a, b = [], []
        for p in pairs.itertuples(index=False):
            ya = sets_A(p.context, p.unit, r)
            if ya is None:
                continue
            a.append(effect_rsd(ya, ctx_A[p.context], cov_A.get(p.context)))
            b.append(effect_rsd(sets_L(p.context, p.unit, r), ctx_L[p.context], cov_L.get(p.context)))
        r_rand.append(_corr(a, b, method))
        n_used.append(len(a))
    rr = np.array(r_rand, float)
    return dict(table=tab, r_complex=r_c, n_pairs=int(tab[["eff_A", "eff_L"]].notna().all(axis=1).sum()),
                r_random=rr.tolist(), n_pairs_random=n_used, R=R, method=method,
                n_random_ge_complex=int(np.sum(rr[np.isfinite(rr)] >= r_c)) if np.isfinite(r_c) else None,
                random_quantiles=(np.nanquantile(rr, [0.025, 0.5, 0.975]).tolist() if np.isfinite(rr).any()
                                  else None))


def distinct_gstar(pairs: pd.DataFrame, fmap) -> dict:
    """BENCH secondary count: distinct g* genes, and distinct frozen coarse complexes containing at least one g*."""
    g = sorted(set(pairs.gstar.dropna()))
    cx = sorted({u for x in g for u in fmap.units_of_gene(x, 1)})
    return dict(n_pairs=int(len(pairs)), n_distinct_gstar=len(g), n_distinct_complexes_of_gstar=len(cx),
                n_distinct_complexes_in_pairs=int(pairs.unit.nunique()))


def _map_rate(df: pd.DataFrame, m: str) -> float:
    s = df[df["map"] == m]
    return float(s.confirmed.mean()) if len(s) else np.nan


def bridge_unweighted(hits: pd.DataFrame, B: int = 2000, seed: int = 5200) -> dict:
    """BRIDGE-L multiverse 'unweighted confirmation rate': confirmed / discovered per map (no size-bin weights), same
    hit set as bridge_rate (lockbox-testable hits in a size bin); difference zavg - GLS, with the bridge_rate
    resampling (within map and bin, seed) when B > 0."""
    h = hits.assign(bin=[size_bin(int(s)) for s in hits["size"]]).dropna(subset=["bin"])
    out = dict(rate_zavg=_map_rate(h, "zavg"), rate_gls=_map_rate(h, "gls"))
    out["diff"] = out["rate_zavg"] - out["rate_gls"]
    if B > 0:
        rng = np.random.default_rng(seed)
        cells = [g for _, g in h.groupby(["map", "bin"])]
        diffs = []
        for _ in range(B):
            bh = pd.concat([g.iloc[rng.integers(len(g), size=len(g))] for g in cells])
            diffs.append(_map_rate(bh, "zavg") - _map_rate(bh, "gls"))
        out["ci95"] = np.nanquantile(diffs, [0.025, 0.975]).tolist()
    return out
