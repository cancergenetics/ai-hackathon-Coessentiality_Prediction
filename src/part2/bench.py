"""Benchmark and bridge-L estimands (DESIGN v3.1 §2.3, §2.5; BENCH_complex_vs_gene, BRIDGE-L_confirmation_rate).
Both read the lockbox only after every confirmatory entry is closed (enforced by assoc.read_lockbox).

BENCH: Delta(C, c) = z_L(S_C) - z_L(g*), z_L = HC3 t of the context coefficient in L signed so that positive = the
A direction; g* = member of C (after removals) with the largest A t in the A direction (ties: larger |beta|), chosen
in A only. Breadth terciles fixed on the primary pairs' A breadth. CI by PatientID-blocked bootstrap over L models
(B = 2000, seed 5100; pairs and g* fixed, t recomputed per resample).

BRIDGE-L: per map, confirmation rate (confirmed / discovered, lockbox-testable contexts) within size bins
{3-4, 5-8, 9-16, 17-60}, averaged with the pooled (both maps) discovered-cluster size distribution; difference
zavg - GLS with a bootstrap over discovered hits within map and bin (B = 2000, seed 5200).
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
                blocks_L: pd.Series | None = None, B: int = 2000, seed: int = 5100) -> dict:
    """pairs: columns context, unit, beta_A, gstar, breadth (A). Returns per-pair table, mean Delta and per-tercile
    means with 95% bootstrap CIs."""
    ys = {(r.context, r.unit): (score_L(E_L, r.context, r.unit), E_L[r.gstar].astype(float).clip(-3, 3))
          for r in pairs.itertuples(index=False)}
    cov_of = (lambda c: covariates_L.get(c)) if isinstance(covariates_L, dict) else (lambda c: covariates_L)
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

    obs = {m: rate(h[h.map == m]) for m in ("zavg", "gls")}
    rng = np.random.default_rng(seed)
    diffs = []
    cells = {k: g for k, g in h.groupby(["map", "bin"])}
    for _ in range(B):
        res = []
        for (m, b), g in cells.items():
            res.append(g.iloc[rng.integers(len(g), size=len(g))])
        bh = pd.concat(res)
        diffs.append(rate(bh[bh.map == "zavg"]) - rate(bh[bh.map == "gls"]))
    ci = np.nanquantile(diffs, [0.025, 0.975]).tolist()
    return dict(rate_zavg=obs["zavg"], rate_gls=obs["gls"], diff=obs["zavg"] - obs["gls"], ci95=ci,
                supports_bridge=bool(ci[0] > 0), bin_weights=w.to_dict())
