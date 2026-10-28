"""Bootstrap resampling schemes (docs/TEST_PROTOCOL.md §6). Every scheme yields per-pair multiplicity weights, so
every statistic is computed on the same replicates for all methods (paired), on top of design weights ``base_w``.

* ``cluster_weights`` (PRIMARY): Part 1's paired cluster bootstrap (src/metrics.cluster_bootstrap) re-expressed as a
  weight generator with the identical random stream: clusters = np.unique(r_[ca, cb]); per replicate
  m = bincount(rng.integers(0, C, C)); within-cluster pair weight m_c, between-cluster pair weight m_c1 * m_c2.
  Equality of AP replicates with Part 1's function is tested.
* ``protein_weights`` (secondary): the same with proteins as units (pair weight m_a * m_b).
* ``pair_weights`` (secondary, anti-conservative under dependence): Part 1's stratified pair bootstrap
  (src/metrics._boot_weights: positives and negatives resampled separately).
Intervals: percentile; p_two_sided = 2 min(P(d <= 0), P(d >= 0)) (Part 1 rule).
"""
from __future__ import annotations

import numpy as np

from testeval.part1 import metrics as p1m


def cluster_weights(ca, cb, n_boot: int, seed: int):
    ca, cb = np.asarray(ca), np.asarray(cb)
    uniq, inv = np.unique(np.r_[ca, cb], return_inverse=True)
    ia, ib = inv[:len(ca)], inv[len(ca):]
    same = ia == ib
    C = len(uniq)
    rng = np.random.default_rng(seed)
    for _ in range(n_boot):
        m = np.bincount(rng.integers(0, C, C), minlength=C).astype(np.float64)
        yield np.where(same, m[ia], m[ia] * m[ib])


def protein_weights(pa, pb, n_boot: int, seed: int):
    pa, pb = np.asarray(pa), np.asarray(pb)
    uniq, inv = np.unique(np.r_[pa, pb], return_inverse=True)
    ia, ib = inv[:len(pa)], inv[len(pa):]
    C = len(uniq)
    rng = np.random.default_rng(seed)
    for _ in range(n_boot):
        m = np.bincount(rng.integers(0, C, C), minlength=C).astype(np.float64)
        yield m[ia] * m[ib]


def pair_weights(y, n_boot: int, seed: int):
    yield from p1m._boot_weights(y, n_boot, seed)


def _flat(fns: dict, w) -> dict:
    out = {}
    for name, f in fns.items():
        for k, v in f(w).items():
            out[f"{name}|{k}"] = v
    return out


def run(weight_gen, base_w, fns: dict, diffs=(), alpha: float = 0.05) -> dict:
    """fns: name -> fn(w) -> {statistic: value} (one call per method computes all its statistics from one sort);
    keys are flattened to 'name|statistic'. diffs: pairs of flattened keys (paired differences, same replicates).
    Returns point (w = base_w), percentile CI, bootstrap SD, and paired differences with relative gain."""
    base = np.asarray(base_w, dtype=np.float64)
    point = _flat(fns, base)
    stats = list(point)
    reps = {k: [] for k in stats}
    n = 0
    for m in weight_gen:
        r = _flat(fns, m * base)
        for k in stats:
            reps[k].append(r[k])
        n += 1
    reps = {k: np.asarray(v, dtype=np.float64) for k, v in reps.items()}
    q = [alpha / 2, 1 - alpha / 2]
    out = dict(n_boot=n, alpha=alpha, stat={}, diff={})
    for k in stats:
        r = reps[k]
        lo, hi = np.nanquantile(r, q) if np.isfinite(r).any() else (np.nan, np.nan)
        out["stat"][k] = dict(point=point[k], ci_low=float(lo), ci_high=float(hi),
                              boot_sd=float(np.nanstd(r, ddof=1)), n_nan=int(np.isnan(r).sum()))
    for a, b in diffs:
        d = reps[a] - reps[b]
        with np.errstate(divide="ignore", invalid="ignore"):
            rel = reps[a] / reps[b] - 1
        lo, hi = np.nanquantile(d, q)
        rlo, rhi = np.nanquantile(rel, q)
        p = 2 * min(np.nanmean(d <= 0), np.nanmean(d >= 0))
        out["diff"][f"{a} - {b}"] = dict(
            diff=point[a] - point[b], ci_low=float(lo), ci_high=float(hi),
            rel_gain=(point[a] / point[b] - 1) if point[b] else float("nan"), rel_ci_low=float(rlo),
            rel_ci_high=float(rhi), p_two_sided=float(min(p, 1.0)))
    return out


def holm(pvals: dict) -> dict:
    return p1m.holm(pvals)
