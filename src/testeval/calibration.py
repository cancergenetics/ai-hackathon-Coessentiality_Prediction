"""Score -> probability calibration fitted on TRAIN pairs only and applied unchanged to test pairs
(docs/TEST_PROTOCOL.md §5.3).

* Isotonic (primary): sklearn IsotonicRegression(increasing, y in [0, 1], out_of_bounds='clip') with sample weights,
  stored as its (X_thresholds_, y_thresholds_) knots and applied with np.interp, which is exactly sklearn's predict
  (linear interpolation between knots, clipped outside; tested).
* Platt (secondary): weighted maximum-likelihood logistic regression p = sigmoid(a + b * z) on the standardised score
  z = (s - mu) / sd (mu, sd from the fit set), no penalty, no target smoothing (n_pos in the thousands).
* Unmapped pairs: their floor score is file-specific (min mapped - 1), so before applying a calibrator the score of
  every unmapped pair is replaced by the fit set's floor (unmapped = the lowest score level in both sets).
Metrics: Brier, Brier skill vs a constant forecast, calibration-in-the-large (sum p / sum y), calibration slope and
intercept (logistic recalibration of y on logit p), ECE with 20 equal-mass bins (primary) and 15 equal-width bins
on [0, 1] and 16 log-spaced bins on [1e-4, 1] (secondary), and reliability tables (20 equal-mass bins; 16
log-spaced bins, which resolve the high-probability tail at ~0.5% prevalence and are used for the diagram). All accept design weights.
"""
from __future__ import annotations

import numpy as np
from scipy import optimize
from scipy.special import expit

EPS = 1e-6                        # probability clip for logit (calibration slope / intercept)


def fit_isotonic(s, y, w=None) -> dict:
    from sklearn.isotonic import IsotonicRegression
    ir = IsotonicRegression(y_min=0.0, y_max=1.0, increasing=True, out_of_bounds="clip")
    ir.fit(np.asarray(s, np.float64), np.asarray(y, np.float64), sample_weight=w)
    return dict(kind="isotonic", x=ir.X_thresholds_.tolist(), y=ir.y_thresholds_.tolist())


def fit_platt(s, y, w=None) -> dict:
    s = np.asarray(s, np.float64)
    y = np.asarray(y, np.float64)
    w = np.ones(len(s)) if w is None else np.asarray(w, np.float64)
    mu, sd = float(np.average(s, weights=w)), float(np.sqrt(np.average((s - np.average(s, weights=w)) ** 2, weights=w)))
    z = (s - mu) / sd
    a, b = _logit_fit(z, y, w)
    return dict(kind="platt", mu=mu, sd=sd, a=a, b=b)


def _logit_fit(x, y, w):
    W = w.sum()

    def f(th):
        eta = th[0] + th[1] * x
        # weighted negative log-likelihood (stable): log(1 + e^eta) - y * eta
        nll = np.sum(w * (np.logaddexp(0.0, eta) - y * eta)) / W
        p = expit(eta)
        g = np.array([np.sum(w * (p - y)), np.sum(w * (p - y) * x)]) / W
        return nll, g

    p0 = np.clip(np.sum(w * y) / W, 1e-6, 1 - 1e-6)
    r = optimize.minimize(f, x0=np.array([np.log(p0 / (1 - p0)), 0.0]), jac=True, method="BFGS",
                          options=dict(gtol=1e-10, maxiter=1000))
    return float(r.x[0]), float(r.x[1])


def apply(cal: dict, s) -> np.ndarray:
    s = np.asarray(s, np.float64)
    if cal["kind"] == "isotonic":
        return np.interp(s, np.asarray(cal["x"]), np.asarray(cal["y"]))
    if cal["kind"] == "platt":
        return expit(cal["a"] + cal["b"] * (s - cal["mu"]) / cal["sd"])
    raise ValueError(cal["kind"])


def floor_to_fit(s, mapped, fit_floor: float) -> np.ndarray:
    s = np.asarray(s, np.float64).copy()
    s[~np.asarray(mapped, bool)] = fit_floor
    return s


def brier(p, y, w=None) -> float:
    return float(np.average((np.asarray(p) - np.asarray(y, np.float64)) ** 2, weights=w))


def _bins_equal_mass(p, w, n_bins):
    order = np.argsort(p, kind="mergesort")
    cw = np.cumsum(w[order])
    before = cw - w[order]
    b = np.minimum((before / cw[-1] * n_bins).astype(int), n_bins - 1)
    out = np.empty(len(p), dtype=int)
    out[order] = b
    return out


def reliability(p, y, w=None, n_bins: int = 20, scheme: str = "mass") -> list:
    p = np.asarray(p, np.float64)
    y = np.asarray(y, np.float64)
    w = np.ones(len(p)) if w is None else np.asarray(w, np.float64)
    if scheme == "mass":
        b = _bins_equal_mass(p, w, n_bins)
    elif scheme == "log":            # n_bins log-spaced bins on [1e-4, 1]; p < 1e-4 joins the first bin
        edges = np.logspace(-4, 0, n_bins + 1)
        b = np.clip(np.searchsorted(edges, p, side="right") - 1, 0, n_bins - 1)
    else:
        b = np.minimum((p * n_bins).astype(int), n_bins - 1)
    rows = []
    for k in range(n_bins):
        m = b == k
        if not m.any():
            continue
        ww = w[m].sum()
        rows.append(dict(bin=k, weight=float(ww), n=int(m.sum()), mean_pred=float(np.sum(w[m] * p[m]) / ww),
                         obs_rate=float(np.sum(w[m] * y[m]) / ww), p_min=float(p[m].min()), p_max=float(p[m].max())))
    return rows


def ece(p, y, w=None, n_bins: int = 20, scheme: str = "mass") -> float:
    rows = reliability(p, y, w, n_bins, scheme)
    tot = sum(r["weight"] for r in rows)
    return float(sum(r["weight"] * abs(r["mean_pred"] - r["obs_rate"]) for r in rows) / tot)


def summary(p, y, w=None, ref_prevalence: float | None = None) -> dict:
    p = np.asarray(p, np.float64)
    y = np.asarray(y, np.float64)
    w = np.ones(len(p)) if w is None else np.asarray(w, np.float64)
    obs = float(np.sum(w * y) / w.sum())
    b = brier(p, y, w)
    out = dict(brier=b, observed_prevalence=obs, mean_pred=float(np.sum(w * p) / w.sum()),
               calibration_in_the_large=float(np.sum(w * p) / np.sum(w * y)) if np.sum(w * y) > 0 else float("nan"),
               ece_mass20=ece(p, y, w, 20, "mass"), ece_width15=ece(p, y, w, 15, "width"),
               ece_log16=ece(p, y, w, 16, "log"))
    if ref_prevalence is not None:
        bref = brier(np.full(len(p), ref_prevalence), y, w)
        out.update(brier_ref_fit_prevalence=bref, bss_vs_fit_prevalence=1 - b / bref)
    bobs = obs * (1 - obs)
    out.update(brier_ref_observed_prevalence=bobs, bss_vs_observed_prevalence=1 - b / bobs if bobs > 0 else float("nan"))
    lp = np.log(np.clip(p, EPS, 1 - EPS) / np.clip(1 - p, EPS, 1))
    a, s = _logit_fit(lp, y, w)
    out.update(calibration_intercept=a, calibration_slope=s)
    return out
