"""Association statistics (DESIGN v3.1 §2.5): HC3 OLS with Satterthwaite df, Freedman–Lane permutation restricted to
PatientID blocks within lineage, GPD tail p-values, Benjamini–Bogomolov hierarchical FDR.

HC3:  V_j = sum_i g_ij^2 e_i^2 / (1 - h_i)^2 with G = X (X'X)^-1 (so beta_j = g_j' y), h = diag(hat matrix).
Satterthwaite df (Bell & McCaffrey 2002 style, normal homoskedastic working model): V_j = e' W^2 e with
W = diag(g_ij / (1 - h_i)) and e = M y, so V_j ~ sigma^2 sum_k lambda_k chi2_1 with lambda = eig(M W^2 M) and
df = (tr A)^2 / tr(A^2), A = M W^2 M. With Q the thin-QR basis of X (M = I - QQ') and a_i = w_i^2:
  tr A = sum_i a_i (1 - h_i),  tr A^2 = sum_i a_i^2 (1 - 2 h_i) + ||Q' diag(a) Q||_F^2.
The df depends on X and the tested column only.

Freedman–Lane: y = Z gamma + beta x + e. Reduced-model residuals r = M_Z y are permuted (y* = Z gamma_hat + P r) and
the full model is refitted; the statistic is the HC3 t of x. Because g_j' Z = 0, beta*_j = g_j' P r and
e* = M P r, which vectorises over permutations. Permutations exchange whole PatientID blocks of equal size within a
lineage stratum (within-block order kept; Winkler et al. 2015 'whole-block' exchangeability), so models of the same
patient never get split and lineage composition is preserved.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats as sst


# ------------------------------------------------------------------------------------------------ HC3 OLS
class HC3Design:
    """Pre-computed pieces of a full-rank design X (n x p) for HC3 inference on one column j."""

    def __init__(self, X: np.ndarray, j: int):
        X = np.asarray(X, float)
        n, p = X.shape
        Q, R = np.linalg.qr(X)
        if np.min(np.abs(np.diag(R))) < 1e-10 * max(1.0, np.max(np.abs(np.diag(R)))):
            raise np.linalg.LinAlgError("design matrix is rank deficient")
        self.n, self.p, self.j = n, p, j
        self.Q, self.R = Q, R
        self.h = np.einsum("ij,ij->i", Q, Q)
        if np.any(self.h > 1 - 1e-10):
            raise ValueError("leverage 1 observation (a context or covariate level with a single model)")
        G = np.linalg.solve(R, Q.T).T                       # X (X'X)^-1 = Q R^-T
        self.g = G[:, j]
        self.w = self.g / (1 - self.h)
        a = self.w ** 2
        trA = np.sum(a * (1 - self.h))
        QaQ = Q.T @ (a[:, None] * Q)
        trA2 = np.sum(a ** 2 * (1 - 2 * self.h)) + np.sum(QaQ ** 2)
        self.df = float(trA ** 2 / trA2)

    def resid(self, Y: np.ndarray) -> np.ndarray:
        Y = np.asarray(Y, float)
        return Y - self.Q @ (self.Q.T @ Y)

    def fit(self, Y: np.ndarray):
        """Y: n or n x k. Returns beta_j, se_HC3, t (arrays of length k)."""
        Y = np.asarray(Y, float)
        one = Y.ndim == 1
        Y2 = Y[:, None] if one else Y
        beta = self.g @ Y2
        E = self.resid(Y2)
        var = (self.w[:, None] ** 2 * E ** 2).sum(0)
        se = np.sqrt(var)
        t = beta / se
        return (beta[0], se[0], t[0]) if one else (beta, se, t)


def pvalue_t(t, df, alternative: str = "two-sided"):
    t = np.asarray(t, float)
    if alternative == "two-sided":
        return 2 * sst.t.sf(np.abs(t), df)
    if alternative == "greater":
        return sst.t.sf(t, df)
    if alternative == "less":
        return sst.t.cdf(t, df)
    raise ValueError(alternative)


def ols_hc3(y, X, j: int, alternative: str = "two-sided") -> dict:
    d = HC3Design(X, j)
    b, se, t = d.fit(y)
    return dict(beta=float(b), se=float(se), t=float(t), df=d.df, p=float(pvalue_t(t, d.df, alternative)), n=d.n)


# ------------------------------------------------------------------------------------------------ permutations
def blocked_permutations(n_perm: int, rng: np.random.Generator, blocks=None, strata=None, n: int | None = None):
    """n_perm x n index arrays. Whole blocks (e.g. PatientID) of equal size are exchanged within a stratum (e.g.
    lineage); within-block order (sorted position) is kept. blocks/strata: length-n label arrays (None = singletons /
    one stratum). A block spanning strata is assigned to the stratum of its first model."""
    if n is None:
        n = len(blocks) if blocks is not None else len(strata)
    blocks = np.arange(n).astype(str) if blocks is None else np.asarray(blocks).astype(str)
    strata = np.zeros(n, int).astype(str) if strata is None else np.asarray(strata).astype(str)
    df = pd.DataFrame(dict(b=blocks, s=strata, i=np.arange(n)))
    members = {b: g.i.values for b, g in df.groupby("b", sort=True)}
    bst = df.groupby("b", sort=True).s.first()
    groups = []
    for (s, k), g in pd.DataFrame(dict(b=list(members), s=bst.reindex(list(members)).values,
                                       k=[len(v) for v in members.values()])).groupby(["s", "k"], sort=True):
        if len(g) > 1:
            groups.append(np.stack([members[b] for b in g.b]))       # n_blocks x k (positions)
    out = np.tile(np.arange(n), (n_perm, 1))
    for P in groups:
        nb = P.shape[0]
        for r in range(n_perm):
            pi = rng.permutation(nb)
            out[r, P.ravel()] = P[pi].ravel()
    return out


def freedman_lane(y, X, j: int, perms: np.ndarray, alternative: str = "two-sided", gpd: bool = True) -> dict:
    """Freedman–Lane permutation p-value of the HC3 t of column j (X includes the intercept). perms: B x n."""
    y = np.asarray(y, float)
    X = np.asarray(X, float)
    d = HC3Design(X, j)
    _, _, t_obs = d.fit(y)
    Z = np.delete(X, j, axis=1)
    Qz, _ = np.linalg.qr(Z)
    r = y - Qz @ (Qz.T @ y)
    Rp = r[perms].T                                          # n x B
    beta = d.g @ Rp
    E = d.resid(Rp)
    t_null = beta / np.sqrt((d.w[:, None] ** 2 * E ** 2).sum(0))
    p = perm_pvalue(t_obs, t_null, alternative, gpd=gpd)
    return dict(t=float(t_obs), p_perm=p, n_perm=int(perms.shape[0]), df=d.df)


def perm_pvalue(t_obs: float, t_null: np.ndarray, alternative: str = "two-sided", gpd: bool = True,
                n_exc: int = 250, min_exceed: int = 10, gpd_below: float = 1e-4) -> float:
    """(1 + #exceed) / (B + 1); only when that empirical p is below gpd_below (registered rule: 'GPD tail when
    p_perm < 1e-4', with 10,000 permutations) and there are fewer than min_exceed exceedances is a generalised Pareto
    tail fitted to the n_exc largest null statistics (Knijnenburg et al. 2009): p = (n_exc / B) * (1 - F_GPD(t - u)).
    The empirical p is exactly valid; the GPD estimate is not (its noise makes P(p_hat <= x) > x at small x), so it is
    confined to the range the permutations cannot resolve (results/part2/verify/h3a_fdr_diagnosis.md)."""
    s = {"two-sided": np.abs, "greater": lambda v: v, "less": lambda v: -v}[alternative]
    T, tn = float(s(np.asarray(t_obs))), s(np.asarray(t_null, float))
    B = len(tn)
    k = int(np.sum(tn >= T))
    p_emp = (k + 1) / (B + 1)
    if not gpd or k >= min_exceed or B < 4 * n_exc or p_emp >= gpd_below:
        return float(p_emp)
    srt = np.sort(tn)[::-1]
    u = 0.5 * (srt[n_exc - 1] + srt[n_exc])
    exc = srt[:n_exc] - u
    try:
        c, loc, scale = sst.genpareto.fit(exc, floc=0)
        p = n_exc / B * sst.genpareto.sf(T - u, c, loc=0, scale=scale)
    except Exception:  # noqa: BLE001
        return float(p_emp)
    if not np.isfinite(p) or p <= 0:
        # observed statistic beyond the support of a bounded (shape < 0) GPD fit: fall back to the exponential tail
        # (shape 0, scale = mean excess), which is heavier and never 0 (no arbitrary floor)
        p = n_exc / B * float(np.exp(-(T - u) / max(float(np.mean(exc)), 1e-12)))
    return float(min(p, p_emp))


# ------------------------------------------------------------------------------------------------ FDR
def bh(p, q: float):
    """Benjamini–Hochberg: (reject mask, adjusted p)."""
    p = np.asarray(p, float)
    m = len(p)
    if m == 0:
        return np.zeros(0, bool), np.zeros(0)
    o = np.argsort(p)
    adj = np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1]
    out = np.empty(m)
    out[o] = np.minimum(adj, 1.0)
    return out <= q, out


def simes(p) -> float:
    p = np.sort(np.asarray(p, float))
    m = len(p)
    return float(np.min(p * m / np.arange(1, m + 1))) if m else 1.0


def benjamini_bogomolov(p, family, q: float = 0.05) -> dict:
    """Families selected by BH on Simes family p-values at q (R of M selected); within each selected family BH at
    level q R / M. Returns per-test rejections, per-test adjusted p (within-family BH x M / R, NaN if family not
    selected), family table. Controls the expected average FDP over selected families (Benjamini & Bogomolov 2014)."""
    p = np.asarray(p, float)
    fam = np.asarray(family).astype(str)
    names = np.array(sorted(set(fam)))
    fp = np.array([simes(p[fam == f]) for f in names])
    sel, fadj = bh(fp, q)
    M, R = len(names), int(sel.sum())
    rej = np.zeros(len(p), bool)
    adj = np.full(len(p), np.nan)
    for f, s in zip(names, sel):
        if not s:
            continue
        idx = np.flatnonzero(fam == f)
        r_, a_ = bh(p[idx], q * R / M)
        rej[idx] = r_
        adj[idx] = np.minimum(a_ * M / R, 1.0)
    ftab = pd.DataFrame(dict(family=names, simes_p=fp, family_q=fadj, selected=sel,
                             n_tests=[int((fam == f).sum()) for f in names],
                             n_rejected=[int(rej[fam == f].sum()) for f in names]))
    return dict(reject=rej, adj=adj, families=ftab, R=R, M=M)
