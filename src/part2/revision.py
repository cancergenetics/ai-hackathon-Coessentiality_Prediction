"""Revision analyses R-1..R-4 (round-1 referee M1-M3). POST HOC and DESCRIPTIVE: nothing here is registered, nothing
writes a ledger or trials line, and no function touches lockbox outcome data. Intent:
results/part2/revision/INTENT.md (committed before computation).

* pool_matched_sets     R-1(a): matched random sets drawn from a gene pool outside every frozen unit, with the
                        acceptance rule of calib.matched_random_sets (size; essentiality +/- 0.1; coherence +/- 1.0 z;
                        R sets; 20,000 attempts; tolerances doubled once) and one of three candidate generators.
* discover              R-1(b): the discovery selection of assoc.run_association (BB -> permutation re-test of the
                        q < 0.1 screen -> BB again) on a table of tests, without ledger or outputs.
* symmetric_flag, conditional_auc   R-1(b) summaries conditional on discovery.
* ns_basis, spline_trend            R-2: natural cubic spline (df = 4) trend of beta on unit essentiality per context.
* fit_one_hc3                       R-3: one HC3 / Satterthwaite fit with the engine's design rules.
* tweedie_posterior, power_eb, standardised_rate   R-4: the lockbox results review's power aid on A estimates only.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats as sst

from part2 import stats

GENERATORS = ("nbr", "unif", "ess")


# ------------------------------------------------------------------------------------------------ R-1(a)
def pool_matched_sets(pool: list, gene_ess: pd.Series, pair_scores: pd.DataFrame, unit_members: dict, units: list,
                      R: int, seed: int, generator: str, tol_ess: float = 0.1, tol_coh: float = 1.0,
                      max_attempts: int = 20_000, nbr_min: int = 10) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    """{f'{unit}#r{k}': pool gene keys}, a per-unit status table and a per-set quality table.

    unit_members: {unit: scored member keys} (as calib.matched_random_sets: members that are A-effect columns).
    Target coherence uses the full pair_scores (S x S); candidate coherence uses the pool block. Generators:
      nbr   seed gene uniform from the pool + (k-1) drawn from its top-max(2k, nbr_min) pool neighbours by pair score
      unif  uniform k-subset of the pool, same acceptance
      ess   uniform k-subset, coherence not matched (size + essentiality only)
    Stream: numpy default_rng([seed, unit index, {nbr: 1, unif: 2, ess: 3}])."""
    if generator not in GENERATORS:
        raise ValueError(generator)
    gidx = GENERATORS.index(generator) + 1
    pool = list(pool)
    N = len(pool)
    ess_p = gene_ess.reindex(pool).values.astype(float)
    PSp = pair_scores.reindex(index=pool, columns=pool).values.astype(float)
    order = None
    if generator == "nbr":
        M = PSp.copy()
        np.fill_diagonal(M, -np.inf)
        M[~np.isfinite(M)] = -np.inf
        order = np.argsort(-M, axis=1, kind="stable")
    out, rows, qual = {}, [], []
    for ui, u in enumerate(units):
        mem = list(unit_members[u])
        k = len(mem)
        t_ess = float(np.nanmean(gene_ess.reindex(mem).values.astype(float))) if k else np.nan
        T = pair_scores.reindex(index=mem, columns=mem).values.astype(float)
        iu = np.triu_indices(k, 1)
        t_coh = float(np.nanmean(T[iu])) if k >= 2 else np.nan
        iu_k = np.triu_indices(k, 1)
        rng = np.random.default_rng([seed, ui, gidx])
        status, sets = "unavailable", []
        if 2 <= k <= N:
            for mult in (1.0, 2.0):
                seen, sets = set(), []
                for _ in range(max_attempts):
                    if generator == "nbr":
                        s0 = int(rng.integers(N))
                        m_nb = max(2 * k, nbr_min)
                        nb = order[s0][: min(m_nb, N - 1)]
                        if len(nb) < k - 1:
                            continue
                        s = np.sort(np.concatenate([[s0], rng.choice(nb, k - 1, replace=False)]))
                    else:
                        s = np.sort(rng.choice(N, k, replace=False))
                    key = s.tobytes()
                    if key in seen:
                        continue
                    if abs(float(np.nanmean(ess_p[s])) - t_ess) > tol_ess * mult:
                        continue
                    if generator != "ess" and np.isfinite(t_coh) and \
                            not abs(float(np.nanmean(PSp[np.ix_(s, s)][iu_k])) - t_coh) <= tol_coh * mult:
                        continue
                    seen.add(key)
                    sets.append(s)
                    if len(sets) >= R:
                        break
                if len(sets) >= R:
                    status = "ok" if mult == 1.0 else "tolerances doubled"
                    break
            if status == "unavailable":
                sets = []
        for r, s in enumerate(sets):
            name = f"{u}#n{r}"
            out[name] = [pool[i] for i in s]
            qual.append(dict(unit=u, set=name, size=k, d_ess=float(np.nanmean(ess_p[s])) - t_ess,
                             d_coh=float(np.nanmean(PSp[np.ix_(s, s)][iu_k])) - t_coh if k >= 2 else np.nan))
        rows.append(dict(unit=u, size=k, target_ess=t_ess, target_coh=t_coh, status=status, n_sets=len(sets)))
    return out, pd.DataFrame(rows), pd.DataFrame(qual, columns=["unit", "set", "size", "d_ess", "d_coh"])


def frac_lt(p, a: float = 0.05) -> dict:
    p = np.asarray(p, float)
    p = p[np.isfinite(p)]
    return dict(n_tests=int(len(p)), frac_p_lt_005=float(np.mean(p < a)) if len(p) else None)


# ------------------------------------------------------------------------------------------------ R-1(b)
def discover(res: pd.DataFrame, q: float, perm_fn, screen_q: float = 0.1) -> pd.DataFrame:
    """Discovery selection of assoc.run_association on a test table (columns context, p, evaluable): BB over
    families = contexts at q, permutation re-test of every test with adjusted q < screen_q (perm_fn(index list) ->
    {index: p_perm}), BB again on p_final. Returns res with q_adj, reject, p_perm, p_final, q_adj_perm, discovery."""
    from part2.assoc import _fdr
    ff = dict(method="benjamini_bogomolov", q=q, families="contexts")
    res = res.copy()
    res["evaluable"] = res.evaluable.astype(bool) & np.isfinite(res.p.astype(float))
    res["q_adj"], res["reject"], _ = _fdr(res, "p", ff, None)
    scr = res.index[res.q_adj.fillna(1) < screen_q]
    res["p_perm"] = pd.Series(perm_fn(list(scr)), dtype=float).reindex(res.index)
    res["p_final"] = res.p.where(res.p_perm.isna(), res.p_perm)
    res["q_adj_perm"], res["discovery"], _ = _fdr(res, "p_final", ff, None)
    res["discovery"] = res.discovery.astype(bool)
    return res


def symmetric_flag(real: pd.DataFrame, rand: pd.DataFrame) -> pd.DataFrame:
    """For every random test (context, unit, rep, abs_t): does its |t| exceed every other member of its group, i.e.
    the other matched sets of the unit in that context and the real unit? real: context, unit, abs_t."""
    rt = real.set_index(["context", "unit"]).abs_t
    out = rand.copy()
    mx_other = []
    for (c, u), g in rand.groupby(["context", "unit"], sort=False):
        a = g.abs_t.values.astype(float)
        r_t = float(rt.get((c, u), np.nan))
        for i in range(len(a)):
            others = np.concatenate([np.delete(a, i), [r_t]])
            mx_other.append((g.index[i], np.nanmax(others) if np.isfinite(others).any() else np.nan))
    m = pd.Series(dict(mx_other))
    out["max_abs_t_others"] = m.reindex(out.index)
    out["exceeds_all_others"] = out.abs_t > out.max_abs_t_others
    return out


def prob_greater(x: float, ys: np.ndarray) -> float:
    """P(x > Y) over the sample ys, ties counted 1/2; NaN for an empty sample."""
    ys = np.asarray(ys, float)
    ys = ys[np.isfinite(ys)]
    if not len(ys) or not np.isfinite(x):
        return np.nan
    return float((np.sum(x > ys) + 0.5 * np.sum(x == ys)) / len(ys))


def conditional_auc(real: pd.DataFrame, rand: pd.DataFrame, col: str, by=("context",)) -> pd.Series:
    """Per real row: P(real[col] > random[col]) against the random rows sharing the `by` keys."""
    by = list(by)
    grp = {(k if isinstance(k, tuple) else (k,)): g[col].values for k, g in rand.groupby(by)}
    return pd.Series([prob_greater(v, grp.get(tuple(k), np.array([])))
                      for v, k in zip(real[col].values, real[by].itertuples(index=False, name=None))],
                     index=real.index)


# ------------------------------------------------------------------------------------------------ R-2
def ns_knots(x: np.ndarray, df: int = 4) -> np.ndarray:
    """ns(x, df) knots: boundary = range, df - 1 interior knots at equally spaced quantiles."""
    x = np.asarray(x, float)
    inner = np.quantile(x, np.linspace(0, 1, df + 1)[1:-1])
    return np.concatenate([[x.min()], inner, [x.max()]])


def ns_basis(x: np.ndarray, knots: np.ndarray) -> np.ndarray:
    """Natural cubic spline basis with an intercept (ESL eq. 5.4-5.5): columns 1, x, d_k - d_{K-1} (k = 1..K-2);
    K knots -> K columns. Linear beyond the boundary knots."""
    x = np.asarray(x, float)
    xi = np.asarray(knots, float)
    K = len(xi)

    def d(k):
        return (np.clip(x - xi[k], 0, None) ** 3 - np.clip(x - xi[K - 1], 0, None) ** 3) / (xi[K - 1] - xi[k])
    cols = [np.ones_like(x), x] + [d(k) - d(K - 2) for k in range(K - 2)]
    return np.column_stack(cols)


def spline_trend(tests: pd.DataFrame, ess: pd.Series, df: int = 4) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per context: natural-spline (df) OLS of beta on unit essentiality over all evaluable tested units, plus the
    linear fit for reference. Returns (per-context table, per-test table with fitted_nl, residual_z_nl,
    concordant_nl, fitted_lin, residual_z_lin, concordant_lin)."""
    ev = tests.evaluable.astype(str).str.lower().eq("true") if tests.evaluable.dtype != bool else tests.evaluable
    t = tests[ev][["context", "unit", "beta"]].copy()
    t["beta"] = t.beta.astype(float)
    t["essentiality"] = t.unit.map(ess).astype(float)
    t = t[np.isfinite(t.beta) & np.isfinite(t.essentiality)]
    ctx_rows, out = [], []
    for c, g in t.groupby("context", sort=True):
        x, y = g.essentiality.values, g.beta.values
        kn = ns_knots(x, df)
        B = ns_basis(x, kn)
        coef, *_ = np.linalg.lstsq(B, y, rcond=None)
        fit = B @ coef
        res = y - fit
        sd = float(np.sqrt(np.sum(res ** 2) / (len(y) - B.shape[1])))
        b1, a1 = np.polyfit(x, y, 1)
        fit1 = a1 + b1 * x
        sd1 = float(np.std(y - fit1, ddof=2))
        rss0 = float(np.sum((y - y.mean()) ** 2))
        ctx_rows.append(dict(context=c, n_units=len(g), knots=";".join(f"{k:.4f}" for k in kn),
                             resid_sd_nl=sd, resid_sd_lin=sd1, r2_nl=1 - float(np.sum(res ** 2)) / rss0,
                             r2_lin=1 - float(np.sum((y - fit1) ** 2)) / rss0, slope_lin=float(b1),
                             spearman_rho=float(sst.spearmanr(x, y)[0]),
                             fitted_at_q10=float((ns_basis(np.quantile(x, [0.1]), kn) @ coef)[0]),
                             fitted_at_q50=float((ns_basis(np.quantile(x, [0.5]), kn) @ coef)[0]),
                             fitted_at_q90=float((ns_basis(np.quantile(x, [0.9]), kn) @ coef)[0]),
                             fitted_at_min=float((ns_basis(np.array([x.min()]), kn) @ coef)[0])))
        out.append(g.assign(fitted_nl=fit, residual_beta_nl=res, residual_z_nl=res / sd,
                            concordant_nl=np.sign(y) == np.sign(fit), fitted_lin=fit1,
                            residual_z_lin=(y - fit1) / sd1, concordant_lin=np.sign(y) == np.sign(fit1)))
    return pd.DataFrame(ctx_rows), pd.concat(out, ignore_index=True)


def scale_labels(rz_nl, concordant_nl, tier, control_recovery, direction) -> pd.DataFrame:
    """Uniform R-2 rule: scale_artefact_possible = concordant_nl & |rz_nl| < 2 (any context); corrected lead = H3a,
    T1, not a control readout, |rz_nl| >= 2, not scale_artefact_possible."""
    rz = np.asarray(rz_nl, float)
    conc = np.asarray(concordant_nl, bool)
    art = conc & (np.abs(rz) < 2)
    ctl = pd.Series(control_recovery).fillna("").astype(str).str.startswith("readout").values
    lead = (np.asarray(direction) == "H3a") & (np.asarray(tier) == "T1") & ~ctl & (np.abs(rz) >= 2) & ~art
    return pd.DataFrame(dict(scale_artefact_possible=art, lead_corrected=lead))


# ------------------------------------------------------------------------------------------------ R-3
def fit_one_hc3(y: pd.Series, ctx: pd.Series, cov: pd.DataFrame | None, min_n: int = 3) -> dict:
    """HC3 / Satterthwaite fit of the context column with the engine's design rules (assoc._design) on the models of y;
    NaN unless >= min_n context and reference models remain."""
    from part2.assoc import _design
    idx = y.index
    X, names, ok = _design(ctx.reindex(idx), None if cov is None else cov.reindex(idx))
    yy = y.values[ok].astype(float)
    m = np.isfinite(yy)
    X, yy = X[m], yy[m]
    n_ctx = int(X[:, 1].sum()) if len(X) else 0
    n = int(len(yy))
    out = dict(n=n, n_context=n_ctx, beta=np.nan, se=np.nan, t=np.nan, df=np.nan, n_covariates=len(names) - 2)
    if n_ctx < min_n or n - n_ctx < min_n:
        return out
    # covariates that became constant or collinear on the rows with a score are dropped by the same greedy rule
    keep =[0, 1] + [j for j in range(2, X.shape[1]) if np.unique(X[:, j]).size > 1]
    X = X[:, keep]
    sel = [0, 1]
    for j in range(2, X.shape[1]):
        if np.linalg.matrix_rank(X[:, sel + [j]]) == len(sel) + 1:
            sel.append(j)
    X = X[:, sel]
    try:
        d = stats.HC3Design(X, 1)
    except (np.linalg.LinAlgError, ValueError):
        return out
    b, se, t = d.fit(yy[:, None])
    out.update(beta=float(b[0]), se=float(se[0]), t=float(t[0]), df=float(d.df), n_covariates=X.shape[1] - 2)
    return out


# ------------------------------------------------------------------------------------------------ R-4
def t_to_z(t, df) -> np.ndarray:
    return sst.norm.ppf(sst.t.cdf(np.asarray(t, float), np.asarray(df, float))).clip(-12, 12)


def tweedie_posterior(z_all: np.ndarray, z_eval: np.ndarray, tweedie_fn) -> tuple[np.ndarray, np.ndarray]:
    """EB posterior mean and variance at z_eval from Tweedie's formula fitted on z_all (tweedie_fn = the lockbox
    results review's Lindsey-method `tweedie`), interpolated over the sorted z_all grid as in that review."""
    Em, Vm = tweedie_fn(z_all)
    o = np.argsort(z_all)
    return np.interp(z_eval, z_all[o], Em[o]), np.interp(z_eval, z_all[o], Vm[o])


def power_eb(z_A, shrink, post_var, se_A, se_L, df_L, absb, alpha: float = 0.05) -> dict:
    """Power of the one-sided alpha L test (the review's power aid): naive, EB posterior mean, EB mean + variance."""
    crit = sst.t.ppf(1 - alpha, df_L)
    k = se_A / se_L
    sh = np.clip(shrink, 0, None)
    m = np.abs(z_A) * sh
    return {"naive": sst.nct.sf(crit, df_L, absb / se_L),
            "eb_mean": sst.nct.sf(crit, df_L, absb * sh / se_L),
            "eb_mean_var": sst.norm.sf((crit - k * m) / np.sqrt(1 + k ** 2 * post_var))}


def standardised_rate(values: np.ndarray, bins: np.ndarray, weights: pd.Series) -> tuple[float, float]:
    """(sum_b w_b mean_b(values) / sum w, Poisson-binomial SD of that rate when values are success probabilities)."""
    v = pd.Series(np.asarray(values, float))
    b = pd.Series(np.asarray(bins, object))
    mean_b = v.groupby(b).mean()
    var_b = (v * (1 - v)).groupby(b).sum() / b.value_counts().reindex(mean_b.index) ** 2
    w = weights.reindex(mean_b.index)
    return float((mean_b * w).sum() / w.sum()), float(np.sqrt((var_b * w ** 2).sum()) / w.sum())


# ------------------------------------------------------------------------------------------------ permutation workers
_BLOCKS = None


def init_perm_worker(blocks: pd.Series, path0: str | None = None):
    """Process-pool initializer: the PatientID block labels (assoc.permutation_p `blocks`)."""
    global _BLOCKS
    _BLOCKS = blocks


def perm_task(a: tuple):
    """(key, y, X, ids, alternative, n_perm, seed sequence, strata) -> (key, Freedman-Lane permutation p), exactly as
    assoc.run_association re-tests a screened test (whole PatientID blocks within strata, GPD tail)."""
    from part2.assoc import permutation_p
    key, y, X, ids, alt, n_perm, seedseq, strata = a
    return key, permutation_p(y, X, pd.Index(ids), alt, n_perm, np.random.default_rng(list(seedseq)), _BLOCKS, strata)
