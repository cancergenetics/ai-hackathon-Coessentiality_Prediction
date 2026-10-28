"""Association statistics with known answers: HC3 vs statsmodels, Satterthwaite df vs the simulated variance of the
HC3 variance estimator, Freedman–Lane type-I under unequal variances (n1 = 8 vs n0 = 1000), blocked permutations,
GPD tail, Benjamini–Bogomolov FDR control on simulated families."""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from part2 import stats  # noqa: E402


def _toy(n=200, seed=0, het=True):
    rng = np.random.default_rng(seed)
    x = (rng.random(n) < 0.2).astype(float)
    z1, z2 = rng.normal(size=n), rng.integers(0, 3, n)
    X = np.column_stack([np.ones(n), x, z1, z2 == 1, z2 == 2]).astype(float)
    sd = np.where(x == 1, 3.0, 1.0) if het else 1.0
    y = 0.3 * x + 0.5 * z1 + rng.normal(size=n) * sd
    return y, X


def test_hc3_matches_statsmodels():
    import statsmodels.api as sm
    for seed in range(5):
        y, X = _toy(seed=seed)
        fit = sm.OLS(y, X).fit(cov_type="HC3")
        for j in range(X.shape[1]):
            r = stats.ols_hc3(y, X, j)
            assert np.isclose(r["beta"], fit.params[j], rtol=1e-10, atol=1e-12)
            assert np.isclose(r["se"], fit.bse[j], rtol=1e-10)
            assert np.isclose(r["t"], fit.tvalues[j], rtol=1e-10)
    # vectorised fit over several outcomes equals column-wise fits
    y, X = _toy(seed=9)
    Y = np.column_stack([y, 2 * y + 1, y[::-1]])
    d = stats.HC3Design(X, 1)
    b, se, t = d.fit(Y)
    for k in range(3):
        r = stats.ols_hc3(Y[:, k], X, 1)
        assert np.isclose(b[k], r["beta"]) and np.isclose(se[k], r["se"])


def test_satterthwaite_df_matches_simulated_moments():
    """Under homoskedastic normal errors var(V)/E(V)^2 = 2/df for the HC3 variance V (Satterthwaite moment match)."""
    rng = np.random.default_rng(1)
    n = 60
    x = np.r_[np.ones(6), np.zeros(n - 6)]
    X = np.column_stack([np.ones(n), x, rng.normal(size=n)])
    d = stats.HC3Design(X, 1)
    E = rng.normal(size=(n, 40000))
    V = d.fit(E)[1] ** 2
    df_sim = 2 * V.mean() ** 2 / V.var()
    assert 0 < d.df < n - X.shape[1]
    assert abs(df_sim / d.df - 1) < 0.05, (df_sim, d.df)
    # brute-force eigenvalue formula
    M = np.eye(n) - d.Q @ d.Q.T
    A = M @ np.diag(d.w ** 2) @ M
    lam = np.linalg.eigvalsh(A)
    assert np.isclose(d.df, lam.sum() ** 2 / (lam ** 2).sum(), rtol=1e-8)


@pytest.mark.parametrize("sd1", [3.0, 1 / 3])
def test_freedman_lane_type1_unequal_variance(sd1):
    """n1 = 8 context models vs n0 = 1000, error SD 3 (or 1/3) in the context group, covariate present.
    Studentised (HC3) Freedman–Lane: rejection rate at alpha = 0.05 should be ~0.05."""
    rng = np.random.default_rng(42)
    n1, n0, reps, B = 8, 1000, 1000, 499
    n = n1 + n0
    x = np.r_[np.ones(n1), np.zeros(n0)]
    rej_perm = rej_t = 0
    for r in range(reps):
        z = rng.normal(size=n)
        X = np.column_stack([np.ones(n), x, z])
        y = 0.5 * z + rng.normal(size=n) * np.where(x == 1, sd1, 1.0)
        perms = np.array([rng.permutation(n) for _ in range(B)])
        out = stats.freedman_lane(y, X, 1, perms, "two-sided", gpd=False)
        rej_perm += out["p_perm"] <= 0.05
        rej_t += stats.ols_hc3(y, X, 1)["p"] <= 0.05
    rate, rate_t = rej_perm / reps, rej_t / reps
    print(f"sd1={sd1}: FL type-I {rate:.3f}; HC3+Satterthwaite type-I {rate_t:.3f}")
    assert 0.03 <= rate <= 0.07, rate
    assert rate_t <= 0.10, rate_t


def test_blocked_permutations_respect_blocks_and_strata():
    rng = np.random.default_rng(0)
    blocks = np.array(["p1", "p1", "p2", "p3", "p3", "p4", "p5", "p6", "p6", "p7"])
    strata = np.array(["a", "a", "a", "a", "a", "a", "b", "b", "b", "b"])
    P = stats.blocked_permutations(200, rng, blocks, strata)
    for pi in P:
        assert sorted(pi) == list(range(10))
        assert all(strata[pi] == strata)                       # lineage preserved
        # the two members of a 2-model block map to the two members of one 2-model block (same stratum)
        for b in ("p1", "p3", "p6"):
            idx = np.flatnonzero(blocks == b)
            src = set(blocks[pi[idx]])
            assert len(src) == 1 and (blocks == src.pop()).sum() == 2
    moved = (P != np.arange(10)).any(0)
    assert moved[[0, 1, 3, 4, 2, 5]].all()                    # stratum a blocks are exchanged
    assert not moved[[7, 8]].any()                            # p6: the only 2-model block in stratum b stays
    assert set(np.unique(P[:, 6])) == {6, 9}                  # p5 / p7 (size 1, stratum b) swap only with each other


def test_gpd_tail():
    rng = np.random.default_rng(3)
    null = rng.normal(size=20000)
    from scipy.stats import norm
    for t in (3.8, 4.2):
        p = stats.perm_pvalue(t, null, "greater", gpd=True)
        true = norm.sf(t)
        assert true / 5 < p < true * 5, (t, p, true)
    assert stats.perm_pvalue(1.0, null, "greater") == pytest.approx((np.sum(null >= 1.0) + 1) / 20001)


def test_gpd_only_below_registered_resolution():
    """Registered rule (H3a/H3b/H2 `statistic`): GPD tail only when p_perm < 1e-4. Above that the exact permutation p
    (k + 1) / (B + 1) is returned, even with < 10 exceedances: the GPD extrapolation's noise made the re-entered p
    anti-conservative at BH thresholds (size 1.4-1.7x at 5e-4..2e-3 with B = 1000; results/part2/verify/
    h3a_fdr_diagnosis.md)."""
    rng = np.random.default_rng(5)
    null = rng.normal(size=1000)
    for t in (2.8, 3.1, 6.0):                                   # 0-5 exceedances: still the exact p at B = 1000
        assert stats.perm_pvalue(t, null, "greater") == (np.sum(null >= t) + 1) / 1001
    null = rng.normal(size=10000)
    t = 3.0                                                   # p ~ 1.3e-3: exact permutation p at B = 10,000
    assert stats.perm_pvalue(t, null, "greater") == (np.sum(null >= t) + 1) / 10001
    assert stats.perm_pvalue(5.0, null, "greater") < 1e-4      # beyond every permutation: GPD tail


def test_bh_known_answer():
    p = np.array([0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205, 0.212, 0.216])
    rej, adj = stats.bh(p, 0.05)
    assert rej.tolist() == [True, True] + [False] * 8
    from statsmodels.stats.multitest import multipletests
    assert np.allclose(adj, multipletests(p, method="fdr_bh")[1])


def test_benjamini_bogomolov_controls_average_fdr():
    """M = 30 families x 40 tests; 6 families carry 10 signals each. BB's error measure (expected average FDP over
    selected families) must be <= q; global null: P(any family selected) <= q."""
    rng = np.random.default_rng(7)
    q, M, m, reps = 0.1, 30, 40, 400
    fam = np.repeat(np.arange(M), m)
    from scipy.stats import norm
    avg_fdp, sel_null = [], 0
    for r in range(reps):
        mu = np.zeros(M * m)
        for f in range(6):
            mu[f * m: f * m + 10] = 3.0
        z = rng.normal(size=M * m) + mu
        p = norm.sf(z)
        bb = stats.benjamini_bogomolov(p, fam, q)
        null = mu == 0
        fdps = []
        for f in bb["families"].family[bb["families"].selected]:
            i = fam.astype(str) == f
            Rf = bb["reject"][i].sum()
            fdps.append((bb["reject"][i] & null[i]).sum() / max(Rf, 1))
        avg_fdp.append(np.sum(fdps) / max(bb["R"], 1))
        p0 = norm.sf(rng.normal(size=M * m))
        sel_null += stats.benjamini_bogomolov(p0, fam, q)["R"] > 0
    se = np.std(avg_fdp) / np.sqrt(reps)
    print(f"BB average FDP {np.mean(avg_fdp):.4f} (q={q}); global-null selection rate {sel_null / reps:.3f}")
    assert np.mean(avg_fdp) <= q + 2 * se
    assert sel_null / reps <= q + 2 * np.sqrt(q * (1 - q) / reps)
