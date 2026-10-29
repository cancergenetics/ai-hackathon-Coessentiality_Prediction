"""Revision analyses R-1..R-4 (src/part2/revision.py) on synthetic data with known answers: the natural-spline basis
(column space of ns(df = 4); linear beyond the boundary), the spline trend and the uniform scale-artefact / lead rule,
pool-matched sets (tolerances, determinism, unavailable units, generators), the discovery selection, the symmetric
flag and conditional P(real > random), the single HC3 fit and the power / standardised-rate helpers."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import stats as sst

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from part2 import revision as rv, stats  # noqa: E402


# ------------------------------------------------------------------------------------------------ R-2
def test_ns_basis_spans_natural_cubic_splines():
    rng = np.random.default_rng(0)
    x = np.sort(rng.normal(size=300))
    kn = rv.ns_knots(x, 4)
    assert len(kn) == 5 and kn[0] == x.min() and kn[-1] == x.max()
    assert np.allclose(kn[1:-1], np.quantile(x, [0.25, 0.5, 0.75]))
    B = rv.ns_basis(x, kn)
    assert B.shape == (300, 5) and np.linalg.matrix_rank(B) == 5
    # linear and constant functions are reproduced exactly
    for y in (np.ones_like(x), 2 - 3 * x):
        c, *_ = np.linalg.lstsq(B, y, rcond=None)
        assert np.allclose(B @ c, y)
    # every basis function is a cubic spline with knots kn: piecewise cubic with continuous 2nd derivative ...
    xx = np.linspace(kn[0] - 1, kn[-1] + 1, 4001)
    Bx = rv.ns_basis(xx, kn)
    h = xx[1] - xx[0]
    d2 = np.diff(Bx, 2, axis=0) / h ** 2
    assert np.abs(np.diff(d2, axis=0)).max() < 1e-2
    # ... and linear beyond the boundary knots (natural)
    out = (xx[1:-1] < kn[0] - 0.01) | (xx[1:-1] > kn[-1] + 0.01)
    assert np.abs(d2[out]).max() < 1e-6


def test_spline_trend_absorbs_curvature_and_keeps_outliers():
    rng = np.random.default_rng(1)
    rows = []
    for c in ("lineage:A", "lineage:B"):
        e = rng.uniform(-2, 0.2, 300)
        beta = 0.05 * e ** 2 + rng.normal(0, 0.01, 300)            # curved trend: linear fit leaves structure
        beta[7] += 0.2                                              # planted outlier
        rows.append(pd.DataFrame(dict(context=c, unit=[f"U{i}" for i in range(300)], beta=beta, evaluable=True,
                                      e=e)))
    t = pd.concat(rows, ignore_index=True)
    ess = t.drop_duplicates("unit").set_index("unit").e        # one essentiality per unit (context A's)
    mb = t.context == "lineage:B"
    t.loc[mb, "beta"] = 0.05 * t.unit[mb].map(ess).values ** 2 + rng.normal(0, 0.01, int(mb.sum()))
    by, per = rv.spline_trend(t.drop(columns="e"), ess)
    a = by.set_index("context").loc["lineage:A"]
    assert a.r2_nl > a.r2_lin and a.resid_sd_nl < a.resid_sd_lin
    pa = per[per.context == "lineage:A"].set_index("unit")
    assert pa.residual_z_nl["U7"] > 10 and abs(pa.residual_z_nl.drop("U7")).max() < 5
    assert set(per.columns) >= {"fitted_nl", "residual_z_nl", "concordant_nl", "residual_z_lin"}


def test_scale_labels_uniform_rule():
    lab = rv.scale_labels(rz_nl=[2.5, 1.0, 1.0, -3.0, 2.2, 2.2],
                          concordant_nl=[True, True, False, False, True, True],
                          tier=["T1", "T1", "T1", "T2", "T1", "T1"],
                          control_recovery=["", None, "", "", "readout (PC20)", "same unit, other control context (x)"],
                          direction=["H3a"] * 6)
    assert lab.scale_artefact_possible.tolist() == [False, True, False, False, False, False]
    assert lab.lead_corrected.tolist() == [True, False, False, False, False, True]


# ------------------------------------------------------------------------------------------------ R-1(a)
def _pool_world(seed=2, n=200):
    rng = np.random.default_rng(seed)
    genes = [f"G{i}" for i in range(n)]
    ess = pd.Series(rng.uniform(-1.0, 0.0, n), index=genes)
    F = rng.normal(size=(n, 5))
    S = F @ F.T / 2
    np.fill_diagonal(S, np.nan)
    return genes, ess, pd.DataFrame(S, index=genes, columns=genes)


@pytest.mark.parametrize("gen", rv.GENERATORS)
def test_pool_matched_sets_respect_tolerances_and_seed(gen):
    genes, ess, ps = _pool_world()
    pool, units_g = genes[40:], genes[:40]
    unit_members = {"U0": units_g[:4], "U1": units_g[4:10], "U2": units_g[10:13], "BIG": genes[:10] + genes[40:200]}
    units = list(unit_members)
    sets, st, q = rv.pool_matched_sets(pool, ess, ps, unit_members, units, R=5, seed=11, generator=gen,
                                       max_attempts=4000)
    st = st.set_index("unit")
    assert st.loc["BIG", "status"] == "unavailable" and st.loc["BIG", "n_sets"] == 0   # k > pool size
    for u in ("U0", "U1", "U2"):
        if st.loc[u, "status"] == "unavailable":
            continue
        mult = 1.0 if st.loc[u, "status"] == "ok" else 2.0
        mine = {k: v for k, v in sets.items() if k.startswith(f"{u}#")}
        assert len(mine) == 5 and len({tuple(v) for v in mine.values()}) == 5
        for v in mine.values():
            assert len(v) == len(unit_members[u]) and set(v) <= set(pool)
            assert abs(ess[v].mean() - st.loc[u, "target_ess"]) <= 0.1 * mult + 1e-12
            if gen != "ess":
                M = ps.loc[v, v].values
                coh = np.nanmean(M[np.triu_indices(len(v), 1)])
                assert abs(coh - st.loc[u, "target_coh"]) <= 1.0 * mult + 1e-12
    assert (st.n_sets > 0).sum() >= 1
    sets2, _, _ = rv.pool_matched_sets(pool, ess, ps, unit_members, units, R=5, seed=11, generator=gen,
                                       max_attempts=4000)
    assert sets2 == sets
    assert set(q.columns) == {"unit", "set", "size", "d_ess", "d_coh"}


def test_pool_nbr_sets_are_neighbourhoods():
    genes, ess, ps = _pool_world(n=120)
    pool = genes[20:]
    sets, st, _ = rv.pool_matched_sets(pool, ess, ps, {"U": genes[:3]}, ["U"], R=3, seed=1, generator="nbr",
                                       tol_ess=10, tol_coh=100, max_attempts=50)
    P = ps.loc[pool, pool]
    for v in sets.values():                      # every non-seed member is among some member's top-10 neighbours
        tops = {g: set(P.loc[g].drop(g).sort_values(ascending=False).index[:10]) for g in v}
        assert any(all(h in tops[g] for h in v if h != g) for g in v)


# ------------------------------------------------------------------------------------------------ R-1(b)
def test_discover_equals_bb_when_permutation_returns_analytic_p():
    rng = np.random.default_rng(4)
    p = np.concatenate([rng.uniform(size=300), rng.uniform(0, 1e-4, 20)])
    res = pd.DataFrame(dict(context=np.repeat(["a", "b", "c", "d"], 80), p=p, evaluable=True))
    out = rv.discover(res, 0.05, lambda idx: {i: res.p[i] for i in idx})
    bb = stats.benjamini_bogomolov(res.p.values, res.context.values, 0.05)
    assert (out.discovery.values == bb["reject"]).all() and out.discovery.sum() > 0
    out2 = rv.discover(res, 0.05, lambda idx: {i: 1.0 for i in idx})          # permutation removes every hit
    assert out2.discovery.sum() == 0 and out2.p_perm.notna().sum() == (out.q_adj < 0.1).sum()


def test_symmetric_flag_and_conditional_auc():
    real = pd.DataFrame(dict(context=["c", "c"], unit=["U", "V"], abs_t=[5.0, 1.0]))
    rand = pd.DataFrame(dict(context=["c"] * 4, unit=["U", "U", "V", "V"], rep=[0, 1, 0, 1],
                             abs_t=[6.0, 2.0, 0.5, 0.9]))
    f = rv.symmetric_flag(real, rand)
    assert f.exceeds_all_others.tolist() == [True, False, False, False]
    assert f.max_abs_t_others.tolist() == [5.0, 6.0, 1.0, 1.0]
    assert rv.prob_greater(3.0, np.array([1, 3, 5])) == pytest.approx(0.5)
    a = rv.conditional_auc(real, rand, "abs_t", ("context",))
    assert a.tolist() == [0.75, 0.5]
    b = rv.conditional_auc(real, rand, "abs_t", ("context", "unit"))
    assert b.tolist() == [0.5, 1.0]


# ------------------------------------------------------------------------------------------------ R-3 / R-4
def test_fit_one_hc3_matches_engine_and_min_n():
    rng = np.random.default_rng(5)
    idx = pd.Index([f"M{i}" for i in range(80)])
    ctx = pd.Series((np.arange(80) % 4 == 0).astype(float), index=idx)
    cov = pd.DataFrame({"x": rng.normal(size=80), "const": 1.0}, index=idx)
    y = pd.Series(rng.normal(size=80) - 0.8 * ctx, index=idx)
    y.iloc[3] = np.nan
    f = rv.fit_one_hc3(y, ctx, cov)
    m = y.notna().values
    X = np.column_stack([np.ones(m.sum()), ctx.values[m], cov.x.values[m]])
    ref = stats.ols_hc3(y.values[m], X, 1)
    assert f["beta"] == pytest.approx(ref["beta"]) and f["t"] == pytest.approx(ref["t"])
    assert f["n"] == 79 and f["n_context"] == 20 and f["n_covariates"] == 1
    few = rv.fit_one_hc3(y, ctx * (np.arange(80) < 8), cov)                    # 2 context models
    assert np.isnan(few["beta"]) and few["n_context"] == 2


def test_power_and_standardised_rate():
    z = np.array([4.0, 4.0])
    se_A, se_L, df_L, b = np.array([0.02, 0.02]), np.array([0.04, 0.04]), np.array([30.0, 30.0]), np.array([.08, .08])
    p1 = rv.power_eb(z, np.array([1.0, 1.0]), np.array([0.0, 0.0]), se_A, se_L, df_L, b)
    assert np.allclose(p1["naive"], p1["eb_mean"])                    # no shrinkage: EB mean = naive
    ncp = 0.08 / 0.04
    assert p1["naive"][0] == pytest.approx(sst.nct.sf(sst.t.ppf(0.95, 30), 30, ncp))
    p2 = rv.power_eb(z, np.array([0.5, 0.5]), np.array([1.0, 1.0]), se_A, se_L, df_L, b)
    assert (p2["eb_mean"] < p1["eb_mean"]).all()
    w = pd.Series({"a": 0.25, "b": 0.75})
    r, sd = rv.standardised_rate(np.array([1, 0, 1, 1]), np.array(["a", "a", "b", "b"]), w)
    assert r == pytest.approx(0.25 * 0.5 + 0.75 * 1.0) and sd == 0.0
    r, sd = rv.standardised_rate(np.array([0.5, 0.5]), np.array(["a", "b"]), w)
    assert r == pytest.approx(0.5) and sd == pytest.approx(np.sqrt(0.25 * (0.25 ** 2 + 0.75 ** 2)))
