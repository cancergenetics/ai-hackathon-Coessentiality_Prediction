"""H1 (loss calls, event enumeration, member-resolved scores, alignment sign, matched null, classes) and H2
(residualised sub-module score, A-fitted coefficients re-used on lockbox rows) and the controls kill tally, on toy /
synthetic data with known answers."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from part2 import controls, covariates as cv, h1, h2, synth  # noqa: E402
from part2.cmap import FrozenMap  # noqa: E402
from part2.cscore import Annotation, Removal  # noqa: E402


def test_loss_calls_rules():
    idx = [f"m{i}" for i in range(8)]
    g = ["G"]
    lof = pd.DataFrame({"G": [1, 0, 0, 0, 0, 0, np.nan, 1]}, index=idx, dtype=float)
    hd = pd.DataFrame({"G": [0, 1, 0, 0, 0, 0, 0, 0]}, index=idx, dtype=float)
    dmg = pd.DataFrame({"G": [1, 0, 0, 0, 1, 0, 1, 1]}, index=idx, dtype=float)
    ex = pd.DataFrame({"G": [5, 0.5, 0.2, 5, 5, np.nan, 5, 5]}, index=idx, dtype=float)
    lin = pd.Series(["a"] * 8, index=idx)
    loss, kind = h1.loss_calls(lof, hd, dmg, ex, lin)
    # m0 LoF; m1 hom-del; m2 silenced (lineage median 5 >= 3, expr < 1, no hom-del); m3 reference; m4 single
    # damaging hit without LoF call -> NaN; m5 expression not measured -> reference; m6 MSI-H (LoF NaN) -> NaN
    assert loss.G.fillna(-1).tolist() == [1, 1, 1, 0, -1, 0, -1, 1]
    assert kind.G.tolist()[:3] == ["LoF", "homdel", "silencing"]
    loss2, _ = h1.loss_calls(lof, hd, dmg, ex, lin, lof_pass=set())                 # fails functional check
    assert np.isnan(loss2.G["m0"]) and loss2.G["m1"] == 1
    assert h1.event_type(pd.Series(["homdel", "LoF"])) == "CN"                        # tie -> CN
    assert h1.event_type(pd.Series(["LoF", "LoF", "homdel"])) == "LoF"
    assert h1.event_type(pd.Series(["silencing", "LoF", "homdel"])) == "silencing"


def test_enumerate_events():
    keys = {s: f"{s} ({i})" for i, s in enumerate(["G", "PAR", "A", "B", "NEAR", "X", "Y", "Z"])}
    fm = FrozenMap(members={"C1": [keys[s] for s in ("G", "PAR", "A", "B", "NEAR")],
                            "C2": [keys[s] for s in ("X", "Y", "Z")]},
                   level={"C1": 1, "C2": 1}, attached={"C2": [keys["G"]]})
    an = Annotation(paralogs={"G": {"PAR"}, "PAR": {"G"}}, chrom={"G": "1", "NEAR": "1", "A": "2", "B": "3"},
                    pos={"G": 1_000_000, "NEAR": 5_000_000, "A": 1, "B": 1})
    n = 100
    ids = [f"m{i}" for i in range(n)]
    loss = pd.DataFrame({"G": np.r_[np.ones(20), np.zeros(80)]}, index=ids)
    fold = pd.Series(np.where(np.arange(n) % 4 == 0, "L", "A"), index=ids)
    for k, exp_N in (("LoF", ["A (2)", "B (3)", "NEAR (4)"]), ("homdel", ["A (2)", "B (3)"])):
        kind = pd.DataFrame({"G": np.where(loss.G == 1, k, "")}, index=ids)
        ev = h1.enumerate_events(fm, loss, kind, an, fold)
        e1 = next(e for e in ev if e.unit == "C1")
        assert e1.P == ["PAR (1)"] and e1.N == exp_N                    # CN-type: member within 10 Mb dropped
        assert e1.n_A == 15 and e1.n_L == 5 and e1.testable
        e2 = next(e for e in ev if e.unit == "C2")
        assert e2.link == "attached" and e2.P == [] and len(e2.N) == 3
    assert [e.index for e in ev] == [0, 1]
    assert h1.units_of(ev) == sorted(["C1|G|N", "C1|G|P", "C2|G|N"])
    few = loss.copy()
    few.iloc[:10, 0] = 0
    assert h1.enumerate_events(fm, few, kind, an, fold) == []             # < 12 loss models in A


def test_alignment_sign_and_collapse_index():
    """kappa_N > 0 when members' dependency grows in its baseline direction (registry statement); < 0 for collapse."""
    rng = np.random.default_rng(0)
    n = 400
    loss = pd.Series((rng.random(n) < 0.3).astype(float))
    base = np.array([-1.0, -0.8, 0.6])                                   # two negative, one positive baseline
    for grow, sign in ((+0.3, 1), (-0.5, -1)):
        E = pd.DataFrame(base + rng.normal(0, 0.1, (n, 3)), columns=["a (1)", "b (2)", "c (3)"])
        E += np.outer(loss, grow * np.sign(base))                       # magnitude +grow in baseline direction
        b = h1.baselines(E, list(E.columns), loss)
        y = h1.aligned_score(E, list(E.columns), b)
        from part2.robust import fit_one
        k = fit_one(y, loss, None)["beta"]
        assert np.sign(k) == sign and k == pytest.approx(grow, abs=0.03)


def test_matched_null_and_classes_on_synthetic():
    sim = synth.simulate(1)
    z = sim.lof * 0
    loss, kind = h1.loss_calls(sim.lof, z, sim.lof, z + 5.0, sim.meta.OncotreeLineage)
    events = h1.enumerate_events(sim.fmap, loss, kind, sim.annot, sim.fold)
    A = sim.fold.index[sim.fold == "A"]
    E = sim.E.loc[A]
    ev = next(e for e in events if e.gene == "ARID1A")
    c = h1.component_scores(E, ev, loss.loc[A, "ARID1A"])
    ess = E.clip(-3, 3).mean()
    ps = 10 * E[sorted(sim.fmap.all_genes(1))].clip(-3, 3).corr()
    tab = pd.concat([sim.contexts, sim.meta.drop(columns=["TP53_altered", "WGD"])], axis=1).loc[A]
    cov = cv.assemble(ev.covariates, tab)
    kw = dict(annot=sim.annot, essentiality=ess, pair_scores=ps, K=300)
    n1 = h1.matched_null(E, sim.fmap, ev, "P", c, loss.loc[A, "ARID1A"], cov, **kw)
    n2 = h1.matched_null(E, sim.fmap, ev, "P", c, loss.loc[A, "ARID1A"], cov, **kw)
    assert n1["status"] in ("ok", "tolerances doubled") and n1["n_sets"] == 300
    assert np.array_equal(n1["null"], n2["null"])                       # seed 7000 + event index: deterministic
    assert abs(np.median(n1["null"])) < 0.05 and n1["q025"] < 0 < n1["q975"]
    tiny = FrozenMap(members={u: sim.fmap.members[u] for u in [ev.unit, sim.fmap.coarse()[-1]]},
                     level={ev.unit: 1, sim.fmap.coarse()[-1]: 1})
    assert h1.matched_null(E, tiny, ev, "P", c, loss.loc[A, "ARID1A"], cov, max_attempts=2000,
                           **kw)["status"] == "null unavailable"


def test_h2_residualised_and_lockbox_reuse():
    rng = np.random.default_rng(4)
    n = 500
    keys = [f"g{i} ({i})" for i in range(8)]
    f = rng.normal(size=n)
    X = -0.5 + np.outer(f, np.ones(8)) * 0.4 + rng.normal(0, 0.2, (n, 8))
    ctx = (rng.random(n) < 0.3).astype(float)
    X[:, :4] += -0.5 * ctx[:, None]                                       # sub-module S1 = g0..g3 diverges
    E = pd.DataFrame(X, columns=keys)
    fm = FrozenMap(members={"P": keys, "P.S1": keys[:4], "P.S2": keys[4:]}, level={"P": 1, "P.S1": 2, "P.S2": 2},
                   parent={"P": "", "P.S1": "P", "P.S2": "P"})
    assert h2.h2_units(fm) == ["P.S1", "P.S2"]
    y, coef, (Ss, Sb) = h2.residualised(E.iloc[:350], fm, "P.S1")
    from part2.robust import fit_one
    assert fit_one(y, pd.Series(ctx[:350]), None)["beta"] < -0.3
    yL, coefL, _ = h2.residualised(E.iloc[350:], fm, "P.S1", coef=coef)
    assert coefL == coef
    a, b = coef
    assert np.allclose(yL, cscore_mean(E.iloc[350:], keys[:4]) - (a + b * cscore_mean(E.iloc[350:], keys[4:])))
    Y, coefs = h2.h2_frames(E, fm, ["cx"], lambda c: Removal(genes=("g0",), paralogs=False), Annotation())
    assert coefs[("cx", "P.S1")] is not None and Y["cx"].shape == (n, 2)
    f_L = h2.score_L_fn(fm, {("cx", "P.S1"): None}, lambda c: Removal(), Annotation())
    assert f_L(E, "cx", "P.S1").isna().all()                             # never refitted on lockbox rows


def cscore_mean(E, g):
    from part2.cscore import mean_score
    return mean_score(E, g)


def test_kill_tally():
    tab = pd.DataFrame([
        dict(design_id="PC1", direction="H3b", counted=True, level="C", evaluable_kill=True, recovered=True),
        dict(design_id="PC2", direction="H3b", counted=True, level="C", evaluable_kill=True, recovered=False),
        dict(design_id="PC4", direction="H3b", counted=True, level="C", evaluable_kill=True, recovered=False),
        dict(design_id="PC17", direction="H3a", counted=True, level="C", evaluable_kill=True, recovered=True),
        dict(design_id="PC20", direction="H3a", counted=True, level="C", evaluable_kill=False, recovered=False),
        dict(design_id="PC6", direction="H1", counted=True, level="C", evaluable_kill=False, recovered=False),
        dict(design_id="PC3", direction="none", counted=False, level="G", evaluable_kill=False, recovered=False),
    ])
    k = controls.kill_tally(tab)
    assert k["H3b"]["verdict"] == "killed by controls" and k["H3b"]["fraction"] == pytest.approx(1 / 3)
    assert k["H3a"]["verdict"] == "passes" and k["H3a"]["thin_coverage"]
    assert k["H1"]["verdict"].startswith("n/a") and k["H2"]["evaluable"] == 0
