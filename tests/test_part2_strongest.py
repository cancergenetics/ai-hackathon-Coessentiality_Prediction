"""Registered robustness and the strongest-findings rule (src/part2/strongest.py) on synthetic data with known answers:
LOSO by the registered grouping (a hit carried by one group fails, a broad hit passes), leave-top-2-out (two extreme
context models), candidate selection (q, primary call incl. rare route, R7 defect-sensitive, locus, H1 E2), rule
application and ranking (tier class first, then |t_A|)."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from part2 import strongest  # noqa: E402


def _world(seed=3, n=700):
    rng = np.random.default_rng(seed)
    idx = pd.Index([f"M{i}" for i in range(n)])
    lin = pd.Series(np.repeat(list("abcdefg"), n // 7), index=idx)
    disease = lin + pd.Series(rng.integers(0, 3, n).astype(str), index=idx)       # 21 primary diseases
    ctx = pd.Series((rng.random(n) < 0.25).astype(float), index=idx)
    cov = pd.DataFrame({"library:x": (rng.random(n) < 0.5).astype(float)}, index=idx)
    return rng, idx, lin, disease, ctx, cov


def test_loso_and_top2_known_answers():
    rng, idx, lin, disease, ctx, cov = _world()
    broad = pd.Series(rng.normal(size=len(idx)) - 0.6 * ctx, index=idx)
    one = pd.Series(rng.normal(size=len(idx)) - 2.5 * ctx * (lin == "c"), index=idx)
    r_b = strongest.hit_robustness(broad, ctx, cov, lin)
    r_o = strongest.hit_robustness(one, ctx, cov, lin)
    assert r_b["loso_robust"] and r_b["loso_folds"] == 7 and r_b["top2_robust"]
    assert r_o["p_A_refit"] < 0.05 and not r_o["loso_robust"]                  # vanishes when lineage c is left out
    r_d = strongest.hit_robustness(broad, ctx, cov, disease)                   # H3a unit: primary disease
    assert r_d["loso_folds"] == 21 and r_d["loso_robust"]
    # leave-top-2-out drops the two context models most extreme in the effect's direction (and only those)
    out = pd.Series(rng.normal(size=len(idx)), index=idx)
    top = ctx[ctx == 1].index[:2]
    out[top] = -60.0
    r_t = strongest.hit_robustness(out, ctx, cov, lin)
    ref = strongest.robust.fit_one(out.drop(top), ctx, cov)
    assert r_t["beta_A_refit"] < 0 and np.isclose(r_t["top2_beta"], ref["beta"]) and not r_t["top2_robust"]
    pos = (-out).copy()
    assert np.isclose(strongest.hit_robustness(pos, ctx, cov, lin)["top2_beta"], -ref["beta"])
    # a group holding every context model: its fold is skipped and counted, not failed
    ctx2 = pd.Series(((lin == "a") & (pd.Series(rng.random(len(idx)), index=idx) < 0.5)).astype(float))
    y2 = pd.Series(rng.normal(size=len(idx)) - 1.0 * ctx2, index=idx)
    r_s = strongest.hit_robustness(y2, ctx2, cov, lin)
    assert r_s["loso_skipped"] == 1 and r_s["loso_folds"] == 6


def test_candidates_rule_and_ranking():
    disc = pd.DataFrame(dict(context=["c1"] * 6, unit=[f"U{i}" for i in range(6)],
                             beta=[-0.3, 0.2, -0.1, 0.4, -0.2, 0.3], t=[-9, 5, -12, 7, -4, 3],
                             q_final=[1e-5, 0.01, 0.001, 0.07, 0.02, 0.03],
                             tier=["T3-neg", "T1", "T2", "T1", "T3-nc", "T1"],
                             locus=[False, False, False, False, True, False]))
    calls = pd.DataFrame(dict(context=["c1"] * 6, unit=disc.unit,
                              primary_call=["confirmed", "confirmed", "not lockbox-testable: rare-context route passed",
                                            "confirmed", "confirmed", "not confirmed"],
                              defect_sensitive=["False"] * 6))
    calls.loc[1, "defect_sensitive"] = "True"
    x = strongest.h3_candidates(disc, calls)
    assert x.candidate.tolist() == [True, False, True, False, False, False]  # R7, q, locus, not confirmed
    assert x.route[2].startswith("rare-context")
    x["loso_robust"] = [True, True, True, True, True, False]
    r = strongest.apply_rule(x)
    s = r[r.strongest]
    assert s.unit.tolist() == ["U2", "U0"] and s["rank"].tolist() == [1, 2]   # T2 before T3 despite smaller |t|
    x.loc[0, "loso_robust"] = False
    assert strongest.apply_rule(x).strongest.sum() == 1


def test_h1_candidates_respect_e2_and_eligibility():
    calls = pd.DataFrame(dict(event=["Z1|A", "Z2|B", "Z3|C", "Z4|D"], axis=["P", "P", "N", "N"],
                              eligible=[True, True, True, False],
                              primary_call=["confirmed", "not lockbox-testable: rare-context route passed",
                                            "confirmed", "confirmed"],
                              defect_sensitive=[False, False, False, False]))
    disc = pd.DataFrame(dict(unit=["Z1|A|P", "Z2|B|P", "Z3|C|N", "Z4|D|N"], beta=[-0.3, -1, 0.05, 0.1],
                             t=[-5, -6, 4, 3], q_final=[1e-4, 1e-4, 1e-3, 1e-3], tier=["T1", "T3-neg", "T1", "T1"]))
    ev = pd.DataFrame(dict(event=["Z1|A", "Z2|B", "Z3|C", "Z4|D"], locus=[False] * 4))
    head = pd.DataFrame(dict(unit=disc.unit, headline_eligible=[True, True, False, True]))
    x = strongest.h1_candidates(disc, calls, ev, head)
    assert x.unit.tolist() == ["Z1|A|P", "Z2|B|P", "Z3|C|N"]                   # ineligible call dropped
    assert x.candidate.tolist() == [True, True, False]                         # E2: not headline-eligible


def test_post_hoc_annotations():
    ctl = pd.DataFrame(dict(design_id=["PC6", "PC20"], context=["LoF:ARID1A", "heme"], unit=["Z1", "Z2"]))
    assert strongest.control_recovery("H1loss:ARID1A", "Z1|ARID1A|P", ctl) == "readout (PC6)"
    assert strongest.control_recovery("lineage:Myeloid", "Z2", ctl) == "readout (PC20)"
    assert strongest.control_recovery("lineage:Bone", "Z2", ctl).startswith("same unit, other control context (PC20")
    assert strongest.control_recovery("lineage:Bone", "Z9", ctl) == ""
    ess = strongest.unit_essentiality(pd.Series({"a": -2.0, "b": 0.0, "c": -1.0}), {"U": ["a", "b"], "V": ["c", "x"]})
    assert ess["U"] == -1.0 and ess["V"] == -1.0
    rng = np.random.default_rng(5)
    e = pd.Series(rng.uniform(-2, 0, 200), index=[f"U{i}" for i in range(200)])
    beta = -0.1 * e + rng.normal(0, 0.01, 200)                     # essential units look less essential
    beta.iloc[0] = 1.0                                              # one unit far off the line
    tests = pd.DataFrame(dict(context="lineage:x", unit=e.index, beta=beta.values, evaluable=True))
    by, per = strongest.scale_effect(tests, e)
    assert np.isclose(by.slope.iloc[0], -0.1, atol=0.02) and by.spearman_rho.iloc[0] < -0.8
    p = per.set_index("unit")
    assert p.loc["U0", "residual_z"] > 5 and p.residual_z.drop("U0").abs().max() < 4
    assert p.drop("U0").scale_effect_concordant.mean() > 0.9
