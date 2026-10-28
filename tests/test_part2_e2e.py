"""End-to-end synthetic test of the Part 2 association engine (DESIGN v3.1 §3): planted complexes, lineage and
genotype effects (incl. a confounded lineage effect), a sub-module divergence, H1 paralog-compensation / collapse /
sensitisation / neutral events and a PatientID-blocked lockbox split, run through the gated pipeline
(tests/part2_synth_pipeline.py). Checks recovery of planted effects, H1 classes, lockbox order and FDR under the null
(mean FDP over 50 null simulations per scenario). No DepMap data is read."""
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy import stats as sst

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))

import part2_synth_pipeline as P  # noqa: E402
from part2 import registry, synth  # noqa: E402

Q = 0.05


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    return P.run_pipeline(synth.simulate(1), tmp_path_factory.mktemp("e2e") / "repo", n_perm=P.N_PERM, K=1000,
                          B_bench=100), tmp_path_factory


def test_planted_effects_recovered(run):
    r, _ = run
    assert r["kill_rule"] == {"H3a": "passes", "H3b": "passes", "H1": "passes", "H2": "passes"}
    for d in ("H3a", "H3b"):
        assert r[d]["power"] == 1.0, r[d]                 # every planted (context, complex) pair discovered
        assert r[d]["false_pos"] <= 1, r[d]
    assert r["H2"]["power"] == 1.0 and r["H2"]["false_pos"] == 0, r["H2"]
    assert r["H3b"]["true_pos"] == 5


def test_h1_classes(run):
    r, _ = run
    calls = r["H1"]["calls"]
    want = {"PARALOG COMPENSATION", "COLLAPSE", "RESIDUAL SENSITISATION"}
    for k, (got, exp) in calls.items():
        if exp in want:
            assert got == exp, (k, got, exp)                 # every planted non-neutral class called
        else:
            assert got not in want, (k, got, exp)            # no planted-neutral axis gets a non-neutral class
    assert r["H1"]["accuracy"] >= 6 / 7
    ci = r["H1"]["collapse_index"]
    assert 0.4 < ci["C0011|GLOSB"] < 0.8 and ci["C0012|GLOSC"] < -0.3     # planted 60% collapse; sensitisation


def test_lockbox_order_and_confirmation(run):
    r, _ = run
    assert r["bench_gate_before_confirmations"].startswith("benchmark lockbox read before")
    conf = r["confirmation"]
    assert set(conf) == {"CTRL_PC1", "CTRL_PC20", "CTRL_PC21", "CTRL_PC6", "H3a_lineage_atlas", "H3b_driver_family",
                         "H2_submodule_divergence", "H1_within_complex_GI"}
    assert all(v["n_confirmed"] == v["n_hits"] for k, v in conf.items() if k.startswith("CTRL"))
    assert conf["H3b_driver_family"]["n_confirmed"] >= 4
    assert "bench" in r and np.isfinite(r["bench"]["mean"])     # BENCH read only after every entry closed
    assert r["runtime_s"]["total"] < 600


def _binom_ok(x: int, n: int, q: float = Q) -> bool:
    """FDP under a complete null is 1{any discovery}: E[FDP] <= q is rejected only if P(X >= x | n, q) < 0.01."""
    return sst.binom.sf(x - 1, n, q) >= 0.01


def test_fdr_global_null_50():
    df = P.fdr_sims(50, seed0=10_000, null=True, lineage_null=True)
    for d in ("H3a", "H3b"):
        fdp = df[f"{d}_fdp"]
        assert _binom_ok(int((fdp > 0).sum()), len(df)), (d, fdp.mean())
        assert fdp.mean() <= 2 * Q, (d, fdp.mean())


def test_fdr_h3b_null_with_confounded_lineage_effects_50():
    """No genotype effects, but planted lineage effects incl. Pancreas on a complex while KRAS is enriched in
    Pancreas: lineage_FE must keep H3b at the null; H3a must find the lineage effects without false units.
    H3a's registered error rate is Benjamini-Bogomolov's (families = contexts): the expected average, over selected
    context families, of the within-family FDP. The pooled FDP over all H3a discoveries is not controlled by BB (each
    family's BH admits ~a(T_f + 1) false discoveries for T_f true ones; ~1.4 true per family here gives pooled
    ~0.066 with exact p-values, results/part2/verify/h3a_fdr_diagnosis.md), so it is reported, not bounded."""
    df = P.fdr_sims(50, seed0=20_000, null=True, lineage_null=False)
    assert _binom_ok(int((df.H3b_fdp > 0).sum()), len(df)), df.H3b_fdp.mean()
    assert df.H3a_fdp_bb.mean() <= Q + 0.02, (df.H3a_fdp_bb.mean(), df.H3a_fdp.mean())
    assert df.H3a_power.mean() >= 0.95, df.H3a_power.mean()
