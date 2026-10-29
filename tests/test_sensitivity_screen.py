"""scripts/48_sensitivity_screen.py on the synthetic project (no DepMap value is read): the screen-quality covariates
are appended to the registered ones, complete-case dropout is counted, the trial ledger is untouched and the run is
refused a second time."""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))

import part2_synth_pipeline as P  # noqa: E402
from part2 import assoc, covariates as covmod, registry, synth  # noqa: E402

spec = importlib.util.spec_from_file_location("sens_screen", ROOT / "scripts" / "48_sensitivity_screen.py")
S = importlib.util.module_from_spec(spec)
sys.modules["sens_screen"] = S
spec.loader.exec_module(S)
R = S.R

N_PERM = 200


@pytest.fixture(scope="module")
def setup(tmp_path_factory):
    sim = synth.simulate(0)
    rep = P.build_repo(tmp_path_factory.mktemp("sens") / "repo", sim)
    proj = rep["proj"]
    tab = pd.read_parquet(proj.root / P.REL_CTX)
    E_A = assoc.load_A_effects(proj=proj)
    A = E_A.index
    env = R.Env(proj=proj, fmap=sim.fmap, ctx=tab.loc[A], E_A=E_A, cov=tab.loc[A], annot=sim.annot,
                inputs=rep["inputs"], blocks=tab.PatientID.loc[A], strata=tab.OncotreeLineage.loc[A],
                strict=False, write=True)
    prim = R.run_h3(env, "H3a_lineage_atlas", n_perm=N_PERM)
    return env, prim


def test_covariates_are_registered_plus_screen_quality(setup):
    env, _ = setup
    entry, _ = R.entry_of(env, "H3a_lineage_atlas")
    ctxs = R.primary_contexts(entry)
    prim, _ = covmod.per_context(entry, ctxs, env.cov)
    sens, names = S.sensitivity_covariates(entry, ctxs, env.cov)
    for c in ctxs:
        assert names[c] == covmod.covariates_for(entry, c)[0] + ["NNMD", "CasActivity"]
        assert list(sens[c].columns) == list(prim[c].columns) + ["NNMD", "CasActivity"]
        pd.testing.assert_frame_equal(sens[c][prim[c].columns], prim[c])
        assert np.allclose(sens[c]["NNMD"].values, env.cov.ScreenNNMD.astype(float).values)
        assert sens[c]["CasActivity"].isna().equals(env.cov.CasActivity.isna())


def test_run_writes_no_trial_line_and_counts_dropout(setup, tmp_path):
    env, prim = setup
    aid = "H3a_lineage_atlas"
    trials = env.proj.trials.read_bytes()
    pcal = prim["calibration"]
    od = tmp_path / "out"
    out = S.run(env, aid, R.with_ci(prim["res"], env.fmap), pcal, od, n_perm=N_PERM)
    assert env.proj.trials.read_bytes() == trials                          # no trial line
    recs = [json.loads(l) for l in (od / "engine_record.jsonl").read_text().splitlines()]
    assert len(recs) == 1 and recs[0]["tag"] == S.LABEL
    s = json.loads((od / "summary.json").read_text())
    assert s["label"] == S.LABEL and "not a trial" in s["standing"]
    # complete-case dropout = models with CasActivity / NNMD missing (registered covariates are complete here)
    miss = env.cov.ScreenNNMD.isna() | env.cov.CasActivity.isna()
    assert s["dropout_complete_case"]["models_dropped_any_context"] == int(miss.sum()) > 0
    d = out["dropout_context"].set_index("context")
    for c in d.index:
        x = env.ctx[c]
        assert d.loc[c, "dropped_context"] == int((miss & (x == 1)).sum())
        assert d.loc[c, "n_ref_sens"] == int((~miss & (x == 0)).sum())
    # the sensitivity design has the two extra columns and fewer models
    r = out["res"]
    assert (r.n <= len(env.E_A) - miss.sum()).all()
    c = out["comparison"]
    assert int(c.surviving.sum()) == s["n_surviving"] and int(c.disc_primary.sum()) == int(prim["res"].discovery.sum())
    assert s["n_surviving"] + s["n_lost"] == s["n_discovery_primary"]
    assert s["n_surviving"] + s["n_new"] == s["n_discovery"]
    assert s["pseudo_context"]["primary"] == pcal["frac_p_lt_005"]
    assert out["rho"] is not None and len(out["rho"]) == len(c)
    with pytest.raises(registry.RegistryError):                            # runs once
        S.run(env, aid, R.with_ci(prim["res"], env.fmap), pcal, od, n_perm=N_PERM)


def test_complete_covariates_reproduce_registered_fit(setup, tmp_path):
    """With NNMD / CasActivity constant the design builder drops them, so the sensitivity equals the registered run
    (plumbing adds nothing else)."""
    env, prim = setup
    env2 = R.Env(**{**env.__dict__, "cov": env.cov.assign(ScreenNNMD=1.0, CasActivity=2.0)})
    out = S.run(env2, "H3a_lineage_atlas", R.with_ci(prim["res"], env.fmap), prim["calibration"], tmp_path / "o",
                n_perm=N_PERM)
    a = prim["res"].set_index(["context", "unit"])
    b = out["res"].set_index(["context", "unit"])
    assert np.allclose(a.beta, b.beta.reindex(a.index)) and np.allclose(a.p, b.p.reindex(a.index))
    assert (a.discovery == b.discovery.reindex(a.index)).all()
    assert out["summary"]["n_lost"] == 0 and out["summary"]["n_new"] == 0
    assert out["summary"]["pseudo_context"]["change"] == 0
