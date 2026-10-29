"""scripts/55_export_member_profiles.py on synthetic data: the export reuses the engine's per-member fit, the aligned
beta follows ALIGN = +1, and the per-event breadth recomputed from the export equals h1.classify's breadth (the check
that guards the real export against the committed classes.tsv)."""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from part2 import covariates as cv, h1, synth  # noqa: E402

_spec = importlib.util.spec_from_file_location("export55", ROOT / "scripts" / "55_export_member_profiles.py")
exp = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = exp
_spec.loader.exec_module(exp)


@pytest.fixture(scope="module")
def setup():
    sim = synth.simulate(1)
    z = sim.lof * 0
    loss, kind = h1.loss_calls(sim.lof, z, sim.lof, z + 5.0, sim.meta.OncotreeLineage)
    events = h1.enumerate_events(sim.fmap, loss, kind, sim.annot, sim.fold)
    A = sim.fold.index[sim.fold == "A"]
    E = sim.E.loc[A]
    ev = next(e for e in events if e.gene == "ARID1A")
    tab = pd.concat([sim.contexts, sim.meta.drop(columns=["TP53_altered", "WGD"])], axis=1).loc[A]
    cov = {ev.context: cv.assemble(ev.covariates, tab)}
    return sim, ev, E, loss.loc[A], cov


def test_export_matches_engine_and_classes_breadth(setup):
    sim, ev, E, loss_A, cov = setup
    units = [f"{ev.id}|{a}" for a in ev.axes]
    tests = pd.DataFrame(dict(context=ev.context, unit=units, evaluable=True, discovery=False, q_adj=np.nan))
    full = exp.build_profiles(E, [ev], loss_A, cov, tests)
    c = h1.component_scores(E, ev, loss_A[ev.gene])
    ctx = loss_A[ev.gene].rename(ev.context)
    # N axis rows == the engine's own member_profile (identical fit, alignment, baselines)
    N = full[full.axis == "N"].set_index("member")
    ref = h1.member_profile(E, list(ev.N), ctx, cov[ev.context], c.b)
    for col, rcol in (("beta", "beta"), ("se", "se"), ("ci95_lo", "ci95_lo"), ("ci95_hi", "ci95_hi"),
                      ("baseline", "baseline")):
        assert np.allclose(N[col], ref[rcol].reindex(N.index), equal_nan=True), col
    Np = N[N.in_aligned_set]
    assert set(Np.index) == set(c.Nprime)
    assert np.allclose(Np.aligned_beta, np.sign(Np.baseline) * Np.beta)          # ALIGN = +1
    assert (Np.aligned_ci95_lo <= Np.aligned_beta).all() and (Np.aligned_beta <= Np.aligned_ci95_hi).all()
    assert (N.n > 12).all() and (N.n_loss >= 12).all() and (N.n_loss < N.n).all()
    # breadth from the export == classify's breadth (the assertion used on the real runs)
    res = pd.DataFrame(dict(context=ev.context, unit=units, evaluable=True, discovery=False, q_adj=np.nan))
    cls = h1.classify([ev], {ev.id: c}, res, E, loss_A, cov, sim.fmap, sim.annot, E.clip(-3, 3).mean(), None, K=30,
                      max_attempts=3000)
    chk = exp.check_breadth(full, cls)
    assert chk["committed_events"] == 1 and chk["committed_not_evaluable"] == []
    bad = cls.copy()
    bad.loc[:, "N_frac_aligned_pos"] += 0.2
    with pytest.raises(AssertionError):
        exp.check_breadth(full, bad)


def test_unscreenable_axis_not_exported_but_breadth_checked(setup):
    sim, ev, E, loss_A, cov = setup
    c = h1.component_scores(E, ev, loss_A[ev.gene])
    units = [f"{ev.id}|{a}" for a in ev.axes]
    tests = pd.DataFrame(dict(context=ev.context, unit=units, evaluable=False))
    full = exp.build_profiles(E, [ev], loss_A, cov, tests)
    assert not full.counted_evaluable.any()
    pos = float((full[(full.axis == "N") & full.in_aligned_set].aligned_beta > 0).mean())
    classes = pd.DataFrame(dict(event=[ev.id], N_frac_aligned_pos=[pos], N_frac_aligned_neg=[1 - pos]))
    chk = exp.check_breadth(full, classes)        # pre-fix-style classes.tsv: breadth present, counted test not
    assert chk["evaluable_exported"] == 0 and chk["committed_not_evaluable"] == [ev.id]
    del c
