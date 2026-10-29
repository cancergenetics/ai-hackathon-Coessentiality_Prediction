"""Registered calibrations of part2.calib (PI ruling R4) on toy / synthetic data."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from part2 import calib, synth  # noqa: E402
from part2.cmap import FrozenMap, symbol  # noqa: E402
from part2.cscore import Annotation  # noqa: E402


def _sim():
    sim = synth.simulate(2)
    A = sim.fold.index[sim.fold == "A"]
    return sim, sim.E.loc[A]


def test_matched_random_sets_rule_and_determinism():
    sim, E = _sim()
    ess = calib.clip_essentiality(E)
    ps = 10 * E.clip(-3, 3)[sorted(sim.fmap.all_genes(1))].corr()
    units = sim.fmap.coarse()[:6]
    s1, st = calib.matched_random_sets(E, sim.fmap, units, ess, ps, 5, 4600)
    s2, _ = calib.matched_random_sets(E, sim.fmap, units, ess, ps, 5, 4600)
    assert s1 == s2 and len(st) == 6
    for k, genes in s1.items():
        u = k.split("#")[0]
        row = st.set_index("unit").loc[u]
        assert len(genes) == row["size"] and not set(genes) & set(sim.fmap.members[u])
        tol = 0.1 if row.status == "ok" else 0.2
        assert abs(ess[genes].mean() - row.target_ess) <= tol + 1e-12
    assert set(st.status) <= {"ok", "tolerances doubled", "unavailable"}


def test_dispersed_sets_and_never_expressed_pool():
    an = Annotation(chrom={f"g{i}": str(1 + i % 3) for i in range(60)},
                    pos={f"g{i}": (i // 3) * 4_000_000 for i in range(60)})
    keys = [f"g{i} ({i})" for i in range(60)]
    ex = pd.DataFrame(0.2, index=[f"m{j}" for j in range(5)], columns=keys)
    ex.iloc[0, :10] = 3.0                                                 # g0..g9 expressed in one model
    pool = calib.never_expressed_pool(ex, keys, an)
    assert pool == sorted(keys[10:])
    sets = calib.dispersed_sets(pool, [4, 6], an, np.random.default_rng(0))
    for genes in sets.values():
        assert all(not an.near(symbol(a), symbol(b)) for i, a in enumerate(genes) for b in genes[i + 1:])
    assert [len(v) for v in sets.values()] == [4, 6]


def test_segment_pseudo_contexts_filters_event_locus_and_size():
    an = Annotation(chrom={"MTAP": "9", "X1": "9", "Y1": "3"}, pos={"MTAP": 21_000_000, "X1": 25_000_000,
                                                                     "Y1": 50_000_000})
    tab = pd.DataFrame(dict(segment=["s_same", "s_ok", "s_small", "s_amp"], event=["homdel"] * 3 + ["highamp"],
                            chrom=["9", "3", "3", "3"], start_pos=[22_000_000, 49_000_000, 49_000_000, 49_000_000],
                            end_pos=[23_000_000, 51_000_000, 51_000_000, 51_000_000]))
    n = 200
    m = pd.DataFrame({"s_same": [1.0] * 40 + [0.0] * 160, "s_ok": [1.0] * 45 + [0.0] * 155,
                      "s_small": [1.0] * 20 + [0.0] * 180, "s_amp": [1.0] * 40 + [0.0] * 160}, index=range(n))
    pcs, near = calib.segment_pseudo_contexts(tab, m, "MTAP", "homdel", 40, an)
    assert list(pcs.columns) == ["s_ok"] and "Y1" in near["s_ok"] and "MTAP" not in near["s_ok"]


def test_random_splits_and_test_columns():
    rng = np.random.default_rng(1)
    n = 400
    keys = [f"g{i} ({i})" for i in range(10)]
    E = pd.DataFrame(-0.5 + np.outer(rng.normal(size=n), np.ones(10)) * 0.4 + rng.normal(0, 0.2, (n, 10)),
                     columns=keys, index=[f"m{i}" for i in range(n)])
    fm = FrozenMap(members={"P": keys, "P.S1": keys[:4], "P.S2": keys[4:]}, level={"P": 1, "P.S1": 2, "P.S2": 2},
                   parent={"P": "", "P.S1": "P", "P.S2": "P"})
    Y = calib.random_split_scores(E, fm, "P.S1", None, Annotation(), 5, np.random.default_rng(4800))
    assert Y.shape == (n, 5) and np.allclose(Y.mean(), 0, atol=1e-10)    # OLS residuals have mean 0
    ctx = pd.Series((rng.random(n) < 0.02).astype(float), index=E.index)  # < 12 context models
    f = calib.test_columns(Y, ctx, None)
    assert f.p.isna().all()
    ctx2 = pd.Series((rng.random(n) < 0.3).astype(float), index=E.index)
    assert calib.test_columns(Y, ctx2, None).p.notna().all()
