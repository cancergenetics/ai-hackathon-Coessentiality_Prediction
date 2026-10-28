"""Tier tests and bridge endpoint on synthetic networks with known answers (no BioPlex / proteomics file is read)."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from part2 import validation as V  # noqa: E402


def synth(seed=0, n=600, p=0.01, planted=((0, 8), (8, 20))):
    rng = np.random.default_rng(seed)
    entrez = np.arange(1, n + 1) * 10
    keys = set()
    iu = np.triu_indices(n, 1)
    e = rng.random(len(iu[0])) < p
    for a, b in zip(iu[0][e], iu[1][e]):
        keys.add((int(entrez[a]), int(entrez[b])))
    for lo, hi in planted:                                      # dense planted complexes
        for a in range(lo, hi):
            for b in range(a + 1, hi):
                if rng.random() < 0.6:
                    keys.add((int(entrez[a]), int(entrez[b])))
    arr = np.array(sorted(keys))
    deg = pd.Series(np.r_[arr[:, 0], arr[:, 1]]).value_counts().sort_index()
    return entrez, keys, deg


def test_t1_detects_planted_and_is_calibrated_on_random_sets():
    entrez, keys, deg = synth()
    net = V.Network(entrez, keys, deg)
    ess = np.where(np.arange(len(entrez)) % 5 == 0, "pan", "selective")
    rng = np.random.default_rng(1)
    rnd = [rng.choice(np.arange(40, 600), 6, replace=False) for _ in range(200)]
    clusters = [np.arange(0, 8), np.arange(8, 20)] + rnd
    t1 = V.t1_test(net, clusters, ess, B=10_000, seed=7)
    assert (t1.p.iloc[:2] < 0.001).all()
    pr = t1.p.iloc[2:].values
    assert 0.0 <= (pr < 0.05).mean() <= 0.10                 # discrete, conservative null
    t2 = pd.DataFrame(dict(i=t1.i, n_meas=0, mean_r=np.nan, null_mean=np.nan, p=np.nan))
    tiers = V.assign_tiers(t1, t2)
    assert (tiers.tier.iloc[:2] == "T1").all()
    assert (tiers.tier.iloc[2:] == "T1").mean() < 0.05


def test_t2_detects_coabundance():
    rng = np.random.default_rng(3)
    n = 300
    f = rng.normal(size=(100, 1))
    Y = rng.normal(size=(100, n))
    Y[:, :6] += 2 * f
    from part2.scores import pairwise_pearson
    R = pairwise_pearson(Y)
    measured = np.ones(n, bool)
    t2 = V.t2_test(R, measured, [np.arange(6), np.arange(50, 56)], B=2000)
    assert t2.p.iloc[0] < 0.001 and t2.p.iloc[1] > 0.01


def test_bridge_bootstrap_detects_better_map_and_null_difference():
    entrez, keys, deg = synth(planted=((0, 8), (8, 20), (20, 30), (30, 34)))
    net = V.Network(entrez, keys, deg)
    good = [np.arange(0, 8), np.arange(8, 20), np.arange(20, 30), np.arange(30, 34)]
    rng = np.random.default_rng(5)
    bad = [rng.choice(np.arange(100, 600), len(c), replace=False) for c in good]
    w = V.bin_weights([len(c) for c in good], [len(c) for c in bad])
    res = V.bridge_bootstrap(V.pair_table(net, good), V.pair_table(net, bad), len(net.obs_idx), w, B=300)
    assert res["diff"] > 0.3 and res["excludes_0"] and res["ci95"][0] > 0
    res0 = V.bridge_bootstrap(V.pair_table(net, good), V.pair_table(net, good), len(net.obs_idx), w, B=100)
    assert res0["diff"] == 0 and not res0["excludes_0"]


def test_match_coverage_trims_least_stable():
    a = pd.DataFrame(dict(cluster_id=["Z1", "Z2", "Z3"], size=[10, 5, 5], stability=[0.9, 0.2, 0.5],
                          cohesion=[6, 6, 6]))
    b = pd.DataFrame(dict(cluster_id=["G1"], size=[15], stability=[0.9], cohesion=[6]))
    ta, tb, info = V.match_coverage(a, b)
    assert list(ta.cluster_id) == ["Z1", "Z3"] and info["genes_after"] == dict(a=15, b=15)


def test_require_frozen_refuses_uncommitted(tmp_path):
    d = V.paths.RESULTS / "construction" / "__nonexistent_map__"
    with pytest.raises(Exception):
        V.require_frozen([d])


def test_strata_min_size():
    deg = np.r_[np.ones(50), np.full(3, 1000), np.full(30, 8)].astype(int)
    ess = np.array(["selective"] * len(deg))
    st = V.strata(deg, ess, min_n=10)
    u, c = np.unique(st, return_counts=True)
    assert (c >= 10).all()
