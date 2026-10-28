"""Toy simulation through the same pipeline (locked CoessentialityModel with small settings): planted complexes are
recovered, a null without complexes gives no persistent clusters, arm-block co-dependency yields only locus-flagged
clusters and earns no split-half credit; dyad maps score 0 on the selection criterion."""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from part2 import consensus as C  # noqa: E402
from part2 import genefilter, locus, pipeline, selection, toy  # noqa: E402
from part2.graph import build_graph  # noqa: E402
from part2.part1 import Config  # noqa: E402
from part2.scores import zavg_block  # noqa: E402

CFG = Config(k=10, d=64)
CC = C.ClusterConfig(seeds=(0, 1, 2))
Z0, K = 5.0, 10


@pytest.fixture(scope="module", params=["planted", "null", "armblock"])
def scen(request):
    T = toy.simulate(request.param, seed=0)
    S, st = genefilter.dependency_bearing(T.X)
    full, boots = pipeline.matrix_runs(T.X, S, 4, 0.8, 1000, CFG)
    graphs = [build_graph(b, K, Z0) for b in boots]
    cm = C.build_map(graphs, CC, full, Z0, levels=C.build_levels(graphs, CC))
    return request.param, T, S, st, full, cm


def test_scenario(scen):
    name, T, S, st, full, cm = scen
    sidx = np.flatnonzero(S)
    arm, band = T.arm[sidx], T.band[sidx]
    flags = [locus.locus_flag(c["members"], arm, band)["locus"] for c in cm.clusters]
    assert ((T.gene_class == "never") & S).sum() <= 0.02 * (T.gene_class == "never").sum()   # filter works
    if name == "planted":
        loc = np.full(len(S), -1)
        loc[sidx] = np.arange(len(sidx))
        truth = [loc[c][loc[c] >= 0] for c in T.complexes]
        pg = np.unique(np.concatenate(truth))
        assert toy.adjusted_rand(truth, cm.labels(len(sidx), 1), pg) >= 0.8
        f1 = toy.best_match_f1(truth, [c["members"] for c in cm.clusters])
        assert np.median(f1) >= 0.9 and (f1[np.array([len(t) for t in truth]) >= 5] >= 0.75).all()
        assert sum(flags) == 0
    elif name == "null":
        assert len(cm.clusters) == 0
    else:
        assert all(flags), "arm-block co-dependency must only give locus-flagged clusters"


def test_split_half_does_not_reward_arm_blocks_or_dyads():
    T = toy.simulate("armblock", seed=0)
    S, st = genefilter.dependency_bearing(T.X)
    sidx = np.flatnonzero(S)
    perm = np.random.default_rng(3000).permutation(len(T.X))
    halves = [tuple(zavg_block(np.ascontiguousarray(T.X[np.sort(h)]), np.arange(T.X.shape[1]), S, 0, CFG)[0]
                    for h in (perm[:150], perm[150:]))]
    sel = genefilter.selective_mask(st["frac_below"][sidx])
    r = selection.split_half_statistic(halves, Z0, K, sel, CC, arm=T.arm[sidx], band=T.band[sidx], rewire_seed=None)
    assert r["stat"] == 0.0
    dy = C.ClusterMap(clusters=[dict(members=np.array([2 * i, 2 * i + 1]), id=i, level=1) for i in range(50)])
    assert selection.recovered_coverage(dy, dy, np.ones(len(sidx), bool))[0] == 0.0
