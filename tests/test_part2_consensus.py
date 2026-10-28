"""Consensus + persistence clustering: co-clustering frequencies, persistence on planted cliques, size / cohesion
rules, nesting, null graphs, determinism in fresh subprocesses (PYTHONHASHSEED=0)."""
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from part2 import consensus as C  # noqa: E402
from part2.graph import Graph  # noqa: E402

CFG = C.ClusterConfig(seeds=(0, 1, 2))


def test_coclustering_hand():
    M = np.array([[0, 0, 1, 1, 2],
                  [0, 0, 0, 1, 1]])
    r, c, f = C.coclustering(M)
    d = {(int(a), int(b)): float(x) for a, b, x in zip(r, c, f)}
    assert d == {(0, 1): 1.0, (2, 3): 0.5, (0, 2): 0.5, (1, 2): 0.5, (3, 4): 0.5}


def planted_runs(n_runs=6, seed=0, sizes=(4, 6, 10, 15), sub=(8, 7), n_noise=60, p_noise=0.02, p_drop=0.1):
    """Cliques (weights ~3-4) + one nested complex (two dense halves, weaker cross edges) + sparse noise edges."""
    rng = np.random.default_rng(seed)
    groups, start = [], 0
    for s in sizes:
        groups.append(np.arange(start, start + s))
        start += s
    nest = np.arange(start, start + sum(sub))
    halves = [nest[: sub[0]], nest[sub[0]:]]
    start += sum(sub)
    n = start + n_noise
    runs = []
    for t in range(n_runs):
        e = {}
        for g in groups + halves:
            for i in range(len(g)):
                for j in range(i + 1, len(g)):
                    if rng.random() > p_drop:
                        e[(g[i], g[j])] = rng.uniform(3, 4)
        for i in halves[0]:
            for j in halves[1]:
                if rng.random() < 0.35:
                    e[(i, j)] = rng.uniform(1, 2)
        for _ in range(int(p_noise * n * n / 2)):
            i, j = sorted(rng.choice(n, 2, replace=False))
            e.setdefault((int(i), int(j)), rng.uniform(1, 1.5))
        k = sorted(e)
        r = np.array([a for a, _ in k])
        c = np.array([b for _, b in k])
        w = np.array([e[x] for x in k])
        runs.append(Graph(n, r, c, w + 2, w))            # v = w + z0 - 1 with z0 = 3
    return runs, groups, halves, nest


def test_persistence_recovers_planted_cliques_and_nesting():
    runs, groups, halves, nest = planted_runs()
    cm = C.build_map(runs, CFG)
    found = [set(c["members"].tolist()) for c in cm.clusters]
    for g in groups:
        assert max(len(set(g) & f) / len(set(g) | f) for f in found) >= 0.75
    l1 = [set(c["members"].tolist()) for c in cm.level(1)]
    l2 = [c for c in cm.level(2)]
    # the nested complex appears whole at level 1 and its halves as level-2 sub-modules, or its halves at level 1
    whole = max(len(set(nest) & f) / len(set(nest) | f) for f in l1)
    halves_any = [max(len(set(h) & f) / len(set(h) | f) for f in found) for h in halves]
    assert whole >= 0.75 or min(halves_any) >= 0.75
    for c in l2:
        parent = next(p for p in cm.level(1) if p["id"] == c["parent"])
        assert set(c["members"]) <= set(parent["members"])
    labels = np.concatenate([c["members"] for c in cm.level(1)])
    assert len(labels) == len(np.unique(labels))                       # one cluster per gene per level
    assert all(3 <= c["size"] <= CFG.max_size for c in cm.clusters)


def test_null_graph_gives_no_persistent_clusters():
    rng = np.random.default_rng(3)
    n, runs = 150, []
    for t in range(6):
        keys = np.unique(rng.integers(0, n, 300) * n + rng.integers(0, n, 300))
        r, c = keys // n, keys % n
        ok = r < c
        r, c = r[ok], c[ok]
        w = rng.uniform(1, 2, len(r))
        runs.append(Graph(n, r, c, w + 2, w))
    cm = C.build_map(runs, CFG)
    assert len(cm.clusters) == 0, [c["size"] for c in cm.clusters]


def test_size_and_cohesion_rules():
    runs, groups, _, _ = planted_runs(sizes=(4, 6, 10))
    n = runs[0].n
    ref = np.full((n, n), 0.0)
    for g in groups:
        ref[np.ix_(g, g)] = 10.0
    ref[np.ix_(groups[2], groups[2])] = 2.5                  # big clique fails coarse cohesion (z0 = 3)
    cm = C.build_map(runs, CFG, ref_scores=ref, z0=3.0)
    sizes = [c["size"] for c in cm.level(1)]
    assert 10 not in sizes and cm.counts["cohesion_fail"] >= 1
    cm2 = C.build_map(runs, C.ClusterConfig(seeds=(0, 1, 2), max_size=5))
    assert all(c["size"] <= 5 for c in cm2.clusters) and cm2.counts["oversize"] >= 1


DET = r"""
import sys, json
sys.path.insert(0, %r)
sys.path.insert(0, %r)
from test_part2_consensus import planted_runs
from part2 import consensus as C
runs, *_ = planted_runs(seed=5)
cm = C.build_map(runs, C.ClusterConfig(seeds=(0, 1, 2)))
print(json.dumps([[c["level"], c["parent"], round(c["gamma"], 6), round(c["support"], 6), c["members"].tolist()]
                  for c in cm.clusters]))
"""


def _run(hashseed):
    env = dict(os.environ, PYTHONHASHSEED=str(hashseed))
    code = DET % (str(ROOT / "src"), str(ROOT / "tests"))
    out = subprocess.run([sys.executable, "-c", code], env=env, capture_output=True, text=True, check=True)
    return json.loads(out.stdout.strip().splitlines()[-1])


def test_determinism_in_subprocess():
    a, b = _run(0), _run(0)
    assert a == b and len(a) > 0
    assert _run(12345) == a                                             # no hash-order dependence either
