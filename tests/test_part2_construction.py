"""FZ-C construction decisions (results/part2/construction/DECISIONS.md): nested placement (D2), containment support
(D1), attachment (D3), GLS rank-mapped cohesion scale (D4), parallel stage 1 == serial, A-only run sets."""
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))

from part2 import attach, construct  # noqa: E402
from part2 import consensus as C  # noqa: E402
from test_part2_consensus import planted_runs  # noqa: E402

NEST = C.ClusterConfig(seeds=(0, 1, 2), support_mode="containment", placement="nested")


def _check_hierarchy(cm):
    l1 = np.concatenate([c["members"] for c in cm.level(1)]) if cm.level(1) else np.array([])
    assert len(l1) == len(np.unique(l1))                                 # one level-1 cluster per gene
    for p in cm.level(1):
        kids = [c for c in cm.level(2) if c["parent"] == p["id"]]
        if kids:
            km = np.concatenate([c["members"] for c in kids])
            assert len(km) == len(np.unique(km))                         # one sub-module per gene per parent
            assert set(km.tolist()) <= set(np.asarray(p["members"]).tolist())


def test_nested_placement_recovers_planted_and_keeps_hierarchy():
    runs, groups, halves, nest = planted_runs()
    cm = C.build_map(runs, NEST)
    found = [set(c["members"].tolist()) for c in cm.clusters]
    for g in groups + halves:
        assert max(len(set(g) & f) / len(set(g) | f) for f in found) >= 0.75
    _check_hierarchy(cm)
    assert all(3 <= c["size"] <= NEST.max_size for c in cm.clusters)


def test_nested_placement_null_graph_empty():
    rng = np.random.default_rng(3)
    from part2.graph import Graph
    n, runs = 150, []
    for _ in range(6):
        keys = np.unique(rng.integers(0, n, 300) * n + rng.integers(0, n, 300))
        r, c = keys // n, keys % n
        ok = r < c
        w = rng.uniform(1, 2, ok.sum())
        runs.append(Graph(n, r[ok], c[ok], w + 2, w))
    assert len(C.build_map(runs, NEST).clusters) == 0


def test_parent_absorbs_more_persistent_children():
    """A persistent coarse candidate that contains accepted clusters becomes their parent (they move to level 2)."""
    a = dict(members=np.arange(0, 6), n_adjacent=8, support=1.0, j=5, j_start=3, j_end=10, gamma=0.3, stability=1)
    b = dict(members=np.arange(6, 12), n_adjacent=7, support=1.0, j=5, j_start=4, j_end=10, gamma=0.3, stability=1)
    p = dict(members=np.arange(0, 13), n_adjacent=3, support=0.9, j=0, j_start=0, j_end=2, gamma=0.05, stability=1)
    part = dict(members=np.r_[np.arange(10, 12), np.arange(20, 24)], n_adjacent=3, support=0.8, j=1, j_start=1,
                j_end=3, gamma=0.07, stability=1)
    counts = {k: 0 for k in ("oversize", "cohesion_fail", "reclustered", "overlap_dropped", "duplicate",
                              "nested_child", "overlap_conflict", "trimmed", "parent_rule_fail", "parent_absorbed")}
    acc, kids = C._place_nested([p, a, b, part], [], 40, NEST, None, 0.0, counts, None, None)
    sizes = sorted(len(x["members"]) for x in acc)
    assert sizes == [4, 13] and counts["parent_absorbed"] == 2 and counts["trimmed"] == 1
    big = next(i for i, x in enumerate(acc) if len(x["members"]) == 13)
    assert sorted(len(k["members"]) for k in kids[big]) == [6, 6]
    # a candidate mostly overlapping the more persistent cluster is dropped, not trimmed
    clash = dict(members=np.r_[np.arange(0, 4), 30, 31], n_adjacent=3, support=0.7, j=1, j_start=1, j_end=3,
                 gamma=0.07, stability=1)
    counts2 = dict(counts)
    acc2, _ = C._place_nested([a, clash], [], 40, NEST, None, 0.0, counts2, None, None)
    assert len(acc2) == 1 and counts2["overlap_conflict"] == 1


def test_parallel_levels_equal_serial():
    runs, *_ = planted_runs(seed=2)
    cfg = C.ClusterConfig(seeds=(0, 1), gammas=C.gamma_grid(n=4))
    a = C.build_levels(runs, cfg)
    b = construct.build_levels_parallel(runs, cfg, n_proc=2)
    for x, y in zip(a, b):
        assert np.array_equal(x.runs, y.runs) and np.array_equal(x.fr, y.fr) and np.allclose(x.ff, y.ff)


def test_rank_map_is_monotone_onto_reference_distribution():
    rng = np.random.default_rng(0)
    n = 60
    A = rng.normal(size=(n, n))
    A = (A + A.T) / 2
    R = rng.gamma(2.0, size=(n, n))
    R = (R + R.T) / 2
    A[3, 5] = A[5, 3] = np.nan
    M = construct.rank_map(A.astype(np.float32), R.astype(np.float32))
    iu = np.triu_indices(n, 1)
    a, m = A[iu], M[iu]
    ok = np.isfinite(a)
    assert np.isnan(m[~ok]).all() and np.array_equal(M, M.T, equal_nan=True)
    o = np.argsort(a[ok])
    assert (np.diff(m[ok][o]) >= -1e-6).all()
    assert abs(np.median(m[ok]) - np.median(R[iu])) < 0.05


def test_bh_partial_equals_padded_bh():
    rng = np.random.default_rng(1)
    p = rng.uniform(size=50) ** 3
    full = np.r_[p, np.ones(150)]
    o = np.argsort(full, kind="stable")
    q = full[o] * len(full) / np.arange(1, len(full) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    ref = np.empty_like(q)
    ref[o] = np.minimum(q, 1)
    assert np.allclose(attach.bh_partial(p, 200), ref[:50])


def test_attach_planted_member_and_null():
    rng = np.random.default_rng(2)
    nS, nR = 300, 200
    rect = rng.normal(0, 1, size=(nR, nS)).astype(np.float32)
    cl = [np.arange(0, 10), np.arange(10, 30)]
    rect[0, :10] = 8.0                                     # gene 0 belongs to cluster 0
    rect[1, 10:22] = 8.0                                   # gene 1: 12/20 members of cluster 1
    rect[2, 10:18] = 8.0                                   # gene 2: 8/20 (< 50%) -> not attached
    rect[3, :] = 8.0                                       # gene 3: hub (high to everything) -> not significant
    rect[4, :] = np.nan                                    # untested
    df, s = attach.attach(rect, cl, z0=5.0)
    got = set(zip(df.gene.tolist(), df.cluster.tolist()))
    assert got == {(0, 0), (1, 1)}
    assert s["n_genes_tested"] == nR - 1 and s["n_tests"] == (nR - 1) * 2


def test_runsets_are_subsets_of_A_and_patient_blocked():
    pytest.importorskip("pandas")
    from part2 import data, paths, runsets
    if not paths.MODEL_CSV.exists():
        pytest.skip("no DepMap Model.csv")
    A = set(runsets.A())
    assert len(A) == 1007
    pid = data.load_models().PatientID
    for s in (1000, 1019):
        b = runsets.models(f"boot{s}")
        assert set(b) <= A and len(b) == round(0.8 * len(A))
    for s in runsets.HALF_SEEDS:
        h1, h2 = runsets.half(s)
        assert not set(h1) & set(h2) and set(h1) | set(h2) == A
        assert not set(pid.reindex(h1).dropna()) & set(pid.reindex(h2).dropna())
    for bad in ("all", "U", "L", "boot999", "half1999a"):
        with pytest.raises(ValueError):
            runsets.models(bad)


DET = r"""
import sys, json
sys.path.insert(0, %r)
sys.path.insert(0, %r)
from test_part2_consensus import planted_runs
from part2 import consensus as C
runs, *_ = planted_runs(seed=5)
cm = C.build_map(runs, C.ClusterConfig(seeds=(0, 1, 2), support_mode="containment", placement="nested"))
print(json.dumps([[c["level"], c["parent"], round(c["gamma"], 6), round(c["support"], 6), c["members"].tolist()]
                  for c in cm.clusters]))
"""


def test_nested_determinism_in_subprocess():
    def run(hs):
        env = dict(os.environ, PYTHONHASHSEED=str(hs))
        code = DET % (str(ROOT / "src"), str(ROOT / "tests"))
        out = subprocess.run([sys.executable, "-c", code], env=env, capture_output=True, text=True, check=True)
        return json.loads(out.stdout.strip().splitlines()[-1])
    a = run(0)
    assert a and run(0) == a and run(12345) == a
