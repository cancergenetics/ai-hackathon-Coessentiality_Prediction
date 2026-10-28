"""Graph construction: mutual top-k, threshold, capped weights, edge-count matching, degree-preserving rewire,
locus flag."""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from part2 import graph, locus  # noqa: E402


def _block():
    # 0-1-2 strong triangle, 3 points to 0 but 0 does not point back with k=2
    B = np.full((5, 5), 0.0)
    for a, b, v in [(0, 1, 9), (0, 2, 8), (1, 2, 7), (0, 3, 6.5), (3, 4, 4.5), (2, 4, 1)]:
        B[a, b] = B[b, a] = v
    np.fill_diagonal(B, np.nan)
    return B


def test_mutual_knn_and_threshold():
    g = graph.build_graph(_block(), k=2, z0=4.0)
    assert sorted(zip(g.r.tolist(), g.c.tolist())) == [(0, 1), (0, 2), (1, 2), (3, 4)]
    assert np.allclose(g.w, graph.cap_weight(g.v, 4.0))
    assert g.w.max() <= 4.0 and g.w.min() >= 1.0
    g5 = graph.build_graph(_block(), k=2, z0=5.0)
    assert (3, 4) not in set(zip(g5.r.tolist(), g5.c.tolist()))


def test_matched_graph_count_and_weights():
    rng = np.random.default_rng(0)
    A = rng.standard_normal((60, 60))
    B = (A + A.T) / 2
    np.fill_diagonal(B, np.nan)
    ref = rng.uniform(1, 4, 500)
    g, thr, ok = graph.matched_graph(B, 5, 40, ref)
    assert g.m == 40 and ok and (g.v >= thr).all()
    assert g.w.min() >= ref.min() and g.w.max() <= ref.max()
    assert np.corrcoef(np.argsort(np.argsort(g.v)), np.argsort(np.argsort(g.w)))[0, 1] > 0.99
    g2, _, ok2 = graph.matched_graph(B, 2, 10 ** 6, ref)
    assert not ok2 and g2.m > 0


def test_rewire_preserves_degree_and_simplicity():
    rng = np.random.default_rng(1)
    n = 80
    keys = np.unique(rng.integers(0, n, 600) * n + rng.integers(0, n, 600))
    r, c = keys // n, keys % n
    ok = r < c
    r, c = r[ok], c[ok]
    nr, nc = graph.degree_preserving_rewire(r, c, n, 10, 0)
    assert np.array_equal(np.bincount(np.r_[r, c], minlength=n), np.bincount(np.r_[nr, nc], minlength=n))
    assert (nr < nc).all() and len(np.unique(nr * n + nc)) == len(nr)
    assert (nr * n + nc != r * n + c).mean() > 0.5


def test_locus_flag():
    assert locus.parse_location("19q13.43") == ("19q", "19q13")
    assert locus.parse_location("Xp22.33") == ("Xp", "Xp22")
    assert locus.parse_location(np.nan) == ("", "")
    arm = np.array(["9p", "9p", "1q", "2q"])
    band = np.array(["9p21", "9p21", "1q21", "2q11"])
    assert locus.locus_flag([0, 1, 2, 3], arm, band)["locus"]          # 50% on 9p
    assert not locus.locus_flag([0, 2, 3], arm, band)["locus"]


def test_locus_window_rule():
    arm = np.array(["1p", "1q", "2q", "3q"])
    band = np.array(["1p36", "1q21", "2q11", "3q11"])
    chrom = np.array(["1", "1", "2", "3"])
    pos = np.array([120_000_000, 128_000_000, 5, 5])       # 1p/1q across the centromere, 8 Mb apart
    f = locus.locus_flag([0, 1, 2, 3], arm, band, chrom=chrom, pos=pos)
    assert f["window_share"] == 0.5 and f["arm_share"] == 0.25 and f["locus"]
    pos2 = np.array([100_000_000, 128_000_000, 5, 5])
    assert not locus.locus_flag([0, 1, 2, 3], arm, band, chrom=chrom, pos=pos2)["locus"]
    assert locus.max_in_window(np.array(["1", "1", ""]), np.array([0, 10_000_000, -1])) == 2


def test_gene_coordinates_join():
    from part2 import paths
    if not paths.GENE_COORDS.exists():
        return
    L = locus.gene_locations()
    assert (L.pos >= 0).mean() > 0.99
    both = (L.chrom != "") & (L.arm != "")
    assert (L.chrom[both] == L.arm[both].str[:-1]).all()        # Ensembl chromosome agrees with HGNC location
