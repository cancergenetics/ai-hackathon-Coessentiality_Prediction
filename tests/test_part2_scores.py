"""Lean locked scorer: bit-exact reproduction of Part 1's stored zavg_BEF matrix on all models, null moments on all
models equal Part 1's manifest, GLS equal to Part 1's predict_pairs GLS (gls_signed), top-K helpers.
The two real-data tests use the cached all-model blocks (scripts/31_score_blocks.py --set all --S-from all, ~2 min
each if absent). No hu.MAP file is read."""
import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from part2 import data, genefilter, paths, scores  # noqa: E402
from part2.part1 import Config, CoessentialityModel, p1_baselines  # noqa: E402

REAL = paths.GENE_EFFECT_CSV.exists() and paths.ZAVG_NPY.exists()


@pytest.fixture(scope="module")
def all_models_S():
    ge = data.load_gene_effect()
    ids = list(ge.models)
    return ids, genefilter.dependency_bearing_for_models(ids)


# ------------------------------------------------------------------------------------------ helpers
def test_select_topk_ties_and_nan():
    S = np.array([[np.nan, 1.0, 1.0, 0.5, -np.inf],
                  [2.0, np.nan, 2.0, 2.0, 2.0]])
    I, V = scores.select_topk(S, 3)
    assert I[0].tolist() == [1, 2, 3] and V[0].tolist() == [1.0, 1.0, 0.5]
    assert I[1].tolist() == [0, 2, 3]                       # ties by ascending index
    I, V = scores.select_topk(S[:1], 6)
    assert I[0].tolist() == [1, 2, 3, -1, -1, -1] and np.isnan(V[0, 3:]).all()


def test_mirror_upper():
    B = np.arange(16, dtype=np.float32).reshape(4, 4)
    scores.mirror_upper(B)
    assert np.allclose(B, B.T, equal_nan=True) and np.isnan(np.diag(B)).all() and B[2, 0] == 2.0


def _toy_X(seed=0, L=80, G=400, nan_genes=()):
    rng = np.random.default_rng(seed)
    f = rng.standard_normal((L, 5))
    X = rng.standard_normal((L, G)) * 0.3 + f @ rng.standard_normal((5, G)) * 0.2
    X[:, :10] += np.outer(rng.standard_normal(L), np.ones(10))      # planted module
    for g in nan_genes:
        X[rng.random(L) < 0.3, g] = np.nan
    return X.astype(np.float32)


def test_zavg_block_matches_pair_api_on_toy():
    X = _toy_X()
    cfg = Config(k=3, d=16, q=5.0, n_null=5000)
    S = np.zeros(X.shape[1], bool)
    S[::3] = True
    B, (TI, TV), _, info = scores.zavg_block(X, np.arange(X.shape[1]), S, K_all=7, config=cfg)
    m = CoessentialityModel(X, cfg)
    s = np.flatnonzero(S)
    ia, ib = np.triu_indices(len(s), 1)
    z, _ = m.zavg(s[ia], s[ib])
    assert np.abs(B[ia, ib] - z).max() < 1e-4                       # float32 storage vs pair API (float64)
    assert np.array_equal(B, B.T, equal_nan=True)
    assert set(info["null_moments"]) == {"B", "E", "F"}
    assert (TI[:, 0] >= 0).all() and (TI != s[:, None]).all()      # all-gene top-K excludes self


def test_zavg_rect_equals_full_block_and_leaves_SxS_unchanged():
    """Attachment rectangle (non-S x S) equals the all-gene block restricted to (non-S, S); S x S block and top-K
    are identical with and without rect."""
    X = _toy_X(2)
    cfg = Config(k=3, d=16, q=5.0, n_null=5000)
    G = X.shape[1]
    S = np.zeros(G, bool)
    S[::3] = True
    gidx = np.arange(G)
    B0, (TI0, TV0), _, _ = scores.zavg_block(X, gidx, S, K_all=7, config=cfg)
    B1, (TI1, TV1), _, info = scores.zavg_block(X, gidx, S, K_all=7, config=cfg, rect=True)
    assert np.array_equal(B0, B1, equal_nan=True) and np.array_equal(TI0, TI1)
    Rv, rpos = info["_rect"]
    Bfull = scores.zavg_block(X, gidx, np.ones(G, bool), K_all=0, config=cfg)[0]
    assert np.array_equal(rpos, np.flatnonzero(~S))
    assert np.array_equal(Rv, Bfull[np.ix_(~S, S)], equal_nan=True)


def test_pearson_block_matches_part1_pearson_pairs():
    X = _toy_X(3, nan_genes=(1, 2, 7)).astype(np.float64)
    S = np.zeros(X.shape[1], bool)
    S[:40] = True
    B = scores.pearson_block(X, np.arange(X.shape[1]), S)[0]
    ia, ib = np.triu_indices(40, 1)
    ref = p1_baselines.pearson_pairs(X, ia, ib)["signed"]
    assert np.nanmax(np.abs(B[ia, ib] - ref)) < 1e-5


def test_gls_block_matches_part1_gls_pairs_with_nan_patterns():
    X = _toy_X(1, nan_genes=(3, 4, 5, 50, 51)).astype(np.float64)
    X[:, 200:] = np.nan_to_num(X[:, 200:])
    S = np.zeros(X.shape[1], bool)
    S[:120] = True
    B, _, _, info = scores.gls_block(X, np.arange(X.shape[1]), S)
    assert info["n_patterns"] >= 3
    ia, ib = np.triu_indices(120, 1)
    ref = p1_baselines.gls_pairs(X, ia, ib)["signed"]
    ok = np.isfinite(ref)
    assert np.abs(B[ia, ib][ok] - ref[ok]).max() < 1e-4 * np.abs(ref[ok]).max()   # float32 storage
    assert np.allclose(B[ia, ib][ok], ref[ok], rtol=1e-5, atol=1e-4)


# ------------------------------------------------------------------------------------------ real data
@pytest.mark.skipif(not REAL, reason="DepMap / Part 1 outputs not available")
def test_zavg_all_models_bit_exact_vs_part1_matrix(all_models_S):
    ids, S = all_models_S
    blk = scores.score_block(ids, "zavg_BEF", S)
    Z = np.load(paths.ZAVG_NPY, mmap_mode="r")
    s = blk.s_index
    rows = np.sort(np.random.default_rng(0).choice(len(s), 200, replace=False))
    ref = np.asarray(Z[s[rows]])
    got = blk.vals[rows]
    same = (ref[:, s] == got) | (np.isnan(ref[:, s]) & np.isnan(got))
    assert same.all(), f"max |diff| {np.nanmax(np.abs(ref[:, s].astype(float) - got))}"
    # all-gene top-50 lists over all protein genes equal those of the stored rows
    full = ref.astype(np.float64)
    full[np.arange(len(rows)), s[rows]] = -np.inf
    I, V = scores.select_topk(full, 50)
    assert np.array_equal(I, blk.topk_all[0][rows]) and np.array_equal(V, blk.topk_all[1][rows])


@pytest.mark.skipif(not REAL, reason="DepMap / Part 1 outputs not available")
def test_null_moments_all_models_equal_part1_manifest(all_models_S):
    ids, S = all_models_S
    blk = scores.score_block(ids, "zavg_BEF", S)
    man = json.loads((paths.PART1_ALL_PAIRS / "manifest.json").read_text())["method"]["null_moments"]
    for c in "BEF":
        assert abs(blk.meta["null_moments"][c]["mean"] - man[c]["mean"]) <= 1e-12
        assert abs(blk.meta["null_moments"][c]["sd"] - man[c]["sd"]) <= 1e-12
    assert blk.meta["config"] == {k: getattr(Config(), k) for k in Config.__dataclass_fields__}


@pytest.mark.skipif(not REAL, reason="DepMap / Part 1 outputs not available")
def test_gls_all_models_vs_part1_predict_baseline_engine(all_models_S):
    ids, S = all_models_S
    blk = scores.score_block(ids, "gls", S)
    from src import predict as p1predict                                     # Part 1 src/predict.py (locked)
    ge = data.load_gene_effect()
    genes = data.protein_genes()
    eng = p1predict.BaselineEngine(ge.X)
    rng = np.random.default_rng(1)
    a = rng.integers(0, len(blk.s_index), 3000)
    b = rng.integers(0, len(blk.s_index), 3000)
    keep = a != b
    a, b = a[keep], b[keep]
    gi = genes.depmap_index.values
    ref = eng.compute("gls", gi[blk.s_index[a]], gi[blk.s_index[b]])
    got = blk.vals[a, b].astype(np.float64)
    assert np.isfinite(ref).all()
    rel = np.abs(got - ref) / np.maximum(1.0, np.abs(ref))
    assert rel.max() < 1e-5, rel.max()                                        # float32 storage
