"""Context-table building blocks on toy data with known answers (three-valued logic, LoY antimode, segment merging,
fusion pairing), plus sanity checks on the built table when it exists (no dependency data involved)."""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from part2 import contexts as C, paths  # noqa: E402

nan = np.nan


def test_kleene_logic():
    a = np.array([1, 1, 1, 0, 0, 0, nan, nan, nan])
    b = np.array([1, 0, nan, 1, 0, nan, 1, 0, nan])
    np.testing.assert_array_equal(C.k_or(a, b), [1, 1, 1, 1, 0, nan, 1, nan, nan])
    np.testing.assert_array_equal(C.k_and(a, b), [1, 0, nan, 0, 0, 0, nan, 0, nan])
    np.testing.assert_array_equal(C.cmp(np.array([0.1, 0.6, nan]), "<", 0.5), [1, 0, nan])


def test_lof_rule_toy():
    """damaging==2 | (damaging>=1 & a < P-0.75) | homdel | (male & X & damaging>=1), three-valued."""
    D = np.array([2, 1, 1, 0, 0, 1, 1, 0])
    a = np.array([2, 1.0, 2.0, 0.2, 2, nan, 2.0, nan])
    P = np.array([2, 2, 2, 2, 2, 2, 2, nan])
    hd = np.array([0, 0, 0, 1, 0, nan, 0, 0])
    maleX = np.array([0, 0, 0, 0, 0, 0, 1, 0])
    lof = C.k_or(C.cmp(D, "==", 2), C.k_and(C.cmp(D, ">=", 1), C.cmp(a - (P - 0.75), "<", 0.0)), hd,
                 C.k_and(maleX, C.cmp(D, ">=", 1)))
    np.testing.assert_array_equal(lof, [1, 1, 0, 1, 0, nan, 1, 0])


def test_loy_threshold_bimodal():
    rng = np.random.default_rng(0)
    v = np.r_[rng.normal(2, 1, 300), rng.normal(35, 5, 400)]
    t = C.loy_threshold(v)
    assert 8 < t["threshold"] < 25 and abs(t["low_mode"] - 2) < 2 and abs(t["high_mode"] - 35) < 3
    with pytest.raises(ValueError):
        C.loy_threshold(rng.normal(0, 1, 500))


def test_merge_segments_toy():
    models = [f"m{i}" for i in range(20)]
    cols = [1, 2, 3, 4, 5]
    E = pd.DataFrame(0.0, index=models, columns=cols)
    E.loc[models[:10], [1, 2, 3]] = 1.0          # genes 1-3 identical -> one segment
    E.loc[models[:3], 4] = 1.0                   # gene 4 differs (Jaccard 0.3) -> own segment
    E.loc[models[0], 2] = nan                    # NaN ignored in Jaccard (complete pairs)
    ann = pd.DataFrame(dict(symbol=["A", "B", "C", "D", "E"], chrom=["9"] * 4 + ["10"], pos=[10, 20, 30, 40, 5],
                            location=["9p21.3"] * 3 + ["9p21.2", "10q1"], band=["9p21"] * 4 + ["10q1"]), index=cols)
    tab, mat = C.merge_segments(E, ann, "homdel", jaccard=0.9)
    assert list(tab.n_genes) == [3, 1]
    assert tab.segment.iloc[0] == "9p21.3 homdel A-C" and tab.segment.iloc[1] == "9p21.2 homdel D"
    assert mat.iloc[:, 0].sum() == 10 and mat.iloc[:, 1].sum() == 3


def test_fusion_list_is_preregistered():
    assert C.FUSIONS["KMT2A-r"] == [("KMT2A", p) for p in ("MLLT3", "AFF1", "MLLT1", "MLLT10", "AFDN", "ELL")]
    assert set(C.FUSION_MIN_U) <= set(C.FUSIONS) and min(C.FUSION_MIN_U.values()) >= 1


TABLE = paths.CONTEXT_CACHE / "context_table.parquet"


@pytest.mark.skipif(not TABLE.exists(), reason="context table not built")
def test_built_table_sanity():
    T = pd.read_parquet(TABLE)
    log = json.loads((paths.CONTEXT_RESULTS / "build_log.json").read_text())
    assert len(T) == 1345 and T.index.is_unique
    for name, mn in C.FUSION_MIN_U.items():
        assert T[f"fusion:{name}"].sum() >= mn
    # MSI-H exclusion: primary non-MMR LoF contexts are NaN in every MSI-H line, MMR genes are not
    msi = T.MSI_H == 1
    assert T.loc[msi, "LoF:ARID1A"].isna().all() and T.loc[msi, "LoF:MLH1"].notna().any()
    # LoY only in males; female_or_LoY consistent
    assert T.loc[T.sex_female == 1, "LoY"].isna().all()
    assert (T.female_or_LoY[T.sex_female == 1] == 1).all()
    # CN scale is ploidy-relative (housekeeping median r ~ 1 in WGD and non-WGD lines)
    cs = log["cn_scale_check"]
    assert 0.9 < cs["median_r_housekeeping_WGD"] < 1.15 and 0.9 < cs["median_r_housekeeping_nonWGD"] < 1.15
    # named contexts carry the expected rough sizes (round-2 biology review, U)
    assert 70 <= T.MSI_H.sum() <= 90 and 150 <= T["9p21_loss"].sum() <= 230
    assert T["melanoma"].sum() == 69
