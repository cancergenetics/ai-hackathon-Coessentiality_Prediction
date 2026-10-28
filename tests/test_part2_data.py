"""Data layer: universe, gene keys vs Part 1 genes.tsv, subset invariance, PatientID-blocked sampling, gene filter."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from part2 import data, genefilter, paths  # noqa: E402

pytestmark = pytest.mark.skipif(not paths.GENE_EFFECT_CSV.exists(), reason="DepMap not available")


def test_universe_excludes_engineered():
    U = data.universe_models()
    ge = data.load_gene_effect()
    assert len(ge.models) == 1361 and len(U) == 1345
    m = data.load_models().loc[U]
    assert m.EngineeredModel.isna().all() and m.CulturedResistanceDrug.isna().all()


def test_gene_keys_consistent_with_part1():
    ge = data.load_gene_effect()
    g = data.protein_genes()
    assert (ge.columns[g.depmap_index.values] == g.depmap_column.values).all()
    assert (ge.entrez[g.depmap_index.values] == g.entrez.values).all()
    from part2.part1 import p1_all_pairs
    from part2.part1 import CoessentialityModel  # noqa: F401
    nan_frac = np.isnan(ge.X).mean(0)
    tab, _ = p1_all_pairs.protein_genes(ge.columns, paths.HGNC_FILE, nan_frac, 0.10)
    assert (tab.depmap_index.values == g.depmap_index.values).all()
    assert (tab.low_coverage.values == g.low_coverage.values).all()


def test_subset_is_function_of_the_set():
    ge = data.load_gene_effect()
    ids = list(ge.models[[5, 1, 9, 3]])
    A = ge.subset(ids)
    B = ge.subset(ids[::-1])
    assert np.array_equal(A, B, equal_nan=True)
    assert data.model_set_hash(ids) == data.model_set_hash(ids[::-1])
    with pytest.raises(KeyError):
        ge.subset(["ACH-NOPE"])


def test_blocked_subsample_keeps_patients_together():
    U = data.universe_models()
    pid = data.load_models().PatientID
    sub = data.blocked_subsample(U, 0.8, 1000)
    assert sub == data.blocked_subsample(U, 0.8, 1000)
    assert sub != data.blocked_subsample(U, 0.8, 1001)
    assert abs(len(sub) - round(0.8 * len(U))) <= 2
    rest = sorted(set(U) - set(sub))
    assert not set(pid.reindex(sub)) & set(pid.reindex(rest))


def test_blocked_split_stratified():
    U = data.universe_models()
    m = data.load_models().reindex(U)
    a, b = data.blocked_split(U, 0.5, 7, strata=m.OncotreeLineage)
    assert sorted(a + b) == sorted(U) and not set(a) & set(b)
    assert not set(m.PatientID.reindex(a)) & set(m.PatientID.reindex(b))
    ca, cb = m.loc[a].OncotreeLineage.value_counts(), m.loc[b].OncotreeLineage.value_counts()
    big = ca.add(cb, fill_value=0)
    big = big[big >= 30].index
    assert ((ca.reindex(big) - cb.reindex(big)).abs() <= 0.25 * big.map(lambda s: ca.get(s, 0) + cb.get(s, 0))).all()


def test_filter_rule_on_toy_and_real_size():
    X = np.zeros((200, 4))
    X[:2, 0] = -1          # 1% below -0.5 -> in
    X[:1, 1] = -1          # 0.5% -> out
    X[:2, 2] = -1
    X[:11, 2] = np.nan     # 94.5% observed -> out
    X[:, 3] = -2           # pan-essential -> in
    m, st = genefilter.dependency_bearing(X)
    assert m.tolist() == [True, False, False, True]
    m2, _ = genefilter.dependency_bearing(X, low_coverage=np.array([0, 0, 0, 1]))
    assert m2.tolist() == [True, False, False, False]
    S = genefilter.dependency_bearing_for_models(list(data.load_gene_effect().models))
    assert S.sum() == 6163                                       # round-2 algorithms review, all models
