"""Complex scores (DESIGN v3.1 §2.4), member removals, covariate assembly and the covariate-of-context rule, frozen map
round trip -- all on toy data with known answers; the machine-readable removal / control rules are checked against
the registry text. No dependency data is read."""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from part2 import controls, covariates as cv, cscore, registry  # noqa: E402
from part2.cmap import FrozenMap, symbol  # noqa: E402
from part2.cscore import Annotation, Removal  # noqa: E402

REG = registry.load_all(ROOT / "docs" / "registry")


def test_mean_score_clip_and_half_rule():
    E = pd.DataFrame({"a (1)": [-5.0, 1.0, np.nan, np.nan], "b (2)": [-1.0, 2.0, 4.0, np.nan],
                      "c (3)": [0.0, np.nan, np.nan, np.nan], "d (4)": [1.0, np.nan, 0.0, 1.0]})
    s = cscore.mean_score(E, list(E.columns))
    assert s.iloc[0] == pytest.approx((-3 - 1 + 0 + 1) / 4)            # clipped at -3
    assert s.iloc[1] == pytest.approx(1.5)                             # 2 of 4 observed = 50% -> kept
    assert s.iloc[2] == pytest.approx((3 + 0) / 2)                     # clipped at +3
    assert np.isnan(s.iloc[3])                                         # 1 of 4 < 50%


def annot():
    return Annotation(paralogs={"ARID1A": {"ARID1B"}, "ARID1B": {"ARID1A"}, "DDX3X": {"DDX3Y"}, "DDX3Y": {"DDX3X"}},
                      chrom={"ARID1A": "1", "ARID1B": "6", "NEAR": "1", "FAR": "1", "DDX3X": "X", "DDX3Y": "Y",
                             "XOTHER": "X"},
                      pos={"ARID1A": 26_700_000, "ARID1B": 157_000_000, "NEAR": 30_000_000, "FAR": 60_000_000,
                           "DDX3X": 41_000_000, "DDX3Y": 12_000_000, "XOTHER": 1},
                      arm={"NOPOS": "1p", "ARID1A": "1p"})


def test_removals():
    a = annot()
    mem = ["ARID1A (1)", "ARID1B (2)", "NEAR (3)", "FAR (4)", "NOPOS (5)", "DDX3X (6)", "DDX3Y (7)", "XOTHER (8)"]
    assert Removal(genes=("ARID1A",)).removed(mem, a) == ["ARID1A (1)", "ARID1B (2)"]           # + in-complex paralog
    assert Removal(genes=("ARID1A",), paralogs=False).removed(mem, a) == ["ARID1A (1)"]
    # +/-10 Mb window around a CN gene (NEAR at 3.3 Mb; FAR at 33 Mb); arm fallback for a gene without coordinates
    assert Removal(genes=("ARID1A",), cn_genes=("ARID1A",)).removed(mem, a) == \
        ["ARID1A (1)", "ARID1B (2)", "NEAR (3)", "NOPOS (5)"]
    assert Removal(genes=("ARID1A",), retain=("ARID1B",)).removed(mem, a) == ["ARID1A (1)"]
    assert Removal(chrY=True).removed(mem, a) == ["DDX3X (6)", "DDX3Y (7)"]                  # Y + X gametolog
    assert Removal(ygenes=True).removed(mem, a) == ["DDX3Y (7)"]


def test_removal_rules_match_registry_text():
    """Every H3a/H3b/H2/BENCH/BRIDGE context has a rule, and every gene of the rule is named in its registry text."""
    n = 0
    for k, e in REG.items():
        cs = e["context_family"]["contexts"]
        if e["direction"] in ("controls", "H1", "H1S") or cs == "FZ-X":   # H1 / H1S: component rules (h1.py)
            continue
        for c in cs:
            r = cscore.removal_for_context(c["id"], homdel_share=1.0)
            txt = c.get("removed_members", "")
            for g in r.genes + r.cn_genes:
                assert g in txt, (k, c["id"], g, txt)
            if r.cn_genes:
                assert "10 Mb" in txt, (k, c["id"])
            if r.chrY:
                assert "chrY" in txt and "paralog" in txt
            if not (r.genes or r.cn_genes or r.chrY) and txt:
                assert txt.startswith("none"), (k, c["id"], txt)
            n += 1
    assert n > 100
    assert cscore.removal_for_context("LoF:ARID1A", homdel_share=0.2).cn_genes == ()
    with pytest.raises(KeyError):
        cscore.removal_for_context("unregistered_context")


def test_control_specs_match_registry_text():
    for k, e in REG.items():
        if e["direction"] != "controls" or e["score"] == "gene":
            continue
        spec = controls.CONTROL_SPECS[e["control"]["design_id"]]
        txt = e["control"]["removed_members"] + " " + e["control"].get("readout_C", "")
        rem = spec.get("removal")
        genes = (rem.genes + rem.cn_genes + rem.retain) if rem else ()
        lost = spec.get("lost", ())
        for g in genes + (tuple(lost) if lost != "chrY" else ()):
            assert g in txt, (k, g, txt)
        if rem is not None and rem.cn_genes:
            assert "10 Mb" in txt
        if rem is not None and rem.retain:
            assert "RETAINED" in txt or "retained" in txt, (k, txt)
        if "target" in spec:
            assert spec["target"] in e["control"]["targets"], k


def test_multiverse_scores():
    rng = np.random.default_rng(0)
    n, k = 400, 6
    f = rng.normal(size=n)
    load = np.linspace(0.5, 1.0, k)
    X = -1 + np.outer(f, load) + rng.normal(0, 0.3, (n, k))
    cols = [f"g{i} ({i})" for i in range(k)]
    E = pd.DataFrame(X, columns=cols)
    eg, fit = cscore.score(E, cols, "eigengene")
    assert abs(np.corrcoef(eg, f)[0, 1]) > 0.95
    assert np.corrcoef(eg, cscore.mean_score(E, cols))[0, 1] > 0.95              # same sign as the mean
    eg2, _ = cscore.score(E.iloc[:50], cols, "eigengene", fit=fit)                 # A-fitted params re-used
    assert np.allclose(eg2, eg.iloc[:50])
    full = pd.concat([E, pd.DataFrame(rng.normal(0, 1, (n, 20)), columns=[f"o{i} (9{i})" for i in range(20)])], axis=1)
    r, _ = cscore.score(full, cols, "within_line_rank", ranks=full.rank(axis=1, pct=True))
    assert r.between(0, 1).all() and r.mean() < 0.3                                  # members are the dependent ones
    expr = pd.DataFrame(5.0, index=E.index, columns=cols)
    expr.iloc[:, 0] = 0.0
    m, _ = cscore.score(E, cols, "expression_masked", expr=expr)
    assert np.allclose(m, cscore.mean_score(E, cols[1:]))
    pan = pd.DataFrame(-2.0, index=E.index, columns=["p1 (1001)", "p2 (1002)"])
    pr, _ = cscore.score(pd.concat([E, pan], axis=1), cols, "pan_essential_referenced", pan_genes=list(pan.columns))
    assert np.allclose(pr, cscore.mean_score(E, cols) / 2.0)


def test_unit_scores_evaluability():
    fm = FrozenMap(members={"C0": ["ARID1A (1)", "ARID1B (2)", "FAR (4)"], "C1": ["NEAR (3)", "FAR (4)"]},
                   level={"C0": 1, "C1": 1})
    E = pd.DataFrame(np.ones((5, 4)), columns=["ARID1A (1)", "ARID1B (2)", "NEAR (3)", "FAR (4)"])
    Y, info, _ = cscore.unit_scores(E, fm, ["C0", "C1"], Removal(genes=("ARID1A",)), annot())
    assert info.set_index("unit").evaluable.to_dict() == {"C0": False, "C1": True}   # C0 keeps only FAR
    assert Y["C0"].isna().all() and (Y["C1"] == 1).all()


def test_covariate_rule_and_assembly():
    assert "lineage_FE" in cv.context_functions("lineage:Bone")
    assert cv.context_functions("melanoma") >= {"lineage_FE"}
    assert cv.context_functions("MSI_H") >= {"MSIScore", "log_burden"}
    assert "Ploidy" in cv.context_functions("WGD")
    assert "TP53_altered" in cv.context_functions("TP53_altered")
    assert "TP53_altered" in cv.context_functions("H1loss:TP53")
    assert not cv.context_functions("LoF_inclMSI:ARID1A")
    # the registered entries already follow the rule: the guard drops nothing
    for k, e in REG.items():
        if e["context_family"]["contexts"] == "FZ-X":
            continue
        for c in e["context_family"]["contexts"]:
            kept, dropped = cv.covariates_for(e, c["id"])
            assert not dropped, (k, c["id"], dropped)
    tab = pd.DataFrame(dict(OncotreeLineage=["A", "B", "B", None], LibraryAvana=[1, 0, 1, 1], LibraryKY=[0, 1, 1, 0],
                            ScreenType3DO=[0, 0, 1, 0], ModelType=["x", "x", "y", "x"], ScreenNNMD=[1.0, 2, 3, 4]))
    X = cv.assemble(["lineage_FE", "library", "ScreenType3DO", "ModelType", "NNMD"], tab)
    assert list(X.columns) == ["lineage_FE:B", "library:LibraryAvana", "library:LibraryKY", "ScreenType3DO",
                               "ModelType:y", "NNMD"]
    assert X.iloc[3].isna().iloc[0] and X.iloc[1, 0] == 1


def test_design_drops_collinear_and_rare_levels():
    from part2.assoc import _design
    n = 40
    ctx = pd.Series(np.r_[np.ones(15), np.zeros(25)])
    cov = pd.DataFrame({"ScreenType3DO": np.r_[np.zeros(30), np.ones(10)]})
    cov["ModelType:y"] = cov.ScreenType3DO                                             # identical -> dropped
    cov["lineage_FE:rare"] = np.r_[1, np.zeros(n - 1)]                               # 1 model -> reference level
    cov["const"] = 1.0
    X, names, ok, dropped = _design(ctx, cov, return_dropped=True)
    assert names == ["intercept", "context", "ScreenType3DO"]
    assert set(dropped) == {"ModelType:y", "lineage_FE:rare", "const"}


def test_frozen_map_round_trip_and_cluster_map():
    from part2.consensus import ClusterMap
    cm = ClusterMap(clusters=[dict(id=0, level=1, parent=-1, members=np.array([0, 1, 2, 3, 4, 5])),
                              dict(id=1, level=1, parent=-1, members=np.array([6, 7, 8])),
                              dict(id=2, level=2, parent=0, members=np.array([0, 1, 2])),
                              dict(id=3, level=2, parent=0, members=np.array([3, 4]))])
    keys = [f"G{i} ({i})" for i in range(9)]
    fm = FrozenMap.from_cluster_map(cm, keys, locus={"C0001": True})
    assert fm.coarse() == ["C0000", "C0001"] and fm.submodules("C0000") == ["C0000.S1", "C0000.S2"]
    assert fm.sibling_complement("C0000.S1") == ["G3 (3)", "G4 (4)", "G5 (5)"]
    fm.attached["C0001"] = ["X (99)"]
    back = FrozenMap.from_membership(fm.to_membership())
    assert back.members == fm.members and back.parent == fm.parent and back.attached == fm.attached
    assert back.locus["C0001"] and not back.locus["C0000"]
    assert fm.units_of_symbol("G7") == ["C0001"] and symbol("ABC (12)") == "ABC"


@pytest.mark.skipif(not (ROOT / "datasets" / "Paralogs").exists(), reason="annotation files not available")
def test_load_annotation_known_pairs():
    """Ensembl 116 annotation (labels only): the registry's canonical in-complex paralog pairs pass the H1 rule."""
    keys = ["ARID1A (8289)", "ARID1B (57492)", "SMARCA4 (6597)", "SMARCA2 (6595)", "VPS4B (9525)", "VPS4A (27183)",
            "DDX3X (1654)", "DDX3Y (8653)", "STAG2 (10735)", "STAG1 (10274)", "MTAP (4507)", "CDKN2A (1029)"]
    a = cscore.load_annotation(keys)
    for x, y in [("ARID1A", "ARID1B"), ("SMARCA4", "SMARCA2"), ("VPS4B", "VPS4A"), ("STAG2", "STAG1")]:
        assert a.is_paralog(x, y), (x, y)
    assert a.chrom["DDX3Y"] == "Y" and a.x_gametolog("DDX3X")
    assert a.near("MTAP", "CDKN2A") and a.arm["MTAP"] == "9p"
