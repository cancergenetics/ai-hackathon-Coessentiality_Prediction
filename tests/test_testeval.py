"""Tests for the test-evaluation harness (src/testeval). Toy data with known answers, positive / negative controls,
equivalence with sklearn and with the locked Part 1 functions. No hu.MAP test file is read here."""
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from sklearn.isotonic import IsotonicRegression  # noqa: E402
from sklearn.metrics import average_precision_score, roc_auc_score  # noqa: E402

from testeval import bootstrap as B  # noqa: E402
from testeval import calibration as C  # noqa: E402
from testeval import evaluate as E  # noqa: E402
from testeval import pairs as P  # noqa: E402
from testeval.metrics import Ranked, top_cluster_ids  # noqa: E402
from testeval.part1 import metrics as p1m  # noqa: E402

Y5 = np.array([1, 0, 1, 0, 0])
S5 = np.array([0.9, 0.8, 0.7, 0.6, 0.5])


# ------------------------------------------------------------------------------------------- AP / AUROC / ties
def test_ap_hand_computed():
    # ranks: + - + - -  -> AP = 1/2 * 1 + 1/2 * 2/3
    want = 0.5 * 1 + 0.5 * (2 / 3)
    assert Ranked(Y5, S5).ap() == pytest.approx(want, abs=1e-15)
    assert p1m.weighted_ap(Y5, S5) == pytest.approx(want, abs=1e-15)
    assert average_precision_score(Y5, S5) == pytest.approx(want, abs=1e-15)


def test_ap_ties_are_blocks():
    y, s = np.array([1, 0, 1, 0]), np.array([1.0, 1.0, 0.0, 0.0])
    # block 1: P = 1/2 at R = 1/2; block 2: P = 2/4 at R = 1  -> AP = 0.5
    assert Ranked(y, s).ap() == pytest.approx(0.5)
    assert average_precision_score(y, s) == pytest.approx(0.5)
    yc = np.r_[np.ones(7), np.zeros(93)]
    assert Ranked(yc, np.zeros(100)).ap() == pytest.approx(0.07)            # constant scorer = prevalence


def test_auroc_hand_and_ties():
    assert Ranked(Y5, S5).auroc() == pytest.approx(5 / 6)
    assert Ranked(np.array([1, 0]), np.array([1.0, 1.0])).auroc() == pytest.approx(0.5)


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_weighted_ap_auroc_match_sklearn_with_ties(seed):
    rng = np.random.default_rng(seed)
    n = 3000
    y = rng.random(n) < 0.05
    s = np.round(rng.normal(size=n) + 1.5 * y, 1)                         # many ties
    w = rng.integers(0, 4, n).astype(float)
    r = Ranked(y, s)
    assert r.ap(w) == pytest.approx(average_precision_score(y, s, sample_weight=w), abs=1e-12)
    assert r.ap(w) == pytest.approx(p1m.PresortedScores(y, s).ap(w), abs=1e-14)
    assert r.auroc(w) == pytest.approx(roc_auc_score(y, s, sample_weight=w), abs=1e-12)
    # integer weights == replication
    rep = np.repeat(np.arange(n), w.astype(int))
    assert r.ap(w) == pytest.approx(average_precision_score(y[rep], s[rep]), abs=1e-12)


# ------------------------------------------------------------------------------------------- decision rules
def test_confusion_at_threshold():
    c = Ranked(Y5, S5).confusion(0.7)
    assert (c["tp"], c["fp"], c["fn"], c["tn"]) == (2, 1, 0, 2)
    assert c["precision"] == pytest.approx(2 / 3) and c["recall"] == 1 and c["f1"] == pytest.approx(0.8)
    c = Ranked(Y5, S5).confusion(0.95)                                      # above the maximum: nothing predicted
    assert c["n_pred_pos"] == 0 and c["recall"] == 0


def test_f1_max_threshold_and_tie_rule():
    f = Ranked(Y5, S5).f1_max()
    assert f["threshold"] == 0.7 and f["f1"] == pytest.approx(0.8)
    # F1 tie: thresholds 0.9 (TP1, FP0: 2/3) and 0.6 (TP2, FP2: 4/6 = 2/3) -> the higher threshold wins
    y, s = np.array([1, 0, 0, 1, 0]), np.array([0.9, 0.8, 0.7, 0.6, 0.5])
    f = Ranked(y, s).f1_max()
    assert f["threshold"] == 0.9 and f["f1"] == pytest.approx(2 / 3)


def test_topk_rule():
    r = Ranked(Y5, S5)
    t, K = E.topk_threshold(r, np.ones(5), 0.4)                             # K = round(0.4 * 5) = 2
    assert K == 2 and t == 0.8
    t, K = E.topk_threshold(Ranked(Y5, np.array([0.9, 0.8, 0.8, 0.6, 0.5])), np.ones(5), 0.4)
    assert t == 0.8                                                         # whole tie block at the boundary


def test_top_cluster_ids():
    y = np.array([1, 1, 1, 1, 0, 1])
    ca, cb = np.array([3, 3, 3, 5, 5, 7]), np.array([3, 3, 3, 5, 9, 8])
    assert top_cluster_ids(y, ca, cb, 2) == [3, 5]


# ------------------------------------------------------------------------------------------- calibration
def test_isotonic_toy_and_sklearn_equivalence():
    s, y = np.array([1.0, 2, 3, 4]), np.array([0.0, 1, 0, 1])
    cal = C.fit_isotonic(s, y)
    assert np.allclose(C.apply(cal, s), [0, 0.5, 0.5, 1])                   # pool-adjacent-violators
    rng = np.random.default_rng(0)
    s = rng.normal(size=500)
    y = (rng.random(500) < 1 / (1 + np.exp(-2 * s))).astype(float)
    w = rng.random(500)
    ir = IsotonicRegression(y_min=0, y_max=1, out_of_bounds="clip").fit(s, y, sample_weight=w)
    new = np.r_[rng.normal(size=200) * 2, -10, 10]
    assert np.allclose(C.apply(C.fit_isotonic(s, y, w), new), ir.predict(new), atol=1e-12)


def test_platt_recovers_parameters():
    rng = np.random.default_rng(1)
    s = rng.normal(2.0, 3.0, 200_000)
    p = 1 / (1 + np.exp(-(-4.0 + 0.8 * s)))
    y = (rng.random(len(s)) < p).astype(float)
    cal = C.fit_platt(s, y)
    a_raw = cal["a"] - cal["b"] * cal["mu"] / cal["sd"]
    b_raw = cal["b"] / cal["sd"]
    assert a_raw == pytest.approx(-4.0, abs=0.08) and b_raw == pytest.approx(0.8, abs=0.02)


def test_calibration_metrics_known_answers():
    p, y = np.array([0.1, 0.9]), np.array([0.0, 1.0])
    assert C.brier(p, y) == pytest.approx(0.01)
    rng = np.random.default_rng(2)
    p = rng.beta(0.5, 20, 400_000)
    y = (rng.random(len(p)) < p).astype(float)                              # perfectly calibrated forecasts
    sm = C.summary(p, y, ref_prevalence=float(y.mean()))
    assert sm["ece_mass20"] < 0.002 and sm["calibration_slope"] == pytest.approx(1, abs=0.03)
    assert sm["calibration_in_the_large"] == pytest.approx(1, abs=0.02)
    const = C.summary(np.full(len(y), y.mean()), y, ref_prevalence=float(y.mean()))
    assert const["bss_vs_fit_prevalence"] == pytest.approx(0, abs=1e-12)
    # ECE hand example: 2 equal-mass bins {0.1, 0.2} (obs 0) and {0.8, 0.9} (obs 1)
    e = C.ece(np.array([0.1, 0.2, 0.8, 0.9]), np.array([0, 0, 1, 1.0]), n_bins=2)
    assert e == pytest.approx(0.5 * 0.15 + 0.5 * 0.15)


def test_floor_to_fit():
    s = C.floor_to_fit(np.array([3.0, -7.0, 1.0]), np.array([True, False, True]), -4.5)
    assert s.tolist() == [3.0, -4.5, 1.0]


# ------------------------------------------------------------------------------------------- bootstraps
def _toy_clustered(seed, n_cl=40, size=6, between=True):
    rng = np.random.default_rng(seed)
    prot = np.arange(n_cl * size)
    cl = prot // size
    eff = rng.normal(0, 1.0, n_cl)                                          # cluster-level score effect
    a, b = np.triu_indices(len(prot), 1)
    same = cl[a] == cl[b]
    keep = same | (between & (rng.random(len(a)) < 0.05))
    a, b, same = a[keep], b[keep], same[keep]
    y = same & (rng.random(len(a)) < 0.6)
    s = 1.2 * y + np.where(same, eff[cl[a]], 0) + rng.normal(size=len(a))
    return y, s, cl[a], cl[b], a, b


def test_cluster_generator_equals_part1_cluster_bootstrap():
    y, s, ca, cb, _, _ = _toy_clustered(0)
    s2 = s + np.random.default_rng(9).normal(size=len(s))
    f = 0.3
    w = p1m.ht_weights(ca, cb, f)
    ref = p1m.cluster_bootstrap(y, {"x": s, "z": s2}, ca, cb, f=f, n_boot=300, seed=4, diffs=[("x", "z")])
    fns = {"x": E._ap_fn(y, s), "z": E._ap_fn(y, s2)}
    got = B.run(B.cluster_weights(ca, cb, 300, 4), w, fns, [("x|ap", "z|ap")])
    for k in ("x", "z"):
        assert got["stat"][f"{k}|ap"]["point"] == pytest.approx(ref["ap"][k]["ap"], abs=1e-14)
        assert got["stat"][f"{k}|ap"]["ci_low"] == pytest.approx(ref["ap"][k]["ci_low"], abs=1e-14)
        assert got["stat"][f"{k}|ap"]["ci_high"] == pytest.approx(ref["ap"][k]["ci_high"], abs=1e-14)
    d, rd = got["diff"]["x|ap - z|ap"], ref["diff"]["x-z"]
    for k in ("ci_low", "ci_high", "p_two_sided", "rel_ci_low"):
        assert d[k] == pytest.approx(rd[k], abs=1e-14)


def test_protein_and_pair_weights():
    pa, pb = np.array(["A", "A", "B"]), np.array(["B", "C", "C"])
    for w in B.protein_weights(pa, pb, 50, 0):
        assert w.shape == (3,) and (w >= 0).all()
    rng = np.random.default_rng(0)
    m = np.bincount(rng.integers(0, 3, 3), minlength=3).astype(float)       # units sorted: A, B, C
    w0 = next(B.protein_weights(pa, pb, 1, 0))
    assert np.allclose(w0, [m[0] * m[1], m[0] * m[2], m[1] * m[2]])
    y = np.array([1, 0, 0, 1, 0, 0, 0])
    for w in B.pair_weights(y, 20, 0):                                      # stratified: class totals fixed
        assert w[y == 1].sum() == 2 and w[y == 0].sum() == 5


def test_bootstrap_ci_coverage_sanity():
    """Clustered simulation: true population values from a 30x larger draw; the 95% cluster-bootstrap CI must cover
    the truth in roughly 95% of 120 replicates (accept 0.85-1.0) for a mean statistic (Brier-type) and for AP."""
    rng = np.random.default_rng(123)

    def draw(n_cl, rng):
        size = rng.integers(2, 9, n_cl)
        cl = np.repeat(np.arange(n_cl), size * (size - 1) // 2)
        eff = rng.normal(0, 0.8, n_cl)[cl]
        y = rng.random(len(cl)) < 0.3
        s = y * 1.0 + eff + rng.normal(size=len(cl))
        return y, s, cl

    yb, sb, _ = draw(6000, rng)
    p_true = 1 / (1 + np.exp(-sb))
    truth_mean = float(np.mean((p_true - yb) ** 2))
    truth_ap = Ranked(yb, sb).ap()
    cov_m = cov_ap = 0
    n_rep = 120
    for r in range(n_rep):
        y, s, cl = draw(80, rng)
        p = 1 / (1 + np.exp(-s))
        fn = {"m": lambda w, p=p, y=y: dict(brier=float(np.sum(w * (p - y) ** 2) / np.sum(w))), "a": E._ap_fn(y, s)}
        out = B.run(B.cluster_weights(cl, cl, 300, r), np.ones(len(y)), fn)
        cov_m += out["stat"]["m|brier"]["ci_low"] <= truth_mean <= out["stat"]["m|brier"]["ci_high"]
        cov_ap += out["stat"]["a|ap"]["ci_low"] <= truth_ap <= out["stat"]["a|ap"]["ci_high"]
    assert 0.85 <= cov_m / n_rep <= 1.0
    assert 0.85 <= cov_ap / n_rep <= 1.0


# ------------------------------------------------------------------------------------------- controls
def test_positive_control_perfect_scorer():
    rng = np.random.default_rng(3)
    y = rng.random(5000) < 0.01
    s = y + 0.1 * rng.random(5000)
    r = Ranked(y, s)
    assert r.ap() == 1.0 and r.auroc() == 1.0
    assert r.f1_max()["f1"] == 1.0


def test_negative_controls_random_and_label_shuffle():
    rng = np.random.default_rng(4)
    n = 400_000
    y = rng.random(n) < 0.005
    assert Ranked(y, rng.normal(size=n)).ap() == pytest.approx(y.mean(), rel=0.12)
    informative = 3 * y + rng.normal(size=n)
    assert Ranked(y, informative).ap() > 0.5
    ys = rng.permutation(y)                                                 # label shuffle
    assert Ranked(ys, informative).ap() == pytest.approx(ys.mean(), rel=0.12)
    assert Ranked(ys, informative).auroc() == pytest.approx(0.5, abs=0.02)


# ------------------------------------------------------------------------------------------- pair handling
def _write(path, rows):
    path.write_text("".join(f"{a}\t{b}\n" for a, b in rows))
    return path


def test_labelled_table_rules(tmp_path):
    pos = _write(tmp_path / "pos.txt", [("B", "A"), ("A", "B"), ("C", "D"), ("E", "E"), ("F", "G")])
    neg = _write(tmp_path / "neg.txt", [("A", "C"), ("G", "F"), ("C", "A"), ("H", "I")])
    trp = _write(tmp_path / "trp.txt", [("D", "C")])
    trn = _write(tmp_path / "trn.txt", [("A", "X")])
    t, au = P.labelled_table(pos, neg, (trp, trn))
    assert au["pos_lines"] == 5 and au["pos_unique_pairs"] == 4 and au["pos_reversed_duplicates"] == 2
    assert au["self_pairs"] == 1 and au["conflicts"] == 1                   # E-E; F-G in both files
    an = t[t.analysed].set_index(["a", "b"])
    assert set(an.index) == {("A", "B"), ("C", "D"), ("A", "C"), ("H", "I")}
    assert an.loc[("A", "B"), "label"] == 1 and an.loc[("A", "C"), "label"] == 0
    assert an.loc[("C", "D"), "in_train_pos"] and not an.loc[("A", "B"), "in_train_pos"]
    assert au["train_overlap"]["pos_in_train_pos"] == 1
    # train proteins: C, D, A, X -> A-B has one unseen protein; only H-I (neg) has both unseen
    assert au["pairs_both_proteins_unseen_in_train"] == dict(pos=0, neg=1)
    assert au["proteins_analysed"] == 6 and au["proteins_also_in_train"] == 3


def test_id_file_has_no_labels_no_self_pairs(tmp_path):
    pos = _write(tmp_path / "pos.txt", [("B", "A"), ("E", "E")])
    neg = _write(tmp_path / "neg.txt", [("A", "B"), ("C", "D")])
    t, _ = P.labelled_table(pos, neg)
    rec = P.write_id_file([t], tmp_path / "ids.tsv")
    lines = (tmp_path / "ids.tsv").read_text().splitlines()
    assert lines == ["A\tB", "C\tD"] and rec["n_pairs"] == 2


def test_read_and_attach_scores(tmp_path):
    f = tmp_path / "s.tsv"
    f.write_text("uniprot_a\tuniprot_b\tscore\tmapped\nB\tA\t0.5\t1\nC\tD\t-2.0\t0\n")
    sc = P.read_scores(f)
    t = pd.DataFrame(dict(a=["A", "C"], b=["B", "D"]))
    s, m = P.attach_scores(t, sc, "x")
    assert s.tolist() == [0.5, -2.0] and m.tolist() == [True, False]
    with pytest.raises(ValueError):
        P.attach_scores(pd.DataFrame(dict(a=["A", "Q"], b=["B", "R"])), sc, "x")
    f.write_text("uniprot_a\tuniprot_b\tscore\tmapped\nB\tA\t0.5\t1\nA\tB\t0.6\t1\n")
    with pytest.raises(ValueError):
        P.read_scores(f)


# ------------------------------------------------------------------------------------------- pipeline on a toy
def _toy_set(tmp_path, seed=0, n_cl=30):
    rng = np.random.default_rng(seed)
    prots = [f"P{i:04d}" for i in range(n_cl * 5)]
    cl = np.repeat(np.arange(n_cl), 5)
    a, b = np.triu_indices(len(prots), 1)
    y = (cl[a] == cl[b]) & (rng.random(len(a)) < 0.8)
    pos = [(prots[i], prots[j]) for i, j in zip(a[y], b[y])]
    neg = [(prots[i], prots[j]) for i, j in zip(a[~y], b[~y])]
    _write(tmp_path / "pos.txt", pos)
    _write(tmp_path / "neg.txt", neg)
    unm = np.isin(a, [0, 7]) | np.isin(b, [0, 7])                           # two "unmapped" proteins
    sd = tmp_path / "scores"
    sd.mkdir()
    for k, m in enumerate(E.METHODS):
        s = (1.5 - 0.2 * k) * y + rng.normal(size=len(a))
        floor = s[~unm].min() - 1
        s = np.where(unm, floor, s)
        pd.DataFrame(dict(uniprot_a=[prots[i] for i in a], uniprot_b=[prots[j] for j in b], score=s,
                          mapped=(~unm).astype(int))).to_csv(sd / f"{m}.tsv", sep="\t", index=False)
    return tmp_path / "pos.txt", tmp_path / "neg.txt", sd


def test_pipeline_end_to_end_toy(tmp_path, monkeypatch):
    monkeypatch.setenv("PYTHONHASHSEED", "0")
    pos, neg, sd = _toy_set(tmp_path)
    es = E.load_set("toy", pos, neg, sd)
    rules = E.fit_rules(es)
    for m, r in rules["methods"].items():
        assert r["fit_floor"] < es.scores[m][es.mapped[m]].min()
    cl = E.annotate_clusters(es)
    res = E.evaluate(es, rules, cl, None, n_boot=50, mc_seeds=(1,))
    pm = res["metrics"]
    assert res["n_unmapped_any"] == int((~es.mapped_all).sum()) > 0
    assert pm["zavg_BEF"]["ap"] == pytest.approx(average_precision_score(es.y, es.scores["zavg_BEF"]))
    assert pm["zavg_BEF"]["ap"] > pm["pearson_signed"]["ap"]                # planted ordering of effect sizes
    assert set(pm["zavg_BEF"]["decision_rules"]) == {"rule_a_z3", "rule_b_f1max_dev", "rule_c_topk"}
    assert "rule_a_z3" not in pm["gls_signed"]["decision_rules"]
    assert res["bootstrap"]["cluster"]["n_boot"] == 50
    assert "zavg_BEF|ap - gls_signed|ap" in res["bootstrap"]["cluster"]["diff"]
    v = res["verdict"]
    assert v["candidates"]["zavg_BEF"]["status"] == "confirmatory"
    # mapped-only subset drops the unmapped pairs
    sub = es.subset(es.mapped_all, "m")
    assert len(sub.table) == int(es.mapped_all.sum())


def test_cluster_construction_requires_hashseed(tmp_path, monkeypatch):
    pos, neg, sd = _toy_set(tmp_path, n_cl=4)
    es = E.load_set("toy", pos, neg, sd)
    monkeypatch.setenv("PYTHONHASHSEED", "1")
    with pytest.raises(RuntimeError):
        E.annotate_clusters(es)


def test_verdict_logic_fixed_sequence():
    def pm(ap, macro, nl):
        return dict(ap=ap, macro_cluster_ap=macro, ap_no_largest_component=nl)
    metrics = {"zavg_BEF": pm(0.40, 0.30, 0.40), "B": pm(0.36, 0.20, 0.38),
               "gls_signed": pm(0.30, 0.31, 0.30), "spearman_signed": pm(0.20, 0.25, 0.2),
               "pearson_signed": pm(0.19, 0.24, 0.2)}
    diff = {}
    for c in ("zavg_BEF", "B"):
        for b in E.BASELINES:
            diff[f"{c}|ap - {b}|ap"] = dict(diff=0.1, ci_low=0.01, ci_high=0.2, p_two_sided=0.01)
    v = E.verdict(metrics, dict(cluster=dict(diff=diff)))
    assert v["pairs"]["zavg_BEF vs gls_signed"]["passes"]                   # macro 0.30 >= 0.95 * 0.31
    assert v["candidates"]["zavg_BEF"]["clearly_beats"]
    assert not v["pairs"]["B vs gls_signed"]["passes"]                      # macro 0.20 < 0.2945
    assert v["candidates"]["B"]["status"] == "confirmatory" and not v["candidates"]["B"]["clearly_beats"]
    diff["zavg_BEF|ap - gls_signed|ap"]["ci_low"] = -0.01                   # CI includes 0 -> primary fails
    v = E.verdict(metrics, dict(cluster=dict(diff=diff)))
    assert not v["candidates"]["zavg_BEF"]["clearly_beats"]
    assert v["candidates"]["B"]["status"] == "exploratory"


# ------------------------------------------------------------------------------------------- process guards
def test_only_entry_points_touch_test_files():
    """Test-file names appear only in testeval/paths.py; test_files() is called only by the two entry points."""
    hits, calls = [], []
    for p in list((ROOT / "src").rglob("*.py")) + list((ROOT / "scripts").glob("*")) + list((ROOT / "tests").glob("*.py")):
        if p.is_dir() or p.name == "test_testeval.py":
            continue
        txt = p.read_text(errors="ignore")
        if re.search(r"Test_(positive|negative)", txt):
            hits.append(p.name)
        if "test_files(" in txt:
            calls.append(p.name)
    assert hits == ["paths.py"]
    assert sorted(set(calls) - {"paths.py"}) == ["19_make_pair_files.py", "21_evaluate_test.py"]


def test_part1_imports_are_the_locked_clone():
    from testeval import part1
    st = part1.lock_state()
    assert st["matches"], st
    import src
    assert str(Path(src.__file__).resolve()).startswith(str((ROOT / "external" / "part1_locked").resolve()))


# ------------------------------------------------------------------------------------------- the look, on a toy
def test_mode_test_end_to_end_on_toy(tmp_path, monkeypatch):
    """Runs scripts/21_evaluate_test.py's look (mode_test) on toy files with the git gate stubbed: the look writes
    a 'started' record before reading labels, a result with the protocol commit, a 'done' record, and refuses a
    second look."""
    import importlib.util
    import json
    monkeypatch.setenv("PYTHONHASHSEED", "0")                              # no re-exec on import
    spec = importlib.util.spec_from_file_location("eval21", ROOT / "scripts" / "21_evaluate_test.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    from testeval import paths
    (tmp_path / "tr").mkdir()
    (tmp_path / "te").mkdir()
    trp, trn, trs = _toy_set(tmp_path / "tr", seed=1, n_cl=10)
    tep, ten, tes = _toy_set(tmp_path / "te", seed=2, n_cl=10)
    res = tmp_path / "results"
    monkeypatch.setattr(paths, "TRAIN_POS", trp)
    monkeypatch.setattr(paths, "TRAIN_NEG", trn)
    monkeypatch.setattr(paths, "SCORES", tmp_path / "scores")
    (tmp_path / "scores").mkdir()
    (tmp_path / "scores" / "test").symlink_to(tes)
    monkeypatch.setattr(paths, "RESULTS", res)
    monkeypatch.setattr(paths, "FITTED_RULES", res / "fitted_rules.json")
    monkeypatch.setattr(paths, "LOOKS", res / "looks.jsonl")
    monkeypatch.setattr(paths, "PROJECT", tmp_path)
    monkeypatch.setattr(paths, "test_files", lambda: (tep, ten))
    frozen = dict(head="h" * 40, protocol_last_commit="p" * 40, protocol_sha256="s" * 64, dirty=[])
    monkeypatch.setattr(mod, "require_frozen", lambda: frozen)
    monkeypatch.setattr(mod, "git_state", lambda paths=None: dict(dirty=[]))
    monkeypatch.setattr(mod, "load_matrix", lambda: (np.zeros((5, 1), np.float32), pd.Index(["AAA (1)"]), {}))
    fit = E.load_set("toy train", trp, trn, trs)
    res.mkdir()
    (res / "fitted_rules.json").write_text(json.dumps(E.fit_rules(fit), default=float))
    monkeypatch.setattr(E, "evaluate", lambda *a, **k: _orig_evaluate(*a, **dict(k, n_boot=40)))   # speed only
    mod.mode_test()
    looks = [json.loads(l) for l in (res / "looks.jsonl").read_text().splitlines()]
    assert [r["status"] for r in looks] == ["started", "done"]
    out = json.loads((res / "test_results.json").read_text())
    assert out["protocol"]["commit"] == "p" * 40 and out["confirmatory"] is True
    assert "verdict" in out["results"]["itt"] and "mapped_only" in out["results"]
    assert (res / "figures" / "reliability_test_itt.png").exists()
    with pytest.raises(SystemExit):
        mod.mode_test()                                                     # single look


_orig_evaluate = E.evaluate


def test_log_bins():
    p = np.array([1e-6, 2e-4, 0.05, 0.9])
    rows = C.reliability(p, np.array([0, 0, 1, 1.0]), n_bins=4, scheme="log")   # edges 1e-4, 1e-3, 1e-2, 1e-1, 1
    assert [r["bin"] for r in rows] == [0, 2, 3] and rows[0]["n"] == 2
