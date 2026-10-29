"""Descriptive export of the test look (scripts/23_export_test_curves.py): toy checks of the curve / down-sampling /
per-cluster functions, and that the committed export reproduces the look's AP and macro AP. Reads only committed
result files; no label file is read here."""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from sklearn.metrics import average_precision_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
RES = ROOT / "results" / "test_eval"
EXPORT = RES / "export"

spec = importlib.util.spec_from_file_location("export23", ROOT / "scripts" / "23_export_test_curves.py")
X = importlib.util.module_from_spec(spec)
spec.loader.exec_module(X)


def _toy(seed=0, n=5000):
    rng = np.random.default_rng(seed)
    y = rng.random(n) < 0.1
    s = np.round(rng.normal(size=n) + 1.5 * y, 1)               # heavy ties
    return y, s


def test_downsampled_curve_integrates_to_exact_ap():
    y, s = _toy()
    full = X.pr_curve(y, s)
    ref = average_precision_score(y, s)
    assert abs(full.inc.sum() - ref) < 1e-12
    ds = X.downsample(full, max_points=50)
    assert len(ds) <= 50 and len(full) > 50
    assert abs(ds.ap_increment.sum() - ref) < 1e-12
    assert ds.recall.iloc[-1] == 1.0 and np.all(np.diff(ds.recall.values) >= 0)
    assert np.all(np.diff(ds.threshold.values) < 0)
    assert abs(X.naive_step_ap(full.recall.values, full.precision.values) - ref) < 1e-12


def test_per_cluster_mean_is_part1_macro_ap():
    from testeval.part1 import metrics as p1m
    y, s = _toy(1, 2000)
    rng = np.random.default_rng(2)
    ca, cb = rng.integers(0, 30, len(y)), rng.integers(0, 30, len(y))
    tab = pd.DataFrame(dict(protein=np.arange(60), cluster=np.repeat(np.arange(30), 2)))
    pc = X.per_cluster(y, {"m": s}, ca, cb, np.ones(len(y)), tab)
    ref = p1m.macro_cluster_ap(y, s, ca, cb, np.ones(len(y)))
    assert len(pc) == ref["n_clusters_macro"] and abs(pc.m.mean() - ref["macro_cluster_ap"]) < 1e-14


@pytest.mark.skipif(not (EXPORT / "pr_curves.tsv.gz").exists(), reason="export not run")
def test_committed_export_reproduces_look_ap_and_macro_ap():
    tr = json.loads((RES / "test_results.json").read_text())["results"]
    pr = pd.read_csv(EXPORT / "pr_curves.tsv.gz", sep="\t")
    sz = pd.read_csv(EXPORT / "ap_by_size.tsv", sep="\t")
    for u in ("itt", "mapped_only"):
        for m, ref in tr[u]["metrics"].items():
            c = pr[(pr.universe == u) & (pr.method == m)]
            assert 0 < len(c) <= X.MAX_POINTS
            assert abs(c.ap_increment.sum() - ref["ap"]) <= 1e-9
            for by in X.BINS:
                a = sz[(sz.universe == u) & (sz.method == m) & (sz.bin_by == by)]
                assert abs(a[a.bin == "all"].macro_ap.item() - ref["macro_cluster_ap"]) <= 1e-9
                b = a[(a.bin != "all") & (a.n_clusters > 0)]
                assert b.n_clusters.sum() == ref["n_clusters_macro"]
    recs = [json.loads(l) for l in (RES / "looks.jsonl").read_text().splitlines() if l.strip()]
    ex = [r for r in recs if r.get("type") == X.EXPORT_TYPE]
    assert ex and all(r["non_confirmatory"] is True and "status" not in r for r in ex)
    assert [r["status"] for r in recs if "status" in r] == ["intent", "started", "done"]
