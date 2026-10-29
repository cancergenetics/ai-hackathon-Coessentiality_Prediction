"""Descriptive export of the single confirmatory test look (NON-CONFIRMATORY; not a relook).

The look (scripts/21_evaluate_test.py --mode test --run, results/test_eval/test_results.json) computed PR curves and
per-cluster APs internally but did not export them. This script re-derives two figure inputs from exactly the same
inputs and definitions, and changes nothing that the look produced:

  results/test_eval/export/pr_curves.tsv.gz  PR curve of every method on ITT and S1 (mapped-only), <= MAX_POINTS
                                             points per method x universe; column ap_increment carries the exact
                                             step-integral mass accrued since the previous exported point, so
                                             sum(ap_increment) = AP of the look
  results/test_eval/export/ap_by_size.tsv    macro per-cluster AP (Part 1 macro_cluster_ap, per size bin) and pooled
                                             AP over the pairs touching the bin's clusters, per method x universe,
                                             with a cheap descriptive within-bin cluster-resampling CI
  results/test_eval/export/README.md

Inputs and definitions are the harness's own: labels through testeval.paths.test_files() + evaluate.load_set (as
in mode_test), cached score TSVs (SHA-256 checked against results/test_eval/score_records/test and test_results.json),
clusters through evaluate.annotate_clusters on ITT and on the S1 subset (protocol §5.2, as run_protocol does).

Sanity checks (asserted; the export is refused otherwise): audit counts and cluster statistics equal the look's;
AP integrated from the exported curves and macro AP equal test_results.json to 1e-9 for every method and universe.

No verdict, decision or pre-specified statistic depends on these outputs. Size bins are fixed below, before the
export is run. The script refuses a second export unless --rerun-reason is given (then logged as such).
Run: PART1_DIR=<Part 1 folder> .venv/bin/python scripts/23_export_test_curves.py
"""
import argparse
import gzip
import io
import json
import os
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

MAX_POINTS = 2000
TOL = 1e-9
N_BOOT = 2000                     # descriptive within-bin cluster resampling (not the protocol's bootstrap)
BOOT_SEED = 0
REF_BASELINE = "gls_signed"
UNIVERSES = {"itt": "ITT (all analysed pairs)", "mapped_only": "S1 (pairs mapped under every method)"}
# Size bins fixed before the export was run. A cluster's "touching positive pairs" are the positive pairs with >= 1
# endpoint in it (= its within-cluster positives when no component is Louvain-split, as on the test set).
BINS = {
    "n_pos_pairs": [("1", 1, 1), ("2-5", 2, 5), ("6-20", 6, 20), ("21-100", 21, 100), ("101+", 101, 10 ** 9)],
    "n_proteins": [("2", 2, 2), ("3-4", 3, 4), ("5-8", 5, 8), ("9-16", 9, 16), ("17+", 17, 10 ** 9)],
}
EXPORT_TYPE = "descriptive_export"
STATEMENT = ("descriptive presentation of the single confirmatory look; no decision or verdict depends on it; "
             "not a relook")


# --------------------------------------------------------------------------------------------- pure functions
def pr_curve(y, s):
    """Full PR curve at sklearn tie blocks (predicted positive iff score >= threshold), descending threshold.
    Returns a DataFrame with threshold, tp, fp, recall, precision and inc (= d recall x precision, the step-integral
    term whose sum is AP, identical to Ranked.ap / sklearn average_precision_score)."""
    from testeval.metrics import Ranked
    r = Ranked(y, s)
    tp, fp = r._cum(None)
    P = tp[-1]
    rec = tp / P
    prec = tp / np.maximum(tp + fp, 1e-300)
    inc = np.diff(np.r_[0.0, rec]) * prec
    return pd.DataFrame(dict(threshold=r.block_score, tp=tp, fp=fp, recall=rec, precision=prec, inc=inc))


def downsample(curve: pd.DataFrame, max_points: int = MAX_POINTS) -> pd.DataFrame:
    """Keep <= max_points points: the first tie block, the first block reaching each of (max_points - 2) evenly
    spaced recall levels, and the last block. ap_increment = sum of inc over the blocks since the previous kept point
    (exact; sum(ap_increment) = AP)."""
    n = len(curve)
    rec = curve.recall.values
    grid = np.linspace(0.0, 1.0, max_points - 2)
    keep = np.unique(np.r_[0, np.minimum(np.searchsorted(rec, grid - 1e-12, side="left"), n - 1), n - 1])
    starts = np.r_[0, keep[:-1] + 1]
    out = curve.iloc[keep].drop(columns="inc").reset_index(drop=True)
    out["ap_increment"] = np.add.reduceat(curve.inc.values, starts)
    out.insert(1, "n_pred_pos", out.tp + out.fp)
    assert len(out) <= max_points
    return out


def naive_step_ap(rec, prec) -> float:
    """AP from the exported points alone (rectangles between kept points) - an approximation, reported only."""
    return float(np.sum(np.diff(np.r_[0.0, rec]) * prec))


def per_cluster(y, scores: dict, ca, cb, w, cl_table) -> pd.DataFrame:
    """Per-cluster AP exactly as Part 1 macro_cluster_ap (same groups, order, eligibility and weighted_ap), plus the
    cluster's size measures. Rows: eligible clusters (>= 1 positive and >= 1 negative among touching pairs)."""
    from testeval.part1 import metrics as p1m
    y = np.asarray(y).astype(bool)
    nprot = cl_table.groupby("cluster").size()
    rows = []
    for c, ix in p1m._groups_touching(np.asarray(ca), np.asarray(cb)).items():
        if y[ix].any() and (~y[ix]).any():
            d = dict(cluster=int(c), n_proteins=int(nprot[c]), n_pos_pairs=int(y[ix].sum()),
                     n_neg_pairs=int((~y[ix]).sum()), ix=ix)
            for m, s in scores.items():
                d[m] = p1m.weighted_ap(y[ix], s[ix], w[ix])
            rows.append(d)
    return pd.DataFrame(rows)


def bin_label(v, bins):
    for name, lo, hi in bins:
        if lo <= v <= hi:
            return name
    raise ValueError(f"value {v} outside every bin")


def ap_by_size(universe, y, scores: dict, w, pc: pd.DataFrame, n_boot=N_BOOT, seed=BOOT_SEED) -> list:
    from testeval.part1 import metrics as p1m
    methods = list(scores)
    rows = []
    for by, bins in BINS.items():
        lab = pc[by].map(lambda v: bin_label(v, bins))
        for name, _, _ in [("all", 0, 0)] + bins:
            sub = pc if name == "all" else pc[lab == name]
            if not len(sub):
                for m in methods:
                    rows.append(dict(universe=universe, bin_by=by, bin=name, n_clusters=0, method=m))
                continue
            ix = np.unique(np.concatenate(sub.ix.values))
            A = sub[methods].values                                   # clusters x methods
            rng = np.random.default_rng(seed)
            draws = rng.integers(0, len(sub), (n_boot, len(sub)))
            bm = A[draws].mean(axis=1)                                # n_boot x methods
            g = methods.index(REF_BASELINE)
            for j, m in enumerate(methods):
                d = bm[:, j] - bm[:, g]
                rows.append(dict(
                    universe=universe, bin_by=by, bin=name, n_clusters=int(len(sub)),
                    n_proteins_median=float(sub.n_proteins.median()), n_proteins_min=int(sub.n_proteins.min()),
                    n_proteins_max=int(sub.n_proteins.max()), n_pos_pairs_median=float(sub.n_pos_pairs.median()),
                    pooled_n_pos=int(y[ix].sum()), pooled_n_neg=int((~y[ix]).sum()), method=m,
                    macro_ap=float(A[:, j].mean()),
                    macro_ap_ci_low=float(np.quantile(bm[:, j], 0.025)),
                    macro_ap_ci_high=float(np.quantile(bm[:, j], 0.975)),
                    macro_ap_minus_gls=float(A[:, j].mean() - A[:, g].mean()),
                    macro_ap_minus_gls_ci_low=float(np.quantile(d, 0.025)) if m != REF_BASELINE else 0.0,
                    macro_ap_minus_gls_ci_high=float(np.quantile(d, 0.975)) if m != REF_BASELINE else 0.0,
                    pooled_ap=p1m.weighted_ap(y[ix], scores[m][ix], w[ix])))
    return rows


# --------------------------------------------------------------------------------------------- the export
def sha256_bytes(b: bytes) -> str:
    import hashlib
    return hashlib.sha256(b).hexdigest()


def _git(*a) -> str:
    return subprocess.run(["git", "-C", str(ROOT), *a], capture_output=True, text=True).stdout.strip()


def check_scores(E, paths, tr: dict) -> dict:
    """Every cached score TSV must hash to its score record (written at scoring) and to the look's input record."""
    from testeval.protocol import sha256_file
    out = {}
    for m in E.METHODS:
        f = paths.SCORES / "test" / f"{m}.tsv"
        h = sha256_file(f)
        rec = json.loads((paths.RESULTS / "score_records" / "test" / f"{m}.run.json").read_text())["output"]["sha256"]
        look = tr["inputs"]["scores"][m]["sha256"]
        if not (h == rec == look):
            raise SystemExit(f"{m}: score TSV {h} != score record {rec} / look input {look}: refused")
        out[m] = h
    return out


def main():
    if os.environ.get("PYTHONHASHSEED") != "0":           # Part 1 Louvain tie-breaking, as in scripts/21
        os.execve(sys.executable, [sys.executable, *sys.argv], dict(os.environ, PYTHONHASHSEED="0"))
    ap = argparse.ArgumentParser()
    ap.add_argument("--rerun-reason", default=None, help="only if an earlier export record exists (logged)")
    a = ap.parse_args()
    t_start = time.perf_counter()

    from testeval import evaluate as E
    from testeval import part1, paths
    from testeval.protocol import append_jsonl, read_jsonl, sha256_file

    out_dir = paths.RESULTS / "export"
    tr_path = paths.RESULTS / "test_results.json"
    tr_sha = sha256_file(tr_path)
    tr = json.loads(tr_path.read_text())
    prior = [r for r in read_jsonl(paths.LOOKS) if r.get("type") == EXPORT_TYPE]
    if prior and not a.rerun_reason:
        raise SystemExit(f"a descriptive export is already recorded ({prior[-1].get('time')}); "
                         "a repeat needs --rerun-reason")

    # ---- inputs, exactly as mode_test loads them
    score_sha = check_scores(E, paths, tr)
    paths.verify_train_files()
    p1 = part1.assert_locked()
    pos, neg = paths.test_files()
    if {"test_pos": sha256_file(pos), "test_neg": sha256_file(neg)} != {k: tr["inputs"][k] for k in ("test_pos", "test_neg")}:
        raise SystemExit("label files differ from the look's inputs: refused")
    t0 = time.perf_counter()
    es = E.load_set("hu.MAP test", pos, neg, paths.SCORES / "test", train_paths=(paths.TRAIN_POS, paths.TRAIN_NEG))
    if json.loads(json.dumps(es.audit, default=int)) != tr["audit"]:
        raise SystemExit("audit counts differ from the look: refused")
    # ---- universes and clusters, exactly as run_protocol builds them
    cl = E.annotate_clusters(es)
    s1 = es.subset(es.mapped_all, es.name + " mapped-only")
    cl1 = E.annotate_clusters(s1)
    sets = {"itt": (es, cl), "mapped_only": (s1, cl1)}
    t_load = time.perf_counter() - t0
    for u, (_, c) in sets.items():
        if json.loads(json.dumps(c["stats"], default=int)) != tr["results"][u]["cluster_stats"]:
            raise SystemExit(f"{u}: cluster statistics differ from the look: refused")

    # ---- PR curves + AP check
    curves, checks = [], {}
    for u, (s, _) in sets.items():
        for m in E.METHODS:
            full = pr_curve(s.y, s.scores[m])
            ds = downsample(full)
            ap_look = tr["results"][u]["metrics"][m]["ap"]
            ap_full = float(full.inc.sum())
            ap_export = float(ds.ap_increment.sum())
            naive = naive_step_ap(ds.recall.values, ds.precision.values)
            assert abs(ap_full - ap_look) <= TOL, (u, m, ap_full, ap_look)
            assert abs(ap_export - ap_look) <= TOL, (u, m, ap_export, ap_look)
            checks.setdefault(u, {})[m] = dict(ap_look=ap_look, ap_export=ap_export,
                                               abs_diff=abs(ap_export - ap_look), n_points=len(ds),
                                               n_tie_blocks=len(full), naive_step_ap_from_points=naive,
                                               naive_abs_diff=abs(naive - ap_look))
            ds.insert(0, "method", m)
            ds.insert(0, "universe", u)
            curves.append(ds)

    # ---- AP by cluster size + macro AP check
    size_rows = []
    for u, (s, c) in sets.items():
        pc = per_cluster(s.y, s.scores, c["ca"], c["cb"], s.w, c["table"])
        for m in E.METHODS:
            ref = tr["results"][u]["metrics"][m]
            macro = float(pc[m].mean())
            assert len(pc) == ref["n_clusters_macro"], (u, m, len(pc), ref["n_clusters_macro"])
            assert abs(macro - ref["macro_cluster_ap"]) <= TOL, (u, m, macro, ref["macro_cluster_ap"])
            checks[u][m].update(macro_look=ref["macro_cluster_ap"], macro_export=macro,
                                macro_abs_diff=abs(macro - ref["macro_cluster_ap"]))
        rows = ap_by_size(u, s.y, s.scores, s.w, pc)
        for by in BINS:                                       # bins recombine to the overall macro AP
            r = pd.DataFrame([x for x in rows if x["bin_by"] == by and x["bin"] != "all" and x["n_clusters"] > 0])
            for m in E.METHODS:
                rm = r[r.method == m]
                assert abs((rm.macro_ap * rm.n_clusters).sum() / rm.n_clusters.sum()
                           - tr["results"][u]["metrics"][m]["macro_cluster_ap"]) <= TOL
        size_rows += rows

    # ---- write
    out_dir.mkdir(parents=True, exist_ok=True)
    buf = io.StringIO()
    pd.concat(curves, ignore_index=True).to_csv(buf, sep="\t", index=False, float_format="%.17g")
    with gzip.GzipFile(out_dir / "pr_curves.tsv.gz", "wb", mtime=0) as fh:
        fh.write(buf.getvalue().encode())
    sz = pd.DataFrame(size_rows)
    for k in ("n_proteins_min", "n_proteins_max", "pooled_n_pos", "pooled_n_neg"):
        sz[k] = sz[k].astype("Int64")                         # empty bins: NA, not 3.0-style floats
    sz.to_csv(out_dir / "ap_by_size.tsv", sep="\t", index=False, float_format="%.17g")
    runtime = round(time.perf_counter() - t_start, 1)
    head = _git("rev-parse", "HEAD")
    dirty = [l for l in _git("status", "--porcelain").splitlines() if l.strip()]
    (out_dir / "README.md").write_text(readme(checks, runtime, head))
    if sha256_file(tr_path) != tr_sha:
        raise SystemExit("test_results.json changed during the export (must never happen)")
    outputs = {f: sha256_file(out_dir / f) for f in ("pr_curves.tsv.gz", "ap_by_size.tsv", "README.md")}
    rec = dict(type=EXPORT_TYPE, non_confirmatory=True, statement=STATEMENT,
               commit=head, worktree_dirty=dirty, script="scripts/23_export_test_curves.py",
               script_sha256=sha256_file(Path(__file__)), rerun_reason=a.rerun_reason,
               reads=dict(labels="paths.test_files() via evaluate.load_set (as mode_test)", scores=score_sha),
               test_results_sha256=tr_sha, test_results_unchanged=True, part1_commit=p1["commit"],
               outputs={f"results/test_eval/export/{k}": v for k, v in outputs.items()},
               checks=dict(max_ap_abs_diff=max(v["abs_diff"] for d in checks.values() for v in d.values()),
                           max_macro_abs_diff=max(v["macro_abs_diff"] for d in checks.values() for v in d.values()),
                           tolerance=TOL),
               runtime_s=dict(load_and_clusters=round(t_load, 1), total=runtime))
    append_jsonl(paths.LOOKS, rec)
    print(json.dumps(dict(rec, checks=checks), indent=1, default=str))


def readme(checks, runtime, head) -> str:
    rows = "\n".join(f"| {u} | {m} | {v['ap_look']:.6f} | {v['abs_diff']:.1e} | {v['macro_look']:.6f} | "
                     f"{v['macro_abs_diff']:.1e} | {v['n_points']} | {v['naive_abs_diff']:.1e} |"
                     for u, d in checks.items() for m, v in d.items())
    bins = "; ".join(f"{by}: " + ", ".join(b[0] for b in bs) for by, bs in BINS.items())
    return f"""# Test look: descriptive export (non-confirmatory)

**This is a {STATEMENT}.** It re-derives two figure inputs that `scripts/21_evaluate_test.py --mode test --run`
computed internally but did not export. `test_results.json`, the protocol and the harness are unchanged; nothing
here is a pre-specified statistic of docs/TEST_PROTOCOL.md, and no claim may be read as a new test.

Produced by `scripts/23_export_test_curves.py` at commit `{head}` (the script and these files are committed together
in the next commit); logged in `results/test_eval/looks.jsonl` as a `"type": "{EXPORT_TYPE}"` record with
`non_confirmatory: true` (no `status` field, so the ledger's look logic, which counts `started` / `done` records,
is unaffected). Runtime {runtime} s.

## Inputs and definitions (the harness's own)
- Labels: `testeval.paths.test_files()` through `evaluate.load_set` with the train files, exactly as `mode_test`;
  label-file SHA-256 and audit counts asserted equal to `test_results.json`.
- Scores: `data_cache/test_eval/scores/test/<method>.tsv`; SHA-256 asserted equal to
  `results/test_eval/score_records/test/<method>.run.json` and to `test_results.json` `inputs.scores`.
- Universes: ITT = all analysed pairs; S1 = pairs mapped under every method (`EvalSet.mapped_all`).
- Clusters: `evaluate.annotate_clusters` (Part 1 `build_clusters`, PYTHONHASHSEED=0) on each universe separately,
  as `run_protocol`; cluster statistics asserted equal to the look's.

## `pr_curves.tsv.gz`
One row per exported point; columns `universe, method, threshold, n_pred_pos, tp, fp, recall, precision,
ap_increment`. Points are sklearn tie blocks (predicted positive iff score >= threshold), descending threshold. At most
{MAX_POINTS} points per method and universe: the first block, the first block reaching each of {MAX_POINTS - 2} evenly
spaced recall levels, and the last block (recall 1, precision = prevalence). Thresholds are raw method scores
(unmapped pairs carry the method's floor score in ITT).

**Integrating AP.** `ap_increment` is the exact step-integral mass (delta recall x precision over all tie blocks)
accrued since the previous exported point, so **AP = sum(ap_increment)** per method x universe. That equals
`test_results.json` AP to <= {TOL:g} (asserted; table below). Rectangles over the exported points alone
(`sum(diff(recall) * precision)`) only approximate AP, because the down-sampling drops tie blocks; that error is
listed for information (`naive |diff|`).

## `ap_by_size.tsv`
Per universe x size measure x bin x method:
- `macro_ap`: mean per-cluster AP of the bin's clusters. Per-cluster AP is Part 1 `macro_cluster_ap`'s own term:
  AP over the pairs touching the cluster (>= 1 endpoint in it), for clusters with >= 1 positive and >= 1 negative
  among them. Bin `all` = the look's `macro_cluster_ap` (asserted), and the cluster-count-weighted mean of the bins
  recombines to it (asserted).
- `pooled_ap`: AP over the union of pairs touching any cluster in the bin (a between-cluster pair can enter two bins).
- Size measures (bins fixed before the export was run): {bins}. `n_pos_pairs` = positive pairs touching the cluster
  (its within-cluster positives; no component was Louvain-split on the test set).
- CIs (`macro_ap_ci_*`, `macro_ap_minus_gls_ci_*`): descriptive percentile intervals from resampling the bin's
  clusters with replacement (B = {N_BOOT}, seed {BOOT_SEED}), per-cluster APs held fixed, paired across methods. They
  are **not** the protocol's paired cluster bootstrap (§6.1) and ignore between-cluster dependence; no pooled-AP CI.

## Reproduction of the look
| universe | method | AP (look) | AP |diff| | macro AP (look) | macro |diff| | points | naive |diff| |
|---|---|---|---|---|---|---|---|
{rows}

## Note on the frozen marker
`tests/test_testeval.py` (a hashed harness file in `docs/TEST_PROTOCOL.FROZEN`) gained this script on its grep-guard
allowlist in the same commit. It already differed from the marker after the look (8fc91f6, 8f63b6e), so
`require_frozen()` was already refusing, which only blocks a further look. No file of the look's code path
(`src/testeval`, scripts 19-22) changed, and every number in `test_results.json` is unchanged.
"""


if __name__ == "__main__":
    main()
