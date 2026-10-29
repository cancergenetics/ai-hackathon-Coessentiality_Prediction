# Test look: descriptive export (non-confirmatory)

**This is a descriptive presentation of the single confirmatory look; no decision or verdict depends on it; not a relook.** It re-derives two figure inputs that `scripts/21_evaluate_test.py --mode test --run`
computed internally but did not export. `test_results.json`, the protocol and the harness are unchanged; nothing
here is a pre-specified statistic of docs/TEST_PROTOCOL.md, and no claim may be read as a new test.

Produced by `scripts/23_export_test_curves.py` at commit `bff303087bf50840bb56b39e8f6bd0dc30447982` (the script and these files are committed together
in the next commit); logged in `results/test_eval/looks.jsonl` as a `"type": "descriptive_export"` record with
`non_confirmatory: true` (no `status` field, so the ledger's look logic, which counts `started` / `done` records,
is unaffected). Runtime 17.1 s.

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
2000 points per method and universe: the first block, the first block reaching each of 1998 evenly
spaced recall levels, and the last block (recall 1, precision = prevalence). Thresholds are raw method scores
(unmapped pairs carry the method's floor score in ITT).

**Integrating AP.** `ap_increment` is the exact step-integral mass (delta recall x precision over all tie blocks)
accrued since the previous exported point, so **AP = sum(ap_increment)** per method x universe. That equals
`test_results.json` AP to <= 1e-09 (asserted; table below). Rectangles over the exported points alone
(`sum(diff(recall) * precision)`) only approximate AP, because the down-sampling drops tie blocks; that error is
listed for information (`naive |diff|`).

## `ap_by_size.tsv`
Per universe x size measure x bin x method:
- `macro_ap`: mean per-cluster AP of the bin's clusters. Per-cluster AP is Part 1 `macro_cluster_ap`'s own term:
  AP over the pairs touching the cluster (>= 1 endpoint in it), for clusters with >= 1 positive and >= 1 negative
  among them. Bin `all` = the look's `macro_cluster_ap` (asserted), and the cluster-count-weighted mean of the bins
  recombines to it (asserted).
- `pooled_ap`: AP over the union of pairs touching any cluster in the bin (a between-cluster pair can enter two bins).
- Size measures (bins fixed before the export was run): n_pos_pairs: 1, 2-5, 6-20, 21-100, 101+; n_proteins: 2, 3-4, 5-8, 9-16, 17+. `n_pos_pairs` = positive pairs touching the cluster
  (its within-cluster positives; no component was Louvain-split on the test set).
- CIs (`macro_ap_ci_*`, `macro_ap_minus_gls_ci_*`): descriptive percentile intervals from resampling the bin's
  clusters with replacement (B = 2000, seed 0), per-cluster APs held fixed, paired across methods. They
  are **not** the protocol's paired cluster bootstrap (§6.1) and ignore between-cluster dependence; no pooled-AP CI.

## Reproduction of the look
| universe | method | AP (look) | AP |diff| | macro AP (look) | macro |diff| | points | naive |diff| |
|---|---|---|---|---|---|---|---|
| itt | zavg_BEF | 0.332901 | 0.0e+00 | 0.296778 | 0.0e+00 | 1906 | 1.5e-04 |
| itt | B | 0.310456 | 5.6e-17 | 0.254717 | 0.0e+00 | 1906 | 1.5e-04 |
| itt | gls_signed | 0.211548 | 8.3e-17 | 0.291753 | 0.0e+00 | 1906 | 1.2e-04 |
| itt | spearman_signed | 0.204309 | 2.8e-17 | 0.175750 | 0.0e+00 | 1906 | 1.1e-04 |
| itt | pearson_signed | 0.196604 | 2.8e-17 | 0.173491 | 0.0e+00 | 1906 | 1.1e-04 |
| itt | spec_k20_tinf_std | 0.239019 | 2.8e-17 | 0.264966 | 0.0e+00 | 1906 | 1.2e-04 |
| mapped_only | zavg_BEF | 0.348758 | 5.6e-17 | 0.309849 | 0.0e+00 | 1999 | 1.1e-04 |
| mapped_only | B | 0.325225 | 0.0e+00 | 0.266275 | 0.0e+00 | 1999 | 8.1e-05 |
| mapped_only | gls_signed | 0.221527 | 0.0e+00 | 0.304287 | 0.0e+00 | 1999 | 1.2e-04 |
| mapped_only | spearman_signed | 0.213938 | 0.0e+00 | 0.184837 | 0.0e+00 | 1998 | 1.0e-04 |
| mapped_only | pearson_signed | 0.205859 | 0.0e+00 | 0.181466 | 0.0e+00 | 1999 | 1.2e-04 |
| mapped_only | spec_k20_tinf_std | 0.250329 | 5.6e-17 | 0.275802 | 0.0e+00 | 1998 | 1.2e-04 |

## Note on the frozen marker
`tests/test_testeval.py` (a hashed harness file in `docs/TEST_PROTOCOL.FROZEN`) gained this script on its grep-guard
allowlist in the same commit. It already differed from the marker after the look (8fc91f6, 8f63b6e), so
`require_frozen()` was already refusing, which only blocks a further look. No file of the look's code path
(`src/testeval`, scripts 19-22) changed, and every number in `test_results.json` is unchanged.
