# Dry run of the test-evaluation harness on development data

**Data used.** No hu.MAP test pair was read, scored or counted. Inputs:
- the hu.MAP **train** pairs (read-only Part 1 folder);
- Part 1's dev/IFT protein assignment (SHA-256 checked against `external/part1_locked/results/splits/protein_assignment.sha256`);
- the locked Part 1 predictor (commit c7f95ee).

The protocol is in `docs/TEST_PROTOCOL.md`; §10 there summarises these results.

## Commands
```bash
scripts/20_score_test_pairs.sh dryrun                               # 333 s: dev (344,055 pairs) + IFT (43,994) x 6 methods
scripts/20_score_test_pairs.sh train                                # 183 s: all train pairs (634,779) x 6 methods
.venv/bin/python scripts/21_evaluate_test.py --mode fit             # 6 s  -> ../test_eval/fitted_rules.json
.venv/bin/python scripts/21_evaluate_test.py --mode dryrun          # parts A, B, C, I (run here as --parts AI / B / C)
```

## Files

| file | content |
|---|---|
| `scoring_log.txt`, `score_records/` | the predictor's run records: clone commit, input SHA-256s, null-moment checks (difference 0), mapping counts, runtimes |
| `A_ift_reproduction.json` | **A.** Part 1's IFT look 1 recomputed with Part 1's definitions from freshly scored IFT pairs. Each field is compared with `look01_2026-10-27.json` and the `2-posthoc` look |
| `B_ift_protocol.json` | **B.** The full protocol pipeline (ITT plus mapped-only, rules, calibration, bootstraps, verdict) with the IFT as the evaluation set and Part 1 dev pairs as the fit set; HT design weights |
| `C_train_as_test.json`, `C_log.txt` | **C.** The exact test-mode code path (weights 1, clusters rebuilt from the file's positives) on all train pairs. It is in-sample, so it checks mechanics, runtime and the expected CI width |
| `summary.json` | headline numbers, the batch-invariance check (**I**) and runtimes per step |
| `figures/` | reliability diagrams (isotonic and Platt, 16 log-spaced bins) |

## Results

### A. Reproduction of IFT look 1

| method | harness HT-AP | Part 1 HT-AP |
|---|---|---|
| zavg_BEF | 0.345139 | 0.345139 |
| GLS | 0.242318 | 0.242318 |
| Spearman | 0.146273 | 0.146273 |
| Pearson | 0.134064 | 0.134064 |
| B | 0.317398 | 0.317398 |
| PC-removed Pearson | 0.266779 | 0.266779 |

- **Point metrics:** 29 metric-suite fields per method agree to a maximum absolute difference of 1.5e-8 (B AUROC, tie rounding).
- **Bootstrap:** the 1,000-resample cluster-bootstrap CIs agree to 9e-16.
- **Counts:** they match Part 1: 847 / 40,786 mapped and 9 / 2,352 unmapped, f = 0.26481.

### B. Protocol on the IFT

ITT, cluster bootstrap with B = 2,000:

| comparison | ΔAP [95% CI] |
|---|---|
| zavg_BEF vs GLS | +0.099 [+0.034, +0.201] |
| zavg_BEF vs Spearman | +0.191 [+0.094, +0.302] |
| zavg_BEF vs Pearson | +0.203 [+0.098, +0.324] |

- **Verdict:** zavg_BEF clearly beats all three baselines. B fails against GLS on macro-AP non-inferiority.

### C. All train pairs

| | zavg_BEF | GLS | Spearman | Pearson |
|---|---|---|---|---|
| AP | 0.389 | 0.245 | 0.214 | 0.195 |

- **ΔAP vs GLS:** +0.144 [+0.078, +0.217].
- **Clusters:** Kish n_eff = 25.8.

### I. Batch invariance

Scores are identical (difference 0) whether a pair is scored in the IFT file or in the all-train file.

## Runtime (8-core arm64 laptop)

| step | time |
|---|---|
| A | 9 s |
| B | 66 s |
| C, all ITT bootstraps on 635k pairs | 830 s |
| C, mapped-only | 144 s |

The C split of the ITT bootstraps:

| bootstrap | time |
|---|---|
| cluster, all statistics | 252 s |
| protein | 112 s |
| pair | 119 s |
| Monte-Carlo seeds 1–4 | 330 s |
