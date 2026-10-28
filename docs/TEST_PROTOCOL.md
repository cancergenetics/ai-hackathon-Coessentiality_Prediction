# TEST_PROTOCOL.md: pre-registered evaluation of the locked Part 1 method on the withheld hu.MAP test pairs

**Version 1.0 (draft for PI freeze), 2026-10-28.**

**How it is frozen.** The protocol is frozen by committing this file together with `src/testeval/`, `scripts/19_*`, `scripts/20*`, `scripts/21_*`, `tests/test_testeval.py` and `results/test_eval/fitted_rules.json`.

**The single look.** `scripts/21_evaluate_test.py --mode test` refuses to run unless:
- this file and the fitted rules are committed, and
- `src/` and `scripts/` are clean.

It then records the commit hash of this file (`git log -1 -- docs/TEST_PROTOCOL.md`), HEAD and the SHA-256 of this file in `results/test_eval/test_results.json` and in `results/test_eval/looks.jsonl`.

**What may not change.** Nothing in this document may be changed after the look. Any deviation is reported as a deviation (§8.3).

---

## 0. Transition from development to test, and what has been seen

| item | status at freeze |
|---|---|
| Part 1 lock | commit `c7f95ee007576be9311c0a5f08c6ac63ca23a009` (2026-10-28 12:55:43 UTC), cloned read-only to `external/part1_locked/` |
| Test files | `datasets/huMAP/Test_positive_pairs.txt` (SHA-256 `c4e3e68e…2d96`), `Test_negative_pairs.txt` (`c4a3e546…e9de`); they arrived 13:55 UTC, about 1 h after the Part 1 lock (DATA.md) |
| Test content seen before freeze | the PI read the line counts (3,312 / 597,426) and the first 3 lines of each file. The harness author printed the first ~200 bytes (6 lines) of the positive file as an octal dump to check the delimiter and line endings. **Not computed:** scores, predictions, metrics, de-duplication or overlap counts, protein lists, mapping rates. |
| Development data used here | hu.MAP **train** pairs (Part 1 folder, read-only), Part 1's dev/IFT split (`data_cache/splits/protein_assignment.csv`, SHA-256 checked against Part 1's committed hash), DepMap 26Q3 CRISPRGeneEffect, HGNC |
| Part 1 never accessed test labels | verified separately in `docs/part1_audit.md` (pending; recommended before the look, §9) |

Everything below was fixed from Part 1's documents and code and from development data only.

## 1. Locked method and inputs

**Predictor**
- Run as `external/part1_locked/scripts/predict_pairs.py`.
- It reads the raw `datasets/DepMAP/CRISPRGeneEffect.csv` (no cache), SHA-256-checks the parsed matrix against the Part 1 layer manifests, and asserts the null moments for B/E/F equal the manifests.
- The ID mapping is Part 1's (`src.data.build_uniprot_map`: HGNC Approved, UniProt → Entrez → DepMap column).

**Locking checks**
- `scripts/20_score_test_pairs.sh` checks that the clone is at `c7f95ee` and has no modified tracked files.
- Every run record (`results/test_eval/score_records/test/*.run.json`) stores:
  - the clone commit;
  - input SHA-256s;
  - null-moment checks;
  - mapping counts and runtimes.

**Failure rule.** If the locked predictor fails, the look is not run and the failure is reported. No Part 1 code is modified (decision D0).

**Labels.** No label reaches the predictor. It receives one ID-only file per evaluated set (§3.2).

## 2. Methods scored

All methods are scored on exactly the same pairs (§3).

| name in results | role | command (locked predictor unless noted) | Part 1 source |
|---|---|---|---|
| `zavg_BEF` | **primary (final method)** | `--method zavg:B,E,F` | D3/D4, REPORT §3 |
| `B` | pre-declared secondary candidate (D3 fixed sequence, second) | `--method single` | `results/phaseB_v2/decision.json` |
| `gls_signed` | baseline: GLS (Wainberg 2021), dev-selected variant: signed t, no olfactory-receptor PC removal, exact marginal-Σ NaN protocol | `--method gls` | `results/baselines_v2/best_variants.json` |
| `spearman_signed` | baseline: Spearman ρ, pairwise complete, signed | `--method spearman` | same |
| `pearson_signed` | baseline: Pearson r, pairwise complete, signed | `--method pearson` | same |
| `spec_k20_tinf_std` | **post hoc reference, secondary, not a candidate**: PC-removed Pearson (gene-standardised, top 20 PCs removed, no whitening) | `scripts/20b_score_posthoc_reference.py` | Part 1 review B1, IFT look `2-posthoc`, `final_eval.py --posthoc-reference` |

**Baselines.** These are exactly the variants of Part 1's IFT look 1; the selection rule was max dev HT-AP, with the same tuning budget for every method.

**Post hoc reference.**
- `predict_pairs.py` does not expose it. `src/testeval/posthoc.py` computes it with Part 1's own `LayerEngine`, ID mapping and unmapped rule.
- It asserts that the frozen config equals manifest column 64 (screen_v1) and that the null moments equal the manifest's (observed difference 0).
- The dry run reproduces its IFT HT-AP (0.26678) to 1e-16.
- It enters no verdict.

## 3. Pair set

### 3.1 Rules

These rules are fixed a priori and applied to every method identically (`src/testeval/pairs.py`).

| rule | decision | justification |
|---|---|---|
| R1 unordered pair, de-duplication | canonical `(a, b)` with `a < b` on the UniProt accession string (Part 1 `data.load_pairs` rule); a pair listed several times (same or reversed orientation, within or across lines of a file) counts once | the task is an unordered pair (Part 1 DESIGN §1); duplicates would give a pair extra weight in AP |
| R2 self pairs (`a == b`) | **excluded** from every analysis; counted | not a pair of two proteins; outside the predictor's domain (Part 1's loader rejects them; the predictor marks them unmapped) |
| R3 pair in both the positive and the negative file | **excluded** from every analysis; counted | label undefined; any assignment would be arbitrary |
| R4 pair also in a train file (either label) | **kept** in the primary analysis; flagged; dropped in sensitivity set S2 | the test files as delivered define the test population. All scored methods are label-free, so train labels cannot leak into their scores; only the train-fitted thresholds and calibrators could profit |
| R5 unmapped pair (a protein without a unique DepMap gene, both proteins on the same gene, or an undefined score) | **kept** in the primary analysis at the locked predictor's floor score: min(mapped scores in the scored file) − 1 (Part 1 DESIGN §3.10). This is intention-to-score (ITT): every labelled pair is scored. The mapped-only analysis is secondary (S1) | the method's real-world coverage is part of its performance. All unmapped pairs form one tied block at the bottom for every method, so the rule treats all methods identically. Part 1 reported this variant on the IFT (0.332 vs 0.345 mapped-only) |

### 3.2 Scoring input

One ID-only file per set: the unique canonical non-self pairs of positives ∪ negatives, sorted, no labels (`scripts/19_make_pair_files.py`). Positives and negatives therefore share one floor.

**Primary pair set (ITT).** All canonical pairs remaining after R2 and R3.

**Pre-specified secondary sets.** Every metric is reported for every method on each:
- **S1 mapped-only:** pairs mapped under every method. This is Part 1's IFT definition ("scorable").
- **S2 no train overlap:** drop pairs that occur in either train file.
- **S3 unseen proteins:** pairs whose two proteins both occur in no train pair. This is protein-disjoint, as in Part 1's development evaluation.

**Counts reported.** Lines, unique pairs, duplicates, reversed duplicates, self pairs, conflicts, overlap with train by label, proteins also in train, and unmapped pairs by reason and by label.

**Expected prevalence.** 3,312 / 600,738 = 0.551% before de-duplication, against 0.522% for train.

## 4. Primary metric: AUPRC (average precision) at natural prevalence

**Definition.** AP is computed on the ITT set with every pair weighted 1:

$\text{AP} = \sum_n (R_n - R_{n-1}) P_n$

- The sum runs over distinct score thresholds; it is non-interpolated.
- Ties are one block (sklearn `average_precision_score`).
- It is implemented by calling Part 1's `src.metrics.metric_suite` / `weighted_ap` with w = 1. The harness's fast O(n) re-implementation is tested equal to sklearn and Part 1.

**Relation to Part 1's HT-weighted AP.**
- Part 1's development folds and IFT were *cluster subsamples* of the train pairs.
- A pair within a cluster enters a subsample with probability f. A pair between two clusters enters with probability ≈ f², so between-cluster pairs, nearly all of them negatives, were under-represented by 1/f.
- Horvitz–Thompson weights restored the natural pair mix of the full labelled universe: 1 within a cluster, 1/f between clusters (f = 0.265 on the IFT). This restored the natural prevalence (0.53%) and the natural share of hard negatives.
- The test file is not a subsample: it is the full labelled test population as delivered. Every pair therefore has inclusion probability 1 and HT weight 1, so HT-AP ≡ plain AP at the file's natural prevalence.
- The definitions coincide; only the weights differ (all 1). The dry run confirms the HT machinery: with Part 1's weights it reproduces the IFT HT-AP exactly (§10).

**Reference values (not targets).**

| method | dev HT-AP | IFT HT-AP, mapped | IFT HT-AP, ITT |
|---|---|---|---|
| zavg_BEF | 0.470 | 0.345 | 0.332 |
| GLS | 0.312 | 0.242 | 0.233 |

## 5. Secondary metrics

All are computed on the ITT set and on S1–S3, for every method.

### 5.1 Ranking metrics
- **AUROC:** weighted, ties count ½; equals sklearn `roc_auc_score`.
- **Precision at recall 0.10 / 0.25 / 0.50:** Part 1 `weighted_precision_at_recall`, using the first tie block that reaches the recall level.
- **AP with only hard negatives:** negatives in the same positive-graph component (Part 1).

### 5.2 Cluster-structured metrics

**Test clusters** are defined exactly as in Part 1 §3.1, applied to the test positives (Part 1 `src.clusters.build_clusters`, imported unchanged):
- universe: all proteins of the ITT pair set, including unmapped ones;
- graph: the ITT positive pairs;
- clusters: connected components;
- every component with more than 40 proteins is split by `networkx.community.louvain_communities(seed=0)` (resolution 1, integer node labels sorted by accession), run under `PYTHONHASHSEED=0`, which the script enforces;
- proteins without a positive pair are singleton clusters;
- cluster ids are deterministic.

Part 1 used networkx Louvain, not Leiden, and the protocol keeps Louvain.

The cluster-structured metrics on these clusters:
- **Macro per-cluster AP (Part 1 `macro_cluster_ap`):** for each cluster whose touching pairs (≥ 1 endpoint in the cluster) include at least one positive and one negative, compute the AP over those pairs, then take the unweighted mean over clusters.
- **AP without the largest component:** the largest component is the one with most proteins (first by accession on ties). Drop every pair with an endpoint in it. This is Part 1's non-inferiority metric. On the IFT, Part 1's "largest component" was the train giant component; on the test it is the largest *test* component (dry-run note §10).
- **AP without the 10 largest clusters:** the 10 clusters with most within-cluster positives (Part 1 `21_phaseB_v2.top_clusters`); drop every pair touching them.

### 5.3 Essentiality strata
- **Common-essential gene:** the fraction of CRISPRGeneEffect models with effect < −0.5 is ≥ 0.5 (Part 1 `universe.common_essential`, all models).
- **Strata:** AP for pairs with both genes non-essential, mixed, and both essential.
- **Pair set:** mapped pairs only (S1), because unmapped proteins have no gene.

### 5.4 Precision, recall and F1 at pre-specified decision rules

These rules are fixed from development data only and stored in `results/test_eval/fitted_rules.json`. Predicted positive means score ≥ t.

| rule | definition | applies to | status |
|---|---|---|---|
| **(b) F1-max on train** | t = the score threshold that maximises F1 over **all hu.MAP train pairs** at their natural prevalence (unweighted; unmapped train pairs at the train file's floor); ties go to the highest such t | every method | **primary operating point** |
| (c) top-k | K = round(π_train × N), with π_train = 3,314 / 634,779 = 0.522% and N = number of pairs in the analysed set; t = the K-th ranked score (whole tie block) | every method | secondary: an equal prediction budget for all methods; it uses only the count N of the test set, never its labels |
| (a) z ≥ 3 | t = 3.0 on the zavg_BEF score | zavg_BEF only | secondary: the cutoff of Part 1's released genome-wide `top_pairs.parquet` (results/predictions/README.md). It is not label-tuned, so it shows what a user of the released table gets |

**Fitted thresholds** (from `fitted_rules.json`). Each threshold is followed by the train-set F1, which is in-sample.

| method | t, rule (b) | train F1 |
|---|---|---|
| zavg_BEF | 6.567 | 0.478 |
| B | 11.78 | 0.463 |
| GLS | 3.253 | 0.339 |
| Spearman | 0.197 | 0.297 |
| Pearson | 0.202 | 0.277 |
| PC-removed Pearson | 0.135 | 0.397 |

For zavg_BEF, rule (a) on train gives P 0.285 and R 0.496.

**Why rule (b) is primary.**
- It is the only operating point that is defined for every method on an equal footing from labels already used in development.
- All scores are label-free functions of the same DepMap matrix, and train and test prevalences are close (0.52% vs 0.55%).
- Rule (c) compares methods at a fixed budget, without any threshold transfer.
- Part 1 pre-specified no label-tuned threshold. Its only fixed cutoff is z ≥ 3, and Part 1 reported precision at recall levels (kept in §5.1).

### 5.5 Calibration

The scores are not probabilities, so each method is mapped to probabilities by a calibrator fitted on the hu.MAP train pairs only (all 634,779, unweighted, natural prevalence 0.522%). The calibrator is applied unchanged to the test pairs.

**Calibrators**
- **Isotonic (primary):** sklearn `IsotonicRegression(increasing, y ∈ [0, 1], out_of_bounds='clip')`. It is stored as knots in `fitted_rules.json` and applied by linear interpolation, which is identical to sklearn's `predict` (tested).
- **Platt (secondary):** maximum-likelihood logistic regression on the standardised score, with no penalty and no target smoothing.
- **Unmapped pairs:** the floor score is specific to each file, so before calibration the score of every unmapped test pair is set to the train file's floor. Unmapped is the lowest score level in both sets.

**Reported** (`src/testeval/calibration.py`):
- **Brier score;**
- **Brier skill score:** against the constant forecast π_train (primary) and against the observed test prevalence (secondary);
- **calibration-in-the-large:** Σp / Σy;
- **calibration intercept and slope:** logistic regression of y on logit p, with p clipped to [1e-6, 1 − 1e-6];
- **ECE:** **20 equal-mass bins** (primary; bins of equal pair count after a stable sort on p), plus 15 equal-width bins on [0, 1] and 16 log-spaced bins on [10⁻⁴, 1] (secondary);
- **reliability tables:** for both the equal-mass and the log-spaced bins;
- **reliability diagram:** log-spaced bins on log–log axes (`results/test_eval/figures/reliability_test_itt.png`). At ~0.5% prevalence, 19 of 20 equal-mass bins lie below p ≈ 0.003, so only log-spaced bins resolve the high-probability tail.

At a prevalence of about 0.5%, ECE in absolute probability units is small for every method. Brier skill and the reliability diagram are therefore the informative quantities.

**Prevalence shift.** Test prevalence is about 0.551% against 0.522% for train (+5.6% relative, from line counts). A calibrator fitted on train is expected to under-predict slightly: calibration-in-the-large is around 0.95 if discrimination is unchanged. The protocol makes no recalibration or prior-shift correction with test information; the shift is reported and discussed.

## 6. Uncertainty

### 6.1 Primary: paired cluster bootstrap
This is Part 1's scheme (`src/metrics.cluster_bootstrap`), re-expressed as a weight generator with the identical random stream. The harness's replicates are tested bit-equal to Part 1's.
- **Units:** the test clusters of §5.2.
- **Replicate weights:** per replicate, m = bincount(rng.integers(0, C, C)) with `rng = numpy.random.default_rng(0)`. A pair within a cluster gets weight m_c; a pair between clusters gets m_c1 · m_c2.
- **Resamples and seed:** **B = 2,000**, seed `SEEDS[0] = 0`.
- **Paired:** all methods and statistics share the same replicates.
- **Intervals:** 95% percentile. Two-sided p-value = 2 · min(P(Δ ≤ 0), P(Δ ≥ 0)) (Part 1 rule).
- **Statistics bootstrapped:**
  - per method: AP, AUROC, AP without the largest component, F1 / precision / recall at rule (b), and isotonic Brier;
  - paired Δ:
    - zavg_BEF − each baseline (AP, AUROC, AP without the largest component);
    - B − each baseline (AP);
    - B − zavg_BEF (AP);
    - PC-removed Pearson − zavg_BEF (AP).
- **Fixed during resampling:** thresholds and calibrators are fitted on train and not refitted, so their sampling error is not propagated (stated as a limitation).
- **Macro per-cluster AP** is reported as a point estimate, as in Part 1's non-inferiority checks.

### 6.2 Secondary bootstraps (AP and ΔAP only)
- **Protein (node) bootstrap:** proteins are resampled; pair weight m_a · m_b; B = 2,000, seed 0.
- **Stratified pair bootstrap:** Part 1's `_boot_weights`; positives and negatives resampled separately; B = 2,000, seed 0. It ignores dependence between pairs that share proteins, so it is reported as anti-conservative.

### 6.3 Monte-Carlo stability
The cluster bootstrap of ΔAP (zavg_BEF − each baseline) is repeated with seeds 1–4 (B = 2,000 each). The verdict uses seed 0 only. The other seeds show how much the CI endpoints depend on the seed.

## 7. Comparison criterion and multiplicity

This is the test-set version of Part 1's "clearly beats" (DESIGN v3.1 §3.7, `final_eval.decide`). Candidate c clearly beats baseline b iff all of the following hold on the **ITT set**:
1. **Significance:** the 95% paired cluster-bootstrap CI of ΔAP = AP_c − AP_b has lower bound > 0 (seed 0, B = 2,000).
2. **Effect size:** relative gain AP_c / AP_b − 1 ≥ 0.10 (point estimate).
3. **Non-inferiority:** macro per-cluster AP_c ≥ 0.95 × macro AP_b **and** AP-without-largest-component_c ≥ 0.95 × that of b.

**Combination rules**
- **Intersection–union:** c clearly beats the baselines iff (1)–(3) hold against **each** of GLS, Spearman and Pearson. The joint claim "beats all three" is then a test at one-sided level 2.5% with no multiplicity adjustment.
- **Fixed sequence:** zavg_BEF (primary) first, which is confirmatory. B is confirmatory only if zavg_BEF clearly beats all three baselines; otherwise its result is exploratory. zavg_BEF stays the final method whatever the outcome.
- **Holm-adjusted p-values** of the three primary ΔAP comparisons are reported as information.
- **Descriptive only:** everything else (secondary metrics, S1–S3, the reference, other bootstraps) has unadjusted 95% CIs.

**Why criterion (1) is added relative to Part 1's IFT gate.**
- Part 1 dropped the CI requirement on the IFT only because of power: n_eff ≈ 7.7 clusters, pass probability ≤ 0.22 at a 20% true gain.
- The test set is hu.MAP's other half of the complexes, so its size is expected to be similar to the full train set.
- On the full train set at natural prevalence (dry run C, §10: 3,314 positives, 245 clusters, Kish n_eff = 25.8 clusters weighted by within-cluster positives, against 7.8 on the IFT), the ΔAP(zavg_BEF − GLS) cluster-bootstrap CI was [+0.078, +0.217], about half the IFT width. A test set of similar size therefore has adequate power at Part 1's observed effect sizes (+0.10 to +0.19). The test effect will probably be smaller than this in-sample value, which is why criterion (1) can still fail.
- Criteria (2) and (3) are Part 1's IFT gate unchanged.

**Secondary verdict.** The same verdict computed on S1 (mapped-only), reported as secondary.

## 8. Execution

### 8.1 Before the freeze (development data only; already run)
```bash
scripts/20_score_test_pairs.sh train                        # locked predictor on all train pairs, ~3 min
.venv/bin/python scripts/21_evaluate_test.py --mode fit     # -> results/test_eval/fitted_rules.json, ~6 s
scripts/20_score_test_pairs.sh dryrun                       # Part 1 dev + IFT pairs, ~5.5 min
.venv/bin/python scripts/21_evaluate_test.py --mode dryrun  # -> results/test_eval_dryrun/
.venv/bin/python -m pytest tests/test_testeval.py -q        # 29 tests, ~10 s
```

### 8.2 The look (once, after the freeze commit)
```bash
git add docs/TEST_PROTOCOL.md src/testeval scripts tests results/test_eval/fitted_rules.json results/test_eval/score_records/train results/test_eval_dryrun
git commit -m "Freeze test protocol v1.0"                                       # = protocol commit
scripts/20_score_test_pairs.sh test                         # label-free scoring of the test pairs (~3 min); refuses unless frozen
.venv/bin/python scripts/21_evaluate_test.py --mode test    # THE single look (~20-25 min; dry run C: 17 min for 635k pairs)
git add results/test_eval && git commit -m "Test look: results (unchanged)"
```

**Outputs**
- `data_cache/test_eval/scores/test/<method>.tsv` (gitignored);
- `results/test_eval/score_records/test/<method>.run.json`;
- `results/test_eval/looks.jsonl`: a `started` record written **before** any test label is read, and a `done` record;
- `results/test_eval/test_results.json`, which holds:
  - the protocol commit, HEAD and the SHA-256 of this file;
  - the SHA-256 of `fitted_rules.json`;
  - the Part 1 lock state;
  - the SHA-256 of the input files and of every score TSV;
  - the audit counts and all results;
  - runtimes;
- `results/test_eval/figures/reliability_test_itt.png`.

### 8.3 Single look and deviations

**Single look**
- `--mode test` refuses to run if any test look is already logged.
- A re-run is possible only with `--relook-reason "..."`, and only for a crash or I/O failure. It is logged. If a completed (`done`) look exists, its result is marked `confirmatory: false`.
- A changed protocol commit after a logged look is refused.
- Results are committed as produced.

**Fixes after the look**
- A bug found after the look in harness code that does not change any definition here is fixed in a new commit. The look-1 numbers and the corrected numbers are both reported, as a deviation in `results/test_eval/DEVIATIONS.md`.
- A Part 1 bug is reported, not fixed (D0).

**Independent verification.** An agent that did not write the harness re-derives AP, AUROC and ΔAP from the score TSVs and the label files without importing `testeval` (Part 1 practice).

## 9. Decisions for the PI before freezing

1. **CI requirement (§7.1)** in the confirmatory criterion: keep (recommended, given dry-run C power), or revert to Part 1's IFT gate (§7.2 and §7.3 only).
2. **R2 / R4:**
   - R2: self pairs excluded (recommended) or kept at the floor.
   - R4: train-overlap pairs kept in the primary analysis with S2 as sensitivity (recommended), or excluded from the primary.
3. **Primary operating point:** rule (b) F1-max on train (recommended) vs rule (c) top-k.
4. **Fitted rules:** freeze `results/test_eval/fitted_rules.json` (train-only fit, already produced) in the same commit (recommended), so that no fitting happens at look time.
5. **Part 1 audit:** finish `docs/part1_audit.md` (no test access during Part 1) before the look, so that the audit cannot be influenced by test results.
6. **Disclosure:** the harness author printed the first 6 lines of the positive test file (format check, §0). Record this in the freeze note.

## 10. Dry run on development data (`results/test_eval_dryrun/`, no test data)

Details are in `results/test_eval_dryrun/README.md` and `summary.json`.

### A. Reproduction of Part 1's IFT look 1

The IFT pairs were re-derived from the train files and Part 1's hash-checked assignment, giving 847 / 40,786 mapped pairs and 9 / 2,352 unmapped pairs, identical to Part 1. They were scored by the locked predictor and evaluated with Part 1's definitions (HT weights, f = 0.26481).

| method | HT-AP, harness | HT-AP, Part 1 look 1 |
|---|---|---|
| zavg_BEF | 0.345139 | 0.345139 |
| B | 0.317398 | 0.317398 |
| GLS | 0.242318 | 0.242318 |
| Spearman | 0.146273 | 0.146273 |
| Pearson | 0.134064 | 0.134064 |
| PC-removed Pearson (look `2-posthoc`) | 0.266779 | 0.266779 |

**Agreement**
- All 29 fields of Part 1's `metric_suite` agree to ≤ 1.5 × 10⁻⁸, including macro AP, AP without the largest component, the essentiality strata and AP with unmapped pairs at the minimum score (0.332 / 0.233).
- The largest difference is B's AUROC, from predictor-vs-cache float rounding at tied scores.
- The 1,000-resample paired cluster-bootstrap CIs agree to 9 × 10⁻¹⁶, for example zavg_BEF − GLS: +0.103 [+0.035, +0.209].
- The harness's own weight generator gives bit-identical bootstrap CIs.

### B. Full protocol on the IFT (rules fitted on Part 1 dev pairs)

**Design.** Evaluation used HT design weights; the bootstrap units were clusters rebuilt from IFT positives.

**ITT**

| method | AP | AUROC | macro AP |
|---|---|---|---|
| zavg_BEF | 0.332 | 0.754 | 0.355 |
| GLS | 0.233 | 0.714 | 0.359 |
| Spearman | 0.141 | — | — |
| Pearson | 0.129 | — | — |

**ΔAP, cluster bootstrap B = 2,000**

| comparison | ΔAP [95% CI] |
|---|---|
| zavg_BEF vs GLS | +0.099 [+0.034, +0.201] |
| zavg_BEF vs Spearman | +0.191 [+0.094, +0.302] |
| zavg_BEF vs Pearson | +0.203 [+0.098, +0.324] |

**Bootstrap sensitivity.** CI endpoints vary by ≤ 0.011 across seeds 1–4. The protein-bootstrap CI is narrower ([+0.046, +0.162] vs GLS) and the pair-bootstrap CI narrower still ([+0.075, +0.124]), as expected.

**Verdict.** zavg_BEF "clearly beats" all three baselines under §7. B fails against GLS on macro-AP non-inferiority (0.288 < 0.95 × 0.359), as in Part 1.

**Rule (b) on the IFT**

| method | P | R | F1 |
|---|---|---|---|
| zavg_BEF | 0.60 | 0.32 | 0.41 |
| GLS | 0.69 | 0.21 | 0.32 |

**Isotonic calibration (fit on dev)**

| method | Brier skill | calibration slope | calibration-in-the-large |
|---|---|---|---|
| zavg_BEF | 0.22 | 0.97 | 0.87 |
| GLS | 0.16 | — | — |

Platt is clearly worse than isotonic: its single-slope logistic misfits the upper tail.

**Mapped-only.** It equals Part 1's IFT: AP 0.345 / 0.242.

**Definitional difference.** "AP without the largest component" on the IFT is 0.374 here against 0.470 in Part 1. Part 1 removed the *train* giant component; the protocol removes the largest component of the evaluated set. The test uses the protocol definition.

### C. Test-mode code path on all hu.MAP train pairs

**Design.** All 3,314 / 631,465 train pairs, 38,182 of them unmapped, with weights 1. Rules were fitted on the same pairs, so the result is in-sample and only checks mechanics and power.

**Clusters.** 245 clusters, with the 144-protein giant component split into 7 communities. Kish n_eff = 25.8.

**ITT AP**

| method | AP |
|---|---|
| zavg_BEF | 0.389 |
| B | 0.366 |
| GLS | 0.245 |
| Spearman | 0.214 |
| Pearson | 0.195 |
| PC-removed Pearson | 0.295 |

**ΔAP(zavg_BEF − GLS).** +0.144 [+0.078, +0.217].

**Macro AP.** zavg_BEF 0.298 vs GLS 0.291 (non-inferiority passes with little margin).

**Runtime** (8-core laptop)

| step | time |
|---|---|
| scoring 635k pairs, all 6 methods | 183 s |
| evaluation, ITT with all bootstraps | 830 s |
| mapped-only | 144 s |
| total | ~17 min |

### I. Batch invariance

Every IFT pair has the same score (max |Δ| = 0) and the same mapped flag whether it is scored in the IFT file or in the all-train file, for all 6 methods.

### Tests

`tests/test_testeval.py` has 29 tests, all passing. They cover:
- hand-computed AP, AUROC and F1;
- tie blocks;
- equality with sklearn and with Part 1;
- isotonic and Platt fits;
- ECE and Brier;
- equality with Part 1's cluster bootstrap;
- bootstrap CI coverage on a clustered simulation (0.85–1.0);
- a positive control (perfect scorer → AP 1);
- negative controls (random scores and label shuffle → AP ≈ prevalence);
- the pair rules R1–R4;
- the full look on toy files, including refusal of a second look;
- guards that the test files are referenced only by the entry points and that the imported `src` is the locked clone.
