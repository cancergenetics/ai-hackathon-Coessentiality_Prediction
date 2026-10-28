# TEST_PROTOCOL.md: pre-registered evaluation of the locked Part 1 method on the withheld hu.MAP test pairs

**DRAFT v1.2** (2026-10-28; revised after design review rounds 1 and 2, see `docs/reviews/round1_protocol_response.md` and §9).

**Freezing.** At the freeze the PI replaces the line above with `**FROZEN v1.2** (date)`. The pre-freeze status word (the one in the line above) must then occur nowhere in this file. The PI then runs `scripts/22_freeze_protocol.py` and commits this file together with `docs/TEST_PROTOCOL.FROZEN`. Agents never make the freeze commit.

**The frozen marker.** `docs/TEST_PROTOCOL.FROZEN` lists the SHA-256 of:
- this file;
- `results/test_eval/fitted_rules.json`;
- every harness file (`src/testeval/*.py`, `scripts/19_*`, `scripts/20_*`, `scripts/20b_*`, `scripts/21_*`, `scripts/22_*`, `tests/test_testeval.py`);
- the venv package list (`pip freeze`);
- the two test files, which must equal the DATA.md values.

**What every test-touching entry point checks.** `scripts/19_make_pair_files.py --set test`, `scripts/20_score_test_pairs.sh test` and `scripts/21_evaluate_test.py --mode test` call `require_frozen()`. It verifies:
- the marker is committed;
- every listed hash equals the current file;
- the header reads FROZEN and the file contains no pre-freeze status word;
- `src/`, `scripts/`, `tests/`, this file, the marker and the fitted rules are committed and unmodified.

**What may not change.** Nothing in this document changes after the freeze. Any later change is a deviation (§8.3).

---

## 0. Transition from development to test: everything seen or computed on the test files before the freeze

| # | who / when | what was seen or computed | consequence |
|---|---|---|---|
| T0 | PI, 2026-10-28 (M0) | line counts (3,312 / 597,426) and the first 3 lines of each file (format check) | none |
| T1 | Part 1 audit agent, commit 7d1fbdd (`docs/part1_audit.md` §4, `scripts/verify/audit_part1.py`) | **descriptive counts without scores** (see the next table) | these counts were known before this protocol was written (v1.0 in 7a842f2). R4, R5 and S2 were therefore set with them available. The choices follow Part 1 DESIGN §3.10 and are kept (§3.1). S2 is redundant (overlap 0) and is kept only as an automatic check |
| **D1** | Part 1 audit agent, intermediate run (RESEARCH_TREE D1) | counted, per test file, how many test pairs occur in Part 1's label-free genome-wide top-pair table (zavg_BEF ≥ 3) | this is **equivalent to precision and recall of rule (a) z ≥ 3 on the test set**. The counts were deleted from all outputs and the PI never saw them; z ≥ 3 is Part 1's release cutoff and predates the incident. Rule (a) is therefore **non-blind**. It is reported for completeness only, labelled "non-blind (D1)", is secondary, and enters no verdict (§5.4) |
| T2 | harness author (this protocol), 2026-10-28 | an octal dump of the first ~200 bytes (6 lines) of the positive file, to check delimiter and line endings | none |
| T3 | everyone | genome-wide label-free scores for all gene pairs exist (Part 1 `outputs/all_pairs`, `results/predictions/genome_wide_top_pairs.tsv.gz`) and therefore cover the test pairs; no one has joined them with test labels except in D1 | stated as a limit of the single-look guarantee (§8.4) |

**Counts from T1 (no scores).** "Train-positive proteins" are the proteins of the train positive file.

| quantity | value |
|---|---|
| unique test pairs, pos / neg | 3,312 / 597,426 |
| duplicates, self pairs, pairs in both files | 0 / 0 / 0 |
| test pairs also in a train file | 0 |
| test proteins | 1,121 |
| of which also in train | 329 |
| of which in train-positive proteins | 322 |
| test positives in C1 / C2 / C3 (both / one / no protein in train-positive proteins) | 265 / 738 / 2,309 |
| test proteins mapped to DepMap | 1,087 / 1,121 |
| scorable pairs (both proteins mapped, distinct genes), pos | 3,159 / 3,312 (95.4%) |
| scorable pairs, neg | 563,648 / 597,426 (94.3%) |

**Other facts at the freeze**

| item | status |
|---|---|
| Part 1 lock | commit `c7f95ee007576be9311c0a5f08c6ac63ca23a009` (2026-10-28 12:55:43 UTC), cloned read-only to `external/part1_locked/` |
| Test files | SHA-256 in DATA.md (`c4e3e68e…`, `c4a3e546…`); created 13:55 UTC, 59.6 min after the Part 1 lock |
| Part 1 never accessed test labels | **audit complete: PASS** (`docs/part1_audit.md`): timeline, code and document references, content leakage, hyper-parameter selection |
| Development data used by this protocol | hu.MAP **train** pairs (Part 1 folder, read-only); Part 1's dev/IFT split (hash-checked); DepMap 26Q3 CRISPRGeneEffect; HGNC |

Apart from T0–T3 and D1, everything below was fixed from Part 1's documents and code and from development data only.

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
- **Fixed reporting language.** Whatever its test AP, it is reported as "a post hoc PC-removed Pearson reference, not a candidate". If it exceeds zavg_BEF, the text reports that fact with its CI and does not state or imply that a better method exists.

## 3. Pair set

### 3.1 Rules

These rules are fixed a priori and applied to every method identically (`src/testeval/pairs.py`).

| rule | decision | justification |
|---|---|---|
| R1 unordered pair, de-duplication | canonical `(a, b)` with `a < b` on the UniProt accession string (Part 1 `data.load_pairs` rule); a pair listed several times counts once | the task is an unordered pair (Part 1 DESIGN §1). T1 found no duplicates; the rule is kept as a check |
| R2 self pairs (`a == b`) | **excluded**; counted | not a pair of two proteins; Part 1's loader rejects them. T1 found none |
| R3 pair in both the positive and the negative file | **excluded**; counted | label undefined. T1 found none |
| R4 pair also in a train file | **kept** in the primary analysis; flagged | the delivered test files define the population, and all scored methods are label-free. T1 found 0 such pairs |
| R5 unmapped pair (a protein without a unique DepMap gene, both proteins on one gene, or an undefined score) | **kept** in the primary analysis at the locked predictor's floor, min(mapped scores in the scored file) − 1 (Part 1 DESIGN §3.10). This is intention-to-score (ITT). The mapped-only set S1 is Part 1's own evaluation universe | coverage is part of real-world performance. All unmapped pairs form one tied bottom block for every method, so the rule treats all methods identically. Part 1 reported both variants on the IFT (0.332 ITT vs 0.345 mapped). T1 showed 4.6% of positives and 5.7% of negatives are unscorable; the rule was not changed in response |

### 3.2 Scoring input and analysis sets

**Scoring input.** One ID-only file per set: the unique canonical non-self pairs of positives ∪ negatives, sorted, no labels (`scripts/19_make_pair_files.py`). Positives and negatives therefore share one floor.

| set | definition | role |
|---|---|---|
| **ITT** | all canonical pairs remaining after R2 and R3 | **primary metric (§4)**; locked verdict also reported (§7) |
| **S1 mapped-only** | pairs mapped under every method. Clusters are rebuilt from S1's own positives, exactly as Part 1 built them on scorable pairs | **set of the confirmatory locked verdict (§7.1)**, because it is Part 1's evaluation universe |
| S2 no train overlap | drop pairs that occur in a train file | automatic check only: analysed iff any overlap exists (T1: 0, so S2 = ITT and is reported as such) |
| **C1 / C2 / C3** | pairs with both / one / neither protein among train-positive proteins (audit definition); negatives assigned by the same rule | secondary strata (generalisation to new proteins). They replace the v1.0 set "S3, both proteins unseen in any train pair", which C3 approximates (7 proteins differ) |

**Expected composition (T1).** Prevalence 0.551% (train 0.522%). Positives split C1 / C2 / C3 = 265 / 738 / 2,309, so the test is mostly protein-new, more so than Part 1's IFT.

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

All are computed on the ITT set and on S1, for every method. AP and ΔAP (with cluster-bootstrap CIs) are computed for C1/C2/C3.

### 5.1 Ranking metrics
- **AUROC:** weighted, ties count ½; equals sklearn `roc_auc_score`.
- **Precision at recall 0.10 / 0.25 / 0.50:** Part 1 `weighted_precision_at_recall`, using the first tie block that reaches the recall level.
- **AP with only hard negatives:** negatives in the same positive-graph component (Part 1).

### 5.2 Cluster-structured metrics

**Test clusters** are defined exactly as in Part 1 §3.1, applied to the test positives (Part 1 `src.clusters.build_clusters`, imported unchanged):
- universe: all proteins of the analysed set: the ITT set (including unmapped proteins) for ITT, C1–C3 and S2; the scorable proteins of S1 for S1, as Part 1 did;
- graph: that set's positive pairs;
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
| (a) z ≥ 3 | t = 3.0 on the zavg_BEF score | zavg_BEF only | secondary, **non-blind (incident D1, §0)**: the cutoff of Part 1's released genome-wide `top_pairs.parquet`. It is reported only with the label "non-blind (D1)", never as headline P/R/F1, and enters no verdict |

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
- Part 1 pre-specified no label-tuned threshold. Its only fixed cutoff is z ≥ 3, which is non-blind after D1. Part 1 reported precision at recall levels, which are kept in §5.1.

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
- **Statistics bootstrapped (ITT and S1):**
  - per method: AP, AUROC, AP without the largest component, macro per-cluster AP, F1 / precision / recall at rule (b), and isotonic Brier;
  - **macro AP in a replicate:** a cluster drawn m_c times counts m_c times; its per-cluster AP uses the replicate pair weights. At multiplicity 1 this equals Part 1's value, which is tested;
  - paired Δ and ratios:
    - zavg_BEF − each baseline: AP and AUROC;
    - zavg_BEF and B − each baseline: AP, macro AP, AP without the largest component. The ratios of the last two are the descriptive non-degradation CIs (§7.2);
    - B − zavg_BEF (AP);
    - PC-removed Pearson − zavg_BEF (AP).
- **Strata C1–C3:** AP and ΔAP only.
- **Fixed during resampling:** thresholds and calibrators are fitted on train and not refitted, so their sampling error is not propagated (stated as a limitation).

### 6.2 Secondary bootstraps (AP and ΔAP only)
- **Protein (node) bootstrap:** proteins are resampled; pair weight m_a · m_b; B = 2,000, seed 0.
- **Stratified pair bootstrap:** Part 1's `_boot_weights`; positives and negatives resampled separately; B = 2,000, seed 0. It ignores dependence between pairs that share proteins, so it is reported as anti-conservative.

### 6.3 Monte-Carlo stability
The cluster bootstrap of ΔAP (zavg_BEF − each baseline) is repeated with seeds 1–4 (B = 2,000 each). The verdict uses seed 0 only. The other seeds show how much the CI endpoints depend on the seed.

## 7. Decision rule, additional statistic and multiplicity

### 7.1 Locked decision rule (Part 1, unchanged)

The decision is Part 1's locked IFT "clearly beats" rule, as Part 1 defined and coded it (DESIGN v3.1 §3.7, IFT parts). It is computed by Part 1's own function `external/part1_locked/scripts/final_eval.py::decide`, called unchanged.

For candidate c and baseline b:

| criterion | Part 1 name | condition |
|---|---|---|
| (ii) effect size | `gate_ii_ift` | relative AP gain AP_c / AP_b − 1 ≥ 0.10 (point estimate) |
| (iii) non-inferiority | — | macro per-cluster AP_c ≥ 0.95 × that of b, **and** AP-without-largest-component_c ≥ 0.95 × that of b |

**Combination rules**
- c clearly beats the baselines iff both criteria hold against **each** of GLS, Spearman and Pearson (intersection–union).
- **Fixed sequence:** zavg_BEF first, which is confirmatory. B is confirmatory only if zavg_BEF passes, otherwise exploratory.
- Bootstrap CIs are information only, as in Part 1.
- Part 1's dev parts (i) and (ii-dev) were met before the lock and are cited by reference (Part 1 `results/phaseB_v2/`).

**Pair set**
- **Confirmatory verdict:** computed on **S1 (mapped-only)**. That is Part 1's evaluation universe (scorable pairs; Part 1 DESIGN §3.10), with clusters and the largest component rebuilt on S1's positives.
- **Also reported:** the same rule on ITT.
- **Plain AP:** `ap_ht` in Part 1's suite equals plain AP, because all weights are 1.

**One unavoidable adaptation.** Part 1's "largest component" was the train giant component; on the test it is the largest component of the test positives (§5.2).

**zavg_BEF stays the final method whatever the outcome.**

### 7.2 Additional pre-specified statistic (reported regardless; not part of the decision)

For each candidate vs baseline, on ITT and on S1:
- **CI-excludes-0:** whether the 95% paired cluster-bootstrap CI of ΔAP (seed 0, B = 2,000) has lower bound > 0.
- **Joint statement:** "zavg_BEF CI-excludes-0 against all three baselines" (intersection–union, one-sided level 2.5%).
- **Borderline flag:** set if |CI lower bound| ≤ 0.011, the largest seed-to-seed change of the CI endpoints in dry run B. Seed 0 stays decisive.
- **Non-degradation ratios:** bootstrap 95% CIs of the ratios macro-AP_c / macro-AP_b and AP-without-largest_c / AP-without-largest_b, descriptive. The point-estimate gate of §7.1 can be decided by a small margin; it passed by about 3% (macro) in dry run C.
- **Holm-adjusted p-values** of the three zavg_BEF ΔAP comparisons.

**Why this statistic exists.**
- Part 1 dropped the CI requirement on the IFT only because of power (n_eff ≈ 7.7 clusters).
- The test is about the size of the train set. In dry run C (all train pairs; 245 clusters; Kish n_eff 25.8), ΔAP(zavg_BEF − GLS) was +0.144 with 95% CI [+0.078, +0.217].
- The CI statistic is therefore informative here, but by PI decision it does not change the locked rule.

**Multiplicity.** The decision is one intersection–union family plus a fixed sequence, so no adjustment is needed. Everything else (§5, S2, C1–C3, the reference, the secondary bootstraps, the ITT verdict) is descriptive, with unadjusted 95% CIs.

### 7.3 Pre-written readings of the possible outcomes

| outcome | reading |
|---|---|
| Locked rule passes **and** the CI excludes 0 against all baselines | "The final method, locked before the test files existed, outperforms Pearson, Spearman and GLS on held-out hu.MAP complexes (AP …, ΔAP vs GLS … [CI])." |
| Locked rule passes, CI includes 0 against at least one baseline | "The pre-registered decision rule was met; the improvement over <baseline> is not statistically resolved at 95% (ΔAP … [CI])." No claim of significance against that baseline |
| Locked rule fails | "On held-out complexes the final method did not meet the pre-registered criterion against <baseline> (failed: <criterion>)." The test numbers are reported in full. The method is not changed. Part 2 still uses zavg_BEF (locked), and its value is then argued only from Part 2's own label-free comparisons |

**Which set governs the reading.** The reading is selected by the **S1** locked rule (§7.1) together with the **S1** CI statistic (§7.2).
- The quoted AP and ΔAP are given for ITT first, then for S1.
- If the ITT and S1 verdicts (or CI statistics) disagree, S1 governs, and the disagreement is reported explicitly in the text, with both results.

In every case:
- the C1–C3 strata are reported, showing whether gains hold on protein-new pairs;
- rule (a) is reported as non-blind.

**Firewall.** No Part 2 parameter or claim may cite or be tuned on test-look results, except the final reporting of this evaluation.

## 8. Execution

### 8.1 Before the freeze (development data only; already run by agents)
```bash
scripts/20_score_test_pairs.sh train                        # locked predictor on all train pairs, ~3 min
.venv/bin/python scripts/21_evaluate_test.py --mode fit     # -> results/test_eval/fitted_rules.json, ~6 s
scripts/20_score_test_pairs.sh dryrun                       # Part 1 dev + IFT pairs, ~5.5 min
.venv/bin/python scripts/21_evaluate_test.py --mode dryrun  # -> results/test_eval_dryrun/
.venv/bin/python -m pytest tests/test_testeval.py -q
```

### 8.2 Freeze and the look (PI only, in this order)
```bash
# (1) freeze: edit the header line to "**FROZEN v1.2** (date)", then
.venv/bin/python scripts/22_freeze_protocol.py              # writes docs/TEST_PROTOCOL.FROZEN (refuses on the pre-freeze status word or a test-file hash mismatch)
git add docs/TEST_PROTOCOL.md docs/TEST_PROTOCOL.FROZEN && git commit -m "Freeze test protocol v1.2"
#     -> report the commit hash and the marker SHA-256 to the user in chat (notarisation)
# (2) label-free scoring of the test pairs (refuses unless frozen), ~3 min
scripts/20_score_test_pairs.sh test
# (3) look intent: no label is read; writes an 'intent' record to the git-tracked results/test_eval/looks.jsonl
.venv/bin/python scripts/21_evaluate_test.py --mode test --intent
git add results/test_eval/looks.jsonl results/test_eval/score_records/test && git commit -m "Test look intent"
# (4) THE look (refuses without a committed intent, or if any earlier look exists on disk, in any ref or in the reflog);
#     it commits the 'started' record itself before reading labels and prints that hash -> report it to the user. ~30-35 min
.venv/bin/python scripts/21_evaluate_test.py --mode test --run
git add results/test_eval && git commit -m "Test look: results (unchanged)"
#     -> report the results commit hash and the SHA-256 of test_results.json to the user in chat
```

**Timing and environment.** Steps (2)–(4) run **immediately after** the freeze commit, in the same session.
- No package may be installed, removed or upgraded in `.venv` between the freeze and the end of the look; the marker's `pip freeze` hash enforces this.
- Part 2 work that needs new packages waits until the look is committed, or uses a separate environment.

**Outputs**
- `results/test_eval/looks.jsonl` (git-tracked): `intent` (committed by the PI), then `started` (committed by `--run` **before** any label is read), then `done`.
- `results/test_eval/test_results.json`: protocol and marker commits, marker SHA-256, the SHA-256 of `fitted_rules.json`, the Part 1 lock state, the input and score-TSV SHA-256s, audit counts, all results, the verdict summary and runtimes.
- `results/test_eval/figures/reliability_test_itt.png`.
- Score TSVs: `data_cache/test_eval/scores/test/` (gitignored); their SHA-256 are recorded.

### 8.3 Single look, re-runs and deviations

**Refusals**
- `--run` refuses if a `started` or `done` record exists in the ledger on disk, **in any commit reachable from any ref, or in any commit recorded in the reflog** (`git log --all` and `git log --reflog`). It also refuses if `test_results.json` exists on disk or in that history.
- `--run` **commits the `started` record itself** (git add + commit of `looks.jsonl` only) **before any label is read**. It prints the commit hash, and the PI reports it to the user. If the commit fails, the look stops before reading labels.
- This closes the partial-look gap (a crash after results are written but before `done`) and the erase-and-re-run gap (reset, deletion, a deleted branch). The latter is tested in a real repository.

**Re-runs**
- A re-run needs `--relook-reason`. It is logged and **automatically `confirmatory: false`**. The first look's results are reported as the confirmatory ones, if they exist.
- If the first run crashed before writing any results, the re-run is still marked non-confirmatory, and the crash is reported.

**Fixes after the look.** A bug found after the look in harness code that does not change any definition here is fixed in a new commit. The original and corrected numbers are both reported, in `results/test_eval/DEVIATIONS.md`. A Part 1 bug is reported, not fixed (D0).

**Independent verification (pre-registered)**
- **Verifier A** is an agent that did not write the harness. It receives only this frozen protocol, the label files and the score TSVs, and implements independently, without importing `testeval`:
  - R1–R5 and the floor rule;
  - the test clusters (Part 1 `build_clusters` may be imported);
  - AP, AUROC, macro AP and AP without the largest component;
  - the seed-0 cluster bootstrap;
  - the §7.1 verdict and the §7.2 statistic.
- Verifier A reports **before** opening `test_results.json`.
- **Tolerances:** point metrics ≤ 1e-9; bootstrap CI endpoints ≤ 1e-9 (same RNG stream); the verdict identical.
- **Verifier B** re-scores a 1% random sample of the test pairs (seed 0) from the raw CSV through the locked predictor, in a separate file, and checks the scores are identical (batch invariance, §10.I).

### 8.4 What the single-look guarantee is (and is not)

**The guarantee is procedural, not technical.**
- The test labels are plain files that any process with file access can read.
- Genome-wide label-free scores covering all test pairs already exist (T3), so test AP could be computed without this harness.
- Incident D1 shows this risk is real.

**Mitigations.** There is no escrow of the test files (PI decision); they stay in place. The mitigations are:
1. Every agent brief forbids reading test files and computing test metrics.
2. The freeze marker hashes the harness, the protocol, the fitted rules, the venv package list and the test files (checked against DATA.md).
3. The look ledger is git-tracked. The intent is committed before the look, and the `started` record is committed by `--run` before any label is read.
4. Refusal on any earlier look found on disk, in any ref or in the reflog. Any re-run is automatically non-confirmatory.
5. The PI reports the freeze commit and marker hash, the intent commit, the `started` commit and the results commit with the SHA-256 of `test_results.json` to the user in chat, outside the agents' reach.
6. A recursive grep test limits references to the test files to the entry points (plus the completed audit script and grep-guard pattern lists).
7. The train files are pinned by SHA-256 (Part 1 DATA.md), so `PART1_DIR` can relocate them but never substitute them. `--mode dryrun` refuses inputs that are not pinned-train pairs with the same label.
8. All prior contact is disclosed in §0.

**Residual limits.** A process that reads the label files directly, or that expires the reflog, is not technically prevented. The guarantee rests on items 1, 5 and 8.

## 9. Choices fixed in v1.1 (no open options) and changelog

**Fixed choices**
- **Pair rules:** R2 excludes self pairs; R4 keeps train-overlap pairs (none exist).
- **Primary operating point:** rule (b), F1-max on train.
- **Rule (a) z ≥ 3:** non-blind, secondary.
- **Fitted rules:** `fitted_rules.json` is frozen in the marker; no fitting happens at look time.
- **Decision:** the Part 1 locked rule on S1, plus the CI statistic, both reported regardless.

**v1.0 → v1.1 (design review round 1; the item-by-item table is in `docs/reviews/round1_protocol_response.md`)**

| area | change |
|---|---|
| freeze | marker file with hashes and a header check; `require_frozen` re-verifies every hash and a clean `tests/` |
| look | two steps, intent then run; git-tracked ledger; refusal on prior looks found in git history; re-runs automatically non-confirmatory |
| §0 | lists T1 (audit counts) and D1; audit marked PASS; rule (a) non-blind |
| §7 | the decision is Part 1's locked rule (Part 1's own `decide`) on S1; the CI-excludes-0 criterion becomes an additional statistic, with a borderline flag and non-degradation CIs; outcome readings written in advance |
| strata | S3 replaced by the audit strata C1–C3; S2 reduced to an automatic check |
| bootstrap | macro AP added |
| verification | verifier specification and tolerances |
| §8.4 | honest statement of the procedural guarantee |
| post hoc reference | reporting language fixed |

**v1.1 → v1.2 (design review round 2)**

| item | change |
|---|---|
| M-B1 / P-M1 | §10 wording fixed; a test applies the documented header edit and passes the text check |
| M-B3 | `--run` commits the `started` record before labels are read; the prior-look scan includes all refs and the reflog; real-git test |
| M-B2 | escrow dropped (PI decision); guarantee and mitigations restated (§8.4) |
| P-M2 | S1 governs the §7.3 reading; disagreements with ITT are reported |
| P-S1 | verifier CI tolerance ≤ 1e-9 |
| S-B4 | train files pinned by SHA-256 |
| S-B5 | no `.venv` changes between freeze and look; the look runs immediately after the freeze |

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

**Verdict**
- **Locked rule (§7.1), S1 and ITT:** zavg_BEF passes against all three baselines. B fails against GLS on macro-AP non-inferiority (S1 0.297 < 0.95 × 0.369), as in Part 1 look 1.
- **CI statistic (§7.2):** the CI excludes 0 against all baselines for both candidates; no comparison is borderline.
- **Non-degradation ratio, zavg_BEF / GLS macro AP:** 0.99 [0.84, 1.31], wide on the IFT.
- **Strata:** all IFT pairs are C3 relative to the dev positives (protein-disjoint split), so the C3 AP equals the ITT AP.

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

**Macro AP.** zavg_BEF 0.298 vs GLS 0.291. Ratio 1.03 [0.96, 1.15]. The locked rule passes for zavg_BEF; B fails against GLS on macro AP.

**Runtime** (8-core laptop)

| step | time |
|---|---|
| scoring 635k pairs, all 6 methods | 183 s |
| evaluation, ITT with all bootstraps, including macro AP | 1,017 s |
| S1 (rebuilt clusters, full cluster bootstrap) | 492 s |
| stratum C1 | 191 s |
| total, about 28 min of compute | ≈ 1,700 s |

The look on ~600k test pairs is expected to take about 30–35 min, including C1–C3.

### I. Batch invariance

Every IFT pair has the same score (max |Δ| = 0) and the same mapped flag whether it is scored in the IFT file or in the all-train file, for all 6 methods.

### Tests

`tests/test_testeval.py` has 34 tests, all passing. They cover:
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
- the full two-step look on toy files: refusal without an intent, refusal without a committed intent, refusal of a second look even after the ledger is wiped, and an automatically non-confirmatory re-run;
- the freeze gate: a pre-freeze header and changed hashes are refused; the documented header edit to FROZEN passes the text check; the current repository is refused;
- pinned train files: substitutes are refused;
- in a real git repository, a committed `started` record is still found after `reset --hard`, file deletion or deletion of its branch (reflog scan);
- the locked verdict through Part 1's `decide`;
- the bootstrap macro AP equals Part 1's value;
- guards that the test files are referenced only by the entry points and that the imported `src` is the locked clone.
