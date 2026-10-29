# REPORT: co-essential protein complexes, held-out test and cancer-context dependencies

<!-- Source convention (as in manuscript/main.md): every number is followed by an HTML comment naming the committed
file it comes from, with JSON key paths, TSV rows/columns, or "calc:" for arithmetic on cited numbers. Where a number
is quoted from the manuscript, the comment names the manuscript's own source (main text checked 421/421 by
scripts/check_manuscript_numbers.py, RESEARCH_TREE.md D21). Aliases:
TR   = results/test_eval/test_results.json
EXP  = results/test_eval/export/ap_by_size.tsv
SR   = results/test_eval/score_records/test/<method>.run.json
TV   = results/test_eval/verify/README.md
RT   = RESEARCH_TREE.md
MS   = manuscript/main.md
SUP  = manuscript/supplementary.md
PR1  = manuscript/review/peer_review_round1.md
RR1  = manuscript/review/response_round1.md
PR2  = manuscript/review/peer_review_round2.md
RR2  = manuscript/review/response_round2.md
MAN  = results/part2/lockbox/manifest.json
CON  = results/part2/construction/README.md
SUMC = results/part2/construction/summary/summary.json
TIERS= results/part2/construction/tiers/summary.json
SEL  = results/part2/construction/selection/selection.json
VAR  = results/part2/construction/variability/variability.json
BR   = results/part2/construction/bridge/bridge.json
CTRL = results/part2/controls
DR   = results/part2/directions
LB   = results/part2/lockbox_confirm
RATES= LB/rates_primary_R7.tsv
CMP  = LB/comparison_primary_vs_corrected_R7.tsv
BLR7 = LB/bridge_primary_vs_corrected_R7.json
BENCH= LB/BENCH_complex_vs_gene/bench.json
LBR  = LB/REPORT.json
LD   = manuscript/tables/lockbox_derived.json
PD   = manuscript/tables/part2_derived.json
RV   = results/part2/verify/lockbox_results_review.md
ROB  = results/part2/robustness/summary.json
ANN  = results/part2/robustness/annotations_summary.json
ROBR = results/part2/robustness/README.md
SF   = results/part2/robustness/strongest_findings.tsv
SD   = manuscript/tables/strongest_derived.json
REV  = results/part2/revision
REVR = results/part2/revision/README.md
R1   = REV/r1_summary.json
R1A  = REV/r1a_rates.tsv
R2   = REV/r2_summary.json
R3   = REV/r3_summary.json
R4   = REV/r4_summary.json
R5   = REV/r5_summary.json
RVW  = results/part2/verify/revision_analyses_review.md
AMD  = docs/registry/AMENDMENTS.md
RETRO= docs/retros/retro_M6_final.md
-->

Status: final report after two rounds of peer review. Round 1 asked for major revision and round 2 for minor revision; every round-2 item has been addressed (Section 7). Submission now waits only on placeholders that the human authors must fill in (Section 10). The manuscript is `manuscript/main.md` (with `manuscript/supplementary.md`; clean copy in `manuscript/submission/`). This report summarises the project for the record and adds no analysis of its own. <!-- RT: D18, D21; RR2 -->

---

## 1. Executive summary

The locked Part 1 co-complex score, zavg_BEF, was evaluated once on the withheld hu.MAP 3.0 test pairs. It reached average precision (AP) 0.333 [95% CI 0.256, 0.410], against 0.212 for GLS, 0.204 for Spearman and 0.197 for Pearson (ΔAP vs GLS +0.121 [+0.079, +0.159]). <!-- TR: results.itt.bootstrap.cluster.stat.*|ap; results.itt.bootstrap.cluster.diff."zavg_BEF|ap - gls_signed|ap" --> The pre-registered verdict passed. The gain comes from pair-rich test clusters: macro per-cluster AP is no better than GLS's (0.297 vs 0.292). <!-- TR: results.itt.metrics.*.macro_cluster_ap -->

In Part 2 the score built a map of 367 complexes without using any curated complex. <!-- SUMC: zavg.n_level1 --> Neither bridge from pair accuracy to complex-level value held. The primary (BioPlex) bridge is null (−0.025 [−0.075, +0.031]). The secondary bridge (BRIDGE-L) is negative: zavg-map context discoveries confirmed in held-out lines less often than GLS-map ones (−0.065 [−0.107, −0.025]). <!-- BR: primary.diff, primary.ci95; BLR7: primary_evaluable.primary.diff, primary_evaluable.primary.ci95 --> Post hoc, power explains only 26–41% of that gap. <!-- REVR: "R-4", "26–41%" -->

Gene sets from other co-dependency modules, matched for size, essentiality and coherence, associated with lineages as often as complexes (28.1% vs 28.1%). Sets matched for essentiality alone did so less often (16.6%, post hoc). <!-- DR/H3a_lineage_atlas/calibration_D8.json: matched_random_sets.pooled.frac_p_lt_005, matched_random_sets.real_pooled.frac_p_lt_005; R1A row H3a_lineage_atlas all_units nonmember_ess pooled: frac_p_lt_005 --> Lockbox confirmation (H3a 52.5%) matched what power predicts. <!-- RATES row H3a_lineage_atlas overall: primary_rate_evaluable; LD: power_band --> The complex-versus-best-gene benchmark (BENCH) is inconclusive.

A registered rule leaves 254 strongest findings. <!-- ROB: per_direction.*.strongest --> Of the 37 post hoc lead examples, 26 lie beyond the lineage-essentiality trend in its own direction, so the uniform scale rule cannot exclude an artefact for them. DEMETER2 RNAi agrees in sign for 29 of 35 evaluable leads. <!-- R2: lead_corrected, lead_corrected_concordant_nl; R3: leads.n_evaluable, leads.same_sign --> Within complexes, 11 of 24 genetic losses leave the complex-wide dependency within ±0.1 Chronos units, and none within ±0.05. <!-- PD: H1.n_neutral, H1.n_events; R5: A1_h1_margin_005.n_inside_005 --> A driver defect, found after the lockbox read, strips claims from 53 calls. <!-- LB/README.md "PI ruling R7"; LD: defect_sensitive_total -->

The revision for peer review narrowed claims and added post hoc analyses (R-1 to R-5) on construction set A or committed tables. No registered result, call or metric changed, and no lockbox value or test label was opened. <!-- RR1 intro; RR2 "What was new in this round" -->

---

## 2. Approach, Part 1 → test transition

**What was locked.** Part 1 was locked at commit c7f95ee (28 Oct 2026, 12:55:43 UTC) and used from a read-only clone (`external/part1_locked`). <!-- DATA.md "Part 1 (locked method)"; TR: part1.commit, part1.matches --> zavg_BEF is the equal-weight mean of three component z-scores (correlation F, neighbourhood cosine B, Laplacian-eigenmap cosine E) on gene-standardised CRISPRGeneEffect profiles with the top 20 principal components removed and ridge whitening. It is unsupervised (no parameters fitted to labels), but its configuration was selected on hu.MAP training complexes. <!-- MS: Methods "The locked method"; RR1 M5.2 --> Baselines: Part 1's signed Pearson, Spearman and GLS; B alone was a pre-declared secondary candidate; a PC-removed Pearson reference was added post hoc. Nothing learned here changed the method.

**Transition and audit.** The hu.MAP test files reached the project at 13:55 UTC, 59.6 min after the final Part 1 commit; the split itself had been published earlier. <!-- DATA.md "Notes"; SUP Note 1 (docs/part1_audit.md "Transition timeline"); RR1 M5.1 --> An independent audit (`docs/part1_audit.md`) found no evidence that test pairs, labels or derived information entered Part 1 (PASS on timeline, references, content and hyper-parameters). The test set is pair-disjoint from training (0 overlaps in all four positive/negative combinations). It has 1,121 proteins, 329 of them also in training. Of the 3,312 positives, 2,309 have neither protein among the training-positive proteins (stratum C3). <!-- RT: Tree "Part 1 audit"; TR: audit.train_overlap, audit.proteins_analysed, audit.proteins_also_in_train, audit.c_strata.C3.pos, audit.analysed_pos -->

**Disclosed contacts before the freeze.** <!-- SUP Note 1 table (docs/TEST_PROTOCOL.md §0); RT: D0, D1 -->
- **D0 / T0.** The PI read the line counts (3,312 positive and 597,426 negative lines) and the first 3 lines of each file as a format check. <!-- RT: D0 -->
- **D1.** In an intermediate run, the audit agent counted test pairs in Part 1's unsupervised top-pair table (zavg_BEF ≥ 3), equivalent to precision and recall at the z ≥ 3 rule. The counts were deleted unseen; rule (a) is reported as non-blind and enters no verdict. <!-- RT: D1 -->
- **T1–T3.** Descriptive counts that used no labels; an octal dump of the first ~200 bytes of the positive file; genome-wide scores already covered all test pairs, so the single-look guarantee is procedural, not cryptographic. <!-- SUP Note 1 table -->

**Protocol** (`docs/TEST_PROTOCOL.md` v1.2, three design-review rounds). Primary metric: AP at natural prevalence on ITT (all 600,738 canonical pairs; 33,931 unmapped pairs at a tied floor score). The locked "clearly beats" verdict is computed on S1, the pairs mapped under every method, with ITT alongside (D3). Uncertainty: paired cluster bootstrap over 238 test clusters (B = 2,000, seed 0). Operating point rule (b) and calibrators were fitted on the 634,779 training pairs before the freeze. A dry run reproduced Part 1's internal final test to ≤ 1.5 × 10⁻⁸. <!-- RT: "Previous status (M2)", D3; TR: audit.unique_canonical_pairs, results.itt.n_unmapped_any, results.itt.cluster_stats.n_clusters, results.itt.bootstrap.cluster.n_boot, .seed; MS: Methods (results/test_eval/fitted_rules.json; docs/TEST_PROTOCOL.md §10) -->

**The single confirmatory look.** Freeze commit 039a843 (marker SHA-256 524c4d68… over protocol, rules, code, tests, packages and test files); "started" record 7afcf06 before labels were read; results 0861e59; `prior_look_evidence.any` false. <!-- TR: protocol.commit, protocol.marker_sha256, started_commit, prior_look_evidence.any --> A verifier that had not written the harness re-implemented the protocol and committed its numbers (c69d828) first: every verdict-relevant number matched to ≤ 2.6 × 10⁻¹⁵. <!-- TV: "Verdict", "Order of work" --> A later descriptive export is logged as non-confirmatory. <!-- results/test_eval/looks.jsonl (descriptive_export) -->

---

## 3. Test-set results

### 3.1 Primary and secondary metrics

Point estimates with 95% paired cluster-bootstrap CIs (B = 2,000, seed 0). P/R/F1 are at the primary operating point, rule (b). Calibration is for the isotonic calibrator fitted on training pairs. The protocol pre-specified a bootstrap CI for Brier only; slope and calibration-in-the-large (CITL = Σp/Σy) have none. <!-- TR: results.{itt,mapped_only}.bootstrap.cluster.stat; results.{itt,mapped_only}.metrics.*.calibration.isotonic; SUP Note 10 "Calibration" -->

**ITT (primary metric set; 3,312 positives, 597,426 negatives, prevalence 0.551%)** <!-- TR: results.itt.n_pos, results.itt.n_neg, audit.analysed_prevalence -->

| Method | AP | AUROC | Precision | Recall | F1 | Brier (×10⁻³) | Cal. slope | CITL | Brier skill |
|---|---|---|---|---|---|---|---|---|---|
| **zavg_BEF** | **0.333** [0.256, 0.410] | 0.789 [0.735, 0.829] | 0.553 [0.455, 0.645] | 0.328 [0.259, 0.403] | 0.411 [0.346, 0.472] | 4.27 [3.26, 5.35] | 1.007 | 0.893 | 0.222 |
| GLS | 0.212 [0.156, 0.273] | 0.730 [0.685, 0.777] | 0.567 [0.475, 0.645] | 0.205 [0.155, 0.266] | 0.301 [0.239, 0.365] | 4.72 [3.55, 5.94] | 0.996 | 0.915 | 0.139 |
| Spearman | 0.204 [0.106, 0.352] | 0.756 [0.698, 0.819] | 0.306 [0.214, 0.440] | 0.272 [0.136, 0.448] | 0.288 [0.179, 0.408] | 5.00 [3.89, 6.27] | 0.841 | 1.195 | 0.089 |
| Pearson | 0.197 [0.103, 0.338] | 0.749 [0.692, 0.813] | 0.273 [0.196, 0.380] | 0.282 [0.148, 0.453] | 0.277 [0.178, 0.390] | 4.95 [3.88, 6.18] | 0.870 | 1.166 | 0.097 |

<!-- TR: results.itt.bootstrap.cluster.stat."<m>|ap", "<m>|auroc", "<m>|precision_rule_b", "<m>|recall_rule_b", "<m>|f1_rule_b", "<m>|brier_isotonic"; results.itt.metrics.<m>.calibration.isotonic.calibration_slope, .calibration_in_the_large, .bss_vs_observed_prevalence; m ∈ {zavg_BEF, gls_signed, spearman_signed, pearson_signed} -->

**S1** (locked-verdict set; 3,159 positives, 563,648 negatives): AP 0.349 [0.274, 0.435] for zavg_BEF against 0.222 (GLS), 0.214 (Spearman) and 0.206 (Pearson); AUROC 0.815 against 0.750 / 0.779 / 0.771. Full S1 metrics: `manuscript/tables/table1_test_metrics.md`. <!-- TR: results.mapped_only.n_pos, .n_neg, results.mapped_only.bootstrap.cluster.stat.*|ap, *|auroc -->

Reading:
- zavg_BEF has the highest AP, AUROC, recall and F1; at rule (b) its precision matches GLS's (0.553 vs 0.567) with higher recall. <!-- TR: results.itt.bootstrap.cluster.stat.*_rule_b -->
- Its isotonic calibrator kept its slope (1.007) but under-predicted total positive mass (CITL 0.893; no recalibration). <!-- TR: results.itt.metrics.*.calibration.isotonic -->
- Secondary candidate B: AP 0.310 [0.234, 0.392]. It failed the locked rule against GLS (macro-AP ratio 0.875 < 0.95 on S1) and was below zavg_BEF (ΔAP −0.022 [−0.032, −0.012]). <!-- TR: results.itt.bootstrap.cluster.stat."B|ap", results.itt.bootstrap.cluster.diff."B|ap - zavg_BEF|ap", summary.locked_verdict_S1.B; TV: "Locked verdict" -->
- The post hoc PC-removed Pearson reference reached AP 0.239, below zavg_BEF (−0.094 [−0.126, −0.061]). It closes part of the gap to Pearson, so preprocessing probably contributes to the gain; a preprocessing-matched GLS could not be tested without a second look. <!-- TR: results.itt.metrics.spec_k20_tinf_std.ap, results.itt.bootstrap.cluster.diff."spec_k20_tinf_std|ap - zavg_BEF|ap"; MS: Discussion -->

### 3.2 Paired differences (zavg_BEF − baseline)

| Contrast | ITT ΔAP [95% CI] | Relative gain | S1 ΔAP [95% CI] | ITT ΔAUROC [95% CI] |
|---|---|---|---|---|
| vs GLS | **+0.121** [+0.079, +0.159] | +57% | +0.127 [+0.089, +0.167] | +0.059 [+0.030, +0.090] |
| vs Spearman | +0.129 [+0.023, +0.216] | +63% | +0.135 [+0.024, +0.229] | +0.033 [−0.003, +0.067] |
| vs Pearson | +0.136 [+0.039, +0.220] | +69% | +0.143 [+0.038, +0.233] | +0.040 [+0.002, +0.076] |

<!-- TR: results.itt.bootstrap.cluster.diff."zavg_BEF|ap - *|ap" (diff, ci_low, ci_high, rel_gain); results.mapped_only.bootstrap.cluster.diff."zavg_BEF|ap - *|ap"; results.itt.bootstrap.cluster.diff."zavg_BEF|auroc - *|auroc" -->

- **Locked verdict:** "clearly beats" holds against all three baselines on S1 and on ITT. <!-- TR: summary.locked_verdict_S1.zavg_BEF, summary.locked_verdict_ITT.zavg_BEF -->
- Holm-adjusted p: < 0.003 (GLS; no replicate ≤ 0, corrected in round 1 from a one-sided "< 0.0005"), 0.022 (Spearman), 0.018 (Pearson). <!-- SUP Note 10 "Inference" (TR: holm_adjusted_p_primary); RR1 minor 17 -->
- The ΔAUROC CI against Spearman includes 0. CIs are wide because the 238 clusters have a Kish effective number of 21.3. <!-- TR: results.itt.bootstrap.cluster.diff."zavg_BEF|auroc - spearman_signed|auroc"; results.itt.cluster_stats.kish_neff_clusters -->
- On protein-new pairs (C3), AP was 0.420 against 0.272 for GLS (+0.148 [+0.093, +0.199]). <!-- TR: results.stratum_C3.bootstrap.cluster.stat, results.stratum_C3.bootstrap.cluster.diff."zavg_BEF|ap - gls_signed|ap" -->

### 3.3 Macro-AP and AP by cluster size

Macro per-cluster AP weights each test cluster equally. It does not separate zavg_BEF from GLS: 0.297 vs 0.292 (+0.005 [−0.014, +0.047]). <!-- TR: results.itt.metrics.*.macro_cluster_ap; results.itt.bootstrap.cluster.diff."zavg_BEF|macro_cluster_ap - gls_signed|macro_cluster_ap" --> Without the 10 clusters with most positives, AP was 0.345 vs 0.264 (not bootstrapped). <!-- TR: results.itt.metrics.*.ap_no_top10_clusters -->

Descriptive, post-look breakdown by positive pairs per cluster (ITT). The bins were fixed before the export was run. The CIs resample each bin's clusters and are **not** the protocol bootstrap. <!-- results/test_eval/export/README.md "ap_by_size.tsv" -->

| Positive pairs per cluster | Clusters | Macro AP zavg_BEF | Macro AP GLS | Δ vs GLS [95% CI] |
|---|---|---|---|---|
| 1 | 91 | 0.256 | 0.289 | −0.033 [−0.090, +0.023] |
| 2–5 | 67 | 0.259 | 0.283 | −0.024 [−0.061, +0.014] |
| 6–20 | 43 | 0.380 | 0.339 | +0.041 [+0.008, +0.080] |
| 21–100 | 25 | 0.388 | 0.261 | +0.126 [+0.075, +0.182] |
| 101+ | 6 | 0.348 | 0.206 | +0.142 [+0.078, +0.198] |

<!-- EXP rows universe=itt, bin_by=n_pos_pairs, method ∈ {zavg_BEF, gls_signed}: n_clusters, macro_ap, macro_ap_minus_gls, macro_ap_minus_gls_ci_low, macro_ap_minus_gls_ci_high -->

The gain grows with the number of positive pairs per cluster. For clusters with five or fewer positive pairs the locked score is not better than GLS.

### 3.4 Runtime (laptop)

| Step | Time | Peak RSS |
|---|---|---|
| Scoring all 600,738 test pairs, zavg_BEF (F 4.6 s, B 16.9 s, E 86.0 s) | 112.6 s | 4.0 GB |
| Scoring, GLS / Spearman / Pearson | 17.4 / 7.7 / 7.5 s | 2.3–2.5 GB |
| Scoring, B (secondary) / PC-removed Pearson (post hoc) | 24.7 / 9.8 s | 3.1 GB / – |
| Scoring, all six methods | 179.7 s | – |
| Evaluation: load 8.1 s; ITT metrics + all bootstraps 1,344.1 s; S1 483.8 s; strata C1–C3 161.8 s | 1,998 s (33.3 min) | – |
| of which ITT cluster bootstrap / protein / pair / MC seeds 1–4 | 758.4 / 129.7 / 127.7 / 310.1 s | – |
| Descriptive export | 17.1 s | – |

<!-- SR (zavg_BEF, gls_signed, spearman_signed, pearson_signed, B, spec_k20_tinf_std): runtime.total_s, runtime.scoring.components_s, runtime.peak_rss_gb; RT: Tree "TEST_PROTOCOL v1.2 frozen" (179.7 s, corrected from "3.3 min"; calc: 112.6+17.4+7.7+7.5+24.7+9.8 = 179.7); TR: runtime_s (calc: 8.1+0.6+1344.1+483.8+6.4+65.7+89.7 = 1998.4; 6.4+65.7+89.7 = 161.8), results.itt.bootstrap.runtime_s; results/test_eval/export/README.md -->

Accuracy for cost: zavg_BEF costs about 6.5× GLS's scoring time (112.6 vs 17.4 s) for +0.121 AP on ITT. The evaluation's cost is the bootstraps, not the scoring. <!-- calc: 112.59/17.42; SR; TR -->

---

## 4. Part 2 approach

### 4.1 Construction set A and lockbox L

Before any dependency value was read, the 1,345 non-engineered DepMap lines (U) were split into a lockbox L of 338 lines (25.1%) and a construction set A of 1,007. The split is PatientID-blocked iterative stratification on model labels only: one label per lineage plus 18 genotype labels, seed 3000. No patient straddles A and L. Runtime 0.6 s. <!-- MAN: n_U, n_L, n_A, seed, n_patients_straddling, runtime_s; results/part2/lockbox/README.md; calc: 338/1345 --> All construction, discovery and calibration used A only; L was read once per entry through a ledgered loader. <!-- MS: Methods "Firewalls, lockbox and contexts" --> L lines were in Part 1's development data and were processed with A in one release, so confirmation in L shows reproducibility within that release. <!-- SUP Note 15 item 26b --> The context table (MSI-H, functional LoF, hotspots, copy-number segments, fusions, sex, loss of Y, WGD, TP53, NRF2, lineage) was built in 86.5 s. <!-- results/part2/contexts/build_log.json: runtime_s; MS: Methods -->

### 4.2 Complex maps, selection and runtime

- **Score and graph.** zavg_BEF was recomputed on A over S, the 6,022 genes with effect < −0.5 in ≥ 1% of A lines. Mutual top-k edges above z0 were clustered by consensus Leiden-CPM over 20 PatientID-blocked 80% subsamples × 5 seeds, with persistence and sub-module splitting; rules were fixed on a toy simulation first (86 of 105 planted complexes recovered; null data gave no cluster). <!-- SUMC: zavg.genes_S; MS: Methods "Complex construction" (results/part2/construction/DECISIONS.md) -->
- **Selection.** A registered split-half statistic, which uses no labels, chose z0 = 6, k = 10 (0.343), but it is flat across z0 (0.341 at z0 = 4, 0.338 at z0 = 5). <!-- SEL: zavg_BEF.selected, zavg_BEF.grid -->
- **Comparator.** A GLS map with the same S, k and edge count per run. <!-- CON "Selection" -->

| | zavg_BEF map | GLS map |
|---|---|---|
| Level-1 complexes (genes) | 367 (3,642) | 379 (2,862) |
| Sub-modules | 132 | 108 |
| Median size | 9 | 5 |
| Locus-flagged, level 1 | 28.3% | 22.2% |
| T1 / T2 / T3-neg / T3-nc | 159 / 74 / 124 / 10 | 159 / 81 / 131 / 8 |
| Map-to-map stability (median best-match Jaccard; random sets 0.06 / 0.10) | 0.69 | 0.67 |
| Pairwise ARI between subsample maps | 0.754 | 0.747 |

<!-- SUMC: zavg.n_level1, zavg.genes_level1, zavg.n_level2, zavg.size_level1.median, zavg.locus_rate_level1, gls.*; CON "Maps", "Tiers"; TIERS: tiers.zavg.level1; VAR: zavg_BEF.size_bins.all.observed_pooled.median, gls.size_bins.all.observed_pooled.median, *.size_bins.all.random_S_pooled.median, *.pairwise.ari_both.median -->

Locus clusters (≥ 50% of members within 10 Mb or on one chromosome arm) are flagged and excluded from biological claims. <!-- MS: Results "A complex map, and a null primary bridge" --> "Complex" here means a supported co-dependency module, not a curated complex. The variability analysis shares subsamples with the frozen map, and its random reference is a floor, not a null. <!-- CON "Map-to-map variability"; MS: Discussion "Limitations" -->

**Runtime of construction** (laptop, under load from concurrent agents): <!-- CON "Runtimes" -->

| Step | Time |
|---|---|
| Toy decision validation | 2,529 s |
| 59 score runs: 29 zavg_BEF (A 159 s at 3.99 GB, 20 subsamples, 8 halves), 29 GLS (1–3 s each), 1 Pearson | 4,218 s (70.3 min) |
| Split-half selection (48 jobs, 3 processes) | 4,115 s wall |
| Map build zavg / GLS; tiers + bridge | 250.3 / 241.8 s; 52.0 s |
| Map variability (40 subsample maps) / FZ-X freeze | 325.1 s / 27.6 s |

<!-- CON "Runtimes" ("All 59 score runs"); results/part2/construction/runs/timings.jsonl (calc: sum runtime_s = 4218.1 over 59 records); results/part2/construction/{zavg_map,gls_map}/runtime.json: total_s; BR: timings.total_s; results/part2/construction/variability/maps_runtime.json: maps_wall_s; results/part2/fzx/FZX_manifest.json: runtime_s -->

### 4.3 Tiers

BioPlex 3.0 and CCLE proteomics were read only after the map commit (code guard). T1: within-cluster BioPlex edges above a degree- and essentiality-preserving null; T2: not T1, but protein co-abundance above a size-matched null; T3-neg: testable, not supported; T3-nc: not testable. <!-- CON intro; results/part2/verify/fzc_review.md §1 --> Among non-locus complexes, 155 of 263 (59%) are T1, against 151 of 295 (51%) for GLS. <!-- MS: Results (SUMC: zavg.tiers_nonlocus_level1, gls.tiers_nonlocus_level1) --> Complex names are post hoc HGNC gene-group annotations. Some level-1 complexes merge two complexes or are pathway modules.

### 4.4 The primary bridge (BioPlex)

The registered primary methods endpoint is size-standardised within-cluster BioPlex precision at matched gene coverage, with a protein-node bootstrap (B = 2,000, seed 5000). <!-- CON "Bridge" --> The zavg map dropped its 115 least stable complexes to match coverage. Precision was 0.157 for zavg and 0.182 for GLS: **−0.025 [−0.075, +0.031], null**. <!-- BR: primary.coverage_matching.trimmed.n_clusters_dropped, primary.zavg, primary.gls, primary.diff, primary.ci95 --> The interval allows an advantage of up to 0.031 for zavg or 0.075 for GLS: absence of evidence of a difference, not evidence of equivalence. <!-- SUP Note 15 item 5; RR1 M3.2 -->

| Sensitivity / secondary | zavg − GLS [95% CI] |
|---|---|
| No coverage matching | −0.020 [−0.054, +0.017] |
| Locus clusters excluded (pre-declared) | +0.003 [−0.060, +0.074] |
| Protein co-abundance (secondary) | −0.061 [−0.101, −0.018] (favours GLS) |
| Pair-level AP against BioPlex over S × S (no CI) | 0.0345 vs GLS 0.0215, Pearson 0.0223 |

<!-- CON "Bridge" sensitivities and secondary tables; BR: sensitivity.*, secondary.protein_coabundance, secondary.pair_level_AP_bioplex -->

An independent reviewer rebuilt both maps byte-identically and reproduced the bridge exactly (ACCEPT-WITH-NOTES). <!-- results/part2/verify/fzc_review.md "Verdict"; RT: Tree "FZ-C independent review" -->

### 4.5 Directions

Every direction is a registry entry (`docs/registry/*.yaml`); the engine refuses unregistered or modified entries and logs every run in a trial ledger. Each ran once on A, in registry order, after a committed intent. <!-- MS: Methods "Registry and reproducibility" --> Complex score: mean clipped CRISPR effect over fixed members, minus the context gene, its in-complex paralogues and nearby copy-number genes. Model: OLS, HC3, Satterthwaite df, PatientID-blocked Freedman–Lane permutation for screening hits (10,000; GPD tail only below 10⁻⁴). <!-- MS: Methods "Association statistics" -->

| Direction | Question | Tests (m) | Discoveries | Error control | Wall time |
|---|---|---|---|---|---|
| H3a lineage atlas | One-vs-rest lineage × complex | 8,441 | 1,145 in 21 of 23 families (942 non-locus) | Benjamini–Bogomolov selective FDR, q = 0.05 | 450.7 s |
| H3b driver/genome state | 20 contexts × complex | 7,340 | 161 in 13 families (121 non-locus) | BB selective FDR | 199.5 s |
| H1 within-complex GI | Does genetic loss of g change dependency on the rest of its complex? | 24 evaluable (registered 28) | 5 | BH | 53.4 s |
| H1S silencing (exploratory) | As H1, for expression-silencing events | 293 (registered 338) | 31 | BH | 437.4 s |
| H2 sub-module divergence (exploratory) | Does one sub-module carry the dependency beyond its siblings? | 452 in 210 families | 109 | Stage-wise | 57.1 s |

<!-- DR/H3a_lineage_atlas/summary.json: m_evaluable, n_discovery, n_discovery_nonlocus, n_families; MS: Results (trial_line.n_families_selected_perm); DR/H3b_driver_family/summary.json; DR/H1_within_complex_GI/runtime.json note; MS: Results H1 (DR/H1_within_complex_GI/summary.json); DR/H1S_silencing_GI/runtime.json note; RT: Tree "H1S"; DR/H2_submodule_divergence/summary.json (MS: m_evaluable, n_eligible_families, n_discovery); DR/*/runtime.json: wall_s -->

- **BB error control.** BB controls the within-family false-discovery proportion averaged over selected families. In 1,000 simulations its error rate was 0.046; the pooled proportion, 0.068, is reported but not claimed (D6). The guarantee assumes independence across families, which overlapping contexts do not satisfy. <!-- results/part2/verify/h3a_fdr_diagnosis.md "Validation"; RT: D6; SUP Note 15 item 11b -->
- **H1 classes.** Paralog compensation, residual sensitisation, collapse, or neutral (90% CI within ±0.1 Chronos units by TOST). A class also needs the effect outside a matched null of random member sets. <!-- MS: Results H1 -->

### 4.6 Controls and kill rule

The design registered 22 numbered positive controls (PC1–PC22). PC19 is registered as two entries (PC19a, PC19b), so there are 23 control entries, each mapped to a direction. They ran once on A, after the map and context freezes, in 28.0 s. <!-- docs/registry/README.md; RT: "Previous status (M2)"; CTRL/runtime.json: wall_s --> **Kill rule:** a direction is dropped if fewer than half of its evaluable complex-level controls are recovered. All 15 evaluable complex-level controls were recovered (one-sided p 3.3 × 10⁻⁸⁸ to 0.023): H3a 2 of 2, H3b 10 of 10, H1 3 of 3; H2 has none evaluable, so it is exploratory. <!-- MS: Results "Positive controls" (CTRL/kill_rule.json: *.evaluable, *.recovered, H2.verdict; CTRL/controls_summary.tsv) --> Coverage is thin (H3a rests on 2 controls; one multi-gene H1 readout, PC15), and the rule is lenient: it shows that the machinery recovers known biology, not that it has power for new effects. <!-- CTRL/kill_rule.json: H3a.thin_coverage; RT: Tree "Controls independent review"; RR1 minor 10 -->

### 4.7 Registered calibrations

These were implemented and run once each **after** the discovery runs and before any lockbox read (ruling R4). <!-- RT: D8 R4 --> Runtimes: H3a 35.7 s, H3b 86.3 s, H2 33.0 s. <!-- DR/*/calibration_D8.json: runtime_s -->

| Calibration | Result |
|---|---|
| Pseudo-contexts matched in size and lineage composition (fraction p < 0.05) | H3a 4.9%, H3b 4.5%; response threshold 0.075 not reached |
| Matched sets drawn from other complexes' members (10 per complex; size, essentiality and coherence; fraction p < 0.05, matched vs real) | H3a 28.1% vs 28.1%; H3b 11.6% vs 12.0% |
| Discoveries whose \|t\| beats all their matched sets | H3a 502 / 1,145; H3b 103 / 161 |
| Never-expressed gene sets, fraction p < 0.05 (expected 5%) | H3a 13.5% (melanoma, Skin, Lymphoid); H3b 6.2% |
| CN segment null | 6.0% (≤ 0.075; no response) |
| H2 random member splits | 78.3 rejections per replicate vs 109 real discoveries |

<!-- DR/H3a_lineage_atlas/calibration.json: frac_p_lt_005; DR/H3b_driver_family/README.md; DR/H3a_lineage_atlas/calibration_D8.json and DR/H3b_driver_family/calibration_D8.json: matched_random_sets.pooled.frac_p_lt_005, matched_random_sets.real_pooled.frac_p_lt_005, matched_random_sets.discoveries_exceeding_all_matched_sets, never_expressed_sets.pooled.frac_p_lt_005; RT: Tree "Registered calibrations" (H3b never-expressed 6.2%, CN segment 6.0%); DR/H2_submodule_divergence/calibration_D8.json: random_splits.stagewise_rejections_per_replicate.mean, random_splits.real.n_discovery -->

Pre-declared sensitivities (not trials): screen-quality covariates kept 843 of 1,145 H3a and 115 of 161 H3b discoveries (366.7 s and 151.2 s), with Esophagus/Stomach, WGD and MSI-H fragile; conservative p = max(p_perm, analytic) kept 1,113 and 130. In Lymphoid a complex's β correlates with its essentiality (ρ −0.240; −0.228 with screen-quality covariates), the scale effect that R-2 later addressed. <!-- DR/SENSITIVITY_SCREEN.md; DR/*/sensitivity_screen_quality/summary.json: runtime_s; DR/gpd_sensitivity_summary.json: *.n_conservative; DR/H3a_lineage_atlas/sensitivity_screen_quality/scale_effect_rho.tsv rows lineage:Lymphoid -->

### 4.8 Lockbox confirmation

The read came after a dry check that stopped on four uncoded steps (D9), amendments A3–A5, a pre-read code review (GO-WITH-NOTES), a 2.7 h dress rehearsal with a dummy L loader, and a committed intent. <!-- RT: D9, D10; LB/REHEARSAL.md "Result" (9,734 s) --> The 30 registered entries were then read once each, in registry order, with no re-read or resume (2.6 h, peak 3.67 GB). <!-- RT: D12; RT: Tree "Lockbox confirmations" --> **Confirmation rule:** same sign and one-sided p < 0.05 in L, for contexts with ≥ 5 context and ≥ 5 reference L lines. Rarer contexts take a labelled rare-context route (leave-one-lineage-out, DEMETER2 on A lines, leave-one-line-out), which is **not** an independent replication. <!-- MS: Methods "Lockbox confirmation, BENCH and BRIDGE-L" -->

| Direction (primary = read as executed, R7) | Confirmed / in lockbox | NOT EVALUABLE | Rate, evaluable | Rate, all-in |
|---|---|---|---|---|
| H3a | 601 / 1,145 | 0 | 0.525 | 0.525 |
| H3a, beats all matched sets / not | 287 / 502 vs 284 / 581 | 0 | 0.572 vs 0.489 | — |
| H3b | 60 / 150 | 28 (26 in WGD) | 0.492 | 0.400 (corrected relook 0.46) |
| H2 (exploratory) | 76 / 109 | 1 | 0.704 | 0.697 |
| H1S (exploratory) | 9 / 18 | 1 | 0.529 | 0.500 |
| Controls (23 entries; 44 readouts) | 29 / 44 | 13 | 0.936 | 0.659 |

<!-- RATES rows overall, beats_all_matched_sets, not_beats_all: primary_n_lockbox, primary_n_not_evaluable, primary_n_confirmed, primary_rate_evaluable, primary_rate_all_in; MS: Results "The lockbox read, and a defect" (CMP; LD: H3b.confirmed_per_context.WGD); RT: D14 (corrected 0.46) -->

- **Power aid** (post hoc, by the reviewer; not a test): 0.47–0.52 expected for H3a, 0.43–0.44 for H3b. **Confirmation is at the rate power predicts**, and so is the beats-matched-sets gap (0.52–0.53 vs 0.41–0.43). Gene-set-generic and scale effects would replicate too, so this cannot separate biology from reproducible artefacts. <!-- MS: Results "About half the associations confirm" (LD: power_band) -->
- **H1 in L.** ARID1A → ARID1B was confirmed (it is control PC6). EP300 → CREBBP passed only the rare route (4 L loss lines). FOCAD is defect-sensitive. TP53 was confirmed but cannot be promoted (E2). PTEN was not confirmed. <!-- MS: Results "H1 in the lockbox"; RT: D12 -->

### 4.9 BENCH and BRIDGE-L

- **BENCH (complex vs its best gene g\*)**. Over 1,295 discovery pairs, Δ = t_L(complex) − t_L(g\*).
  - High-breadth tercile (registered criterion): +0.240 [−0.033, +0.362]; 60% of pairs favour the complex.
  - The bootstrap is compromised: 16.5% of pair × replicate fits failed silently, 1,032 of 2,000 replicates used a reduced pair set, and with the covariate design assembled per replicate the lower bound is −0.122.
  - **Revised reading: inconclusive, with a point estimate favouring the complex; the registered interval is unreliable.** The first version's conclusion against the complex is withdrawn.
  - Low tercile: the best gene replicated better (−0.338 [−0.515, −0.201]), largely expected from how breadth is defined.
  - Runtime 8,192.7 s (B = 2,000).

  <!-- LBR: BENCH_complex_vs_gene.distinct.n_pairs; BENCH: primary.mean_high, primary.ci95_high, primary.mean_low, primary.ci95_low; RV: "60%"; MS (LD: BENCH.frac_pair_rep_nan, BENCH.n_rep_with_any_nan, BENCH.ci95_high_assembled); RR1 M4; LBR: BENCH_complex_vs_gene.runtime_s -->
- **BRIDGE-L** (the registered secondary bridge) asks whether zavg-map context discoveries confirm in L more often than discoveries from identical machinery on the GLS map. The answer is no.
  - R7 primary, evaluable hits (zavg 1,267, GLS 1,326): 0.514 vs 0.578, **−0.065 [−0.107, −0.025]**.
  - As executed: −0.072 [−0.111, −0.033]. Corrected relook: −0.069 [−0.108, −0.029].
  - Matched gene coverage: −0.074 [−0.122, −0.026]. T1/T2 clusters only: −0.074 [−0.120, −0.026].
  - Power explains 26–41% of the difference (post hoc, R-4; Section 7.2).
  - Runtime 733.8 s.

  <!-- BLR7: primary_evaluable.n_hits, primary_evaluable.primary.rate_zavg, .rate_gls, .diff, .ci95; as_executed.primary.diff, .ci95; corrected_nonconfirmatory.primary.diff, .ci95; primary_evaluable.multiverse.matched_gene_coverage, primary_evaluable.multiverse.T1_T2_clusters; REVR "R-4"; LBR: BRIDGE-L_confirmation_rate.runtime_s -->

### 4.10 Part 2 runtime summary

| Step | Wall time |
|---|---|
| Context table / lockbox draw / FZ-X | 86.5 s / 0.6 s / 27.6 s |
| Construction (score runs + selection + maps + validation) | ≈ 8,880 s (2.5 h) |
| Controls | 28.0 s |
| Discovery: H3a + H3b + H1 + H1S + H2 | 1,198 s (20 min) |
| Calibrations | 155 s |
| Screen-quality sensitivity | 518 s |
| Lockbox read (30 entries; dominated by the BENCH bootstrap) | 2.6 h |
| Per-hit robustness (strongest findings) | 117.5 s |
| Revision analyses R-1 to R-4 (Section 7.5) | ≈ 1,815 s |

<!-- REVR wall-time table (calc: 1797+1.5+14+2.7 = 1815.2); results/part2/contexts/build_log.json: runtime_s; MAN: runtime_s; results/part2/fzx/FZX_manifest.json: runtime_s; calc: 4218.1 + 4115 + 250.3 + 241.8 + 52.0 = 8877.2 (CON "Runtimes"); CTRL/runtime.json: wall_s; calc: 450.67+199.45+53.35+437.42+57.11 = 1198.0 (DR/*/runtime.json); calc: 35.7+86.3+33.0 = 155.0 (DR/*/calibration_D8.json); calc: 366.7+151.2 = 517.9 (DR/*/sensitivity_screen_quality/summary.json); RT: D12; ROB: runtime_s -->

The cost is dominated by the zavg_BEF score runs (the ARPACK embedding takes 118.6 of the 159 s per run on A) and the BENCH bootstrap. <!-- results/part2/construction/runs/timings.jsonl row set=A, method=zavg_BEF: stages_s.embedding_arpack, runtime_s; LBR -->

---

## 5. Best methods and strongest findings

### 5.1 Best method: the locked score at the pair level

zavg_BEF is the best pair-level co-complex ranker tested. It beat GLS, Spearman and Pearson on held-out hu.MAP 3.0 complexes under a rule frozen before the look, also on protein-new pairs, and ranked BioPlex pairs better than GLS and Pearson (AP 0.0345 vs 0.0215 / 0.0223; no CI). <!-- TR: summary; results.stratum_C3; BR: secondary.pair_level_AP_bioplex -->

Its scope is limited in three ways. The gain is carried by pair-rich complexes (Section 3.3). Its configuration was selected on hu.MAP training labels derived from Complex Portal, so agreement with curated complexes partly re-measures that selection. Preprocessing probably contributes part of the gain. <!-- MS: Discussion; DESIGN.md §0 F3 --> For complex construction, the zavg map's advantage is coverage (3,642 vs 2,862 genes), not better-supported complexes (Section 6). <!-- SUMC: zavg.genes_level1, gls.genes_level1 -->

### 5.2 The registered strongest findings (254)

**Rule** (DESIGN §2.5 and the registry): discovery q < 0.05; primary lockbox confirmation or a passed rare route; not defect-sensitive and not a locus cluster; leave-one-out robust (H3a by primary disease; H3b and H1 by lineage). <!-- ROBR: "Rule (applied mechanically)" --> The per-hit checks were registered but computed only after the lockbox read, on A only, against a candidate list committed first (fd03b03). <!-- ROBR: "Disclosure" -->

| Direction | Candidates | Strongest | Units | T1 / T2 / T3-neg |
|---|---|---|---|---|
| H3a | 524 | 207 | 137 | 124 / 56 / 27 |
| H3b | 46 | 45 | 41 | 27 / 9 / 9 |
| H1 | 2 | 2 (ARID1A → ARID1B, control PC6; EP300 → CREBBP, rare route) | 2 | 1 / – / 1 |
| **Total** | | **254** | | **152 / 65 / 37** |

<!-- ROB: per_direction.*.candidates, .strongest, .strongest_units, .strongest_per_tier; ANN: n, tier -->

- **Testability.** Leave-one-out applied fully to 199 findings and partly to 49; 6 passed the rare route instead (EP300 → CREBBP, five EWSR1-fusion findings). It is strict where one disease dominates a lineage (Bowel kept 0 of 29) and vacuous where the remainder has < 2 lines (Pancreas 17 of 17). <!-- MS: Results "The registered strongest findings" (ANN: loso_partially_testable; ROB); ROBR: "H3a LOSO behaviour" -->
- **Alternative registry readings** (post hoc, R-5). Only the two H1 calls could change, so the list would hold 252–254 findings; alternative rankings change only the order (Spearman 0.37 for |t| alone, 0.89 for |t| × tier weight). <!-- R5: A2_alternative_readings.loso_h1.total_alternative_bound, A2_alternative_readings.ranking.alt1_abs_t.spearman_vs_registered, A2_alternative_readings.ranking.alt2_abs_t_x_tierweight.spearman_vs_registered -->

**Post hoc annotations and labels** (none selects or ranks). 107 of 207 H3a and 34 of 45 H3b findings beat all their matched sets, a flag that no longer means "complex-specific"; 7 are readouts of registered controls. <!-- ANN --> The revised presentation labels (R-2, Section 7.2) replace the first version's Lymphoid/Myeloid-only artefact label (56) with one spline-trend rule in every lineage: 95 H3a findings are possible scale artefacts, and 37 are lead examples (T1, not a control readout, ≥ 2 residual SDs from the trend). For 26 of the 37 leads a scale artefact cannot be excluded; 9 are in less reliable lineages (Lymphoid, melanoma); DEMETER2 agrees in sign for 29 of 35. <!-- SD: possible_lineage_scale_artefact, revision.rev_scale_artefact_possible, revision.rev_lead, revision.rev_lead_concordant_nl, revision.rev_lead_less_reliable_lineage, revision.rev_lead_D2_same_sign, revision.rev_lead_D2_evaluable -->

Examples highlighted in the revised manuscript (A β; L β [95% CI], L lines; never-expressed rate; one-sided DEMETER2 p; negative = more dependent). The first four lie beyond the trend in its own direction (spline residual z −2.79 to −4.36), so DEMETER2 is their main independent support. <!-- SF rows as cited per line; SD: revision.text_highlighted_leads -->
- Ovary/Fallopian Tube, 17-member unit including CDH1, CTNNA1, PARD3, HNF1B and PAX8 (Z0035): −0.140; −0.202 [−0.298, −0.107], 16 lines; 0.22; p = 0.00041. The complex replicated better than its best gene, PAX8 (Δ +1.05). <!-- SF row Z0035; SD: revision.text_highlighted_leads."lineage:Ovary/Fallopian Tube|Z0035"; SUP Note 17 g* table -->
- Bone, SNAPC1–SNAPC5 (Z0303): −0.210; −0.251 [−0.370, −0.132], 15 lines; 0.19; p = 0.0089. <!-- SF row Z0303; SD: revision.text_highlighted_leads."lineage:Bone|Z0303" -->
- Head and Neck, 16-member unit including EGFR, GRB2, GAB1 and PTPN11 (Z0051): −0.133; −0.107 [−0.155, −0.059], 19 lines; 0.19; p = 0.00027. GRB2 alone replicated better (Δ −1.23). <!-- SF row Z0051; SD: revision.text_highlighted_leads."lineage:Head and Neck|Z0051"; SUP Note 17 g* table -->
- Myeloid, 12-member AP2 adaptor unit (Z0254): −0.189; −0.239 [−0.380, −0.098], 11 lines; 0.01; p = 0.00079. DNM2 alone replicated better (Δ −1.69). <!-- SF row Z0254; SD: revision.text_highlighted_leads."lineage:Myeloid|Z0254"; SUP Note 17 g* table -->
- Outside H3a, where no trend is defined: RB1-loss lines depend more on a 14-member cyclin E–CDK2/E2F unit (A −0.201; L −0.156 [−0.235, −0.076], 17 lines). <!-- SF row Z0069 -->
- Set aside: the lymphoid proteasome-associated finding (Z0190; residual z +5.92) lies in the trend's direction on a highly essential unit, in a less reliable lineage, and is weak in DEMETER2 (p = 0.14). <!-- MS: Results (S9 rows Z0190: rev_residual_z_nl, rev_D2_p_dir) -->

All functional readings are post hoc, from gene symbols, and the examples are a post hoc selection for presentation. These findings are replicated, robustness-checked module-level associations. Neither the matched-set comparison nor the lockbox shows that the complex, rather than one or two members, is the operative unit, and none has been tested outside DepMap. <!-- MS: Discussion "What survives" -->

---

## 6. What did not work

- **The primary (BioPlex) bridge is null** (−0.025 [−0.075, +0.031]), and protein co-abundance favoured GLS. Locus clusters explain the deficit only post hoc (+0.003 with them excluded). <!-- BR: primary.diff, primary.ci95, sensitivity.locus_excluded.diff -->
- **The secondary bridge (BRIDGE-L) is negative** (−0.065 [−0.107, −0.025]), also at matched coverage and among T1/T2 clusters, so locus clusters do not explain it. Power explains 26–41%; the remainder (−0.042) has no tested explanation. <!-- BLR7: primary_evaluable.primary, primary_evaluable.multiverse; R4: headline.obs_minus_exp; REVR "R-4" -->
- **Complex-level associations are shared by other coherent co-dependency modules** (28.1% vs 28.1% for H3a), including coherent non-member modules where they could be matched; essentiality-only sets associate less often (16.6%). A complex is therefore not more context-associated than a coherent, essentiality-matched module, though it is more so than an essentiality-matched gene set. Never-expressed sets show a lineage excess (13.5% vs 5%), pointing to a lineage shift of the CRISPR effect scale; most H2 divergence recurs in random splits (78.3 vs 109). <!-- DR/H3a_lineage_atlas/calibration_D8.json: matched_random_sets.*, never_expressed_sets.pooled.frac_p_lt_005; R1A rows H3a_lineage_atlas common_units_nbr, all_units nonmember_ess; DR/H2_submodule_divergence/calibration_D8.json -->
- **The "beats all matched sets" flag is mostly selection** (random-atlas discoveries 32.9% vs real 46.4%); across discoveries real complexes are not stronger, and the lockbox gap between flagged and unflagged discoveries is what power predicts. <!-- R1: directions.H3a_lineage_atlas.b.summary.random_symmetric_flag_rate, real_beats_all_rate; RV §4 -->
- **BENCH is inconclusive** (Section 4.9).
- **H1 recovers known biology.** ARID1A → ARID1B is control PC6; FOCAD → PELO–HBS1L overlaps control PC10; EP300 → CREBBP is the reciprocal of a published dependency. The neutral result is narrower than first stated: 11 of 24 losses changed the complex-wide dependency by less than ±0.1 Chronos units, none would be neutral at ±0.05 (minimal containing margin 0.059–0.097), and L cannot re-establish it (2 of 11 L 90% CIs inside ±0.1). <!-- RT: D11; PD: H1.n_neutral, H1.n_events; R5: A1_h1_margin_005.n_inside_005, A1_h1_margin_005.min_margin_range; LD: H1_lockbox.neutral_L_ci90_inside_margin -->
- **Not run** (declared before the read): BENCH within-line rank and the A–L random-set correlation; alternative per-hit scores; per-hit robustness for H1S and H2. <!-- AMD A3–A5; SUP Note 14 "(not run)" -->

---

## 7. Peer review and revision

### 7.1 History

| Round | Recommendation | Main points | Response |
|---|---|---|---|
| 1 | **Major revision** (numbers reproduce; issues interpretive) | M1 matched sets are other complexes' fragments, flag lacks a selection-aware reference; M2 scale artefact not carried to highlights; M3 bridges overstated, power one-sided; M4 BENCH bootstrap; M5 Part 1 wording; M6 H1 neutrality; M7 rule readings | R-1 to R-4 on A (intent 4157f1e → 740f439); claims narrowed; new title |
| 2 (77c86fe) | **Minor revision** (M1, M3, M5 resolved; M4 for reporting; M2, M6 partly; M7 mostly) | MUST: N1 H1 at ±0.05; N2 "reliable lineage"; N3 scope of the independent check. SHOULD: alternative readings, coherence attribution, atlas spread | R-5 from committed tables (4a7c715 → 5fac192); all MUST and SHOULD items addressed |

<!-- PR1 §2 "Recommendation", §3; PR2 header, §1, §2; RR1; RR2; RT: D18–D21; git log 4157f1e, 740f439, 4a7c715, 5fac192 -->

The title is now "A co-essentiality score selected on training complexes wins a single held-out test but yields neither better-supported complexes nor more replicable complex-level dependencies in DepMap"; "better-supported" because the primary bridge is null, not negative. <!-- MS title; RR1 M3.4 --> In neither round did a registered result, call, label or metric change, and no lockbox value, persisted read or test label was opened. <!-- RR1 intro; RR2 "Not changed" -->

### 7.2 Revision analyses (post hoc, descriptive; A only or committed tables)

| Analysis | Result | Reading |
|---|---|---|
| R-1(a) matched sets from the 2,380 genes outside every complex | Coherent non-member sets could be matched for only 90 of 367 units, the less coherent ones; there real 0.230, non-member 0.251 (H3a). Essentiality-only sets 0.166 vs real 0.281 (H3a), 0.088 vs 0.119 (H3b) | **Parity holds with coherent modules, not with any gene set.** Attributing the excess to coherence was dropped: coherence is measured on A, so matching on it partly conditions on the outcome |
| R-1(b) selection-aware "beats all" reference: 10 random atlases through the full discovery procedure | Real 0.464 vs random-atlas 0.329 (H3a; per atlas 0.291–0.364, 0 of 10 reach the real rate); 0.669 vs 0.483 (H3b; 0.448–0.517). P(\|t\| real > random discovery) 0.48 | **Selection alone lifts the rate from about 9% to 33%; the remaining excess is modest, within-group, and not an effect-size advantage.** The flag no longer means "complex-specific" |
| R-2 spline trend, one scale-artefact rule in all lineages | 95 possible artefacts (was 56); 37 leads (35 kept, 2 lost, 2 new); 26 of 37 beyond the trend in its own direction; 4 within 0.05 of the threshold | **The uniform rule cannot exclude an artefact for most leads.** Lymphoid (0.27) and melanoma (0.33) leads less reliable; others not thereby reliable (0.01–0.22) |
| R-3 DEMETER2 RNAi, 423 A lines | Leads: 29 of 35 same sign, 20 one-sided p < 0.05, ρ 0.86; all 207 H3a findings: 157 of 194 same sign | **Concordant for most leads**; six discordant leads named. Orthogonal in perturbation, not in lines |
| R-4 power aid applied to both maps' BRIDGE-L hits | Expected Δ −0.023 (SD 0.019) vs observed −0.065; share 26–41% across three power versions; remainder −0.042 (z −2.2, approximate: hits treated as independent); model under-predicts both levels (zavg by 0.055, GLS by 0.097) | **Power explains a quarter to two fifths of BRIDGE-L**; the rest is unexplained |
| R-5 (round 2) counts from committed tables | H1: 0 of 11 neutral events inside ±0.05 (minimal margin 0.059–0.097). Alternative readings: 252–254 findings; ranking changes order only. Atlas spread as in R-1(b) | **H1 neutrality rests entirely on the registered ±0.1**; the 254 are stable up to the two H1 calls |

<!-- R1A pooled rows; R1: directions.*.b.summary (real_beats_all_rate, random_symmetric_flag_rate[_all_tests], mean_p_beats_rand_disc_t_ctx); R2: scale_artefact_possible, artefact_D16, lead_corrected, transitions_lead, lead_corrected_concordant_nl, per_context; R3: n_A_models_with_D2, leads, all_207; R4: headline.*; R5: A1_h1_margin_005, A2_alternative_readings, A3_atlas_flag_range; REVR "four of the 37"; RVW "26%", "41%"; calc: 0.5137-0.4588, 0.5785-0.4817 -->

**Independent check of R-1 to R-4: PASS-WITH-NOTES, no MUST (19ad2ae).** Its scope is now stated in the Methods: the agent read the results summary first (not blind), reran scripts 59–61 and recomputed the numbers from committed tables with its own code, but did not rerun script 58 (random atlases, non-member test table) and reused rather than re-implemented the empirical-Bayes step of the power aid. Its SHOULD items were applied to the revision README. <!-- RVW "Verdict", "Independence caveat"; RT: D19; RR2 N3 -->

### 7.3 Main text changes

The BioPlex bridge is called primary (null) and BRIDGE-L secondary (negative) throughout; BENCH is inconclusive; "label-free" and "locked before the test data existed" are gone; the Holm bound is < 0.003; the H1 Discussion heading says "11 of 24 … and none within ±0.05"; "reliable lineage" is removed and the named leads carry their never-expressed rates. Fig. 6 and Supplementary Notes 17 (revision) and 18 (glossary) were added. The number check passes 421 of 421 in the main text and 320 of 320 in Notes 16–17 (952 of 1,208 over the whole supplement; Notes 1–15 predate the convention). Abstract 149 words; Introduction + Results + Discussion 4,805. <!-- RR1 M3, M4, M5, minor 17; RR2 N1, N2, N9, "Number checker on the supplement"; RT: D20, D21 -->

### 7.4 Requests not done, and why

| Request | Reason |
|---|---|
| Preprocessing-matched baseline; ΔAP leaving out each large test cluster (M5) | Needs the test labels again: a second look, which the protocol forbids |
| BENCH bootstrap with an explicit failure policy (M4.1) | Possible from the persisted, hash-checked L reads, but **we chose** not to make a second, post hoc look at L (round 1 wrongly called it impossible) |
| H1 margin in residual-SD units (M6) | Residual SDs were not exported (the ±0.05 part was done in R-5; round 1 wrongly called it not committed) |
| Exact H1 outcome under the primary-disease leave-one-out unit (M7) | Needs an A refit; bounded at 252–254 |
| Rescaling by never-expressed inflation (M2(d)) | Not computed; reliability column instead (round 2: "weak but tolerable") |
| Minimal detectable difference, primary bridge (M3.2); L effect of g\* alone (minor 9) | Not in committed outputs; CI bounds and Δ interpreted instead |
| Number check of supplement Notes 1–15 and the cover letter; some figure layout fixes | Notes predate the alias convention; layout handled in legends |

<!-- SUP Note 17 "Referee requests not carried out, and why"; RR1 "Summary of requests not carried out"; RR2 "M4.1", N1, SHOULD M7; PR2 §3 -->

### 7.5 Revision runtime (laptop)

Script 58 (R-1; 7 workers; it reran the full discovery procedure on ten random atlases) took 1,797 s at 3.36 GB; scripts 59, 60 and 61 (R-2 to R-4) took 1.5, 14 and 2.7 s; script 62 (R-5) reads committed tables and its runtime was not logged. <!-- REVR wall-time table; R5 (no runtime field) -->

---

## 8. Process and integrity

### 8.1 Amendments, errata and rulings, with timing

Outcome access (UTC): test labels 28 Oct 18:18; first A-side association outcome 29 Oct 00:20; H1 outcome 00:35; calibration outcomes from 01:03; lockbox read 06:50:05–09:23:05. <!-- SUP Note 14 header -->

| ID | What changed | When (UTC, 28–29 Oct) | Relative to outcome access |
|---|---|---|---|
| A1, A2 | Silencing events moved to exploratory H1S; genetic-loss rule (SMARCB1 into H1) | 22:54:50, 23:28:31 | Before any H1 run; labels only |
| R1/E1, R2 | H1 sign typo corrected to the meaning the engine used; H1S crash fixed, H1 not re-run | 00:41:42, 00:47 | After the H1 run |
| R3, R4 | GPD p stays primary (p_cons sensitivity); registered calibrations implemented and run once | 00:59, 01:02–01:06 | After discovery; before L |
| R5/E2 | "Detectable, within equivalence margin" label (TP53, PTEN); not promotable | 01:07 | After the H1 run; post hoc, conservative |
| A3–A5 | Lockbox seeds; BENCH secondaries not run; BRIDGE-L multiverse CIs; literal readings | 02:03:27–03:25:03 | Before L |
| R6, R7 | Failed fits NOT EVALUABLE; read as executed stays primary; 560-fit corrected relook alongside; defect-sensitive calls carry no claim | 10:19:44, 10:41:49 | After the read |
| D14–D16 | Reporting rules; robustness checks; 254 fixed; annotations and labels | After the read | A only |
| R-1 to R-4; R-5 | Revision analyses (round 1); round-2 counts | 13:43–14:50; 15:41–15:43 | After peer review; A only or committed tables |

<!-- SUP Note 14 table; SUP Note 7; AMD; LB/README.md (R6, R7 verbatim); MS: Discussion Table 3; git log 4157f1e, 740f439, 4a7c715, 5fac192 -->

The full record is in `manuscript/supplementary.md` Note 14 and `docs/registry/AMENDMENTS.md`. Rules fixed after outcomes were seen are listed as such: E2, the defect-sensitive label, the power aid, the robustness readings, the lead and artefact labels, and the revision analyses. Each makes the reading more conservative or is descriptive. <!-- MS: Discussion "Amendments and post hoc decisions" -->

### 8.2 The lockbox driver defect

- **The defect.** The driver built the L covariate design over all models and then subset it to L, so single-line lineages could form a spurious group. <!-- MS: Results "The lockbox read, and a defect" -->
- **Footprint.** 53 fits failed (50 from the defect; 3 genuinely: CTRL_PC21 ×2, H1S HNF1B), and 510 successful fits carried a spurious column (median |Δβ| 2.4 × 10⁻⁴, max 0.046). H3a is unaffected; 26 of the 28 failed H3b fits are in WGD. <!-- LB/README.md "PI ruling R7"; MS -->
- **Handling.** R7 keeps the read as executed as primary, with a corrected, non-confirmatory relook of all 560 affected fits (persisted matrices, resume path only). Calls that differ are defect-sensitive and support no claim: 53 (H3b 30, BRIDGE-L 10, controls 11, H1 FOCAD 1, H2 1). WGD and controls PC9, PC10, PC12 and PC13 lose their lockbox evidence. <!-- LB/README.md R7; MS (LD: defect_sensitive_total, defect_sensitive_by_entry) -->
- **How it was caught.** The dress rehearsal had shown the same failures, but nobody checked fit status. The 53 failed fits were found first, after the read (D12). The R6 reproduction check (68f0d3f) then found the second footprint, the 510 shifted fits, so the scope grew from 53 to 560. The results reviewer's label-only span audit predicted the 510/50/3 split exactly. <!-- RT: D12, D13, D14; RETRO "Stop" 2 (dated correction), 3 -->

### 8.3 Registry ↔ code gaps

Four registered steps had no code when the step they gate ran. <!-- RETRO "Stop" 1; LESSONS.md rule 15 -->
1. **Calibrations**: implemented after the discovery runs (R4), found by review 2 and the PI. <!-- RT: D8 R4; results/part2/verify/discovery_review.md "Review 2" -->
2. **Lockbox machinery**: the driver's own dry check stopped before any L read (D9); implemented as registered (A3–A5). <!-- RT: D9; LB/BLOCKED.md -->
3. **L covariate design**: built differently from the A side (Section 8.2); the pre-read review never executed the L path against the A path. <!-- RETRO "Stop" 2 -->
4. **Per-hit robustness**: found at manuscript time (D15) and computed after the read, on A only. <!-- RT: D15; results/part2/robustness/INTENT.md -->

An earlier mismatch (GPD tail rule with 1,000 permutations inflating H3a test size) was caught on synthetic data and fixed before any real association. <!-- results/part2/verify/h3a_fdr_diagnosis.md "Trigger", "Verdict" -->

During the revision, eight stale numbers in source records were found and corrected at the source (4db7dc7): the test scoring time, the score-run count, the control count, the FZ-X and context-build times, the retro's defect-discovery order, the robustness README's control-readout count, and script names in the code review. This report uses the corrected values. <!-- RT: D20 "Stale sources fixed (4db7dc7)"; git show 4db7dc7 -->

### 8.4 Independent reviews (agents that did not produce the work)

| Review | Verdict | What it caught or confirmed |
|---|---|---|
| Design rounds 1–3 (6 + 5 + 1 reviewers) | DESIGN v3.1 | Bridge CI by protein bootstrap; genotype lockbox strata; read order |
| Part 1 audit; test verifier | PASS; REPRODUCED (≤ 2.6 × 10⁻¹⁵) | No leakage (its own run caused D1); numbers committed before comparison |
| FZ-C; FZ-X; controls reviews | ACCEPT(-WITH-NOTES) | Byte-identical rebuild; SMARCB1 typing (A2); lineage-merge convention |
| Discovery reviews 1 and 2 | ACCEPT-WITH-NOTES | 15,781 tests reproduced; GPD-tail dependence; Lymphoid scale effect |
| Pre-read lockbox code review | GO-WITH-NOTES | Missed the covariate defect |
| Lockbox results review | PASS | 54/54 calls recomputed; 510/50/3 predicted; 16.5% silent BENCH failures; power aid |
| Reference check | E1–E15 applied | 9 metadata errors, 2 claim mismatches |
| Peer review rounds 1 and 2 | MAJOR, then MINOR REVISION | Section 7.1 |
| Revision-analyses check | PASS-WITH-NOTES | Not blind; script 58 and the EB step not rerun |

<!-- RT: "Previous status (M2)", D1, D5, D14, D17, D18, D19, D21, Tree entries; docs/reviews/round3_final.md (DESIGN.md header); TV; results/part2/verify/{fzc_review,fzx_review,controls_review,discovery_review,lockbox_code_review,lockbox_results_review}.md "Verdict"; RETRO "Stop" 5 (a303165); PR1 §2; RVW "Verdict"; PR2 §1 -->

Reviewers, the peer reviewer included, and producers share a model family. Reviews reproduced every confirmatory number, but they are not independent in the human sense. <!-- MS: Discussion "Limitations"; PROMPT.md ("peer reviewer (Opus 5.5)") -->

---

## 9. Limitations

1. **Procedural single look.** The test labels were plain files, and unsupervised scores already covered all test pairs. D1 shows the risk was real. <!-- SUP Note 15 item 1 -->
2. **Narrow test win.** The gain over GLS is carried by pair-rich clusters, with an effective cluster count of 21.3. The training labels derive from Complex Portal, so the test partly re-measures selection. Preprocessing probably contributes, and a matched baseline can no longer be tested on this test set. <!-- TR: results.itt.cluster_stats.kish_neff_clusters; SUP Note 15 items 2–3; MS: Discussion -->
3. **Calibration.** CITL was 0.893, with no pre-specified CI. Threshold and calibrator sampling error is not propagated. <!-- TR: results.itt.metrics.zavg_BEF.calibration.isotonic.calibration_in_the_large; SUP Note 15 item 4 -->
4. **Map choices.** Selection is flat across z0; 28% of complexes are locus clusters; some merge two complexes or pathways; about 40% of S genes are in no complex. <!-- SEL; SUMC: zavg.locus_rate_level1; SUP Note 15 items 6–10 -->
5. **Late checks.** Calibrations and robustness checks ran late (R4, D15). <!-- RT: D8, D15 -->
6. **Lockbox.** L (338 lines) has been read; non-confirmation of small effects is often uninformative (325 H3a discoveries have |β| < 0.05); the rare route is not independent replication; L shares Part 1's development data and the release's processing with A; 53 calls are claimless. **No held-out lines remain for a new confirmatory test.** <!-- MAN: n_L; MS (PD: H3a.n_abs_beta_lt_0.05); SUP Note 15 item 26b; LD: defect_sensitive_total -->
7. **Lineage-scale artefacts.** The rule cannot flag effects beyond the trend in its own direction (26 of 37 leads); never-expressed sets reach p < 0.05 at up to 33% in some lineages, and statistics were not rescaled; the |z| ≥ 2 threshold is unregistered. <!-- R2: lead_corrected_concordant_nl, per_context.melanoma.never_expressed_frac_p_lt_005; SUP Note 15 item 16 -->
8. **Comparators.** Coherent non-member sets matched only a quarter of the units, and coherence is measured on A, so matching on it partly conditions on the outcome. <!-- REVR "R-1"; RR2 N4 -->
9. **Partial answers.** BENCH is inconclusive; H1 neutrality rests on the loose ±0.1 margin; most of BRIDGE-L is unexplained by power. <!-- BENCH; R5: A1_h1_margin_005; REVR "R-4" -->
10. **One resource.** All dependencies come from DepMap 26Q3 CRISPR (Chronos); DEMETER2 was used on A lines only; nothing was tested in another collection or in vivo. Controls are thin (H3a 2, H2 none). <!-- SUP Note 15 item 27; R3; CTRL/kill_rule.json -->
11. **Unverified metadata.** Licences (hu.MAP test files, CCLE proteomics, BioPlex), the HGNC download date and the DepMap release need confirmation (Section 10). <!-- DATA.md; RR1 "Reproducibility" -->
12. **Agent-run study** under a human PI. The lockbox loader's run log is permanently lost (loader count inferred from the ledger). Only the main text and Supplementary Notes 16–17 are fully machine-checked. <!-- RV §1; SUP Note 15 item 26; RT: D21 -->

---

## 10. Follow-ups, remaining ideas and open items

**Placeholders the human authors must fill in before submission** (agents cannot resolve these): <!-- RR1 "Items blocking submission"; RR2 "Not changed"; RT: D18 -->
- [AUTHORS], [CORRESPONDING_AUTHOR], [COMPETING_INTERESTS];
- [ARCHIVE_DOI] for this repository (https://github.com/cancergenetics/ai-hackathon-Coessentiality_Prediction). Part 1 is public at https://github.com/cancergenetics/ai-hackathon-Coessentiality; its commit 3c684c9 is content-identical to the locked c7f95ee for every file used here (`external/README.md`);
- [DEPMAP_RELEASE_CONFIRMATION];
- licence checks: [LICENCE_HUMAP], [LICENCE_CCLE_PROTEOMICS], [LICENCE_BIOPLEX];
- [HGNC_DOWNLOAD_DATE].

**Scientific follow-ups.** L has been read and the hu.MAP test set used, so any new confirmatory claim needs its own registration and new held-out data; work on A or L from here is exploratory. Ordered by expected gain in confidence, then novelty:

1. **Lineage-scale normalisation (highest priority).** Re-express complex scores relative to the lineage's essentiality trend before testing, or use never-expressed sets as an empirical null to rescale statistics per context (the M2(d) request). The matched-set test needs far more than 10 sets per complex to give a usable empirical p-value. <!-- ANN; R2; SUP Note 15 item 11 -->
2. **Orthogonal replication with data already provided.** R-3 checked DEMETER2 on A lines only; a registered RNAi replication of the 254 findings and the 11 neutral H1 events on L lines with RNAi data would be the cheapest independent perturbation test. CCLE proteomics could test whether contexts change complex protein abundance. <!-- DATA.md; R3: n_A_models_with_D2 -->
3. **External CRISPR replication.** A second screen collection, requested from the user (PROMPT.md), not downloaded.
4. **Explain the bridge results.** Locus-aware construction (copy-number-corrected profiles, or no arm-block edges), then both bridges again; the unexplained BRIDGE-L remainder (−0.042) and the larger under-prediction of GLS's level; a BioPlex-integrated map validated on held-out BioPlex or proteomics (D4); and the untested idea that any clustering captures the gain on pair-rich complexes while extra coverage adds weaker modules. <!-- RT: D4; MS: Discussion; R4 -->
5. **Declined revision requests that remain possible:** the BENCH bootstrap with a failure policy, the H1 residual-SD margin and the H1 primary-disease leave-one-out refit (Section 7.4), each as a disclosed post hoc analysis.
6. **Directions cut at D4:** H4 copy-number dosage (CYCLOPS-like, hard because of locus artefacts) and rewiring (context-specific coherence; BioPlex 293T vs HCT116). <!-- RT: D4; DESIGN.md §2.6 "Cut" -->
7. **Declared not-run items:** alternative per-hit scores; BENCH within-line rank and the A–L random-set correlation; per-hit robustness for H1S and H2. <!-- SUP Note 14 "(not run)" -->
8. **Pair-level method.** zavg_BEF does not beat GLS on small complexes (≤ 5 positive pairs); a preprocessing-matched GLS is the first comparison, on a new held-out benchmark. <!-- EXP -->
9. **Process:** a registry↔code traceability table that fails on unmapped steps; rehearsals that assert fit-status counts and A-path = L-path; committed run logs for irreversible runs; the alias convention and number check in every document from the start; checking committed files before writing "not done" reasons (round 1 gave two wrong ones). <!-- RETRO "Proposed pruned LESSONS.md" 11, 13–15; RT: D21; RR1 M4.1, M6 annotations -->

---

## 11. Reproducibility

**Environment.** Run `scripts/00_setup_part1.sh` to clone the locked Part 1 repository (c7f95ee) into `external/part1_locked` and build `.venv`. Set `PART1_DIR` and `PYTHONHASHSEED=0` for every step. `scripts/check_data.sh` verifies every dataset checksum against DATA.md. <!-- scripts/00_setup_part1.sh; scripts/check_data.sh; SUP Note 9 -->

**Pipeline, in order** (each step's exact command is in its intent note):

| Stage | Scripts |
|---|---|
| Part 1 test evaluation | `19`–`23` (pair files without labels → scoring → fit → `22_freeze_protocol.py` → single look → descriptive export); verifier `results/test_eval/verify/verify_independent.py` |
| Part 2 setup and construction | `05`, `32`, `33` (contexts, lockbox) → `30`, `35` (toy) → `34_run_A_blocks.sh` → `36`–`39` (maps, validation after the map commit, variability) |
| Freeze, A-side analyses | `41`–`43` (FZ-X, A1, A2) → `34_run_controls.py` → `44_run_direction.py` (H3a, H3b, H1, H1S, H2) → `45`–`49` |
| Lockbox | `52`–`54` (A3–A5) → `51_lockbox_confirm.py` (single read; ledger blocks re-reads; `--relook R6` resume-only) → `55_lockbox_r7_report.py` |
| Strongest findings; revision | `56`, `57`; `58_rev_matched_sets.py` → `59`–`61` → `62_rev_round2.py` (committed tables only) |
| Figures, tables, checks | `tab_*.py`, `fig_1`–`fig_6_*.py`, `check_manuscript_numbers.py` (`--start-at`), `make_submission_copy.py` |

<!-- scripts/ listing and docstrings; CON "Pipeline and files"; ROBR "Reproduce"; REVR "Reproduce"; SUP Note 9; LB/CORRECTED_REFIT_INTENT.md "Commands"; RR2 "Number checker" -->

**Seeds**: test bootstrap 0 (MC 1–4, spot-check 7); lockbox 3000; Leiden 0–4; subsamples 1000–1019; split halves 2000–2003; permutations 4000+ (controls 4400 + index); BRIDGE-L GLS discovery 4500; revision sets 4600 / 4700 and atlas permutations 4000 / 4100; bridge bootstrap 5000, co-abundance 5001, BENCH 5100, BRIDGE-L 5200; variability reference 7000. <!-- DESIGN.md §2.1; MAN: seed; CON; BR; MS: Methods; AMD A3; TV; REVR "Seeds" -->

**Data.** DATA.md lists source, version, SHA-256 and licence for every input; the DepMap and HGNC checksums equal Part 1's. **Never committed:** `datasets/` (including the hu.MAP test files), `data_cache/` (score blocks, context tables and the persisted lockbox reads, whose hashes are in `results/part2/lockbox_reads.jsonl`), `external/`, `.venv/` and binary arrays. The ledger blocks any re-run of the lockbox read. <!-- DATA.md; .gitignore; LB/README.md "Retention rule"; SUP Note 9 -->
