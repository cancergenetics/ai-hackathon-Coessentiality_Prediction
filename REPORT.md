# REPORT: co-essential protein complexes, held-out test and cancer-context dependencies

<!-- Source convention (as in manuscript/main.md): every number is followed by an HTML comment naming the committed
file it comes from, with JSON key paths, TSV rows/columns, or "calc:" for arithmetic on cited numbers. Where a number
is quoted from the manuscript, the comment names the manuscript's own source (checked 579/579 by
scripts/check_manuscript_numbers.py, RESEARCH_TREE.md D17). Aliases:
TR   = results/test_eval/test_results.json
EXP  = results/test_eval/export/ap_by_size.tsv
SR   = results/test_eval/score_records/test/<method>.run.json
TV   = results/test_eval/verify/README.md
RT   = RESEARCH_TREE.md
MS   = manuscript/main.md
SUP  = manuscript/supplementary.md
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
AMD  = docs/registry/AMENDMENTS.md
RETRO= docs/retros/retro_M6_final.md
-->

Status: draft for PI review, not committed. The submission manuscript is `manuscript/main.md` (with `manuscript/supplementary.md`); this report summarises the project for the record and does not add any analysis.

---

## 1. Executive summary

The locked Part 1 co-complex score, zavg_BEF, was evaluated once on the withheld hu.MAP 3.0 test pairs. It reached average precision (AP) 0.333 [95% CI 0.256, 0.410], against 0.212 for GLS, 0.204 for Spearman and 0.197 for Pearson (ΔAP vs GLS +0.121 [+0.079, +0.159]). <!-- TR: results.itt.bootstrap.cluster.stat.*|ap; results.itt.bootstrap.cluster.diff."zavg_BEF|ap - gls_signed|ap" --> The pre-registered verdict passed. The gain comes from pair-rich test clusters: macro per-cluster AP is no better than GLS's (0.297 vs 0.292). <!-- TR: results.itt.metrics.*.macro_cluster_ap -->

In Part 2 the score built a curation-free complex map of 367 complexes. <!-- SUMC: zavg.n_level1 --> Neither pre-registered bridge held: BioPlex agreement was null (−0.025 [−0.075, +0.031]), and zavg-map context discoveries confirmed in held-out lines less often than GLS-map ones (−0.065 [−0.107, −0.025]). <!-- BR: primary.diff, primary.ci95; BLR7: primary_evaluable.primary.diff, primary_evaluable.primary.ci95 --> Matched random gene sets associated with lineages as often as real complexes (28.1% vs 28.1%). <!-- DR/H3a_lineage_atlas/calibration_D8.json: matched_random_sets.pooled.frac_p_lt_005, matched_random_sets.real_pooled.frac_p_lt_005 --> Lockbox confirmation (H3a 52.5%) matched what power predicts. <!-- RATES row H3a_lineage_atlas overall: primary_rate_evaluable; LD: power_band --> Complex scores did not out-replicate their best gene. A registered rule leaves 254 strongest findings. <!-- ROB: per_direction.*.strongest --> A driver defect, found after the lockbox read, strips claims from 53 calls. <!-- LB/README.md "PI ruling R7"; LD: defect_sensitive_total -->

---

## 2. Approach, Part 1 → test transition

**What was locked.** The Part 1 repository was locked at commit c7f95ee (28 Oct 2026, 12:55:43 UTC) and used from a read-only clone (`external/part1_locked`, created by `scripts/00_setup_part1.sh`). <!-- DATA.md "Part 1 (locked method)"; TR: part1.commit, part1.matches --> zavg_BEF is the equal-weight mean of three component z-scores (first-order correlation F, second-order neighbourhood cosine B, Laplacian-eigenmap cosine E) on gene-standardised CRISPRGeneEffect profiles with the top 20 principal components removed and ridge whitening. It has no parameters fitted to labels. <!-- MS: Methods "The locked method" (LITERATURE.md §0) --> The baselines are Part 1's development-selected signed Pearson, Spearman and GLS. B alone was a pre-declared secondary candidate. A PC-removed Pearson reference (`spec_k20_tinf_std`) was added post hoc and is not a candidate. Nothing learned in this project changed the method.

**Transition and audit.** The hu.MAP test files appeared at 13:55 UTC, 59.6 min after the final Part 1 commit. <!-- DATA.md "Notes"; SUP Note 1 (docs/part1_audit.md "Transition timeline") --> An independent audit (`docs/part1_audit.md`) found no evidence that test pairs, labels or derived information entered Part 1 (PASS on timeline, references, content and hyper-parameters). The test set is pair-disjoint from training (0 overlaps in all four positive/negative combinations). It has 1,121 proteins, 329 of them also in training. Of the 3,312 positives, 2,309 have neither protein among the training-positive proteins (stratum C3). <!-- RT: Tree "Part 1 audit"; TR: audit.train_overlap, audit.proteins_analysed, audit.proteins_also_in_train, audit.c_strata.C3.pos, audit.analysed_pos -->

**Disclosed contacts before the freeze.** <!-- SUP Note 1 table (docs/TEST_PROTOCOL.md §0); RT: D0, D1 -->
- **D0 / T0.** The PI read the line counts (3,312 positive and 597,426 negative lines) and the first 3 lines of each file as a format check. No scoring. <!-- RT: D0; DATA.md "Notes" -->
- **D1.** In an intermediate run, the Part 1 audit agent counted how many test pairs appear in Part 1's label-free top-pair table (zavg_BEF ≥ 3). That is equivalent to precision and recall at the z ≥ 3 rule. The counts were deleted and the PI never saw them. Rule (a) (z ≥ 3) is therefore reported as non-blind and enters no verdict. <!-- RT: D1; MS: Results "The locked score wins" -->
- **T1–T3.** Label-free descriptive counts by the audit agent (T1); an octal dump of the first ~200 bytes of the positive file (T2). Genome-wide label-free scores already covered all test pairs (T3), so the single-look guarantee is procedural, not cryptographic. <!-- SUP Note 1 table -->

**Protocol.** `docs/TEST_PROTOCOL.md` v1.2 went through three design-review rounds (6, 5 and 1 reviewers). <!-- RT: "Previous status (M2)" -->
- Primary metric: non-interpolated AP at natural prevalence on ITT (all 600,738 canonical pairs; the 33,931 unmapped pairs sit at a tied floor score for every method). <!-- TR: audit.unique_canonical_pairs, results.itt.n_unmapped_any -->
- The locked "clearly beats" verdict (Part 1's rule) is computed on S1, the pairs mapped under every method (3,159 positives), with ITT quoted alongside (D3). <!-- TR: results.mapped_only.n_pos; RT: D3 -->
- Uncertainty: paired cluster bootstrap over 238 test clusters, B = 2,000, seed 0; protein and pair bootstraps are secondary. <!-- TR: results.itt.cluster_stats.n_clusters, results.itt.bootstrap.cluster.n_boot, results.itt.bootstrap.cluster.seed -->
- Operating points (rule b = F1-maximising threshold on all hu.MAP training pairs) and isotonic/Platt calibrators were fitted on the 634,779 training pairs before the freeze and applied unchanged. <!-- MS: Methods "Calibration"; results/test_eval/fitted_rules.json -->
- A dry run before the freeze reproduced Part 1's internal final test to ≤ 1.5 × 10⁻⁸. <!-- MS: Methods (docs/TEST_PROTOCOL.md §10) -->

**The single confirmatory look.** Freeze commit 039a843 with a marker hashing the protocol, fitted rules, evaluation code, tests, package list and test files (marker SHA-256 524c4d68…). Intent ab674b7; "started" record 7afcf06, committed before labels were read; results 0861e59. <!-- TR: protocol.commit, protocol.marker_sha256, protocol.head, started_commit; RT: Tree "TEST_PROTOCOL v1.2 frozen" --> The harness refuses to run if any earlier look exists; `prior_look_evidence.any` is false and the run is marked confirmatory. <!-- TR: prior_look_evidence.any, confirmatory --> A verifier that had not written the harness re-implemented the protocol and committed its numbers (c69d828) before opening the harness output. Every verdict-relevant number matched to ≤ 2.6 × 10⁻¹⁵. 36 secondary calibration fields differ by ≤ 4.1 × 10⁻⁷, explained by the ECE bin-remainder rule and optimiser tolerance. <!-- TV: "Verdict", "Order of work", "Discrepancies" --> A later descriptive export (curves and size bins) is logged as non-confirmatory and reproduces the look's AP to ≤ 10⁻¹⁶. <!-- results/test_eval/export/README.md "Reproduction of the look"; results/test_eval/looks.jsonl (descriptive_export) -->

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

**S1 (locked-verdict set; pairs mapped under every method; 3,159 positives, 563,648 negatives)** <!-- TR: results.mapped_only.n_pos, results.mapped_only.n_neg -->

| Method | AP | AUROC | Precision | Recall | F1 | Brier (×10⁻³) | Cal. slope | CITL |
|---|---|---|---|---|---|---|---|---|
| **zavg_BEF** | **0.349** [0.274, 0.435] | 0.815 [0.768, 0.853] | 0.553 [0.460, 0.647] | 0.343 [0.269, 0.429] | 0.424 [0.358, 0.492] | 4.25 [3.25, 5.43] | 1.017 | 0.914 |
| GLS | 0.222 [0.163, 0.291] | 0.750 [0.702, 0.800] | 0.567 [0.480, 0.653] | 0.215 [0.163, 0.284] | 0.311 [0.250, 0.381] | 4.73 [3.60, 6.06] | 1.006 | 0.934 |
| Spearman | 0.214 [0.113, 0.369] | 0.779 [0.722, 0.844] | 0.306 [0.212, 0.443] | 0.286 [0.139, 0.474] | 0.295 [0.184, 0.418] | 5.03 [3.97, 6.41] | 0.850 | 1.229 |
| Pearson | 0.206 [0.108, 0.352] | 0.771 [0.714, 0.839] | 0.273 [0.194, 0.387] | 0.296 [0.152, 0.482] | 0.284 [0.184, 0.399] | 4.98 [3.97, 6.28] | 0.879 | 1.198 |

<!-- TR: results.mapped_only.bootstrap.cluster.stat.*; results.mapped_only.metrics.*.calibration.isotonic -->

Reading:
- zavg_BEF has the highest AP, AUROC, recall and F1 at rule (b). Its precision at rule (b) is about the same as GLS's (0.553 vs 0.567), but it recalls more positives at that threshold. <!-- TR: results.itt.bootstrap.cluster.stat.*_rule_b -->
- zavg_BEF's isotonic calibrator kept its slope (1.007) but under-predicted the total positive mass (CITL 0.893). No recalibration was made. The Spearman and Pearson calibrators over-predicted (CITL 1.195, 1.166). <!-- TR: results.itt.metrics.*.calibration.isotonic -->
- Secondary candidate B: AP 0.310 [0.234, 0.392]. It failed the locked rule against GLS (macro-AP ratio 0.875 < 0.95 on S1) and was below zavg_BEF (ΔAP −0.022 [−0.032, −0.012]). <!-- TR: results.itt.bootstrap.cluster.stat."B|ap", results.itt.bootstrap.cluster.diff."B|ap - zavg_BEF|ap", summary.locked_verdict_S1.B; TV: "Locked verdict" -->
- The post hoc PC-removed Pearson reference reached AP 0.239, below zavg_BEF (−0.094 [−0.126, −0.061]). <!-- TR: results.itt.metrics.spec_k20_tinf_std.ap, results.itt.bootstrap.cluster.diff."spec_k20_tinf_std|ap - zavg_BEF|ap" -->

### 3.2 Paired differences (zavg_BEF − baseline)

| Contrast | ITT ΔAP [95% CI] | Relative gain | S1 ΔAP [95% CI] | ITT ΔAUROC [95% CI] |
|---|---|---|---|---|
| vs GLS | **+0.121** [+0.079, +0.159] | +57% | +0.127 [+0.089, +0.167] | +0.059 [+0.030, +0.090] |
| vs Spearman | +0.129 [+0.023, +0.216] | +63% | +0.135 [+0.024, +0.229] | +0.033 [−0.003, +0.067] |
| vs Pearson | +0.136 [+0.039, +0.220] | +69% | +0.143 [+0.038, +0.233] | +0.040 [+0.002, +0.076] |

<!-- TR: results.itt.bootstrap.cluster.diff."zavg_BEF|ap - *|ap" (diff, ci_low, ci_high, rel_gain); results.mapped_only.bootstrap.cluster.diff."zavg_BEF|ap - *|ap"; results.itt.bootstrap.cluster.diff."zavg_BEF|auroc - *|auroc" -->

- **Locked verdict:** "clearly beats" holds against all three baselines on S1 and on ITT. <!-- TR: summary.locked_verdict_S1.zavg_BEF, summary.locked_verdict_ITT.zavg_BEF -->
- Holm-adjusted p for the three primary ΔAP comparisons: < 0.0005 (GLS), 0.022 (Spearman), 0.018 (Pearson). <!-- SUP Note 10 "Inference" (TR: holm_adjusted_p_primary) -->
- The ΔAUROC CI against Spearman includes 0. <!-- TR: results.itt.bootstrap.cluster.diff."zavg_BEF|auroc - spearman_signed|auroc" -->
- The CIs are wide because the 238 clusters have a Kish effective number of only 21.3. <!-- TR: results.itt.cluster_stats.kish_neff_clusters -->
- On protein-new pairs (C3), AP was 0.420 against 0.272 for GLS (+0.148 [+0.093, +0.199]). <!-- TR: results.stratum_C3.bootstrap.cluster.stat, results.stratum_C3.bootstrap.cluster.diff."zavg_BEF|ap - gls_signed|ap" -->

### 3.3 Macro-AP and AP by cluster size

Macro per-cluster AP weights each test cluster equally. It does not separate zavg_BEF from GLS: 0.297 vs 0.292 (+0.005 [−0.014, +0.047]). <!-- TR: results.itt.metrics.*.macro_cluster_ap; results.itt.bootstrap.cluster.diff."zavg_BEF|macro_cluster_ap - gls_signed|macro_cluster_ap" --> Without the 10 clusters with most positives, AP was 0.345 vs 0.264. <!-- TR: results.itt.metrics.*.ap_no_top10_clusters -->

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
| Evaluation: load 8.1 s; ITT metrics + all bootstraps 1,344.1 s; S1 483.8 s; strata C1–C3 161.8 s | 1,998 s (33.3 min) | – |
| of which ITT cluster bootstrap / protein / pair / MC seeds 1–4 | 758.4 / 129.7 / 127.7 / 310.1 s | – |
| Descriptive export | 17.1 s | – |

<!-- SR (zavg_BEF, gls_signed, spearman_signed, pearson_signed, B, spec_k20_tinf_std): runtime.total_s, runtime.scoring.components_s, runtime.peak_rss_gb; TR: runtime_s (calc: 8.1+0.6+1344.1+483.8+6.4+65.7+89.7 = 1998.4; 6.4+65.7+89.7 = 161.8), results.itt.bootstrap.runtime_s; results/test_eval/export/README.md -->

Accuracy for cost: zavg_BEF costs about 6.5× GLS's scoring time (112.6 vs 17.4 s) for +0.121 AP on ITT. The evaluation's cost is the bootstraps, not the scoring. <!-- calc: 112.59/17.42; SR; TR -->

---

## 4. Part 2 approach

### 4.1 Construction set A and lockbox L

Before any dependency value was read, the 1,345 non-engineered DepMap lines (U) were split into a lockbox L of 338 lines (25.1%) and a construction set A of 1,007. The split is PatientID-blocked iterative stratification on model labels only: one label per lineage plus 18 genotype labels, seed 3000. No patient straddles A and L. Runtime 0.6 s. <!-- MAN: n_U, n_L, n_A, seed, n_patients_straddling, runtime_s; results/part2/lockbox/README.md; calc: 338/1345 --> All construction, discovery and calibration used A only. L was read once per registry entry through a ledgered loader. <!-- MS: Methods "Firewalls, lockbox and contexts" --> The context table (MSI-H, LoF with a functional check, hotspots, copy-number segments, fusions, sex, loss of Y, WGD, TP53, NRF2, lineage) was built from DepMap Omics and Model.csv in 86.5 s. <!-- results/part2/contexts/build_log.json: runtime_s; MS: Methods -->

### 4.2 Complex maps, selection and runtime

- **Score.** zavg_BEF was recomputed on A over S, the 6,022 genes with effect < −0.5 in ≥ 1% of A lines. On all U lines the recomputation reproduces Part 1 to ≤ 10⁻¹². <!-- SUMC: zavg.genes_S; MS: Methods "Complex construction" -->
- **Graph and clustering.** Mutual top-k edges above z0 were clustered by consensus Leiden-CPM over 20 PatientID-blocked 80% subsamples × 5 Leiden seeds, with multi-resolution persistence and top-down sub-module splitting. The rules were fixed on a toy simulation first: 86 of 105 planted complexes and all 39 planted sub-modules were recovered, and null data gave no cluster. <!-- MS: Methods "Graph and clustering" (results/part2/construction/DECISIONS.md) -->
- **Selection.** A label-free split-half statistic chose z0 = 6, k = 10 (0.343). The criterion is flat across z0 at k = 10 (0.341 at z0 = 4, 0.338 at z0 = 5), so z0 is not well determined. A degree-rewired graph scores 0. <!-- SEL: zavg_BEF.selected, zavg_BEF.grid; CON "Selection" -->
- **Comparator.** A GLS map with the same S, the same k and the same edge count per run. <!-- CON "Selection" -->

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

Locus clusters (≥ 50% of members within 10 Mb or on one chromosome arm) are flagged and excluded from biological claims. <!-- MS: Results "A curation-free complex map" --> The variability analysis shares subsamples with the frozen map, and its random reference is a floor, not a null. <!-- CON "Map-to-map variability" -->

**Runtime of construction** (laptop, under load from concurrent agents): <!-- CON "Runtimes" -->

| Step | Time |
|---|---|
| Toy decision validation | 2,529 s |
| Score runs: zavg_BEF on A 159 s (3.99 GB); 20 subsamples, 8 halves, 29 GLS runs (1–3 s each), Pearson on A | 4,218 s total (70.3 min) |
| Split-half selection (48 jobs, 3 processes) | 4,115 s wall |
| Map build zavg / GLS | 250.3 / 241.8 s |
| Tiers + bridge + pair-level AP | 52.0 s |
| Map variability: 20 subsample maps × 2 methods / analysis | 325.1 s / 32 s |
| FZ-X freeze (families, events, lockbox-testable list) | 27.6 s |

<!-- CON "Runtimes"; results/part2/construction/runs/timings.jsonl (calc: sum runtime_s = 4218.1 over 59 records); results/part2/construction/zavg_map/runtime.json: total_s; results/part2/construction/gls_map/runtime.json: total_s; BR: timings.total_s; results/part2/construction/variability/maps_runtime.json: maps_wall_s; CON "Map-to-map variability" (analysis 32 s); results/part2/fzx/FZX_manifest.json: runtime_s -->

### 4.3 Tiers

BioPlex 3.0 and CCLE proteomics were read only after the map commit, which a code guard enforces. <!-- CON intro; results/part2/verify/fzc_review.md §1 -->
- **T1**: within-cluster BioPlex edges above a degree- and essentiality-preserving null.
- **T2**: not T1, but protein co-abundance above a size-matched null (A lines only).
- **T3-neg**: testable, not supported. **T3-nc**: not testable.

Among non-locus complexes, 155 of 263 (59%) are T1, against 151 of 295 (51%) for GLS. <!-- MS: Results (SUMC: zavg.tiers_nonlocus_level1, gls.tiers_nonlocus_level1) --> Complex names are post hoc HGNC gene-group annotations. Some level-1 complexes merge two complexes or are pathway modules.

### 4.4 The bridge (primary methods endpoint)

The endpoint is size-standardised within-cluster BioPlex precision at matched gene coverage, with a protein-node bootstrap (B = 2,000, seed 5000). <!-- CON "Bridge" --> The zavg map dropped its 115 least stable complexes to match coverage. Precision was 0.157 for zavg and 0.182 for GLS: **−0.025 [−0.075, +0.031], null**. <!-- BR: primary.coverage_matching.trimmed.n_clusters_dropped, primary.zavg, primary.gls, primary.diff, primary.ci95 -->

| Sensitivity / secondary | zavg − GLS [95% CI] |
|---|---|
| No coverage matching | −0.020 [−0.054, +0.017] |
| Locus clusters excluded (pre-declared) | +0.003 [−0.060, +0.074] |
| Protein co-abundance (secondary) | −0.061 [−0.101, −0.018] (favours GLS) |
| Pair-level AP against BioPlex over S × S (no CI) | 0.0345 vs GLS 0.0215, Pearson 0.0223 |

<!-- CON "Bridge" sensitivities and secondary tables; BR: sensitivity.*, secondary.protein_coabundance, secondary.pair_level_AP_bioplex -->

An independent reviewer rebuilt both maps byte-identically and reproduced the bridge exactly (ACCEPT-WITH-NOTES). <!-- results/part2/verify/fzc_review.md "Verdict"; RT: Tree "FZ-C independent review" -->

### 4.5 Directions

Every direction is a machine-readable registry entry (`docs/registry/*.yaml`). The engine refuses unregistered or modified entries and appends every run to a trial ledger. Each direction ran once on A, in registry order, after a committed intent. <!-- MS: Methods "Registry, amendments and reproducibility"; MS: Results "Most complex-level context associations" --> Complex score: the mean clipped CRISPR effect over fixed members, with the context gene, its in-complex paralogues and nearby copy-number genes removed. Model: OLS with HC3 errors, Satterthwaite df, and PatientID-blocked Freedman–Lane permutation for screening hits (10,000 permutations; GPD tail only below 10⁻⁴). <!-- MS: Methods "Association statistics" -->

| Direction | Question | Tests (m) | Discoveries | Error control | Wall time |
|---|---|---|---|---|---|
| H3a lineage atlas | One-vs-rest lineage × complex | 8,441 | 1,145 in 21 of 23 families (942 non-locus) | Benjamini–Bogomolov selective FDR, q = 0.05 | 450.7 s |
| H3b driver/genome state | 20 contexts × complex | 7,340 | 161 in 13 families (121 non-locus) | BB selective FDR | 199.5 s |
| H1 within-complex GI | Does genetic loss of g change dependency on the rest of its complex? | 24 evaluable (registered 28) | 5 | BH | 53.4 s |
| H1S silencing (exploratory) | As H1, for expression-silencing events | 293 (registered 338) | 31 | BH | 437.4 s |
| H2 sub-module divergence (exploratory) | Does one sub-module carry the dependency beyond its siblings? | 452 in 210 families | 109 | Stage-wise | 57.1 s |

<!-- DR/H3a_lineage_atlas/summary.json: m_evaluable, n_discovery, n_discovery_nonlocus, n_families; MS: Results (trial_line.n_families_selected_perm); DR/H3b_driver_family/summary.json; DR/H1_within_complex_GI/runtime.json note; MS: Results H1 (DR/H1_within_complex_GI/summary.json); DR/H1S_silencing_GI/runtime.json note; RT: Tree "H1S"; DR/H2_submodule_divergence/summary.json (MS: m_evaluable, n_eligible_families, n_discovery); DR/*/runtime.json: wall_s -->

- **BB error control.** BB controls the within-family false-discovery proportion averaged over selected families. In 1,000 simulations its error rate was 0.046; the pooled proportion, 0.068, is reported but not claimed (D6). <!-- results/part2/verify/h3a_fdr_diagnosis.md "Validation"; RT: D6 -->
- **H1 classes.** Paralog compensation, residual sensitisation, collapse, or neutral (90% CI within ±0.1 Chronos units by TOST). A class also needs the effect outside a matched null of random member sets. <!-- MS: Results H1 -->

### 4.6 Controls and kill rule

There are 23 pre-registered positive controls, each mapped to a direction. They ran once on A, after the map and context freezes, in 28.0 s. <!-- docs/registry/README.md; CTRL/runtime.json: wall_s --> **Kill rule:** a direction is dropped if fewer than half of its evaluable complex-level controls are recovered.
- All 15 evaluable complex-level controls were recovered (one-sided p from 3.3 × 10⁻⁸⁸ to 0.023). <!-- MS: Results "Positive controls" (CTRL/kill_rule.json; CTRL/controls_summary.tsv) -->
- H3a passes 2 of 2, H3b 10 of 10, H1 3 of 3. <!-- CTRL/kill_rule.json: *.evaluable, *.recovered -->
- H2 has no evaluable control, so it is exploratory. <!-- CTRL/kill_rule.json: H2.verdict -->
- Coverage is thin: H3a rests on 2 controls, and only one H1 recovery (PC15) is a multi-gene readout. <!-- CTRL/kill_rule.json: H3a.thin_coverage; RT: Tree "Controls independent review" -->

### 4.7 Registered calibrations

These were implemented and run once each **after** the discovery runs and before any lockbox read (ruling R4). <!-- RT: D8 R4 --> Runtimes: H3a 35.7 s, H3b 86.3 s, H2 33.0 s. <!-- DR/*/calibration_D8.json: runtime_s -->

| Calibration | Result |
|---|---|
| Pseudo-contexts matched in size and lineage composition (fraction p < 0.05) | H3a 4.9%, H3b 4.5%; response threshold 0.075 not reached |
| Matched random gene sets (10 per complex; size, essentiality and coherence; fraction p < 0.05, random vs real) | H3a 28.1% vs 28.1%; H3b 11.6% vs 12.0% |
| Discoveries whose \|t\| beats all their matched sets | H3a 502 / 1,145; H3b 103 / 161 |
| Never-expressed gene sets, fraction p < 0.05 (expected 5%) | H3a 13.5% (melanoma, Skin, Lymphoid); H3b 6.2% |
| CN segment null | 6.0% (≤ 0.075; no response) |
| H2 random member splits | 78.3 rejections per replicate vs 109 real discoveries |

<!-- DR/H3a_lineage_atlas/calibration.json: frac_p_lt_005; DR/H3b_driver_family/README.md; DR/H3a_lineage_atlas/calibration_D8.json and DR/H3b_driver_family/calibration_D8.json: matched_random_sets.pooled.frac_p_lt_005, matched_random_sets.real_pooled.frac_p_lt_005, matched_random_sets.discoveries_exceeding_all_matched_sets, never_expressed_sets.pooled.frac_p_lt_005; RT: Tree "Registered calibrations" (H3b never-expressed 6.2%, CN segment 6.0%); DR/H2_submodule_divergence/calibration_D8.json: random_splits.stagewise_rejections_per_replicate.mean, random_splits.real.n_discovery -->

Pre-declared sensitivities (none is a trial): <!-- MS: Results "Sensitivity analyses" -->
- **Screen-quality covariates** kept 843 of 1,145 H3a and 115 of 161 H3b discoveries (366.7 s and 151.2 s). Esophagus/Stomach (2 of 30), WGD (12 of 26) and MSI-H (7 of 14) were fragile. <!-- DR/SENSITIVITY_SCREEN.md; DR/*/sensitivity_screen_quality/summary.json: runtime_s -->
- **Conservative p = max(p_perm, analytic)** kept 1,113 and 130. <!-- DR/gpd_sensitivity_summary.json: H3a_lineage_atlas.n_conservative, H3b_driver_family.n_conservative -->
- **Lymphoid lineage-scale effect.** In Lymphoid, a complex's β correlates with its essentiality (ρ −0.240 over 367 complexes; −0.228 with screen-quality covariates). <!-- DR/H3a_lineage_atlas/sensitivity_screen_quality/scale_effect_rho.tsv rows lineage:Lymphoid (MS) -->

### 4.8 Lockbox confirmation

The single read happened only after several gates:
- A dry check stopped the read because four registered steps had no code (D9).
- Amendments A3–A5 were made, and a pre-read code review returned GO-WITH-NOTES.
- A dress rehearsal on A with a dummy L loader took 2.7 h.
- An intent with the hashed hypothesis list was committed.

<!-- RT: D9, D10; LB/REHEARSAL.md "Result" (9,734 s) --> The 30 registered entries were then read once each, in registry order, with no re-read and no resume. The read took 2.6 h, with a peak of 3.67 GB. <!-- RT: D12; RT: Tree "Lockbox confirmations" --> **Confirmation rule:** same sign and one-sided p < 0.05 in L, for contexts with ≥ 5 context and ≥ 5 reference L lines. Rarer contexts take a labelled rare-context route (leave-one-lineage-out, DEMETER2 on A lines, leave-one-line-out), which is **not** an independent replication. <!-- MS: Methods "Lockbox confirmation" -->

| Direction (primary = read as executed, R7) | Confirmed / in lockbox | NOT EVALUABLE | Rate, evaluable | Rate, all-in |
|---|---|---|---|---|
| H3a | 601 / 1,145 | 0 | 0.525 | 0.525 |
| H3a, beats all matched sets / not | 287 / 502 vs 284 / 581 | 0 | 0.572 vs 0.489 | — |
| H3b | 60 / 150 | 28 (26 in WGD) | 0.492 | 0.400 (corrected relook 0.46) |
| H2 (exploratory) | 76 / 109 | 1 | 0.704 | 0.697 |
| H1S (exploratory) | 9 / 18 | 1 | 0.529 | 0.500 |
| Controls (23 entries; 44 readouts) | 29 / 44 | 13 | 0.936 | 0.659 |

<!-- RATES rows overall, beats_all_matched_sets, not_beats_all: primary_n_lockbox, primary_n_not_evaluable, primary_n_confirmed, primary_rate_evaluable, primary_rate_all_in; MS: Results "Of the 28 failed H3b fits, 26 are in WGD" (CMP; LD: H3b.confirmed_per_context.WGD); RT: D14 (corrected 0.46) -->

- **Power aid.** A post hoc aid by the independent reviewer (not registered, not a test) predicts 0.47–0.52 for H3a and 0.43–0.44 for H3b. **Confirmation is at the rate power predicts.** The beats-matched-sets gap is about what effect size and power predict (0.52–0.53 vs 0.41–0.43). <!-- MS: Results "About half the associations confirm" (LD: power_band) -->
- **H1 in L.** ARID1A → ARID1B was confirmed (it is control PC6). EP300 → CREBBP passed the rare route (4 L loss lines). FOCAD is defect-sensitive. TP53 was confirmed but cannot be promoted (E2). PTEN was not confirmed. <!-- MS: Results "H1 in the lockbox"; RT: D12 -->

### 4.9 BENCH and BRIDGE-L

- **BENCH (complex vs its best gene g\*)**. Over 1,295 discovery pairs, Δ = t_L(complex) − t_L(g\*).
  - High-breadth tercile (registered criterion): +0.240 [−0.033, +0.362]. The CI includes 0, so **complex scores did not out-replicate their best gene.**
  - Low-breadth tercile: the best gene replicated clearly better (−0.338 [−0.515, −0.201]).
  - Overall mean: −0.066 [−0.155, +0.031].
  - In the bootstrap, 16.5% of pair × replicate fits failed silently. With the covariate design assembled per replicate, the high-tercile lower bound is −0.122. The conclusion holds either way.
  - Runtime 8,192.7 s (B = 2,000).

  <!-- LBR: BENCH_complex_vs_gene.distinct.n_pairs; BENCH: primary.mean_high, primary.ci95_high, primary.mean_low, primary.ci95_low, primary.mean, primary.ci95; MS (LD: BENCH.frac_pair_rep_nan, BENCH.ci95_high_assembled); LBR: BENCH_complex_vs_gene.runtime_s -->
- **BRIDGE-L** asks whether zavg-map context discoveries confirm in L more often than discoveries from identical machinery on the GLS map. The answer is no.
  - R7 primary, evaluable hits (zavg 1,267, GLS 1,326): 0.514 vs 0.578, **−0.065 [−0.107, −0.025]**.
  - As executed: −0.072 [−0.111, −0.033]. Corrected relook: −0.069 [−0.108, −0.029].
  - Matched gene coverage: −0.074 [−0.122, −0.026]. T1/T2 clusters only: −0.074 [−0.120, −0.026].
  - Runtime 733.8 s.

  <!-- BLR7: primary_evaluable.n_hits, primary_evaluable.primary.rate_zavg, .rate_gls, .diff, .ci95; as_executed.primary.diff, .ci95; corrected_nonconfirmatory.primary.diff, .ci95; primary_evaluable.multiverse.matched_gene_coverage, primary_evaluable.multiverse.T1_T2_clusters; LBR: BRIDGE-L_confirmation_rate.runtime_s -->

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

<!-- results/part2/contexts/build_log.json: runtime_s; MAN: runtime_s; results/part2/fzx/FZX_manifest.json: runtime_s; calc: 4218.1 + 4115 + 250.3 + 241.8 + 52.0 = 8877.2 (CON "Runtimes"); CTRL/runtime.json: wall_s; calc: 450.67+199.45+53.35+437.42+57.11 = 1198.0 (DR/*/runtime.json); calc: 35.7+86.3+33.0 = 155.0 (DR/*/calibration_D8.json); calc: 366.7+151.2 = 517.9 (DR/*/sensitivity_screen_quality/summary.json); RT: D12; ROB: runtime_s -->

The cost is dominated by the zavg_BEF score runs (the ARPACK embedding takes 118.6 of the 159 s per run on A) and the BENCH bootstrap. <!-- results/part2/construction/runs/timings.jsonl row set=A, method=zavg_BEF: stages_s.embedding_arpack, runtime_s; LBR -->

---

## 5. Best methods and strongest findings

### 5.1 Best method: the locked score at the pair level

zavg_BEF is the best pair-level co-complex ranker tested:
- It beat GLS, Spearman and Pearson on held-out hu.MAP 3.0 complexes under a rule frozen before the look.
- It also beat them on protein-new pairs.
- It ranked BioPlex pairs better than GLS and Pearson (AP 0.0345 vs 0.0215 / 0.0223; no CI).

<!-- TR: summary; results.stratum_C3; BR: secondary.pair_level_AP_bioplex -->

Its scope is limited in two ways. The gain is carried by pair-rich complexes (Section 3.3). The method was also selected on hu.MAP training labels derived from Complex Portal, so agreement with curated complexes partly re-measures that selection. <!-- MS: Discussion; DESIGN.md §0 F3 --> For complex construction, the zavg map's advantage is coverage (3,642 vs 2,862 genes), not better complexes (Section 6). <!-- SUMC: zavg.genes_level1, gls.genes_level1 -->

### 5.2 The registered strongest findings (254)

**Rule.** The rule was fixed in DESIGN §2.5 and the registry. A finding qualifies if all of these hold:
- discovery q < 0.05;
- primary lockbox confirmation (or a passed rare route);
- not defect-sensitive and not a locus cluster;
- leave-one-out robust (H3a by primary disease; H3b and H1 by lineage).

<!-- ROBR: "Rule (applied mechanically)" --> The per-hit checks were registered but computed only after the lockbox read, on A only, against a candidate list committed first (fd03b03). <!-- ROBR: "Disclosure" -->

| Direction | Candidates | Strongest | Units | T1 / T2 / T3-neg |
|---|---|---|---|---|
| H3a | 524 | 207 | 137 | 124 / 56 / 27 |
| H3b | 46 | 45 | 41 | 27 / 9 / 9 |
| H1 | 2 | 2 (ARID1A → ARID1B, control PC6; EP300 → CREBBP, rare route) | 2 | 1 / – / 1 |
| **Total** | | **254** | | **152 / 65 / 37** |

<!-- ROB: per_direction.*.candidates, .strongest, .strongest_units, .strongest_per_tier; ANN: n, tier -->

- **The leave-one-out rule is strict for single-disease lineages.** It removed 317 of 524 H3a candidates, 308 of them through exactly one fold. Bowel kept 0 of 29.
- **It is vacuous where the remainder has fewer than 2 lines.** Pancreas passed 17 of 17. Leave-one-out is partly untestable for 21 H3a and 33 H3b findings.

<!-- ROB: per_direction.H3a.loso_fail_any_fold; ROBR: "H3a LOSO behaviour"; ANN: loso_partially_testable -->

**Post hoc descriptive annotations** (none selects or ranks): <!-- ANN -->
- **Matched sets.** 107 of 207 H3a and 34 of 45 H3b findings beat all their matched sets.
- **Control readouts.** 7 findings are readouts of registered controls; 18 share a unit with a control read in another context.
- **Lineage-scale trend.** 142 of 207 H3a findings have the sign the lineage-scale trend predicts. 79 depart from it by |residual z| ≥ 2.

**Reporting labels** (post hoc, D16/D17): <!-- SD: lead_example, possible_lineage_scale_artefact, possible_lineage_scale_artefact_per_context, median_unit_essentiality_artefact, median_unit_essentiality_all_units -->
- **Lead set** (37): H3a, T1, not a control readout, |residual z| ≥ 2.
- **Possible lineage-scale artefact** (56: 52 Lymphoid, 4 Myeloid): scale-concordant with |residual z| < 2, on units more essential than average (median unit essentiality −0.63 vs −0.30).

Examples from the lead set (A β, then L β [95% CI], L lines; negative = more dependent): <!-- SF rows as cited per line (via MS Results "The registered strongest findings") -->
- Ovary/Fallopian Tube, unit with CDH1, CDH2, CTNNA1, CTNND1, PARD3, PARD6B, PRKCI, HNF1B, PAX8: −0.140; −0.202 [−0.298, −0.107], 16 lines; residual z −4.47. <!-- SF row Z0035 -->
- Bone, SNAPC1–SNAPC5: −0.210; −0.251 [−0.370, −0.132], 15 lines; residual z −3.74. <!-- SF row Z0303 -->
- Lung, APC/AXIN1/CSNK1A1/CTNNBIP1 unit: +0.086; +0.077 [+0.032, +0.122], 31 lines; residual z +3.01. <!-- SF row Z0174 -->
- Head and Neck, EGFR/GRB2/GAB1/PTPN11/SHC1 unit: −0.133; −0.107 [−0.155, −0.059], 19 lines; residual z −2.58. <!-- SF row Z0051 -->
- Myeloid, AP2 adaptor unit: −0.189; −0.239 [−0.380, −0.098], 11 lines; residual z −3.09. <!-- SF row Z0254 -->
- Pleura, heparan-sulfate-synthesis-like unit (EXT1/EXT2 …): −0.265; −0.182 [−0.408, +0.044], 5 lines. This one rests on a small lineage, and its leave-one-out check was only partly testable. <!-- SF row Z0067 -->
- Outside H3a: RB1-LoF lines depend more on a cyclin E–CDK2/E2F unit (A −0.201; L −0.156 [−0.235, −0.076], 17 lines). <!-- SF row Z0069 -->

All functional readings above are post hoc, from gene symbols. These findings are replicated, robustness-checked module-level associations. Neither the matched-set comparison nor the lockbox shows that the complex, rather than one or two members, is the operative unit, and none has been tested outside DepMap. <!-- MS: Discussion "What survives" -->

---

## 6. What did not work

- **The BioPlex bridge is null.** The pre-registered thesis was that the zavg map agrees better with physical evidence than the GLS-module standard. It was not supported: −0.025 [−0.075, +0.031]. Protein co-abundance favoured GLS (−0.061 [−0.101, −0.018]). <!-- BR: primary.diff, primary.ci95, secondary.protein_coabundance --> Locus clusters explain the deficit only post hoc (+0.003 with them excluded). <!-- BR: sensitivity.locus_excluded.diff -->
- **BRIDGE-L favours GLS.** zavg-map discoveries confirmed less often (−0.065 [−0.107, −0.025]). This held at matched coverage and among T1/T2 clusters, so the locus explanation of the first bridge does not explain this one. We have no tested explanation. <!-- BLR7: primary_evaluable.primary, primary_evaluable.multiverse -->
- **Complex-level associations are mostly not complex-specific.** Matched random gene sets reach p < 0.05 as often as real complexes (H3a 28.1% vs 28.1%). <!-- DR/H3a_lineage_atlas/calibration_D8.json --> Never-expressed sets show a lineage excess (13.5% vs 5% expected), which points to a lineage shift of the CRISPR effect scale. <!-- DR/H3a_lineage_atlas/calibration_D8.json: never_expressed_sets.pooled.frac_p_lt_005 --> Most H2 divergence is also seen in random splits (78.3 vs 109). <!-- DR/H2_submodule_divergence/calibration_D8.json --> The lockbox adds no independent evidence that "beats matched sets" marks complex biology, because the gap is what power predicts. <!-- MS: Results; RV §4 -->
- **The complex does not beat its best gene.** In the registered high-breadth tercile the CI includes 0. Where few members carry the association, the single gene replicates clearly better. <!-- BENCH: primary.ci95_high, primary.ci95_low -->
- **H1 recovers known biology.** All three headline-eligible H1 events are registered controls or published:
  - ARID1A → ARID1B is control PC6;
  - FOCAD → PELO–HBS1L overlaps control PC10;
  - EP300 → CREBBP is the reciprocal of the published CREBBP-deficient → EP300 dependency.

  The informative H1 result is that 11 of 24 genetic losses left their complex's dependency unchanged within ±0.1 Chronos units. L cannot re-establish that neutrality: only 2 of the 11 L 90% CIs lie inside the margin. <!-- RT: D11; MS: Results H1 (PD: H1.n_neutral, H1.n_events; LD: H1_lockbox.neutral_L_ci90_inside_margin) -->
- **Not run** (declared before the read): BENCH within-line-rank and the A–L random-set correlation; alternative per-hit scores (eigengene, within-line rank, expression-masked, pan-essential-referenced); per-hit robustness for H1S and H2. <!-- AMD A3–A5; SUP Note 14 "(not run)" -->

---

## 7. Process and integrity

### 7.1 Amendments, errata and rulings, with timing

Outcome access points (UTC):
- test labels: 28 Oct 18:18;
- first A-side association outcome: 29 Oct 00:20 (H3a);
- H1 outcome: 00:35;
- calibration outcomes: from 01:03;
- lockbox read: 06:50:05 to 09:23:05.

<!-- SUP Note 14 header -->

| ID | What changed | When (UTC) | Relative to outcome access |
|---|---|---|---|
| A1 | 318 silencing events moved from H1 to exploratory H1S | 28 Oct 22:54:50 | Before any H1 run; labels only |
| A2 | Genetic-loss rule (≥ 50% LoF or hom-del); SMARCB1 moved into H1 | 28 Oct 23:28:31 | Before any H1 run; labels only |
| R1 / E1 | Sign typo in the registered H1 aligned score corrected to the stated meaning (the engine already used it) | 29 Oct 00:41:42 | After the H1 run |
| R2 | Engine crash on leverage-1 H1S designs fixed; H1 not re-run | 29 Oct 00:47 | After H1, before H1S |
| R3 | GPD-extrapolated p stays primary; p_cons sensitivity pre-declared | 29 Oct 00:59 | After discovery, before the sensitivity |
| R4 | Registered calibrations implemented and run once | 29 Oct 01:02–01:06 | After all discovery runs; before L |
| R5 / E2 | "Detectable, within equivalence margin" label (TP53, PTEN); not promotable | 29 Oct 01:07 | After the H1 run and its review; post hoc, conservative |
| A3 | Lockbox permutation seeds; BENCH within-line rank not run | 29 Oct 02:03:27 | Before L |
| A4 | BENCH secondaries not run; multiverse point estimates only | 29 Oct 02:42:37 | Before L |
| A5 | A4 erratum; BRIDGE-L multiverse CIs; literal readings | 29 Oct 03:25:03 | Before L |
| R6 | Failed L fits relabelled NOT EVALUABLE; corrected relook authorised (non-confirmatory) | 29 Oct 10:19:44 | After the read |
| R7 | Read as executed stays primary; 560-fit corrected relook alongside; defect-sensitive calls carry no claim | 29 Oct 10:41:49 | After the read and the failed R6 check |
| D14–D16 | Reporting rules; robustness checks computed; 254 fixed; post hoc annotations and labels | 29 Oct, after the read | After the read; A only |

<!-- SUP Note 14 table; SUP Note 7; AMD; LB/README.md (R6, R7 verbatim) -->

The full record, including D5–D16, is in `manuscript/supplementary.md` Note 14 and `docs/registry/AMENDMENTS.md`. Rules fixed after outcomes were seen are listed as such: E2, the defect-sensitive label, the power aid, the robustness readings, and the lead and artefact labels. Each makes the reading more conservative or is descriptive. <!-- MS: Discussion "Limitations" -->

### 7.2 The lockbox driver defect

- **The defect.** The driver built the L covariate design over all models and then subset it to L. When a context's reference lineage had no L line, single-line lineages formed a spurious group. <!-- MS: Results "The lockbox read, and a defect" -->
- **Footprint.**
  - 53 fits failed: 50 because of the defect and 3 genuinely (CTRL_PC21 ×2, H1S HNF1B).
  - 510 successful fits carried a spurious column (median |Δβ| 2.4 × 10⁻⁴, max 0.046).
  - H3a is unaffected. Of the 28 failed H3b fits, 26 are in WGD. <!-- LB/README.md "PI ruling R7"; MS -->
- **Handling.**
  - R7 keeps the single read as executed as the primary result.
  - A corrected, non-confirmatory relook of all 560 affected fits is reported alongside. It reads persisted matrices through the resume path only.
  - Any call that differs between the two is defect-sensitive and supports no claim: 53 calls (H3b 30, BRIDGE-L 10, controls 11, H1 FOCAD 1, H2 1).
  - Cost: WGD and four controls (PC9, PC10, PC12, PC13) lose their lockbox evidence.

  <!-- LB/README.md R7; MS (LD: defect_sensitive_total, defect_sensitive_by_entry); MS: Results "H1 in the lockbox" -->
- **How it was caught.** The dress rehearsal had shown the same failures, but nobody checked fit status. The failed fits were found after the read (D12). The R6 reproduction check (68f0d3f) then showed the 510 shifted fits, so the scope grew from 53 to 560. The results reviewer's label-only span audit predicted the 510/50/3 split exactly. <!-- RT: D12, D13, D14; RETRO "Stop" 2–3 -->

### 7.3 Registry ↔ code gaps

Four registered steps had no code when the step they gate ran. <!-- RETRO "Stop" 1; LESSONS.md rule 15 -->
1. **Calibrations.** Missing at the discovery runs; implemented after them (R4). Found when review 2 and the PI compared the registry text with the engine. <!-- RT: D8 R4; results/part2/verify/discovery_review.md "Review 2" -->
2. **Lockbox machinery.** The driver, the BENCH secondary, BRIDGE-L GLS discovery and lockbox seeds were missing. The driver's own dry check stopped before any L read (D9); the gaps were implemented as registered (A3–A5). <!-- RT: D9; LB/BLOCKED.md -->
3. **L covariate design.** Built differently from the A side (Section 7.2). The pre-read code review compared the driver with the registry in tables but did not execute the L path against the A path. <!-- RETRO "Stop" 2 -->
4. **Per-hit robustness** (LOSO, leave-top-2-out). Registered in DESIGN, never computed at discovery. Found at manuscript time (D15) and computed after the read on A only. <!-- RT: D15; results/part2/robustness/INTENT.md -->

There was also one earlier code↔registry mismatch, caught on synthetic data before any real association: the GPD tail rule plus 1,000 permutations inflated the H3a test size. A multi-simulation FDR test failed, and the diagnosis fixed it (BB rate 0.046 over 1,000 simulations). <!-- results/part2/verify/h3a_fdr_diagnosis.md "Trigger", "Verdict"; RT: Tree "Association engine" -->

### 7.4 Independent reviews (agents that did not produce the work)

| Review | Verdict | What it caught or confirmed |
|---|---|---|
| Design rounds 1–3 (6 + 5 + 1 reviewers) | DESIGN v3.1; protocol ready to freeze | Bridge CI by protein bootstrap; genotype lockbox strata; lockbox read order |
| Part 1 audit | PASS | No test leakage into Part 1; its intermediate run caused, and disclosed, D1 |
| Test verifier | REPRODUCED (≤ 2.6 × 10⁻¹⁵) | Numbers committed before comparison |
| FZ-C review | ACCEPT-WITH-NOTES | Byte-identical rebuild; bridge exact; disclosed the pre-freeze proteomics read for the LoF check (D5) |
| FZ-X review (+ delta) | ACCEPT-WITH-NOTES → ACCEPT | SMARCB1 was typed silencing; it became amendment A2 |
| Controls review | ACCEPT-WITH-NOTES | 15 tests reproduced; undocumented lineage-merge convention (≤ 2% change in p) |
| Discovery reviews 1 and 2 | ACCEPT-WITH-NOTES | 15,781 tests reproduced; GPD-tail dependence; fragile PTEN call; Lymphoid scale effect; specificity per hit |
| Pre-read lockbox code review | GO-WITH-NOTES | Persist/hash/resume guard; amendment record. Missed the covariate defect |
| Lockbox results review | PASS | 54/54 calls recomputed under both designs; 510/50/3 predicted from labels; 16.5% silent BENCH failures; power aid |
| Reference check | E1–E15 applied | 9 metadata errors, 2 claim mismatches |

<!-- RT: "Previous status (M2)", D1, D5, D14, D17, Tree entries; docs/reviews/round3_final.md (DESIGN.md header); TV; results/part2/verify/{fzc_review,fzx_review,controls_review,discovery_review,lockbox_code_review,lockbox_results_review}.md "Verdict"; RETRO "Stop" 5 (a303165) -->

Reviewers and producers share a model family. Reviews reproduced every confirmatory number, but they are not independent in the human sense. <!-- MS: Discussion "Limitations" -->

---

## 8. Limitations

1. **Procedural single look.** The test labels were plain files, and label-free scores already covered all test pairs. D1 shows the risk was real. <!-- SUP Note 15 item 1 -->
2. **Narrow test win.** The gain over GLS is carried by pair-rich clusters, with an effective cluster count of 21.3. The training labels derive from Complex Portal, so the test partly re-measures selection. <!-- TR: results.itt.cluster_stats.kish_neff_clusters; SUP Note 15 items 2–3 -->
3. **Calibration.** CITL was 0.893, with no pre-specified CI. Threshold and calibrator sampling error is not propagated. <!-- TR: results.itt.metrics.zavg_BEF.calibration.isotonic.calibration_in_the_large; MS: Methods -->
4. **Map choices.** Selection is flat across z0. 28% of complexes are locus clusters. Some complexes merge two complexes or pathways. About 40% of S genes are in no complex. <!-- SEL; SUMC: zavg.locus_rate_level1; SUP Note 15 items 6–10 -->
5. **Calibrations and robustness checks ran late** (R4, D15). The calibration response scope's literal reading was confirmed after results were seen. <!-- RT: D8, D15 -->
6. **Lockbox.** L has 338 lines and has now been read. Non-confirmation of small effects is often uninformative (325 of 1,145 H3a discoveries have |β| < 0.05). The rare route is not independent replication. The driver defect leaves 53 claimless calls. **No held-out cell lines remain for a new confirmatory test.** <!-- MAN: n_L; MS (PD: H3a.n_abs_beta_lt_0.05); LD: defect_sensitive_total -->
7. **Lineage-scale artefacts.** The Lymphoid ρ −0.24 and the never-expressed excess mean some lineage findings may reflect the Chronos dynamic range. The residual-z labels use an unregistered |z| ≥ 2 threshold. <!-- ANN: scale_effect_H3a.per_context.lineage:Lymphoid.rho; SUP Note 15 item 16 -->
8. **One resource.** All dependencies come from DepMap 26Q3 CRISPR (Chronos). Nothing was tested in another screen collection or in vivo. <!-- SUP Note 15 item 27 -->
9. **Thin controls.** H3a rests on 2 controls, and H2 has none. <!-- CTRL/kill_rule.json -->
10. **Unverified licences.** Licences are unverified for the hu.MAP test files, CCLE proteomics and BioPlex. The HGNC download date is unknown. <!-- DATA.md -->
11. **Agent-run study** under a human PI. The run log of the lockbox read was not committed, so the loader count (30) is inferred from ledger invariants. <!-- RV §1 -->

---

## 9. Follow-ups and remaining ideas (prioritised)

L has been read. Any new confirmatory claim therefore needs its own registration and new held-out data; work on A or L from here is exploratory. Ordered by expected gain in confidence, then novelty:

1. **Lineage-scale normalisation (highest priority for confidence).** Re-express complex scores relative to the lineage's essentiality trend before testing. Options are residualising on unit essentiality per lineage, using a rank or quantile-normalised effect, or testing against matched random sets as the null rather than a flag. This targets the two largest threats to the context findings: the Lymphoid scale effect (ρ −0.24) and the never-expressed excess (13.5%). The matched-set test needs far more than 10 sets per complex to give a usable empirical p-value. <!-- ANN; DR/H3a_lineage_atlas/calibration_D8.json; SUP Note 15 item 11 -->
2. **Orthogonal replication with data already provided.**
   - **DEMETER2 v6 (RNAi).** It is in hand but used only in the rare-context route. A registered RNAi replication of the 254 strongest findings and the 11 neutral H1 events, in lines with RNAi data, would be the cheapest independent perturbation test.
   - **CCLE proteomics.** Used only for T2 tiers and the LoF check. Testing whether contexts change complex protein abundance would add a non-CRISPR readout.

   <!-- DATA.md; MS: Methods "Lockbox confirmation", "Tiers"; SD: n -->
3. **External CRISPR replication.** A second screen collection would be needed. Per PROMPT.md it must be requested from the user, not downloaded.
4. **Explain the bridge failures.**
   - Locus-aware construction (copy-number-corrected profiles, or excluding arm-block edges before clustering), then the bridge comparison again.
   - A BioPlex-integrated map (D4 remaining idea). It would need a held-out BioPlex split or proteomics as the validator, because BioPlex could no longer validate itself.
   - The untested hypothesis to examine: any reasonable clustering already captures the pair-level gain on pair-rich complexes, while the zavg map's extra coverage adds weaker modules. <!-- RT: D4; MS: Discussion -->
5. **H4 copy-number dosage direction** (CYCLOPS-like partial-loss vulnerabilities at complex level; cut at D4). Novel, but locus artefacts make it hard. <!-- RT: D4; DESIGN.md §2.6 "Cut" -->
6. **Rewiring.** Context-specific change of complex coherence (differential co-essentiality), and BioPlex 293T vs HCT116 interactome differences against lineage co-essentiality. Both were cut at D4 and remain open. <!-- RT: D4 -->
7. **Close the declared not-run items.** Alternative per-hit scores (eigengene, within-line rank, expression-masked); BENCH within-line rank and the A–L random-set correlation; per-hit robustness for H1S and H2. <!-- SUP Note 14 "(not run)" -->
8. **Pair-level method.** The test points to small complexes (≤ 5 positive pairs) as where zavg_BEF does not beat GLS. Any new pair-level method needs a new held-out benchmark, because the hu.MAP test set has now been used. <!-- EXP -->
9. **Process** (from the M6 retrospective):
   - a mechanical registry↔code traceability table that fails on unmapped steps;
   - rehearsals that assert fit-status counts and A-path = L-path on identical rows;
   - committed run logs for irreversible runs;
   - word budgets and the number check at draft time.

   <!-- RETRO "Proposed pruned LESSONS.md" 11, 13–15 -->

---

## 10. Reproducibility

**Environment.** Run `scripts/00_setup_part1.sh` to clone the locked Part 1 repository (c7f95ee) into `external/part1_locked` and build `.venv`. Set `PART1_DIR` and `PYTHONHASHSEED=0` for every step. `scripts/check_data.sh` verifies every dataset checksum against DATA.md. <!-- scripts/00_setup_part1.sh; scripts/check_data.sh; SUP Note 9 -->

**Pipeline, in order** (each step's exact command is in its intent note):

| Stage | Scripts |
|---|---|
| Part 1 test evaluation | `19_make_pair_files.py` → `20_score_test_pairs.sh` (train, test; label-free ID files) → `21_evaluate_test.py --mode fit` → `22_freeze_protocol.py` (PI; freeze commit) → `21_evaluate_test.py --mode test --intent`, then `--run` → `23_export_test_curves.py` (descriptive). Verifier: `results/test_eval/verify/verify_independent.py` |
| Part 2 setup | `05_profile_part2.py` → `32_build_contexts.py` → `33_make_lockbox.py` |
| Construction | `30_toy_construction.py`, `35_toy_decisions.py` → `34_run_A_blocks.sh` (→ `31_score_blocks.py`) → `36_build_maps.py` → `37_validate_maps.py` (after the map commit) → `38_describe_maps.py` → `39_map_variability.py` |
| Freeze and amendments | `41_freeze_fzx.py` → `42_amend_A1.py` → `43_amend_A2.py` |
| A-side analyses | `34_run_controls.py` → `44_run_direction.py --id` H3a, H3b, H1, H1S, H2 → `45_gpd_sensitivity.py` → `46_run_calibrations.py` → `47_h1_label_final.py` → `48_sensitivity_screen.py` → `49_specificity_counts.py` |
| Lockbox | `52_amend_A3.py` → `53_amend_A4.py` → `54_amend_A5.py` → `51_lockbox_confirm.py` (single read; ledger blocks any re-read; `--relook R6` resume-only) → `55_lockbox_r7_report.py`; `55_export_member_profiles.py` (A only) |
| Strongest findings | `56_strongest_findings.py --candidates`, then `--run` → `57_annotate_strongest.py` |
| Figures, tables, checks | `tab_manuscript.py`, `tab_part2_discoveries.py`, `tab_lockbox.py`, `fig_1`–`fig_5_*.py`, `tab_strongest.py`, `check_manuscript_numbers.py manuscript/main.md` (read committed results only; byte-identical on re-run) |

<!-- scripts/ listing and docstrings; CON "Pipeline and files"; ROBR "Reproduce"; SUP Note 9; LB/CORRECTED_REFIT_INTENT.md "Commands" -->

**Seeds** (DESIGN §2.1, the registry and amendments): <!-- DESIGN.md §2.1; MAN: seed; CON; BR; MS: Methods; AMD A3; TV -->
- test bootstrap seed 0, MC seeds 1–4, spot-check seed 7;
- lockbox 3000;
- Leiden 0–4; subsamples 1000–1019; split halves 2000–2003;
- permutations 4000+ (controls 4400 + index, A3); BRIDGE-L GLS discovery 4500;
- bridge bootstrap 5000, co-abundance 5001, BENCH 5100, BRIDGE-L 5200;
- variability random reference 7000.

**Data.** DATA.md lists source, version, SHA-256 and licence for every input. The DepMap and HGNC checksums equal Part 1's. <!-- DATA.md -->

**Never committed** (`.gitignore`, retention rules):
- `datasets/` (all raw data, including the hu.MAP test files);
- `data_cache/` (score blocks, context tables, and the persisted lockbox reads, whose hashes alone are committed in `results/part2/lockbox_reads.jsonl` and which only the driver's resume path may open);
- `external/` (the Part 1 clone) and `.venv/`;
- binary arrays (`*.npy`, `*.npz`, `*.parquet`).

Re-running the lockbox read is blocked by the ledger by design. <!-- .gitignore; LB/README.md "Retention rule"; SUP Note 9 -->
