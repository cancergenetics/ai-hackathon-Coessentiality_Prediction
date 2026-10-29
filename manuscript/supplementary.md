# Supplementary Information

**A label-free co-essentiality score wins a pre-registered held-out test, but its advantage does not carry over to complexes or to complex-level cancer dependencies**

<!-- Source convention as in main.md: every number carries an HTML comment naming its committed source. "LB" =
results/part2/lockbox_confirm; "CMP" = LB/comparison_primary_vs_corrected_R7.tsv; "LD" =
manuscript/tables/lockbox_derived.json; "RV" = results/part2/verify/lockbox_results_review.md; "RVD" = its output
folder; "DR" = results/part2/directions; "PD" = manuscript/tables/part2_derived.json. -->

## Contents

- Supplementary Note 1. Transition from development to test: all prior contact with the test files
- Supplementary Note 2. Methods detail: construction, association engine, lockbox driver
- Supplementary Note 3. The lockbox read, the driver defect and the corrected relook
- Supplementary Note 4. Post hoc power aid for confirmation rates
- Supplementary Note 5. BENCH and BRIDGE-L estimands
- Supplementary Note 6. Registry summary (30 entries) with the lockbox outcome per entry
- Supplementary Note 7. Amendments, errata and rulings
- Supplementary Note 8. Member-resolved H1 profiles of the headline-eligible events
- Supplementary Note 9. Reproduction
- Supplementary Tables S1–S8 (files)

---

## Supplementary Note 1. Transition from development to test

The Part 1 method was locked at commit c7f95ee (28 October 2026, 12:55:43 UTC). The hu.MAP 3.0 test files appeared in the project folder at 13:55 UTC, 59.6 min later. <!-- docs/part1_audit.md "Transition timeline"; DATA.md --> The only contacts with the test files before the protocol freeze were the following. <!-- docs/TEST_PROTOCOL.md §0 -->

| ID | Who | What | Consequence |
|---|---|---|---|
| T0 / D0 | PI | Line counts (3,312 / 597,426) and the first 3 lines of each file, as a format check | None |
| T1 | Part 1 audit agent | Label-free descriptive counts: overlaps, strata, mappability | Known when the protocol was written; the choices follow Part 1's design |
| D1 | Part 1 audit agent, intermediate run | Counted test pairs in Part 1's label-free top-pair table (zavg_BEF ≥ 3), which is equivalent to precision and recall at z ≥ 3 | Counts deleted; PI never saw them; rule (a) reported only as "non-blind (D1)" and enters no verdict |
| T2 | Protocol author | Octal dump of the first ~200 bytes of the positive file | None |
| T3 | Everyone | Genome-wide label-free scores covering all test pairs existed | The single-look guarantee is procedural |

<!-- docs/TEST_PROTOCOL.md §0 rows T0–T3, D1 -->

Decision notes D2–D4 concern data scope and are not contacts. D2: no CORUM or Complex Portal was obtained, so any naming of complexes is post hoc. D3: the locked verdict is computed on S1, and ITT is quoted alongside. D4: the scope of Part 2 was cut, and BioPlex is used for validation only. <!-- RESEARCH_TREE.md D2–D4 --> An independent audit found no evidence that test information entered Part 1 (docs/part1_audit.md).

## Supplementary Note 2. Methods detail

*Construction.* The complete parameters of both maps are in `results/part2/construction/{zavg_map,gls_map}/parameters.json`, and the rationale for each construction rule is in `results/part2/construction/DECISIONS.md` (D1–D6):
- containment support for small complexes (D1);
- placement and trimming of overlapping candidates (D2);
- attachment of genes outside S (D3);
- edge-count matching of the GLS comparator (D4);
- tiers and the bridge endpoint (D5);
- split-half selection (D6).

A toy simulation fixed the rules before any real map was built: 86 of 105 planted complexes and all 39 planted sub-modules were recovered, and null data gave no cluster. <!-- results/part2/construction/DECISIONS.md "Toy results" -->

*Association engine.* Complex scores, the OLS model with HC3 errors and Satterthwaite df, the PatientID-blocked Freedman–Lane permutation, Benjamini–Bogomolov selective FDR and stage-wise H2 testing are implemented in `src/part2/` and specified per entry in `docs/registry/*.yaml`. The engine's lineage-merge convention merges indicator levels with fewer than two lines into the reference. DESIGN does not state it; it changes p by at most 2% (relative) in the controls. <!-- results/part2/verify/controls_review.md §2 --> In 1,000 simulations with planted lineage effects, the selective error rate was 0.046 and the pooled false-discovery proportion 0.068. <!-- results/part2/verify/h3a_fdr_diagnosis.md "Validation" -->

*Registered calibrations (D8, R4).* The matched random sets, never-expressed sets, copy-number segment null and H2 random splits were specified operationally in `results/part2/directions/CALIBRATION_D8_INTENT.md` before they were run. They were run once each, after the discovery runs and before any lockbox read. <!-- DR/CALIBRATION_D8_INTENT.md; git 1e674bd → fd582a6, 4fb83d7, 1a2a445 -->

*Lockbox driver.* `scripts/51_lockbox_confirm.py` reads each registered entry once through `assoc.read_lockbox`, which commits the ledger line before it loads and refuses a second read. It persists the read with its hash, then applies the registered confirmation rule. Hits in contexts with fewer than 5 context or 5 reference L lines take the rare-context route, which uses only A lines and DEMETER2. <!-- LB/README.md; LB/INTENT.md --> Lockbox permutation seeds are registered by A3, and the BENCH and BRIDGE-L bootstrap seeds by the entries (5100, 5200). <!-- docs/registry/AMENDMENTS.md A3 -->

## Supplementary Note 3. The lockbox read, the driver defect and the corrected relook

*Timeline (29 October 2026, UTC).*

| Time | Event | Commit |
|---|---|---|
| 01:45:43 | Dry check stops before any read: four registered steps have no code | 21fe1be |
| 02:03–03:25 | Amendments A3, A4, A5 (no L value read) | b349e98, d9a42dd, 5082cc8 |
| 03:17 | Pre-read code review: GO-WITH-NOTES | 2bd58c0 |
| 06:16 | Dress rehearsal on A with a dummy L loader (2.7 h); its failed fits went unchecked | fb54005 |
| 06:47:06 | Intent note with the hashed list of all hypotheses | ed065f4 |
| 06:47:56 | Decision D11 (framing of H1; novelty criteria) | cb3ce32 |
| 06:50:05 – 09:23:05 | Single read of 30 entries in registry order | ledger |
| 06:53:55 | Member-resolved H1 profiles exported from A (no L access) | a37e938 |
| 10:10:31 | Results committed | 5bc9444 |
| 10:11 | D12: defect found; uncorrected results recorded | 4af8e60 |
| 10:19:44 | R6 note, corrected-refit intent and fix | 4d48039, 38e301a |
| 10:26:22 | R6 reproduction check: 510 of 3,014 successful fits differ; refit stopped | 68f0d3f |
| 10:41:49 | R7 note and intent addendum | b0ae055 |
| 10:49:06 | Corrected relook output (560 fits) | 94b42c3 |
| 11:37:54 | Independent results review: PASS | 63e0e62 |

<!-- git log; RESEARCH_TREE.md D9–D14; RV §1 -->

*Diagnosis.* The driver assembled the one-hot covariate design over all models and then subset it to L. The engine instead assembles the design over the rows analysed and merges levels with fewer than two lines into the reference. The reviewer reproduced both outcomes of the defect:
- *Failed fit* (H3b, WGD, Z0026): the all-model reference lineage (Adrenal Gland) has no L line, and Testis has a single complete-case L line. Testis therefore became its own reference group with leverage 1, and HC3 refused the fit. With the design assembled over L, the fit gives β = −0.0241, p = 0.880 (not confirmed).
- *Shifted fit* (H3b, TP53_altered, Z0096): two singleton lineages form a two-line group, which adds one design column (32 instead of 31). The primary gives p = 0.0506 (not confirmed); the corrected design gives p = 0.0490 (confirmed). This is one of the three defect-sensitive flips.

<!-- RV §2 "Defect diagnosis, reproduced under both designs" -->

Over all 3,126 hits of the 30 entries, a span audit using labels and the L row index only predicts exactly 510 shifted fits, 50 failures due to the defect and 3 genuine failures. 2,504 hits have an identical span. <!-- RV §2 "Full-scale label audit"; RVD/span_audit.tsv -->

*Affected fits by entry* (span differs under the two constructions): all 122 successful H3b fits, 17 H1 axes (16 N axes plus EP300 P), 220 H1S, 7 H2, 120 BRIDGE-L GLS hits and 11 controls. H3a is unaffected. <!-- LB/CORRECTED_REFIT_INTENT.md; RV N1 --> The relook's intent note named ARID1A P instead of EP300 P; the count was right (erratum, D14).

*Defect-sensitive calls (53).* These comprise H3b 30 (28 refitted, 2 flips), BRIDGE-L 10 (9 refitted, 1 flip), controls 11 (PC9 ×3, PC10 ×3, PC12 ×2, PC13 ×2, PC18), H1 1 (FOCAD N) and H2 1 (Z0026.2 in WGD). <!-- LB/defect_sensitive_calls_R7.tsv; LD: defect_sensitive_by_entry --> In the corrected relook, the 50 refitted calls go as follows: H3b 9 confirmed and 19 not; BRIDGE-L 6 and 3; controls 10 and 1 (PC18); H1 FOCAD confirmed; H2 confirmed. <!-- RV §3 "Defect-sensitive list" --> None of these calls supports a claim.

*Reported rates.* Table S6 gives every rate in three versions: primary evaluable, primary all-in and corrected. The H3b NOT EVALUABLE fits are informative, not missing at random (26 of 28 in WGD), so the evaluable rate is biased upward. <!-- RV MUST-1 -->

## Supplementary Note 4. Post hoc power aid

This analysis is not registered, was computed after the read by the independent reviewer, and is not a test. <!-- RV §4 --> For each eligible hit on the lockbox route, it computes the power of the one-sided α = 0.05 L test (noncentral t, L Satterthwaite df). The expected effect is the A estimate under four versions: unshrunk, halved, the empirical-Bayes posterior mean (Tweedie's formula on all evaluable A z-scores of the direction: 8,441 for H3a, 7,340 for H3b), and that posterior integrated over its variance. The expected SE is either the A SE scaled by the L/A line counts or the observed L SE (median ratio 0.99 and 1.04). <!-- RV §4 "Method"; RVD/tables_review.json: power -->

| Direction | Observed | Unshrunk | 50% shrink | EB mean | EB mean + variance |
|---|---|---|---|---|---|
| H3a (1,145) | 0.525 | 0.69 | 0.32–0.33 | 0.51–0.52 | 0.47 |
| H3b (150) | 0.49 evaluable / 0.46 corrected / 0.40 all-in | 0.68–0.71 | 0.32–0.33 | 0.43–0.44 | 0.43–0.44 |

<!-- RVD/tables_review.json: power.*.expected_rate (rounded); CMP -->

The median EB shrinkage factor is 0.74 (H3a) and 0.59 (H3b). Split by the complex-specificity flag, the expected rates are 0.52–0.53 against 0.41–0.43 for H3a (observed 0.572 vs 0.489), and 0.50–0.51 against 0.29–0.30 for H3b (observed 0.583 / 0.511 vs 0.341 / 0.360). <!-- RVD/tables_review.json: power_by_beats_all; LD: power_band --> The caveats are that the tests are correlated, the prior is pooled across contexts, the t→z conversion is approximate, and permutation p-values are approximated by t-tests.

## Supplementary Note 5. BENCH and BRIDGE-L estimands

*BENCH.* For each primary pair (an H3a or H3b discovery in a lockbox-testable context), g* is the member with the largest A |t| in the complex's direction, chosen in A after the context's removals (ties broken by |β|). Δ = t_L(S_C) − t_L(g*), signed toward the A direction. Breadth is the fraction of all members, removed ones included, with one-sided p < 0.05 in A (literal reading, A5). Terciles are formed on the primary pairs' breadth with ties broken by pair order (A5, reading 11). The 95% CI is a PatientID-blocked bootstrap over L lines with pairs and g* fixed (B = 2,000, seed 5100). The registered criterion is the high-tercile CI. <!-- docs/registry/BENCH_complex_vs_gene.yaml: estimand, success_criterion; AMENDMENTS.md A5 --> The pairs cover 1,295 (C, c) combinations in 353 complexes with 920 distinct g*. <!-- LD: BENCH --> The bootstrap failure rate (16.5%) and the sensitivity of the CI to covariate assembly are in main text §7.6. Not run: the within-line-rank multiverse score (A3) and the A–L random-set correlation (A4, A5).

*BRIDGE-L.* GLS-map discovery used the H3a/H3b machinery on A: per-context covariates, 12/12 evaluability, a screen at q < 0.1, 10,000 Freedman–Lane permutations, seed 4500 and BB at q = 0.05 in two calls. It gave 1,206 lineage and 139 driver discoveries. <!-- AMENDMENTS.md A4 reading 4; LB/BRIDGE-L_confirmation_rate/gls_discoveries.tsv --> The zavg statuses come from the H3a and H3b confirmation tables, so no new zavg read was made. Per map, the confirmation rate is computed within the size bins {3–4, 5–8, 9–16, 17–60} and averaged with the pooled bin weights of the discovered clusters. Weights are renormalised over the bins present (A5, N12), and rare-context hits are excluded. The CI is a bootstrap over hits within map and bin (B = 2,000, seed 5200). <!-- docs/registry/BRIDGE-L_confirmation_rate.yaml; AMENDMENTS.md A4, A5 -->

## Supplementary Note 6. Registry summary

Thirty entries, all closed. "m" is the registered test count or the size of the target set. A-side and lockbox results are summarised from the committed outputs. <!-- docs/registry/*.yaml: id, direction, family_m, status; CMP; LB/primary_relabelled_R7.tsv.gz; results/part2/controls/kill_rule.json -->

| Entry | Direction | m | Standing | A side | Lockbox (primary; corrected in brackets) |
|---|---|---|---|---|---|
| H3a_lineage_atlas | H3a | 8,441 | confirmatory | 1,145 discoveries | 601/1,145 confirmed [601/1,145] |
| H3b_driver_family | H3b | 7,340 | confirmatory | 161 discoveries | 60/122 evaluable, 60/150 all-in [69/150]; rare route 5 pass, 6 NE |
| H1_within_complex_GI | H1 | 28 | confirmatory | 5 discoveries; 3 headline-eligible | ARID1A P confirmed (= PC6); EP300 P rare route passed; FOCAD N NE (defect) [confirmed; no claim]; TP53 N confirmed, not promotable; PTEN N not confirmed |
| H1S_silencing_GI | H1S | 338 | exploratory | 31 discoveries | 9/17 evaluable [9/17]; rare route 6/6 pass |
| H2_submodule_divergence | H2 | set by H3 discoveries | exploratory (no evaluable control) | 109 discoveries | 76/108 evaluable [77/109] |
| BENCH_complex_vs_gene | benchmark | 1,295 pairs | confirmatory (estimation) | — | high tercile Δ +0.240 (−0.033 to +0.362): complex did not out-replicate best gene |
| BRIDGE-L_confirmation_rate | bridge | 16,294 (GLS map) | secondary bridge | 1,345 GLS-map discoveries | Δ −0.065 (−0.107 to −0.025): favours GLS |
| CTRL_PC1 (WRN, MSI-H) | H3b | 2 | control | recovered | 2/2 |
| CTRL_PC2 (MDM2/PPM1D/USP7, TP53 WT) | H3b | 4 | control | recovered | 4/4 |
| CTRL_PC3 (TP53/CDKN1A/USP28) | gene level | 3 | control, not counted | gene level | 3/3 |
| CTRL_PC4 (SHOC2, RAS) | H3b | 2 | control | recovered | 2/2 |
| CTRL_PC5 (SHOC2, BRAF vs RAS) | H3b | 2 | control | recovered | 1/2 (unit p_L 0.053) |
| CTRL_PC6 (ARID1B, ARID1A-LoF) | H1 | 2 | control | recovered | 2/2 |
| CTRL_PC7 (SMARCA2, SMARCA4-LoF) | H1 | 1 | control | gene only | 1/1 |
| CTRL_PC8 (STAG1, STAG2-LoF) | H1 | 1 | control | gene only | 1/1 |
| CTRL_PC9 (PRMT5/WDR77, 9p21) | H3b | 3 | control | recovered | 3 NE (defect) [3/3; no claim] |
| CTRL_PC10 (PELO/HBS1L, FOCAD or MSI-H) | H3b | 3 | control | recovered | 3 NE (defect) [3/3; no claim] |
| CTRL_PC11 (NFE2L2, NRF2) | H3b | 2 | control | recovered | 2/2 |
| CTRL_PC12 (VPS4A, VPS4B loss) | H1 | 2 | control | recovered | 2 NE (defect) [2/2; no claim] |
| CTRL_PC13 (DDX3X, female or LoY) | H3b | 2 | control | recovered | 2 NE (defect) [2/2; no claim] |
| CTRL_PC14 (CDK4/6, CCND1, RB1-LoF) | gene level | 3 | control, not counted | gene level | 2/3 (CDK6 p_L 0.062) |
| CTRL_PC15 (BRD9, SMARCB1-LoF or SS18::SSX) | H1 | 2 | control | recovered | rare route: unit passed; BRD9 NE (DEMETER2) |
| CTRL_PC16 (MEN1/MLLT1, KMT2A-r) | gene level | 2 | control, not counted | gene level, mixed | rare route: MEN1 passed; MLLT1 failed |
| CTRL_PC17 (SOX10/MITF, melanoma) | H3a | 3 | control | recovered | 3/3 |
| CTRL_PC18 (KIF18A, WGD) | H3b | 1 | control | not evaluable at complex level | 1 NE (defect) [not confirmed; no claim] |
| CTRL_PC19a (PIK3CB, PTEN-LoF) | H3b | 2 | control | recovered | 2/2 |
| CTRL_PC19b (PIK3CA hotspot) | H3b | 2 | control | recovered | 2/2 |
| CTRL_PC20 (MYB, heme) | H3a | 2 | control | recovered | 2/2 |
| CTRL_PC21 (CDK6/CCND3 sub-module, heme vs breast) | H2 | 2 | control | not evaluable | 2 NE (genuine) |
| CTRL_PC22 (ncBAF sub-module) | H2 | 1 | control | not evaluable | rare route NE (DEMETER2) |

The kill rule, which is defined on A, passes for H3a (2/2), H3b (10/10) and H1 (3/3). H2 has no evaluable control and is exploratory. <!-- results/part2/controls/kill_rule.json -->

## Supplementary Note 7. Amendments, errata and rulings

The full records are in `docs/registry/AMENDMENTS.md`. Each touched registry entry carries an `amendments` item, and a test checks that every changed field is declared. <!-- docs/registry/AMENDMENTS.md header --> The main text gives the deviations table with timing relative to outcome access. In brief:
- **A1** (28 Oct 22:54:50; before any H1 run; labels only): 318 silencing events moved to the exploratory H1S.
- **A2** (23:28:31; before any H1 run): genetic-loss rule (≥ 50% LoF or hom-del); SMARCB1 moved to H1.
- **E1** (29 Oct 00:41:42; after the H1 run): sign typo in the aligned score corrected to the stated semantics, which the engine had used since 28 Oct 20:09.
- **E2** (01:07; after the H1 run and its review; post hoc, conservative): "detectable, within equivalence margin" label; such events cannot be promoted.
- **A3** (02:03:27; before any L read): lockbox permutation seeds; BENCH within-line rank not run.
- **A4** (02:42:37; before any L read): BENCH secondaries not run; multiverse point estimates; H1S rare-route N axes not evaluable (later withdrawn).
- **A5** (03:25:03; before any L read): A4 erratum; distinct-complex count reported; H1S rare route on A alignment; BRIDGE-L multiverse CIs; readings 8–11.
- **R1–R5** (29 Oct ~00:41–01:07; after the A-side discovery runs, before any L read): E1; the engine crash fix (R2); GPD p primary with a pre-declared conservative sensitivity (R3); calibrations implemented after discovery (R4); E2 (R5).
- **R6, R7** (10:19:44, 10:41:49; after the read): see Supplementary Note 3.

<!-- docs/registry/AMENDMENTS.md; RESEARCH_TREE.md D7–D14; LB/README.md -->

## Supplementary Note 8. Member-resolved H1 profiles (A, descriptive)

Loss coefficient β_j with 95% CI for every member of the three headline-eligible events (Chronos units; < 0 = member more essential in loss lines). Role: P = in-complex paralog; N' = counted member with a baseline direction; nb = no baseline direction (not counted). The full export for all events is `results/part2/directions/H1_within_complex_GI/member_profiles.tsv`. <!-- DR/H1_within_complex_GI/member_profiles.tsv -->

| Event (loss lines) | Member | Role | β_j [95% CI] |
|---|---|---|---|
| ARID1A (37) | ARID1B | P | −0.374 [−0.512, −0.237] |
| | SMARCB1 | N' | +0.109 [−0.037, +0.254] |
| | SMARCC1 | N' | +0.087 [+0.006, +0.167] |
| | SMARCA4 | N' | −0.035 [−0.150, +0.081] |
| | SMARCC2 | N' | +0.032 [−0.022, +0.086] |
| | SMARCE1 | N' | +0.025 [−0.083, +0.132] |
| | SMARCD2 | N' | +0.011 [−0.050, +0.073] |
| | DPF2 | N' | +0.011 [−0.077, +0.099] |
| | IER5L | nb | −0.076 [−0.147, −0.005] |
| | SS18 | nb | −0.015 [−0.085, +0.055] |
| EP300 (12) | CREBBP | P | −1.093 [−1.462, −0.725] |
| FOCAD (30) | PELO | N' | −0.808 [−1.084, −0.531] |
| | HBS1L | N' | −0.334 [−0.467, −0.202] |
| | AVEN | N' | +0.091 [+0.023, +0.159] |
| | SKIC2 | N' | −0.089 [−0.156, −0.021] |
| | RNF14 | N' | +0.053 [−0.005, +0.112] |
| | VDAC1 | N' | +0.053 [−0.050, +0.155] |
| | DIS3L | N' | −0.045 [−0.107, +0.018] |
| | SKIC3 | N' | −0.039 [−0.106, +0.028] |
| | TUT7 | nb | −0.175 [−0.318, −0.032] |

<!-- DR/H1_within_complex_GI/member_profiles.tsv rows Z0109|ARID1A, Z0186|EP300, Z0286|FOCAD: beta, ci95_lo, ci95_hi, in_aligned_set -->

EP300's complex-wide axis is not evaluable (only SMARCD1 has a baseline direction), so only its paralog is exported. The export was committed at 06:53:55 UTC on 29 October, during the lockbox read window. It reads A only and changes no call (main text, Methods).

## Supplementary Note 9. Reproduction

From the repository root, with `PART1_DIR` set and `PYTHONHASHSEED=0`:

```
.venv/bin/python scripts/tab_manuscript.py
.venv/bin/python scripts/tab_part2_discoveries.py
.venv/bin/python scripts/tab_lockbox.py
.venv/bin/python scripts/fig_1_test_evaluation.py
.venv/bin/python scripts/fig_2_construction.py
.venv/bin/python scripts/fig_3_context_atlas.py
.venv/bin/python scripts/fig_4_h1_panel.py
.venv/bin/python scripts/fig_5_lockbox.py
```

These scripts read committed result files only, never datasets, data caches or persisted lockbox reads. A second run reproduces every figure file byte for byte. The full pipeline (scores, maps, contexts, directions, calibrations and the lockbox read) is in `scripts/` in numbered order, with each step's command in its intent note. Re-running the lockbox read is blocked by the ledger by design.

## Supplementary Tables

| Table | Content | File |
|---|---|---|
| S1 | Pre-registered controls (A) and kill-rule tally | `manuscript/tables/tableS1_controls.md`, `tableS1b_kill_rule.tsv` |
| S2 | Largest T1 non-locus complexes with post hoc HGNC gene-group labels (for example Z0001 iron–sulfur assembly, Z0002 mitochondrial tRNA synthetases with mtDNA replication and transcription, Z0003 cytochrome c oxidase assembly, Z0004 centriole duplication, Z0008 spindle checkpoint with MIS12, Z0010 TRAPP/COPII, Z0012 complex I; labels are post hoc) | `results/part2/construction/summary/zavg_top_T1_complexes.tsv` |
| S3 | All 1,306 H3a/H3b discoveries with specificity, locus, GPD and screen-quality flags | `manuscript/tables/tableS3_context_discoveries.md` (full `.tsv`) |
| S4 | H1 panel, 24 events | `manuscript/tables/tableS4_h1_panel.md` (full `.tsv`) |
| S5 | H2 discoveries with rank against own random splits | `manuscript/tables/tableS5_h2_split_ranks.md` (full `.tsv`) |
| S6 | Lockbox confirmation rates per direction: evaluable / all-in / corrected, defect-sensitive counts, power expectation | `manuscript/tables/tableS6_lockbox_rates.md` |
| S7 | Lockbox results for all H3a/H3b discoveries; strongest specific findings | `manuscript/tables/tableS7_lockbox_h3.md` (full `tableS7_lockbox_h3_full.tsv`) |
| S8 | Lockbox results for H1, H2, BENCH and BRIDGE-L | `manuscript/tables/tableS8_lockbox_h1_h2_bench_bridge.md` |

<!-- S2 labels: results/part2/construction/summary/zavg_top_T1_complexes.tsv: POSTHOC_hgnc_gene_group; main.md §2.3 -->
