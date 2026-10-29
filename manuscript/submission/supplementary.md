# Supplementary Information

**A co-essentiality score selected on training complexes wins a single held-out test but yields neither better-supported complexes nor more replicable complex-level dependencies in DepMap**

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
- Supplementary Note 10. Extended test-evaluation results and calibration details
- Supplementary Note 11. Extended construction results
- Supplementary Note 12. Per-direction results (A side and lockbox)
- Supplementary Note 13. Sensitivity analyses
- Supplementary Note 14. Deviations from pre-registration (full table)
- Supplementary Note 15. Extended limitations
- Supplementary Note 16. The registered strongest findings
- Supplementary Note 17. Revision analyses R-1 to R-5 (post hoc, descriptive) and requests not carried out
- Supplementary Note 18. Glossary of identifiers
- Supplementary Tables S1–S9 (files)

---

## Supplementary Note 1. Transition from development to test

The Part 1 method was locked at commit c7f95ee (28 October 2026, 12:55:43 UTC). The hu.MAP 3.0 test files appeared in the project folder at 13:55 UTC, 59.6 min later. The only contacts with the test files before the protocol freeze were the following.

| ID | Who | What | Consequence |
|---|---|---|---|
| T0 / D0 | PI | Line counts (3,312 / 597,426) and the first 3 lines of each file, as a format check | None |
| T1 | Part 1 audit agent | Descriptive counts without labels: overlaps, strata, mappability | Known when the protocol was written; the choices follow Part 1's design |
| D1 | Part 1 audit agent, intermediate run | Counted test pairs in Part 1's unsupervised top-pair table (zavg_BEF ≥ 3), which is equivalent to precision and recall at z ≥ 3 | Counts deleted; PI never saw them; rule (a) reported only as "non-blind (D1)" and enters no verdict |
| T2 | Protocol author | Octal dump of the first ~200 bytes of the positive file | None |
| T3 | Everyone | Genome-wide unsupervised scores covering all test pairs existed | The single-look guarantee is procedural |

Decision notes D2–D4 concern data scope and are not contacts. D2: no CORUM or Complex Portal was obtained, so any naming of complexes is post hoc. D3: the locked verdict is computed on S1, and ITT is quoted alongside. D4: the scope of Part 2 was cut, and BioPlex is used for validation only. An independent audit found no evidence that test information entered Part 1 (docs/part1_audit.md).

*Wording.* The first version said the score was "locked before the test data existed". The hu.MAP 3.0 test split was published earlier [Fischer 2025]; what is true is that the method was locked 59.6 min before the test files reached this project and was never run on them before the single look.

*Order of Part 2 design relative to the test (28 October 2026, UTC).*

| Time | Event | Commit |
|---|---|---|
| 12:55:43 | Part 1 method locked | c7f95ee |
| 13:55 | Test files reach the project | – |
| 15:08 | DESIGN v1 (Part 2) | d6bf744 |
| 17:37 | DESIGN v3.1, final after three review rounds | bdb8816 |
| 18:14 | Test protocol frozen | 039a843 |
| 18:18 | Test labels first read (the single look) | looks.jsonl |
| 18:52 | Test look results committed | 0861e59 |
| 19:27 | Part 2 registry draft (29 entries) | ee45596 |

The Part 2 design was written after the test files arrived, and the Part 2 registry after the Part 1 result was known (zavg_BEF had won). Nothing indicates that Part 2 used test information, and Part 2 code cannot read the test files (firewall tests); Part 2 uses zavg_BEF regardless of the test outcome (DESIGN §1).

## Supplementary Note 2. Methods detail

*Construction.* The complete parameters of both maps are in `results/part2/construction/{zavg_map,gls_map}/parameters.json`, and the rationale for each construction rule is in `results/part2/construction/DECISIONS.md` (D1–D6):
- containment support for small complexes (D1);
- placement and trimming of overlapping candidates (D2);
- attachment of genes outside S (D3);
- edge-count matching of the GLS comparator (D4);
- tiers and the bridge endpoint (D5);
- split-half selection (D6).

A toy simulation fixed the rules before any real map was built: 86 of 105 planted complexes and all 39 planted sub-modules were recovered, and null data gave no cluster.

*Association engine.* Complex scores, the OLS model with HC3 errors and Satterthwaite df, the PatientID-blocked Freedman–Lane permutation, Benjamini–Bogomolov selective FDR and stage-wise H2 testing are implemented in `src/part2/` and specified per entry in `docs/registry/*.yaml`. The engine's lineage-merge convention merges indicator levels with fewer than two lines into the reference. DESIGN does not state it; it changes p by at most 2% (relative) in the controls. In 1,000 simulations with planted lineage effects, the selective error rate was 0.046 and the pooled false-discovery proportion 0.068. References for these methods are given in the main-text Methods.

*Registered calibrations (D8, R4).* The matched random sets, never-expressed sets, copy-number segment null and H2 random splits were specified operationally in `results/part2/directions/CALIBRATION_D8_INTENT.md` before they were run. They were run once each, after the discovery runs and before any lockbox read.

*Lockbox driver.* `scripts/51_lockbox_confirm.py` reads each registered entry once through `assoc.read_lockbox`, which commits the ledger line before it loads and refuses a second read. It persists the read with its hash, then applies the registered confirmation rule. Hits in contexts with fewer than 5 context or 5 reference L lines take the rare-context route, which uses only A lines and DEMETER2. Lockbox permutation seeds are registered by A3, and the BENCH and BRIDGE-L bootstrap seeds by the entries (5100, 5200).

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

*Diagnosis.* The driver assembled the one-hot covariate design over all models and then subset it to L. The engine instead assembles the design over the rows analysed and merges levels with fewer than two lines into the reference. The reviewer reproduced both outcomes of the defect:
- *Failed fit* (H3b, WGD, Z0026): the all-model reference lineage (Adrenal Gland) has no L line, and Testis has a single complete-case L line. Testis therefore became its own reference group with leverage 1, and HC3 refused the fit. With the design assembled over L, the fit gives β = −0.0241, p = 0.880 (not confirmed).
- *Shifted fit* (H3b, TP53_altered, Z0096): two singleton lineages form a two-line group, which adds one design column (32 instead of 31). The primary gives p = 0.0506 (not confirmed); the corrected design gives p = 0.0490 (confirmed). This is one of the three defect-sensitive flips.

Over all 3,126 hits of the 30 entries, a span audit using labels and the L row index only predicts exactly 510 shifted fits, 50 failures due to the defect and 3 genuine failures. 2,504 hits have an identical span.

*Affected fits by entry* (span differs under the two constructions): all 122 successful H3b fits, 17 H1 axes (16 N axes plus EP300 P), 220 H1S, 7 H2, 120 BRIDGE-L GLS hits and 11 controls. H3a is unaffected. The relook's intent note named ARID1A P instead of EP300 P; the count was right (erratum, D14).

*Defect-sensitive calls (53).* These comprise H3b 30 (28 refitted, 2 flips), BRIDGE-L 10 (9 refitted, 1 flip), controls 11 (PC9 ×3, PC10 ×3, PC12 ×2, PC13 ×2, PC18), H1 1 (FOCAD N) and H2 1 (Z0026.2 in WGD). In the corrected relook, the 50 refitted calls go as follows: H3b 9 confirmed and 19 not; BRIDGE-L 6 and 3; controls 10 and 1 (PC18); H1 FOCAD confirmed; H2 confirmed. None of these calls supports a claim.

*Reported rates.* Table S6 gives every rate in three versions: primary evaluable, primary all-in and corrected. The H3b NOT EVALUABLE fits are informative, not missing at random (26 of 28 in WGD), so the evaluable rate is biased upward.

## Supplementary Note 4. Post hoc power aid

This analysis is not registered, was computed after the read by the independent reviewer, and is not a test. For each eligible hit on the lockbox route, it computes the power of the one-sided α = 0.05 L test (noncentral t, L Satterthwaite df). The expected effect is the A estimate under four versions: unshrunk, halved, the empirical-Bayes posterior mean (Tweedie's formula on all evaluable A z-scores of the direction: 8,441 for H3a, 7,340 for H3b), and that posterior integrated over its variance. The expected SE is either the A SE scaled by the L/A line counts or the observed L SE (median ratio 0.99 and 1.04).

| Direction | Observed | Unshrunk | 50% shrink | EB mean | EB mean + variance |
|---|---|---|---|---|---|
| H3a (1,145) | 0.525 | 0.69 | 0.32–0.33 | 0.51–0.52 | 0.47 |
| H3b (150) | 0.49 evaluable / 0.46 corrected / 0.40 all-in | 0.68–0.71 | 0.32–0.33 | 0.43–0.44 | 0.43–0.44 |

The median EB shrinkage factor is 0.74 (H3a) and 0.59 (H3b). Split by the complex-specificity flag, the expected rates are 0.52–0.53 against 0.41–0.43 for H3a (observed 0.572 vs 0.489), and 0.50–0.51 against 0.29–0.30 for H3b (observed 0.583 / 0.511 vs 0.341 / 0.360). The caveats are that the tests are correlated, the prior is pooled across contexts, the t→z conversion is approximate, and permutation p-values are approximated by t-tests. References for these methods are given in the main-text Methods.

## Supplementary Note 5. BENCH and BRIDGE-L estimands

*BENCH.* For each primary pair (an H3a or H3b discovery in a lockbox-testable context), g* is the member with the largest A |t| in the complex's direction, chosen in A after the context's removals (ties broken by |β|). Δ = t_L(S_C) − t_L(g*), signed toward the A direction. Breadth is the fraction of all members, removed ones included, with one-sided p < 0.05 in A (literal reading, A5). Terciles are formed on the primary pairs' breadth with ties broken by pair order (A5, reading 11). The 95% CI is a PatientID-blocked bootstrap over L lines with pairs and g* fixed (B = 2,000, seed 5100). The registered criterion is the high-tercile CI. The pairs cover 1,295 (C, c) combinations in 353 complexes with 920 distinct g*. The bootstrap failure rate (16.5%) and the sensitivity of the CI to covariate assembly are in the main text (BENCH) and Supplementary Note 12. Not run: the within-line-rank multiverse score (A3) and the A–L random-set correlation (A4, A5). Because 16.5% of bootstrap fits failed silently and the lower bound moves from −0.033 to −0.122 with the design assembled within replicates, the revised manuscript treats the registered interval as unreliable and the BENCH result as inconclusive, with a point estimate favouring the complex in the high tercile; the first version's conclusion that the complex "does not beat its best gene" is withdrawn.

*BRIDGE-L.* GLS-map discovery used the H3a/H3b machinery on A: per-context covariates, 12/12 evaluability, a screen at q < 0.1, 10,000 Freedman–Lane permutations, seed 4500 and BB at q = 0.05 in two calls. It gave 1,206 lineage and 139 driver discoveries. The zavg statuses come from the H3a and H3b confirmation tables, so no new zavg read was made. Per map, the confirmation rate is computed within the size bins {3–4, 5–8, 9–16, 17–60} and averaged with the pooled bin weights of the discovered clusters. Weights are renormalised over the bins present (A5, N12), and rare-context hits are excluded. The CI is a bootstrap over hits within map and bin (B = 2,000, seed 5200).

## Supplementary Note 6. Registry summary

Thirty entries, all closed. "m" is the registered test count or the size of the target set. A-side and lockbox results are summarised from the committed outputs.

| Entry | Direction | m | Standing | A side | Lockbox (primary; corrected in brackets) |
|---|---|---|---|---|---|
| H3a_lineage_atlas | H3a | 8,441 | confirmatory | 1,145 discoveries | 601/1,145 confirmed [601/1,145] |
| H3b_driver_family | H3b | 7,340 | confirmatory | 161 discoveries | 60/122 evaluable, 60/150 all-in [69/150]; rare route 5 pass, 6 NE |
| H1_within_complex_GI | H1 | 28 | confirmatory | 5 discoveries; 3 headline-eligible | ARID1A P confirmed (= PC6); EP300 P rare route passed; FOCAD N NE (defect) [confirmed; no claim]; TP53 N confirmed, not promotable; PTEN N not confirmed |
| H1S_silencing_GI | H1S | 338 | exploratory | 31 discoveries | 9/17 evaluable [9/17]; rare route 6/6 pass |
| H2_submodule_divergence | H2 | set by H3 discoveries | exploratory (no evaluable control) | 109 discoveries | 76/108 evaluable [77/109] |
| BENCH_complex_vs_gene | benchmark | 1,295 pairs | confirmatory (estimation) | — | high tercile Δ +0.240 (−0.033 to +0.362): inconclusive (bootstrap partly failed; interval unreliable) |
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

The kill rule, which is defined on A, passes for H3a (2/2), H3b (10/10) and H1 (3/3). H2 has no evaluable control and is exploratory.

## Supplementary Note 7. Amendments, errata and rulings

The full records are in `docs/registry/AMENDMENTS.md`. Each touched registry entry carries an `amendments` item, and a test checks that every changed field is declared. The main text gives a short table (Table 3), and Supplementary Note 14 the full deviations table with timing relative to outcome access. In brief:
- **A1** (28 Oct 22:54:50; before any H1 run; labels only): 318 silencing events moved to the exploratory H1S.
- **A2** (23:28:31; before any H1 run): genetic-loss rule (≥ 50% LoF or hom-del); SMARCB1 moved to H1.
- **E1** (29 Oct 00:41:42; after the H1 run): sign typo in the aligned score corrected to the stated semantics, which the engine had used since 28 Oct 20:09.
- **E2** (01:07; after the H1 run and its review; post hoc, conservative): "detectable, within equivalence margin" label; such events cannot be promoted.
- **A3** (02:03:27; before any L read): lockbox permutation seeds; BENCH within-line rank not run.
- **A4** (02:42:37; before any L read): BENCH secondaries not run; multiverse point estimates; H1S rare-route N axes not evaluable (later withdrawn).
- **A5** (03:25:03; before any L read): A4 erratum; distinct-complex count reported; H1S rare route on A alignment; BRIDGE-L multiverse CIs; readings 8–11.
- **R1–R5** (29 Oct ~00:41–01:07; after the A-side discovery runs, before any L read): E1; the engine crash fix (R2); GPD p primary with a pre-declared conservative sensitivity (R3); calibrations implemented after discovery (R4); E2 (R5).
- **R6, R7** (10:19:44, 10:41:49; after the read): see Supplementary Note 3.
- **D15, D16** (after the read; A only): the registered per-hit robustness checks were computed and the strongest findings fixed by the registered rule, with post hoc descriptive annotations (Supplementary Note 16).

## Supplementary Note 8. Member-resolved H1 profiles (A, descriptive)

Loss coefficient β_j with 95% CI for every member of the three headline-eligible events (Chronos units; < 0 = member more essential in loss lines). Role: P = in-complex paralog; N' = counted member with a baseline direction; nb = no baseline direction (not counted). The full export for all events is `results/part2/directions/H1_within_complex_GI/member_profiles.tsv`.

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
.venv/bin/python scripts/tab_strongest.py
.venv/bin/python scripts/check_manuscript_numbers.py manuscript/main.md
```

These scripts read committed result files only, never datasets, data caches or persisted lockbox reads. A second run reproduces every figure file byte for byte. The full pipeline (scores, maps, contexts, directions, calibrations and the lockbox read) is in `scripts/` in numbered order, with each step's command in its intent note. Re-running the lockbox read is blocked by the ledger by design.

## Supplementary Note 10. Extended test-evaluation results and calibration details

*Inference.* The Holm-adjusted p-values for the three primary ΔAP comparisons were < 0.003 (GLS; no replicate ≤ 0 of 2,000, so the two-sided bootstrap p is < 0.001 and its Holm adjustment < 0.003; the first version's "< 0.0005" gave the one-sided bound without the factor 3), 0.022 (Spearman) and 0.018 (Pearson). The intervals were stable across Monte-Carlo seeds 1–4 (lower bound of the GLS comparison 0.078–0.082). The secondary protein bootstrap (+0.093 to +0.154) and the anti-conservative pair bootstrap (+0.111 to +0.132) gave narrower intervals for the same comparison. AUROC was 0.789 (0.735–0.829); the gain over GLS was +0.059 (+0.030 to +0.090), while the interval for the gain over Spearman included 0 (+0.033, −0.003 to +0.067).

*Locked verdict on S1.* S1 has 3,159 positives and 563,648 negatives. On S1, AP was 0.349 (0.274–0.435), and the gains were +0.127 (+0.089 to +0.167), +0.135 (+0.024 to +0.229) and +0.143 (+0.038 to +0.233) against GLS, Spearman and Pearson. The macro per-cluster AP ratio against GLS was 1.017 (0.955–1.161), which passes the 0.95 non-inferiority margin but shows no per-cluster advantage.

*Secondary candidate and post hoc reference.* The pre-declared secondary candidate B (the second-order component alone) did not meet the rule: its macro per-cluster AP was 0.875× that of GLS on S1, below the 0.95 margin, as in the Part 1 internal test, and its AP was below zavg_BEF's (−0.022, −0.032 to −0.012). A post hoc PC-removed Pearson reference, which is not a candidate, reached AP 0.239, below zavg_BEF (−0.094, −0.126 to −0.061).

*Precision–recall curves and operating points.* The curves (Fig. 1a) come from a descriptive export of the single look; the export is not a relook, and its integrated AP reproduces the look's AP to ≤ 10⁻¹⁶. The curves coincide up to recall 0.05. Above recall 0.1, zavg_BEF has the highest precision of the four methods; at recall 0.3–0.4 interpolated precision was 0.60 and 0.44 for zavg_BEF, against at most 0.24 and 0.09 for the baselines. zavg_BEF's in-sample training F1 at the rule (b) threshold was 0.478, against 0.411 on test. At an equal prediction budget (rule c, the top 0.522% of pairs), F1 was 0.415 for zavg_BEF and 0.268 for GLS. At z ≥ 3 (rule a, non-blind owing to D1), precision was 0.286 and recall 0.460.

*Calibration.* For zavg_BEF the Brier skill score was 0.222 (GLS 0.139, Spearman 0.089, Pearson 0.097), the calibration slope 1.007 and calibration-in-the-large 0.893. Calibration-in-the-large is Σp/Σy; the protocol pre-specified no bootstrap CI for it, so none is reported, and no recalibration was made. The Spearman and Pearson calibrators over-predicted (calibration-in-the-large 1.195 and 1.166; slopes 0.84 and 0.87).

*Size breakdown and strata.* By positive pairs per cluster, macro AP differed from GLS by −0.033 (1 pair; 91 clusters), −0.024 (2–5; 67), +0.041 (+0.008 to +0.080; 6–20), +0.126 (+0.075 to +0.182; 21–100) and +0.142 (+0.078 to +0.198; 101+, 6 clusters). These CIs resample the bin's clusters and are not the protocol bootstrap. Without the largest component (39 proteins) AP was 0.330 against 0.209. In C1 and C2 the gains over GLS were +0.081 (+0.026 to +0.140) and +0.050 (+0.017 to +0.092). Among S1 pairs of two common-essential genes AP was 0.511 against 0.363; among pairs of two non-essential genes, 0.197 against 0.169 (no CIs).

*Runtime.* Scoring all 600,738 test pairs took 112.6 s for zavg_BEF and 7.5–17.4 s for the baselines on a laptop; the full evaluation with all bootstraps took 1,998 s.

## Supplementary Note 11. Extended construction results

*Coverage.* At the same edge count the zavg map clusters 27% more genes than the GLS map, with higher coverage of pan-essential (62% vs 42%), selective (76% vs 66%) and other S genes (48% vs 34%) (Fig. 2b). 2,380 of the 6,022 S genes are in no level-1 complex, including SMARCA2, STAG1 and KIF18A. Both maps were byte-identical when rebuilt in a fresh process by an independent reviewer.

*Selection.* At k = 10, z0 = 4 gave 0.341 (SD 0.010), z0 = 5 0.338 (SD 0.019) and z0 = 6 0.343 (SD 0.005); k = 20 gave 0.204–0.256. A degree-rewired graph scored 0 at every grid point. GLS's own split-half preference, (4, 10) at 0.315, is reported only.

*Attachments.* Genes outside S were attached when at least half the complex's members exceed z0 with the gene and a hypergeometric test passes BH q < 0.01; 2,563 attachments covering 2,000 genes were made among 5.7 million tests (about 26 expected false). Attached genes are not members; 485 non-locus level-1 attachments may define H1 events.

*Content.* The largest non-locus complexes supported by BioPlex include, by post hoc reading of their member symbols, cytosolic iron–sulfur assembly, mitochondrial aminoacyl-tRNA synthetases with mtDNA replication and transcription factors, cytochrome c oxidase assembly, centriole duplication, the spindle assembly checkpoint with the MIS12 kinetochore complex, TRAPP and COPII trafficking, and complex I (Supplementary Table S2). Some level-1 complexes merge two complexes through shared or hub genes (for example PRC1 with PRC2, Z0061), and some T1 modules are pathways (for example N-glycosylation, Z0276).

*Tiers.* BioPlex 3.0 (293T ∪ HCT116; 167,374 edges, observed for 5,267 of the 6,022 S genes) and CCLE proteomics (239 A lines) were read only after the map commit. About 8 T1 calls are expected false. T2 is not locus-aware: 14 zavg locus clusters are T2 because copy-number dosage makes neighbouring proteins co-abundant. An independent implementation assigned the same tier to 494 of 499 zavg clusters.

*Variability.* 73% of complex-by-map Jaccard values were ≥ 0.5, and random size-matched sets reached 0.5 once in 734,000 draws. Recovery rose with size, from a median of 0.50 for 3–4-gene complexes (52% ≥ 0.5) to 0.82 for 17–60-gene complexes (93%). The GLS map had 71% ≥ 0.5. A separate stability statistic, best-match Jaccard against the stage-1 runs, was higher for zavg in every size bin, but it has no chance reference (a rewired graph gives 0.84–0.90 for 3–4-gene clusters) and is descriptive only.

*Bridge sensitivities.* Without coverage matching the difference was −0.020 (−0.054 to +0.017). A reviewer probe that dropped random zavg complexes to matched coverage gave a similar point estimate (−0.020 over 50 draws), so matching by stability does not drive the result. After matching, 32% of the zavg complexes kept were locus-flagged. Pair-level AUROC against BioPlex over S × S was 0.606 for zavg_BEF, 0.579 for GLS and 0.595 for Pearson (13.9 million pairs, prevalence 0.0026; no CI).

## Supplementary Note 12. Per-direction results

### H3a and H3b on A

In H3a, 1,666 tests passed the screen and were re-tested by permutation; Biliary Tract and Prostate had no discovery. In H3b, 190 tests were re-tested. Discovery counts per lineage track the number of lines (Spearman 0.70). Registered-control readouts among the H3b discoveries are TP53 alteration with the MDM2/MDM4/PPM1D unit, MSI-H with the WRN unit, MTAP/9p21 loss with the PRMT5/WDR77 unit (a locus cluster) and NRF2 activation with the NFE2L2-containing unit (Supplementary Table S3, "control readout").

*Matched sets.* BB applied to the random-set tests selected 21 H3a families and made 11,861 rejections from 80,500 tests (14.7%), against 1,145 of 8,441 (13.6%) for the complexes. Seventeen complexes had no matched sets and 12 had both tolerances doubled (H3a). Of the 1,145 H3a discoveries, 62 have no matched sets, and of the 161 H3b discoveries 7. A post hoc, descriptive filter requested after these results were seen, β outside the central 95% of the matched sets' β, gives 699 of 1,083 (H3a) and 124 of 154 (H3b); with 10 sets per complex the quantiles are coarse, and the filter neither defines nor removes a discovery. Specificity by lineage: Lymphoid 42 of 197 discoveries beat all matched sets (21%), against 22 of 30 for Esophagus/Stomach and 28 of 48 for Lung.

*Flagged A-side examples (first-version "beats all matched sets" flag, which Supplementary Note 17 shows is largely produced by selection; not the registered strongest list).* Ranked by |t|, the strongest complex-specific lineage associations on A are lymphoid lines' lower dependency on a 22-member unit including ATRN, BCAR1, CRK, CRKL, ILK, ITGAV, ITGB5 and FERMT2 (β = +0.182) and their higher dependency on an 18-member unit including ADSL, ADSS2, ATIC, CAD, CTPS1, DHODH, GART, GMPS and IMPDH2 (−0.291). Others include lymphoid ALAD–UROS (−0.186), lymphoid EEFSEC/GPX4/PSTK/SECISBP2/SEPHS2/SEPSECS/TXNRD1 (−0.207), bowel CTNNB1/TCF7L2/HNF1A/HNF4A/GATA6/FOXA2 (−0.262) and CNS/Brain SLC3A2/SLC7A1/SLC7A5/TRPM7 (+0.176). The melanoma BRAF/MAP2K1/MAPK1/DUSP4/MITF association is the readout of control PC17. Among driver contexts: TP53-altered lines with ATM/TP53BP1/SHLD1/SHLD2/CCDC6/FBXO42 (−0.100), RB1-LoF lines with CCNE1/CCNE2/CDK2/E2F1/E2F3/TFDP1/SKP2 (−0.201), PTEN-LoF lines with IGF1R/INSR/IRS1/IRS2/FURIN (+0.096) and EWSR1-fusion lines with the EWSR1/FLI1-containing unit (−0.274; rare-context route). The bowel unit fails the registered leave-one-primary-disease-out check (Supplementary Note 16) and is therefore not among the strongest findings.

*Never-expressed sets and response scope.* The H3a never-expressed excess (13.5%; 13 BB rejections) is concentrated in melanoma (33%), Skin (28%) and Lymphoid (27%); in H3b the rate was 6.2% (1 rejection). The registered response ("if that fraction exceeds 0.075") refers to the pseudo-context fraction, which was 0.049, so it did not fire. Under a broader reading it would have fired, but its consequence (permutation re-testing of every real test with p < 0.01, then BB again) can add discoveries and cannot remove any. The literal reading was fixed before the calibrations ran but after the discovery outcomes were known, and confirmed after the calibration results were seen. The copy-number segment null gave 6.0% overall (9p21 loss 4.5%; MYCN amplification 6.1%), within the 0.075 threshold; no qualifying segment existed for MYC amplification.

*Scale effect.* Lymphoid ρ = −0.240 (p = 3 × 10⁻⁶); Myeloid −0.135; Head and Neck +0.20, Skin +0.21 and Bone +0.19 show the opposite sign. Removing the two most supportive lines lifts nine discoveries, all lymphoid, above p = 0.05.

### H1 and H1S

*Matched nulls and contrasts.* ARID1A → ARID1B: q = 2.5 × 10⁻⁶; matched null −0.072 to +0.076; paralog-versus-rest contrast −0.391 (p = 1.6 × 10⁻⁸); κ_N = −0.034 (90% CI −0.072 to +0.004). EP300 → CREBBP: q = 2.5 × 10⁻⁶; null −0.097 to +0.098; contrast −1.044 (p = 1.9 × 10⁻⁴); analytic p 3.0 × 10⁻⁵; the complex-wide axis is not evaluable because only SMARCD1 has a baseline direction. FOCAD: q = 2.3 × 10⁻⁴; null −0.062 to +0.040; 5 of 8 counted members moved in the aligned direction. Removing the three most supportive lines left ARID1A, EP300 and FOCAD at p ≤ 7 × 10⁻⁴.

*TP53 and PTEN.* TP53 loss (605 lines): κ_N = +0.056 (90% CI +0.038 to +0.075; q = 3.1 × 10⁻⁵), resting on four members near the baseline threshold (CCDC6, PABIR1, SHLD1, SHLD2), while TP53BP1 and ATM are excluded by that threshold. PTEN loss (94 lines): κ_N = −0.057 (−0.092 to −0.021; q = 0.044). At the registered m = 28 its q would be 0.051; removing one to three supportive lines raised its p to 0.017–0.038, and its counted members are only SZT2, TSC2 and ZFP36L2.

*Other labels.* Inconclusive: MED12, PLAA, SMARCA4, KEAP1, VHL, BRCA1 and PBRM1 (matched null unavailable). Not classed: EP300 and SMAD4 (complex-wide axis) and MAP2K4 (11 complete-case loss lines; not evaluable under ruling R2). SMARCB1 copy-number sensitivity: κ_N = +0.025 (90% CI −0.038 to +0.089; 15 loss lines among 732 complete cases), against +0.018 (−0.043 to +0.080; 17 of 960) in the primary analysis.

*H1S (exploratory).* Amendment A1 moved 318 silencing events to H1S; A2 then moved one of them, SMARCB1, to H1 as a genetic loss, leaving 317 events (registered m 338). They gave 321 tests (317 complex-wide axes and 4 paralog axes); m = 293 were evaluable (the remainder have a leverage-1 design); 31 BH discoveries, 26 outside locus clusters. Labels: 20 collapse, 4 detectable-within-margin (E2), 191 neutral, 45 inconclusive, 13 with no matched null and 44 not evaluable; two paralog-relief calls (SOX10, BTK). Collapse calls are dominated by lineage-restricted non-expression of B-cell, melanocyte and renal genes, which is the confounding A1 anticipated. In the lockbox, 9 of 17 evaluable calls were confirmed and all six rare-route calls passed; none is promoted.

### H2

Over 210 eligible (context, parent) families, 134 tests were re-tested by permutation; the 109 discoveries (97 outside locus clusters) involve 59 sub-modules of 40 parent complexes, and 31 are lymphoid. The H3 pseudo-contexts passed through the same machinery gave 6 stage-wise rejections against about 0 expected. 66 discoveries rank first or second against their splits (median rank 2). The post hoc β filter flags 78 of 109, and the same filter on t flags 51. The strongest rank-1 candidates were a lymphoid IPMK/IPPK/ITPK1/NAMPT/NMNAT1/PPCDC/PPCS/RFK sub-module (A β −0.290; L −0.202, confirmed), the GAB1/GRB2/PTPN11/SHC1 sub-module in lymphoid lines (L +0.149, confirmed), a lymphoid proteasome sub-module (L −0.131, confirmed) and the Head and Neck ATM/SHLD1/SHLD2/TP53BP1 sub-module (L −0.056, −0.138 to +0.026; not confirmed).

### Lockbox, per direction

*H3a.* Confirmed per context (Table S6): Lymphoid 137 of 197, Head and Neck 72 of 126, CNS/Brain 64 of 107, melanoma 45 of 91, Esophagus/Stomach 6 of 30.

*H3b.* MSI-H 12 of 14, TP53 alteration 39 of 92, WGD 0 of 26 (all 26 NOT EVALUABLE in the primary read). The eleven EWSR1-fusion discoveries, with fewer than five fusion lines in L, went by the rare-context route: five passed and six could not be evaluated because DEMETER2 lacks the genes or lines. Split by the specificity flag, the observed H3b rates were 0.583 against 0.341 among evaluable fits and 0.511 against 0.360 in the corrected relook, against an expected 0.50–0.51 and 0.29–0.30; 30 of these calls are defect-sensitive.

*BRIDGE-L multiverse.* Unweighted: −0.054 (−0.094 to −0.017). As executed: −0.072 (−0.111 to −0.033). Corrected relook: −0.069 (−0.108 to −0.029). The zavg rate was lower in every size bin, with the smallest gap for complexes of 17–60 genes.

*H1 per event.* ARID1A P: β_P −0.294 (95% HC3 CI −0.593 to +0.005; permutation p 0.0062; 12 L lines); also control PC6 in L (gene-level ARID1B −0.293, p 0.0063). EP300 P: L estimate −0.374 (−1.820 to +1.072; 4 L lines); the rare route passed (leave-one-lineage-out same sign in 30 folds, leave-one-line-out same sign over 12 A loss lines, DEMETER2 β −0.204 on 7 lines). FOCAD N: NOT EVALUABLE (defect); the corrected relook would confirm it (+0.129, +0.046 to +0.213), but it carries no claim. TP53 N: +0.104 (+0.064 to +0.144; 195 L lines), not promotable under E2. PTEN N: −0.046 (−0.135 to +0.043; p 0.15; 31 L lines), not confirmed. Neutral events in L: point estimates −0.129 to +0.006.

*Controls.* 29 of 31 evaluable control units confirmed (29 of 44 all-in; 39 of 42 corrected). The two primary failures were PC14's CDK6 readout (p_L = 0.062) and the PC5 complex readout (p_L = 0.053). PC9, PC10, PC12 and PC13 failed only because of the defect and would confirm in the corrected relook; PC18 (KIF18A with WGD) failed and would not. PC21 is not evaluable for a genuine reason. Among rare-route controls, the PC15 unit and PC16's MEN1 passed, PC16's MLLT1 failed and BRD9 could not be evaluated in DEMETER2.

*BENCH.* 1,295 pairs (1,145 H3a and 150 H3b) in 353 complexes with 920 distinct g*. Mid tercile Δ −0.100 (−0.226 to +0.155); overall −0.066 (−0.155 to +0.031); 60% of high-tercile pairs have Δ > 0; Δ rises with breadth (Spearman ρ = 0.20; no p-value, because pairs share complexes and lines). In the bootstrap, most undefined fits arise because resampled small contexts leave a line with leverage 1 (13,010 cells) or fewer than two context lines (179); the code drops them silently, so 1,032 of the 2,000 replicates average over a reduced set of pairs. The registered multiverse, reported as point estimates (A4), points the same way: eigengene −0.125 overall (+0.220 high), g* by |β| −0.076 (+0.273), T1/T2 pairs −0.084 (+0.221; 949 pairs).

*Independent review.* The reviewer recomputed 54 calls (20 H3a, 17 H3b, the 5 eligible H1 calls and 12 control units) with separate code under both designs, re-hashed all 30 persisted reads, and returned PASS with four required reporting fixes, all applied (D14). The relook's intent note named ARID1A P as an affected H1 fit; the affected paralog-axis fit is EP300 P (erratum).

## Supplementary Note 13. Sensitivity analyses

All three sensitivities were declared before they were computed; none writes a trial line or changes a discovery.

*Screen quality (H3a, H3b).* ScreenNNMD and CasActivity were added as covariates with complete-case analysis, which drops 119 A lines (all 93 KY-only screens among them). 843 of 1,145 H3a discoveries were kept (73.6%; 96 new) and 115 of 161 H3b discoveries (71.4%; 4 new); the pseudo-context rate did not change (0.049 and 0.044). Lymphoid was robust (188 of 197) despite its low specificity; Esophagus/Stomach (2 of 30), WGD (12 of 26) and MSI-H (7 of 14) were fragile. Complex-specific discoveries survived more often: 390 of 502 against 453 of 643 in H3a, and 77 of 103 against 38 of 58 in H3b. The sensitivity changes covariates and sample at once, so a lost discovery may reflect adjustment, lost power or the changed sample.

*Permutation tail.* With 10,000 permutations, the registered generalised Pareto tail extrapolates p-values below 10⁻⁴, sometimes far below the analytic p (for example 10⁻²⁴ against 10⁻¹⁰). Using p = max(p_perm, analytic p) kept 1,113 of 1,145 H3a, 130 of 161 H3b, 5 of 5 H1, 29 of 31 H1S and 107 of 109 H2 discoveries; no discovery was gained. The H3b losses concentrate in WGD (16 of 26), and three single-discovery H3b families (NRAS hotspot, female sex, SMARCA4 LoF) are selected only through the extrapolated tail. An independent reviewer obtained the same counts from separate code.

*Locus clusters.* 203 of the 1,145 H3a discoveries (17.7%) and 40 of the 161 H3b discoveries (24.8%) are in locus-flagged complexes, which make up 28.3% of complexes; locus clusters are excluded from all claims. In MSI-H, 7 of 14 discoveries are locus-flagged.

*SMARCB1 copy-number covariates (H1).* See Supplementary Note 12; the event stayed neutral.

## Supplementary Note 14. Deviations from pre-registration (full table)

Every amendment (A), erratum (E), ruling (R) and decision note (D5–D19) after the design freeze, with its timing relative to access to outcomes. Outcome access points (UTC): test labels 28 Oct 18:18; first A-side association outcome 29 Oct 00:20 (H3a run); H1 outcome 00:35; registered calibration outcomes from 01:03; first lockbox (L) value 06:50:05; last L read 09:23:05. Part 1 disclosures D0–D4 are in Supplementary Note 1.

| ID | What changed | When (UTC) | Relative to outcome access | Effect on reported results |
|---|---|---|---|---|
| D5 | FZ-C accepted with notes; disclosed that the context builder read proteomics for the LoF check before the map commit; variability required; bridge null reported with the pre-written paragraph | 28 Oct, after 22:01 | After the BioPlex bridge outcome; before any association outcome | None on construction; "zavg more stable" kept descriptive |
| D6 | H3a/H3b error control claimed as BB selective FDR, not pooled FDP; GPD/permutation implementation fixed to the registered rule | 28 Oct, ~21:50 | Synthetic data only; before any real association | Scope of the error-rate claim |
| D7 | H1 amendments decided from event counts and labels; SMARCB1 CN sensitivity pre-declared | 28 Oct, 22:54–23:28 | Before any H1 outcome | See A1, A2 |
| A1 | 318 silencing events moved from H1 to exploratory H1S | 28 Oct 22:54:50 | Before any H1 run; labels only | H1S never headlined |
| A2 | Genetic-loss rule (≥ 50% LoF or hom-del); SMARCB1 moved to H1 | 28 Oct 23:28:31 | Before any H1 run; labels only | SMARCB1 in H1 (neutral) |
| R1 / E1 | Sign typo in the registered H1 aligned score corrected to the stated semantics | 29 Oct 00:41:42 | **After the H1 run**; the engine had used the stated semantics since 28 Oct 20:09 | None (the run already used it) |
| R2 | Engine crash on leverage-1 H1S designs fixed; MAP2K4 labelled not evaluable | 29 Oct 00:47 | After the H1 run, before the H1S run | H1 not re-run; synthetic outputs identical by hash |
| R3 | GPD-extrapolated p kept primary; p_cons = max(p_perm, p_analytic) pre-declared as sensitivity; ranking by \|t\| | 29 Oct 00:59 | After the discovery outcomes, before the sensitivity was computed | Sensitivity reported (Note 13) |
| R4 | Registered calibrations (matched sets, never-expressed sets, segment null, H2 splits), absent from the engine, implemented and run once | 29 Oct 01:02–01:06 | **After all discovery runs**, before any L read | Main reading of the context results; responses applied mechanically; none fired |
| R5 / E2 | "Detectable, within equivalence margin" label when a non-neutral class and TOST-neutral both hold | 29 Oct 01:07 | **After the H1 run and its review**; post hoc, conservative | TP53, PTEN (and four H1S events) not headline, not promotable |
| D8 | Calibration response scope (literal reading), screen-quality sensitivity pre-declared, "beats all matched sets" reporting flag | 29 Oct 01:35 | After discovery and calibration outcomes; before any L read | Response did not fire; flag defines "complex-specific" |
| D9 | Lockbox dry check stopped; missing registered machinery implemented as registered | 29 Oct 01:45–01:54 | Before any L read | Enables BENCH, BRIDGE-L, rare route |
| A3 | Lockbox permutation seeds registered; BENCH within-line-rank declared not run | 29 Oct 02:03:27 | Before any L read | Seeds fixed; one multiverse score missing |
| A4 | BENCH secondary components not run; multiverse point estimates only; H1S rare-route N axes not evaluable | 29 Oct 02:42:37 | Before any L read | Partly reversed by A5 |
| A5 | A4 erratum; distinct-complex count reported; H1S rare route on A alignment; BRIDGE-L multiverse with CIs; literal readings 8–11 | 29 Oct 03:25:03 | Before any L read | As listed; A–L random-set correlation not run |
| D10 | Crash policy (persist and hash each read; resume only); dress rehearsal before the read | 29 Oct 03:18 | Before any L read | Rehearsal showed the defect's failures, which went unchecked |
| D11 | H1 framed as recovery of controls and known biology; novelty to rest on lockbox-confirmed, complex-specific, non-locus H3a/H3b hits; member profiles to be exported | 29 Oct 06:47:56 | After the intent (06:47:06), before the first L read | Framing of H1 and the Discussion |
| (export) | Member-resolved H1 profiles exported from A | 29 Oct 06:53:55 | During the read window, after H1's L read; reads A only | Fig. 4c; no call changed |
| D12 / R6 | Defect found; failed fits NOT EVALUABLE; corrected relook authorised (non-confirmatory, resume path only) | 29 Oct 10:11 / 10:19:44 | **After the read** | Evaluable and all-in rates reported |
| D13 / R7 | Primary = read as executed; corrected relook of 560 fits alongside; 53 defect-sensitive calls carry no claim | 29 Oct 10:40 / 10:41:49 | **After the read and the failed R6 check** | FOCAD, WGD and four controls lose lockbox evidence |
| D14 | Results review PASS; reporting rules (H3b three rates; BRIDGE-L R7 headline; BENCH on the top-tercile CI with disclosure, no iid p; power aid; beats-all gap not evidence); intent erratum (EP300 P, not ARID1A P) | 29 Oct 11:38 | **After the read** | Reporting only; numbers unchanged |
| D15 | Registered per-hit robustness (leave-one-out, leave-top-2-out) computed on A, against an intent and candidate list committed first (fd03b03) | 29 Oct, after the read | **After the read**; A only, no L value opened | The strongest-findings list (254) replaces the earlier provisional set |
| D16 | Strongest findings fixed by the registered rule (readings: H1 leave-one-out unit is lineage; ranking is tier then \|t_A\|); post hoc descriptive annotations (matched sets, control recovery, partial testability, scale-effect concordance); reporting rule: lead with T1, non-control, residual-robust hits; label scale-concordant lymphoid/myeloid hits as a possible lineage-scale artefact | 29 Oct, after D15 | **After the read**; post hoc readings and labels | Supplementary Table S9; labels select and rank nothing |
| D17 | Submission version; post hoc reporting labels (lead example; possible lineage-scale artefact) | 29 Oct, after D16 | **After the read** | Superseded by the revised labels (D19) |
| D18 | Peer review round 1 (major revision): response plan | 29 Oct | **After the read** | Interpretation and reporting |
| D19 (R-1 to R-4) | Revision analyses: non-member matched sets and selection-aware flag reference; spline trend with one artefact rule in all lineages; DEMETER2 on A; BRIDGE-L power for both maps; intent 4157f1e committed before computation | 29 Oct | **After the read and peer review**; A only; L sample sizes and primary call labels only | Narrower matched-set reading; flag no longer defines "complex-specific"; revised labels in Table S9; BENCH called inconclusive |
| (not run) | Alternative scores (eigengene, within-line rank, expression-masked, pan-essential-referenced) per hit for H3a/H3b; per-hit robustness for H1S and H2 | — | Never run | H1S and H2 are exploratory and outside the strongest-findings rule |

## Supplementary Note 15. Extended limitations

*Part 1 (test evaluation)*
1. **The single-look guarantee is procedural.** The test labels were plain files, and unsupervised genome-wide scores covering the test pairs already existed. Incident D1 shows that the risk was real; one operating point is non-blind as a result.
2. **The gain is carried by pair-rich clusters.** Macro per-cluster AP does not differ from GLS's; in the descriptive size breakdown, made after the look, zavg_BEF is not better than GLS for clusters with ≤ 5 positive pairs. The test has 238 clusters with a Kish effective number of 21.3, so CIs are wide for Spearman and Pearson.
3. **Training-label re-measurement.** zavg_BEF was selected on hu.MAP training labels derived from Complex Portal [Fischer 2025], so agreement with physical complexes partly re-measures that selection.
4. **Calibration.** Calibration-in-the-large of the train-fitted isotonic calibrator was 0.893 on the test set, with no CI pre-specified; no recalibration was made, and the calibrators' own sampling error is not propagated.

*Complex map and bridges*
5. **Neither bridge supports the map.** The primary (BioPlex) bridge is null, with an interval that allows an advantage of up to 0.031 for zavg or 0.075 for GLS; protein co-abundance favours GLS; the secondary bridge (BRIDGE-L) is negative, and power accounts for roughly a quarter to two fifths of its difference (Supplementary Note 17). The locus-cluster explanation of the BioPlex deficit is post hoc and does not explain BRIDGE-L.
6. **Selection is flat across z0** at k = 10, so a different z0 could have been selected on the same evidence.
7. **Locus clusters are frequent** (28% at level 1); they are flagged, not removed, and T2 is not locus-aware.
8. **Stability has no chance baseline.** The variability analysis shares subsamples with the frozen map, and its random reference is a floor.
9. **Some complexes merge two complexes or are pathway modules**, so "complex" means physically supported module.
10. **Coverage.** About 40% of S genes are in no level-1 complex, including SMARCA2, STAG1 and KIF18A, which limits paralog-based analyses.

*Context associations (H3a, H3b, H2)*
11. **Context associations are shared by coherent co-dependency modules.** The registered comparators are fragments of other complexes; coherent non-member modules could be matched for only 90 of 367 (the less coherent) units; essentiality-only sets are associated less often. The "beats all matched sets" flag rests on 10 sets per complex, is met by a third of random-atlas discoveries after selection, and no longer defines complex-specificity; 69 discoveries have no matched sets, and the 95% filter is post hoc.
11b. **Selective FDR.** BB controls the average within-family error over selected families, not the pooled error (0.068 in a simulation that planted lineage effects only). Its guarantee assumes independence across families and PRDS within them; lineage contexts share reference lines and H3b contexts overlap, and no simulation covered the H3b dependence. Permutation p-values for screened tests re-enter the procedure (a hybrid p). With 21 of 23 H3a families selected, BB is close to within-family BH at q ≈ 0.046.
12. **Calibrations implemented after the discovery runs** (R4). The literal reading of the response scope was fixed before the calibrations ran and confirmed after their results were seen; a broader reading would have added re-tests but removed none.
13. **GPD tail.** Extrapolated permutation p-values can lie many orders of magnitude below the analytic p; using the larger drops 32 H3a, 31 H3b, 2 H1S and 2 H2 discoveries.
14. **Screen quality and scale.** About a quarter of H3a/H3b discoveries are lost with screen-quality covariates (with 119 lines dropped). Esophagus/Stomach, WGD and MSI-H are fragile. The lymphoid association between complex essentiality and β (ρ −0.24) persists after adjustment and suggests a Chronos dynamic-range effect. Never-expressed sets show a residual lineage artefact (13.5% at p < 0.05); we did not rescale statistics by a per-context inflation factor and instead mark leads in melanoma and Lymphoid as less reliable. That flag does not certify the other lineages: among the other lead contexts the rate is 0.01 (Myeloid) to 0.22, about four times the nominal 0.05 for Ovary/Fallopian Tube, Bone and Head and Neck.
15. **Per-hit robustness was computed after the lockbox read.** The registered leave-one-out and leave-top-2-out checks were deferred at the discovery runs and run after the read, on A only, against a committed candidate list. Alternative scores (eigengene, within-line rank, expression-masked, pan-essential-referenced) remain not run. The leave-one-primary-disease-out rule is strict where one disease dominates a lineage (Bowel kept 0 of 29 candidates) and partly untestable where the remainder has fewer than two lines (21 H3a and 33 H3b findings).
16. **The scale-effect labels are descriptive and incomplete.** The revised rule (spline trend, all lineages) cannot flag findings beyond the trend in its own direction (26 of 37 leads), and four leads lie within 0.05 of the unregistered |z| = 2 threshold. The labels were fixed after the robustness results were known. They select the examples shown in the text, so they are a post hoc selection for presentation, not purely descriptive annotations; BioPlex tiers (validation data) also steer which findings are highlighted (T1 leads). All four leads named in the main text lie beyond the trend in its own direction, so for them the DEMETER2 check, not the trend rule, is the main independent support.
17. **H2 is exploratory and mostly reproduced by random splits** (78.3 vs 109 stage-wise rejections). Lockbox confirmation (0.70) shows that the divergence recurs, not that it is sub-module-specific; the rank-1 subset was not confirmed more often.

*H1*
18. **The H1 headline events re-find known biology.** ARID1A → ARID1B is control PC6; FOCAD → PELO–HBS1L overlaps control PC10. EP300 → CREBBP rests on 12 A loss lines at the evaluability floor, a module without physical support (T3-neg) and the rare-context route (4 L loss lines). The FOCAD lockbox call is defect-sensitive.
19. **Neutral H1 events rest on A, under a loose margin.** L has a median of 8 loss lines per neutral event, and only 2 of 11 L 90% CIs lie inside ±0.1. The ±0.1 margin is larger than many H3a effects called real (325 of 1,145 have |β| < 0.05). At a ±0.05 margin none of the 11 would be neutral: no committed 90% CI lies inside it (post hoc; Supplementary Note 17, R-5), so the neutral call depends entirely on the registered ±0.1 margin. The registered margin in residual-SD units was not reported, because the per-fit residual SD was not exported. One of the 11 (ARID1A) is neutral only on the complex-wide axis.
20. **Thin control coverage.** H3a rests on two controls, the H1 pass on one multi-gene readout, and H2 has none; five controls were not evaluable at complex level; in L, four controls are defect-sensitive.
21. **Amendments and post hoc rules.** A1 and A2 changed the H1 panel before any H1 run, from labels only. E1 corrected a sign typo after the run, matching the code that ran. E2 and the 95% specificity filter are post hoc. The PTEN call depends on the evaluability-adjusted m (q = 0.044 at m = 24, 0.051 at the registered 28). The member-resolved profiles were exported after the discovery run, from A only, during the lockbox read window.

*Lockbox*
22. **A defect in the confirmation driver**, found only after the read although the dress rehearsal had shown its failures. The handling rules (R6, R7) were made after the outcomes were seen; they keep the read as executed as primary and strip claims from all 53 defect-sensitive calls, at the cost of H3b's WGD context and four controls. The H3b rates depend on how failures are counted (0.40–0.49).
23. **Lockbox size and power.** L has 338 lines, about a third of A. The rare-context route re-uses A lines and DEMETER2 and is not an independent replication. Non-confirmation of small effects is often uninformative (325 of 1,145 H3a discoveries have |β| < 0.05).
24. **The power aid is post hoc.** It was computed after the read by the reviewer, is not registered, pools an empirical-Bayes prior across contexts, treats correlated tests as independent and uses normal approximations.
25. **BENCH is inconclusive.** 16.5% of bootstrap pair × replicate fits fail and are dropped; the high-tercile CI depends on how the covariate design is assembled within a replicate (lower bound −0.033 or −0.122), so the registered interval is unreliable. The persisted, hash-checked L reads would allow a re-run with an explicit failure policy; we chose not to make this second, post hoc look at L, to keep the single-read design. Within-line-rank scoring and the A–L random-set correlation were not run, and the multiverse is point estimates only.
26. **The loader count is inferred** from ledger invariants and file times, because the run log of the lockbox loader was not committed; it is permanently lost.
26b. **Lockbox independence.** L lines were part of the data on which zavg_BEF was developed in Part 1, and their Chronos values were inferred jointly with A in the same release. Confirmation establishes reproducibility within DepMap 26Q3 processing, not independence from shared artefacts. The power aid's prior includes gene-set-generic and scale effects, which would also replicate, so agreement with power cannot separate biology from reproducible artefacts.

*General*
27. **Single cell-line resource.** All dependencies come from one CRISPR resource (DepMap 26Q3, Chronos); DEMETER2 enters only the rare-context route. No result has been tested in another screen collection, under different culture conditions or in vivo.
28. **Agent-run study.** Analyses, reviews and rulings were carried out by AI agents under a human principal investigator. Independent agent reviews reproduced every confirmatory number, but reviewers and producers share a model family.

## Supplementary Note 16. The registered strongest findings

*Rule and computation.* The rule is stated in DESIGN v3.1 §2.5 and the success criterion of H3a, H3b and H1; it is implemented in `src/part2/strongest.py` and run by `scripts/56_strongest_findings.py`. The intent note and candidate list were committed (fd03b03) before any robustness value was computed, and the checks ran after the lockbox read on A only (118 s). A full A refit reproduced every committed discovery β and t for all 575 hits. The PI accepted two readings of ambiguous registry text after the run: the H1 leave-one-out unit is the lineage, and "A statistic × tier" means ranking by tier class, then |t_A|.

| Direction | Candidates | Strongest (leave-one-out robust) | Units | T1 / T2 / T3-neg |
|---|---|---|---|---|
| H3a | 524 | 207 | 137 | 124 / 56 / 27 |
| H3b | 46 (41 lockbox, 5 rare route) | 45 | 41 | 27 / 9 / 9 |
| H1 | 2 | 2 (ARID1A → ARID1B P; EP300 → CREBBP P, rare route) | 2 | 1 / – / 1 |

*Leave-one-out behaviour.* Of the 317 H3a failures, 308 fail exactly one of the 81–82 primary-disease folds. Failures cluster in lineages dominated by one primary disease: Bowel 0 of 29 robust (colorectal adenocarcinoma is 60 of 64 Bowel A lines), Head and Neck 1 of 54, Skin 1 of 22 and Peripheral Nervous System 5 of 44. Where fewer than 2 context lines remain, a fold is skipped, as the engine rule requires, which is why Pancreas passes 17 of 17. The registered rule therefore tests "not carried by one primary disease"; it is strict for lineages with a small non-dominant remainder and vacuous when that remainder has fewer than 2 lines. Lymphoid (92 of 119) and Myeloid (21 of 37) have several sizeable diseases and pass more often. The one H3b failure is MSI-H Z0334 (1 of 30 lineage folds).

*Relation to the earlier provisional set.* The draft's provisional "strongest specific" set (261 H3a and 29 H3b) used two unregistered filters (beats all matched sets; not a control readout). Of it, 107 H3a and 28 H3b findings are in the registered list. Those filters are now descriptive annotation columns.

*Annotations (post hoc, descriptive).* Over the 254: tier T1 152, T2 65, T3-neg 37; beats all matched sets H3a 107 of 207, H3b 34 of 45; control readouts 7 (H1 ARID1A P, PC6; H3a Lymphoid and Myeloid Z0113, PC20; H3b MSI-H Z0286, PC10; TP53-altered Z0166, PC2; NRAS Z0132, PC4/PC5; MSI-H Z0157, PC1); same unit as a control in another context 18; leave-one-out partly untestable H3a 21 (Pancreas 17, Ovary 3, Pleura 1) and H3b 33 (TP53-altered 28, EWSR1 fusion 5).

*Scale-effect annotation.* For each H3a context, β was regressed on unit essentiality over all 367 tested units (slope, intercept, Spearman ρ and residual SD in `results/part2/robustness/scale_effect_by_context.tsv`). Concordant means the hit's β has the sign of its fitted value; where the slope is near 0, the fitted sign is set mostly by the context's intercept, so "concordant" there means "same sign as the lineage-wide shift". Lymphoid ρ is −0.24 and Myeloid −0.135 (slope −0.029); the other lineages have |ρ| ≤ 0.21, mostly positive. Overall, 142 of 207 H3a findings are concordant, 79 have |residual z| ≥ 2 and 55 are both.

*Testability.* The leave-one-out check was fully applicable to 199 findings (H3a 186, H3b 12, H1 ARID1A → ARID1B), partly applicable to 49 lockbox-route findings (H3a 21, H3b 28 TP53-altered) and replaced by the rare-context route for 6 (five EWSR1-fusion H3b findings and EP300 → CREBBP). Under the alternative readings of the two ambiguous registry items (post hoc, committed tables only; Supplementary Note 17, R-5), the H1 leave-one-out unit can affect only the two H1 calls (H3a 207 and H3b 45 are unchanged; 252–254 in all; whether either call changes cannot be determined without a refit), and ranking by |t_A| alone or by |t_A| × tier weight changes only the order (Spearman correlation with the registered rank 0.37 and 0.89).

*Revised reporting labels (Table S9, `rev_*` columns; post hoc, Supplementary Note 17).* The first-version labels below are kept in Table S9 for traceability; the revised manuscript uses the revised labels: 95 possible scale artefacts in any lineage, 37 leads (26 beyond the trend in its direction; 9 in less reliable lineages; 29 of 35 evaluable concordant in DEMETER2).

*First-version reporting labels (Table S9; post hoc, D16).* 37 H3a findings are lead examples (T1, no control annotation, |residual z| ≥ 2); 10 of them are lymphoid or myeloid. 56 findings (52 lymphoid, 4 myeloid) are labelled a possible lineage-scale artefact (concordant, |residual z| < 2); 38 of them have positive β, and their median unit essentiality is −0.63, against −0.41 for all H3a findings and −0.30 over all 367 units. Of the 45 H3b findings, 25 are T1 and not control readouts.

## Supplementary Note 17. Revision analyses (referee rounds 1 and 2; post hoc, descriptive)

Four analyses were run after peer review, on the construction set A only, with an intent committed before any value was computed (`results/part2/revision/INTENT.md`, 4157f1e; results 740f439). None is registered, none is a test with error control, and none changes a registry entry, discovery, confirmation or ledger. No lockbox value was read; R-4 uses only committed L sample sizes and the committed primary call labels. An independent agent checked them and returned PASS-WITH-NOTES; its notes are reflected below. Scope of that check: it was not blind (the agent read the results README before computing); it reran scripts 59–61 (byte-identical apart from runtime fields) and recomputed the numbers of all four analyses from the committed tables with its own code; it did not rerun script 58, i.e. the ten random atlases with their permutation p-values and the non-member test table, and it reused, rather than re-implemented, the empirical-Bayes (Tweedie) step of the R-4 power aid. Every number it could recompute matched.

*R-1(a). Matched sets drawn from genes outside every complex.* The pool is the 2,380 S genes in no unit of the map. Three generators were declared in advance: coherent neighbourhoods matched on size, essentiality and coherence (N-nbr); uniform draws with the same matching (N-unif); and uniform draws matched on size and essentiality only (N-ess). Uniform draws never reached complex-level coherence (0 units). N-nbr could be matched for 90 of 367 units in H3a (87 in H3b), and these are the less coherent units (median target coherence 8.9 z against 12.3 z for the others). N-ess could be matched for 355 (H3a) and 354 (H3b).

| Direction | Units compared | Complexes | Sets from other complexes | Coherent non-member sets (N-nbr) | Essentiality-only sets (N-ess) |
|---|---|---|---|---|---|
| H3a | all 367 | 0.281 | 0.281 | – | – |
| H3a | the 90 with N-nbr sets | 0.230 | 0.248 | 0.251 | 0.164 |
| H3a | the 355 with N-ess sets | 0.281 | 0.282 | – | 0.166 |
| H3b | all 367 | 0.120 | 0.116 | – | – |
| H3b | the 87 with N-nbr sets | 0.113 | 0.112 | 0.108 | 0.084 |
| H3b | the 354 with N-ess sets | 0.119 | 0.116 | – | 0.088 |

Selective FDR over the non-member tests made 1,666 of 20,700 (N-nbr) and 2,459 of 81,650 (N-ess) rejections in H3a, and 44 of 17,400 and 337 of 70,800 in H3b. Over the complex-member matched sets it made 11,861 of 80,500 in H3a, against 1,145 of 8,441 for the complexes. Reading: a complex is not more context-associated than a coherent, essentiality-matched co-dependency module; sets matched on essentiality alone are associated less often, though still well above 5%. The coherent non-member comparison rests on the less coherent quarter of the map.

*R-1(b). Selection-aware reference for "beats all matched sets".* Ten random atlases, each replacing every unit by its r-th matched set, went through the full discovery procedure (screen, 10,000-permutation re-test, selective FDR; full mode). They gave 1,016 to 1,277 H3a discoveries per atlas (mean 1,109; real 1,145) and a mean of 138 in H3b (real 161).

| | H3a | H3b |
|---|---|---|
| Real discoveries beating all 10 matched sets | 502 / 1,083 = 0.464 | 103 / 154 = 0.669 |
| Random-atlas discoveries beating the 9 sibling sets and the real unit | 3,649 / 11,087 = 0.329 | 667 / 1,381 = 0.483 |
| Same flag over all random tests (no selection; 1/11 = 0.091) | 0.087 | 0.087 |
| Mean P(\|t\| real discovery > \|t\| random discovery), same context | 0.482 | 0.476 |
| Same for \|β\| | 0.493 | 0.494 |

Across the ten atlases, the random-discovery rate ranged from 0.291 to 0.364 (H3a) and from 0.448 to 0.517 (H3b); no atlas reached the real rate (R-5; a spread across atlases, not a confidence interval). Real discoveries beat their own matched sets more often than random discoveries beat theirs, but across discoveries their |t| and |β| are not larger. The excess is relative to each complex's own matched sets, not an effect-size advantage. The first version's unconditional contrast (13.3% of all H3a tests against 9.1% under exchangeability; no CI, because tests share complexes and lines) is superseded by this reference. Per-context values: `REV/r1b_per_context.tsv` (Fig. 6b).

*R-2. Lineage scale effect with a non-linear trend.* Per H3a context, a natural cubic spline (4 df) of β on unit essentiality was fitted over all 367 units. It raises R² only slightly (Lymphoid 0.068 linear, 0.120 spline; every other context ≤ 0.083). One rule was applied in every lineage: *possible scale artefact* = sign as the fitted trend predicts and |spline residual z| < 2; *lead* = T1, not a control readout, |spline residual z| ≥ 2 and not a possible scale artefact. Results: 95 of 207 H3a findings are possible scale artefacts (56 under the first version's Lymphoid/Myeloid-only rule); 37 are leads (35 kept, Bone Z0130 and Myeloid Z0278 lost, CNS/Brain Z0146 and Breast Z0253 new). 26 of the 37 lie beyond the trend in the trend's direction, which the rule cannot flag by construction, and four lie within 0.05 of |z| = 2 (Bone Z0001, Breast Z0253, Lung Z0087, Kidney Z0362), so the lead set is sensitive to the threshold. Never-expressed gene sets reach p < 0.05 in 27% of Lymphoid and 33% of melanoma tests (Skin 28%); the 9 leads in these lineages (Lymphoid 5, melanoma 4) are marked less reliable in Table S9 (`rev_lineage_less_reliable`). The flag marks only the three highest rates; the other lead contexts range from 0.01 (Myeloid) to 0.22 (Ovary/Fallopian Tube), and Bone and Head and Neck are at 0.19. Not done: excluding scale-concordant leads (referee M2(c)), because the rule as specified keeps them and we report them instead; and rescaling statistics by a per-context inflation factor from the never-expressed sets (M2(d)), which would need new computation; the never-expressed rate is reported as a reliability column instead.

*R-3. DEMETER2 on A.* The same complex score, removals and HC3 design were refitted on DEMETER2 for the 423 A models with RNAi data (no FDR).

| Hits | Evaluable | Same sign | One-sided p < 0.05 | Spearman β_A vs β_D2 |
|---|---|---|---|---|
| Revised leads | 35 of 37 | 29 | 20 | 0.86 |
| Leads with ≥ 12 context models | 28 | 24 | 18 | 0.87 |
| All 207 H3a strongest | 194 | 157 | 88 | 0.65 |

Not evaluable: Peripheral Nervous System Z0190 (no PNS model has DEMETER2 data) and Kidney Z0362 (fit failed). Discordant leads: Soft Tissue Z0080 and Z0336, Kidney Z0323, Pancreas Z0283, Lung Z0001 and Bone Z0001. Same sign but weak (one-sided p 0.10–0.14): Lymphoid Z0190, melanoma Z0190, Pleura Z0067 (6 context models) and Lung Z0150. Per-hit values are in Table S9 (`rev_D2_*`).

*R-4. BRIDGE-L power comparison.* The lockbox review's power aid was applied unchanged to the zavg-map and GLS-map hits, with an empirical-Bayes prior per map and family.

| Power version | Expected zavg | Expected GLS | Expected Δ (SD) | Observed Δ | Share of observed Δ |
|---|---|---|---|---|---|
| Naive | 0.686 | 0.703 | −0.017 (0.018) | −0.065 | 26% |
| EB posterior mean | 0.500 | 0.526 | −0.026 (0.018) | −0.065 | 41% |
| EB mean + variance (headline) | 0.459 | 0.482 | −0.023 (0.019) | −0.065 | 35% |

GLS-map discoveries have larger A effects (median |β| 0.081 against 0.066) and similar |z| and L counts. The headline remainder (−0.042, 2.2 SDs) uses an SD that treats hits as independent and ignores uncertainty in the prior, so it is approximate. The headline model also under-predicts both maps' observed levels (zavg 0.459 against 0.514; GLS 0.482 against 0.578), more so for GLS; GLS hits exceed their expected rate in every size bin (0.57–0.59 against 0.47–0.50), zavg hits less so (0.48–0.57 against 0.42–0.49). Reading: power accounts for roughly a quarter to two fifths of the BRIDGE-L difference.

*Referee requests not carried out, and why.*
- Preprocessing-matched GLS or Spearman baseline, and ΔAP leaving out each of the six largest test clusters (M5): both need the test labels again, i.e. a second look at the test set, which the protocol forbids and this revision did not open.
- BENCH bootstrap with an explicit failure policy (M4.1): the persisted L reads would allow it, but we chose not to make a second, post hoc look at L, to keep the single-read design. The registered interval is reported as unreliable and the result as inconclusive. (Round 1 said this "needs the lockbox values again", which overstated it.)
- H1 margin in residual-SD units (M6): needs each fit's residual SD, which was not exported; not computed. The ±0.05 sensitivity was computed in round 2 from the committed 90% CIs (R-5); round 1 wrongly said it was not in committed outputs.
- Counts of the strongest findings under the alternative registry readings (M7.2): bounded in round 2 from committed tables (R-5); the exact H1 outcome under the primary-disease unit would need an A refit and was not computed.
- Per-context rescaling by never-expressed inflation (M2(d)): not computed (above).

*R-5. Referee round 2 (post hoc; committed tables only).* An intent note (`results/part2/revision/INTENT_round2.md`) was committed before `scripts/62_rev_round2.py` computed three counts. No statistic was refitted, and no dataset, lockbox value or test label was opened.
- *H1 at ±0.05.* Of the 11 events labelled neutral at the registered ±0.1, none has a committed 90% CI inside ±0.05 (none on the boundary either). The smallest symmetric margin containing each CI is 0.059 (STK11) to 0.097 (KDM6A), median 0.076. The registered labels are unchanged.
- *Alternative registry readings.* H1 leave-one-out by primary disease instead of lineage: only ARID1A P and EP300 P can change; the other three eligible H1 calls fail other criteria (defect-sensitive, not confirmed, or not headline-eligible). The total is 252–254, and the per-fold fits needed to decide were not computed. Both calls passed all 30 lineage folds. Ranking by |t_A| alone changes every rank (top-10 overlap 5, top-30 overlap 17); |t_A| × tier weight (3, 2, 1; illustrative) changes 225 ranks (overlaps 9 and 28). No finding enters or leaves.
- *Atlas range.* See R-1(b) above (`REV/r5_atlas_flag_range.tsv`).

*Best gene g* of the text-highlighted findings (BENCH, committed lockbox results).* Δ = t_L(complex) − t_L(g*), signed toward the A direction (positive: the complex replicated better than its best gene). The L estimate of g* alone is not in the committed outputs.

| Context | Unit (members) | g* | Breadth | Δ |
|---|---|---|---|---|
| Ovary/Fallopian Tube | Z0035 (17) | PAX8 | 0.59 | +1.05 |
| Bone | Z0303 (5) | SNAPC5 | 0.80 | +1.29 |
| Lung | Z0174 (6) | CSNK1A1 | 0.83 | +0.17 |
| Head and Neck | Z0051 (16) | GRB2 | 0.50 | −1.23 |
| Pleura | Z0067 (14) | FGFR1 | 0.93 | +0.63 |
| Myeloid | Z0254 (12) | DNM2 | 0.83 | −1.69 |
| Lymphoid | Z0030 (18) | ADSS2 | 0.83 | +0.47 |
| RB1 loss | Z0069 (14) | SKP2 | 0.79 | +0.83 |

## Supplementary Note 18. Glossary of identifiers

| Identifier | Meaning |
|---|---|
| zavg_BEF | The locked Part 1 co-complex score (mean of three z-scored similarity components) |
| ITT, S1 | All canonical test pairs (intention-to-score); pairs mapped under every method |
| C1–C3 | Test pairs with both, one or neither protein among training-positive proteins |
| A, L, U | Construction set (1,007 lines), lockbox (338 lines), their union (1,345 lines) |
| S | The 6,022 genes with a dependency in ≥ 1% of A lines |
| Z0001–Z0367 | Level-1 complex (unit) identifiers of the frozen zavg map; Z0xxx.n are sub-modules |
| T1, T2, T3-neg, T3-nc | Tiers: BioPlex-supported; co-abundant (not T1); tested, unsupported; not testable |
| H3a, H3b, H1, H1S, H2 | Lineage atlas; driver/genome-state family; within-complex genetic interactions; exploratory silencing events; sub-module divergence |
| BENCH | Complex score against its best member in the lockbox |
| BRIDGE-L | Registered secondary bridge: lockbox confirmation of zavg-map against GLS-map discoveries |
| PC1–PC22 (PC19a/b) | Pre-registered positive controls (23 entries) |
| Rare-context route | For contexts with < 5 L lines: leave-one-lineage-out, DEMETER2 on A lines and leave-one-line-out must agree in sign; re-uses A lines, not an independent replication |
| FZ-C, FZ-X, FZ-LB | Freeze points: complex maps; contexts, families and events; lockbox split |
| A1–A5 | Registry amendments (H1 panel; lockbox seeds and readings) |
| E1, E2 | Errata: H1 sign typo; "detectable, within equivalence margin" label |
| R1–R7 | PI rulings; R4 = calibrations implemented after discovery; R6, R7 = handling of the lockbox driver defect |
| D0–D21 | Decision notes in RESEARCH_TREE.md (D1 = the non-blind z ≥ 3 operating point) |
| R-1 to R-5 | Revision analyses, rounds 1 (R-1 to R-4) and 2 (R-5) (Supplementary Note 17) |
| Defect-sensitive | A lockbox call that differs between the read as executed and the corrected relook; carries no claim |

## Supplementary Tables

| Table | Content | File |
|---|---|---|
| S1 | Pre-registered controls (A) and kill-rule tally | `manuscript/tables/tableS1_controls.md`, `tableS1b_kill_rule.tsv` |
| S2 | Largest T1 non-locus complexes with post hoc HGNC gene-group labels (for example Z0001 iron–sulfur assembly, Z0002 mitochondrial tRNA synthetases with mtDNA replication and transcription, Z0003 cytochrome c oxidase assembly, Z0004 centriole duplication, Z0008 spindle checkpoint with MIS12, Z0010 TRAPP/COPII, Z0012 complex I; labels are post hoc) | `results/part2/construction/summary/zavg_top_T1_complexes.tsv` |
| S3 | All 1,306 H3a/H3b discoveries with specificity, locus, GPD and screen-quality flags | `manuscript/tables/tableS3_context_discoveries.md` (full `.tsv`) |
| S4 | H1 panel, 24 events | `manuscript/tables/tableS4_h1_panel.md` (full `.tsv`) |
| S5 | H2 discoveries with rank against own random splits | `manuscript/tables/tableS5_h2_split_ranks.md` (full `.tsv`) |
| S6 | Lockbox confirmation rates per direction: evaluable / all-in / corrected, defect-sensitive counts, power expectation | `manuscript/tables/tableS6_lockbox_rates.md` |
| S7 | Lockbox results for all H3a/H3b discoveries; the Markdown rows show the registered strongest findings (S9) in rank order; the earlier unregistered subset is only a column of the full TSV | `manuscript/tables/tableS7_lockbox_h3.md` (full `tableS7_lockbox_h3_full.tsv`) |
| S8 | Lockbox results for H1, H2, BENCH and BRIDGE-L | `manuscript/tables/tableS8_lockbox_h1_h2_bench_bridge.md` |
| S9 | The 254 registered strongest findings with post hoc annotations, first-version labels and the revised labels (spline trend, one rule in all lineages, never-expressed reliability, DEMETER2) | `manuscript/tables/tableS9_strongest_findings.md` (all rows: `tableS9_strongest_findings.tsv`; counts: `strongest_derived.json`) |
