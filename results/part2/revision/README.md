# results/part2/revision: referee round 1 (analyses R-1 to R-4) and round 2 (R-5)

**Every analysis in this folder is POST HOC and DESCRIPTIVE.** None is registered, none is a trial, none carries error
control, and none changes a registry entry, discovery, confirmation, primary result, ledger or `trials.jsonl`
(no `run_association` call is made). The intent note `INTENT.md` was committed (4157f1e) before any value was computed.

- **Data.** Construction set A only (1,007 models; DEMETER2: the 423 A models with RNAi data). No lockbox L value,
  persisted read or L row was opened. R-4 reads only the committed L sample sizes (`n_L`, `n_context_L`, by column) and
  the committed R7 primary call labels. No hu.MAP, STRING or curated complex resource.
- **Code.** `src/part2/revision.py` (tested in `tests/test_part2_revision.py`), `scripts/58_rev_matched_sets.py` (R-1),
  `59_rev_scale_effect.py` (R-2), `60_rev_demeter2.py` (R-3), `61_rev_bridge_power.py` (R-4),
  `62_rev_round2.py` (R-5, round 2; committed tables only, intent `INTENT_round2.md`).
- **Seeds.** Fixed as in `INTENT.md`: entry calibration seeds 4600 / 4700 for set generation, entry permutation seeds
  4000 / 4100 for the random-atlas permutations, `PYTHONHASHSEED=0`.

| Script | Wall time | Peak RSS |
|---|---|---|
| 58 (R-1; 7 worker processes) | 1,797 s (set generation 242 s; random atlases H3a ≈ 16 min, H3b ≈ 7 min) | 3.36 GB main, 0.71 GB per worker |
| 59 (R-2) | 1.5 s | 0.19 GB |
| 60 (R-3) | 14 s | 2.93 GB |
| 61 (R-4) | 2.7 s | 0.25 GB |
| 62 (R-5) | 1.6 s (PI rerun, outputs byte-identical) | 0.17 GB |

Wall times are whole-process (interpreter start and imports included). The `runtime_s` fields in the summary JSONs
time the computation only (58: 1,795.6 s; 59: 0.11 s; 60: 12.6 s; 61: 1.2 s).

**Reproduction checks (all passed).**
- The regenerated complex-member matched sets reproduce every committed calibration t (H3a 80,500, H3b 70,200 tests)
  to the stored 6 significant digits (max relative difference 5.0e-6). `INTENT.md` said |Δt| < 1e-6, but the committed
  table stores t as `%.6g`, so agreement was checked to that precision. This is a disclosed deviation.
- The discovery-selection code reproduces the committed discovery flags exactly when given the committed permutation
  p: H3a 1,145, with the 1,666 screened tests; H3b 161, with 190 screened.
- 4 committed permutation p per direction were recomputed bit-identically.
- R-2's linear residual z reproduces the D16 annotation, and its labels reproduce D16's 37 leads and 56 artefacts.
- R-4's observed rates reproduce the committed R7 primary-evaluable BRIDGE-L: 0.5137 against 0.5785, Δ −0.0648.

---

## R-1 Matched-set reference (referee M1)

### (a) Matched sets drawn from genes outside every frozen complex

The pool is the S genes that are A-effect columns and appear in no unit of the zavg membership table: 2,380 genes. The
acceptance rule is the one the registered calibration uses (size; essentiality ±0.1; coherence ±1.0 z; R = 10;
tolerances doubled once). The three candidate generators are pre-declared in `INTENT.md`.

**Matching quality** (`r1a_matching_quality.tsv`):

| Generator | Units with sets (H3a / H3b) | Median \|Δess\| | Median \|Δcoh\| (z) | Note |
|---|---|---|---|---|
| N-nbr (neighbourhood draws, full matching) | 90 / 87 of 367 | 0.061 | 0.72 | Matched units have lower coherence (median target 8.9 z, against 12.3 z for unmatched units) |
| N-unif (uniform draws, full matching) | 0 / 0 | – | – | Uniform non-member sets never reach complex-level coherence |
| N-ess (uniform draws, no coherence matching) | 355 / 354 | 0.051 | 11.4 | Size and essentiality matched only |

**p < 0.05 rates** (`r1a_rates.tsv`; pooled over contexts):

| Direction | Scope | Real complexes | Complex-member sets (existing) | Non-member sets |
|---|---|---|---|---|
| H3a | all units | 0.281 (8,441) | 0.281 (80,500) | – |
| H3a | units with N-nbr sets (90) | 0.230 | 0.248 | **N-nbr 0.251** |
| H3a | units with N-ess sets (355) | 0.281 | 0.282 | **N-ess 0.166** |
| H3b | all units | 0.120 (7,340) | 0.116 (70,200) | – |
| H3b | units with N-nbr sets (87) | 0.113 | 0.112 | **N-nbr 0.108** |
| H3b | units with N-ess sets (354) | 0.119 | 0.116 | **N-ess 0.088** |

BB rejections over the non-member tests (families = contexts):
- H3a: N-nbr 1,666 of 20,700; N-ess 2,459 of 81,650.
- H3b: N-nbr 44 of 17,400; N-ess 337 of 70,800.

**Reading.**
- Coherent co-dependency neighbourhoods of genes that belong to **no** frozen complex are as context-associated as
  real complexes and as fragments of other complexes. This holds on the units where such neighbourhoods can be
  matched, which are the less coherent units.
- Essentiality-matched sets **without** coherence are less often associated: 16.6% against 28.1% (H3a) and 8.8%
  against 11.9% (H3b). That is still well above 5%.
- So the parity result is "a complex is not more context-associated than a coherent, essentiality-matched gene module".
  It is not "than any essentiality-matched gene set". Sets matched for coherence behave like complexes; sets matched
  for essentiality alone associate less often. (Round 2: an earlier sentence attributed the excess to coherence. That
  attribution is dropped, because coherence is itself measured on A as co-variation across lines, so matching on it
  partly conditions on the outcome.)
- Caveat: N-nbr could be matched for only about a quarter of the units, and these are the less coherent ones.
- Like-for-like on the same 90 / 87 units: real 0.230, N-nbr 0.251, N-ess 0.164 (H3a); real 0.113, N-nbr 0.108,
  N-ess 0.084 (H3b) (`r1a_rates.tsv`; added after the independent review, 19ad2ae).

### (b) Selection-aware reference for "beats all matched sets"

**Method.**
- Each of the 10 random atlases replaces every unit by its r-th complex-member matched set.
- Each atlas goes through the full discovery selection: BB, then a 10,000-permutation Freedman–Lane re-test of the
  q < 0.1 screen, then BB again.
- This was full mode: the cost rule did not trigger, so no approximation was used.

**Atlas discoveries.**
- H3a: 1,016 to 1,277 random discoveries per atlas (mean 1,109), against 1,145 real.
- H3b: mean 138 per atlas, against 161 real.

**Results** (`r1b_*`):

| | H3a | H3b |
|---|---|---|
| Real discoveries that beat all 10 matched sets | 502 / 1,083 = **0.464** | 103 / 154 = **0.669** |
| Random-atlas discoveries that beat all 10 others in their group (9 sibling sets + the real unit) | 3,649 / 11,087 = **0.329** | 667 / 1,381 = **0.483** |
| Unconditional reference (all random tests; 1/11 = 0.091) | 0.087 | 0.087 |
| Mean P(\|t\| real discovery > \|t\| random discovery), same context | **0.482** | **0.476** |
| Same, \|β\| | 0.493 | 0.494 |
| Same, \|t\|, against discoveries of the same unit's sets (n real) | 0.505 (837) | 0.517 (81) |
| Median \|t\|: real against random discoveries | 3.74 / 3.77 | 3.84 / 3.91 |
| Median \|β\|: real against random discoveries | 0.068 / 0.069 | 0.046 / 0.047 |

**Reading.**
- Selection alone lifts the "beats all" rate from about 9% to 33% (H3a) and 48% (H3b).
- Against that selection-aware reference, real discoveries beat their matched sets somewhat more often: 46% against
  33%, and 67% against 48%. This is a modest excess, far smaller than the 46% against 9% contrast implied by an
  unselected reference.
- Conditional on discovery, the effect-size distributions of real complexes and of other-complex fragments are the
  same: P(real > random) ≈ 0.48 for |t| and ≈ 0.49 for |β|.
- The two readings are consistent: the excess is within a group (a real unit against its own 10 matched sets), not
  a marginal effect-size advantage; across discoveries, real complexes are not stronger.
- Per-context values are in `r1b_per_context.tsv`.

## R-2 Lineage scale effect with a non-linear trend (referee M2)

**Method.**
- Per H3a context, a natural cubic spline with 4 df (ns() knots, fixed in `INTENT.md`) of β on unit essentiality,
  fitted over all 367 tested units.
- Uniform rule in every context: `scale_artefact_possible` = sign(β) = sign(fitted β) AND |nonlinear residual z| < 2.
- Corrected lead = T1, not a control readout, |nonlinear residual z| ≥ 2, not `scale_artefact_possible`.

**Results.**
- The spline raises R² only slightly. Lymphoid goes from 0.068 (linear) to 0.120. Every other context stays ≤ 0.083 (CNS/Brain).
- 79 of the 207 H3a strongest findings have |nonlinear residual z| ≥ 2, the same number as with the linear fit.
- 95 are `scale_artefact_possible`, against 56 under D16's Lymphoid/Myeloid-only rule.
- **Corrected leads: 37**, the same count as D16: 35 are kept, 2 are lost and 2 are new.
  - Lost: Bone Z0130 and Myeloid Z0278.
  - New: CNS/Brain Z0146 and Breast Z0253.
- **26 of the 37 corrected leads are scale-concordant.** They are more extreme than the trend in the trend's direction.
  The rule as specified keeps them, because |z| ≥ 2 excludes the artefact flag by construction.
- The lead set is boundary-sensitive: four of the 37 leads lie within 0.05 of |z_nl| = 2 (Bone Z0001, Breast Z0253,
  Lung Z0087, Kidney Z0362), so "35 kept, 2 lost, 2 new" is a count identity, not stability.
- Not done: the referee's M2(c) (exclude scale-concordant leads) is not applied, since the registered-style rule keeps
  them by construction (disclosed as 26/37); M2(d) (rescale statistics by a per-context inflation factor from the
  never-expressed sets) was not carried out; the never-expressed rate is reported as a reliability column instead.
- The flagged Z0190 (Lymphoid, z +5.9) and Z0196, Z0363 and Z0001 stay leads under the non-linear fit. A flexible fit
  does not absorb them.
- Tables: `r2_trend_by_context.tsv`, `r2_strongest_h3a.tsv`, `r2_leads_corrected.tsv`, `r2_summary.json`.

**Per context** (contexts with strongest findings; never-expressed = p < 0.05 rate of the 100 never-expressed sets,
expected 0.05):

| Context | Strongest | \|z_nl\| ≥ 2 | Scale-artefact-possible | Corrected leads | Never-expressed p<0.05 |
|---|---|---|---|---|---|
| Lymphoid | 92 | 19 | 53 | 5 | **0.27** |
| Myeloid | 21 | 10 | 10 | 4 | 0.01 |
| Pancreas | 17 | 6 | 7 | 2 | 0.09 |
| CNS/Brain | 14 | 7 | 7 | 3 | 0.17 |
| Lung | 12 | 6 | 4 | 4 | 0.12 |
| Bone | 11 | 6 | 5 | 3 | 0.19 |
| melanoma | 11 | 6 | 3 | 4 | **0.33** |
| Esophagus/Stomach | 5 | 1 | 2 | 0 | 0.16 |
| Kidney | 5 | 3 | 2 | 3 | 0.19 |
| Peripheral Nervous System | 5 | 3 | 1 | 1 | 0.11 |
| Soft Tissue | 5 | 5 | 0 | 3 | 0.12 |
| Ovary/Fallopian Tube | 3 | 3 | 0 | 1 | 0.22 |
| Breast | 2 | 2 | 0 | 2 | 0.10 |
| Bladder/Urinary Tract | 1 | 0 | 0 | 0 | 0.08 |
| Head and Neck | 1 | 1 | 0 | 1 | 0.19 |
| Pleura | 1 | 1 | 0 | 1 | 0.13 |
| Skin | 1 | 0 | 1 | 0 | **0.28** |

Lineage reliability: 9 of the 37 leads (Lymphoid 5, melanoma 4) are in contexts where the never-expressed sets reach
p < 0.05 at 27–33%. Myeloid leads (4) sit in a context with a well-calibrated never-expressed null (1%).

## R-3 DEMETER2 orthogonal check (referee M2(e))

**Method.**
- The same complex score with removals, and the same HC3 design with the registered covariates available.
- DEMETER2 on the 423 A models with RNAi data. No FDR.

**Results** (`r3_demeter2_per_hit.tsv`, `r3_summary.json`):

| Hits | Evaluable | Same sign | One-sided p < 0.05 (A direction) | Spearman β_A vs β_D2 | Median n_context |
|---|---|---|---|---|---|
| Corrected leads | 35 / 37 | 29 (0.83) | 20 (0.57) | 0.86 | 17 |
| Leads with ≥ 12 context models | 28 | 24 (0.86) | 18 (0.64) | 0.87 | 19 |
| All 207 H3a strongest | 194 | 157 (0.81) | 88 (0.45) | 0.65 | 20 |

**Not evaluable.**
- PNS Z0190: no PNS model has DEMETER2 data.
- Kidney Z0362: the fit failed on the RNAi subset.

**Leads discordant in DEMETER2:** Soft Tissue Z0080 (n_ctx 6), Kidney Z0323, Pancreas Z0283, Soft Tissue Z0336,
Lung Z0001 and Bone Z0001.

**Same sign but weak in DEMETER2** (one-sided p 0.10–0.14):
- Lymphoid Z0190 (β_D2 +0.06, p 0.14);
- melanoma Z0190 (p 0.10);
- Pleura Z0067 (n_ctx 6, p 0.13);
- Lung Z0150 (n_ctx 77, p 0.142) (added after the independent review, 19ad2ae).

## R-4 BRIDGE-L power comparison (referee M3.3)

**Method.**
- The lockbox results review's power aid, reused unchanged (`tweedie` imported from
  `scripts/verify/lockbox_results_review_tables.py`).
- The EB prior is fitted per map and family on all evaluable A tests.
- The L SE is predicted from the A SE and the committed L counts. The L df is the A df scaled by the context counts;
  a sensitivity uses n_ctx,L − 1.
- Rates are standardised by the bridge's size bins.

| Power version | Expected zavg | Expected GLS | Expected Δ (SD) | Observed Δ | Observed − expected (z) |
|---|---|---|---|---|---|
| naive | 0.686 | 0.703 | −0.017 (0.018) | −0.065 | −0.048 (−2.7) |
| EB posterior mean | 0.500 | 0.526 | −0.026 (0.018) | −0.065 | −0.038 (−2.1) |
| **EB mean + variance (headline)** | **0.459** | **0.482** | **−0.023 (0.019)** | **−0.065** | **−0.042 (−2.2)** |

Observed rates: zavg 0.514 (1,267 hits), GLS 0.578 (1,326 hits). The df sensitivity changes the expected values by
≤ 0.004 (`r4_bridge_power.tsv`).

**Reading.**
- Power explains about a third of the BRIDGE-L difference under the headline version; across the three power
  versions the share is 26–41% (naive −0.017, EB mean −0.026, EB mean + variance −0.023, against −0.065). GLS-map discoveries are expected to confirm slightly more
  often, because they have larger |β_A| (median 0.081 against 0.066) and similar |z_A| and L counts
  (`r4_by_map_family.tsv`).
- The remaining −0.042 is about 2.2 SD beyond the power expectation. This SD treats hits as independent; hits share
  units and contexts and the EB prior is estimated, so the z is approximate. The headline model under-predicts the
  observed level of both maps (zavg by 0.055, GLS by 0.097; 0.5785 − 0.4817 = 0.0968, corrected from
  0.096 in referee round 2). GLS hits out-confirm their expected rate in every
  size bin (`r4_by_size_bin.tsv`): observed 0.57–0.59 against expected 0.47–0.50. zavg hits do so less: 0.48–0.57
  against 0.42–0.49.
- So the BRIDGE-L difference is **not fully explained by power**. This is a descriptive reading of one post hoc
  comparison.

## R-5 Referee round 2: three counts from committed tables (POST HOC)

Intent `INTENT_round2.md` (4a7c715) was committed before any value was computed. Script `scripts/62_rev_round2.py`
reads committed tables only. It refits nothing, opens no dataset, lockbox value or test label, and changes no label,
call or rank. Output: `r5_summary.json` plus the three TSVs below.

### A1. H1 neutrality at ±0.05 (referee N1)

Source: `results/part2/directions/H1_within_complex_GI/classes.tsv`, the 11 events labelled NEUTRAL, all TOST-neutral
at the registered ±0.1. **0 of 11** have a committed 90% CI inside ±0.05, and the count is the same with closed
bounds, so no event sits on the boundary. The smallest symmetric margin that would contain each CI ranges from 0.059
(STK11) to 0.097 (KDM6A), median 0.076 (`r5_h1_margin_005.tsv`). This confirms the referee's count. The neutral call
rests entirely on the registered ±0.1 margin. The registered label is unchanged.

### A2. Strongest findings under the alternative registry readings (referee round 1 M7.2)

**H1 LOSO unit** (registered reading: lineage; alternative: primary disease, the H3a unit).
- The unit applies only to H1 rows. H3a (207) and H3b (45) are unchanged.
- Of the 5 H1 eligible calls, only the 2 candidates can change: ARID1A P and EP300 P. The other three fail on
  non-LOSO criteria. FOCAD N is defect-sensitive, PTEN N was not confirmed, and TP53 N is not headline-eligible
  under E2.
- Bound under the alternative reading: H1 0–2, total **252–254**. Whether either call changes is **not determined**,
  because the per-fold H1 fits by primary disease are not committed and no A refit was made in this round. Under the
  registered unit, both calls passed all 30 lineage folds with none skipped.
- If an H1 call dropped, every finding ranked below it would move up one place.

**Ranking** ("A statistic × tier"; registered reading: tier class, then |t_A|). Ranking selects nothing, so all 254
findings stay the same under any reading. Only the order changes (`r5_ranking_alternatives.tsv`):

| Reading | Rows whose rank changes | Spearman vs registered | Top-10 / top-30 overlap | Median / max rank shift |
|---|---|---|---|---|
| alt1: \|t_A\| alone | 254 | 0.37 | 5 / 17 | 58.5 / 217 |
| alt2: \|t_A\| × tier weight (T1 3, T2 2, T3 1; illustrative, not registered) | 225 | 0.89 | 9 / 28 | 10 / 174 |

The rank of each highlighted finding under registered / alt1 / alt2:
- Ovary Z0035: 85 / 143 / 98
- Bone Z0303: 53 / 100 / 60
- Head and Neck Z0051: 19 / 34 / 21
- Myeloid Z0254: 49 / 94 / 54
- RB1 loss Z0069: 7 / 13 / 8
- ARID1A P: 60 / 110 / 67
- EP300 P: 229 / 66 / 224

Reading: the H1 LOSO reading can change at most the two H1 calls, and the ranking reading changes only the order. Of
the alternatives, ignoring tier reorders the list the most. No count in the text depends on rank.

### A3. Spread of the R-1(b) selection-aware reference across the 10 random atlases (referee N5)

Per atlas, the share of random-atlas discoveries that beat all ten others in their group is given below
(`r5_atlas_flag_range.tsv`). The pooled rates reproduce `r1_summary.json`: 3,649 / 11,087 and 667 / 1,381.

| | Per-atlas range | Median | Pooled | Real discoveries | Atlases at or above the real rate |
|---|---|---|---|---|---|
| H3a | 0.291–0.364 | 0.333 | 0.329 | 0.464 | 0 of 10 |
| H3b | 0.448–0.517 | 0.482 | 0.483 | 0.669 | 0 of 10 |

This is a spread across atlases, not a confidence interval for the real rate. Real discoveries and atlas
discoveries share units, contexts and lines.

## Files

| File | Contents |
|---|---|
| `INTENT.md` | Pre-computation intent (committed first) |
| `r1_summary.json` | R-1 checks, matching status, BB, replicate counts, conditional summaries, runtime, memory |
| `r1a_matching_quality.tsv`, `r1a_rates.tsv`, `r1a_nonmember_tests.tsv.gz` | R-1(a) |
| `r1b_per_context.tsv`, `r1b_real_discoveries.tsv`, `r1b_random_discoveries.tsv.gz` | R-1(b) |
| `r2_trend_by_context.tsv`, `r2_strongest_h3a.tsv`, `r2_leads_corrected.tsv`, `r2_summary.json` | R-2 |
| `r3_demeter2_per_hit.tsv`, `r3_summary.json` | R-3 |
| `r4_bridge_power.tsv`, `r4_by_map_family.tsv`, `r4_by_size_bin.tsv`, `r4_per_hit.tsv.gz`, `r4_summary.json` | R-4 |
| `INTENT_round2.md` | Round-2 intent for R-5 (committed first, 4a7c715) |
| `r5_h1_margin_005.tsv`, `r5_ranking_alternatives.tsv`, `r5_atlas_flag_range.tsv`, `r5_summary.json` | R-5 (round 2) |

## Reproduce

```
export PART1_DIR=".../coessentiality_protein complexes"; export PYTHONHASHSEED=0
.venv/bin/python scripts/58_rev_matched_sets.py --workers 7
.venv/bin/python scripts/59_rev_scale_effect.py
.venv/bin/python scripts/60_rev_demeter2.py        # needs r2_strongest_h3a.tsv
.venv/bin/python scripts/61_rev_bridge_power.py
.venv/bin/python scripts/62_rev_round2.py          # R-5; committed tables only
```
