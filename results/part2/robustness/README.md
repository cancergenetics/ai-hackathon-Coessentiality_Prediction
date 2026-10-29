# results/part2/robustness

This folder holds the registered per-hit robustness checks (LOSO, leave-top-2-out, breadth) and the final
strongest-findings list. The rule is in DESIGN v3.1 §2.5 and the `success_criterion` of H3a, H3b and H1. It is
implemented in `src/part2/strongest.py` and run by `scripts/56_strongest_findings.py`. The intent note
(`INTENT.md`) and the candidate list (`candidates.tsv`) were committed in fd03b03, before any robustness value was
computed.

**PI ruling (after the run).** The 254 findings that pass the registered rule are the primary strongest-findings
list. The PI accepted two readings, which are disclosed here:
- the H1 LOSO unit is the lineage;
- findings are ranked by tier first, then by |t_A|.

## Disclosure: computed after the lockbox read

**These checks were computed after the single lockbox read (5bc9444).** They were registered before it: the
strongest-findings rule, LOSO and leave-top-2-out are in DESIGN §2.4–2.5 and in each entry's `multiverse.robustness`.
They were not run at the discovery runs ("deferred for a later step"), which is a registry-to-code gap (LESSONS rule
15). Running them now does not change any discovery, confirmation, registry entry or ledger, and it is not a trial.

- **Data used:** the checks use A models only. No lockbox value, persisted read or loader was opened.
- **Columns copied from committed tables:** the L columns (`beta_L`, CI, `p_L`, `n_context_L`) come from
  `results/part2/lockbox_confirm/<ID>/confirmation.tsv`.
- **Remaining forking path:** the rule, the LOSO units, the candidate list and the ranking were fixed in `INTENT.md`
  before computation. What remains is a forking path in *reading* ambiguous registry text, not in choosing analyses.
  The two ambiguous points are the H1 LOSO unit and the meaning of "A statistic × tier", both described below.

## Rule (applied mechanically)

A finding is strongest if all of these hold:

- discovery q < 0.05;
- the PRIMARY lockbox call is "confirmed" or "rare-context route passed" (labelled in `route`);
- it is not defect-sensitive (R7);
- it is not a locus cluster;
- it is **LOSO-robust**: every evaluable fold has the same sign as the full A fit and p < 0.05 under the discovery
  test's alternative. A fold that leaves fewer than 2 context or 2 reference models is skipped and counted, not
  failed.

For H1, the call must also be non-neutral (eligible) and headline-eligible under E2.

The LOSO unit for each direction:

| Direction | LOSO unit |
|---|---|
| H3a | OncotreePrimaryDisease (registry text) |
| H3b | OncotreeLineage ("leave-one-lineage-out (LOSO)") |
| H1 | OncotreeLineage |

The H1 entry writes "LOSO" without expanding it. Every entry that expands the term reads it as leave-one-lineage-out,
and the H1 rare route used the same unit.

- **Leave-top-2-out and breadth:** reported for every hit, not used to select hits. The registry `use` field says
  "reported per hit, never selected", and the success criterion names LOSO only.
- **Ranking:** "A statistic × tier" has no registered numeric tier weights. Hits are ranked by tier class first (T1,
  then T2, then T3-neg/T3-nc), then by |t_A| in descending order.
- **Reproduction check:** the full A refit reproduced every committed discovery β (to 1e-6) and t (to 1e-4) for all
  575 hits.

## Results

| Direction | Candidates | LOSO-robust (strongest) | Units | Leave-top-2 robust among strongest | Manuscript "specific" subset surviving |
|---|---|---|---|---|---|
| H3a | 524 | **207** | 137 | 207 | 107 / 261 |
| H3b | 46 (41 lockbox, 5 rare route) | **45** | 41 | 45 | 28 / 29 |
| H1 | 2 (5 eligible calls computed) | **2** (ARID1A→ARID1B P; EP300→CREBBP P, rare route) | 2 | 2 | – |

The strongest findings by tier are:

| Direction | T1 | T2 | T3-neg |
|---|---|---|---|
| H3a | 124 | 56 | 27 |
| H3b | 27 | 9 | 9 |
| H1 | 1 | – | 1 |

The one H3b failure is MSI_H Z0334, which failed 1 of 30 lineage folds.

**H3a LOSO behaviour.**

- Of the 317 H3a failures, 308 fail exactly one of the 81–82 primary-disease folds.
- Failures cluster in lineages dominated by one primary disease. Examples:
  - Bowel: 0 of 29 robust. Colorectal adenocarcinoma is 60 of the 64 Bowel A models.
  - Head and Neck: 1 of 54 robust. HNSCC is 56 of 58.
  - Skin: 1 of 22 robust. Melanoma is 44 of 51.
  - Peripheral Nervous System: 5 of 44 robust. Neuroblastoma is 31 of 35.
- In these lineages, leaving out the dominant disease leaves only a handful of context models.
- When fewer than 2 context models remain, the fold is skipped instead of failed, as the engine rule requires. This is
  why Pancreas passes 17 of 17: 52 of its 53 models are one disease.
- The registered rule therefore tests "not carried by one primary disease". It is strict for lineages with a small
  non-dominant remainder and vacuous when the remainder has fewer than 2 models. This is a property of the registered
  rule, reported here, not altered.
- Lymphoid (92 of 119 robust) and Myeloid (21 of 37) have several sizeable diseases and pass more often.

Three caveats from the manuscript carry over to this list:

- **Lymphoid scale effect.** 55 of the 92 lymphoid strongest findings have positive β (lower lymphoid dependency on
  the unit). The manuscript's lineage-scale caveat applies to them (main text, Results, "Possible lineage-scale
  artefacts", and Supplementary Note 16; this README earlier cited "§4.5" of a superseded draft). The revision analysis
  R-2 (`results/part2/revision/`, post hoc) replaces the Lymphoid/Myeloid-only label with one rule for all lineages.
- **Control readouts.** The literal `control` column marks 2 rows "control readout" (H3b MSI_H Z0157, PC1; TP53_altered
  Z0166, PC2). The post hoc
  `control_recovery` annotation below, which also counts by-definition equivalent contexts (e.g. H1loss:ARID1A =
  LoF:ARID1A, PC6; Lymphoid/Myeloid within heme, PC20), marks **7** readouts. The manuscript uses the annotation count
  (7). *Corrected 2026-10-29: this caveat gave only the literal count (2), which did not match the annotation count (7).*
- **Control units in other contexts.** The literal `control` column marks 21 rows "control unit, other context"; the
  `control_recovery` annotation, after moving equivalent contexts to "readout", gives 18.

## Post hoc descriptive annotations (PI ruling; not filters, not tests)

These annotations were added after the rule was applied, at the PI's request. **They are post hoc and descriptive.
They select nothing, change no rank and are not hypothesis tests.** They are computed by
`scripts/57_annotate_strongest.py` from committed A-side files: discoveries, calibration_D8 tests, H3a `tests.tsv.gz`
and `controls_summary.tsv`.

The one exception is unit essentiality, which needs the A gene-effect matrix. It is computed once on A models only
(`--essentiality`) and committed as `unit_essentiality_A.tsv`. Every annotation can therefore be re-derived from
committed files. No L value is used.

| Column | Meaning |
|---|---|
| `beats_all_matched_sets` | H3: \|t\| exceeds every matched random set of the unit in that context (calibration_D8, D8 flag). H1: NA, because H1 uses its own matched null (`classes.tsv`). |
| `control_recovery` | `readout (PCx)`: the (context, complex unit) pair is a registered complex-level control readout. By-definition equivalent contexts count too, e.g. H1loss:ARID1A = LoF:ARID1A (PC6) and lineage:Lymphoid/Myeloid ⊂ heme (PC20); the map is `strongest.CONTROL_CONTEXT_OVERLAP`. `same unit, other control context (...)`: the unit is a control readout in a context that does not overlap. Empty otherwise. |
| `loso_partially_testable`, `loso_folds_tested`, `loso_folds_skipped` | At least one LOSO fold was skipped because fewer than 2 context or reference models were left. |
| `unit_essentiality_A` | Mean over members of the gene's mean clipped A effect (`calib.clip_essentiality`), as in the calibration matching. |
| `scale_effect_context_slope`, `scale_effect_context_rho` | H3a only. OLS slope and Spearman ρ of β on unit essentiality over all 367 evaluable tested units of that context (`tests.tsv.gz`). Full table: `scale_effect_by_context.tsv`. |
| `scale_effect_fitted_beta`, `scale_effect_concordant` | Fitted β = intercept + slope × essentiality. Concordant = sign(β) equals sign(fitted β). Where the slope is near 0, the fitted sign is set mostly by the context's intercept, so "concordant" there means "same sign as the lineage-wide shift". |
| `scale_effect_residual_beta`, `scale_effect_residual_z` | β minus fitted β, and that residual divided by the residual SD of the context. |
| `locus_flag`, `tier` | For completeness. All 254 have `locus_flag` = False. |

### Counts over the 254

| Annotation | Count |
|---|---|
| Tier | T1 152, T2 65, T3-neg 37 |
| `beats_all_matched_sets` | H3a 107/207, H3b 34/45 |
| `control_recovery` = readout | 7: H1 ARID1A P (PC6); H3a Lymphoid and Myeloid Z0113 (PC20); H3b MSI_H Z0286 (PC10), TP53_altered Z0166 (PC2), NRAS Z0132 (PC4, PC5), MSI_H Z0157 (PC1) |
| `control_recovery` = same unit, other control context | 18 |
| `loso_partially_testable` | H3a 21 (Pancreas 17, Ovary 3, Pleura 1); H3b 33 (TP53_altered 28, EWSR1 fusion 5); H1 0 |

### Scale effect

The Lymphoid ρ of −0.24 reproduces the review's figure. Myeloid shows the same direction, weaker (ρ −0.135, slope
−0.029). The other lineages have |ρ| ≤ 0.21, mostly positive: Bone +0.19, CNS +0.11, Head and Neck +0.20,
melanoma +0.18, Skin +0.21.

Overall, 142 of the 207 H3a findings are scale-effect-concordant and 79 have |residual z| ≥ 2. 55 findings are both.

| Context | Concordant / n | \|residual z\| ≥ 2 |
|---|---|---|
| Lymphoid | 64/92 | 18 |
| Myeloid | 11/21 | 11 |
| Pancreas | 7/17 | 6 |
| CNS/Brain | 12/14 | 6 |
| Lung | 11/12 | 6 |
| Bone | 8/11 | 7 |
| melanoma | 7/11 | 6 |
| Kidney | 5/5 | 4 |
| Esophagus/Stomach | 2/5 | 1 |
| PNS | 4/5 | 3 |
| Soft Tissue | 3/5 | 5 |
| Ovary | 3/3 | 3 |
| Breast | 2/2 | 1 |
| Bladder | 0/1 | 0 |
| Head and Neck | 1/1 | 1 |
| Pleura | 1/1 | 1 |
| Skin | 1/1 | 0 |

A concordant hit with a large |residual z| has a β larger than the lineage's essentiality trend predicts. An example
is Lymphoid Z0190, the proteasome β-ring/E4F1 unit, with residual z +5.9.

### Top 15 by tier, then |residual z|

This is a **post hoc, descriptive view of the H3a findings. It is not the registered ranking.** The registered rank is
given in the first column.

| Reg. rank | Context | Unit | Members | Tier | A β | L β | Residual z | Concordant |
|---|---|---|---|---|---|---|---|---|
| 2 | Lymphoid | Z0190 | E4F1, PSMB5, PSMB6, PSMB7, THAP1 | T1 | +0.484 | +0.565 | +5.88 | yes |
| 46 | Pleura | Z0067 | B3GALT6, B3GAT3, B4GALT7, EXT1, EXT2… | T1 | −0.265 | −0.182 | −4.95 | yes |
| 85 | Ovary/Fallopian Tube | Z0035 | AFDN, ARHGAP29, CDC42, CDH1, CDH2… | T1 | −0.140 | −0.202 | −4.47 | yes |
| 81 | melanoma | Z0190 | E4F1, PSMB5, PSMB6, PSMB7, THAP1 | T1 | −0.218 | −0.274 | −4.00 | yes |
| 41 | Soft Tissue | Z0080 | CCT4, MSTO1, PDCD5, PDCL3, PFDN1… | T1 | −0.156 | −0.091 | −3.84 | yes |
| 53 | Bone | Z0303 | SNAPC1, SNAPC2, SNAPC3, SNAPC4, SNAPC5 | T1 | −0.210 | −0.251 | −3.74 | yes |
| 49 | Myeloid | Z0254 | AGFG1, AP2A1, AP2B1, AP2M1, AP2S1… | T1 | −0.189 | −0.239 | −3.09 | yes |
| 15 | CNS/Brain | Z0254 | AGFG1, AP2A1, AP2B1, AP2M1, AP2S1… | T1 | +0.138 | +0.091 | +3.08 | yes |
| 76 | Lung | Z0174 | ANKRD49, APC, AXIN1, CSNK1A1, CTNNBIP1… | T1 | +0.086 | +0.077 | +3.01 | yes |
| 83 | Peripheral Nervous System | Z0190 | E4F1, PSMB5, PSMB6, PSMB7, THAP1 | T1 | −0.171 | −0.184 | −2.89 | no |
| 17 | Lymphoid | Z0014 | BTAF1, CCNC, DRAP1, KAT2A, KDM5C… | T1 | −0.200 | −0.174 | −2.63 | yes |
| 55 | CNS/Brain | Z0081 | CHMP5, CHMP6, HGS, PTPN23, SNF8… | T1 | +0.106 | +0.103 | +2.60 | no |
| 127 | Kidney | Z0323 | ARL8B, SNAPIN, SPG21, VPS16, VPS18… | T1 | −0.143 | −0.086 | −2.58 | yes |
| 19 | Head and Neck | Z0051 | CALM1, CCNJ, CHCHD2, EGFR, FKBP9… | T1 | −0.133 | −0.107 | −2.58 | yes |
| 1 | Lymphoid | Z0016 | ATRN, BCAR1, CRK, CRKL, DOCK5… | T1 | +0.182 | +0.173 | +2.53 | no |

## Files

| File | Contents |
|---|---|
| `INTENT.md` | Rule, LOSO units, candidate list (committed before computation) |
| `candidates.tsv` | Every H3a/H3b discovery and every H1 eligible call, with the non-robustness criteria; `writer_specific_set` reproduces the manuscript's 261/29 |
| `robustness_per_hit.tsv` | Per hit: `loso_folds`, `loso_pass`, `loso_skipped`, `loso_robust`, `top2_beta`, `top2_p`, `top2_same_sign`, `top2_robust`, `breadth` (H3), `strongest`, `rank` |
| `strongest_findings.tsv` | The 254 strongest findings, ranked; A β/t, L β with 95% CI (committed), LOSO and leave-top-2 outcomes, members; plus the post hoc descriptive annotation columns |
| `summary.json` | Counts, runtime (118 s wall, 2.5 GB peak RSS) |
| `unit_essentiality_A.tsv` | Post hoc: unit essentiality on A (mean clipped A effect over members), 367 coarse units |
| `scale_effect_by_context.tsv` | Post hoc: per H3a context, slope / intercept / Spearman ρ of β on essentiality over all tested units, residual SD |
| `annotations_summary.json` | Post hoc annotation counts (per direction, per context) |

## Reproduce

Set `PART1_DIR` and `PYTHONHASHSEED=0`, then run:

```
.venv/bin/python scripts/56_strongest_findings.py --candidates   # committed files only
.venv/bin/python scripts/56_strongest_findings.py --run          # needs INTENT.md + candidates.tsv committed
```

## Not covered

- H1S_silencing_GI (exploratory by registry standing) and H2_submodule_divergence (exploratory: its kill rule is n/a)
  are outside the strongest-findings rule.
- The other multiverse items (alternative scores, covariate sensitivity) remain not run, as before.
