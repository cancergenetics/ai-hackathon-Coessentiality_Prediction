# Complex-specificity counts (POST HOC, DESCRIPTIVE)

**This is a post hoc descriptive filter. It is not part of the registered error control.** The registered discoveries,
q-values, thresholds and registry are unchanged. This filter does not define, promote or demote any discovery, and
lockbox confirmation does not use it. It was requested after the discovery and calibration outcomes were seen. It
writes no trial line and does not read the lockbox.

## Question

A discovery may reflect something specific to the complex, or only the background that any similar gene set shows in
that context. For each primary discovery we ask: does its β lie outside the central 95% of the β values of its own
matched random sets, tested in the same context? "Outside" means strictly below the 2.5% quantile or above the 97.5%
quantile, using numpy linear quantiles over the finite values.

The random sets are the registered D8 calibrations (`CALIBRATION_D8_INTENT.md`, commits `fd582a6`, `4fb83d7`,
`1a2a445`):
- **H3a and H3b.** Matched random gene sets: R = 10 per coarse unit, matched on size, essentiality and coherence.
  Each is scored with the context's removals and tested with the context's covariates. With 10 sets, the 2.5% and
  97.5% quantiles lie close to the minimum and maximum.
- **H2.** Random member splits of the same parent: R = 20 per tested (context, sub-module), residualised on the
  complement as for the real sub-module.

## Source of β (disclosed)

The committed `calibration_D8_tests.tsv.gz` stores t and p for each random set, but not β. Two versions are
therefore given:
- **β (requested).** The same sets were regenerated deterministically with the committed code and seeds
  (`calib.matched_random_sets` / `calib.random_split_scores`, sub-streams as in `scripts/46_run_calibrations.py`)
  and re-fitted on A. The regenerated t equals the committed t in every row: 80,500 rows for H3a, 70,200 for H3b and
  2,180 for H2, with max relative difference 5.0e-6, which is the `%.6g` rounding. So the β quantiles come from
  exactly the committed sets.
- **t, from the committed files only.** This is the same filter on the HC3 t statistic.

Script: `scripts/49_specificity_counts.py`. Outputs: `specificity.json` and `specificity_discoveries.tsv.gz`, one row
per discovery. Runtime 115 s, peak RSS 3.4 GB.

## Totals

| Direction | Discoveries | With matched sets | β outside 95% | (same side as β) | t outside 95% |
|---|---|---|---|---|---|
| H3a | 1,145 | 1,083 | **699** | 698 | 718 |
| H3b | 161 | 154 | **124** | 123 | 128 |
| H2 | 109 | 109 | **78** | 78 | 51 |

Units without matched sets are the 17 H3a and 16 H3b units marked "unavailable" in the D8 calibration. They account
for 62 H3a discoveries and 7 H3b discoveries, which are not counted.

## H3a, per lineage

| Context | Discoveries | With sets | β outside | t outside |
|---|---|---|---|---|
| Lymphoid | 197 | 187 | 79 | 79 |
| Head and Neck | 126 | 121 | 65 | 64 |
| CNS/Brain | 107 | 99 | 60 | 60 |
| melanoma | 91 | 87 | 60 | 61 |
| Peripheral Nervous System | 86 | 79 | 54 | 58 |
| Myeloid | 81 | 77 | 43 | 50 |
| Bone | 69 | 65 | 49 | 50 |
| Bowel | 62 | 62 | 44 | 44 |
| Skin | 57 | 54 | 36 | 37 |
| Breast | 55 | 50 | 40 | 39 |
| Lung | 48 | 46 | 41 | 40 |
| Pancreas | 48 | 43 | 38 | 38 |
| Kidney | 38 | 37 | 29 | 28 |
| Esophagus/Stomach | 30 | 30 | 22 | 25 |
| Uterus | 11 | 11 | 9 | 11 |
| Pleura | 9 | 8 | 7 | 8 |
| Soft Tissue | 9 | 7 | 7 | 7 |
| Ovary/Fallopian Tube | 8 | 8 | 7 | 8 |
| Bladder/Urinary Tract | 6 | 5 | 4 | 5 |
| Cervix | 4 | 4 | 3 | 3 |
| Liver | 3 | 3 | 2 | 3 |

Biliary Tract and Prostate have no discoveries.

Lymphoid has the most discoveries but the lowest β-specific fraction: 79 of 187, or 42%. Other lineages with at least 20 discoveries are at
54–89%. This fits the background-shift reading of the D8 calibration (lineage contrasts move matched random sets too).

## H3b, per context

| Context | Discoveries | With sets | β outside | t outside |
|---|---|---|---|---|
| TP53_altered | 92 | 86 | 63 | 63 |
| WGD | 26 | 26 | 22 | 24 |
| MSI_H | 14 | 14 | 13 | 13 |
| fusion:EWSR1::FLI1\|ERG | 11 | 10 | 10 | 10 |
| hotspot:BRAF | 4 | 4 | 3 | 4 |
| hotspot:KRAS | 3 | 3 | 3 | 3 |
| 9p21_loss | 2 | 2 | 2 | 2 |
| LoF:PTEN | 2 | 2 | 2 | 2 |
| LoF:RB1 | 2 | 2 | 2 | 2 |
| NRF2_active | 2 | 2 | 2 | 2 |
| LoF:SMARCA4 | 1 | 1 | 0 | 1 |
| hotspot:NRAS | 1 | 1 | 1 | 1 |
| sex_female | 1 | 1 | 1 | 1 |

## H2 (exploratory), per context

| Context | Discoveries | β outside | t outside |
|---|---|---|---|
| Lymphoid | 31 | 20 | 13 |
| CNS/Brain | 11 | 9 | 4 |
| Head and Neck | 9 | 6 | 4 |
| Myeloid | 8 | 7 | 5 |
| melanoma | 8 | 5 | 2 |
| Bone | 6 | 3 | 2 |
| Kidney | 6 | 4 | 3 |
| TP53_altered | 5 | 3 | 3 |
| Peripheral Nervous System | 5 | 5 | 3 |
| Bowel | 4 | 3 | 2 |
| Pancreas | 4 | 3 | 3 |
| Skin | 4 | 3 | 2 |
| WGD, hotspot:BRAF, Esophagus/Stomach, Pleura, Soft Tissue | 1 each | 1 each | 1 each |
| hotspot:KRAS, Breast | 1 each | 1 each | 0 |
| Lung | 1 | 0 | 0 |

For H2 the β and t filters disagree more (78 against 51). Random splits have different residual variance from the
real sub-module, so their standard errors differ. The β filter compares effect sizes; the t filter compares evidence.

## Reading

These counts describe how many registered discoveries also stand out against similar gene sets in the same context.
They are not an error rate. With 10 sets per unit (H3a/H3b), the quantiles are coarse. A discovery inside the band is
not shown to be non-specific. A discovery outside the band is not confirmed.
