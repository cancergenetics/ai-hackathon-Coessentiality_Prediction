# Table 2: complex maps on the construction set A

| quantity | zavg_BEF map | GLS map |
|---|---|---|
| level-1 complexes | 367 | 379 |
| genes in level-1 complexes (of S) | 3,642 / 6,022 | 2,862 / 6,022 |
| sub-modules (parents with any) | 132 (74) | 108 (65) |
| level-1 size, median (IQR; max) | 9 (5–13; 29) | 5 (4–9; 37) |
| level-1 size bins 3–4 / 5–8 / 9–16 / 17–60 | 85 / 85 / 149 / 48 | 152 / 121 / 78 / 28 |
| coverage of pan-essential S genes | 594/965 (61.6%) | 410/965 (42.5%) |
| coverage of selective S genes | 1,722/2,269 (75.9%) | 1,494/2,269 (65.8%) |
| coverage of other S genes | 1,326/2,788 (47.6%) | 958/2,788 (34.4%) |
| locus-flagged, level 1 / level 2 | 28.3% / 23.5% | 22.2% / 25.0% |
| tiers level 1: T1 / T2 / T3-neg / T3-nc | 159 / 74 / 124 / 10 | 159 / 81 / 131 / 8 |
| expected false T1 (level 1) | 8.0 | 8.0 |
| tiers level 2: T1 / T2 / T3-neg / T3-nc | 59 / 28 / 42 / 3 | 42 / 24 / 38 / 4 |
| expected false T1 (level 2) | 3.0 | 2.1 |
| non-locus level 1: T1 / T2 / T3-neg / T3-nc | 155 / 60 / 47 / 1 | 151 / 69 / 69 / 6 |
| best-match Jaccard in bootstrap maps, pooled median (fraction ≥ 0.5) | 0.69 (0.73) | 0.67 (0.71) |
| random size-matched sets: median (P95) | 0.06 (0.12) | 0.10 (0.17) |
| complexes with median Jaccard > random P95 | 96.5% | 97.9% |
| pairwise ARI of bootstrap maps, median (IQR) | 0.754 (0.745–0.764) | 0.747 (0.735–0.763) |
| bridge: size-standardised BioPlex precision (matched coverage) | 0.157 | 0.182 |

## Notes

- Sources: results/part2/construction/summary/summary.json, results/part2/construction/tiers/summary.json, results/part2/construction/variability/variability.json, results/part2/construction/bridge/bridge.json. Maps frozen at commit 27fa9b2 (FZ-C).
- Bridge primary difference zavg − GLS: -0.025 [-0.075, +0.031] (protein-node bootstrap, B = 2000, seed 5000); decision: no difference claimed.
- Coverage classes by A-essentiality: pan-essential = effect < −0.5 in ≥ 90% of A models, selective 5–90%.
- Tiers use BioPlex 3.0 and CCLE proteomics (A models only), read only after the maps were committed. T2 is not locus-aware (copy-number dosage makes neighbouring proteins co-abundant).
- Variability: 20 single-subsample maps (PatientID-blocked 80% subsamples of A). The frozen map is built from the same subsamples, so agreement with it is not independent; the random baseline is a floor, not a null.
