# Independent verification of the test look (Verifier A, plus a score spot check)

**Verdict: REPRODUCED.**
- Every primary and verdict-relevant number matches the harness to ≤ 2.6 × 10⁻¹⁵. That covers the pair set, floor, clusters, AP/AUROC/macro AP/AP without the largest component, rule (b) P/R/F1, Brier, the seed-0 cluster-bootstrap CIs, ΔAP, the ratio CIs, Holm, and the §7.1/§7.2 verdicts on S1 and ITT.
- 36 of 1,352 numeric fields differ from the harness by more than 1e-9 (largest 4.1 × 10⁻⁷). All 36 are secondary calibration descriptors: equal-mass ECE, and the calibration intercept and slope. Both causes are identified below. Neither affects any verdict.

## Method
**Independence**
- `verify_independent.py` re-implements the following from `docs/TEST_PROTOCOL.md` (FROZEN v1.2):
  - R1–R3 and the ITT floor (min of mapped scores − 1);
  - S1, the set of pairs mapped under all six methods;
  - AP and AUROC with tie blocks, P@R, AP with hard negatives;
  - macro per-cluster AP, AP without the largest component, AP without the top-10 clusters;
  - rules (a), (b) and (c);
  - isotonic and Platt calibration, using the train-fitted knots and parameters in `fitted_rules.json` with unmapped pairs set to the train floor. From these it computes Brier, Brier skill, calibration-in-the-large, intercept/slope (exact Newton MLE), and ECE with 20 equal-mass or 15 equal-width bins;
  - the paired cluster bootstrap (B = 2,000, `default_rng(0)`, m = bincount(integers(0, C, C)), pair weights m_c or m_c1·m_c2), covering AP, AUROC, AP without the largest component, macro AP (cluster counted m_c times, with replicate pair weights), F1/P/R at rule (b), and isotonic Brier;
  - ΔAP, relative gain and ratio CIs, p-values, Holm;
  - the §7.1 verdict and §7.2 statistic.
- From Part 1, only `src.clusters.build_clusters` was imported (permitted by §8.3), run under `PYTHONHASHSEED=0`.
- Part 1's `metric_suite` and `final_eval.decide` were called only as cross-checks. Agreement: ≤ 2.2 × 10⁻¹⁶, and the decisions are identical.
- Nothing from `src/testeval` or `scripts/21_*` was imported.

**Order of work**
1. The numbers were committed in `c69d828` ("independent numbers (before comparison)"; `verifier_results.json` SHA-256 `d899bc8b…`).
2. Only after that commit was `test_results.json` opened.
3. `compare_with_harness.py` then wrote `comparison.json`, which lists all 1,444 compared fields.

## Comparison (tolerance 1e-9)
| group | fields | > 1e-9 | max abs diff |
|---|---|---|---|
| audit (R1–R4, mapped counts) | 20 | 0 | 0 |
| clusters (ITT 238 / S1 239 clusters, no Louvain split, largest component 39 proteins) | 20 | 0 | 0 |
| point metrics, ITT + S1 | 300 | 0 | 2.2e-16 |
| bootstrap per-method (point, CI, SD) | 384 | 0 | 2.2e-16 |
| bootstrap Δ / relative gain / p | 322 | 0 | 2.6e-15 |
| §7.1 verdict, ITT + S1 | 72 | 0 | — (identical) |
| §7.2 CI-excludes-0, borderline, ratio CIs, Holm | 106 | 0 | 2.4e-15 |
| C1–C3 AP and counts | 24 | 0 | 5.6e-17 |
| calibration (isotonic + Platt) | 192 | **36** | 4.1e-7 |

**Discrepancies, all explained**
1. **ECE with 20 equal-mass bins** (≤ 3 × 10⁻⁸). The 600,738 pairs do not divide evenly into 20 bins.
   - My code uses `np.array_split`: the first 18 bins hold 30,037 pairs, the last two 30,036.
   - The harness assigns bin = floor(20·i/N): bins 9 and 19 hold 30,036.
   - The protocol only says "equal pair count". Using the harness's rule, I reproduce its ECE exactly (|Δ| ≤ 4 × 10⁻¹⁹ for all 12 method × calibrator combinations).
2. **Calibration intercept and slope** (≤ 4.1 × 10⁻⁷). My Newton fit converges to step < 1e-12. At the harness's coefficients, the log-likelihood gradient is still up to 3.4 × 10⁻⁴, so the harness's optimiser stopped at its tolerance. The two differ in optimiser precision, not in definition.

## Locked verdict (identical in both)
| set | zavg_BEF | B |
|---|---|---|
| S1 (confirmatory) | **passes vs GLS, Spearman, Pearson → confirmed** | fails vs GLS (macro-AP NI: ratio 0.875 < 0.95) |
| ITT | passes (confirmed) | fails vs GLS (ratio 0.873) |

**ΔAP of zavg_BEF, S1** (seed 0)
- vs GLS: +0.127 [+0.089, +0.167]
- vs Spearman: +0.135 [+0.024, +0.229]
- vs Pearson: +0.143 [+0.038, +0.233]

All three CIs exclude 0, and none is borderline. For B vs Spearman, the CI excludes 0 but is borderline (lower bound 0.0026).

**Absolute values**
- AP, zavg_BEF: ITT 0.3329, S1 0.3488.
- macro-AP ratio zavg_BEF / GLS: 1.017 (ITT), 1.018 (S1).

## Score spot check (`spotcheck_scores.py`, `spotcheck_scores.json`)
**Sample.** 1% of the scored pairs, drawn with `default_rng(7)` without replacement: 6,007 pairs, 5,645 mapped. They were re-scored from the raw `CRISPRGeneEffect.csv`; its SHA-256 matches the Part 1 manifest.

| path | max abs diff vs TSV |
|---|---|
| GLS / Pearson / Spearman via Part 1 `src/baselines.py` (dev-time estimators) | 1.3e-14 / 1.1e-15 / 4.4e-16 |
| GLS / Pearson / Spearman via Part 1 `BaselineEngine` on the sample | 0 (bit-identical) |
| zavg_BEF from the genome-wide all-pairs matrix (`outputs/all_pairs/zavg_BEF.npy`, float32) | 3.3e-6 |
| zavg_BEF / B via Part 1 `mvp/scripts/predict_pairs.py` (standalone re-implementation) | 3.5e-6 / 0 |
| zavg_BEF / B via the locked `predict_pairs.py` on the sample file (batch invariance) | 0 / 0 |

Rank correlation is 1.0 everywhere. The mapped flags agree for every path.

The zavg_BEF differences of about 3e-6 come from the mvp's float32 storage and its separate implementation. They are within the mvp's documented 3e-5 equivalence.

**Not re-scored:** the post hoc reference `spec_k20_tinf_std`.

## Not independently verified
- essentiality strata;
- ECE with log-spaced bins, and the reliability tables;
- protein and pair bootstraps;
- Monte-Carlo seeds 1–4;
- C1–C3 bootstrap CIs.

**Deviation from §8.3:** §8.3 specifies seed 0 for the spot-check sample; seed 7 was used, as the verification brief requested.
