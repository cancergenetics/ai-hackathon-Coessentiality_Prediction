# H1_within_complex_GI: A-side discovery (single look)

**Commits.** Intent `ef41ddb`, then results `d12c627`. The SMARCB1 copy-number-covariate sensitivity is `de9bc3f`.
The kill rule passes: 3 of 3 controls, which is thin coverage.

**Discovery.**
- 24 events and 26 component tests. m = 24 evaluable tests; the registered m is 28.
- 5 tests passed the screen and were re-tested by permutation.
- 5 BH discoveries.

**PTEN depends on the evaluability-adjusted m.** The engine counts only evaluable tests in m, as the registry allows
("non-evaluable tests are listed and excluded from m"). Here m = 24:
- 2 N axes have |N'| = 1 after the baseline rule (EP300, SMAD4);
- MAP2K4 has 11 loss models.

Z0130|PTEN|N is a discovery at m = 24 (q = 0.0437). At the registered m = 28 it would be lost (q = 0.051). The other
four discoveries are kept either way (`../h1_label_final.json`).

**Classes and reported labels.** In `classes.tsv`, the engine classes are in `class_complex_wide` and
`class_paralog`. The reported labels are in `label_final` and `label_final_paralog`.
- **Paralog compensation:** Z0109|ARID1A, with β_P = −0.37, and Z0186|EP300, with β_P = −1.09.
- **Residual sensitisation:** Z0286|FOCAD, with κ_N = 0.14.
- **E2 (ruling R5, post hoc, conservative).** Z0130|PTEN (engine class COLLAPSE, κ_N = −0.057) and Z0158|TP53
  (engine class RESIDUAL SENSITISATION, κ_N = 0.056) also meet TOST-neutral: their 90% CIs lie inside ±0.1.
  - Both are reported as "DETECTABLE, WITHIN EQUIVALENCE MARGIN (negligible)".
  - They are not headline findings.
  - They still go to lockbox confirmation on κ_N, but a confirmation cannot promote them.
- **Neutral:** 11 events. **Inconclusive:** 6 events.
- **Not evaluable.** MAP2K4 is labelled NOT EVALUABLE in the reporting layer. This follows ruling R2; the engine run
  predates the fix and H1 was not re-run. SMAD4 has no evaluable axis.

**Notes.**
- **E1:** the sign typo in the score_rule text is corrected; the engine semantics (ALIGN = +1, in place since
  `201d8bf`) are unchanged.
- **D8:** not applicable. The H1 calibration is the per-event matched null, which ran with discovery.
- **SMARCB1 sensitivity.** Adding the copy-number covariates gives κ_N = 0.025, with 90% CI [−0.038, 0.089], n 732.
  The primary is 0.018 [−0.043, 0.080], n 960. Both are NEUTRAL.
- **GPD sensitivity (R3):** all 5 discoveries are kept.
