# H3b_driver_family: A-side discovery (single look)

**Commits.** Intent `388924a`, then results `5e1c1d1`. The kill rule passes: 10 of 10 evaluable controls.

**Discovery.**
- m = 7,340 tests in 20 families.
- 190 tests passed the screen and were re-tested by permutation.
- 161 BB discoveries, 121 of them outside locus clusters.

**Calibration run with the discovery.** Pseudo-contexts: fraction p < 0.05 = 0.045. There were 0 BB rejections,
and the response did not fire.

**D8: the other registered calibrations were implemented and run after the discovery run.** PI ruling R4; commit
`4fb83d7`; result `calibration_D8.json`. The results are below; nothing in discovery changed.
- **Matched random sets.** Fraction p < 0.05 = 0.116, against 0.120 for the real tests. 103 of the 161 discoveries
  exceed all of their matched sets.
- **Never-expressed sets.** Fraction = 0.062, against 0.05 expected, with 1 BB rejection.
- **CN segment-matched null.** This is the registered pseudo-context null for the CN contexts, and it carries the
  registered response.
  - 9p21_loss: 0.045, from 3 segments.
  - MYCN_amp: 0.061, from 37 segments.
  - MYC_amp: unavailable. No `highamp` segment at another locus has an A count within ±25% of 174.
  - Pooled: 0.060, which is at or below 0.075, so the response did not fire.

**GPD sensitivity (R3).** This is a sensitivity analysis: 130 of the 161 discoveries are kept (31 lost).
