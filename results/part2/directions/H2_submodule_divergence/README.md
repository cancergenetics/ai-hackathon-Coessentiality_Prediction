# H2_submodule_divergence: A-side discovery (single look, EXPLORATORY)

**Commits.** Intent `5ea25ca`, then results `1274e40`. H2 is exploratory because neither of its kill-rule controls
(PC21, PC22) is evaluable.

**Discovery.**
- 210 eligible (context, parent) families from the committed H3a and H3b discoveries.
- m = 452 tests; 134 passed the screen and were re-tested by permutation.
- 109 stage-wise discoveries, 97 of them outside locus clusters.

**D8: the registered calibrations were implemented and run after the discovery run.** PI ruling R4; commit
`1a2a445`; result `calibration_D8.json`. They are reported only; no response is registered.
- **Random member splits of the same parents.** These keep the size of each sub-module (20 replicates per tested
  pair).
  - Fraction p < 0.05: 0.375 for the splits, 0.460 for the real tests.
  - The stage-wise machinery rejects on average 78.3 tests per replicate (maximum 90), against 109 real discoveries.
  - **Most H2 "divergence" is reproduced by arbitrary splits of the same parent.** A sub-module-specific claim needs
    the real sub-module to beat its own random splits. This is a major caveat on H2.
- **H3 pseudo-contexts through the stage-wise machinery.** Fraction p < 0.05 = 0.050 with 6 rejections, against
  roughly 0 expected.

**GPD sensitivity (R3).** 107 of the 109 discoveries are kept.
