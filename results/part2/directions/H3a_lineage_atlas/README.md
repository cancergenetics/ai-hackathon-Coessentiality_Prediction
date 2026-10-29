# H3a_lineage_atlas: A-side discovery (single look)

**Commits.** Intent `228077a`, then results `b0dbd83`. The kill rule passes: 2 of 2 controls.

**Discovery.**
- m = 8,441 tests in 23 families.
- 1,666 tests passed the screen and were re-tested by permutation (10,000 permutations).
- 1,145 BB discoveries at q = 0.05, 942 of them outside locus clusters.

**Calibration run with the discovery.** Pseudo-contexts: fraction p < 0.05 = 0.049, against 0.05 expected. There
were 0 BB rejections, so the pre-declared response did not fire.

**D8: the other registered calibrations were implemented and run after the discovery run.** PI ruling R4;
`../CALIBRATION_D8_INTENT.md`; result `calibration_D8.json`, commit `fd582a6`. They were not in the engine when
discovery ran. They are reported only and change no discovery.
- **Matched random sets.** 3,500 sets covering 350 of 367 units. Fraction p < 0.05 = 0.281, the same as the real
  tests (0.281): one-vs-rest lineage contrasts move most co-dependency modules. 502 of the 1,145 discoveries exceed
  every one of their unit's matched sets.
- **Never-expressed sets.** Fraction = 0.135, against 0.05 expected, with 13 BB rejections. Lineage contrasts shift
  even non-expressed gene sets (screen-level effects). Read lineage hits with this background in mind.

**GPD sensitivity (R3).** This is a sensitivity analysis. With p_cons = max(p_final, analytic p), 1,113 of the 1,145
discoveries are kept (32 lost, 0 gained). `discoveries_report.tsv` ranks discoveries by |t| and shows p_perm as
"<1e-4" wherever the GPD tail was used.
