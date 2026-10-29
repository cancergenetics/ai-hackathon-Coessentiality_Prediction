# H1S_silencing_GI: A-side discovery (single look, EXPLORATORY; never headlined)

**Commits.**
- `BLOCKED.md` records the engine crash found by the first dry check.
- R2 fix `b3e8043`; synthetic outputs were identical before and after by SHA-256.
- Intent `b865595`, then results `06c8cd1`.

**Discovery.**
- 317 events and 321 tests. m = 293 evaluable tests; the registered m is 338.
- 49 tests passed the screen and were re-tested by permutation.
- 31 BH discoveries, 26 of them outside locus clusters.

**Classes.**
- **Complex-wide axis:** 23 COLLAPSE and 1 RESIDUAL SENSITISATION. Under E2 (`label_final`), 4 of these are
  "DETECTABLE, WITHIN EQUIVALENCE MARGIN (negligible)": ERBB3, PRR15L, B3GNT5 and CST6. Also 191 NEUTRAL, 45
  INCONCLUSIVE, 13 with no matched null, and 27 NOT EVALUABLE under the R2 fix.
- **Paralog axis:** 2 cases of paralog relief (exploratory): SOX10 and BTK.

The collapse calls are dominated by lineage-restricted non-expression (B-cell, melanocyte and renal genes), the
confounding A1 anticipated.

**Other notes.**
- **D8:** not applicable. The H1S calibration is the matched null, which ran with discovery.
- **GPD sensitivity (R3):** 29 of the 31 discoveries are kept.
