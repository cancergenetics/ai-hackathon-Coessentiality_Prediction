# Registered calibrations run after the discovery runs: intent (D8; PI ruling R4)

**Disclosure D8.** The discovery runs (H3a `b0dbd83`, H3b `5e1c1d1`, H2 `1274e40`) ran the pseudo-context calibration
only. The other negative controls named in the entries' `calibration` fields had no engine implementation and were
not run. That is a deviation. Under PI ruling R4 they are implemented (`src/part2/calib.py`, `scripts/46_run_calibrations.py`,
tests `tests/test_part2_calib.py`; commit `44f30caa4f3ffdcc9ad7e982be33372778ba71f4`) and run **once each on A, after
the discovery outcomes were seen**. Each pre-declared response is applied mechanically. No lockbox read.

Commands (one run each, in this order; stdout/stderr to `<ID>/run_log_calibration.txt`, `/usr/bin/time -l`):

```
export PART1_DIR="<PART1_DIR>"; export PYTHONHASHSEED=0
<main repo>/.venv/bin/python scripts/46_run_calibrations.py --id H3a_lineage_atlas
<main repo>/.venv/bin/python scripts/46_run_calibrations.py --id H3b_driver_family
<main repo>/.venv/bin/python scripts/46_run_calibrations.py --id H2_submodule_divergence
```

## Specification (registry text → operational choices fixed here, before running)

Seeds come from each entry's `calibration.seed` (H3a 4600, H3b 4700, H2 4800); sub-streams are
`numpy.default_rng([seed, …])` as listed.

1. **Matched random gene sets (H3a, H3b).** For each of the 367 frozen coarse units, R = 10 random sets with the same
   size, matched on essentiality (mean clipped A effect, ±0.1) and coherence (mean pairwise zavg_BEF on A from the FZ-C
   A block, ±1.0 z). Sets are drawn from the members of one other frozen coarse unit with ≥ k members (else the union
   of two). If fewer than R sets are accepted in 20,000 attempts, both tolerances are doubled once, else the unit is
   'unavailable'. This is the registered H1 matched-null rule, with R = 10 instead of K = 1,000. Sub-stream
   `[seed, unit index]`. Each set is scored per context with the context's removals (primary score, ≥ 2 members) and
   tested with the context's covariates by the same HC3 / Satterthwaite machinery (≥ 12 / ≥ 12).
   Reported:
   - the fraction of random-set tests with p < 0.05, per context and pooled, against the real tests' fraction;
   - BB rejections over the random-set tests (families = contexts);
   - for each primary discovery, whether its |t| exceeds the largest |t| among its unit's matched sets in that
     context.

   These sets carry real biology (lineage / genotype effects on other complexes), so they are a background, not a
   null. **No response is registered: report only.**
2. **Dispersed never-expressed sets (H3a, H3b).** 100 sets drawn from screened genes with log1p TPM < 1 in every
   RNA-profiled A model and known coordinates (958 profiled models; pool of 558 genes, from expression only). Set
   sizes are drawn with replacement from the coarse-unit size distribution. Members are drawn in random order, and a
   gene within 10 Mb of an already chosen member is skipped. Sub-stream `[seed, 1000000]`. Scored and tested as in 1.
   Reported: the pooled fraction p < 0.05 (expected 0.05) and BB rejections (expected ~0). **Report only.**
3. **CN segment-matched null (H3b CN contexts 9p21_loss, MYC_amp, MYCN_amp).** Pseudo-contexts are segments of
   `results/part2/contexts/segments_ge8.tsv` that meet all of these conditions:
   - the same event type (hom-del for 9p21; the registered file holds only `highamp` amplification segments, which
     are used for MYC/MYCN amp);
   - at another locus (not within 10 Mb of MTAP / MYC / MYCN);
   - an A count within ±25% of the context's n_A (148 / 174 / 28).

   If more than 100 qualify, 100 are taken with sub-stream `[seed, 2000000]`. Label-only count: 3 / 0 / 37
   segments, so for MYC_amp the null is unavailable and is reported as such. Members within 10 Mb of the segment
   span are removed, and the test uses the context's covariates.

   **Response (registered, applied mechanically).** The segment pseudo-contexts are the CN contexts' pseudo-context
   null. If their pooled fraction p < 0.05 exceeds 0.075:
   - every real test of those CN contexts with raw p < 0.01 that was not yet permutation re-tested is re-tested
     (Freedman–Lane, 10,000, PatientID within lineage, seed 4100 + row, as `run_association`);
   - BB is re-run on the whole H3b table, giving `discovery_calibrated` (`discoveries_calibrated.tsv.gz`); the
     changes are reported.
4. **H2 random member splits.** For every tested (eligible, evaluable) (context, sub-module), R = 20 random subsets
   of the parent's members (after the context's removals) of the sub-module's size, residualised on their complement
   exactly like `h2.residualised`. Sub-stream `[4800, sub-module index, context index]`. Each replicate is run
   through the stage-wise machinery (BH within each eligible (context, parent) family at its registered level).
   Reported: the fraction p < 0.05 and the number of stage-wise rejections per replicate, against the real 109.
   **Report only.**
5. **H2: H3 pseudo-contexts through the stage-wise machinery.** These are the same pseudo-contexts the H3a / H3b
   discovery runs generated (seeds 4600 / 4700 + context index, pseudo j matched to context j mod M). Each is tested
   against the eligible sub-modules of its matched context with the H2 covariates, then BH within families at the
   level. Reported: the fraction p < 0.05 and the stage-wise rejections (expected ~0). **Report only.**

Outputs: `<ID>/calibration_D8.json`, `<ID>/calibration_D8_tests.tsv.gz`, run logs, plus a README per direction
recording D8. One commit per calibration run.
