# H3a FDR under planted lineage effects: diagnosis

**Trigger.** `tests/test_part2_e2e.py::test_fdr_h3b_null_with_confounded_lineage_effects_50` failed. It runs
`fdr_sims(50, seed0=20000, null=True, lineage_null=False)`, and the H3a mean FDP was 0.0748, above the bound Q + 0.02 = 0.07.
Only synthetic data were used (`part2.synth.simulate`). No DepMap, huMAP or test-eval file was read.

**Verdict.** There are two causes, and both are needed to explain the failure.

1. **Method bug (candidate c): the permutation re-entry was anti-conservative.** The code used a GPD tail whenever there
   were fewer than 10 permutation exceedances. It then returned min(p_GPD, p_emp), and ran 1,000 permutations. The
   registry `statistic` for H3a, H3b and H2 says "10,000 permutations; GPD tail when p_perm < 1e-4". So the code used the
   GPD for p ≲ 1e-2, which is exactly the range of BH thresholds that decides rejections. The GPD extrapolation is a
   noisy estimate of p. By Jensen's inequality, P(p̂ ≤ x) = x·E[p/p̂] > x. Fixed.
2. **Mis-defined error metric (a cousin of candidate a).** The test bounded the *pooled* FDP over all H3a discoveries.
   The pre-registered procedure (DESIGN v3.1 §2.5 and `docs/registry/H3a_lineage_atlas.yaml`: Benjamini–Bogomolov,
   families = contexts, Simes + BH selection, within-family BH at q·R/M) controls something else: the expected
   *average, over selected families, of the within-family FDP*. It does not control the pooled FDP. Even with exact
   p-values, BB here gives a pooled FDP of 0.066 and a BB error rate of 0.046.

The other two candidates are ruled out:

- **(b) HC3/df is not the cause.** Null-unit analytic p-values are calibrated in every lineage.
- **(d) Monte Carlo noise is not the cause.** Over 1,000 sims the old code's pooled FDP is 0.086 ± 0.003, so the first
  50 seeds (0.075) were actually lucky.

The truth definition (candidate a proper) was checked and is correct. A false positive is counted only on `null_units`,
which are complexes that carry no planted signal. For those, the one-vs-rest difference in expected score is exactly 0,
because no lineage term enters their genes. Signal-carrying units that leak through the "rest" mean are in `truth["H3a"]`
or `other`, and are never counted as false. In this scenario `other` = 0.

## Evidence

All runs use scenario `null=True, lineage_null=False, effect=0.3`, seeds 20000–20999 (n = 1000, superset of the test's
50). Values are shown as mean ± SE.

### 1. Old code (1,000 permutations, GPD when k < 10, min(p_GPD, p_emp))

| quantity | value |
|---|---|
| pooled FDP, first 50 seeds (the test) | 0.0748 ± 0.0114 |
| pooled FDP, 1000 sims (permutation re-entry = `discovery`) | **0.0863 ± 0.0027** |
| pooled FDP, analytic BB only (`reject`, before permutation) | 0.0663 ± 0.0024 |
| BB error rate (average over selected families), analytic | 0.0457 ± 0.0018 |
| BB error rate, after permutation re-entry | **0.0582 ± 0.0020** (> q) |
| E[false discoveries] per sim, analytic → permutation | 0.815 → 1.040 |
| false discoveries that the analytic BB had **not** rejected | 32% |
| false discoveries: median p_perm / p_analytic | 0.45 (81% have p_perm < p) |
| screened null tests: median p_perm / p | 0.76 (65% lower) |

**Analytic p-values are calibrated.** Over the 688,000 null-unit tests, P(p < 0.05, 0.01, 0.001) = 0.0481, 0.0094,
0.00095. Every lineage (n_context ≈ 84, Satterthwaite df ≈ 109) has P(p < 0.001) between 0.00084 and 0.00108. The KS
rejection comes from slight conservativeness of HC3, not from excess small p-values. This rules out (b).

**False discoveries are spread evenly over the 8 lineages (91–163 each).** No context is aberrant, which argues against
a small-group artefact.

**Direct size check of the permutation p rule.** The script `gpd_validity2.py` uses t_null ~ t(110) with the true p
integrated over U(0, 0.03), and 10,000 draws per rule. Each cell is P(p̂ ≤ x)/x, with Monte Carlo SE ≈ 8% at x = 5e-4:

| B | rule | x = 5e-4 | 1e-3 | 2e-3 | 5e-3 |
|---|---|---|---|---|---|
| 1,000 | old (GPD at k < 10, min with p_emp) | **1.82** | **1.52** | **1.30** | **1.15** |
| 1,000 | exact (k+1)/(B+1) | 0 (floor 1/1001) | 1.01 | 1.03 | 1.03 |
| 10,000 | old | 1.13 | 1.00 | 1.00 | 1.00 |
| 10,000 | registered (GPD only when p_perm < 1e-4) | 1.04 | 1.04 | 0.99 | 0.99 |
| 10,000 | exact | 1.08 | 0.98 | 1.01 | 1.00 |

With m ≈ 90 units per family, the BB thresholds are q·R/M·k/m ≈ 5e-4·k. The old rule inflates the size there by 30–80%.

### 2. Why the pooled FDP exceeds q even with exact p-values

The numbers below are per selected family, analytic stage, from 7,180 families. Here a = q·(R/M)·(m0/m) ≈ 0.043:

| true discoveries T_f | families | mean false V_f | predicted a·(T_f+1) | mean family FDP |
|---|---|---|---|---|
| 0 | 32 | 1.13 | (family selected only by nulls) | 1.00 |
| 1 | 4358 | 0.088 | 0.086 | 0.042 |
| 2 | 2488 | 0.135 | 0.129 | 0.043 |
| 3 | 292 | 0.192 | 0.172 | 0.046 |

Within each family, BH behaves exactly as theory says: family FDP ≈ a whatever T_f is. But a BH family with T_f true
rejections also admits about a·(T_f + 1) false ones. Here the 7.2 selected families carry only 10.3 true discoveries
between them, about 1.4 each. Pooled, that gives about a·(1 + F/ΣT) = 0.073, against an observed 0.066. The
BB-controlled quantity is 0.046 ≤ q.

So a pooled bound of q + 0.02 tests a property the registered procedure does not have. It would fail about half the
time at 50 sims even with an exact method, because the pooled FDP after the fix is 0.068. Pooled atlas-wide FDR control
would need BH over all (context, complex) pairs. That is a design change to a frozen registry entry, and it was **not**
made here.

## Fix

1. `src/part2/stats.py::perm_pvalue` now implements the registered rule. It returns the exact (k+1)/(B+1) unless that
   is below `gpd_below` (default 1e-4, the registered value; it also needs < 10 exceedances and B ≥ 1000). Only beyond
   the permutation resolution is the GPD tail used. This affects every caller (H3a, H3b, H2, H1, lockbox
   confirmation), and all of them register the same rule. New unit test:
   `tests/test_part2_stats.py::test_gpd_only_below_registered_resolution`.
2. `tests/part2_synth_pipeline.py` now uses the registered `N_PERM = 10_000` for discovery, the full pipeline and
   `fdr_sims`, replacing 1,000 or 2,000. With B = 1,000 under the registered rule, p cannot go below 1/1001. That is
   above every BB threshold 0.05·k/90 for k = 1, so planted power falls to 0 (checked: seed 20003 has power 0 at
   B = 1,000 and 1 at B = 10,000).
3. `_metrics` now also reports `fdp_bb`, the BB error rate: the mean over families with ≥ 1 discovery of V_f/R_f.
   Every BB-selected family has ≥ 1 rejection, since a Simes p ≤ q·R/M implies that BH at q·R/M rejects. The H3a
   assertion now bounds `fdp_bb` at the **same** Q + 0.02, plus power ≥ 0.95. The pooled FDP is still reported in the
   `fdr_sims` output and in the `__main__` JSON. The global-null test is unchanged: under the complete null, the pooled
   FDP equals the BB error rate, which equals 1{any discovery}. The H3b assertion is unchanged.

The bound was not loosened, and the number of sims and the seeds were not changed.

## Validation (fixed code, 10,000 permutations)

| scenario | n sims | metric | value |
|---|---|---|---|
| planted lineage, null genotype (the failing test's scenario) | 1000 | **H3a BB error rate** | **0.0463 ± 0.0018** |
| | | H3a pooled FDP (reported) | 0.0682 ± 0.0024, analytic-only 0.0663 |
| | | H3a planted power | 1.000 |
| | | permutation re-entry: median p_perm/p on screened null tests | 1.008 (49% lower) |
| | | H3b mean FDP | 0.024 ± 0.005 |
| global null (`null=True, lineage_null=True`), seeds 10000–10299 | 300 | P(any H3a / H3b discovery) | 0.020 ± 0.008 / 0.023 ± 0.009 |
| planted genotype + lineage, effect 0.15, seeds 30000–30199 | 200 | H3a / H3b BB error rate | 0.047 ± 0.004 / 0.046 ± 0.005 |
| | | H3a / H3b pooled FDP | 0.063 ± 0.005 / 0.069 ± 0.007 |
| | | H3a / H3b planted power | 0.984 / 0.909 |

On the test's 50 seeds, the fixed code gives a pooled FDP of 0.060. The full 1,000-sim BB error rate is 0.046 ± 0.002.
The probability that the 50-sim mean exceeds 0.07 is about 0.002 (SE ≈ 0.008 at n = 50, z ≈ 2.9).

**Open item for the authors.** If the paper's atlas claim is meant to be "FDR ≤ 5% over all reported H3a pairs", the
registered BB procedure does not deliver it: it gives about 6.6–6.8% pooled in this design. The claim should be worded
as BB's selective (average-over-selected-lineages) FDR, or the registry should be amended before FZ-X.

Diagnostic scripts (session scratchpad, not committed): `run_sims.py` / `run_sims2.py` (save full H3a/H3b tables per
sim), `analyze.py`, `decomp.py`, `gpd_validity2.py`.

## Side effect: `test_planted_effects_recovered` (single seed, left unchanged, now fails)

The e2e fixture (`synth.simulate(1)`, now run with the registered 10,000 permutations) asserts `false_pos <= 1` for both
H3a and H3b. Under the fixed code, H3a for seed 1 has 2 false positives out of 19 discoveries: Breast:C0010
(p = 0.0011) and Skin:C0073 (p = 0.0017). Both are **also rejected by the analytic BB before any permutation** (q_adj
0.032 / 0.039). The fixed permutation re-entry now reproduces the analytic BB exactly (19 = 19).

The old code passed only by chance. Its noisy B = 1,000 permutation p-values knocked the whole Breast family out of
selection, which dropped C0010 together with 4 true discoveries (14 discoveries instead of 19).

Across 140 planted-effect seeds (effect 0.3, fixed code), P(H3a false_pos ≤ 1) = 0.79 and P(H3b ≤ 1) = 0.88. P(both)
is only 0.69, with E[false_pos] = 0.91 and 0.54. Planted power is 1.0 in every seed, and the H3a BB error rate is 0.042.

So this assertion is a single-realisation check that holds in only about 69% of seeds. It is not a property that BB
guarantees. It was **not** weakened here, and the decision is left to the authors. Options are to bound a
per-realisation quantity that is justified by the sims, or to rely on the FDR simulation tests for error control and
keep the power and H1 checks for this seed.

## PI decision (2026-10-28)

The per-seed assertion `false_pos <= 1` for H3a/H3b in `test_planted_effects_recovered` is removed. It was an
arbitrary single-realisation threshold, not a property of the registered BB procedure. Error control stays asserted,
over many simulations, by `test_fdr_global_null_50` and `test_fdr_h3b_null_with_confounded_lineage_effects_50` (BB rate
≤ Q + 0.02) and by the 1000-sim evidence above. Power (= 1.0) and H3b true_pos (= 5) stay asserted on seed 1. Neither
the method nor the registered procedure was changed to pass a test.

Consequence for reporting: H3a/H3b claims are worded as Benjamini–Bogomolov selective FDR (average over selected
context families) at q = 0.05. The pooled FDP over all reported pairs is not claimed. In this design it runs at about
6.8% (synthetic), and the manuscript states this.
