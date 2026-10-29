# Lockbox machinery before the read: the gaps that remain

> **Resolved by amendment A4 (PI rulings, before any lockbox read; no L values read).**
> - G1: BENCH secondary is NOT RUN.
> - G2: the multiverses report point estimates only.
> - G3: the six H1S N-axis rare-route calls are NOT EVALUABLE in L.
> - The literal readings in section 2 are accepted.
>
> With these rulings the driver's preflight passes without reading any L value.

Date: 2026-10-29.

**No lockbox read.** The ledger does not exist, no target list has been written, and no L value has been loaded.

This note follows the PI rulings on `BLOCKED.md`. The machinery is in place:
- `scripts/48_lockbox_confirm.py` (the driver);
- `src/part2/bench.py`: BENCH secondary, the BRIDGE-L unweighted rate, and a `bridge_rate` bug fix;
- the BRIDGE-L GLS-map discovery runner (`gls_discovery` in the driver);
- amendment A3 (`scripts/49_amend_A3.py`).

It is tested on synthetic data in `tests/test_part2_lockbox_driver.py`.

The driver's preflight runs on the A side only and refuses to start while any item in section 1 is open, so nothing is
read until every gap is ruled on. Run on real A data with the A3 registry, the preflight reports:
- 28 confirmatory entries in registry order;
- 1,295 BENCH primary pairs and 15,047 secondary pairs;
- GLS inputs that are frozen.

## 1. The registry is silent: rulings needed (the driver refuses until they exist)

| # | Step | What is open | Registry field the driver reads |
|---|---|---|---|
| G1 | BENCH secondary, "A–L correlation … vs random gene sets matched on size, essentiality and coherence (seed 7100)" | Five choices are open (see below). | `estimand_secondary: {correlation, R, pair_set, matching, summary}` in `BENCH_complex_vs_gene` |
| G2 | BENCH and BRIDGE-L multiverse ("reported, never selected") | Whether the multiverse gets the primary's bootstrap CIs or point estimates only. Bootstrap CIs for the three BENCH variants add about 11 h. | `multiverse.ci`: `bootstrap as primary` or `point estimate only`, in both entries |
| G3 | Rare-context route for H1S, "DEMETER2 on A models same sign" | Six non-lockbox-testable H1S calls on the N axis (Z0018\|IGSF11, Z0018\|SOX10, Z0035\|PAX8, Z0084\|NFATC2, Z0149\|SPDEF, Z0316\|IGF1R). Their aligned score κ_N needs baselines b_j and the N′ set: either the CRISPR A parameters or a DEMETER2 A refit. The other rare-route hits need no fitted parameters: H1 EP300 (P axis: a mean), H3b EWSR1-fusion (11 hits, S_C), and controls PC15, PC16 and PC22. | `lockbox.rare_route_fitted_parameters`: `A CRISPR parameters` or `DEMETER2 A refit`, in `H1S_silencing_GI` |

The five open choices for G1:
- the correlation statistic (Pearson or Spearman);
- R, the number of matched sets per complex;
- the matching rule. The only implementation is `calib.matched_random_sets`, the D8 operationalisation (±0.1 essentiality, ±1.0 z coherence, sets drawn from the members of other frozen complexes, tolerances doubled once);
- the pair set: primary, or "all evaluable (C, c) pairs of those contexts";
- the summary of r_complex against the R replicate correlations. The function returns r_complex, the R values r_random, the count of replicates with r_random ≥ r_complex, and quantiles.

## 2. Interpretations taken literally from the registry text (flagged for the review, not refused)

- **Breadth:** taken over all of C's members, as the registry says "C's members". g* is chosen among the members after
  the context's removals, as the registry says "after the same removals".
- **g\* by largest |beta| (multiverse):** the member with the largest beta in the complex's A direction (ties go to t).
  This mirrors the primary rule as implemented in `bench.choose_gstar`: largest t in the A direction.
- **"Number of distinct complexes among the g\*":** reported as three counts: distinct g\* genes, distinct frozen coarse
  complexes that contain a g\*, and distinct complexes in the pairs. On real A data these are 920, 353 and 353.
- **BENCH secondary pair set "those contexts":** read as BENCH's 41 lockbox-testable H3a/H3b contexts. The single BENCH
  read covers the members of every complex in the primary and secondary sets. The matched-set genes are members of
  frozen complexes, so the read covers them too.
- **BRIDGE-L GLS pass:**
  - It uses the H3a/H3b machinery: removals, per-context covariates, 12/12 evaluability, a screen at q < 0.1,
    10,000 Freedman-Lane permutations, and seed 4500 for both calls.
  - The BRIDGE-L entry has no `calibration` field or pseudo-context seed. So no pseudo-context calibration is run and
    no calibration response applies.
  - The pass runs after every confirmatory entry is closed, as `run_after` requires, and is committed before its read.
- **BRIDGE-L rates:**
  - Cluster size is the number of frozen members, which equals `clusters.tsv` size for both maps (checked).
  - The zavg confirmation status comes from the H3a/H3b confirmation tables; there is no new zavg read.
  - Hits on the rare-context route are excluded, because the registry says "lockbox-testable contexts only".
  - The "matched gene coverage subset" is `validation.match_coverage` (FZ-C decision D5) applied to the level-1 clusters.
- **Testable in the frozen list but under 5/5 L models after complete cases:** marked not confirmed and flagged. The
  context is lockbox-testable, so it does not take the rare-context route.
- **Hit order and permutation seed:** hits are ordered by (context, unit). The permutation seed is A3's
  `lockbox_permutation` + hit index.

## 3. Bug found and fixed before the read

`bench.bridge_rate` selected hits with `h.map == m`. In pandas ≥ 2.1, `DataFrame.map` is a method, so the
BRIDGE-L estimate would have failed *after* its single read. It now uses `h["map"]`, with a regression test. No other
engine code changed.

## 4. Projected runtime of the single read

Measured on synthetic data: `bench_delta` takes 5.4 ms per pair per resample. The A side of BENCH takes 103 s on real
A data.

| Step | Time |
|---|---|
| Confirmations: 28 reads; about 800 L permutations at about 0.3 s each; rare routes on A | about 25 min |
| BRIDGE-L GLS discovery on A | about 11 min |
| BENCH primary bootstrap (B 2000) | about 3.9 h |
| BENCH secondary (R = 10, secondary pair set) | about 15 min |
| BENCH multiverse | about 11 h with bootstrap CIs; about 1 min with point estimates only |
| BRIDGE-L estimate | under 1 min |

Total: about 5 h with point-estimate multiverse, or about 16 h with bootstrapped multiverse. Peak RSS is about
3–4 GB. The real A-side preflight peaked at 2.8 GB.
