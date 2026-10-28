# Round 1 review: statistics

Reviewer role: biostatistics. Objects: DESIGN.md v1 and docs/TEST_PROTOCOL.md v1.0 draft.

## What I ran
No hu.MAP test file was opened, and no association between a complex and a context was computed.
- **Synthetic simulations** (`/private/tmp/claude-502/-Users-metinyazar--claude-spares-spare-gqTg9n/5a84c62a-d8a2-48ad-b26c-5dd43a72db34/scratchpad/sim_stats.py`, seed 0, 63 s):
  - type-I error of the permutation tests;
  - power of the full-data test vs the split-half rule;
  - power of the complex mean vs the best single gene.
- **Label-free counts:** sizes of mutation contexts, mutation burden, duplicated PatientIDs, and the scale of the CN file.

## Part 1 (TEST_PROTOCOL.md)

The protocol is sound:
- intention-to-score (ITT) set with one shared floor score;
- paired cluster bootstrap, bit-equal to Part 1;
- thresholds and calibrators fitted on train only;
- intersection–union verdict with a fixed testing sequence.

One change blocks the freeze.

**M1. The disclosure in §0 and §9 is out of date.**
- *Problem:* §0 says de-duplication counts, overlap counts and mapping rates were "not computed", and §9.5 says the audit is "pending". docs/part1_audit.md §4 has since computed:
  - unique counts;
  - zero duplicates, self pairs and conflicts;
  - zero overlap with train pairs;
  - protein overlap with train;
  - C1/C2/C3 counts of the positives;
  - scorable counts by label.
- *Change:*
  - List these quantities as seen before the freeze, and mark the audit as PASS.
  - State that S2 equals the primary set.
  - Add the audit's C1/C2/C3 strata, which are defined on train-*positive* proteins. S3 uses a different definition (train pairs), so report the strata alongside S3.
- *Why:* the freeze note must be true, and the strata must match the published audit.

**S-P1 (SHOULD). Report uncertainty for the criterion-3 gate.**
- *Problem:* criterion 3 is a gate on point estimates. In dry run C it passed by only about 8%, so it may decide the verdict.
- *Change:* call it "non-degradation", and report bootstrap CIs for the macro-AP ratio and for the ratio of AP without the largest component. The CIs are descriptive only.

**C-P1 (COULD). Pre-declare a borderline result and check coverage.**
- *Problem:* the coverage unit test accepts 0.85 with an n_eff of about 26 clusters.
- *Change:*
  - Pre-declare "borderline" when the lower bound of the ΔAP CI is within the Monte-Carlo spread (≤0.011) of 0. The verdict still uses seed 0.
  - Add a synthetic coverage check that uses the train cluster-size distribution.

## Part 2 (DESIGN §2)

### MUST

**M2. Mutation contexts are mostly a proxy for mutation burden.**
- *Problem:*
  - 2,886 genes have a damaging mutation in ≥8 models.
  - In 79% of these contexts, most mutant lines fall in the top decile of mutation burden. The median share is 0.62.
  - The model in §2.4 (lineage + library) would attach every hypermutator dependency, for example WRN, to hundreds of passenger contexts.
- *Change:*
  - Add log(burden) and MSIScore as covariates to every M and D3 model.
  - Add ploidy and WGD to every CN model.
  - Report results with and without these covariates.

**M3. The permutation design is invalid, and it is too coarse for BH.**
- *Problem 1:* the design permutes a pooled-variance OLS t.
  - Simulation: 8 lines vs 1,000, t(5) noise, variance ratio 3.
  - Rejection at α = 0.05 is **0.21** (0.10 at α = 0.01).
  - A Welch-type studentised permutation gives 0.04.
- *Problem 2:* with about 10⁶ tests per family, BH needs p ≈ 10⁻⁷. 10⁴ permutations cannot resolve that.
- *Change:*
  - Use an HC3 or Welch t.
  - Use Freedman–Lane residual permutation within lineage strata, with patient blocks.
  - For screening, use robust-t p-values with Satterthwaite df, calibrated on permutation and negative-control QQ plots.
  - For contexts with n1 < 20, confirm by exact within-stratum permutation, with a GPD tail approximation for p < 10⁻⁴.

**M4. The FDR families are not defined.**
- *Problem:*
  - m is not stated for any family; the M family alone could reach about 3×10⁶ tests.
  - Neighbouring CN genes give near-identical context vectors.
- *Change:* enumerate every family and its m before any test.
  - Collapse CN contexts whose vectors have Jaccard ≥ 0.9 into segments.
  - Confirmatory M family: disruption of complex members (D3), hotspot genes, and the pre-registered control genes.
  - All other mutations: a separate exploratory family.
  - Sub-modules: tested only within complexes that pass, stage-wise.

**M5. Forking paths, and selection based on outcomes.**
- *Problem:*
  - §2.6 selects the pipeline variant (BioPlex, γ, score) by the replication rate of context associations, which tunes on the reported outcome.
  - Three score definitions, two γ levels, masking and covariate sets multiply the specifications further.
- *Change:*
  - Select by the label-free criteria of §2.2 only.
  - Freeze one primary specification in DESIGN v2: mean Chronos, coarse level, fixed membership, primary covariates.
  - For every hit, report a multiverse table: the fraction of specifications with the same sign and p < 0.05.
  - Never report the best specification for a hit.

**M6. Self-association.**
- *Problem:* in D3, D5 and EXP, a context gene that is itself a member of the complex lowers its own dependency, which trivially drives the complex score.
- *Change:* compute the score without the context gene(s).

### SHOULD

**S1. The split-half rule.**
- *Problem:*
  - Power at d = 1.5 (discovery at α = 2×10⁻⁴ on half the data, validation at one-sided 0.025):

    | n1 | power |
    |---|---|
    | 8 | 0.20 |
    | 12 | 0.45 |
    | 20 | 0.83 |

  - "Median over 10 splits passes" is undefined for a binary outcome.
  - Random halves replicate confounders rather than testing them.
- *Change:*
  - Define the rule as "passes in ≥6 of 10 splits, with BH recomputed within the half".
  - Make it a gate only for n1 ≥ 20. Below that, report the observed replication rate next to the expected rate.
  - For non-lineage contexts, add cross-lineage consistency: same sign in ≥2 lineages that each have ≥3 context lines.

**S2. Minimum context size.**
- *Problem:* at α = 10⁻⁶, power is 0.25 (n1 = 8, d = 1.5), 0.77 (n1 = 8, d = 2) and 0.61 (n1 = 12, d = 1.5).
- *Change:*
  - Confirmatory analyses need n1 ≥ 12. Contexts with 8–11 lines are exploratory, except pre-registered controls.
  - State the minimum detectable effect for each family.

**S3. Pseudo-replication.**
- *Problem:* 104 screened models come from 48 multi-model patients.
- *Change:* split, permute and bootstrap by PatientID.

**S4. The complex vs gene comparison is not like for like.**
- *Problem:* counts of hits depend on family size (about 1,000 vs 18,000 tests) and on redundancy.
- *Change:*
  - Use the same gene universe, contexts, covariates and family rules for both.
  - Primary comparison, threshold-free: correlation of the half-A and half-B effect estimates.
  - Map gene hits to unique complexes.
  - Null: random sets matched on size *and* on mean within-set correlation.
  - Choose the "best subunit" on half A and test it on half B.
- *Simulation C:* either level can win, so pre-register which you expect.
  - The complex wins when all members shift (k = 5, ρ = 0.3, d = 1): 0.56 vs ≤0.06.
  - It loses when 20% of members carry the effect (k = 10, d = 1.5): 0.00 vs 0.59.
  - The gain shrinks as residual ρ rises, and ρ is high by construction.

**S5. Screen quality.**
- *Problem:* pan-essential strength correlates with NNMD (ρ = 0.38).
- *Change:* for pan-essential complexes, add a mandatory sensitivity analysis that scores relative to the line's median over a disjoint reference set of pan-essential genes.

**S6. Winner's curse.**
- *Change:* report headline effect sizes from the held-out halves (cross-fitted), with CIs.

**S7. D6 coherence.**
- *Problem:* correlation falls with restricted range.
- *Change:*
  - Fit a one-factor model per lineage, with shared covariance and idiosyncratic variance.
  - Call a complex "fragmented" only if the shared covariance falls while member variances do not.
  - Keep D6 exploratory, and check its power on synthetic data first.

**S8. Masking and LOSO.**
- *Problem:* per-model masking of non-expressed members changes the score's composition by lineage.
- *Change:*
  - Make fixed membership the primary score, and masking a sensitivity analysis.
  - Add breadth (the fraction of members with the same sign) and a leave-top-2-out check to LOSO.

### COULD

**C1. The CN scale.**
- *Problem:* OmicsCNGeneWGS centres at about 1.03, with 0.1% and 99.9% quantiles of 0.22 and 3.79. That fits a ratio or log2(ratio + 1), and neither matches the log2 thresholds in §2.1.
- *Change:* resolve the scale from the release notes before setting the thresholds.
