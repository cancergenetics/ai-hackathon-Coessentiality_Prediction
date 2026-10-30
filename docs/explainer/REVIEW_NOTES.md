# Review notes: discrepancies, open questions, uncertain choices

These are notes for the human reviewer of the explainer. Each item was found while checking the documents against the code, the committed results and the session transcript. Where they disagreed, the explainer follows the code and the results. Nothing listed here changes a registered result.

## Discrepancies between documents and code or results

Most important first.

1. **"LOSO" changed meaning, so no member-level robustness check was ever run.**
   - DESIGN v2 (`061bd1d`) defined LOSO as leave-one-*subunit*-out: does the complex effect survive dropping any one member?
   - The registry and `src/part2/robust.py` implement leave-one-*lineage* (or primary-disease)-out.
   - `results/part2/robustness/INTENT.md` says only "early reviews" used the subunit meaning. In fact DESIGN v2 adopted it.
   - This bears directly on the open question of whether the complex, or one or two members, carries the effect. In the worked example, 4 of 5 SNAPC subunits carry it.
2. **Persistence rule.**
   - DESIGN §2.3 and the `consensus.py` docstring describe Jaccard-based support with greedy placement.
   - The frozen build used *containment* support and *nested* placement (`results/part2/construction/DECISIONS.md` D1, D2; `scripts/36_build_maps.py`).
   - Consequence: the `stability` column in `clusters.tsv` reads 0.347 for SNAPC, although the five genes were together in 100% of runs.
3. **Matched "random" gene sets.** DESIGN calls them random gene sets matched for size, essentiality and coherence. They are drawn from members of other frozen complexes. Peer review caught this and the manuscript was corrected, but DESIGN.md still carries the old wording.
4. **Permutation scheme for lineage contexts.** DESIGN says permutations run within patient blocks within lineage. For lineage contexts they run across lineages, as the registry entry states.
5. **Generalised-Pareto tail p-values.**
   - The tail is fitted only when *no* permutation reaches the observed statistic, so every p below 10⁻⁴ is an extrapolation.
   - Worked example: 1.4 × 10⁻⁶, against an empirical floor of 1/10,001.
   - DESIGN's "GPD tail < 10⁻⁴" hides this.
6. **"Each finding is confirmed once on held-out lines."** 6 of the 254 strongest findings passed only the pre-registered rare-context route (A lines plus DEMETER2 on A), not a lockbox test. The explainer says so.
7. **Lockbox firewall level.**
   - `assoc.load_A_effects` loads the whole cached gene-effect matrix, lockbox rows included, into memory and returns only A rows.
   - No lockbox value is returned or used, but the firewall holds at the code-interface level, not the data level.
8. **Who found the calibration gap.**
   - REPORT §8.3 credits "review 2 and the PI".
   - The transcript shows the discovery agent reported "registered calibrations not implemented" at 00:39 UTC, before review 2 was launched (01:18).
9. **"Label-free".**
   - The pair score runs without labels, but its configuration was chosen on hu.MAP training complexes.
   - Peer review had the term withdrawn from the manuscript. DESIGN.md and your brief for this explainer still use it.
   - The explainer says "built without curated complexes" and states the caveat.
10. **Two different "53"s.** There were 53 failed lockbox fits (50 from the defect, 3 genuine) and 53 defect-sensitive calls (50 refitted plus 3 that flipped). The numbers coincide but count different things. The PI's chat once said 4 calls flipped; the records say 3.
11. **Within-complex multiplicity.** False-discovery control used the 24 evaluable events rather than the registered 28. PTEN's discovery depends on this: q = 0.044, against 0.051 with 28. This is disclosed, but easy to miss.
12. **Script 62 runtime.** REPORT §7.5 says the runtime of the last revision script was not logged; tree note D23 records 1.6 s from a PI rerun. REPORT was not updated.
13. **z0 choice is a tie-break.** Split-half recovery is 0.3407, 0.3381 and 0.3428 for z0 = 4, 5 and 6 (SD 0.005–0.019). "Selected z0 = 6" therefore reflects the tie-break rule (prefer the sparser graph), not a measured optimum.
14. **Freedman & Lane (1983)** is not in LITERATURE.md; the method is cited only through Winkler et al. (2014). The `stats.py` docstring cites Winkler 2015, which LITERATURE.md does not list.
15. **Minor wording.**
    - DESIGN calls the 20 resampled construction runs "bootstraps"; they are 80% subsamples without replacement.
    - DESIGN gives S ≈ 6.2k genes; the actual count is 6,022.
    - The complex score is described as a "mean"; it is the mean of effects clipped to ±3, missing if fewer than half the members are observed.

## Open questions

- What explains the part of the GLS map's held-out-line advantage that power does not explain (−0.042)? Why does the power model under-predict the GLS map more (by 0.097 against 0.055)?
- How many lineage findings survive once statistics are rescaled by each lineage's never-expressed false-positive rate? This rescaling was requested in peer review and not done.
- Is the complex, or one or two members, the operative unit? BENCH is inconclusive and the subunit check was never run (item 1).
- Would a registered DEMETER2 replication in lockbox lines confirm the 254 findings?
- How much of the pair-level gain comes from preprocessing (principal-component removal, whitening) rather than from the score's form? This can no longer be tested on this test set.
- How independent are reviews by the same model family? Every reviewer, including the peer reviewer, was Claude Opus or Sonnet 5.5.

## Places where I was unsure what mattered

- **Worked-example choice.** SNAPC in Bone was chosen as the smallest, cleanest complex: 5 genes, all subunits of one known machine, lockbox-confirmed and supported by DEMETER2. But Bone is a "broad shift" context: never-expressed sets reach p < 0.05 in 19% of Bone tests. Myeloid × AP2 has a cleaner null (1%) but a 12-gene unit.
- **What counts as "the method".** The explainer treats the Part 2 pipeline (steps 1–6) as the method and the Part 1 pair score as its input. The held-out test of that input is the project's one clean positive result, so it gets its own section rather than being a step.
- **Turning points.** Eight were kept. Human input mainly set the brief and supplied data, so it appears in the research-process section rather than as a turning point; no human message changed a scientific decision.
- **Engineering insights.** Several (the waiter loop, stale agent bases, commit volume) are about running agents rather than the science. They are kept in category 3 because they cost hours.
