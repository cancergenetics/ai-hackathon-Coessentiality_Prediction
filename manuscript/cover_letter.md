Dear Editor,

Please consider our revised manuscript, "A co-essentiality score selected on training complexes wins a single held-out test but yields neither better-supported complexes nor more replicable complex-level dependencies in DepMap", for publication as a research article. A point-by-point response to the referee is enclosed.

**What we did.** Co-essentiality methods are often tuned and judged on overlapping curated labels, and complex-level analyses of cancer dependencies often assume that a complex, rather than any similar gene set, is the unit that matters. We tested both points in a pre-registered, two-part study:
- **Part 1.** An unsupervised co-complex score, whose configuration was selected on hu.MAP training complexes and which was locked before the withheld hu.MAP 3.0 test files reached the project, was evaluated once on those test pairs.
- **Part 2.** The score was turned into a complex map on 1,007 DepMap lines without reading any curated complex resource. Every complex-level claim built on the map was tested in a single read of 338 held-out lines. <!-- results/part2/lockbox/manifest.json: n_A, n_L -->

Protocols, decision rules, calibrations and confirmation rules were frozen and committed before the outcomes they govern. Independent re-analyses reproduced every confirmatory number.

**What we found.**
1. **The method wins the test it was built for.** Average precision on held-out complexes was 0.333 against 0.212 for GLS, the strongest baseline (paired gain +0.121, 95% CI +0.079 to +0.159). The gain is concentrated in complexes with many positive pairs; macro per-cluster AP does not differ. <!-- results/test_eval/test_results.json: results.itt.bootstrap.cluster -->
2. **The advantage did not translate into better-supported complexes.** The registered primary bridge (BioPlex agreement) is null. In the registered secondary bridge, context associations discovered on the map were confirmed in held-out lines less often than those discovered on a GLS map (−0.065, −0.107 to −0.025); a post hoc power calculation accounts for roughly a quarter to two fifths of that difference. <!-- results/part2/lockbox_confirm/bridge_primary_vs_corrected_R7.json: primary_evaluable; results/part2/revision/README.md R-4 -->
3. **Complex-level context associations are shared by coherent co-dependency modules.** Gene sets drawn from other co-dependency modules, matched for size, essentiality and coherence, associate with lineages as often as complexes; sets matched only for essentiality do so less often. About half the associations confirm in held-out lines, which is what the lockbox's power predicts (post hoc). Whether a complex replicates better than its best single gene is inconclusive.
4. **Within complexes**, known paralog compensation was recovered (ARID1A → ARID1B, a registered control). For 11 of 24 genetic-loss events the complex-wide dependency changed by less than ±0.1 Chronos units, the registered margin; none of the 11 would be neutral at ±0.05, and for one of them (ARID1A) the paralogue still changed.
5. **What survives.** A registered rule, applied after the held-out read on the construction set only, leaves 254 strongest findings: 207 lineage and 45 driver associations confirmed in held-out lines (or, for five, by a rare-context route that is not an independent replication) and 2 within-complex calls (ARID1A → ARID1B, confirmed; EP300 → CREBBP, which passed only the rare-context route). The examples we name were chosen by post hoc labels; for many of them a lineage-scale artefact cannot be excluded, and we report an RNAi (DEMETER2) check for each. <!-- results/part2/robustness/summary.json: per_direction -->

**Why it matters.** The paper is a largely negative, fully disclosed test of an idea that is widely used: that better pair-level co-essentiality yields better complexes, and that complexes are the natural unit for cancer-context dependencies. Its practical message is that complex-specificity should be tested against matched co-dependency modules, conditioned on discovery, not assumed.

**What we disclose.**
- After the single held-out read, we found a defect in the confirmation code. It left 53 fits unevaluable and shifted 510 others slightly. We kept the read as executed as the primary result, reported a corrected re-analysis alongside as non-confirmatory, and removed any claim from the 53 calls that the correction changes. <!-- results/part2/lockbox_confirm/README.md R7 -->
- Some registered calibrations were implemented only after the discovery runs; the registered per-hit robustness checks and all reporting labels were computed after the held-out read, on the construction set only.
- The revision analyses requested by the referee are post hoc and descriptive, run on the construction set only (the second round's counts on committed result tables only). An independent agent recomputed their numbers; the check was not blind and did not rerun the random atlases (Supplementary Note 17).
- Every amendment, erratum and ruling is listed with its timing relative to outcome access.
- The work was carried out by AI agents under a human principal investigator, with independent agent reviews at each stage. [CORRESPONDING_AUTHOR] takes responsibility for the analysis.

All inputs are public DepMap (including CCLE proteomics), BioPlex, Ensembl, HGNC and hu.MAP releases. No restricted data are shared. The analysis is reproducible from these releases and the committed code, registry and result files (https://github.com/cancergenetics/ai-hackathon-Coessentiality_Prediction; [ARCHIVE_DOI]). The manuscript has not been submitted elsewhere. Competing interests: [COMPETING_INTERESTS].

Sincerely,

[AUTHORS]
