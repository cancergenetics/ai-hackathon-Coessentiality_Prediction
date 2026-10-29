Dear Editor,

Please consider our manuscript, "A label-free co-essentiality score wins a pre-registered held-out test, but its advantage does not carry over to complexes or to complex-level cancer dependencies", for publication as a research article.

**What we did.** Co-essentiality methods are often tuned and judged on overlapping curated labels, and complex-level analyses of cancer dependencies often assume that a complex, rather than any similar gene set, is the unit that matters. We tested both points in a pre-registered, two-part study:
- **Part 1.** A label-free co-complex score, locked before the evaluation data existed, was evaluated once on withheld hu.MAP 3.0 test pairs.
- **Part 2.** The score was turned into a curation-free complex map on 1,007 DepMap lines. Every complex-level claim built on the map was confirmed in a single read of 338 held-out lines. <!-- results/part2/lockbox/manifest.json: n_A, n_L -->

Protocols, decision rules, calibrations and confirmation rules were frozen and committed before the outcomes they govern. Independent re-analyses reproduced every confirmatory number.

**What we found.**
1. **The method wins the test it was built for.** Average precision on held-out complexes was 0.333 against 0.212 for GLS, the strongest baseline (paired gain +0.121, 95% CI +0.079 to +0.159). The gain is concentrated in complexes with many positive pairs. <!-- results/test_eval/test_results.json: results.itt.bootstrap.cluster -->
2. **The advantage does not carry over to complexes.** Both pre-registered bridges failed. The map did not agree better with BioPlex physical interactions than a GLS map (null). Context associations discovered on it were confirmed in held-out lines less often than those discovered on the GLS map (−0.065, −0.107 to −0.025). <!-- results/part2/construction/bridge/bridge.json: primary; results/part2/lockbox_confirm/bridge_primary_vs_corrected_R7.json: primary_evaluable -->
3. **Most complex-level context dependencies are not complex-specific.** Random gene sets matched for size, essentiality and coherence are associated with lineages and drivers as often as real complexes. About half the associations confirm in held-out lines, which is what the lockbox's power predicts. Complex scores did not replicate better than their best single gene.
4. **Within complexes, known paralog compensation was recovered** (ARID1A → ARID1B, also a registered control). Most tumour-suppressor losses left their complex's dependency unchanged.
5. **What survives.** A registered rule leaves 254 strongest findings: 207 lineage and 45 driver associations and 2 within-complex calls that were confirmed in held-out lines, are robust to leaving out each primary disease or lineage, lie outside copy-number locus clusters and do not depend on a code defect described below. We lead with physically supported units whose effects depart from their lineage's essentiality trend, such as ovarian dependency on a CDH1/PAX8 unit, bone dependency on SNAPC1–5 and myeloid dependency on AP2. We label 56 lymphoid and myeloid findings that the lineage trend could explain as a possible lineage-scale artefact. Any functional reading is post hoc. <!-- results/part2/robustness/summary.json: per_direction; manuscript/tables/strongest_derived.json: possible_lineage_scale_artefact -->

**Why it matters.** The paper is a largely negative, fully disclosed test of an idea that is widely used: that better pair-level co-essentiality yields better complexes, and that complexes are the natural unit for cancer-context dependencies. Its practical message is that complex-specificity must be tested against matched gene sets, not assumed.

**What we disclose.**
- After the single held-out read, we found a defect in the confirmation code. It left 53 fits unevaluable and shifted 510 others slightly. <!-- results/part2/lockbox_confirm/README.md R7; defect_sensitive_calls_R7.tsv --> We kept the read as executed as the primary result, reported a corrected re-analysis alongside as non-confirmatory, and removed any claim from the 53 calls that the correction changes.
- Some registered calibrations were implemented only after the discovery runs.
- The registered per-hit robustness checks were computed after the held-out read, on the construction set only, against a candidate list committed first.
- Every amendment, erratum and ruling is listed with its timing relative to outcome access in a deviations table.
- The work was carried out by AI agents under a human principal investigator, with independent agent reviews at each stage.

All inputs are public DepMap (including CCLE proteomics), BioPlex, Ensembl, HGNC and hu.MAP releases. No restricted data are shared. The analysis is reproducible from these releases and the committed code, registry and result files. The manuscript has not been submitted elsewhere. We have no competing interests to declare. [Authors to confirm.]

Sincerely,

[[AUTHORS]]
