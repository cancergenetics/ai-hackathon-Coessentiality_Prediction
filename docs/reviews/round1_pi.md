# Design review round 1 — external PI (role: pi)

No hu.MAP test file opened; no complex–context association computed.

**Feasibility check I ran** (label-free, context frequencies only; script in my scratchpad, 15 s): among the 1,361 CRISPR models (default omics entries), 2,886 genes carry a damaging mutation in ≥8 models, but 47% of all damaging calls fall in the top 5% mutation-burden models (79 models have MSIScore > 20). If those models are excluded, only 859 genes reach ≥8. Genes at ≥15 / ≥25 / ≥50 models: 799 / 215 / 30.

## Overall verdict
Part 1 is in good shape: the protocol is rigorous and the audit is clean. Part 2 is competent, but as written it is an atlas with moderate novelty. It runs eight directions, each with a large analysis budget, and it does not yet have a sentence that ties Part 1 to Part 2. The strongest publishable story is:

> A label-free co-complex predictor, validated on withheld complexes, yields a curation-free complex map that is better than the GLS-module standard. Analysing dependencies at the complex level reveals *how* cancers rewire dependencies within complexes (buffering by paralogs vs collapse of the complex, and sub-module specialisation), beyond lineage atlases already published by Wainberg 2021 and Chen 2023.

To get there: cut scope, add a bridge experiment, and make the "complex" claim defensible.

## MUST

**M1. Disclose incident D1 in the frozen protocol.**
- Problem: RESEARCH_TREE D1 says that an audit agent counted test pairs in the z ≥ 3 top-pair table, which is equivalent to precision and recall at rule (a). The tree claims this is "disclosed in TEST_PROTOCOL". It is not. §0 of the protocol says "Not computed: … overlap counts", and §9.6 discloses only the 6-line format dump.
- Change: before the freeze, add the D1 incident to the §0 table and to §9.6, stating that the counts were deleted, the PI never saw them, and rule (a) was fixed beforehand.
- Why: an undisclosed test-information touch, found later by a referee, would sink the paper's credibility on the very point (a clean held-out look) that makes Part 1 publishable.

**M2. Add the bridge experiment: run the complex pipeline on GLS scores as well as zavg_BEF.**
- Problem: nothing in Part 2 shows that the better pair predictor matters downstream. A referee will ask why Part 1 is in the paper.
- Change: run the frozen §2.2 pipeline (same grid rule, same selection criteria) on GLS genome-wide scores, the Wainberg standard, and secondarily on Pearson. Compare, label-free:
  - stability;
  - within-cluster BioPlex and proteomic enrichment;
  - downstream replication rate of context associations (split-half and RNAi).
- Pre-register this as the primary methodological claim of Part 2.
- Why: this is the single result that makes the paper one story rather than two.

**M3. Cut the scope and set kill criteria.**
- Problem: D1–D8, three score types, overlap, multiscale construction and cross-line BioPlex integration will not all finish, and a half-finished atlas is not publishable.
- Change: prioritise as follows.
  - **D1 (lineage/subtype atlas):** the foundation. It goes in the main text as one figure; its novelty is low.
  - **D3 + D4 merged into one headline direction, "Within-complex genetic interactions".** When subunit *g* is disrupted (LoF, deep deletion or silencing), do the remaining subunits become more essential (buffering, e.g. a paralog takes over) or less essential (collapse/epistasis: the complex is already non-functional)? Classify each complex as buffered or collapsing. Test complex × complex pairs only as a restricted secondary analysis. This is interpretable, mechanistic, and novel at complex scale.
  - **D2 (sub-module divergence):** the second headline.
  - **D5:** COULD, only after the headline directions are done. Its CN-correction confound inside Chronos makes it high-risk.
  - **D6, D7:** cut, because they are underpowered.
  - **D8:** fold into D1 as a supplementary confounder analysis.
  - **Kill rule:** each direction is dropped if fewer than half of its pre-registered complex-level positive controls are recovered once the pipeline is frozen.

**M4. Make "complex" a claim that is earned.**
- Problem: the clusters are co-dependency modules. Many will be pathways or lineage programmes rather than physical complexes. This is the first attack a top-journal referee will make ("these are not complexes").
- Change: give each frozen cluster an evidence tier that is fixed *before* any context test:
  - **T1** within-cluster BioPlex density above a degree-matched null;
  - **T2** CCLE proteomic co-abundance above a null;
  - **T3** co-essentiality only.
- The manuscript says "complex" only for T1/T2 and "co-dependency module" for T3. Report findings by tier.
- Keep BioPlex **out of primary construction**. Use it for validation and tiering, which is option (c) in LITERATURE §2. The 293T/HCT116 cross-integration is weak because the two networks share core complexes; demote it to a COULD sensitivity analysis.
- Why: you cannot both construct with BioPlex and validate with it. Validation is the more valuable use.

**M5. Fix the genotype context definitions.**
- Problem: hypermutated models dominate the damaging calls (see the feasibility check). In D3, "any subunit disrupted" will flag most large complexes in MSI lines, and the ≥8-model rule admits about 2,000 mostly passenger genes.
- Change:
  - Define mutation contexts and complex disruption in non-hypermutated models: exclude models above a burden or MSI cut-off fixed a priori, or require a truncating call plus LOH or a deep deletion.
  - Report MSI/hypermutated models as their own context.
  - Restrict the genotype hypothesis space a priori to (i) subunits of frozen complexes and (ii) pre-registered drivers.
  - State the resulting family sizes before testing.

## SHOULD

**S1. Make the replication rule power-aware.**
- Problem: the split-half rule (discovery FDR < 0.1, validation p < 0.05, plus the median over 10 splits) is near-powerless for contexts with about 4 models per half. By a z-approximation at d = 1.5: discovery power ≈ 0.2.
- Change:
  - Contexts with ≥20 models per half: split-half replication.
  - Rarer contexts: full-data FDR, plus leave-one-subunit-out (LOSO), plus same-sign DEMETER2 replication, plus leave-one-model-out.
  - Drop the "median over 10 splits" criterion; report it descriptively.

**S2. Specify the compute plan for permutations.**
- Problem: about 10³ complexes × 10³ contexts × 10⁴ permutations is not feasible on the laptop.
- Change: use analytic OLS/HC3 p-values for all tests. Check their calibration by within-lineage permutation on a random subset (QQ, null control). Use exact permutation only for candidate hits below FDR 0.1.

**S3. Split the positive controls into two classes.**
- Problem: WRN–MSI and KEAP1–NFE2L2 are gene-level controls. A complex-level pipeline may legitimately miss them.
- Change: pre-register (a) **pipeline/gene-level controls** and (b) **complex-level controls**, the latter mapped to frozen clusters after the freeze and before any test. Report recovery with denominators for each class.

**S4. Pre-specify the complex-vs-gene benchmark.**
- Problem: the benchmark in §2.4 has no fixed metric.
- Change: fix the metric before testing, e.g. held-out replication rate and DEMETER2 concordance at matched numbers of tests and context families. This is the second methodological claim.

**S5. Write the Part 1 outcome paragraph both ways before the look.**
- Change: pre-write how the paper reads if zavg_BEF fails criterion §7(1) against GLS. Part 2 uses zavg_BEF regardless, and M2 then carries the "does it matter" argument.

**S6. Simplify construction to one primary map.**
- Change:
  - **Primary:** Leiden-CPM with consensus over 20 cell-line bootstraps, and two nested γ levels.
  - **Sensitivity only:** soft overlap, the eigengene and singscore scores, and BioPlex-integrated variants.
  - Fix the selection grid now (as in §2.2) and do not extend it.

## COULD
- **C1.** Note the PC-removal trade-off explicitly. Lineage-structured co-dependency, such as SOX10/MITF, is partly removed by the top-20-PC filter, so lineage TF modules may fragment. Report whether the lineage positive controls form clusters; this is a known limitation, not something to tune.
- **C2.** Add one line to firewall F3: hu.MAP **train** labels are also excluded from Part 2 construction and parameter choice. They come from curated complexes.
- **C3.** Endorse the protocol §9 recommendations as drafted: keep the CI criterion, exclude self pairs (R2), keep train-overlap pairs (R4), use rule (b) as the primary operating point, and freeze the fitted rules. Freeze after M1.

## What a top-journal referee attacks first
(1) "Not complexes" → M4. (2) "Wainberg/Chen already did this; controls are rediscoveries" → M2, M3. (3) Lineage, hypermutation and screen-quality confounding → M5, S2. (4) Selective inference on rare contexts → S1.
