# What was learned

**The main point.** The most useful lessons are not about any one complex. A better pair score did not give better complexes. In DepMap, "complex-level" cancer dependencies behave like dependencies of any coherent gene module. And a pre-registered, held-out design is what made these negative answers visible. Most other realisations are engineering details, listed last.

Insights are sorted into three kinds:
1. **Scientific**: what a researcher in the field would want to know, regardless of this method.
2. **Methodological**: what carries over to similar analyses.
3. **Engineering**: what you need only to reproduce this work.

Within each kind they are ranked by how general and how surprising they are. "Shown" means a committed result supports it. "Believed" marks an interpretation the project did not test.

## The short version (half a page)

**Scientific**
1. **Better pair accuracy did not give better complexes.**
   - *Shown:* the pair score beat GLS on held-out complexes (AP 0.333 vs 0.212). Yet complexes built from it were no better supported by physical interaction data than complexes built from GLS (−0.025 [−0.075, +0.031]). They were also *less* often confirmed in held-out cell lines (−0.065 [−0.107, −0.025]).
   - *Why:* the pair-level gain sat in large complexes with many pairs, which add little to a map of small ones.
   - [Turning points 1](JOURNEY.md#tp1) and [4](JOURNEY.md#tp4).
2. **Context dependencies of complexes are mostly module-level, not complex-specific.**
   - *Shown:* gene sets of the same size, essentiality and co-dependency coherence, drawn from other complexes, were associated with lineages exactly as often as real complexes (28.1% vs 28.1% of tests). Sets matched only for size and essentiality were associated less often (16.6%).
   - So what predicts a "complex dependency" is coherence plus essentiality, not complex membership.
   - [Turning point 5](JOURNEY.md#tp5); [Step 5](METHOD.md#step5).
3. **Lineages differ in the scale of CRISPR effects, and this masquerades as lineage-specific dependency.**
   - *Shown:* sets of genes that are never expressed (so should have no real dependency) reach p < 0.05 in 13.5% of lineage tests overall, and in 27–33% of lymphoid and melanoma tests.
   - Highly essential complexes look more lineage-specific than they are.
   - Any lineage atlas built on Chronos gene effects needs a scale null.
   - [Turning point 8](JOURNEY.md#tp8); [Step 5](METHOD.md#step5).
4. **Replication at the expected rate does not prove specificity.**
   - *Shown:* 52.5% of lineage findings replicated in held-out lines, about what statistical power predicted (47–52%).
   - Generic and effect-scale signals replicate just as well as specific ones.
   - [Turning point 6](JOURNEY.md#tp6); [Step 6](METHOD.md#step6).

**Methodological**
1. **Use a selection-aware reference for "beats its matched sets".**
   - *Shown:* among unselected tests, a complex beats all 10 of its matched sets 9% of the time. Among *discoveries*, the rate rises to 33% even for random sets that went through the same selection. Real discoveries reach 46%.
   - Comparing 46% with 9% would have overstated specificity five-fold.
   - [Turning point 8](JOURNEY.md#tp8); [Step 5](METHOD.md#step5).
2. **Say what the matched comparator is.**
   - Sets drawn from other complexes, coherent non-member modules and essentiality-only sets answer three different questions.
   - The project first called its comparator "random gene sets"; it was in fact fragments of other complexes.
   - [Turning point 5](JOURNEY.md#tp5).
3. **A split-half drawn from the same pool used for tuning is not a hold-out.**
   - Seal a patient-blocked, stratified lockbox before reading any outcome, and read it once per registered analysis.
   - [Turning point 2](JOURNEY.md#tp2); [Step 1](METHOD.md#step1).
4. **Dress rehearsals must assert outcomes, not just run.**
   - *Shown:* a 2.7-hour rehearsal passed while 53 model fits silently failed, because nobody checked fit status. Check the fit-status counts and that the held-out code path equals the construction path on identical rows.
   - [Turning point 6](JOURNEY.md#tp6).

## 1. Scientific and domain insights

Ranked by generality and surprise.

1. **Pair accuracy is not complex utility.** Expanded from the short version above.
   - The held-out gain over GLS grows with test-complex size, from −0.033 for complexes with one positive pair to +0.142 for those with more than 100 (`results/test_eval/export/ap_by_size.tsv`).
   - Averaged complex by complex, the two scores tie (macro AP 0.297 vs 0.292).
   - *Believed, not tested:* the extra pairs the score gets right sit inside large complexes, which clustering would find anyway.
   - [Turning point 1](JOURNEY.md#tp1), [Turning point 4](JOURNEY.md#tp4).
2. **Complex-level associations are shared by coherent co-dependency modules.**
   - On the 90 complexes where coherent *non-member* modules could be matched, those modules were associated as often as the real complexes (25.1% vs 23.0% of lineage tests; `results/part2/revision/r1a_rates.tsv`).
   - Caveat: such modules could be found for only the less coherent quarter of complexes.
   - [Turning point 5](JOURNEY.md#tp5).
3. **Lineage effect-scale shifts need their own null.**
   - The lymphoid shift is visible as a correlation between a complex's essentiality and its lymphoid coefficient (ρ = −0.24).
   - 26 of the 37 post hoc "lead" findings lie beyond the lineage trend *in its own direction*. A trend-based rule therefore cannot rule out a scale artefact for them.
   - [Turning point 8](JOURNEY.md#tp8).
4. **In a raw co-essentiality graph, never-essential genes link mainly to neighbours on the same chromosome arm.**
   - *Shown:* 70–81% of graph edges between never-essential genes join genes on the same arm, against 3% for random pairs (design review, `docs/reviews/round1_algorithms.md`). This is the known CRISPR copy-number/proximity artefact.
   - A gene filter (keep genes with some dependency) and a locus flag are needed. 28% of the frozen complexes are still flagged.
   - [Turning point 3](JOURNEY.md#tp3); [Step 2](METHOD.md#step2).
5. **Known paralog buffering is recovered, but "neutral" loss is a statement about the margin.**
   - Loss of ARID1A raises dependency on its complex's paralog ARID1B; EP300 loss does the same for CREBBP. These are recoveries of known biology.
   - Most tumour-suppressor losses leave their complex's dependency within ±0.1 Chronos units (11 of 24), but 0 of those 11 stay within ±0.05 (`results/part2/revision/r5_h1_margin_005.tsv`).
   - [Turning point 8](JOURNEY.md#tp8).
6. **In the worked example, snRNA transcription and snRNP assembly form one co-dependency neighbourhood.**
   - Shown in the example: SNRPB2 (an SMN/Gemin complex member) has strong links to three SNAPC subunits. Consensus clustering still separated the two machines.
   - The known SNAPc partner POU2F1 (Oct-1) links to SNAPC1 but is not included in the complex.
   - One example only.
   - [Worked example](METHOD.md#step2).

## 2. Methodological insights

1. **Selection-aware references.** Run random sets through the *same* selection, then compare rates among discoveries; see the short version. *Generalises to:* any "is it more specific than chance?" claim made after an FDR screen. [Turning point 8](JOURNEY.md#tp8).
2. **Name what a comparator is.** See the short version. [Turning point 5](JOURNEY.md#tp5).
3. **Lockbox, not split-half.** See the short version. [Turning point 2](JOURNEY.md#tp2).
4. **Rehearsals with asserted outcomes, and footprint audits before remediation.**
   - The defect first showed as 53 failed fits; the full audit found 560 affected fits.
   - Scope a fix by a label-free audit before writing its plan.
   - [Turning point 6](JOURNEY.md#tp6).
5. **After a defect found post-read, keep the read as executed as primary.**
   - Report a corrected, non-confirmatory relook alongside it, and strip claims from every call that differs between the two (53 calls).
   - This keeps the single-look guarantee honest without hiding the defect.
   - [Turning point 6](JOURNEY.md#tp6).
6. **Test FDR code against the guarantee the procedure actually gives.**
   - Benjamini–Bogomolov (2014) bounds the average false-discovery proportion *within selected families*, not the pooled one. *Shown in simulation:* 0.046 against 0.068 pooled.
   - Assert it over many simulations, never one seed.
   - [Turning point 5](JOURNEY.md#tp5); [Step 4](METHOD.md#step4).
7. **Map every registered step to code, a test and an output before the gated run.**
   - Four registered steps had no code when their gate ran: calibrations, lockbox driver, covariate design, robustness checks.
   - [Turning point 7](JOURNEY.md#tp7).
8. **Power-check confirmation rates.**
   - An empirical-Bayes power estimate (Efron 2011) turned "52% replicated" from an ambiguous number into "about what power predicts".
   - [Turning point 6](JOURNEY.md#tp6).
9. **Verify the evaluation harness on an old answer before the single look.**
   - The test harness first reproduced Part 1's internal test to 1.5 × 10⁻⁸.
   - Independent verifiers commit their numbers before seeing the harness's.
   - [Turning point 1](JOURNEY.md#tp1).
10. **Machine-check every reported number.**
    - Every manuscript number carries a source comment and a script checks it: 421/421 in the main text, 1,181/1,181 in the supplement.
    - The check found citation errors, not value errors.
    - [Turning point 7](JOURNEY.md#tp7).

## 3. Engineering details (reproduction only)

- **Seeds.** Fix `PYTHONHASHSEED` and all seeds; Leiden consensus then rebuilds byte-identically (444 s). [Step 2](METHOD.md#step2).
- **Permutation p-values.** These are extrapolated with a generalised Pareto tail only when no permutation reaches the observed statistic. Read p < 10⁻⁴ as "below the permutation floor" (worked example: 1.4 × 10⁻⁶). [Step 4](METHOD.md#step4).
- **Lockbox reads.** Persist each read with its hash committed at once; the resume path never re-calls the loader. [Step 6](METHOD.md#step6).
- **Waiting on background jobs.** Wait on process IDs, not on text in a log. A waiter keyed on the word "real" (zsh prints "total") ran for 3 h 43 min. [Research process](JOURNEY.md#the-research-process).
- **Commits.** Batch mechanical commits: five commits per lockbox entry made 153 of 335 commits. Always commit the run log of an irreversible run; the lockbox loader's log was lost.
- **Agent worktrees.** Agents' working copies were once created from a stale commit. Check the base and merge or delete agent branches at each merge.
- **Public Part 1 release.** It has a fresh git history, so the locked Part 1 commit is checked by git object ids at the public commit 3c684c9 (`external/README.md`).
- **Hash-pinned files.** Their SHA-256 is fixed in the frozen test protocol, so they cannot be edited, even to remove local paths.
- **Matched-set seeding.** Matched sets are seeded by complex index, so one complex's sets cannot be regenerated alone. The worked example re-implemented the inner loop.
