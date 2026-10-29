# Response to the referee, round 1

We thank the referee for a careful and constructive report. The referee's checks of our numbers were all correct, and every major issue was about how results should be stated. We agree with almost all of them. The revision narrows several claims, adds four post hoc, descriptive analyses on the construction set A (R-1 to R-4; `results/part2/revision/`, intent 4157f1e committed before computation, independently reproduced in `results/part2/verify/revision_analyses_review.md`; *[round 2: the check was not blind and did not rerun script 58 or re-implement the empirical-Bayes power step; see N3]*), and corrects the wording the referee flagged. No registered result, call or metric changed. No lockbox value or test-set file was opened for the revision. Where a request needed new data access or new computation that is not in committed results, we did not do it and say so below.

Locations refer to the revised `manuscript/main.md` (sections), `manuscript/supplementary.md` (Supplementary Notes, "SN"), figures and tables. A clean copy without source comments is in `manuscript/submission/`.

---

## Major issues

### M1. Matched sets are fragments of other complexes; the "complex-specific" flag has no selection-aware reference

> *"The parity result … shows that a complex is not more context-associated than a coherent fragment of another complex … It does not show that associations are 'not about complexes'."*

We agree, and the central message is rewritten accordingly.

1. **Accurate description.** Everywhere, including the abstract, the comparators are now "gene sets from other co-dependency modules, matched for size, essentiality and coherence". The Results heading is now "Context associations are matched by other coherent co-dependency modules", and the Discussion heading "…matched by other coherent modules, not by any gene set". The phrase "moves almost any coherent, essentiality-matched gene set" is replaced by a statement limited to coherent modules.
2. **Second calibration (R-1(a), post hoc).** Using the existing calibration machinery, we drew sets from the 2,380 S genes outside every complex with three pre-declared generators: coherent neighbourhoods with full matching, uniform draws with full matching, and uniform draws matched on size and essentiality only. Results (Results, "Context associations…"; Fig. 6a; SN 17):
   - Uniform non-member sets never reach complex-level coherence.
   - Coherent non-member sets could be matched for only 90 of 367 units, the less coherent ones. On those 90 units they reached p < 0.05 in 25.1% of H3a tests, against 23.0% for the complexes; essentiality-only sets on the same units reached 16.4% (H3b: 10.8%, 11.3%, 8.4%).
   - Over the 355 units with essentiality-only sets, those sets reached 16.6% against 28.1%.
   - The reading is now: a complex is not more context-associated than a coherent, essentiality-matched co-dependency module; sets matched only for essentiality are associated less often. The 90/367 caveat is stated next to the claim.
3. **Selection-conditional reference (R-1(b)).** Ten random atlases, each replacing every complex by one of its matched sets, went through the full discovery procedure (screen, 10,000-permutation re-test, selective FDR). Among their discoveries, 32.9% beat all ten others in their group, against 46.4% of real discoveries (H3b 48.3% against 66.9%). Conditional on discovery, real complexes' |t| was not larger than random discoveries' (mean probability 0.48; |β| 0.49). We say plainly that the excess is relative to a complex's own matched sets, not an effect-size advantage. The flag no longer defines "complex-specific" anywhere in the text (Results; Fig. 6b; Fig. 3d legend; SN 17).
4. **BB over matched-set tests.** Now in the main text: 11,861 rejections from 80,500 tests (14.7%), against 1,145 from 8,441 (13.6%).

Minor 13 (17 units without sets, 12 with doubled tolerances) is stated where the calibration is described.

### M2. The lineage-scale artefact is not carried through to the highlighted findings

We agree that the first-version labels were uneven. R-2 and R-3 replace them; the first-version columns stay in Table S9 for traceability.

- **(a) Non-linear trend.** A natural cubic spline (4 df) of β on unit essentiality was fitted per context over all 367 units. It raises R² only slightly (Lymphoid 0.068 → 0.120; all others ≤ 0.083). 79 findings depart from it by ≥ 2 residual SDs, as with the linear fit. 35 of the 37 first-version leads are kept, Bone Z0130 and Myeloid Z0278 are lost, and CNS/Brain Z0146 and Breast Z0253 are new. The flagged Z0190, Z0196, Z0363 and Z0001 stay above the threshold under the spline. Four leads lie within 0.05 of the threshold, so the lead set is sensitive to it (SN 17).
- **(b) All lineages.** One rule now applies in every lineage. It labels 95 of the 207 H3a findings a possible scale artefact, against 56 before (Results, "The registered strongest findings"; Table S9 `rev_*`).
- **(c) Scale-concordant leads.** We did not exclude them, because the rule as specified cannot flag a finding beyond the trend. Instead the text states plainly that **26 of the 37 leads lie beyond the trend in its own direction, so a scale artefact cannot be excluded for them.** The lymphoid Z0190 finding is explicitly not interpreted.
- **(d) Never-expressed sets as an empirical null.** We did not rescale statistics by a per-context inflation factor: this needs new computation that is not in the committed results. Instead, Table S9 reports each context's never-expressed rate. The nine leads in Lymphoid (27%) and melanoma (33%) are marked less reliable in the text and table. This is listed as a limitation (Discussion; SN 15, 17).
- **(e) DEMETER2 concordance (R-3).** The same score and design were refitted on DEMETER2 for the 423 A lines with RNAi data. 29 of 35 evaluable leads agree in sign, 20 have one-sided p < 0.05, and the rank correlation is 0.86. Over all 207 H3a findings, 157 of 194 agree in sign. The six discordant leads are named (Soft Tissue Z0080 and Z0336, Kidney Z0323, Pancreas Z0283, Lung and Bone Z0001), as are four that agree in sign but are weak (Lymphoid and melanoma Z0190, Pleura Z0067, Lung Z0150) (Results; Fig. 6c; Table S9; SN 17).
- **Point 2 (positive-ρ lineages).** Skin and Head and Neck now appear in the text, and the Fig. 3b legend points to them. Bone SNAPC and Head and Neck EGFR/GRB2 remain examples: they are beyond the spline trend, outside the less reliable lineages and supported in DEMETER2. Their trend direction is shown in Table S9.
- **Point 3 (never-expressed excess).** The excess now discounts leads in the affected lineages, as described under (d).

The text now leads only with findings that are T1, in a reliable lineage and supported in DEMETER2. *[Round 2: "reliable lineage" was not supported. Three of the four examples are in lineages with never-expressed rates of 0.19–0.22, and all four lie beyond the trend in its own direction; the text now gives the rates and relies on DEMETER2 (N2).]* It states that the labels are a post hoc selection for presentation.

### M3. Bridges framed too strongly; power argument applied asymmetrically

1. **BRIDGE-L is called the registered secondary bridge throughout,** including the abstract ("primary bridge", "secondary bridge") and its own Results heading.
2. **Detectable difference for the BioPlex bridge.** We did not compute a formal minimal detectable difference; the bootstrap SD is not in the committed output. We now state what the registered interval implies: it is compatible with an advantage of up to 0.031 for zavg or 0.075 for GLS, so the result is absence of evidence of a difference, not evidence of equivalence (Results, "A complex map, and a null primary bridge"; SN 15 item 5).
3. **Power aid for BRIDGE-L (R-4).** The Supplementary Note 4 aid was applied unchanged to both maps' hits, using A estimates, a per-map empirical-Bayes prior and only committed L sample sizes and calls. GLS-map discoveries have larger A effects (median |β| 0.081 against 0.066), so they are expected to confirm more often. Across three power versions the expected difference is −0.017 to −0.026 against the observed −0.065. Power therefore accounts for roughly a quarter to two fifths of it (26–41%; headline −0.023, SD 0.019). The remainder (−0.042, about 2.2 SDs) is called approximate, because that SD treats hits as independent. The headline model also under-predicts both maps' observed levels (zavg by 0.055, GLS by 0.097) (Results, "The secondary bridge"; Fig. 6d; SN 17). The Discussion applies the power argument to complex-specificity and to BRIDGE-L alike.
4. **Title.** The new title is "A co-essentiality score selected on training complexes wins a single held-out test but yields neither better-supported complexes nor more replicable complex-level dependencies in DepMap". We use "better-supported", not "better", because the primary bridge is null rather than negative. The Discussion heading now reads "The advantage did not translate into better-supported complexes".

### M4. BENCH: partly failed bootstrap; Discussion overstates the result

1. **Re-run with an explicit failure policy.** Not done. Any bootstrap re-run needs the lockbox values again, and this revision opened no lockbox data. We list it as an open item. *[Round 2: stated too strongly. The persisted, hash-checked L reads exist and could be re-read; we chose not to make a second, post hoc look at L, to keep the single-read design.]*
2. **Reporting.** The result is now "inconclusive, with a point estimate favouring the complex" (+0.240; 60% of high-tercile pairs favour the complex). The registered interval is called unreliable (16.5% failed fits; 1,032 of 2,000 replicates on reduced pair sets; lower bound −0.122 with per-replicate design) (Results, "Complex scores against their best gene: inconclusive"; SN 5, 6, 15).
3. **"The complex is not the right unit"** is removed, and so are "complex scores did not replicate better" and "the conclusion holds either way". The low-tercile result is described as largely expected from the definition of breadth.

### M5. Part 1 scope

1. **"Locked before the test data existed"** is removed everywhere. The text now says the method was locked 59.6 min before the test files reached the project, and that the split was published earlier [Fischer 2025] (abstract, Results, cover letter, SN 1).
2. **"Label-free"** is removed from the title, abstract, text and cover letter. The score is described as unsupervised, with no parameters fitted to labels and a configuration selected on hu.MAP training complexes (Introduction, Methods). The Supplementary no longer uses the word either.
3. **Concentration of the gain.** The abstract now states it and says macro per-cluster AP did not differ. The Results keep the Kish effective number (21.3 of 238) and the size breakdown.
4. **Preprocessing-matched baseline and leave-one-large-cluster-out ΔAP.** Not done. Both need the test labels again, which would be a second look at the test set; the protocol forbids it, and this revision opened no test-set file. The Discussion now says that preprocessing probably contributes to the gain, citing the post hoc PC-removed Pearson reference, and that we could not test a preprocessing-matched GLS.

### M6. H1 neutrality

- **Count.** The text now says: "For 11 of the 24 events, the mean complex-wide (non-paralog) dependency changed by less than ±0.1 Chronos units (90% CI inside the margin)". It adds that for ARID1A the paralogue did change, so 10 events showed no change on any tested axis. The abstract no longer carries an H1 count, and the cover letter is aligned (minor 11).
- **Fig. 4d.** The panel title now reads "H1, 24 events: complex-wide labels + 2 paralog calls", and the legend explains that the bars count 26 labels.
- **Residual-SD units and ±0.05 sensitivity.** Not done; neither is in the committed outputs. The text and SN 15 say so. *[Correction, round 2: this reason was wrong for the ±0.05 sensitivity. The committed `classes.tsv` holds every event's 90% CI, so the sensitivity is a column filter: 0 of the 11 neutral events stay neutral at ±0.05 (R-5; `manuscript/review/response_round2.md`, N1). The reason stands for residual-SD units, which need the per-fit residual SD, which was not exported.]*
- **Margin against H3a effects.** The text now notes that 325 H3a discoveries are smaller than half the margin, so "neutral" means "no change larger than 0.1", not "no change". L re-establishes equivalence for only 2 of 11 events (SN 12, 15).

### M7. Strongest findings

- **Alternative registry readings.** Not done. The counts under the alternative leave-one-out unit and ranking were not computed, and the text says so. *[Round 2: now bounded from committed tables (R-5): 252–254 findings; ranking changes only the order.]*
- **Testability split.** Results and SN 16 now give it. The leave-one-out check applies fully to 199 findings and partly to 49; 6 passed the rare-context route instead (EP300 → CREBBP and five EWSR1-fusion findings). The strict (Bowel) and vacuous (Pancreas) cases are stated.
- **Wording.** The cover letter now says the 254 were "confirmed in held-out lines (or, for five, by a rare-context route that is not an independent replication)", and that EP300 → CREBBP "passed only the rare-context route". The Discussion says the same.
- **Timing.** The Results and cover letter state that the rule was applied after the lockbox read. The abstract no longer mentions the 254, to keep within 150 words.

---

## Minor issues

1. **Methods scope.** The sentence is now limited to Part 2. A second sentence adds that Part 1's validation set derives from curated complexes by design (Methods, "Data and versions").
2. **"Curation-free".** The Introduction now defines it narrowly and states the indirect routes: the score's configuration, and controls that encode literature. The term is dropped from figure and section titles.
3. **Lockbox independence.** Added to the Results (lockbox section), the Limitations and SN 15 item 26b. L lines were in Part 1's development data and were processed jointly with A, so confirmation shows reproducibility within the release.
4. **Error control in the abstract.** The abstract now says "1,145 lineage associations (selective FDR)" and "as post hoc power predicts".
5. **CI for 13.3% vs 9.1%.** Not computed. The unconditional contrast is superseded by the selection-aware reference of R-1(b) and is kept only in SN 17, with a note that its binomial reference is approximate.
6. **H1 BH denominator.** A clause in the Methods ("Permutation and error control") now gives PTEN's q at the registered m (0.051).
7. **H1S counts.** Reconciled in the Methods and SN 12. A1 moved 318 events; A2 moved SMARCB1 to H1, leaving 317 events with 321 tests, 293 of them evaluable (registered m 338).
8. **Member lists.** Units are now given with member counts and "including" (Z0035, 17 members; Z0051, 16; Z0254, 12; Z0069, 14). The lymphoid Z0030 example moved to SN 17 with its member count (18).
9. **Ovary lead and g*.** For each highlighted lead, SN 17 now gives g* and the committed BENCH difference Δ. For Ovary, g* is PAX8, and the complex replicated better than PAX8 (Δ +1.05). The text says the lineage factors could drive the association alone. The L effect of g* by itself is not in the committed outputs and was not computed.
10. **Kill rule.** The text now says that passing the lenient rule shows that the machinery recovers known biology, not that it has power for new effects.
11. **Cover-letter wording.** Aligned with M6.
12. **Wainberg/Chen statement.** We could not check those papers' null models without new downloads, so the claim is removed. The Discussion now says we have not compared our references with their nulls and make no claim about how much of their signal is complex-specific.
13. **Matched-set availability.** Stated in the Results, as described under M1.
14. **Fig. 4c timing.** The legend now says the profiles were exported from A after the H1 lockbox read. The Methods say the same.
15. **Internal identifiers.** The main text now uses plain descriptions for rulings, amendments and decision notes; Table 3 has no identifiers. SN 18 is a glossary. Unit IDs are kept only as Table S9 look-up keys and are defined at first use.
16. **Format.** The abstract has 149 words; the Introduction, Results and Discussion have 4,547; the Methods have 1,604. HTML comments are stripped in `manuscript/submission/` by `scripts/make_submission_copy.py`. We kept the Nature Communications format; a structured abstract for Genome Biology can be prepared if needed.
17. **Holm bound.** Table 1 now shows p < 0.001 and Holm p < 0.003 where no replicate fell at or below 0 (`scripts/tab_manuscript.py`). The bound is 3 × 2/B, because the protocol's p is two-sided. The referee's 0.0015 used the one-sided 1/B. SN 10 is corrected to match.
18. **Timeline.** SN 1 now gives the order: Part 1 locked 12:55; test files arrived 13:55; DESIGN v1 15:08; v3.1 17:37; protocol frozen 18:14; labels first read 18:18; results 18:52; registry draft 19:27. It states that the Part 2 registry was written after the Part 1 result was known.

---

## Statistical review

- **BB assumptions.** The Methods and SN 15 (item 11b) now say that the guarantee assumes independence across families and PRDS within them, which the overlapping contexts here do not satisfy. They also say that the simulation planted lineage effects only, with none for the H3b dependence structure; that the p-values are a hybrid; and that with 21 of 23 families selected, BB is close to within-family BH at q ≈ 0.046.
- **Power aid and reproducible artefacts.** Now stated for the overall rates (Results, lockbox section; Discussion; SN 15 item 26b). We also note that observed rates are close to the predicted bands, with H3a slightly above.
- **Power aid for BRIDGE-L.** Done (R-4); see M3.
- **(i) Timing of the robustness rule and labels.** Stated in the Results, Discussion (Table 3) and cover letter.
- **(ii) Response scope of the never-expressed calibration.** Now in the main text: the literal reading was fixed after the discovery outcomes were known and before the calibrations ran.
- **(iii) Lockbox loader run log.** It is permanently lost. The Methods, SN 15 item 26 and `data_code_availability.md` now say so.

---

## Figures and tables

- **Fig. 1.** The legend now covers candidate B and the PC-removed Pearson reference, says that AP without the 10 largest clusters was not bootstrapped, and explains the under-prediction visible in panel c.
- **Fig. 2.** Panel e's axis now reads "subsample maps" (`scripts/fig_2_construction.py`, regenerated), and the legend says these are 80% subsamples, not bootstraps. The legend gives the colour key for panel e and states what z0 means for GLS in panel c (the edge count the GLS graph matches). Not changed: we did not move the overlapping zero line in panel f.
- **Fig. 3.** The legend now explains the grey dots in panel d, assigns the dotted 0.05 line in panel e to the top axis, and points to the positive-ρ lineages in panel b. The panel d legend also says that its flag is largely produced by selection. Not changed: the layout fixes for the inset legend overlap in panel a and for the crowding between the panel e axis label and panel d. We dealt with these in the legends only; the artwork can be adjusted at the proof stage.
- **Fig. 4.** The panel d title now counts labels per axis (regenerated). The legend explains the n in panel a (defined loss lines; FOCAD has 30 analysed) and the SMARCB1 "silencing" type field. Not fixed: the overlap of the panel letter "c" with the ARID1A sub-panel title.
- **Fig. 5.** The legend now explains that EP300 is absent from the H1 row (it went by the rare-context route), and that a confirmed call such as ARID1A can have an HC3 CI crossing 0 because of the permutation rule. The panel c legend now calls the interval unreliable.
- **New Fig. 6** (`scripts/fig_6_revision.py`, generated deterministically from committed revision files). It shows the matched-set comparison (a), the "beats all" flag conditional on discovery (b), DEMETER2 concordance (c) and BRIDGE-L observed against power-expected rates (d). This brings the display items to nine: six figures and three tables. Fig. 3d stays, with a corrected legend.
- **Table 1/1b.** The Holm bound is fixed (minor 17).
- **Table S7.** The Markdown now lists the registered strongest findings in rank order. The earlier subset survives only as a column of the full TSV (`scripts/tab_lockbox.py`).
- **Table S9.** Revised labels (`rev_*`) are added alongside the first-version labels. Z0190 now appears with its trend direction, the less reliable lineage flag and its weak DEMETER2 result. The text sets it aside explicitly.

---

## Reproducibility

- **Data versions and licences.** We cannot resolve these ourselves. They are marked as placeholders for the authors: [DEPMAP_RELEASE_CONFIRMATION], [HGNC_DOWNLOAD_DATE], [LICENCE_BIOPLEX], [LICENCE_CCLE_PROTEOMICS] and [LICENCE_HUMAP] (Methods table; `data_code_availability.md`).
- **Code location.** Marked as [REPO_URL], [ARCHIVE_DOI] and [PART1_REPO_URL] (commit c7f95ee).
- **Number checker.** The checker passes on the revised main text: 393 of 393 numbers matched (92 by text match). We did not extend it to the supplement or cover letter, and we did not tighten the text matches to cited sections. The supplement does not yet follow the alias convention the checker needs: an attempted run left most numbers unresolved rather than mismatched. This remains open.

## Circularity

We agree with each point. The indirect route of curated information through the score's configuration is stated in the Introduction, Discussion and Methods. The statement about Part 1 validation is scoped (minor 1). SN 15 item 16 now says that BioPlex tiers (validation data) steer which findings are highlighted. It also says that the lead label is a post hoc selection for presentation, not a purely descriptive annotation.

## Title and abstract

Both are rewritten (M3.4; M5). The abstract no longer contains "label-free", "locked before the test data existed", "random gene sets", the H1 count or the 254. It now carries the macro-AP result, "selective FDR", "post hoc" and the primary/secondary labels of the bridges.

## Items blocking submission

These need the authors and remain open: the author list and corresponding author ([AUTHORS], [CORRESPONDING_AUTHOR]); competing interests ([COMPETING_INTERESTS]); repository URL and archive DOI, including public availability of the Part 1 repository; licence checks; the HGNC download date; and confirmation of the DepMap release. HTML comments are stripped in the submission copy, and the lengths are within the journal limits.

---

## Summary of requests not carried out

| Request | Reason |
|---|---|
| M2(d) rescaling by never-expressed inflation | New computation, not in committed results; reliability flag reported instead |
| M3.2 formal minimal detectable difference | Bootstrap SD not committed; the CI bounds are interpreted instead |
| M4.1 BENCH re-run with a failure policy | Needs the lockbox values again *[round 2: a choice, not an impossibility; see M4.1]* |
| M5 preprocessing-matched baseline; leave-one-large-cluster-out ΔAP | Needs the test labels again (a second look) |
| M6 margin in residual-SD units; ±0.05 sensitivity | Not in committed outputs *[round 2: true for residual-SD units only; the ±0.05 sensitivity is now reported, 0 of 11]* |
| M7 counts under the alternative registry readings | Not computed *[round 2: bounded, R-5]* |
| Minor 5 CI for 13.3% vs 9.1% | Superseded by the R-1(b) reference |
| Minor 9 L effect of g* alone | Not in committed outputs; Δ reported |
| Number checker on the supplement and cover letter | The supplement needs the alias convention first |
| Some figure layout fixes (Fig. 2f, 3a, 3e, 4c) | Handled in the legends; artwork left unchanged |
