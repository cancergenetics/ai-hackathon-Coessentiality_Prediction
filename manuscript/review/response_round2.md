# Response to the referee, round 2

We thank the referee again. All the requested points are text corrections or counts taken from committed files, and we have made each of them.

**What was new in this round.** We ran three post hoc computations (R-5), all on committed result tables. The intent note `results/part2/revision/INTENT_round2.md` (4a7c715) was committed first; the computations followed in `scripts/62_rev_round2.py` (5fac192) (outputs `results/part2/revision/r5_*`, README section "R-5").
- None of them refits anything.
- None opens a dataset, a lockbox value, a persisted read or a test label.
- No registered result, label, call or rank changed.

Locations below refer to `manuscript/main.md`, `manuscript/supplementary.md` ("SN") and the tables. The clean copy in `manuscript/submission/` was regenerated. It has an abstract of 149 words and an Introduction, Results and Discussion of 4,805 words.

---

## MUST

**N1. H1 at ±0.05; "mostly"; our earlier reason.**
- **Count.** The referee's count is confirmed. None of the 11 neutral events has a committed 90% CI inside ±0.05, and none sits on the boundary. The smallest margin that would contain each CI runs from 0.059 (STK11) to 0.097 (KDM6A). The neutral call therefore depends entirely on the registered ±0.1 margin. This is stated in the Results (H1 section), SN 15 item 19 and SN 17 (R-5), and the cover letter says it too.
- **Heading.** The Discussion heading now reads "Within complexes, 11 of 24 losses leave the complex-wide dependency within ±0.1 Chronos units, and none within ±0.05". It is followed by "'Neutral' here means no change larger than the registered margin".
- **Our earlier reason.** The round-1 letter said the ±0.05 sensitivity was "not in committed outputs". That was wrong, and we have marked it as corrected in `response_round1.md` (M6 and the summary table). The reason still holds for the margin in residual-SD units, which needs each fit's residual SD. That SD was not exported, and the text says so.

**N2. "Reliable lineage".**
- **Wording.** The phrase is removed.
- **Rates in the text.** Each highlighted lead now carries its lineage's never-expressed p < 0.05 rate: Ovary/Fallopian Tube 0.22, Bone 0.19, Head and Neck 0.19 and Myeloid 0.01. The text then states two things:
  - All four lie beyond the trend in its own direction (spline residual z −2.79 to −4.36), so the scale rule cannot exclude an artefact for them, and DEMETER2 is their main independent support.
  - Only Myeloid has a calibrated never-expressed null.
- **Other lineages.** The Results add that leads outside Lymphoid and melanoma are "not thereby reliable": their contexts' rates range from 0.01 to 0.22.
- **Table S9.** The note says the same. Its `rev_lineage_less_reliable` flag marks only the three highest rates and does not certify the other lineages. `strongest_derived.json` now carries the four leads' rates, trend direction and DEMETER2 p, which are the sources the main-text numbers cite.
- **Round-1 letter.** The corresponding sentence in `response_round1.md` (M2) is marked as corrected.

**N3. Scope of the independent check.**
- **Methods.** The "Revision analyses" paragraph now states that the agent:
  - had read the results summary first, so the check was not blind;
  - reran scripts 59–61 and recomputed the numbers from the committed tables with its own code;
  - did not rerun script 58 (the random atlases and the non-member test table);
  - reused, rather than re-implemented, the empirical-Bayes step of the power aid.
- **Elsewhere.** SN 17 gives the same scope in more detail. The cover letter now says "recomputed … not blind and did not rerun the random atlases" instead of "independently reproduced". The round-1 letter is annotated.

## SHOULD

- **M7, alternative registry readings.**
  - *H1 leave-one-out unit.* The unit enters only the H1 rows, so H3a (207) and H3b (45) are unchanged. Of the five eligible H1 calls, only ARID1A P and EP300 P can change; the other three fail on other criteria. The count is therefore 252–254. Both calls passed all 30 lineage folds. Whether either would pass by primary disease is not in the committed outputs. Deciding it needs an A refit, which we did not make, and we say so.
  - *Ranking.* Ranking selects nothing, so the 254 do not change. Ranking by |t_A| alone reorders the list strongly (Spearman 0.37 with the registered rank; top-10 overlap 5). A product |t_A| × tier weight (3/2/1, illustrative) reorders it mildly (Spearman 0.89; top-10 overlap 9). No count in the text depends on rank.
  - *Where reported.* Results ("The registered strongest findings"), SN 16 and SN 17.
- **N4.** "Coherence, not complex membership, carries most of the excess" is removed from the main text and the revision README. The text now reads: "sets matched for coherence behave like complexes, and sets matched for essentiality alone associate less often". It adds that coherence is measured on A as co-variation across lines, so matching on it partly conditions on the outcome.
- **N5.** Across the ten random atlases, the share of random-atlas discoveries beating all ten others was 0.291–0.364 for H3a (real 0.464) and 0.448–0.517 for H3b (real 0.669). No atlas reached the real rate. The range appears in the Results, the Discussion and SN 17, labelled as a spread across atlases and not a confidence interval (`r5_atlas_flag_range.tsv`).
- **N6. Table S7.** The "notes" column now uses the revised label (R-2 `rev_scale_artefact_possible`, any lineage). The first-version heuristic (every lymphoid β > 0) is kept only as the full-TSV column `scale_caveat_v1`, and the table note explains the change (`scripts/tab_lockbox.py`). The Z0016 row, for example, no longer carries a caveat.
- **N7.** The revision README now gives 0.097 (0.5785 − 0.4817 = 0.0968). The main text and supplement already used 0.097. The D20 tree note recorded the correction earlier, and D19 does not contain the value.
- **N8.** The cover letter now says "was tested in a single read of 338 held-out lines".
- **N9. Fig. 6.** The panel b legend now says that the open symbols are contexts, with area proportional to the number of real discoveries. The submission figure list already names `fig6_revision_references.pdf`, generated by `make_submission_copy.py`.
- **Number checker on the supplement.**
  - *Checker changes.* `scripts/check_manuscript_numbers.py` now reads the supplement's per-note alias comments and accepts `--start-at`, so a range of notes can be checked.
  - *Notes 16 and 17* (strongest findings and the revision analyses) pass at 320 of 320 (`results/verify_manuscript/number_check_supplement.tsv`). Several comments were corrected to make this possible, and one value (13.3%) needed its correct source cited.
  - *The whole supplement* matches 952 of 1,208 numbers. The remainder are 55 unmatched and 201 with no source comment, almost all in Notes 1–15. Those notes predate the convention and cite sources at file level or not at all.
  - *What we did not do.* We did not rewrite the comments of Notes 1–15. The Methods now say that only the main text and Notes 16–17 are fully machine-checked. The cover letter is not checked.
- **M4.1 (BENCH).** "Needs the lockbox values again" is replaced in the Results, SN 15 and SN 17 by this: the persisted reads would allow a re-run, but we chose not to make a second, post hoc look at L, to keep the single-read design, so the registered BENCH result stays inconclusive.
- **Optional (M3.2).** We did not add the implied-SE sentence. The interpretation of the CI bounds already in the text covers it.

## Not changed

- No registered result, call, label or rank.
- The figure artwork, apart from the Fig. 6 legend text.
- The author placeholders the referee lists ([AUTHORS], [CORRESPONDING_AUTHOR], [COMPETING_INTERESTS], [REPO_URL], [ARCHIVE_DOI], [PART1_REPO_URL], [DEPMAP_RELEASE_CONFIRMATION], [HGNC_DOWNLOAD_DATE] and the licence checks). These still need the authors.
