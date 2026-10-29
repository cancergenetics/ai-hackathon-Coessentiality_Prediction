# Referee report, round 2

**Manuscript:** "A co-essentiality score selected on training complexes wins a single held-out test but yields neither better-supported complexes nor more replicable complex-level dependencies in DepMap"
**Version reviewed:** main at 77c86fe. I reviewed manuscript/main.md, manuscript/submission/, supplementary.md (Notes 1–18), cover_letter.md, Fig. 1–6, Tables S7 and S9, response_round1.md, results/part2/revision/ (INTENT, README, r1–r4 outputs), results/part2/verify/revision_analyses_review.md and RESEARCH_TREE.md D18–D20.
**Access.** I read committed result files only. I did not open lockbox values, persisted reads, test labels or datasets/, and I computed no metric on test pairs. My one computation was on the committed A-side H1 table (`results/part2/directions/H1_within_complex_GI/classes.tsv`); see N1 below.

---

## 1. Overall assessment and recommendation

The revision is substantive and handled in good faith.
- The central message has been narrowed to what the matched-set design supports.
- Four post hoc, A-only analyses were run. Their intent (4157f1e, 13:43) was committed before the results (740f439, 14:50), and they were checked by an independent agent.
- BENCH is now reported as inconclusive.
- The title, abstract and cover letter drop "label-free", "locked before the test data existed" and "random gene sets".
- Every major issue from round 1 has been addressed at least in part.

My numerical checks of the revision outputs all match (Section 6). `scripts/check_manuscript_numbers.py` passes 393/393 on the revised main text. The submission copy equals the main text with source comments removed, plus a list of figure files.

Some new or remaining statements still overstate the evidence. Two of them can be checked directly against committed files and are wrong or unsupported. Neither needs new data, and both can be fixed in an afternoon.

**Recommendation: minor revision.** The MUST items are listed in Section 7.

---

## 2. Round-1 major issues

| Item | Verdict | Comment |
|---|---|---|
| **M1** Matched sets are fragments of other complexes; the flag has no selection-aware reference | **Resolved** | Comparators are described accurately in the abstract, Results, Discussion, Methods and the Fig. 3a legend. R-1(a) adds non-member references; R-1(b) adds a selection-aware reference for the flag (46.4% vs 32.9% H3a; 66.9% vs 48.3% H3b). The flag no longer defines "complex-specific". The 11,861/80,500 vs 1,145/8,441 comparison is in the main text. Residual points are N4 and N5 (SHOULD). |
| **M2** Lineage-scale artefact not carried through to the highlights | **Partly resolved** | (a) spline trend: done. (b) one rule in all lineages: done (95 vs 56). (e) DEMETER2: done and informative. (c) not done but disclosed (26/37 leads lie beyond the trend in its own direction). (d) replaced by a reliability column. The text that presents the highlighted leads introduces a new overstatement (N2, MUST). |
| **M3** Bridge framing; power aid applied asymmetrically | **Resolved** | BRIDGE-L is called secondary throughout. The width of the primary bridge's interval is interpreted. R-4 applies the power aid to both maps: power accounts for 26–41% of the gap. The text discloses that the headline model under-predicts both maps' levels. The title's "better-supported" is acceptable; "no evidence of better-supported" would be more exact. |
| **M4** BENCH bootstrap; Discussion overstated | **Resolved for reporting** (re-run declined, acceptably; Section 3) | Reported as inconclusive with a point estimate favouring the complex. "The complex is not the right unit" is removed. |
| **M5** Part 1 wording, concentration of the gain, baselines | **Resolved** (analyses declined, acceptably; Section 3) | Wording corrected. The abstract now carries the macro-AP null and the concentration in pair-rich clusters. The Discussion says preprocessing probably contributes. |
| **M6** H1 neutrality | **Partly resolved** | The ARID1A caveat, the Fig. 4d count and the wording "within ±0.1, not no change" are fixed. The ±0.05 sensitivity was declined as "not in committed outputs", but the committed 90% CIs answer it (N1, MUST). The Discussion heading says "mostly", which 11/24 does not support (N1, MUST). |
| **M7** Strongest-findings list | **Mostly resolved** | Testability split given (199 full, 49 partial, 6 rare route). Rare-route wording corrected in the cover letter and Discussion. Timing stated. The alternative registry readings were not computed and no reason is given (Section 3; SHOULD). |

## 3. Are the "not done" justifications acceptable?

| Request | Authors' reason | Judgement |
|---|---|---|
| M5: preprocessing-matched baseline; leave-one-large-cluster-out ΔAP | A second look at the test labels | **Acceptable.** Any new test-set metric is a second look, and the protocol forbids it. The descriptive export after the look was a precedent, but declining is the conservative choice. |
| M4.1: BENCH re-run with an explicit failure policy | Needs the lockbox values again | **Acceptable, but stated too strongly.** The persisted, hash-checked L matrices exist, and the corrected relook (R7) shows they can be re-read non-confirmatorily through the resume path. Re-reading them is possible, so this is a choice. Say "we chose not to re-open the persisted L reads" rather than implying it is impossible. Reporting the interval as unreliable and the result as inconclusive meets the substance of the request. |
| M3.2: formal minimal detectable difference | Bootstrap SD not committed | **Acceptable.** Interpreting the CI bounds is enough. One optional sentence: the reported CI implies an SE of about 0.027, i.e. (0.031 + 0.075)/3.92, so the design had roughly 80% power only for differences of about 0.075. |
| M2(d): rescaling by never-expressed inflation | New computation not in committed results | **Weak but tolerable.** R-1 to R-4 were all new A-only computations, so "not in committed results" is not a reason. The reliability column is a reasonable substitute, provided it is not presented as more than it is (N2). |
| M6: residual-SD units and ±0.05 sensitivity | Not in committed outputs | **Not acceptable for the ±0.05 sensitivity.** `classes.tsv` commits `kappa_N_ci90_lo`/`hi` for every event, so the sensitivity is a column filter. I ran it: **0 of the 11 neutral events have a 90% CI inside ±0.05** (Section 4, N1). Residual-SD units would need the residual SD per fit, which is not committed. Declining that part is acceptable if stated. |
| M7: counts under the alternative registry readings | "Not computed", no reason | **Not acceptable as written.** It is an A-only rerun of `scripts/56_strongest_findings.py` (118 s), and the stakes are small: the H1 leave-one-out unit affects only the two H1 calls, and ranking affects only order. Either run it or state that it can change at most those two H1 calls and the row order. |
| Minor 5: CI for 13.3% vs 9.1% | Superseded by R-1(b) | **Acceptable.** |
| Minor 9: L effect of g* alone | Not in committed outputs | **Acceptable.** I confirmed that `bench_pairs.tsv` holds Δ but not t_L(g*). |
| Number checker on supplement and cover letter | Supplement lacks the alias convention | **Acceptable for now, still open.** I spot-checked Supplementary Note 17 by hand (Section 6). |
| Figure layout fixes (Fig. 2f, 3a, 3e, 4c) | Handled in legends | **Acceptable** for proofs. |
| Minor 17 (Holm bound) | Referee used one-sided 1/B | **I accept the correction.** The protocol's two-sided p is 2·min(·), so the bound is < 2/B = 0.001 and the Holm bound is < 0.003. |

## 4. New issues introduced by the revision

**N1 (MUST). H1 "neutral" at a tighter margin, and "mostly".**
- I applied the ±0.05 sensitivity to the committed `classes.tsv`. None of the 11 NEUTRAL events has a 90% CI inside ±0.05; for example NF1 is −0.022 to +0.065 and STK11 is −0.059 to +0.009.
- The main text says this sensitivity "was not exported", yet it can be read from committed columns. It matters, because it shows the neutral call depends entirely on the registered ±0.1 margin.
- The Discussion heading "Within complexes, losses mostly leave the complex-wide dependency within ±0.1 Chronos units" is not supported: 11 of 24 is 46%.

*Requested.*
- Report "0 of 11 at ±0.05" in the Results and Supplementary Note 15.
- Change the heading, for example to "Within complexes, 11 of 24 losses leave the complex-wide dependency within ±0.1 Chronos units".
- Correct the response letter's M6 reason.

**N2 (MUST). "Reliable lineage" for the highlighted examples.**
- The Results introduce the four bullet examples as "T1, in a reliable lineage and supported in DEMETER2".
- "Less reliable" is set only for melanoma, Skin and Lymphoid (`scripts/tab_strongest.py` `LESS_RELIABLE`), where never-expressed sets reach 27–33%.
- Three of the four highlighted examples are in lineages whose never-expressed rate is also about four times the nominal 0.05 (Table S9 `rev_never_expressed_frac`):
  - Ovary/Fallopian Tube 0.22;
  - Bone 0.19;
  - Head and Neck 0.19.
- Only Myeloid (0.01) has a calibrated never-expressed null.
- All four highlighted examples are also trend-concordant (`rev_concordant_nl` True, residual z −2.79 to −4.36). By the paper's own statement, the scale rule therefore cannot exclude an artefact for them.
- DEMETER2 support (p 0.0003–0.009) is the real argument for these four, and it is a good one.

*Requested.*
- Replace "reliable lineage" with the numeric never-expressed rate for each example.
- State that all four lie beyond the trend in its own direction, so the DEMETER2 check is what supports them.
- Make the same correction in Table S9's note.

**N3 (MUST, wording). Scope of the independent review.**
- Methods ("Revision analyses") and the response letter say "an independent review reproduced every number".
- The review (`revision_analyses_review.md`) says "Every number I could recompute matches". It also says:
  - the reviewer read the README before computing;
  - the R-1(b) random atlases and their permutation p-values, and the R-1(a) non-member test table (script 58), were not rerun;
  - the Tweedie step in R-4 was reused, not re-implemented.

*Requested.* State this scope in one clause.

**N4 (SHOULD). "Coherence, not complex membership, carries most of the excess."**
- This new causal-sounding sentence rests on 90 of 367 units, the less coherent quarter. The caveat follows, which is good.
- Coherence is itself measured on A as co-variation across lines. A set that co-varies across lines is more likely to shift with lineage, so matching on coherence partly conditions on the outcome.
- Suggest "sets matched for coherence behave like complexes; sets matched for essentiality alone associate less often", and drop the attribution.

**N5 (SHOULD). Uncertainty for the R-1(b) excess.**
- The Discussion says flagged complexes beat their matched sets "more often than selection alone predicts" (46.4% vs 32.9%), but gives no uncertainty.
- The ten random atlases give a natural range: per-atlas discoveries are in `r1_summary.json` and per-discovery flags in `r1b_random_discoveries.tsv.gz`.
- Report the range of the per-atlas flag rate.

**N6 (minor). Table S7 caveat is inconsistent with the revised rule.**
- Table S7's Markdown "notes" column still marks every Lymphoid β > 0 row "scale caveat". This is the first-version heuristic.
- Of the 55 Lymphoid β > 0 strongest findings, 21 are not `rev_scale_artefact_possible` under the revised rule. Z0016 is one of them, a non-concordant lead that is still flagged.
- Use the revised label, or say that the S7 note is the first-version heuristic.

**N7 (minor). Rounding mismatch.** The GLS under-prediction is 0.097 in the main text and response letter, and 0.096 in the revision README and RESEARCH_TREE D19. The value 0.5785 − 0.4817 = 0.0968 rounds to 0.097, so correct the README.

**N8 (minor). Cover letter: "was confirmed".** The cover letter says "Every complex-level claim built on the map was confirmed in a single read of 338 held-out lines". This sentence was carried over from version 1 and I missed it in round 1. About half were not confirmed. Say "was tested".

**N9 (minor). Fig. 6.**
- In panel b, marker size varies with context, but the legend explains only "large symbols, pooled". Add what size encodes.
- The figure file is `fig6_revision_references.*`, while the legend cites the script `fig_6_revision.py`. That is consistent, but give the file name in the submission figure list as for the other figures.
- Panel c is otherwise clear.

## 5. Round-1 minor issues, statistics, figures, reproducibility, circularity, title/abstract

**Minor issues 1–18.**
- **Resolved:** 1, 2, 3, 4, 6, 7, 8, 10, 11, 12, 13, 14, 15, 17, 18.
- **Resolved by an acceptable substitution:** 5 and 9.
- **Partly resolved:** 16. The submission copy strips comments, and the word counts are within NC limits (abstract 149; main 4,547). A structured abstract for Genome Biology is deferred, which is fine.

**Statistical review.**
- BB assumptions (independence across families, PRDS within), hybrid p-values, a simulation that planted lineage effects only, and q ≈ 0.046 with 21/23 families selected: **resolved**.
- Power aid cannot separate biology from reproducible artefacts: **resolved**.
- Lockbox shares processing with A and was part of Part 1's development data: **resolved**.
- Timing of the robustness rule and labels; scope of the never-expressed response; lost run log: **resolved**.

**Figures and tables.**
- Fig. 1, 2, 4 and 5 legend corrections: **resolved**.
- Fig. 3: legend fixes **resolved**; artwork overlaps deferred to proofs, which is acceptable.
- Fig. 4: panel d now counts 26 labels; the n in panel a and the SMARCB1 type are explained. **Resolved.**
- Table 1 Holm bound: **resolved**.
- Table S7 now lists the registered findings: **resolved**, apart from N6.
- Table S9 revised columns: **resolved**, apart from N2.
- New Fig. 6 is legible and matches its source files (minor N9).

**Reproducibility.**
- Checker: 393/393 on main (92 by text match), and no tracked file changed when I ran it.
- The supplement is still not machine-checked (open, SHOULD).
- Seeds for R-1 to R-4 are documented (4600/4700, 4000/4100, `PYTHONHASHSEED=0`).
- Intent before results: verified in git (4157f1e at 13:43, then 740f439 at 14:50).
- The independent review is PASS-WITH-NOTES, with its SHOULD items applied (bc61c46).

**Circularity.**
- **Resolved.** The Part 2 scope sentence, the indirect route through the score's configuration, BioPlex tiers steering the highlights, and leads as a post hoc selection for presentation are all stated.
- The revision analyses read no hu.MAP, STRING or curated resource (README; firewall checked by the reviewer for scripts 58–61). I agree with the reviewer's recommendation to widen the firewall test glob to `scripts/[3-6]*_*.py`.

**Title and abstract.**
- **Resolved.** The abstract now matches the evidence: macro-AP null, primary/secondary bridge labels, "post hoc" power, an accurate description of the comparators and the essentiality-only reference.
- One nuance remains in the title: "yields neither better-supported complexes" rests on a null primary bridge. That is acceptable given the explicit "absence of evidence" sentence in the Results.

## 6. Spot-checks of revised numbers against committed files

| Claim (where) | Source | Value | Match |
|---|---|---|---|
| Essentiality-only sets 16.6% vs 28.1% (abstract, Results) | `r1a_rates.tsv`; README table | 0.166 vs 0.281 (355 units) | yes |
| N-nbr 25.1% vs complexes 23.0%, N-ess 16.4% on 90 units | README "Like-for-like" | 0.251 / 0.230 / 0.164 | yes |
| Flag 46.4% vs 32.9% (H3a); 66.9% vs 48.3% (H3b); P(\|t\|) 0.48 | `r1_summary.json` `b.summary` | 0.4635 / 0.3291; 0.669 / 0.483; 0.482 | yes |
| 95 scale-artefact-possible; 37 leads; 26 concordant | `r2_summary.json` | 95; 37; 26 | yes |
| DEMETER2: 423 lines; leads 29/35 same sign, 20 p < 0.05, ρ 0.86; all 157/194 | `r3_summary.json` | identical | yes |
| BRIDGE-L power: expected −0.023 (SD 0.019), remainder −0.042, z −2.2; zavg 0.459 vs 0.514 | `r4_summary.json` `headline` | −0.0229 (0.0189); −0.0419; −2.21; 0.4588 / 0.5137 | yes |
| Ovary Z0035 never-expressed 0.22, concordant, D2 p 0.00041 | Table S9 `rev_*` | 0.22; True; 0.00041 | yes (see N2) |

## 7. Items remaining before submission

**MUST (scientific or wording; no new data needed):**
1. **N1.** Report that 0/11 H1 neutral events hold at ±0.05, from committed `classes.tsv`. Replace "mostly" in the Discussion heading. Correct the response letter's M6 reason.
2. **N2.** Drop "reliable lineage" for the highlighted examples. Give their never-expressed rates (0.19–0.22 for Ovary, Bone and Head and Neck), and state that all four are trend-concordant and rely on DEMETER2.
3. **N3.** State the scope of the independent review of R-1 to R-4: not blind, atlases and the Tweedie step not rerun.

**SHOULD:**
- M7: run the alternative readings or state their bounded impact.
- N4: reword the coherence attribution.
- N5: give the atlas range for the R-1(b) excess.
- N6: Table S7 caveat column.
- Extend the number checker to the supplement.
- M4.1: reword "needs the lockbox values" as a choice.

**Known author placeholders (not scientific blockers, but required before submission):** [AUTHORS], [CORRESPONDING_AUTHOR], [COMPETING_INTERESTS], [REPO_URL], [ARCHIVE_DOI], [PART1_REPO_URL], [DEPMAP_RELEASE_CONFIRMATION], [HGNC_DOWNLOAD_DATE] and the three licence checks.
