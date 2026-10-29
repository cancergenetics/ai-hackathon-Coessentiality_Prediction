# Referee report, round 1

**Manuscript:** "A label-free co-essentiality score wins a pre-registered held-out test, but its advantage does not carry over to complexes or to complex-level cancer dependencies"
**Files reviewed:** manuscript/main.md, supplementary.md, figures/fig1–fig5 (PNG), tables/ (Tables 1, 2, S1–S9 and the derived JSON files), cover_letter.md, data_code_availability.md. I checked claims against results/, docs/registry/, docs/TEST_PROTOCOL.md, DESIGN.md, RESEARCH_TREE.md and the git history. I did not open datasets/ or any hu.MAP files.
**Expertise:** cancer functional genomics (DepMap/CRISPR), protein-complex biology, network methods, statistics.

---

## 1. Summary of the claims

The study has two linked parts.

**Part 1.** A co-essentiality score, zavg_BEF, was developed and locked in earlier work. It is the mean of three z-scored similarity components computed on PC-removed, ridge-whitened DepMap Chronos profiles. It was evaluated once on the withheld hu.MAP 3.0 test split (600,738 pairs, 3,312 positives) against signed Pearson, Spearman and GLS. Pooled AP was 0.333, against 0.212 for GLS (paired cluster-bootstrap Δ +0.121, 95% CI +0.079 to +0.159). The gain held for pairs whose proteins were never training positives (C3). The locked "clearly beats" rule was met. Macro per-cluster AP did not differ from GLS (+0.005), and the gain is concentrated in clusters with many positive pairs.

**Part 2.** The locked score was recomputed on a construction set A of 1,007 lines, with a PatientID-blocked lockbox L of 338 lines held back. It was turned into a consensus-Leiden complex map: 367 level-1 complexes covering 3,642 genes. An edge-matched GLS map served as comparator. The authors claim the following.
- **Bridges.** The registered BioPlex bridge is null (Δ −0.025, −0.075 to +0.031). The lockbox bridge BRIDGE-L is negative: zavg-map discoveries confirm less often than GLS-map discoveries (Δ −0.065, −0.107 to −0.025).
- **Context associations (H3a lineage, H3b driver/genome state).** There were 1,145 and 161 discoveries under Benjamini–Bogomolov selective FDR. Matched gene sets reach p < 0.05 as often as real complexes (28.1% vs 28.1% in H3a). Never-expressed gene sets show a lineage excess (13.5%), and a "lineage-scale" (Chronos dynamic-range) effect is seen in Lymphoid.
- **Within-complex genetic interactions (H1).** The panel recovers known paralog and partner dependencies (ARID1A→ARID1B, EP300→CREBBP, FOCAD→PELO–HBS1L). 11 of 24 losses are "neutral" by TOST at ±0.1 Chronos units.
- **Sub-module divergence (H2, exploratory).** It is largely reproduced by random member splits.
- **Lockbox.** A defect in the confirmation driver was found after the single read and handled by post-read rulings (R6, R7). 53 calls are labelled defect-sensitive and carry no claim. About half of discoveries confirm (H3a 52.5%; H3b 40–49%), which a post hoc empirical-Bayes power aid predicts. Complex scores do not out-replicate their best single gene (BENCH; registered CI includes 0).
- **Strongest findings.** A registered robustness rule, computed after the read on A only, leaves 254 "strongest findings". Post hoc labels pick 37 "lead examples" and 56 "possible lineage-scale artefacts".

The overall message is that pair-level gains in co-complex scoring do not carry over to better complexes. It is also that complex-level context dependencies in DepMap are largely gene-set-generic, so complex-specificity must be tested against matched gene sets.

## 2. Overall assessment

**Novelty.** Moderate to high for the field's methodology. I know of no other DepMap complex-dependency study that combines a single blind test look, a PatientID-blocked lockbox, a machine-enforced registry with a ledger, and a full public record of deviations. The two substantive negative results are useful and will be cited: the bridge failure, and matched-set parity for context associations. The individual biological findings are mostly re-discoveries, as the authors say. The new candidates (the 254 findings) are descriptive.

**Significance.** The practical lesson is that complex-level dependency atlases need matched-set calibration. This is important if it is framed correctly. As written, the strongest form of the lesson ("most complex-level context associations are not about complexes") goes beyond what the matched-set design can show (Major 1).

**Rigour.** The process is unusually rigorous: freeze markers, hashed intents, single-read ledgers, independent re-implementations, a number checker and a candid deviations table. The main weaknesses are inferential, not procedural:
- what the matched sets are (they are fragments of other complexes);
- how the lineage-scale artefact is propagated into the highlighted findings;
- the symmetric use of the power argument;
- a compromised BENCH bootstrap;
- several abstract and title phrasings that are stronger than the evidence.

**Clarity.** The text is accurate but dense. It relies heavily on internal identifiers (R4, R7, D14, E2, A3–A5, FZ-C, "rare-context route", Z-unit IDs), which a general reader cannot follow without the supplement. The main text plus Conclusions runs to about 6,300 words, and the abstract to 189 words. Both exceed Nature Communications norms (about 5,000 and 150 words). Genome Biology would need a structured abstract.

**Recommendation: major revision.** The core results are sound and reproducible from committed files. The requested changes need no new data, but several change how the paper's main messages should be stated.

---

## 3. Major issues

### M1. The "matched random gene sets" are fragments of other complexes. This limits what the parity result can show, and the "complex-specific" flag has no selection-adjusted reference.

*Concern.* The abstract and Results describe "Random gene sets matched for size, essentiality and coherence". The operational specification (results/part2/directions/CALIBRATION_D8_INTENT.md, item 1) draws each set "from the members of one other frozen coarse unit … (else the union of two)", matched on coherence within ±1.0 z. The comparators are therefore pieces of other co-essential modules, with real complex-level biology of their own. The intent note says so ("background, not a null"), and so does main.md l. 86. The parity result (28.1% vs 28.1%) therefore shows that a complex is not more context-associated than a coherent fragment of another complex of similar essentiality. It does not show that associations are "not about complexes" (Discussion heading, l. 167) or that a contrast "moves almost any coherent, essentiality-matched gene set" (l. 86). The registry text (H3a_lineage_atlas.yaml, calibration) says only "random gene sets matched for size, essentiality and coherence". Drawing the sets from other units was an operational choice fixed after the discovery outcomes were known (R4, D8). That choice shapes the headline reading.

Second, the "beats all matched sets" flag defines "complex-specific" (l. 88). It compares discoveries, which are selected for large |t|, with unselected matched sets. The unconditional reference (13.3% vs 9.1%, l. 88) is appropriate. The conditional rate among discoveries (502 of 1,083 = 46%) has no reference and is inflated by selection alone.

*Why it matters.* This is the paper's central biological message and the main take-home in the abstract, conclusions and cover letter.

*Requested changes (A only, no lockbox access).*
1. Describe the comparators accurately everywhere, including the abstract: "gene sets drawn from other co-dependency modules, matched for size, essentiality and coherence". Reword the Discussion heading and l. 86 accordingly.
2. Add a second matched-set calibration with the existing `src/part2/calib.py` machinery, labelled post hoc. Draw sets from S genes outside any level-1 complex, or from all of S, matched on essentiality, with coherence matching (a) relaxed or (b) removed. This separates "any coherent module" from "any essentiality-matched gene set".
3. Give a selection-conditional reference for the flag. Among matched-set tests that are themselves BB rejections (11,861 in H3a per calibration_D8.json `matched_random_sets.bb`), report the fraction whose |t| exceeds the other 9 sets and the real complex in that context. Compare 502/1,083 with this, not with 1/11.
4. Report in the main text that BB over the matched-set tests made 11,861 rejections from 80,500 tests (14.7%), against 1,145 from 8,441 (13.6%) for real complexes. This is the most direct statement of parity and is currently buried in Supplementary Note 12.

### M2. The lineage-scale artefact is not consistently carried through to the highlighted findings

*Concern.* The paper documents two lineage-level artefacts: never-expressed sets at 13.5% p < 0.05 (concentrated in melanoma, Skin and Lymphoid), and a Lymphoid ρ = −0.240 between unit essentiality and β. Yet the "lead example" and "possible lineage-scale artefact" labels (D16, post hoc) handle these artefacts unevenly.

1. **Linear trend, concordant extremes promoted.** The artefact label requires |residual z| < 2 from a *linear* per-context trend of β on unit essentiality. Concordant findings that are *more* extreme than the line become "lead examples". A dynamic-range compression is expected to be stronger for the most essential units, so a linear fit will leave exactly these as positive residuals. Table S9 lists Lymphoid/Myeloid leads that are concordant, positive-β and on highly essential units: Z0196 (unit essentiality −0.69, z +2.23), Z0363 (−0.77, +2.45), Z0001 (Myeloid, −0.85, +2.19) and Z0190 (−1.38, +5.88). The main text (l. 155) declines to interpret Z0190 for exactly this reason. Table S9 still lists it as a lead example, and Z0196, Z0363 and Z0001 are treated as leads. Overall, 26 of the 37 lead examples are scale-concordant (tableS9_strongest_findings.tsv).
2. **Asymmetry across lineages.** Skin, Head and Neck and Bone show ρ of about +0.19 to +0.21 (Fig. 3b; Supplementary Note 12), similar in size to Myeloid (−0.135). Only Lymphoid and Myeloid findings can receive the artefact label. Two of the six text-highlighted leads are in lineages with positive ρ: Bone SNAPC1–5 and Head and Neck EGFR/GRB2.
3. **The never-expressed excess is not propagated.** The registered response did not fire on the literal reading, and the excess is "reported as a finding". It is never used to discount the strongest findings in the affected lineages (melanoma 11, Skin, Lymphoid 92 of 207 H3a findings).

*Why it matters.* The six biological vignettes and the cover letter's "what survives" rest on these labels.

*Requested changes (A only).*
- (a) Fit a monotone (isotonic) or smooth (LOESS or spline) trend of β on unit essentiality per context. Re-derive both labels and report how many of the 37 leads and 56 artefacts change.
- (b) Apply the artefact rule to every context with |ρ| above a stated threshold, not only to Lymphoid and Myeloid.
- (c) Exclude scale-concordant findings from "lead example", or explicitly justify including them.
- (d) Use the 100 never-expressed sets as a per-context empirical null. For example, compute a per-context inflation factor from their t-distribution and re-evaluate the strongest findings with rescaled statistics, reported descriptively.
- (e) Report DEMETER2 concordance (sign and correlation) on A lines for all 254 findings. DEMETER2 is already in the pipeline for the rare-context route and is not subject to Chronos's dynamic-range compression. This is the most direct within-data check on whether the lead examples are platform artefacts.

### M3. The bridges are framed more strongly than the registry and the evidence support, and the power argument is applied asymmetrically

*Concern.*
- The title and abstract state that the advantage "does not carry over to complexes". The BioPlex bridge, which is the *primary* bridge endpoint (BRIDGE-L_confirmation_rate.yaml: "the primary bridge endpoint is BioPlex precision"), is a null with a CI of −0.075 to +0.031. That is absence of evidence, and the minimal detectable difference is not stated.
- BRIDGE-L is registered as a *secondary* bridge endpoint (same file, and Supplementary Note 6: "secondary bridge"). The abstract and Discussion ("Both registered bridges failed") present the two as co-equal.
- For complex-specificity (l. 125), the manuscript argues that a 0.083 confirmation gap "is about the size that effect size and power alone predict". For BRIDGE-L it says "We have no tested explanation" (l. 129). The zavg map has larger units (median 9 vs 5) and covers more genes, so its discoveries may have different A effect sizes and SEs from GLS-map discoveries. The same post hoc power aid could show whether the −0.065 is expected from those differences.

*Requested changes.*
1. Call BRIDGE-L the registered secondary bridge throughout, including the abstract.
2. For the BioPlex bridge, report the difference the design could detect (for example from the bootstrap SD) so that "null" is interpretable.
3. Apply the Supplementary Note 4 power aid to the zavg-map and GLS-map BRIDGE-L hits, using A estimates, the same EB prior construction and L line counts. Report the predicted rate difference next to the observed −0.065. This needs only committed A-side outputs and the existing L confirmation statuses, not a new lockbox read.
4. Soften the title, for example "…does not yield better complexes or more replicable complex-level dependencies". Make "does not carry over" conditional on the measured endpoints.

### M4. BENCH: the registered interval comes from a partly failed bootstrap, and the Discussion overstates the result

*Concern.* 16.5% of pair × replicate fits failed silently, and 1,032 of 2,000 replicates averaged over a reduced set of pairs (Supplementary Note 12; lockbox_derived.json `BENCH.frac_pair_rep_nan`, `n_rep_with_any_nan`). The pairs dropped are not random: they are small contexts with leverage-1 lines. The resulting percentile CI for the registered high tercile (−0.033 to +0.362) is therefore of uncertain coverage. Its lower bound moves to −0.122 under a different covariate assembly. The point estimate in the registered tercile *favours the complex* (+0.240; 60% of pairs Δ > 0). The low-tercile result (best gene wins) is largely expected by construction, because breadth is defined as few members associated. Even so, the Discussion concludes "For most context dependencies, the complex is not the right unit" (l. 169), and the abstract and Conclusions repeat that complex scores "did not replicate better".

*Requested changes.*
1. Re-run the BENCH bootstrap with an explicit failure policy and report both CIs. Either keep the pair set fixed and resample within strata so that leverage-1 designs cannot arise, or drop failing pairs consistently across all replicates.
2. Report the high-tercile estimate as "inconclusive, point estimate favouring the complex" rather than as a negative.
3. Remove or qualify "the complex is not the right unit". The registered test does not support it, and the low-tercile comparison is partly tautological.

### M5. Part 1 scope: wording, concentration of the gain and the baseline set

*Concern.*
1. **"Locked before the test data existed"** (abstract, Conclusions, cover letter). The hu.MAP 3.0 test split was published in 2025. The score was locked 59.6 min before the test *files reached the project* (Supplementary Note 1). Say this.
2. **"Label-free."** zavg_BEF has no label-fitted parameters, but its configuration was selected by complex-grouped CV on hu.MAP training labels derived from Complex Portal (Methods l. 221; Discussion l. 163). "Unsupervised score, configuration selected on training labels" is accurate. "Label-free" in the title is not.
3. **Concentration of the gain.** The Kish effective number of clusters is 21.3 of 238. The six clusters with more than 100 positive pairs contain 1,454 of 3,312 positives (44%; Fig. 1f). Macro per-cluster AP shows no gain (+0.005, −0.014 to +0.047). The abstract reports none of this, although the main text does.
4. **Baselines.** The post hoc PC-removed Pearson closes a large part of the gap to Pearson (0.239 vs 0.197; Table 1). This suggests preprocessing contributes a substantial share of the gain. No baseline receives the same preprocessing as zavg_BEF (top-20 PC removal plus ridge whitening) with GLS or Spearman similarity.

*Requested changes.*
- Correct items 1 and 2 in the title, abstract and cover letter.
- Add one abstract sentence on macro AP and the concentration in pair-rich clusters.
- As a clearly labelled post hoc, non-confirmatory reference (like PC-rm Pearson), report a preprocessing-matched GLS or Spearman baseline. Also report the ΔAP against GLS after leaving out each of the six largest clusters in turn. Both use existing data.

### M6. H1 "neutrality": counting, margin and abstract claim

*Concern.*
1. The abstract says "11 of 24 genetic losses left their complex's dependency unchanged". One of the 11 is ARID1A, whose complex-wide axis is neutral but whose paralogue ARID1B, a member of the same unit Z0109, becomes markedly more essential. The complex's dependency was therefore *not* unchanged.
2. Fig. 4d is titled "H1 (24 events)", but its bars sum to 26. They count tests (axes), not events.
3. The ±0.1 Chronos margin is registered. However, 325 of the 1,145 H3a discoveries have |β| < 0.05 and are reported as real context effects. TP53 (κ_N +0.056) is labelled "within margin". The paper thus applies one effect-size standard to call H3a effects real and a looser one to call H1 effects absent. The registry (H1_within_complex_GI.yaml l. 414) also requires the margin to be "reported in residual-SD units". I did not find this in the manuscript.
4. L re-establishes equivalence for only 2 of 11 events.

*Requested changes.*
- Rephrase as "for 11 of 24 events the mean complex-wide (non-paralog) dependency changed by less than ±0.1 Chronos units", and exclude or footnote ARID1A.
- Retitle Fig. 4d as "26 tests in 24 events".
- Report the margin in residual-SD units as registered, and add a ±0.05 sensitivity.
- Discuss the margin relative to the H3a effect-size distribution.

### M7. The strongest-findings list: post-read readings, uneven robustness and rare-route items

*Concern.*
1. The registered per-hit robustness checks were computed after the lockbox read. This is disclosed, and the candidate list was committed first (fd03b03).
2. After the run, the PI accepted two readings of ambiguous registry text: the H1 leave-one-out unit is the lineage, and ranking is by tier, then |t_A| (Supplementary Note 16). The outcome under the alternative readings is not reported.
3. The leave-one-out check is vacuous for some contexts (Pancreas 17/17 pass because folds are skipped) and very strict for others (Bowel 0/29). It is partly untestable for 54 of the 254 findings, including all 28 TP53-altered findings. "Robust to leaving out each primary disease or lineage" (abstract, cover letter) is therefore not uniformly true.
4. The cover letter says the 254 include "2 within-complex calls that were confirmed in held-out lines". EP300→CREBBP was not confirmed in L (4 L lines; L CI −1.82 to +1.07). It passed the rare-context route, which the paper itself says "is not an independent replication". The same applies to the five EWSR1-fusion H3b findings.

*Requested changes.*
- Report the counts under the alternative registry readings.
- Split Table S9 tallies into "leave-one-out fully testable", "partly testable" and "rare-context route".
- Correct the cover letter and abstract wording ("confirmed in held-out lines or passed the rare-context route").
- State in the abstract that the robustness rule was applied after the lockbox read.

---

## 4. Minor issues

1. **Methods l. 217.** "STRING, hu.MAP predictions and curated complexes were used neither in construction nor in validation" is literally false for Part 1, whose validation set (hu.MAP test pairs) derives from curated complexes. Scope the sentence to Part 2.
2. **"Curation-free map".** The map is built from a score whose configuration was selected on Complex Portal-derived labels (M5.2). The 23 positive controls and the kill rule also encode literature knowledge. Qualify the term once in the Introduction.
3. **Lockbox independence.** L lines were part of the data on which zavg_BEF was developed in Part 1 ("on all U lines the recomputation reproduces Part 1 scores"). L Chronos values were also inferred jointly with A in the same release. Lockbox confirmation therefore establishes reproducibility within DepMap 26Q3 processing, not independence from shared artefacts. State this in the Limitations; currently only "one CRISPR resource" is mentioned.
4. **Error control in the abstract.** "1,145 … discoveries" and "about half … confirmed" carry no qualifier. Add "selective FDR" and "post hoc power calculation" respectively.
5. **Specificity exchangeability (l. 88).** Give a CI or permutation p for 13.3% vs 9.1%. Tests share complexes and lines, so the binomial reference is approximate.
6. **H1 BH denominator.** BH used the evaluable m (24) instead of the registered 28. PTEN's q changes from 0.044 to 0.051 (Supplementary Note 12). This is disclosed in the supplement only; add one clause in the Methods.
7. **H1S counts.** A1 moved 318 events, the main text says 317 events, Supplementary Note 12 says 321 tests, and the registered m is 338. Reconcile these in one sentence.
8. **Incomplete member lists given as complete.**
   - l. 151: AP2 unit "(AP2A1, AP2B1, AP2M1, AP2S1, PICALM, DNM2, TFRC)"; Z0254 has 12 members, including AGFG1, FCHO2, CCDC32, C10orf88 and SLC11A2.
   - l. 153: the Z0030 lymphoid unit is listed with 8 genes; it has 18, including ADSS2.

   Use "including" consistently, or give the member counts.
9. **Ovary lead example (Z0035).** The unit combines junction proteins with the lineage transcription factors PAX8 and HNF1B. These are known lineage dependencies that could drive a one-vs-rest association alone. Report g* and its L effect for each text-highlighted lead; these are available from the BENCH outputs.
10. **Kill rule.** The kill rule is lenient: a direction survives with half of its controls, and the controls are strong, well-known dependencies. H3a rests on 2 controls and H1 on one multi-gene readout. This is stated, but the claim that the "kill rule passes" should not be read as validation of discovery power.
11. **"Most tumour-suppressor losses left their complex's dependency unchanged"** (cover letter). The supported statement is 11 of 24 events, with the ARID1A caveat. Align the wording with the main text.
12. **Discussion l. 167.** The statement that Wainberg 2021 and Chen 2023 "did not calibrate against such matched sets" should be checked against those papers' methods (e.g. Wainberg's module-level enrichment nulls) and phrased precisely.
13. **Matched-set availability.** 17 units had no matched sets and 12 had tolerances doubled (calibration_D8.json `status`). Mention this where the flag is defined.
14. **Member-resolved H1 profiles.** These were exported after H1's L read (Supplementary Note 14). They are A-only and descriptive, but state the timing in the Fig. 4c legend.
15. **Internal identifiers.** Replace them in the main text with plain descriptions, or define each at first use with a glossary table: R4, R6, R7, D1, D8, D14–D16, E1, E2, A1–A5, FZ-C, Z-IDs, "rare-context route".
16. **Format.** The HTML source comments must be stripped from the submitted text. The abstract (189 words) and main text (about 6,300 words) exceed Nature Communications limits; Genome Biology requires a structured abstract.
17. **Holm p for GLS.** The table reports "0" and Supplementary Note 10 reports "< 0.0005". With B = 2,000 and Holm ×3, the bound is p_Holm < 0.0015. Report the bound consistently.
18. **Timeline.** Part 2 design (DESIGN v1 at 15:08, v3.1 at 17:37) and the registry draft (19:27) came after the test files arrived (13:55). The registry draft also came after the Part 1 test results (18:52). Nothing suggests Part 2 used test information, but give the order explicitly in Supplementary Note 1, since the Part 2 registry was written knowing that zavg_BEF had won.

---

## 5. Statistical review

**Is BB selective FDR described correctly?** Mostly yes. Methods l. 267 describes the Benjamini–Bogomolov procedure correctly:
- family p-values by Simes;
- family selection by BH at q over M families;
- BH within each of R selected families at q·R/M;
- the controlled quantity is the expected within-family FDP averaged over selected families.

It correctly declines to claim pooled FDR, reporting 0.068 in simulation. Caveats to add:
- The BB guarantee assumes independence across families, or at least a selection rule compatible with dependence, and PRDS within families for Simes and BH. Here families are heavily dependent. One-vs-rest lineage contexts share their reference lines, and H3b contexts overlap (TP53, WGD, MSI-H, lineage covariates).
- The p-values are a hybrid. A screen at BB q < 0.1 on analytic p is followed by permutation p with GPD extrapolation for screened tests only, and these re-enter the procedure. The validity rests on the simulation in h3a_fdr_diagnosis.md, which planted *lineage* effects only. State that no simulation covered the H3b driver-family dependence structure.
- With 21 of 23 H3a families selected, BB is nearly BH within family at q ≈ 0.046. The selective layer adds little there, which is worth one sentence.

**Do the lockbox and power interpretations follow?** Partly.
- The observed H3a rate (0.525) sits at or just above the post hoc band (0.467–0.519). H3b (0.40–0.49) brackets 0.43–0.44. "About the rate power predicts" is a fair description, and the authors correctly say this is "no evidence" of excess non-replication.
- The power aid assumes that the A effects shrunk by Tweedie's formula are the true L effects. The prior is fitted to all A z-scores, which include the gene-set-generic and scale effects. Those effects would also replicate in L. Agreement with the power prediction therefore cannot separate biology from reproducible artefacts. The manuscript says this for complex-specificity (l. 125) but should say it for the overall rates too.
- In every subset the observed rate exceeds the predicted one by about 0.05–0.07: H3a beats-all 0.572 vs 0.52–0.53, and H3b beats-all 0.583 vs 0.50–0.51. This is worth noting.
- The same tool should be applied to BRIDGE-L (M3).

**Are the amendments and the driver defect handled adequately?** Yes, with minor additions.
- The defect diagnosis is convincing. The reviewer predicted the 510/50/3 split from labels alone.
- R7 keeps the read as executed as primary and strips claims from every call that the correction flips. That is the conservative choice, and reporting all three H3b rates is right.
- Table 3 and Supplementary Note 14 are exemplary in giving timing relative to outcome access.

Additional requests:
- (i) Note in the abstract or main-text Limitations that the robustness rule and both post hoc labels were applied after the read (M7).
- (ii) The response scope of the never-expressed calibration (literal vs broad reading) was "fixed before the calibrations ran but after the discovery outcomes were known". Say this in the main text, not only in Supplementary Note 12.
- (iii) Commit the run log of the lockbox loader for the next version, or state that it is permanently lost. Supplementary Note 15 item 26 currently says the call count is "inferred".

**Claims stronger than the evidence.**
- "Locked before the test data existed" (M5).
- "Label-free" in the title (M5).
- "Its advantage does not carry over to complexes", based on one null and one secondary negative bridge (M3).
- "Random gene sets … associated … as often as real complexes" and "Most complex-level context associations are not about complexes" (M1).
- "For most context dependencies, the complex is not the right unit" (M4).
- "11 of 24 genetic losses left their complex's dependency unchanged" (M6).
- "Confirmed in held-out lines" for rare-route items (M7).
- "The rate that the held-out set's power predicts" in the abstract, without "post hoc" (minor 4).

---

## 6. Figures and tables

**Fig. 1.**
- Panels a, b and e include the secondary candidate B and the post hoc PC-removed Pearson. The legend mentions neither, and panel b's B rows are unexplained. Add both to the legend.
- Panel d shows "AP without 10 largest" with no CI. The legend's "where bootstrapped" covers this, but say explicitly that this metric was not bootstrapped.
- Panel c shows zavg_BEF above the diagonal in the mid range, i.e. under-prediction, consistent with CITL 0.893. One sentence in the legend would help.

**Fig. 2.**
- Panel e's y-axis reads "20 bootstrap maps". The legend and Methods describe 80% PatientID-blocked subsamples, not bootstraps. Make these consistent.
- Panel e has no colour key for zavg vs GLS.
- Panel f: the zero line runs through the blue pair-level annotation text.
- Panel c plots GLS curves against "Edge threshold z0", but the GLS map is edge-matched rather than thresholded at z0. Clarify the x-axis meaning for GLS.

**Fig. 3.**
- Panel d has grey dots that are not explained in the legend.
- Panel e's top-axis label ("Never-expressed sets…") crowds panel d's title. Its dotted reference line at 0.05 belongs to the top axis, but visually it reads against the bottom "Discoveries" axis.
- Panel a's inset legend overlaps the data.
- Panel b shows positive ρ for Skin, Head and Neck and Bone that the text does not discuss (M2).

**Fig. 4.**
- Panel a labels FOCAD n = 46. The text, panel c and discoveries.tsv use 30 analysed loss lines (46 defined, 30 after complete-case CN covariates). Label which n is shown.
- SMARCB1 is labelled "silencing", although A2 moved it into H1 as a genetic loss. The type field is silencing, but explain this in the legend.
- Panel d sums to 26, not 24 (M6).
- The panel letter "c" overlaps the ARID1A sub-panel title.

**Fig. 5.**
- Panel a: the "H1 eligible calls 2/3 +1 NE" row omits EP300 (rare route) without saying so in the legend.
- Panel b shows ARID1A P "confirmed" with a 95% HC3 CI crossing 0. The text explains the permutation rule, but the legend should too.
- Panel c and d legends match their content.

**Tables.**
- Table 1/1b: Holm p "0" (minor 17).
- Table 2 matches the text.
- Table S7's Markdown rows list the earlier unregistered "strongest specific" set (Supplementary Tables list). This risks confusion with S9. Regenerate S7's Markdown with the registered list, or drop the earlier set.
- Table S9 labels Z0190 a lead example while the text sets it aside (M2).

**Missing panels.** A figure showing the matched-set comparison conditional on discovery (M1.3) and the DEMETER2 concordance of the 254 findings (M2e) would carry the paper's main message better than Fig. 3d. Fig. 3d is an A-side heat map that is not the registered list.

---

## 7. Reproducibility

**Data.** All inputs are public, with SHA-256 checksums in DATA.md and a check script. Weaknesses:
- The DepMap release is "26Q3 (user-stated)", with no version stamp.
- The HGNC download date is unknown.
- The licences for BioPlex, CCLE proteomics and the hu.MAP files are "unverified".

These must be resolved before submission.

**Code.**
- The repository URL, archive DOI and the public location of the Part 1 repository (commit c7f95ee) are "to be added". **This blocks submission.**
- Persisted lockbox matrices are not distributed. Because L membership is committed, readers can re-derive L values from the public release, which is acceptable.

**Seeds.** Seeds are documented for every stochastic step:
- bootstrap seed 0 for the test look;
- 4000/4100/4500 for permutations;
- 4600/4700/4800 for calibrations;
- 5000/5100/5200 for the bridge, BENCH and BRIDGE-L;
- A3 for lockbox permutations;
- `PYTHONHASHSEED=0`.

This is sufficient.

**Number checker.** I ran `scripts/check_manuscript_numbers.py manuscript/main.md`: 579/579 matched, and no tracked file changed. However, 140 of 579 were matched only by "text" (the number appears anywhere in a cited Markdown or YAML file). That is a weak check. The supplement and cover letter are not checked. Extend the checker to supplementary.md and cover_letter.md, and tighten text matches to cited sections.

**Independent spot-checks (my own, from results files).**
1. AP 0.333 vs 0.212; Δ +0.121 (+0.079, +0.159); relative gain 57%; Kish neff 21.3; 238 clusters; C3 positives 2,309. `results/test_eval/test_results.json` gives 0.33290, 0.21155, Δ 0.12135 [0.0794, 0.1586], rel_gain 0.574, kish 21.287, n_clusters 238 and C3 pos 2,309. **Match.**
2. n_U 1,345, n_L 338, n_A 1,007, no straddling patients. `results/part2/lockbox/manifest.json` agrees. **Match.**
3. BioPlex bridge 0.157 vs 0.182, Δ −0.025 (−0.075, +0.031). `bridge.json` `primary` gives 0.15703, 0.18223, −0.02520 [−0.0745, +0.0308]. **Match.**
4. BRIDGE-L −0.065 (−0.107, −0.025), hits 1,267/1,326, matched coverage −0.074 (−0.122, −0.026). `bridge_primary_vs_corrected_R7.json` agrees: −0.06481 [−0.1069, −0.0246]; −0.07422 [−0.1216, −0.0256]. **Match.**
5. Matched sets 28.1% vs 28.1% (H3a), 11.6% vs 12.0% (H3b); 502/1,145 and 103/161 beat all. `calibration_D8.json` gives 0.28140/0.28077, 0.11600/0.11989, 502 and 103. **Match.**
6. Additional checks:
   - H3b rates 60/122 = 0.4918, 60/150 = 0.40, corrected 69/150 = 0.46 (CMP);
   - BENCH high tercile +0.240 (−0.033, +0.362);
   - H1 11/24;
   - strongest 207/45/2;
   - 317 LOSO failures;
   - power bands.

   **All match.**

---

## 8. Circularity check

- **STRING.** Not used. No STRING file is in DATA.md, and a code search of `src/part2` and the Part 2 scripts found no STRING, CORUM, Complex Portal or hu.MAP access. The firewall F2 is stated in `src/part2/__init__.py`, and construction scripts declare their inputs (`scripts/31_score_blocks.py`, `36_build_maps.py`).
- **hu.MAP.** Part 2 reads neither hu.MAP files nor test results. **Indirect route:** the locked zavg_BEF configuration was selected on hu.MAP training labels derived from Complex Portal, so curated complex information enters the map via method selection (M5, minor 1–2). Part 1 validation uses hu.MAP test pairs (curated-derived) by design, so l. 217 must be scoped to Part 2.
- **Curated complexes.** CORUM and Complex Portal were not obtained. HGNC gene groups, which partly encode curated complexes, are read only in `scripts/38_describe_maps.py` for a column prefixed `POSTHOC_`. Positive controls use literature knowledge, which is appropriate for controls.
- **BioPlex and proteomics.** Both are used only for validation and tiers after the map commit. The exception is the disclosed read of proteomics by the context builder for the LoF functional check (gene-level, before the map commit). BioPlex tiers are then used to *rank* strongest findings and to choose lead examples (T1). This is not circular, but validation data steer which findings are highlighted. Say so.
- **Post hoc annotations.** They are clearly marked in Table S9 (column names; D16), and Supplementary Note 14 states that "labels select and rank nothing". In practice the post hoc "lead example" label (unregistered |z| ≥ 2 threshold) selects the findings shown in the Results, Discussion and cover letter. Describe it as a post hoc selection for presentation, not as purely descriptive.

---

## 9. Do the title and abstract match the evidence?

Partly. The abstract is admirably balanced in reporting negative results, but:
- "Label-free" and "locked before the test data existed" are inaccurate (M5).
- The abstract omits the null macro AP and the concentration of the gain in pair-rich clusters (M5).
- "Random gene sets" misdescribes the comparators (M1).
- "The rate that the held-out set's power predicts" omits "post hoc" (minor 4).
- "11 of 24 genetic losses left their complex's dependency unchanged" counts ARID1A (M6).
- The 254 strongest findings are presented without saying that the rule was applied after the lockbox read, and that 56 of them are labelled possible artefacts (M7).

The title's second clause rests on one null primary bridge and one negative secondary bridge (M3).

A title such as the following would be accurate: "A co-essentiality score selected on training complexes wins a single held-out test, but yields neither better complexes nor more replicable complex-level dependencies in DepMap".

---

## 10. Items blocking submission (independent of scientific revision)

1. Author list (`[[AUTHORS]]`), competing interests and the cover-letter "[Authors to confirm]" placeholder. A human corresponding author must take responsibility for an agent-run analysis.
2. Repository URL, archive DOI and the public availability of the Part 1 repository at c7f95ee.
3. Unverified licences (BioPlex, CCLE proteomics, hu.MAP files), the unknown HGNC download date and the "user-stated" DepMap version.
4. HTML source comments stripped; abstract and main-text length brought within journal limits.
