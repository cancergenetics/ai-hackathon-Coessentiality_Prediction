# Discovery runs H3a / H3b / H1: independent review

Reviewer: independent agent (did not produce the runs). Reviewed state: `446395d` (branch `review-discovery`, cut from
the producer branch; the producer's later commits `c3c4422` (erratum E1) and `b3e8043` (h1 fix R2) were read but not
reviewed). Scripts: `scripts/verify/discovery_lib.py`, `discovery_review.py`, `discovery_full.py`. They do not import
`src/part2`.

## Verdict: ACCEPT-WITH-NOTES

The process is clean. Every H3a and H3b analytic statistic (15,781 tests) and the 5 H1 discovered tests reproduce
from independent code. The BB / BH discovery sets match exactly. Two items must be disclosed before any claim: the
GPD-tail dependence (M1) and the fragile H1 PTEN call (M2). Neither changes the registered discovery sets.

## MUST (disclose; no re-run)

**M1. GPD-tail sensitivity must sit next to the discoveries, and extreme p_perm values must not be reported as numbers.**

| | H3a | H3b | H1 |
|---|---|---|---|
| Discoveries (registered rule) | 1,145 | 161 | 5 |
| p = max(p_perm, p_analytic), BB / BH re-run | 1,113 (−32) | 130 (−31) | 5 (−0) |
| p = max(p_perm, 1/10,001), i.e. empirical only, no GPD | 1,145 (−0) | 136 (−25; R 13 → 10) | 5 (−0) |
| Screened tests with p_perm < 1e-4 / below analytic by ≥ 10× | 346 / 73 | 65 / 18 | 4 / 1 |
| Smallest p_perm among discoveries | 1.4e-24 | 4.6e-57 | 1.9e-7 |

- **Losses under max(p_perm, p_analytic).**
  - H3a: Myeloid 7, Esophagus/Stomach 6, Uterus 6, Soft Tissue 4, Liver 2, Pleura 2, Skin 2, Head and Neck 1,
    Ovary 1, PNS 1.
  - H3b: WGD 16 of 26, TP53_altered 7, EWSR1 fusion 5, KRAS 2, LoF:PTEN 1.
- **Family selection rests on the GPD tail.** Under the empirical floor, the families hotspot:NRAS, sex_female and
  LoF:SMARCA4 are deselected. Each of them holds one discovery.
- **Inflated precision.** Examples: Cervix (n = 13) Z0037 has analytic p 4.5e-8 and p_perm 1.3e-23. Top H3a hits
  have analytic p of 1e-8 to 1e-17 against p_perm near 1e-22.
- **Required.**
  - Report p_perm < 1e-4 as "< 1e-4 (GPD)".
  - Carry both sensitivity columns in the discovery tables.
  - Label the 31 H3b hits and the 3 GPD-dependent families. The rule is pre-registered, so this is a disclosure
    item, not a protocol breach.

**M2. The H1 Z0130|PTEN COLLAPSE call is fragile and must be labelled as such.**
- It is rejected at the run's m = 24 (p_final 0.0091, q 0.044). At the registered `family_m` 28, BH rejects 4 and
  PTEN is lost. m fell to 24 because EP300|N and SMAD4|N have |N'| = 1 and MAP2K4 has 11 loss models. That follows
  the registry rules, and runtime.json discloses it, but the call depends on it.
- κ_N = −0.057 with 90% CI [−0.092, −0.021]. That CI is also inside the TOST margin, so the event is NEUTRAL as well.
- Z0158|TP53 RESIDUAL SENSITISATION is in the same position: 90% CI [0.038, 0.075].
- The registry does not set the precedence between a non-neutral class and NEUTRAL. State this, and report both
  conditions for TP53 and PTEN.

## Checks

**1. Process: pass.**
- **Commit order.** Each intent commit precedes its run, and each run precedes its results commit:

  | Direction | Intent commit | Run start | Run end | Results commit |
  |---|---|---|---|---|
  | H3a | 228077a, 00:20:31 | 00:20:39 | 00:28:10 | b0dbd83, 00:28:51 |
  | H3b | 388924a, 00:29:56 | 00:30:03 | 00:33:22 | 5e1c1d1, 00:33:49 |
  | H1 | ef41ddb, 00:34:45 | 00:34:53 | 00:35:46 | d12c627, 00:36:25 |

  The intent commits add only INTENT.md and dry_check.json.
- **Commands.** The runtime.json command equals the intent command, wrapped in `/usr/bin/time -l` as declared.
- **Trial ledger.** trials.jsonl has exactly one discovery line per direction, with times 00:28:09, 00:33:22 and
  00:35:10.
- **No rerun traces.** For each `data_cache/part2/assoc/<ID>/discovery.parquet`, the birth time equals its
  modification time and the trial time. The branch reflog is linear (no reset or amend). No other worktree or main
  holds a discovery line. Runner, engine and registry are unchanged from `beda09f` to `446395d`.
- **No L model used.** A ∩ L = ∅ (1,007 / 338). max n is 1,007 (H3a, H3b) and 961 (H1). My A-only recompute
  reproduces n exactly for every H3a and H3b test and for the H1 tests.
- **Runner logic.** The runner is thin, with minor analysis-adjacent code:
  - the H1 matched-null essentiality vector (`E_A.clip(±3).mean()`);
  - a re-implementation of the pseudo-context p-values for the calibration BB report. It is reporting only and is
    guarded by an equality check against the engine.
  - CI and tier/locus annotation;
  - covariate assembly for the SMARCB1 sensitivity.

  None of these feeds discovery.

**2. Independent recompute: pass.**
- **30-test sample.** 10 random H3a discoveries, 10 random H3b discoveries, 5 random non-discoveries and the 5 H1
  discovered tests all match. Maximum relative differences are 3.0e-6 (β), 3.6e-6 (SE), 3.2e-6 (df) and 3.0e-6
  (p), the 6-significant-digit rounding of the TSV. n matches exactly.
- **H1 alignment.** The H1 N' sizes match (8, 3, 4), confirming ALIGN = +1 (y_N = mean sign(b_j)·effect_j), as
  E1 later corrected in the text.
- **Full recompute.** I recomputed all 8,441 H3a and all 7,340 H3b tests with my own HC3 and Bell–McCaffrey df. All
  match, including the registered removals: paralogs, the ±10 Mb windows and X–Y gametologs.
- **Two bugs of mine, fixed.** My first pass mismatched 5 H3b tests. Both causes were bugs in my code, not the
  producer's:
  - an Ensembl gene ID with two Entrez IDs (EIF1AX);
  - patch and PAR coordinates (RNF139, IL9R).
- **BB on my own analytic p.** I get 1,248 / 138 rejections and 1,666 / 190 screened tests, the same sets as
  committed. The screen equals BB run at q = 0.1.
- **BB on the committed p_final.** Discovery sets match exactly: H3a 1,145 (R = 21; Biliary Tract and Prostate not
  selected); H3b 161 (R = 13). H1 pooled BH on 24 tests gives the same 5.
- **TOST.** For STAG2, KDM6A and APC, the 90% CIs and the NEUTRAL calls match.

**3. Sanity.**
- **Lineage pattern.** H3a discoveries per lineage track lineage n (Spearman 0.70). Lymphoid has 197/367 at
  n = 73. Head and Neck has 126 at n = 58; it is homogeneous squamous, against heterogeneous Lung with 48 at n = 93.
  Biliary Tract and Prostate have 0. Top lymphoid hits are biologically coherent:
  - adhesion/integrin units (ILK, ITGAV, BCAR1) are less needed in suspension lines;
  - purine/pyrimidine synthesis is more needed.
- **H3b plausibility.** Hits include TP53_altered → MDM2/MDM4/PPM1D (+0.29), MSI-H → WRN, BRAF hotspot →
  BRAF/MAP2K1/DUSP4 and EWS-FLI1 → its own unit.
- **Note on scale.** In Lymphoid, β correlates negatively with complex essentiality (ρ = −0.24: essential
  complexes look less essential). In Head and Neck, Skin and Bone the correlation is +0.2. That points to a partial
  Chronos dynamic-range component. The pre-registered within_line_rank and NNMD multiverse checks should be read
  first for these lineages. 325 of the 1,145 H3a discoveries have |β| < 0.05 Chronos.
- **Locus share.** 203 of 1,145 H3a discoveries are in locus-flagged units (17.7%), and 40 of 161 in H3b (24.8%).
  Locus units are 104 of 367 (28.3%), so they are not over-represented. The second-ranked H3a hit (Cervix Z0037,
  the 11q13 amplicon) is locus-flagged.
- **Influence.** I removed the k most supportive models by exact dfbeta and refitted with analytic p.

  | k removed | H3a (1,145) p ≥ 0.05 | H3b (161) p ≥ 0.05 |
  |---|---|---|
  | 1 | 0 | 0 |
  | 2 | 9 | 0 |
  | 3 | 47 | 0 |
  | 5 | 247 | 6 |

  - The 9 at k = 2 are all Lymphoid, with analytic p 0.015 to 0.024; after removal p is 0.050 to 0.077.
  - The 6 at k = 5 are all TP53_altered, with |β| < 0.031.
  - No sign flips at any k.
  - The top 5 H3a hits stay at p ≤ 3e-5 after removing 5 models.
  - H1: ARID1A, EP300, FOCAD and TP53 stay at p ≤ 7e-4 after removing 3. PTEN gives 0.017, 0.025 and 0.038 for
    k = 1, 2, 3.

**4. GPD tail.** See M1. Under p = max(p_perm, p_analytic), 32 H3a, 31 H3b and 0 H1 discoveries are lost, 63 in total.

## Notes

- **N1. H1 calls rest on few members near the |b_j| ≥ 0.1 threshold.**
  - TP53 (attached to Z0158): N' is CCDC6, PABIR1, SHLD1 and SHLD2. CCDC6 has b = −0.110. The strongest
    TP53-dependent members are left out by the threshold: TP53BP1 (b = +0.02, β = −0.22) and ATM (b = +0.07,
    β = −0.17).
  - PTEN: N' is SZT2, TSC2 and ZFP36L2, with baselines of −0.109 and −0.111 among them.
  - FOCAD: driven by PELO (aligned +0.81) and HBS1L (+0.33), a coherent signal.
  - Report the member profiles with the class.
- **N2. MAP2K4.** Both MAP2K4 tests are not evaluable (n_context 11), yet classes.tsv labels them INCONCLUSIVE ×2,
  and so do the summary.json class_counts. The producer's later fix R2 (`b3e8043`) addresses classify. The committed
  H1 tables should be relabelled NOT EVALUABLE.
- **N3. EP300.** n_context = 12, exactly at the evaluability floor. p_perm is 2e-7 (GPD) against analytic 3e-5.
- **N4. Design conventions.** The lineage_FE reference level (first sorted), the merging of indicators with fewer
  than 2 models, greedy collinearity dropping and the Bell–McCaffrey df come from engine docstrings, not the
  registry. I mirrored them in independent code; they reproduce exactly.
- **N5. Post-run commit E1.** It corrects the registry text to the semantics actually run. My recompute confirms the
  run used sign(b_j).

## Not verified

- Permutation p-values, including the Freedman–Lane draws and the GPD fits. I used the committed p_perm.
- Pseudo-context calibration.
- The H1 matched null, P-vs-N contrast and breadth conditions.
- The SMARCB1 sensitivity.
- The exact PART1_DIR at run time, which is elided in runtime.json.
- Whether an earlier real run was made and deleted. No trace was found in birth times, the ledger or the reflog.
  Commit times are self-reported but consistent with ledger and file times.

## Reproduce

```
export PART1_DIR=".../coessentiality_protein complexes"; export PYTHONHASHSEED=0
.venv/bin/python scripts/verify/discovery_full.py   --out <scratch>/dr   # all 15,781 H3a/H3b tests, ~30 s
.venv/bin/python scripts/verify/discovery_review.py --out <scratch>/dr   # sample, FDR, TOST, GPD, sanity, ~1 min
```

---

# Review 2: H1S, H2, calibrations, GPD sensitivity, errata

Reviewer: the same independent agent (did not produce the runs). Reviewed state: producer `f4e99dd`, merged into
`review-discovery`. New script: `scripts/verify/discovery_review2.py` (does not import `src/part2`; it reuses
`discovery_lib.py` / `discovery_review.py`). PI rulings R1–R5: RESEARCH_TREE.md D8 (main `3f5e8fd`).

## Verdict: ACCEPT-WITH-NOTES

The process is clean. Every quantity I recomputed reproduces exactly:
- 10 analytic tests;
- the H1S BH and H2 stage-wise discovery sets;
- the GPD sensitivity counts;
- the H3a matched-set and never-expressed rates (RNG-mirrored to |Δt| ≤ 5e-5);
- the H2 random-split rejections.

The R2 fix cannot change an H1 discovery, and E1/E2 touch only text and reporting columns. The substantive finding
comes from the calibrations themselves:
- **H3a.** Lineage signal is about as common in matched random sets as in the real complexes.
- **H2.** Most sub-module divergence is reproduced by random splits.

The READMEs say so, but claims must be scoped accordingly (M1).

## MUST (disclosure and reporting; no re-run)

**M1. Every H3a / H3b / H2 claim must carry the calibration background, and specificity must be shown per hit.**
- **H3a matched random sets.** p < 0.05 in 28.1% of random-set tests against 28.1% of real tests. That makes 1,145
  discoveries "lineage moves this gene set", not "lineage moves this complex specifically".
  - Across all 8,050 evaluable real tests with 10 sets, the real unit beats all of its matched sets 13.3% of the
    time. Exchangeable sets would give 9.1% (1/11). So unit-specificity is real but modest.
  - The same holds for H3b (13.3%).
  - Add a per-discovery `exceeds_all_matched` flag to `discoveries_report.tsv` (H3a 502 of 1,145, 62 without sets;
    H3b 103 of 161). Restrict "complex-specific" wording to flagged hits.
- **H3a never-expressed sets.** 13.5% against 5% expected, with 13 BB rejections. The rate is concentrated in
  melanoma (0.33), Skin (0.28) and Lymphoid (0.27). These are the same lineages where the review-1 scale effect
  appeared.
  - Report the per-lineage never-expressed rate next to each lineage's discovery count.
  - Read those lineages together with the pre-declared ScreenNNMD + CasActivity sensitivity (D8).
- **H2 random splits.** Mean 78.3 stage-wise rejections per replicate (range 66–90), against 109 real. A
  sub-module-specific claim needs each discovery's rank against its own 20 splits; that is not yet reported.
  Without it, H2 discoveries should be described as within-parent heterogeneity.

**M2. State the scope of the registered response explicitly.**
- The registry `calibration` text lists pseudo-contexts, matched random sets and never-expressed sets in one
  sentence, then says "if *that fraction* exceeds 0.075 …". The producer reads "that fraction" as the
  pseudo-context fraction only, plus the CN segment null for H3b CN contexts. That reading is grammatically right,
  and it was fixed in `1e674bd` before the calibration ran, but after the discovery outcomes.
- Under the broader reading, the H3a never-expressed rate (0.135) would have fired the response.
- Firing would not have been conservative: it re-tests p < 0.01 by permutation and re-runs BB. Permuting context
  labels does not remove a screen-level shift that is shared by the gene set, so it would not address what the
  never-expressed sets detect.
- Add one sentence to D8 recording this reading and its consequence.

## Checks

**1. Process: pass.**
- **Commit order** (from `git log` and file birth times in the producer worktree):

  | Item | Intent / pre-declaration | Run start | Outputs born | Results commit |
  |---|---|---|---|---|
  | H1S | `b865595` 00:49:17 | 00:49:24 | 00:56:42 (trial 00:50:16) | `06c8cd1` 00:57:05 |
  | H2 | `5ea25ca` 00:57:57 | 00:58:04 | 00:59:01 (trial 00:59:01) | `1274e40` 00:59:25 |
  | GPD sensitivity | `7b754f8` 00:59:34 | — | 00:59:43–44 | `7e33c6e` 00:59:52 |
  | Calibrations | code `44f30ca` 01:02:23, intent `1e674bd` 01:02:49 | H3a log 01:02:57; H3b 01:03:54; H2 01:05:37 | 01:03:34 / 01:05:22 / 01:06:11 | `fd582a6` / `4fb83d7` / `1a2a445` |

- **Ledger.** trials.jsonl has exactly one discovery line each for H1S and H2, and the same lines appear on every
  branch that holds them. Calibrations and GPD write no trial line; neither intent promises one.
- **No reruns.**
  - `data_cache/part2/assoc/{H1S,H2}/discovery.parquet` have birth = modification = trial time.
  - Each calibration output has birth = modification, and the runner refuses a second run.
  - The producer reflog is linear.
  - `scripts/44_run_direction.py` is unchanged since `beda09f`.
- **A only.** No lockbox model ID appears in any file changed `bff6cad..f4e99dd`. max n is 962 (H1S) and 888 (H2).
  The calibration runner uses the A-hash-checked environment.
- **Minor.** The intents' "planned start" (00:52 / 01:02) is later than the actual start. This is immaterial
  because each intent commit precedes its run.
- **GPD independence.**
  - The producer's `scripts/45_gpd_sensitivity.py` (file born 00:50:14) re-runs the engine `assoc._fdr` on the
    SHA-verified parquet tables. My numbers came from my own BB / BH on `tests.tsv.gz`. No code is shared.
  - The quantity is a deterministic function of committed data, and the PI fixed the rule (R3). Knowing my counts
    could not tune anything.
  - Agreement is exact: H3a 1,113, H3b 130, H1 5. I also recomputed the producer's new numbers: H1S 29 (losing
    ATP8B1|N and ZEB2|N) and H2 107.
  - The "<1e-4" label is equivalent to "the GPD branch was used". The engine fits the GPD only when there are 0
    exceedances and returns min(p_GPD, 1/10,001).
- **E1.**
  - The registry diff (`c3c4422`) is the one `score_rule` line plus the E1 `amendments` record, in both H1 and H1S.
  - The `h1.py` change is docstring only.
  - `ALIGN = +1.0` and the code that uses it are in `201d8bf` (2026-10-28 20:09), an ancestor of the H1 run
    `d12c627`. The sign test is in `040e0df` (20:16), an ancestor of the H1 intent.
  - My review-1 N' sizes already confirmed that the run used sign(b_j).

**2. R2 fix: pass, with a note on scope.**
- **What changed.** `b3e8043` changes two things:
  - `matched_null` returns "null unavailable" when no null fit succeeded. This is the empty case.
  - `classify` classes every axis whose counted test is absent or not evaluable as NOT EVALUABLE, with no fit and no
    null draw. This part is broader than "empty null": it relabels any non-evaluable axis, which previously got
    INCONCLUSIVE from a fit.
- **Effect on H1.** Both changes are confined to tests that never enter BH. Matched-null seeds are per event
  (7000 + event index), so skipping a draw shifts no other event's stream. For H1 the only affected event is
  MAP2K4 (n_context 11), which goes from INCONCLUSIVE to NOT EVALUABLE. This is applied in the reporting layer
  (E2). No discovery, q or other class can change.
- **Hash check.** I re-ran `scripts/verify/45_h1_fix_equivalence.py` on clean `git archive` trees of `c3c4422` and
  `b3e8043`. All 17 hashes (16 runner files + e2e JSON) are identical before and after, and equal the committed
  values. The check is vacuous for the changed branch, however: the synthetic project has no present-but-non-evaluable
  axis and no empty null. The producer's new unit test is what exercises it.
- **Tests.** `test_part2_calib`, `test_part2_h1h2`, `test_registry_entries` and `test_part2_fzx` pass (47).

**3. Calibrations: implemented as registered (operational choices fixed before running); recomputed exactly.**
- **Registry fields.**
  - Sets, matching variables and seeds (4600 / 4700 / 4800) follow the `calibration` fields.
  - R = 10 (H3) and R = 20 (H2), the H1 matched-null tolerances (±0.1 essentiality, ±1.0 z coherence, doubled
    once), 100 never-expressed sets and the ±25% segment rule are the producer's choices. They were declared in
    the intent before running.
  - Matched sets are subsets of other frozen complexes. Coherence matching makes them complex-like, so this is a
    complex-vs-other-complexes background (see M1).
- **Responses.**
  - Only the CN segment null carries a response. It was computed and compared mechanically: pooled 0.060 ≤ 0.075,
    and per context 0.045 / 0.061, so it did not fire.
  - MYC_amp has no qualifying segment, so it has no null.
  - The 9p21 null rests on 3 segments, one of which is Yp11.2 (RPS4Y1) hom-del, a LoY proxy.
  - MYCN contributes 13,499 of the 14,599 pooled tests.
- **H3a matched sets, subsample.** I took 30 random units and mirrored the RNG with my own sampler.
  - Status was 28 ok, 1 tolerances doubled and 1 unavailable, giving 290 sets.
  - All 6,670 committed tests matched (max |Δt| 4.5e-5, max relative Δp 4.7e-6). The p < 0.05 rate is 28.9% in
    both. Real tests in the same units: 25.5%.
  - Exceed-all: 41 of 80 discoveries, both mine and committed.
  - An independent draw on the same units (own stream) gives 30.3% and 40 of 80.
  - The full 502 / 1,145 (62 without sets) reproduces from the committed tables.
- **Never-expressed sets.**
  - My own read of the raw expression CSV gives 958 RNA-profiled A models and a pool of exactly 558 genes.
  - The mirrored sets reproduce all 2,300 tests (max |Δt| 5e-6): 13.48% and 13 BB rejections.
  - An independent draw gives 12.8% and 8 BB rejections.
- **H2 random splits.**
  - The mirrored draw reproduces all 9,040 tests (max |Δt| 4.9e-5): rate 0.375; stage-wise rejections mean 78.3,
    range 66–90.
  - An independent stream gives 0.376 and mean 81.6 (58–103).

**4. H1S / H2: pass.**
- **Analytic tests.** 5 H1S tests (3 discoveries, 2 not) and 5 H2 tests (3 discoveries, 2 not; lineage and
  TP53_altered contexts) match. Maximum relative differences are 2.4e-6 (β), 1.6e-6 (SE), 1.3e-6 (df) and 2.0e-6
  (p), and n matches exactly.
  - H1S uses ALIGN = +1, reference-model baselines and the event covariates.
  - H2 uses Y_s = S_s − (a + b S_sib), with an OLS fit over A after the context's removals, and the context
    covariates + NNMD + CasActivity.
- **H2 stage-wise, from the committed p.** There are 210 (context, parent) families, at levels q R_x / m_x =
  0.0067824 (H3a parents) and 0.0010967 (H3b parents).
  - On p_final: 109 discoveries, identical set.
  - On analytic p: 112 rejections, identical.
- **H1S BH.** m = 293 gives 31 discoveries, identical set. Padding to the registered m = 338 keeps all 31.
- **H1S m.** Of 317 events, all have |N| ≥ 2 and 21 have |P| ≥ 1, so m = 338 as registered.
  - My |N'| equals the committed value for all 317 events. 17 have |N'| < 2, and they are exactly the 17 N units
    absent from the test table, giving 321 tests. This is the rule already applied in H1 (EP300, SMAD4).
  - The remaining 28 are leverage-1 designs, for which HC3 is undefined. They were listed from labels only before
    the run (`degenerate_designs.tsv`, `446395d`). Result: 293.

**5. E2: pass.**
- Both `classes.tsv` files gain three appended columns (`label_final`, `label_final_paralog`, `headline_eligible`).
  All 46 original columns are cell-identical to the run commits (`d12c627`, `06c8cd1`).
- `tests.tsv.gz`, `discoveries.tsv`, `summary.json` and the trial lines are untouched.
- The rule is applied exactly. H1: PTEN and TP53 become DETECTABLE-WITHIN-MARGIN, and MAP2K4 becomes NOT
  EVALUABLE. H1S: 4 N-axis relabels (ERBB3, PRR15L, B3GNT5, CST6), which equals the count of non-neutral axes with
  N_tost; there are 0 P-axis relabels. `headline_eligible` is False for all H1S events.

## Notes

- **N1.** H2 `summary.json` carries `standing: confirmatory` and `headline: true`. These come from the runner
  defaults, because the H2 entry has no `standing` field. This contradicts its exploratory status, which follows
  from its 0 evaluable controls (README, intent, `kill_rule_verdict`). Correct it in the reporting layer so no
  downstream tally reads it as headline.
- **N2.** E2 writes its columns into the run's own `classes.tsv`. The original cells are preserved, but a separate
  reporting file would keep run outputs byte-identical to their commits.
- **N3.** R3 asks tables to show "<1e-4 (GPD)". `discoveries_report.tsv` writes "<1e-4". The meaning is the same;
  align the string.
- **N4.** The H2 pseudo-contexts give 6 stage-wise rejections, against about 0 expected (pooled p < 0.05 = 0.050).
  This is small relative to 109, but it is non-zero and should be reported as is (it is).
- **N5.** Review-1 N2 (MAP2K4) is resolved by E2. M2 (PTEN) is resolved by R5/E2, and the README reports the m
  dependence (q 0.044 at m = 24, 0.051 at m = 28).

## Not verified (part 2)

- The H1S matched null, contrast and breadth calls.
- H1S / H2 permutation p.
- The H3b matched / never-expressed / segment tests: I checked the committed tables and the response arithmetic
  only.
- The H2 pseudo-context replay.

## Reproduce (part 2)

```
export PART1_DIR=".../coessentiality_protein complexes"; export PYTHONHASHSEED=0
.venv/bin/python scripts/verify/discovery_review2.py --out <scratch>/dr2      # ~85 s
# R2: git archive c3c4422 / b3e8043 into clean trees, copy scripts/verify/45_h1_fix_equivalence.py into the first,
# run it in each, compare the JSONs with each other and with results/part2/directions/h1_fix_equivalence.json
```
