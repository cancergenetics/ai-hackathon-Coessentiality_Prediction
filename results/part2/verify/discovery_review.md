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
