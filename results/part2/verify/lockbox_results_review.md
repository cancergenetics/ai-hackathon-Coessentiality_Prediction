# Lockbox read results: independent review

Date: 2026-10-29. Reviewed: producer branch `worktree-agent-a3ec5565842d4e968` at `8cdd176` (intent `ed065f4`, summary
`5bc9444`, R6 `38e301a`/`4d48039`/`68f0d3f`, R7 `b0ae055`/`29ab061`/`94b42c3`/`8cdd176`). Review branch `review-lbres`
starts at `8cdd176` and adds only this note, three scripts under `scripts/verify/` and their outputs under
`results/part2/verify/lockbox_results_review/`. No producer file is modified. No ledger line was written.

My code is new and does not call `src/part2/confirm.py`, `robust.py`, `stats.py`, `assoc._design`, `bench.py` or
`scripts/55_lockbox_r7_report.py`:
- `scripts/verify/lockbox_results_review.py`: covariate one-hot assembly, FE-merge rule, design span, HC3 via
  Frisch–Waugh, Satterthwaite df (exact moments of the HC3 quadratic form), the PatientID-blocked Freedman–Lane
  permutation, the confirmation rule, the rare-context route, and a label-only span audit of every fit.
- `scripts/verify/lockbox_results_review_tables.py`: rates, defect-sensitive list, E2, BRIDGE-L bootstrap and the
  power aid, all from committed per-hit files.
- `scripts/verify/lockbox_results_review_bench.py`: BENCH per-pair Δ and the bootstrap.

Reused from the project, and not re-derived: the registry's covariate-name lists (`covariates.covariates_for`), the
frozen score constructions (`Spec.score_L` / `score_A` / `score_D2`), and the label tables (`World.ctx/loss/cov`).

## Verdict: PASS, with record and interpretation fixes (MUST-1 to MUST-4)

- The read was executed as intended:
  - one ledger line per entry, in registry order;
  - no re-read, and every persisted hash matches;
  - the relook ran only after its intent, and only through the resume path.
- Every confirmation call I recomputed agrees with the producer, in both versions: 54 of 54. These are 20 H3a, 17 H3b
  (15 random plus 2 targeted), 5 H1 and 12 control units from 5 controls.
- Every rate in the comparison table is reproduced exactly.
- The defect-sensitive list is reproduced exactly: 53 calls. E2 labels are never promoted.
- The defect diagnosis is correct, at full scale. From labels alone, my span audit predicts exactly which fits are
  affected: 510 changed, 50 failed by the defect and 3 genuine failures, out of 3,014 fitted and 53 failed fits.
- BRIDGE-L and the BENCH point estimates are reproduced. The BENCH bootstrap CIs are reproduced exactly from L.
- One new finding: 16.5% of BENCH bootstrap fits fail silently (MUST-3).
- What needs fixing is how some numbers are reported and read, not the numbers themselves.

## How I accessed L (exact)

- **L matrices.** Only `assoc.resume_lockbox_read(aid, proj, relook="R6")`. This is the function the driver's R6/R7
  relook calls. It requires the committed ledger line, the committed SHA-256 (it re-hashes the file) and the committed
  R6 note plus intent. It returns the matrix into memory.
  - I called it for all 30 entries: once per entry for the span audit (only the row index was used), plus a second
    time for the 7 entries I refitted. For BENCH, the first call only re-checks the hash; a separate call in the
    BENCH script supplies the values for the Δ and bootstrap recompute.
  - Every call passed the hash check, so **all 30 persisted reads match the committed hashes**.
  - Every call returned 338 rows, none of them A models, and n_genes equals the ledger's n_genes in all 30.
  - Log: `l_access_log.tsv`.
- **Loader guard.**
  - `Project.outcome_loader` was replaced by a guard that refuses any non-A request and logs every call. Log:
    `run_meta.json`.
  - The guard recorded exactly 2 calls: CRISPR for the 1,007 A models, then DEMETER2 for the 1,007 A models (the
    EP300 rare route).
  - **0 loader calls happened inside any resume call.**
- **A matrix.** The A matrix was taken from a read-only memmap of A rows only (`data.load_gene_effect` was patched), so
  no CRISPRGeneEffect L row was materialised.
  - Disclosure: the project's DEMETER2 loader parses the whole DEMETER2 CSV before subsetting to A. That is
    pre-existing behaviour (pre-read review N9). It is not a CRISPR L read.
- **What I did not open.** I never opened `data_cache/part2/lockbox_reads/*` or `lockbox_models.txt` directly. I only
  `stat`ed the persisted files for their mtimes.
- **Written outputs.** These hold only per-hit L statistics for the sampled hits and BENCH pairs, the same fields the
  committed `confirmation.tsv` and `bench_pairs.tsv` already carry.

## 1. Process

| Check | Result |
|---|---|
| The intent precedes every read | `ed065f4` (06:47:06Z) is an ancestor of every target-list commit. Each ledger line's `head` is its target-list commit, which is the direct parent of the read commit, and each target-list commit touches only its own target file. `intent_hits.tsv.gz` SHA-256 = `bba0041e…` (as stated in INTENT). |
| One ledger line per entry; no re-reads | 30 lines and 30 commits touch the ledger, across **all** branches (`git log --all`). Each adds exactly one line, with its own analysis_id. Nothing touches the ledger or `lockbox_reads.jsonl` after `5bc9444`. Every `summary.json` has `"resumed": false`. |
| Registry order | The ledger order equals `confirm_order` (controls by registered seed, then H3a, H3b, H1, H1S, H2), then BENCH, then BRIDGE-L. BRIDGE-L's GLS discovery (`5e5f186`, 09:22:53) precedes its target list and read (09:23:01–05). |
| The loader was called exactly 30 times | **Inferred, not observed.** The run log was not committed. The evidence is mechanical: `read_lockbox` commits the ledger line before it loads and refuses a second line; `_persist_read` refuses an existing file; there are 30 lines, 30 persisted reads with hashes committed 0–1 s after their line, and no resume. File mtimes equal the read times (06:50:05 for CTRL_PC1 … 09:23:06 for BRIDGE-L), and the last mtime is 09:23:06. The 2 h 14 m gap between BENCH's read and its output matches its `runtime_s` 8,193 (B = 2,000 bootstrap). |
| The persisted hashes match the ledger | Yes, for all 30 (re-hashed by the resume path). `n_models` = 338 in every record. |
| The driver code was unchanged during the read | `git diff ed065f4 5bc9444 -- src scripts` is empty. The only changes in the read window are status flips in the registry YAMLs and test edits in `5bc9444` itself. |
| R6/R7 ordering | R6 note and intent at `4d48039` (10:19:44). The fix `38e301a` landed in the same second. The check output `68f0d3f` is at 10:26. The R7 addendum `b0ae055` (10:41:49) precedes the driver's `corrected` mode `29ab061` (10:42:23) and the output `94b42c3` (10:49). The relook path refuses to run without the committed note and intent, so the ordering is enforced by code, not only by timestamps. |
| The relook never called the loader | `resume_lockbox_read` contains no loader call; my guarded run of the same path recorded 0 calls. No ledger or hash-sidecar change after `5bc9444`. |

## 2. Independent recompute (registry rule only)

Rule applied:
- **Confirmed** iff the fit is OK, ≥ 5 context and ≥ 5 reference L models, same sign as A, and one-sided p < 0.05.
  The p-value is from HC3 with Satterthwaite df.
- **Permutation:** exact PatientID-blocked Freedman–Lane, B = 10,000, my own RNG. It replaces the HC3 p when there are
  < 20 context models. It runs within lineage for non-lineage contexts, and blocked but unrestricted for lineage-like
  ones (as in the pre-read review, MUST-3c).
- **Both designs:**
  - *defect*: one-hot assembled over all models, then subset to L;
  - *corrected*: one-hot assembled over the L rows read.
- Output: `recompute_calls.tsv`.

| Set | Calls agree, primary | Calls agree, corrected | β, SE, df vs committed |
|---|---|---|---|
| H3a, 20 random (seed 20261029) | 20/20 | 20/20 (no H3a fit is affected; spans are identical) | \|Δβ\| ≤ 2e-6, rel. ΔSE ≤ 4e-6, \|Δdf\| ≤ 3e-4. This is the 6-significant-digit rounding of the TSV. |
| H3b, 15 random + WGD Z0026 + TP53_altered Z0096 | 17/17 | 17/17 | Same; corrected matches the relook to the printed digits. |
| H1, 5 eligible (ARID1A P, EP300 P, FOCAD N, PTEN N, TP53 N) | 5/5 | 5/5 | Same |
| Controls PC1, PC5, PC9, PC17, PC21 (12 units) | 12/12 | 12/12 | Same |

- **statsmodels cross-check.** statsmodels `OLS(cov_type="HC3")` gives the same SE as my FWL implementation to
  machine precision in every fit.
- **Permutation p-values.**
  - They agree within Monte Carlo error: the largest difference is 0.0165 at p ≈ 0.31. My two seeds differ from each
    other by 0.015 there, too.
  - Where the producer reports a GPD-tail p (e.g. PC17 MITF 6.7e-13), mine is the floor 1/(B+1). That is irrelevant to
    the calls.
- **EP300 P rare route, on A with my fits.** It passes:
  - LOSO same sign in all 30 lineage folds;
  - leave-one-context-model-out same sign over its 12 context models;
  - DEMETER2 on A same sign (β −0.204, 7 context lines).
  This matches `rare_route_pass = True`. (My first DEMETER2 attempt reindexed the scores to all A models; the
  FE-merge then counted models with no DEMETER2 score, and the fit failed. I fixed this to fit on the DEMETER2 rows,
  as the engine does.)
- **Defect diagnosis, reproduced under both designs:**
  - *Failed fit*: H3b WGD Z0026.
    - The all-model reference lineage, Adrenal Gland, has no L model.
    - Testis has one complete-case L model, so it is merged out and becomes a reference group of its own: leverage 1,
      and HC3 refuses the fit.
    - With the L-assembled design, the reference is Ampulla of Vater, which absorbs Testis. The fit gives
      β = −0.0241, p = 0.880 (not confirmed), identical to the relook.
  - *Shifted fit*: H3b TP53_altered Z0096.
    - Vulva/Vagina and Testis are both singletons, so they form a 2-model reference group: one extra design column
      (32 vs 31).
    - The defect design reproduces the primary exactly: β −0.021426, p 0.05064, not confirmed.
    - The corrected design reproduces the relook: β −0.021707, p 0.04904, confirmed. This is one of the 3
      defect-sensitive flips.
- **Full-scale label audit** (`span_audit.tsv`). For all 3,126 hits of the 30 entries, I compared the column spans of
  the two designs on the complete-case L rows, and checked leverage 1. This uses labels and the L row index only.
  - The prediction equals the producer's affected set **exactly**:
    - 510 fitted hits with a different span = the 510 "changed design span";
    - 53 failed fits, all with leverage 1 under the defect design;
    - of those, 50 are fittable under the corrected design (driver defect), and 3 have leverage 1 under both
      (genuine: PC21 ×2, H1S HNF1B);
    - 2,504 hits have an identical span and are not in the relook.
  - No fitted hit has leverage 1 under the defect design.
  - The R6 claim "exactly the predicted set" therefore holds, although the per-hit prediction itself was not committed
    (see N1).

## 3. Rates and labels

- **Rate table.** I rebuilt every row of `comparison_primary_vs_corrected_R7.tsv` from `*/confirmation.tsv` and
  `corrected_refit_nonconfirmatory.tsv`, all 19 direction × subset rows and all 9 numeric fields. There are 0
  mismatches (`rates_independent.tsv`).
  - "Beats all matched sets" uses my own reading of `calibration_D8_tests.tsv.gz`: the discovery's |t| exceeds the
    maximum |t| of its unit's matched random sets (H3a/H3b) or random splits (H2) in that context.
  - Locus status comes from the zavg membership file.
- **Confirmation rule in the committed calls.** Re-applying the rule to the committed L statistics gives 0 mismatches
  in the primary calls. Re-applying it to the relook statistics gives 0 mismatches with `confirmed_relook`.
- **Split counts.** The comparison table omits the "no sets" subset: H3a 62 hits (30 confirmed); H3b 6 (3). Hence
  502 + 581 + 62 = 1,145. `rates.tsv` has these rows.
- **Defect-sensitive list.** Mine has 53 calls, identical to `defect_sensitive_calls_R7.tsv` (no set difference).
  - 50 are NOT EVALUABLE → evaluable:
    - H3b 9 → confirmed, 19 → not confirmed;
    - BRIDGE-L 6 → confirmed, 3 → not confirmed;
    - controls 10 → confirmed, PC18 → not confirmed;
    - H1 FOCAD N → confirmed;
    - H2 Z0026.2 → confirmed.
  - 3 are flips among fitted calls, all with p near 0.05 (H3b TP53_altered Z0096 and Z0284; BRIDGE-L G0347).
- **E2 is never promoted.**
  - H1: 20 rows carry an E2 label or are not headline-eligible. One of them is confirmed (TP53 N), with
    `promotable = False` and status "confirmed (not promotable…)". The only promotable call is ARID1A P.
  - H1S: 293 rows, 9 confirmed, 0 promotable.
  - In `h1_events_primary_vs_corrected_R7.tsv`, `promotable_primary` is True only for ARID1A P. The corrected-only
    FOCAD N confirmation is flagged `defect_sensitive = True`.
- **BRIDGE-L estimate.** My bootstrap uses B = 2,000 and my own RNG, resampling within map × size bin, with pooled bin
  weights. Cluster size = number of `role == member` genes; counting attached genes too does **not** reproduce the
  producer's weights.

| Version | zavg | GLS | Δ (mine) | My 95% CI | Producer |
|---|---|---|---|---|---|
| As executed (failed fits = not confirmed) | 0.5019 | 0.5741 | −0.0722 | [−0.111, −0.031] | −0.0722 [−0.111, −0.033] |
| Primary under R7 (NOT EVALUABLE excluded) | 0.5137 | 0.5785 | −0.0648 | [−0.104, −0.024] | −0.0648 [−0.107, −0.025] |
| Corrected (non-confirmatory) | 0.5094 | 0.5784 | −0.0690 | [−0.108, −0.028] | −0.0690 [−0.108, −0.029] |

All three versions agree: GLS-map discoveries confirm more often, and the CI excludes 0. The bridge is not supported.

## 4. Interpretation aid (POST HOC, DESCRIPTIVE; not registered, not a test)

**Question.** Is about 52% (H3a) and about 46–49% (H3b) confirmation what we would expect, given L's sample sizes and
the winner's curse on the A estimates?

**Method.** Per eligible lockbox-route hit, I computed the power of the one-sided α = 0.05 test in L. The expected
rate is the mean power over hits.
- **Expected effect.** I used four versions of the A estimate:
  - (a) unshrunk;
  - (b) shrunk by a flat 50% factor;
  - (c) the empirical-Bayes posterior mean (Tweedie's formula on the z-scores of **all** evaluable A tests of that
    direction: H3a 8,441, H3b 7,340; log-density by Efron's Lindsey method, degree-7 Poisson GLM);
  - (d) as (c), but integrating over the posterior variance as well.
- **Expected L standard error.** Either the A SE scaled by the L/A context and reference counts (no L outcome used),
  or the observed L SE.
- **Power.** Noncentral t with the L Satterthwaite df; the normal form for (d).
- Output: `power_aid.tsv`, `tables_review.json` → `power`.

| Direction (n) | Observed | (a) naive | (b) 50% | (c) EB mean | (d) EB mean + var | sd of the expected rate |
|---|---|---|---|---|---|---|
| H3a (1,145) | 0.525 | 0.69 | 0.32–0.33 | 0.51–0.52 | 0.47 | 0.013 |
| H3b (150) | 0.492 primary evaluable; 0.460 corrected; 0.400 all-in | 0.68–0.71 | 0.32–0.33 | 0.43–0.44 | 0.43–0.44 | 0.035 |

Ranges span the two SE choices; the observed and A-scaled SEs agree (median ratio 0.99 for H3a, 1.04 for H3b).

**Reading.**
- The observed rates are what the winner's curse plus L power predict. The empirical-Bayes shrinkage factor (median
  0.74 for H3a, 0.59 for H3b) is the main driver.
- H3a is at or slightly above the EB expectation. H3b corrected is within 1 SD of it.
- So roughly half confirming is **not** evidence of non-replication beyond selection. Equally, it is not evidence that
  hits replicate better than their shrunken A effects predict.
- The flat 50% factor is too conservative here, because many H3a hits have very large A t.
- **The "beats all matched sets" split is largely explained by power:**

| Subset | Observed | Expected (d) |
|---|---|---|
| H3a, beats all | 0.572 | 0.52–0.53 |
| H3a, not beats all | 0.489 | 0.41–0.43 |
| H3b, beats all | 0.583 primary evaluable; 0.511 corrected | 0.50–0.51 |
| H3b, not beats all | 0.341 primary; 0.360 corrected | 0.29–0.30 |

  Set-beating hits have larger A effects, so they confirm more. The gap is not specific evidence for complex-level
  signal (MUST-4).
- **Caveats.**
  - The tests are correlated (shared models and units), and the EB prior is pooled across contexts.
  - The t→z conversion and the normal approximation are approximations.
  - The permutation p used for < 20 context models is approximated by the t test.
  - This is an aid for reading the rates, not a result.

## 5. BENCH

**Point estimates from the committed per-pair outputs** (`bench_pairs.tsv`, 1,295 pairs, no NaN):

| Estimate | Recomputed | Reported |
|---|---|---|
| Mean Δ | −0.06597 | −0.06597 |
| Low tercile | −0.338 | −0.338 |
| Mid tercile | −0.100 | −0.100 |
| High tercile | +0.240 | +0.240 |
| Spearman(Δ, breadth) | 0.2014 | 0.2014 |

The tercile labels are the producer's, with 432/431/432 pairs (N6).

**Recomputed from L.** I refitted the per-pair Δ from L (resume path) with my own HC3 fits: S_C, and g\* clipped at
±3, under the same covariates. Every one of the 1,295 Δ matches the committed value (max |diff| 5e-6), and the means
are identical.

**Bootstrap** (`bench_recompute.json`; B = 2,000, the registered seed 5100, the driver's draw order):

| Design used inside each bootstrap replicate | Mean Δ 95% CI | High-tercile 95% CI | Low tercile | Mid tercile |
|---|---|---|---|---|
| Executed (covariates assembled over L, subset to the replicate) | [−0.1549, 0.0313] | [−0.0330, 0.3625] | [−0.515, −0.201] | [−0.226, 0.155] |
| Assembled within the replicate (engine convention) | [−0.1587, 0.0355] | [−0.1224, 0.3623] | [−0.530, −0.207] | [−0.225, 0.180] |

- The executed row reproduces the producer's CIs to 1e-15.
- **Finding: 16.5% of (pair × replicate) Δ values are undefined, and the code silently drops them.**
  - The undefined values are fits refused for leverage 1 (13,010 context × replicate cells) or for fewer than 2
    context models (179).
  - 1,032 of the 2,000 replicates lose at least one pair.
  - `bench_delta` averages each replicate with `nanmean`, so each replicate's mean is over a different pair set.
  - The same failures occur under the replicate-assembled design (16.0%). Most are therefore inherent to resampling
    small contexts, not to the R6 defect. Even so, the two designs differ in 16.8% of cells, and the high-tercile
    lower bound moves from −0.033 to −0.122.
- **Interpretation check against the registered criterion.** The criterion is the top-tercile CI: success iff it
  excludes 0. Otherwise "the paper states that complex scores did not out-replicate the best gene".
  - The top-tercile CI includes 0 under both designs, so the registered statement follows.
  - D12 justifies "the complex does not beat its best gene" with the **overall** mean CI. That is the wrong basis,
    although the conclusion is the same.
  - The result is not evidence of equivalence. In the low-breadth tercile the best gene clearly out-replicates the
    complex (CI excludes 0). In the high tercile the point estimate favours the complex (+0.24; 60% of pairs have
    Δ > 0).
  - Spearman's p = 2.6e-13 treats 1,295 pairs as independent. The pairs share 353 complexes and all L models: the iid
    SE of the mean is 0.038, while the bootstrap half-width is about 0.093. Report rho without that p, or with a
    bootstrap CI.

## MUST

**MUST-1. H3b rate reporting.** The H3b NOT EVALUABLE fits are informative, not missing at random.
- 26 of the 28 are in WGD; in the corrected relook 19 are not confirmed and 9 are confirmed.
- So the primary "evaluable" rate, 60/122 = 0.49, is biased upward relative to the corrected 69/150 = 0.46 and the
  all-in 60/150 = 0.40.
- The H3b "beats all matched sets" contrast is defect-sensitive in 22 of 94 vs 8 of 50 hits: primary evaluable 0.58
  vs 0.34; corrected 0.51 vs 0.36.
- Wherever an H3b rate appears (tree, manuscript, figures), show the three numbers together: evaluable, all-in and
  corrected (non-confirmatory). Do not headline the evaluable rate alone.
- State that the H3b beats-all contrast cannot support a claim (R7: it rests on defect-sensitive calls).
- State that the WGD context is effectively NOT EVALUABLE in the primary.

**MUST-2. BRIDGE-L headline per R7.**
- D12 (`main`) and the `5bc9444` message cite the as-executed estimate, which counts the 37 defect failures as "not
  confirmed": −0.072 [−0.111, −0.033].
- Under R6/R7 the primary excludes NOT EVALUABLE: −0.065 [−0.107, −0.025]. Corrected: −0.069 [−0.108, −0.029].
- Update the tree and manuscript to the R7 primary. Show the as-executed value as a disclosure. The conclusion
  ("negative for the bridge") is unchanged in all three versions.

**MUST-3. BENCH reporting** (§5).
- (a) Justify "complex scores did not out-replicate the best gene" by the registered criterion: the top-tercile CI
  [−0.033, 0.362] includes 0. Do not cite the overall mean.
- (b) Report the low-tercile result (the best gene wins) and the high-tercile point estimate. Do not phrase the result
  as "no better" or equivalent.
- (c) Disclose that 16.5% of bootstrap pair-replicate fits fail and are dropped by `nanmean`, and that the CIs depend
  on the in-replicate covariate assembly: the high-tercile lower bound is −0.033 as executed and −0.122 when
  assembled within the replicate. Do not describe the top tercile as "nearly significant".
- (d) Drop the iid Spearman p-value.

**MUST-4. "Beats matched sets" interpretation.**
- D12 reports "beats matched sets 0.57 vs 0.49" for H3a without context.
- Any text that uses this split must say that the gap is about the size expected from A effect size and L power alone
  (§4).
- It must not be presented as evidence that complex-specific hits replicate better.

## Notes (SHOULD / disclosure)

- **N1. Intent erratum.** `CORRECTED_REFIT_INTENT.md` lists the H1 changed-span set as "17, including ARID1A P and
  TP53 N".
  - ARID1A P (and ARID1A N) have an **unchanged** span: they are not in the relook, and my audit agrees.
  - The 17 are 16 N axes plus EP300 P. The count is right; the example is wrong. Add an erratum.
  - The per-hit prediction was not committed before the check. Only the counts were. My label audit now confirms the
    set independently.
- **N2. Defect-sensitive count in the tree.** D13 says "4 known so far, all at p≈0.05". Those were the 4 status
  differences of the R6 check, including one non-eligible H1S axis. Under the R7 definition the list has 53 calls:
  the 50 refitted fits plus 3 flips. Update D13 to avoid two different numbers.
- **N3. Loader count is inferred.** It rests on the invariants in §1. Commit the run log (`lockbox_run.log`) if it
  still exists in the scratch area.
- **N4. Controls.**
  - Primary: 29/31 evaluable; 29/44 all-in.
  - Corrected: 39/42. 10 defect-failed controls confirm in the relook: PC9 ×3, PC10 ×3, PC12 ×2, PC13 ×2. PC18 KIF18A
    does not.
  - They are defect-sensitive, so the L recovery of PC9, PC10, PC12 and PC13 cannot be claimed from the primary read.
- **N5. H1.** Primary: ARID1A P confirmed and promotable; EP300 P rare route passed; TP53 N confirmed but not
  promotable (E2); PTEN N not confirmed; FOCAD N NOT EVALUABLE (driver defect). FOCAD's corrected confirmation (p 5e-4)
  is defect-sensitive and must not appear as a finding.
- **N6. BENCH terciles** use `pd.qcut` on breadth ranks (432/431/432). A floor(n/3) rule would move one pair. This is
  harmless and follows the registry's "terciles".
- **N7.** The comparison table omits the "no sets" subset (see §3). Add it, or footnote the denominators.

## Not verified

- I did not observe the loader calls directly (N3).
- The H1S, H2 and BRIDGE-L GLS fits are covered only by the label-level span audit and the rate and list
  reconstruction, not by refits.
- The A-side discovery tables (FDR, permutation) were taken as committed. They were reviewed earlier (`fzx_review.md`,
  `discovery_review.md`).
- The score constructions (S_C, H1 components, H2 residualisation) were reused, not re-implemented.

## Reproduce

```
export PART1_DIR="/Users/metinyazar/Desktop/autoresearch/coessentiality_protein complexes"; export PYTHONHASHSEED=0
O=results/part2/verify/lockbox_results_review
.venv/bin/python scripts/verify/lockbox_results_review.py --out $O          # ~1 min; resume path only
.venv/bin/python scripts/verify/lockbox_results_review_tables.py --out $O   # ~10 s; committed files only
.venv/bin/python scripts/verify/lockbox_results_review_bench.py --out $O --B 2000   # resume path only
```
