# Independent review: revision analyses R-1 to R-4 (post hoc, descriptive)

Reviewer: independent agent (Claude Sonnet 5.5). Worktree branch `worktree-agent-a980a8729585608e6`, fast-forwarded to main 11ab5bf
(code and results identical to 740f439). Scope: `results/part2/revision/` (INTENT.md 4157f1e, results 740f439), scripts 58-61,
`src/part2/revision.py`, `tests/test_part2_revision.py`. Referee items M1-M3 in `manuscript/review/peer_review_round1.md`.

## Verdict: PASS-WITH-NOTES

No MUST items. Every number I could recompute matches the README (to the printed precision). The reproduction of scripts 59-61 is
byte-identical except for wall-time and peak-RSS fields. Label discipline (POST HOC / descriptive) and the lockbox/test-set firewall hold.
The notes are about wording and unmatched referee requests, not about numbers.

**Independence caveat (LESSONS rule 4).** I read the README before computing, so my independence is in code and data path, not in blindness.
My numbers were computed from the committed tables with my own code (scratch only, not committed:
`/private/tmp/claude-502/-Users-metinyazar-Desktop-autoresearch-predicting-protein-complexes/5a84c62a-d8a2-48ad-b26c-5dd43a72db34/scratchpad/rv_review/`:
`r1a.py`, `r1b.py`, `coh.py`, `r2.py`, `r3.py`, `r4.py`, `bb.py`; outputs `*.out`). The tables below list my values first, then the README value.
What is and is not independent:

| Analysis | Independent of the author's code | Shared with the author |
|---|---|---|
| R-1(b) | all statistics recomputed from raw committed t tables | the 10 random atlases and their permutation p (not rerun) |
| R-1(a) | rates, unit counts, coherence (own mean of upper-triangle zavg_BEF), BB (own Simes + BH + BB implementation) | the nonmember test table itself (script 58 not rerun) |
| R-2 | whole pipeline (own natural cubic spline basis, OLS, labels, counts) | none |
| R-3 | HC3 beta/SE/t via statsmodels; own p; own aggregation | `cscore.unit_scores`, `assoc._design`, `HC3Design.df`, `load_A_effects` (the project engine, as instructed) |
| R-4 | observed standardised rates, size bins from the membership tables; order-of-magnitude power | the Tweedie EB step (reused review function) was not re-implemented |

## MUST

None.

## SHOULD

1. **R-2 README sentence is wrong at the margin.** "Every other context stays <= 0.08" (spline R-squared). CNS/Brain is 0.0829
   (`r2_trend_by_context.tsv`; Head and Neck 0.0714, Myeloid 0.0723). Say "<= 0.083".
2. **R-4 "about a third" and the 2.2 SD need a caveat in the Reading.**
   - Expected difference -0.023 over observed -0.065 is 35% for the headline (EB mean + variance). It is 26% under naive power and 41% under
     EB mean. Give the range, not one figure.
   - The headline EB model under-predicts the confirmation level of both maps: zavg 0.459 expected vs 0.514 observed, GLS 0.482 vs 0.578
     (naive over-predicts: 0.686 / 0.703). The residual -0.042 is a difference of two miscalibrated levels. GLS is under-predicted by about twice as
     much as zavg (-0.096 vs -0.055). That is the finding, and the README only says "out-confirm" per size bin. State the level miscalibration.
   - The SD 0.019 is the Poisson-binomial SD for independent hits. Hits share units and contexts. "2.2 SD" is not an inferential statement; the README
     already says descriptive, but should say the SD ignores clustering and EB-prior uncertainty.
3. **R-2: referee requests not addressed are not listed.** M2(c) (exclude scale-concordant leads) is not done: the rule keeps them by construction
   (|z| >= 2 excludes the artefact flag), disclosed as 26 of 37 and in the Reading. M2(d) (per-context inflation factor from the 100 never-expressed sets,
   rescaled statistics) is replaced by a reliability column. The README should say in one line that (d) was not carried out and why.
4. **R-2 lead set is fragile at the threshold.** Four of the 37 corrected leads lie within 0.05 of |z_nl| = 2: Bone Z0001 (-2.00), Breast Z0253 (new lead),
   Lung Z0087, Kidney Z0362; three further non-leads sit just inside. The "37 = 37 with 2 lost, 2 new" statement is a count identity, not stability. Add the sentence.
5. **R-3 "same sign but weak" list is incomplete.** README lists Lymphoid Z0190 (0.14), melanoma Z0190 (0.10), Pleura Z0067 (0.13). Lung Z0150 is a lead, same sign,
   p_dir 0.142 (n_ctx 77) and fits "0.10-0.14". Add it (or restate the criterion).
6. **R-1(b) Reading: tension to state.** Real discoveries beat their own matched sets more often (0.464 vs 0.329; 0.669 vs 0.483) while mean P(|t| real > |t| random)
   = 0.482 / 0.476 (< 0.5). Both are true: the excess is within-group (real unit vs its own 10 sets), not a marginal effect-size advantage; marginally
   real discoveries are, if anything, slightly weaker. One sentence prevents "modest excess" being read as "real complexes are stronger".
7. **Firewall coverage (question 3 of the brief): scripts 58-61 are outside the rule 1-3 glob; it does not change the result but should be closed.**
   `tests/test_part2_firewall.py` and the OUTCOME/LOCKBOX rules of `tests/test_part2_lockbox_firewall.py` glob `scripts/[34]*_*.py` plus `src/part2/*.py`.
   `src/part2/revision.py` is covered; scripts 58-61 are not (neither are 51/55-57 from earlier milestones). I applied the same patterns to 58-61 by hand:
   forbidden-input, outcome-access (`load_gene_effect`, `DEMETER2_CSV`, `.subset(`, `universe_models(`), `LOCKBOX_MODELS`, `lockbox_reads`, `resume_lockbox_read`,
   `read_lockbox`: zero hits. The persisted-read rule uses `rglob` and does cover them (passes). Recommend widening the glob to `[3-6]*_*.py` (check
   whether 51/55-57 need an allowlist) so the next script is covered by construction (LESSONS 12).

## COULD

1. R-1(a): the README compares N-nbr (90 units) with N-ess (355 units). The like-for-like comparison is in the committed table and strengthens the
   third bullet: on the same 90 units N-ess is 0.164 (H3a) / 0.084 (H3b) against N-nbr 0.251 / 0.108 and real 0.230 / 0.113. Add it.
2. R-1(a) "Reading" first bullet ("as context-associated as real complexes") rests on 90 and 87 units with median target coherence 8.9 z (others 12.3 z). The
   caveat is present; keep it adjacent to the claim in the manuscript text.
3. Wall times: README gives R-2 1.5 s; the script logs 0.1 s (compute) and my cold run was 4.6 s wall (imports). State what is timed.
4. The TSV/gz tables carry no label of their own; only the JSON summaries and README do. A first-line `# POST HOC, DESCRIPTIVE` is not possible in these TSVs without
   breaking readers; the folder README is the label. Acceptable, noted.
5. `scripts/61` imports `lockbox_results_review_tables.py` (whose other functions read L estimates) for `tweedie` alone. Import is side-effect free (checked: only
   definitions run at import), so no L access, but moving `tweedie` into `src/part2/` would remove the dependency on a script that can read L.
6. Script 58 (30 min) was not rerun. Instead I verified its outputs internally (below).

## Recomputed numbers vs README

| Item | My value | README | Match |
|---|---|---|---|
| R-1b H3a real beats all 10 sets | 502 / 1,083 = 0.4635 | 0.464 | yes |
| R-1b H3b real beats all | 103 / 154 = 0.6688 | 0.669 | yes |
| R-1b H3a random discoveries beating 9 siblings + real unit | 3,649 / 11,087 = 0.3291 (stored flag also 3,649) | 0.329 | yes |
| R-1b H3b random symmetric flag | 667 / 1,381 = 0.4830 | 0.483 | yes |
| R-1b unconditional symmetric (all random tests) | 0.0867 both | 0.087 | yes |
| R-1b mean P(abs t real > abs t random), same context, ties 1/2 | H3a 0.4824; H3b 0.4762 | 0.482 / 0.476 | yes |
| R-1b same with abs beta | 0.4935 / 0.4943 | 0.493 / 0.494 | yes |
| R-1b median abs t real vs random | H3a 3.74 / 3.77; H3b 3.84 / 3.91 | 3.74 / 3.77; 3.84 / 3.91 | yes |
| R-1b median abs beta real vs random | 0.068 / 0.069; 0.046 / 0.047 | same | yes |
| R-1b random discoveries per atlas | H3a mean 1,108.7 (1,016-1,277); H3b 138.1 (111-165) | 1,109; 138 | yes |
| R-1a H3a rates all units: real / member sets | 0.2808 (8,441) / 0.2814 (80,500) | 0.281 / 0.281 | yes |
| R-1a H3b rates all units: real / member sets | 0.1199 (7,340) / 0.1160 (70,200) | 0.120 / 0.116 | yes |
| R-1a H3a on 90 N-nbr units: real / member / N-nbr | 0.2304 / 0.2476 / 0.2509 | 0.230 / 0.248 / 0.251 | yes |
| R-1a H3a on 355 N-ess units: real / member / N-ess | 0.2814 / 0.2824 / 0.1657 | 0.281 / 0.282 / 0.166 | yes |
| R-1a H3b on 87 N-nbr units: real / member / N-nbr | 0.1132 / 0.1116 / 0.1076 | 0.113 / 0.112 / 0.108 | yes |
| R-1a H3b on 354 N-ess units: real / member / N-ess | 0.1192 / 0.1160 / 0.0878 | 0.119 / 0.116 / 0.088 | yes |
| R-1a matched units | 90 of 367 (H3a), 87 (H3b); 0 for N-unif | 90 / 87 | yes |
| R-1a median target coherence: matched vs unmatched | 8.93 vs 12.34 z (H3a); 8.83 vs 12.35 (H3b) | 8.9 vs 12.3 | yes |
| R-1a pool N | 2,380 of 6,022 S genes | 2,380 | yes |
| R-1a BB rejections (own BB implementation, q = 0.05) | N-nbr 1,666 / 44; N-ess 2,459 / 337 (H3a / H3b); member sets 11,861 / 1,278 | 1,666 / 44; 2,459 / 337 | yes |
| R-2 H3a strongest findings | 207 | 207 | yes |
| R-2 abs z_nl >= 2 (and linear) | 79 (79) | 79, "the same" | yes |
| R-2 scale-artefact-possible, uniform rule | 95 (D16 rule 56) | 95 (56) | yes |
| R-2 corrected leads | 37; D16 leads 37 | 37; 37 | yes |
| R-2 kept / lost / new | 35 / 2 (Bone Z0130, Myeloid Z0278) / 2 (CNS/Brain Z0146, Breast Z0253) | same | yes |
| R-2 scale-concordant leads | 26 of 37 | 26 | yes |
| R-2 Z0190 (Lymphoid, +5.92), Z0196, Z0363, Myeloid Z0001 stay leads | all four leads (z 5.92, 2.28, 2.43, 2.13) | stay leads | yes |
| R-2 max abs difference in residual z vs committed | 5e-5; label mismatches 0 / 0 | - | yes |
| R-2 per-context table (17 contexts) | all strongest / z>=2 / artefact / leads / never-expressed rates equal | table | yes |
| R-2 never-expressed p < 0.05 (Lymphoid, melanoma, Skin, Myeloid) | 0.27, 0.33, 0.28, 0.01 | same | yes |
| R-2 R-squared Lymphoid linear to spline | 0.068 to 0.120 (committed table; my residual SD 0.0722 equals it) | same | yes |
| R-2 "every other context <= 0.08" | max 0.0829 (CNS/Brain) | <= 0.08 | NO (SHOULD 1) |
| R-3 leads evaluable | 35 of 37 (PNS Z0190, Kidney Z0362 not) | 35 / 37 | yes |
| R-3 same sign / p_dir < 0.05 / Spearman / median n_ctx | 29 / 20 / 0.86 / 17 | 29 / 20 / 0.86 / 17 | yes |
| R-3 leads with >= 12 context models | 28; 24 same sign; 18 p < 0.05; rho 0.87 | 28 / 24 / 18 / 0.87 | yes |
| R-3 discordant leads | Soft Tissue Z0080 (n 6), Kidney Z0323, Pancreas Z0283, Soft Tissue Z0336, Lung Z0001, Bone Z0001 | same six | yes |
| R-3 per-hit agreement with committed table (35 hits) | max abs diff beta 3e-6, p 5e-6 | - | yes |
| R-3 all 207 (from committed table) | 194 evaluable; 157 same sign; 88 p < 0.05; rho 0.65; median n 20 | same | yes |
| R-4 lockbox-testable hits / evaluable | 2,630 / zavg 1,267, GLS 1,326 | 1,267 / 1,326 | yes |
| R-4 observed standardised rates, difference | zavg 0.5137, GLS 0.5785, -0.0648 | 0.514 / 0.578 / -0.065 | yes |
| R-4 unweighted observed | 0.5217 / 0.5762 | table (0.5217 / 0.5762) | yes |
| R-4 order of magnitude of expected rate (my naive normal power, predicted L SE) | zavg 0.715, GLS 0.722 (README naive t-based 0.686 / 0.703; same size and order) | 0.686 / 0.703 | consistent |
| R-4 EB-shrunk point-effect power (per-hit shrink factors from the committed r4_per_hit table) | zavg 0.528, GLS 0.543 (README EB mean 0.500 / 0.526) | 0.500 / 0.526 | consistent |
| R-4 headline expected / difference / z | not re-derived (Tweedie step reused); committed table: 0.459 / 0.482 / -0.023 / -2.2; shares = 35% | "about a third" | arithmetic yes |
| R-4 median abs beta_A GLS vs zavg | 0.081 vs 0.066 | 0.081 vs 0.066 | yes |
| R-4 per-bin pattern (GLS observed 0.57-0.59 vs expected 0.47-0.50; zavg 0.48-0.57 vs 0.42-0.49) | from `r4_by_size_bin.tsv`: GLS 0.570-0.591 vs 0.470-0.500; zavg 0.481-0.568 vs 0.416-0.488 | same | yes |

## Reproduction and process checks

- **Scripts 59, 60, 61 rerun** with `PART1_DIR` set, `PYTHONHASHSEED=0`, `.venv/bin/python` of the main checkout, `data_cache` and `datasets` symlinked read-only from the main checkout
  (symlinks not committed). `git diff` afterwards: only `runtime_s` and `peak_rss_bytes` lines in `r2_summary.json`, `r3_summary.json`, `r4_summary.json` changed. All TSV, gz and other
  JSON outputs are byte-identical (git reported no change). I restored the three JSON files (`git checkout`) so the revision folder is unmodified.
- **Script 58 not rerun (30 min).** Internal checks: `r1_summary.json` states discovery selection reproduces the committed flags (H3a 1,145 with 1,666 screened; H3b 161 with 190); member-set t reproduce
  to 6 significant digits (max relative difference 5.0e-6; the INTENT.md 1e-6 tolerance was replaced by the stored precision, disclosed in the README); all 12,468 random discoveries have
  `discovery` true and `q_adj_perm` < 0.1 (max 0.04999); 10 replicates per direction; per-replicate counts (H3a 1,016-1,277, H3b 111-165) equal mine. My own BB implementation reproduces the BB
  counts of the nonmember tests and the committed 11,861 member-set rejections (so the "entry q" is 0.05 and the 1,666 for N-nbr is a coincidence with the 1,666 screened real tests, not a copy error).
- **INTENT before computation:** INTENT.md commit 4157f1e (13:43) precedes the results commit 740f439; the only recorded deviation (t precision) is disclosed.
- **No registry / ledger change:** `git show --stat` of the two revision commits lists only `results/part2/revision/`, scripts 58-61, `src/part2/revision.py`, `tests/test_part2_revision.py`;
  `trials.jsonl` and `lockbox_ledger.jsonl` untouched. No `run_association` call in 58-61 or `revision.py`.
- **Firewall.** Grep of scripts 58-61 and `revision.py` for lockbox/outcome/test-set patterns: only comments and the docstring "no lockbox data", the column-restricted read of committed
  `confirmation.tsv` (`n_L`, `n_context_L`, `lockbox_testable`, `eligible`; no L value) and `primary_relabelled_R7.tsv.gz` (`primary_call` only), `strongest_findings.tsv` with an explicit non-L `usecols`
  (script 59), and DEMETER2 through `assoc.load_A_effects(source="demeter2")` (423 rows, subset of A, asserted). I did not read lockbox_reads, any L estimate, or any hu.MAP file.
- **Labels.** `r1_summary.json`, `r2_summary.json`, `r3_summary.json`, `r4_summary.json` carry a `label` "POST HOC, DESCRIPTIVE"; the folder README states it first; INTENT.md says it in the title. Script docstrings say so.
- **Overstatement check of the Readings.**
  - R-1(a): the "parity" reading is carried by 90/87 units of lower coherence and says so; the claim "coherence, not membership, carries most of the excess over essentiality alone" is supported
    like-for-like (COULD 1). OK.
  - R-1(b): "selection alone lifts 9% to 33% / 48%" correct; "modest excess" fair; see SHOULD 6.
  - R-2: "26 of 37 concordant; the rule keeps them" is honest. "A flexible fit does not absorb them" is true for the four named, at 4-df spline.
  - R-3: 29/35 and 20/35 are descriptive; low power at median n_ctx 17 is noted via the n columns. "Same sign" is not a null-calibrated statement (no comparison to random sets); not claimed.
  - R-4: "not fully explained by power" is the right level, subject to SHOULD 2.

## Commands run

```
git merge --ff-only main                                    # base 740f439 -> 11ab5bf (coordinator instruction)
ln -s <main>/{data_cache,datasets,external} .               # read-only symlinks, untracked, not committed
export PART1_DIR="/Users/metinyazar/Desktop/autoresearch/coessentiality_protein complexes" PYTHONHASHSEED=0
.venv/bin/python scripts/59_rev_scale_effect.py             # 4.6 s wall
.venv/bin/python scripts/60_rev_demeter2.py                 # 14.0 s
.venv/bin/python scripts/61_rev_bridge_power.py             # 4.0 s
git diff --stat; git checkout -- results/part2/revision     # only runtime/RSS lines differed; restored
.venv/bin/python <scratch>/{r1a,r1b,coh,r2,r3,r4,bb}.py     # independent re-implementations (scratch only)
.venv/bin/python -m pytest tests/test_part2_revision.py tests/test_part2_lockbox_firewall.py tests/test_part2_firewall.py tests/test_testeval.py -q
    # 53 passed, 1 skipped (test_testeval.py:504 "frozen marker present"), 17.4 s
```

Runtime of the review: about 45 minutes of wall time including environment set-up; compute is under 2 minutes in total.
