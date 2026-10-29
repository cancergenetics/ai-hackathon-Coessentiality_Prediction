# Lockbox machinery: independent pre-read code review

Date: 2026-10-29. Reviewed commit: `d9a42dd` (branch `worktree-agent-a3ec5565842d4e968`; the review branch
`review-lockbox` starts there and adds only this note and `tests/test_lockbox_review_adversarial.py`).
Scope: `scripts/48_lockbox_confirm.py`, the `src/part2/bench.py` changes, the `confirm.*` / `assoc.read_lockbox`
path the driver uses, amendments A3 and A4 (`scripts/49_amend_A3.py`, `scripts/50_amend_A4.py`, registry entries,
`docs/registry/AMENDMENTS.md`), `results/part2/lockbox_confirm/{BLOCKED,GAPS}.md`, and the synthetic tests.

**No L value was read in this review.**
- I did not read any CRISPRGeneEffect or DEMETER2 row of an L model, and I did not open the lockbox model list.
- I did not run `48_lockbox_confirm.py --dry-check` on real data. It calls `assoc.load_A_effects`, which `np.load`s
  the whole cached gene-effect matrix (L rows included) and then subsets it to A, so L values are held in memory even
  though they are not returned (see N9).
- Real files I read: registry entries, committed A-side discovery parquet schemas and counts (A only), `classes.tsv`,
  `controls_summary.tsv` and `calibration.json`. Everything else ran on synthetic data.
- No ledger or target list exists on any branch (`git log --all` on `results/part2/lockbox_ledger.jsonl` and
  `results/part2/lockbox_targets` is empty).

## Verdict: GO-WITH-NOTES, conditional on MUST-1 to MUST-3

The core read-safety design is sound:
- one gated loader;
- the ledger is committed before any load, and rolled back if the commit fails;
- one read per entry;
- BENCH and BRIDGE-L are refused until all 28 confirmatory entries are closed and committed;
- every failure I could provoke is fail-closed.

The confirmation rules match the registry for every entry type, and the hypotheses come from the committed discovery
tables (verified by SHA) with no re-selection.

Two problems remain in an irreversible one-shot step: integrity of the inputs that decide eligibility and route, and
post-read failure handling. Both are cheap to fix before the read. **Do not start `--run` until MUST-1 to MUST-3
are closed, and have the fix diff reviewed.** None of them needs a redesign.

### MUST (before the read)

**MUST-1. Pin the eligibility and route inputs, and record them in the intent commit.** The driver reads these from
the working tree with no committed or clean check:
- `classes.tsv`, which decides H1/H1S eligibility;
- `lockbox_testable.tsv`, which decides the lockbox-or-rare route. Its FZX-manifest SHA is checked only in the
  BLOCKED.md narrative, not in code;
- `controls_summary.tsv`;
- the D8 calibration tables.

The committed target list records `(context, unit, beta_A, alternative)` but not `eligible`, `route`, `label_final`
or `headline_eligible`. `confirm.confirm_entry` also tests whatever hits it is passed and never cross-checks them
against the committed target list. Two tests demonstrate this:
- `test_uncommitted_eligibility_or_route_inputs_refused[classes|testable]`: an uncommitted edit that makes every
  N-axis call eligible, or re-routes the rare contexts to "lockbox", runs the read without any refusal.
- `test_hits_must_match_committed_target_list`: hits with flipped signs are tested after the intent commit.

Fix:
1. Refuse `--run` unless `git status --porcelain` is empty for tracked files. At minimum, `registry.committed_blob`
   every input above, and compare `lockbox_testable.tsv` with its `FZX_manifest.json` SHA.
2. Write `eligible`, `route`, `label_final` and `headline_eligible` per hit into the target list.
3. Have `confirm_entry` refuse if its `hits` differ from the committed target list's `hits`.

**MUST-2. Post-read failure policy and a full-scale dress rehearsal.** A crash after a read is detected, and it never
triggers a silent re-read (`test_crash_after_read_is_fail_closed_and_never_rereads`). It is, however, permanent: the
entry's L matrix lives only in memory, the ledger line forbids a second read, and the driver has no resume path. The
BENCH bootstrap runs for about 4 h after its read on a laptop, with the L matrix only in memory. Several post-read code
paths have never run at real scale:
- H1S, with 317 events (321 event-axis tests);
- `bench_delta` with `gstar` or `gstar_beta` = None. `choose_gstar` can return None, and the result is a KeyError
  after the read;
- the lazy eigengene A fits, computed after the read in `sL_eig`;
- the GLS `check_frozen`, which is first exercised after BENCH has read.

Before the read:
1. Run the whole driver once on real A data in a throw-away clone (not a worktree). Use a dummy outcome loader that
   returns random values for the L model IDs. IDs and context labels only, never values.
2. Move every A-side computation and validity check (g* not null, eigengene fits, GLS frozen inputs, coverage tables)
   ahead of the corresponding target-list commit.
3. Get a PI decision on crash recovery: either persist each read matrix to `data_cache`, with its hash committed right
   after the read, and add a resume mode that never calls the loader; or explicitly accept the loss.
4. Run under `caffeinate -i`.

**MUST-3. Make the amendment record complete and accurate.**
- (a) This branch has no RESEARCH_TREE entry for the A4 rulings: D9 is only on main, and there is no D10. Record them.
- (b) A4 G1's reason says the random-set comparison "needs L values of non-member genes". That is not true for the
  only implemented matching rule. `calib.matched_random_sets` draws from the members of frozen coarse complexes, and
  GAPS §2 itself says the registered read covers them. The other reason, under-specification, stands. Add an erratum.
- (c) Two readings in the code are not listed among A4's accepted readings:
  - Controls in lineage-like contexts (`melanoma`, `heme`, `heme_vs_breast`; `cscore.LINEAGE_LIKE`) get an
    unrestricted, PatientID-blocked permutation, although the controls' rule says "within-lineage". `melanoma` has
    n_L 19 < 20, so the permutation runs and this reading decides PC17.
  - H1/H1S permutations are PatientID-blocked. The H1 rule says only "within-lineage". This matches the discovery
    runs.

  I agree with both readings, but they must be written down before the read.

## 1. Read safety

| Question | Finding |
|---|---|
| Functions that touch `lockbox_models.txt` | Code that reads L values: only `assoc._model_list(proj, "lockbox_models.txt")` inside `assoc.read_lockbox`, after the ledger commit. That list is then passed to `proj.outcome_loader`. Other references: `paths.LOCKBOX_MODELS` (constant); `scripts/33_make_lockbox.py` (its writer); and two pre-existing verifier scripts, `scripts/verify/discovery_lib.py` and `controls_independent.py`, which read the ID list only, for exclusion. These verifiers are outside the firewall scan glob `scripts/[34]*`. The driver never names the list. `test_part2_lockbox_firewall.py` passes. |
| Call sites of `read_lockbox` | `confirm.confirm_entry` (called once per `run_entry`, for the controls, the directions and the BRIDGE-L GLS hits) and `bench_stage`. No other caller. The rare route uses DEMETER2 **on A** through `load_A_effects`. |
| L before the intent commit? | No. `read_lockbox` requires a committed, clean target list and commits the ledger line before it loads anything. `run_entry` computes the rare routes (A side) before it writes the target list. |
| More than once per entry? | The driver: no. `run_entry`, `gate_late` and the preflight all refuse an existing ledger line, and `write_target_list` refuses an existing list. The machinery: `read_lockbox` allows one read per entry **per source**, so a second, DEMETER2 L read of the same entry is accepted (`test_second_source_read_of_same_entry_refused`, xfail). The driver never makes one, and the registry stopping rules say "one lockbox read (crispr)". See S1. |
| BENCH / BRIDGE-L gated? | Yes, twice. `gate_late` checks working-tree status plus the committed ledger. `read_lockbox` checks committed `status: closed` for all 28 entries (`registry.CONFIRMATORY` includes H1S). `gls_discovery` runs only after `gate_late`. Calling it directly has no gate, but it is A-only (N10). |
| BENCH reads only members? | Yes: the members of the complexes in the primary and secondary pair sets. Under A4 the matched sets are empty. The read still covers the secondary-pair members, roughly every coarse complex, although the secondary is NOT RUN. That is registry-literal, so I accept it. Those values have no registered use (N7). |
| Ledger atomicity, crash detection | The ledger is appended, then committed with `git commit -- <ledger>`. If the commit fails, the line is rolled back and nothing is loaded (`test_ledger_commit_failure_rolls_back_and_reads_nothing`). A process killed between append and commit leaves a dirty ledger, which blocks every later read and that entry in the driver (`test_process_killed_between_ledger_append_and_commit_blocks_every_read`). A crash after the read leaves a committed line, an open entry and no `confirmation.tsv`. This is detectable from git alone and never re-read. It is also unrecoverable (MUST-2). Hazard: `confirm.close_entry` would accept closing such an entry by hand, after which BENCH would proceed with a missing confirmation (S6). |
| L into A | `load_A_effects` refuses non-A rows (tested). The driver trusts `World.E_A`; only `real_world()` checks the A hash, through `44.real_env` (`test_driver_refuses_world_with_lockbox_rows_in_A`, xfail; S2). |

## 2. Faithfulness to the registry

| Entry | Hits (source) | Rule in code vs registry | Finding |
|---|---|---|---|
| CTRL ×23 | Every `controls_summary.tsv` unit with an A estimate. The A sign is checked against the registered one-sided direction. | `lockbox_test`: one-sided in the A direction; HC3 with Satterthwaite df; ≥5/5 L models; p < 0.05; exact Freedman-Lane permutation, 10,000 below 20 context models. Seed 6700+index. | Matches, except the lineage-like reading (MUST-3c). Each control has a single context, so `ctxs[0]` for the hom-del share is safe. PC22's complex unit is not evaluable, so no fitted DEMETER2 parameters are needed. If a gene is missing from DEMETER2, the rare route "fails" rather than "not evaluable" (N5). |
| H3a / H3b | `discovery == True` rows of the data_cache parquet, whose SHA equals the committed summary. The columns are bool (checked: 1,145 / 161). | Same sign; p < 0.05; ≥5/5; permutation within lineage for non-lineage contexts. Route from the frozen list. Rare route for KMT2D and the EWSR1 fusion only. 95% CI for every hit. | Matches. A frozen-testable hit with fewer than 5/5 complete cases is not confirmed and is flagged (reading 6; accepted). |
| H1 / H1S | Every evaluable event axis. Eligible = engine class COLLAPSE or RESIDUAL SENSITISATION on N (κ_N), PARALOG COMPENSATION on P (β_P). | Same members; alignment from A baselines (`h1.score_L_fn`); within-lineage permutation; ≥5/5; neutral calls get the 90% CI. | Matches. E2 is honoured: eligibility uses the engine class, as E2 requires ("still go to confirmation"), and `label_final` / `headline_eligible` are carried through unchanged. Non-promotion is not mechanical, though: the `status` reads "confirmed" for PTEN and TP53, and nothing computes the strongest-findings flag (N4). H1S has 6 N-axis calls NOT EVALUABLE (A4 G3). |
| H2 | Committed stage-wise discoveries. | A-fitted coefficients (`h2.score_L_fn`); target list = the parent's members. | Matches. No H2 discovery is in a rare context (the preflight would otherwise raise for a missing ruling). |
| BENCH | H3a+H3b discoveries in BENCH's 41 lockbox-testable contexts. g* by (signed t, \|β\|) among members after removals. Breadth over all members. | Delta = signed z_L(S_C) − z_L(g*); PatientID-blocked bootstrap, B 2000, seed 5100; Spearman(Delta, breadth); `top_tercile_wins` = lower CI > 0. | Matches the estimand. Terciles use `rank(method="first")`, so ties in breadth are broken by (context, unit) order (N11). The T1/T2 multiverse recomputes terciles on its subset instead of reusing the primary's (N8). |
| BRIDGE-L | The registered GLS discovery: two `run_association` calls, seed 4500, q 0.05, screen q 0.1, 12/12 evaluability, release_hits. | `bridge_rate`: pooled bin weights, bootstrap within map and bin, B 2000, seed 5200. zavg statuses from the closed H3a/H3b tables. Rare hits excluded. | Matches. Omitting the pseudo-context calibration is equivalent to the "H3a/H3b machinery unchanged", because the H3a/H3b calibration response did not fire (`calibration.json`: `response_fires: false` for both). `bridge_rate` renormalises the weights over the bins present in each map (N12). |

A3 seeds: `confirm_entry` uses `seed + i`, with i the position in the (context, unit)-sorted hits, as declared. The
offset ranges overlap other bases. That is harmless: each test stays exact, and at worst p-values are correlated
across entries (N6).

## 3. Amendments A3 and A4

Both are pre-read: their commits precede any ledger or target list, which exist on no branch. Both are machine-checked:
`amendment_changes == []`, the declared paths changed, and the `--check` scripts re-derive the text (tests pass). A3 is
limited to the gap: (a) seeds were genuinely missing; (b) `within_line_rank` truly needs genome-wide L ranks. A4 G2
(multiverse point estimates only) is reasonable for BENCH, where bootstrap CIs would add about 11 h.

Rulings or readings I would make differently. These are for the PI and are not MUST items, apart from the record
fixes in MUST-3:
1. **G1, "number of distinct complexes among the g\*"** is A-only. It needs no L value and was already computed
   (920 / 353 / 353). Declaring it NOT RUN loses a registered descriptive at no saving. I would report all three
   counts.
2. **G3.** The same registry sentence prescribes "alignment from A baselines" for the L route, so "A CRISPR
   parameters" is the registry-consistent choice for the DEMETER2 rare route, and it adds no analyst choice. Even
   under NOT EVALUABLE, the code also skips LOSO and LOMO for these 6 calls, although neither needs fitted
   parameters. Also, "NOT EVALUABLE **in L**" is a misnomer: the rare route is an A-side procedure.
3. **G2 for BRIDGE-L.** The multiverse bootstrap costs under 1 min, so the cost argument does not apply. This is mild.
4. **Literal reading 1 (breadth over all members, removed ones included).** g* and S_C are after removals. Including
   a removed, context-defining gene in breadth, such as the LoF gene itself, is inconsistent. I would compute breadth
   after removals, but the literal text supports the current choice. This is mild.

Accepted as written: readings 2, 3, 5, 6 and 7, and reading 4 given the non-firing calibration.

## 4. The `bench.bridge_rate` fix

Correct. In pandas ≥ 2.1 (the venv has 3.1.0), `h.map` is the `DataFrame.map` method. I re-ran the pre-fix code from
`75bf5e6` and it raises `KeyError: False`. The regression test `test_bridge_rate_uses_map_column` exercises exactly
that path and checks the rates. No other `.map` attribute access remains in the driver, `bench.py` or `validation.py`.

## 5. Synthetic tests

The producer's suite (`tests/test_part2_lockbox_driver.py`) passes: 37 tests with the firewall and registry tests,
88 s. It genuinely covers:
- a planted confirm and a planted A-only fail in H3b, controls and H1;
- the rare route firing for KMT2D and the EWSR1 fusion only;
- BENCH and BRIDGE-L refused out of order, by the driver and by `read_lockbox`;
- re-read refusal and the ledger order;
- E2 labels unchanged;
- the BENCH and BRIDGE-L outputs, A3/A4 re-derivation, and the `bridge_rate` regression.

Gaps:
- The synthetic world has **no H1S entry**, the largest confirmatory set.
- The H1 rare route and the G3 NOT EVALUABLE path are only forced in a unit test.
- The controls' rare route (PC15, PC16) is not exercised.
- Crashes and tampering are not tested.

My additions are in `tests/test_lockbox_review_adversarial.py`: 5 pass and 5 are `xfail(strict=True)` gap
demonstrations, 64 s. With `--runxfail`, each xfail fails for exactly the stated reason.

| Test | Result |
|---|---|
| crash after the read: fail-closed, never re-read, BENCH blocked | pass |
| ledger commit fails (pre-commit hook): rolled back, nothing loaded; then exactly one read | pass |
| killed between ledger append and commit: every read blocked, no loader call | pass |
| target list edited after intent: an uncommitted edit is refused; a committed edit is detectable from the ledger `targets_blob` | pass |
| loader leaks an L row into A: refused | pass |
| hits must match the committed target list | xfail (MUST-1) |
| uncommitted `classes.tsv` / `lockbox_testable.tsv` refused (2 params) | xfail (MUST-1) |
| a World with an L row in `E_A` refused by the driver | xfail (S2) |
| a second, DEMETER2 L read of the same entry refused | xfail (S1) |

## SHOULD and NOTES

- **S1.** `read_lockbox`: one read per entry regardless of source, or refuse any source an entry does not register.
- **S2.** `preflight` / `run_all`: check `set_hash(World.E_A.index) == A hash`.
- **S3.** Move to the preflight everything that can fail late: `check_frozen(inputs_gls)`, `_bench_pairs` (g* not
  null), eigengene A fits, `match_coverage` (overlaps MUST-2).
- **S4.** Add an E2 marker to `status` where `label_final` differs from the engine class, for example
  "confirmed (E2: not promotable)".
- **S5.** The driver neither writes nor requires an intent note (BLOCKED.md mentions one). Have `--run` require a
  committed intent file naming HEAD.
- **S6.** `close_entry` called by the driver should require that entry's committed `confirmation.tsv`. Document that
  after a crash an entry must never be closed by hand without a PI ruling.
- **S7.** Add an H1S entry, and controls in rare contexts, to the synthetic world.
- **N4.** The strongest-findings and E2 non-promotion logic lives downstream. Make sure the report reads
  `headline_eligible`.
- **N5.** If a gene is missing from DEMETER2, the rare route counts as a "fail". Label it "DEMETER2 not available"
  instead.
- **N6.** Seeds: prefer `default_rng([base, i])` in future. There is no action now.
- **N7.** BENCH reads the secondary-pair members with no registered use under A4. Registry-literal; disclose.
- **N8.** The T1/T2 multiverse terciles are not fixed from the primary's breadth.
- **N9.** `load_A_effects` loads the full cached matrix, L rows included, into memory. Using `np.load(...,
  mmap_mode="r")` with row indexing would keep L pages unread. The dry-check mode therefore also holds L values in
  memory.
- **N10.** `gls_discovery` has no gate of its own (A-only).
- **N11.** Breadth tercile ties are broken by (context, unit) order. Disclose.
- **N12.** `bridge_rate` renormalises the bin weights per map when a map has no hits in a bin. Disclose.

## Reproduce

```
git -C <main> worktree add <main>/.claude/worktrees/review-lockbox -b review-lockbox d9a42dd   # + symlink .venv, datasets, data_cache, external
export PART1_DIR=".../coessentiality_protein complexes"; export PYTHONHASHSEED=0
.venv/bin/python -m pytest -q tests/test_part2_lockbox_driver.py tests/test_part2_lockbox_firewall.py tests/test_registry_entries.py  # 37 passed
.venv/bin/python -m pytest -q -rxX tests/test_lockbox_review_adversarial.py         # 5 passed, 5 xfailed
.venv/bin/python -m pytest -q --runxfail tests/test_lockbox_review_adversarial.py   # the 5 gaps fail as described
.venv/bin/python -m pytest -q tests/test_part2_assoc.py tests/test_part2_fzx.py     # 24 passed
```
