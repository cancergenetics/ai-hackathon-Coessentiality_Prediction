# Round 2 review: AI-agent failure modes (ai_failure)

**What I ran.** `pytest tests/test_testeval.py` (31 passed). `require_frozen()`, which refuses ("FROZEN missing"). In a scratch git repo: `check_protocol_text()` on the protocol with only its header edited, and `prior_look_evidence()`. I opened no test file, ran no look and computed no metric.

**Round-1 status.**
- M1 (a draft passes the gate): resolved.
- M3 (disclosure): resolved.
- M2 (a deletable ledger): mostly resolved. One gap remains (M-B3).

## B. Protocol and harness: not ready to freeze until M-B1 to M-B3 are fixed (each under an hour)

**M-B1. The documented freeze fails.**
- Problem: `check_protocol_text` rejects "draft" anywhere in the file, and §10 line 540 contains "a draft header". With only the header edited as documented, the check refused, citing line 540. At freeze time an agent will then improvise edits or loosen the check.
- Change: reword line 540. Add a test that applies the real header edit and asserts the check passes.

**M-B2. The escrow cannot follow the protocol's order.**
- Problem: §8.4(9) keeps the test files out of the project until the freeze is reported. But `build_marker()` (script 22) hashes them and fails if they are absent. DESIGN §9 also says escrow was "replaced", which contradicts §8.4.
- Change: either record the DATA.md hashes in the marker and check the files only when scoring and at the look, or drop the escrow. Make DESIGN and the protocol consistent.

**M-B3. A look can be erased and re-run as confirmatory.**
- Problem: `--run` writes the `started` record to disk without committing it. In my simulation, `git checkout -- looks.jsonl` plus deleting the results left no evidence. A look committed on a branch that was then deleted is invisible to `git log --all`, although it is still in the reflog. The toy test stubs `prior_look_evidence`, so the git-history refusal is untested.
- Change:
  - `--run` commits the `started` record itself before reading any label, and prints the hash; the PI reports it to the user.
  - The evidence scan includes `git log --reflog`.
  - Add a test against a real git repo.

**S-B4. `PART1_DIR` bypasses the dry-run guard.**
- Problem: the train paths come from an environment variable. If it points at copies of the test files named `Train_*`, then `20 train`, `21 --mode fit` and `21 --mode dryrun --parts C` run the full test-mode evaluation with no gate. `_check_dryrun_inputs` checks against the same redirected files, so §8.4 mitigation 7 is false.
- Change: assert that the train SHA-256 equals the hash in the frozen `fitted_rules.json` `inputs`, or remove the override.

**S-B5. The pip-freeze hash will block a late look.**
- Problem: installing Infomap or MCL into `.venv` after the freeze makes the gate refuse, which invites edits to `protocol.py`.
- Change: run steps 2–4 directly after the freeze, or give Part 2 its own venv. State which in §8.2.

With these fixed, no further protocol change is needed. §8.4 states the procedural limit honestly.

## C. DESIGN v2 guardrails: the right list, but not yet implementable as written

**M-C1. The lockbox is read outside `read_lockbox`.**
- Problem: bridge (v) and the complex-vs-gene benchmark compute L effects for *all* registered families, for two maps, before any per-hypothesis confirmation. DEMETER2 and proteomics also include L lines.
- Change:
  - At FZ-LB, write an A-only effect matrix that is the only one Part 2 code reads.
  - Add a grep test: CRISPR, DEMETER2 and Omics loaders appear only in the FZ-LB writer and in `read_lockbox`.
  - DEMETER2 and proteomics are A-only until confirmation.
  - Bridge (v) and the benchmark become one ledgered L read after every registry entry closes (a new freeze, FZ-L).

**M-C2. Recomputing the scores on A requires new code that can drift silently.**
- Problem: the locked predictor asserts the manifest null moments (TEST_PROTOCOL §1), so it cannot run on A unchanged.
- Change: before FZ-C, commit a test that the A pipeline (including the lean top-50 variant), run on U, reproduces the locked scores and top-50 lists to ≤1e-12.

**S-C3. The registry must be machine-readable.**
- Problem: `run_association` cannot check prose.
- Change:
  - One schema-checked `docs/registry/<ID>.yaml` per entry: sign, family, hashes of the complex set and the context table, covariates, test, success criterion and stop rule.
  - The wrapper refuses an entry that is uncommitted or has mismatched hashes, and commits `trials.jsonl` after each call.
  - A changed input needs an amended entry.
  - Exploration in A also goes through the wrapper, under `EXPL-*` IDs.

**S-C4. `read_lockbox` needs a pre-committed target list.**
- Problem: H1 may yield many hits, so "one read per hypothesis" is ambiguous.
- Change: commit and hash the exact (complex, context, sign) list before the read; the read returns only those rows, and a second read is refused.

**S-C5. The controls still lack a statistic (round-1 M9).**
- Problem: §2.7 gives signs only. The kill rule does not say which controls belong to which direction.
- Change: give each row a statistic, a one-sided α, a minimum n and its direction. The controls script refuses a second run.

**S-C6. Specify the freeze manifests and notarise them.**
- Name each manifest file and add one `require_freeze(name)` check.
- Report the FZ-LB, FZ-C and FZ-X hashes to the user in chat.
- FZ-LB must come before the "construction infrastructure being built in parallel" touches any model rows.

**C-C7. Fix two loose ends.** Define "selective-gene coverage", the criterion for choosing (z0, k). Fix the dangling reference `docs/reviews/round1_response.md`.
