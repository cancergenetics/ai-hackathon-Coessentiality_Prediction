# Retrospective at milestone M2 (design frozen; test look pending)

Scope: process only. Sources: PROMPT.md, RESEARCH_TREE.md, LESSONS.md, git log (14 commits, 2026-10-28 14:32-17:38 UTC), docs/reviews/*. No huMAP file opened. No metric computed.

## What should we stop doing? (5)

1. Stop hardening the one-shot test harness through three review rounds while design-level flaws go unseen. Protocol/harness took ~28k lines in 7a842f2 and rewrites in 958cfda and 6a40d73. Round 3 still found three blocking DESIGN circularities (lockbox stratified on genes that need complexes built on A; full-L read before confirmations; ill-defined bridge CI): round3_final.md B1-B3.
2. Stop letting agents write to the working tree while a gate needs it clean. The freeze is blocked by untracked src/scripts/tests files from the construction agent (RESEARCH_TREE status; round3_final "untracked Part 2 files"; round2_ai_failure S-B5, pip-freeze hash vs installing Infomap).
3. Stop recording claims in the tree/DESIGN that were not checked against the file. The tree said D1 was "disclosed in TEST_PROTOCOL", but it was not (round1_pi M1, round1_ai_failure MUST-3). DESIGN §9 cited a non-existent round1_response.md (round2_pi, "Gap"). Protocol status said "pending" after the audit passed.
4. Stop re-committing bulky generated outputs. A 5,905-line B_ift_protocol.json diff and 6,264-line C_train_as_test.json diff went in 958cfda, after 10,873 lines in 7a842f2. This is review noise and churn. Commit compact summaries and hashes, and regenerate the rest.
5. Stop hand-editing version labels and "options" into frozen artifacts. Protocol text said "FROZEN v1.1" after v1.2 (round3_final §1; script 22 line 7); "draft" in §10 broke the documented freeze (round2_ai_failure M-B1); the draft passed the gate (round1_ai_failure MUST-1). Generate labels from one source, and test the real freeze edit end to end.

## What worked well? (5)

1. Role-split external reviewers with a fixed brief (BRIEF_round1/2.md): label-free only, MUST/SHOULD/COULD, problem -> change -> why, reply <=200 words. Each round caught real blockers (round1_ai_failure MUST-1/2/6/7; round2_ai_failure M-B1..B3).
2. Reviewers ran label-free checks, not just read: pytest, `require_frozen()`, a scratch git repo to reproduce erase-and-rerun, and PI feasibility counts (round2_ai_failure; round1_pi header; round2_pi table of LoF genes). The counts shrank H1 to an honest panel and drove the scope cut (D4).
3. Explicit response tables mapping every item to made/rejected + location, with reasons for partial acceptance such as escrow (round1_protocol_response.md; DESIGN §9; round1_response.md).
4. Harness validated before any test contact: dry run reproduces Part 1 IFT look 1 to 1.5e-8 and toy tests (34 pass) (7a842f2, 958cfda, round3_final).
5. Early, specific disclosure and audit: Part 1 leakage audit before reuse (7d1fbdd, docs/part1_audit.md); D0/D1 recorded as soon as they happened; protocol §0 lists T0-T3 contact and states the guarantee is procedural (round1_protocol_response.md).

## Proposed pruned LESSONS.md (12 rules)

```
# LESSONS.md - process rules (rule + evidence)
1. Non-evaluating agents and reviewers get no test-label access; every brief says "compute no metric, score or count on test pairs", including top-k/score-threshold tables (they are metrics in disguise). (Part 1 D1; project D1; BRIEF_round1/2.md)
2. Fix PYTHONHASHSEED and all seeds for graph/community code; verify determinism by regenerating in a subprocess. (Part 1 Louvain tie-break bug)
3. Run a power simulation before freezing any statistical gate or CI scheme (and check the CI is well-defined: B, resampling unit, size matching). (Part 1 D2; round3_final B1)
4. Verifiers get only protocol + raw inputs, re-implement from raw data, sweep neighbouring configs, and report before seeing the harness numbers. (Part 1 a2f8622/460fd36; round1_ai_failure SHOULD-4)
5. Pre-register the exact look command; commit the "started" record before reading labels; one look; a relook is non-confirmatory. (Part 1 IFT look 1; round2_ai_failure M-B3)
6. Update the tree in the same commit as each result, and grep any "disclosed/see file X/done" claim against the file before writing it. (Part 1 stale tree; round1_pi M1; round2_pi gap)
7. torch and LightGBM in one pytest process can deadlock on OpenMP; use separate processes. (Part 1 final suite)
8. Request all external data in one batch right after the literature review's dataset section. (day-1 two request rounds; d6bf744)
9. Before any gate that needs a clean tree or pinned env, commit or isolate (worktree/own venv) concurrent agents' work; never let two writers share a gated tree. (round3_final untracked files; round2_ai_failure S-B5)
10. Draw the data-flow order (who is built on which split, who reads the lockbox when) before hardening code; review structure before harness. (round3_final B2/B3; round1_ai_failure MUST-7)
11. Reviews: independent roles, MUST/SHOULD/COULD, require label-free runs (tests, gate calls, scratch-repo simulations) and apply the real freeze edit as a test on the final artifact. (round2_ai_failure M-B1; round1_ai_failure)
12. Commit compact summaries + hashes, not repeated multi-thousand-line generated JSON; one source for version labels. (958cfda 12k-line diff; round3_final version label drift)
```
