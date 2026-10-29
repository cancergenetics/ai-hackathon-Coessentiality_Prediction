# Retrospective at milestone M6 (manuscript submission version; lockbox read and relook complete)

Scope: process only, since M3 (8f63b6e), as run by a PI agent managing coding, writing and review agents. Sources: RESEARCH_TREE.md (D5-D17), LESSONS.md, retro_M2.md, retro_M3.md, `git log --oneline -150` (the window opens mid-lockbox, at 1047a0c CTRL_PC7), results/part2/verify/*.md (MUST items), results/part2/lockbox_confirm/README.md, REHEARSAL.md. No datasets or hu.MAP files opened; no metric computed.

## What should we stop doing? (6)

1. **Stop finding registry-to-code gaps after the step they gate.** This happened four times, although lesson 15 existed from M4.
   - Calibrations were built after the discovery runs (44f30ca; D8 R4, "disclosed deviation").
   - The lockbox driver, BENCH secondary, BRIDGE-L discovery and seeds were missing at the dry check (D9, BLOCKED.md).
   - The L covariate design was built on all models and then subset, which left 53 fits rank-deficient and shifted 510 more (D12/D13; 38e301a, 68f0d3f).
   - The registered LOSO / leave-top-2 robustness was never computed until after the lockbox read (D15; fd03b03 → 01f4c20, "disclosed").
   - Lesson 15 covers only "the engine" and registry entries. The robustness rule lived in DESIGN, so it fell outside the check. The remedy is a mechanical traceability table (new lesson 13), not another reminder.
2. **Stop treating a passing rehearsal and a reading-style code review as evidence the driver is correct.**
   - The dress rehearsal (REHEARSAL.md) used a random dummy loader, ran 2.7 h, reported "no error", and never inspected fit status. It had the same 53 failures (D12).
   - The pre-read review (lockbox_code_review.md, GO-WITH-NOTES) compared the driver with the registry in tables and closed MUST-1 to MUST-3, but did not execute the L path against the A path on the same rows.
   - The defect was found only after the single read, by the R6 reproduction check.
   - Rehearsals and reviews must assert outcomes (fit counts, equality with the A-side path), not completion.
3. **Stop scoping a remediation before measuring the defect's footprint.**
   - R6 (4d48039) authorised a refit of 53 failed fits.
   - The reproduction check (68f0d3f) showed 560 fits were affected (510 changed + 50 failed), so R7 had to extend the intent (b0ae055, 2566ab9). The intent also named the wrong H1 axis (D14).
   - The results reviewer's label-only span audit predicted all of 510/50/3 exactly. Run that audit before writing the intent.
4. **Stop the five-commit-per-entry ceremony for mechanical steps, and stop leaving the one irreversible run unlogged.**
   - 125 of the last 150 commits are "lockbox targets / read / read persisted / Lockbox confirmations / close registry entry" for about 27 entries, most stamped within one minute (06:50).
   - Only the ledger line and the hash-persist commit carry the firewall guarantee. The rest bury the 25 substantive commits.
   - The reviewer still had to write "the loader was called exactly 30 times: inferred, not observed. The run log was not committed" (lockbox_results_review.md §1). Commit the run log, and batch targets, close and report per run.
5. **Stop letting writing agents set scope and selection, and stop checking numbers last.**
   - The first manuscript draft was about 22k words (3efa551, D15). It was condensed to 6.3k in 6d3a11e, so the writing agent's work was largely thrown away.
   - The writer's "261/29" strongest-findings set used two unregistered filters. It was replaced by the registered 254 (D16).
   - `check_manuscript_numbers.py` (579/579) and the reference check (a303165: 9 metadata errors, 2 claim mismatches) ran only after the draft was complete.
   - Give writers a word budget, registered sources only, and run the number check at draft time.
6. **Stop letting descriptive post hoc annotations quietly become selection rules, and stop leaving dead status text in the tree.**
   - D16 says the matched-set, control-recovery, LOSO and scale-concordance columns are "descriptive and none a filter". D17 then defines "lead example" and "possible lineage-scale artefact" from them (37 / 56). Declare each as a reporting rule before the manuscript uses it.
   - The RESEARCH_TREE.md "Status summary" is still headed "M3, 2026-10-28 ~22:30" and lists as "Next" items that are long done. This is the same stale-status failure as M3 stop 5 and lesson 6.

## What worked well? (6)

1. **A single lockbox read enforced by code, not by promises.**
   - Intent ed065f4 precedes every target commit. Each ledger line is committed before its read, and each read is hashed at once (e.g. 747d37b → 8c3962a). Resume mode never calls the loader, and `tests/test_part2_lockbox_firewall.py` covers every module.
   - The results review found 30 ledger lines, no re-read, `git diff ed065f4 5bc9444 -- src scripts` empty, and the relook ordering enforced by code (lockbox_results_review.md §1).
2. **Disciplined handling of the discovered defect.**
   - The R6 reproduction check was written to stop on failure and did stop (68f0d3f, "STOPPED before the refit"). The failure was predicted from labels beforehand.
   - R7 kept the as-executed single read as primary, added a NON-CONFIRMATORY corrected relook over 560 fits, and labelled 53 calls "defect-sensitive" so they cannot carry a claim (D13, README R7).
   - The independent review reproduced 54/54 calls under both designs, plus every rate.
3. **Independent label-free reviewers whose MUST items changed the record.**
   - lockbox_code_review MUST-1 to MUST-3 (pinned inputs, persist/hash/resume with a close guard, complete amendment record) were closed in 5082cc8 and A5.
   - The results review found the 16.5% silent BENCH bootstrap failures and the "beats matched sets" over-reading, both turned into reporting rules (D14).
   - The discovery review's M1 (GPD sensitivity next to the discoveries) was adopted (D8).
4. **Negative and null results reported as such, with pre-declared sensitivities.**
   - The bridge null (D5), "no calibration response fired, so context associations are mostly not complex-specific" (D8), BENCH "complex does not beat its best gene", and BRIDGE-L favouring GLS (D12/D14) were all kept in the headline.
   - The post hoc power aid (D14) stopped the matched-set gap from being sold as complex-specificity.
5. **Intent-before-compute for late-found gaps, plus PI re-runs of agent output.**
   - The robustness intent fd03b03 named the candidates and the rule before any computation, giving 254 findings (D16). The PI re-ran the annotations and reproduced them (dfef529).
   - This is the same discipline as the lockbox, applied cheaply, and the post-lockbox timing was disclosed.
6. **Mechanical verification of the written record.**
   - 579/579 manuscript numbers matched to source comments (6d3a11e). The independent reference check applied edits E1-E15 (a303165).
   - The M3 hygiene complaint was fixed: `git branch -a` shows one branch, where there were six worktree branches at M3.
   - The firewall allowlist for the two review scripts was added as a guard change with a commit (52b61a3), not as prose.

## Proposed pruned LESSONS.md (15 rules; was 16)

```
# LESSONS.md - process rules (rule + evidence). Pruned at M6 (docs/retros/retro_M6_final.md).
1. Non-evaluating agents and reviewers get no test- or lockbox-label access; every brief says "compute no metric, score or count on test pairs", including top-k/score-threshold tables and run logs that print them. (Part 1 D1; project D1; BRIEF_round1/2.md)
2. Fix PYTHONHASHSEED and all seeds for graph/community code; verify determinism by regenerating in a subprocess. (Part 1 Louvain tie-break bug)
3. Before freezing any gate, CI scheme or stochastic test, run a power/size simulation. Check the CI is well-defined (B, resampling unit, size matching). Assert error rates over >=50 sims against the guarantee the registered procedure gives, never a single seed, and match code defaults (B, tail rules) to the registry text. (Part 1 D2; round3_final B1; M3 GPD tail, 02f62ee)
4. Verifiers get only protocol + raw inputs, re-implement from raw data, sweep neighbouring configs, and commit their numbers before seeing the harness numbers. (Part 1 a2f8622/460fd36; c69d828 then 8d009f2)
5. Pre-register the exact look command; commit the "started"/intent record before reading labels; one look; a relook is non-confirmatory and runs only through code that refuses without the committed intent. (Part 1 IFT look 1; ed065f4; 29ab061)
6. Update the tree in the same commit as each result; grep any "disclosed/see file X/done" claim against the file before writing it; delete superseded status text. (Part 1 stale tree; round1_pi M1; M3 stop 5; the 2026-10-28 M3 header still in RESEARCH_TREE.md at M6)
7. Environment gotchas: run torch and LightGBM in separate pytest processes (OpenMP deadlock); wait on background jobs by PID (kill -0), never pgrep -f on a pattern in the waiter's own command line; set PART1_DIR explicitly in worktrees. (Part 1 final suite; M3 watcher loops; fzc_review)
8. Before any gate or merge, run the whole suite including cross-part guards, on a clean tree or pinned env; isolate concurrent agents (worktree/own venv), never share a gated tree, and merge or delete agent branches at each merge. (round3_final untracked files; M3 firewall pattern vs Part 1 guard; 52b61a3; six stale branches at M3)
9. Draw the data-flow order (who is built on which split, who reads the lockbox when) before hardening code. (round3_final B2/B3; round1_ai_failure MUST-7)
10. Reviews: independent roles, MUST/SHOULD/COULD, label-free, and they must execute code paths (tests, gate calls, scratch-repo simulations, the real freeze edit), not only compare code with the registry in tables. (round2_ai_failure M-B1; lockbox_code_review missed the covariate defect)
11. Commit compact summaries + hashes, not repeated generated JSON; one source for version labels; for mechanical per-entry steps commit one batch per run, and always commit the run log of any irreversible run. (958cfda 12k-line diff; 125/150 ceremony commits; "loader calls inferred, not observed")
12. Every firewall is a guard or test, not prose, and the prose is reconciled with DESIGN before the freeze commit. Register the metric, its guarantee and a chance baseline before the test or freeze that uses them. (D5 contexts.py read; D6 pooled FDP vs BB; retro_M3)
13. Registry<->code traceability table before the first gated run: every step in the registry AND every DESIGN rule (calibrations, responses, seeds, secondary estimands, robustness checks, lockbox driver) maps to a code symbol, a test and an output path, and a script fails on any unmapped step. Registry text that is ambiguous or inconsistent gets an amendment before the run. (D8 R4 calibrations; D9 driver gaps, A3/A4; D15 LOSO/leave-top-2; D8 E1 sign typo)
14. A dress rehearsal passes only on asserted outcomes: zero unexplained failed fits, fit-status counts equal to the expected eligible counts, and the lockbox path reproducing the A-side path on identical rows (differential test). After any defect, scope the remediation by a label-free footprint audit before writing its intent. (D12: 53 rank-deficient fits missed by the rehearsal; D13: footprint 560, not 53; 68f0d3f)
15. Writing agents get a word budget, registered sources only, and the number-check script at draft time. Post hoc annotations are descriptive columns; any rule that picks headline examples from them is written as a reporting rule first. (D15: 22k words cut to 6.3k; D16 261/29 vs registered 254; D17 label rules; 579/579 check run late)
```

### Changes and reasons

| Change | Reason |
|---|---|
| 1 sharpened | The lockbox adds a second label surface; run logs and diagnostics are metrics in disguise. |
| 3 = old 3 + 12 merged | Both are "simulate and test against the guarantee before freezing"; frees a slot. |
| 4 sharpened | Names the commit-numbers-first mechanic, which worked (c69d828/8d009f2). |
| 5 sharpened | Adds the code-enforced relook that worked (R6/R7), not only a procedure. |
| 6 sharpened | M3 stop 5 recurred; deleting dead status text is now part of the rule. |
| 8 = old 8 + 13 merged | Whole-suite-before-merge and clean-tree/isolation are one pre-gate checklist. Adds the branch hygiene that was fixed after M3. |
| 10 sharpened | The pre-read review compared tables and missed the covariate defect; reviews must run the path. |
| 11 sharpened | Adds the ceremony-commit batching and the uncommitted run log (D12 review). |
| 12 = old 14 | Unchanged in substance; moved after 11 so the registry rules sit together (12-14). |
| 13 = old 15 rewritten | Lesson 15 failed to prevent the D15 recurrence because it covered only "the engine" and registry entries. It is now a table that a script enforces, and includes DESIGN rules. |
| 14 = old 16 rewritten | The rehearsal missed the driver defect because it checked completion. Adds a differential A-path vs L-path test and the footprint-audit step from R6 to R7. |
| 15 new | Covers M6 failures (22k-word draft, unregistered filters, late number check, annotations becoming selection) that no old rule addressed. |
| 7, 2, 9 unchanged | Still valid and cheap; no new evidence to change them. |
