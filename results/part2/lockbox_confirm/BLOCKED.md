# Lockbox confirmations: STOPPED at the dry check (no read)

> **Update (same day, still no read).** The PI ruled on this note, and the missing machinery has been built:
> - the driver `scripts/51_lockbox_confirm.py`;
> - BENCH secondary and the BRIDGE-L GLS runner and multiverse;
> - amendment A3, which fixes the lockbox seeds and declares BENCH `within_line_rank` not run.
>
> The registry rulings still open are listed in `GAPS.md`. The driver refuses to read until they are recorded.

Date: 2026-10-29. Code: main `75bf5e6`. Agent worktree `worktree-agent-a3ec5565842d4e968`.

**Nothing was read from L.**
- `results/part2/lockbox_ledger.jsonl` does not exist.
- No target list was written, and no entry was closed.
- No intent note was committed.

The rule given for this run was: "if the machinery cannot perform a registered step, STOP and report". The dry
check found registered steps that the machinery cannot perform (section 2), so the run stopped here.

## 1. Dry-check items that pass

| Check | Result |
|---|---|
| Entries load | All 30 entries load through `registry.load_entry`: committed, clean, no placeholders, all `open`. 28 are confirmatory for read order: 23 controls, H3a, H3b, H1, H1S, H2. |
| Read order | Simulated in a scratch git copy of the committed registry, with zero-target lists and an outcome loader that raises if called. BENCH before any read is refused, and the refusal lists all 28 entries. BRIDGE-L after H1 is refused with `missing ['H1S_silencing_GI', 'H2_submodule_divergence']`. A second H3a read is refused (one read per entry). BENCH and BRIDGE-L are accepted only after all 28 entries have a ledger line and a committed `status: closed`. Order *among* the confirmatory entries (controls → H3a → H3b → H1 → H1S → H2) is not enforced by code, and the README does not require it to be. |
| Confirmation rule | `confirm.lockbox_test` implements the rule. The test is one-sided in the A direction (`less` if beta_A < 0). A hit is confirmed on same sign and p < 0.05, with ≥ 5 context and ≥ 5 reference L models. Inference is HC3 with Satterthwaite df. An exact PatientID-blocked Freedman-Lane permutation is used below 20 context models in L: within lineage for non-lineage contexts, unrestricted for lineage contexts (`lineage_contexts`). L 95% CIs are reported for every hit. H1 L scores (`h1.score_L_fn`) use the same members and A baselines. H2 (`h2.score_L_fn`) uses A-fitted coefficients. `confirm.rare_route` implements the rare-context route: LOSO same sign, DEMETER2 on A same sign, and leave-one-model-out same sign. All 49 control rows with an A estimate have A betas in the registered one-sided direction. |
| E2 | `confirm` has no label logic. `classes.tsv` keeps `label_final` "DETECTABLE, WITHIN EQUIVALENCE MARGIN" for PTEN and TP53, and `headline_eligible` = ARID1A, EP300, FOCAD only. A confirmation therefore cannot promote these events, provided the report reads labels from `label_final` and `headline_eligible`. |
| Frozen testable list | The SHA-256 of `results/part2/fzx/lockbox_testable.tsv` equals its entry in `FZX_manifest.json`. The file was last changed in A2 `fb264b8`, before any outcome was read. Not lockbox-testable: H3b/H2/BRIDGE-L LoF:KMT2D and EWSR1 fusion (n_L 4); controls PC15, PC16, PC22; H1 BRCA1 and EP300 (EP300 is a discovery). |
| Matched-set flag | Recomputed from `calibration_D8_tests.tsv.gz`, the flag matches the committed counts: H3a 502/1,145 (62 with no sets), H3b 103/161 (7). H2 beats all of its own splits in 40/109. This H2 count was computed here from the split tests and is not in the committed summary. |

## 2. Registered steps with no implementation (reason for the stop)

1. **BENCH secondary estimand.** The registry specifies an "A–L correlation of context effects (residual-SD units) for complexes vs random gene sets matched on size, essentiality and coherence (seed 7100); number of distinct complexes among the g*". No code computes it. `bench.py` covers only the primary Delta, g* and the bootstrap.
2. **BENCH multiverse `within_line_rank`.** `cscore.score` ranks every gene within a line. The registered BENCH read, however, covers only "all members of the complexes", so this score cannot be computed from the single registered read. The multiverse eigengene score with an A fit, g* by |beta| and the T1/T2 subset are not wired.
3. **BRIDGE-L.** The registered GLS-map discovery pass on A (two `run_association` calls with BB at q 0.05, Freedman-Lane re-entry and seed 4500) has no runner: `44_run_direction.py` accepts only the five direction IDs. The engine itself supports it. The multiverse ("unweighted rate, matched gene coverage subset, T1/T2 clusters only") is not implemented.
4. **No lockbox runner.** No script calls `confirm.*` on real data. The only end-to-end use is `tests/part2_synth_pipeline.py`, and it omits several pieces:
   - the rare-context route;
   - H1S;
   - the controls' CN hom-del removal share;
   - L 90% CIs for neutral H1 calls;
   - routing by the frozen testable list;
   - BRIDGE-L.

   Also, `lockbox_test` decides testability from the L complete-case counts, not from the frozen list, and `confirm_entry` never calls `rare_route`. A driver therefore has to route hits itself.
5. **Lockbox permutation seeds are not registered for the directions.** The `seeds.permutation` values of H3a, H3b, H1, H1S and H2 are their discovery seeds. CTRL 4400–4422 were unused in A (`n_perm = 0`).

Items 1–3 need a PI ruling before any read: implement and review, or declare as not run, as was done for D8 R4. Item 4 is new orchestration code that would run once against L. It should be reviewed, or at least pass an end-to-end synthetic run that covers every route, before the read.

## 3. Runtime estimate (if unblocked)

| Step | Estimate | Basis |
|---|---|---|
| 28 confirmatory reads | about 20–30 min | About 770 hits need L permutation (H3a 701, H2 53, H3b 17, plus controls and H1); 10,000 permutations at n ≈ 338 take about 0.3 s per hit. One CRISPR load per read. |
| GLS discovery pass | about 10–12 min | H3a 451 s and H3b 199 s analogues |
| BENCH primary bootstrap | about 4 h | 1,295 primary pairs (H3a 1,145 + H3b 161 − 11 EWSR1) × 2,001 fits; `bench_delta` measured at 5.4 ms per pair per resample on synthetic data |
| BRIDGE-L | under 1 min | |

Total about 4.5–5 h, with peak RSS of about 3–4 GB, as in the discovery runs.
