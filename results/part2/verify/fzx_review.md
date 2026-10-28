# FZ-X independent review

Reviewed: branch `worktree-agent-a741c9e5631821aca`, commits `408593e` (FZ-X) and `84a0c43` (amendment A1).
I did not produce this work, and I modified no producer file. My own code is `scripts/verify/fzx_independent.py`;
its outputs are in `results/part2/verify/fzx/`.

## Verdict: **ACCEPT-WITH-NOTES**, conditional on M1 and M2 before the H1 / H3b runs

Everything I could recompute independently reproduces exactly:
- the unit sets and the three set hashes;
- every context count;
- m for H3a, H3b and BRIDGE-L;
- all 341 H1 events, with their fields;
- the A1 split, its hashes and its m values.

The firewall holds. A fresh-process re-run is byte-identical. Neither MUST item changes m, a set hash or the FZ-X
event list, so the FZ-X lists can stand. Both items do change what will be run, so each needs a PI ruling or a fix
before the run it affects.

## MUST fix (before the run named)

**M1. SMARCB1 is typed silencing and A1 therefore moves it to the exploratory H1S. Fix before the H1 run.**

`Z0109|SMARCB1` has 17 A loss models:
- 12 carry a LoF call. DESIGN §2.2 defines LoF as "damaging == 2 OR damaging ≥ 1 + copy loss OR homozygous
  deletion", so hom-del lines count here.
- 7 are hom-del calls.
- 5 are silencing only.

So 71% of the loss is genetic. `h1.loss_calls` assigns each model one kind, with hom-del first, then LoF (mutation
only), then silencing. Under that partition, LoF is 5/17 and hom-del is 7/17. Neither reaches 50%, so rule 6 falls
through to silencing-type, and A1 moves SMARCB1 into H1S.

SMARCB1 is one of DESIGN's named LoF genes (it is the PC15/PC22 context). Rule 6 is ambiguous when LoF and hom-del
calls overlap. Under the literal §2.2 definition, SMARCB1 is LoF-type. Under the same reading, MAP2K4 (13/16),
PLAA (22/26) and FOCAD (31/46) would be LoF-type rather than CN-type. That changes their covariates and drops the
±10 Mb N exclusion.

SMARCB1 is the only silencing-typed gene whose loss is at least 50% genetic. Of the other silencing genes, 2 have a
genetic share of at least 30%, and 280 of 314 have none. A1's premise ("silencing = lineage non-expression") does not
hold for SMARCB1.

Needed: a PI ruling on rule 6 for overlapping calls, at minimum for SMARCB1, recorded as an amendment before any H1
run. If SMARCB1 moves into H1 as LoF-type, H1 has 24 events and m = 28, and H1S has 317 events and m = 338.

**M2. The registered gametolog removal misses RBMX. Fix before the H3b run.**

The registered rule for LoY and sex_female removes "every X-linked member that is an Ensembl-116 paralog of a chrY
gene". `cscore.Annotation.x_gametolog` looks up a paralog's chromosome only among DepMap gene-effect columns. The
chrY paralogs of RBMX (RBMY1A1, B, D, E, F and J) are not screened, so RBMX is kept.

This affects four (context, unit) tests:
- LoY and sex_female on `Z0104` (H3b);
- LoY and sex_female on `G0100` (BRIDGE-L).

In those tests, n_used is 12 / 8 instead of 11 / 7. The complete-case n, the evaluable status and m are unchanged.
The frozen `evaluability_tests.tsv.gz` records n_removed for these four rows from the faulty rule.

Fix: take the paralog's chromosome from Ensembl coordinates. A note in the freeze record for those four rows is
enough; m does not change.

## Checks

### 1. Outcome firewall: passes

- **Gene effects.** `39_freeze_fzx.py` calls `assoc.load_A_effects()`, which is A-only and raises on non-A rows. It
  keeps only `.notna()`, deletes the values and then uses column names. `fzx.evaluability` and `control_rows` use
  only the boolean mask, the context labels and the covariate metadata. The HC3 rank check uses X only.
- **Complete-case evaluability.** The complete-case part of evaluability is registered, in every entry's
  `statistic`: "≥ 12 context and ≥ 12 reference A models remain after complete-case". This is missingness only.
  ScreenNNMD and CasActivity enter as registered covariates, through their presence only.
- **L models.** L is used only for label counts.
- **Proteomics.** The LoF functional check reads CCLE proteomics and expression over U. This was pre-registered in
  DESIGN §2.2, and FZ-LB used the same scope.
- **Ledgers.** No `run_association` call and no ledger write. `results/part2/trials.jsonl` and
  `results/part2/lockbox_ledger.jsonl` exist neither on the branch nor on main.

### 2. Independent recompute: all match

My script does not use the producer's fzx, h1, cscore, covariates or contexts code. It reads Ensembl paralogs by
ENSG and midpoints from the coordinates file.

**Hashes**

| Unit set | Units | Set hash |
|---|---|---|
| zavg | 367 | `53f9dbf46370d237` |
| GLS | 379 | `14b6c53a01146a61` |
| H2 | 109 (51 parents) | `141a70b7133ed8a4` |
| H1 (after A1) | 23 events / 27 components | `4a5dcc9d02bc4b23` |
| H1S | 318 events / 339 components | `d6264047856ff8e9` |

**Contexts**
- All registered n_A, ref_A, n_L, ref_L and lockbox_testable values for H3a, H3b, BRIDGE-L and H2 equal my recount
  from the context table with L = U \ A (1,345 = 1,007 + 338).
- H3a has 22 lineages with ≥ 12 A models, plus melanoma. This matches the registry.
- The FZ-R exclusions were recounted: CREBBP 9, KMT2A-r 10, BCR::ABL1 11, PAX3/7::FOXO1 7 and SS18::SSX 5.

**m**

| Family | m |
|---|---|
| H3a | 8,441 |
| H3b | 7,340 |
| BRIDGE-L, H3a contexts | 8,717 |
| BRIDGE-L, H3b contexts | 7,577 |

- All 32,078 (context, unit) rows agree on n, n_context, n_ref and evaluable. n_used differs only in the four M2
  rows.
- The three BRIDGE-L tests that are not evaluable are `MYC_amp × G0137` and LoY / sex_female `× G0196`.
- The smallest complete-case n_context values are KMT2D 12, EWSR1 15 and STAG2 16. No context drops below 12.
- KMT2D and the EWSR1 fusion have n_L = 4 and take the rare-context route.
- No LoF context reaches a hom-del share of 0.5 (the largest is PTEN at 0.17).

**H1**
- Rules 1–8 give the same 341 events: 20 LoF, 3 CN and 318 silencing. Type, n_P, n_N, n_A, ref_A, n_L, ref_L,
  testable, locus and link all match. 124 events are in locus clusters, 123 of them silencing; 25 have a P axis.
- After the A1 split, 21 of 23 H1 events and 253 of 318 H1S events are lockbox-testable.

**Spot checks (10 events)**

| Event | Type | A loss models | Breakdown |
|---|---|---|---|
| RB1 (attached) | LoF | 60 | LoF 49 |
| STK11 | LoF | 43 | LoF 34 |
| BAP1 | LoF | 14 | LoF 12 |
| KEAP1 | LoF | 13 | LoF 12 |
| MED12 | LoF | 12 | LoF 10 |
| BTK (attached) | silencing | 26 | 26 silencing; LoF check failed |
| EPAS1 (attached) | silencing | 55 | 55 silencing |
| NOTCH2 | silencing | 19 | 19 silencing |
| HOMER3 | silencing | 20 | 20 silencing; locus |
| PSD4 | silencing | 15 | 15 silencing |

- All five LoF genes passed the FZ-LB functional check.
- Every attached link is a level-1, non-locus attachment.

### 3. Producer resolutions

- **Registry hash.** `complex_set_hash` is `registry.set_hash` (16 hex characters over the sorted unit IDs), while
  the manifest holds the file SHA-256s. This is consistent with the registry README and with the engine check
  (`assoc.run_association` compares set_hash with the tested unit columns). OK.
- **DepMap column keys, 8 symbol mismatches.**
  - Keying by Entrez → DepMap column is correct for the engine.
  - Note: `Annotation` looks up paralogs by DepMap symbol. For TIPE2 and TSLIG3A, B and C (HGNC TNFAIP8L2 and
    FAM98A, B and C), Ensembl knows only the HGNC names, so their paralogs are missed. This has no effect on any
    FZ-X list, because my ENSG-based recompute matches. Fix it in the same pass as M2, because it matters for the
    matched-null exclusions.
- **Attachments.** Only level-1, non-locus attachments are used (485 of 2,563). This follows FZ-C D3 (locus
  attachments "not used to define H1 events") and rule 4 ("attached to frozen coarse complex C"). Including the
  1,434 locus attachments would give 552 events. OK.
- **Missingness.** The mask is used as registered (see check 1). OK.
- **LoF check on the 893 extra genes.** All 893 have only 1–7 MSS LoF lines. The pre-registered check uses
  `min_lof = 8`, which made them ineligible at FZ-LB. FZ-X ran the same test with `min_lof = 0`, and 130 genes pass.
  - This is the only reading under which "genes not in that table get the same check" has content, so I accept it.
  - Impact if the extra passes are ignored: 337 instead of 341 events. ARHGAP21, SATB1, STRN and DIAPH3 drop out,
    and 5 n_A values shift by 1–2 (ALDOC, CDH2, YAP1, BCR, ACSL4). All of these are silencing events, so they fall in H1S after A1.
- **Script-31 hash.**
  - The FZ-C manifest still records the old SHA-256 of `scripts/31_score_blocks.py`, which is now the only FZ-C
    manifest mismatch. Every FZ-C data artefact still matches its manifest.
  - The change adds a guard and nothing else: `assert_A_only`, called before scoring, is a no-op for A run sets.
    No frozen FZ-C artefact changes.
  - FZX_manifest records the new hash. Anyone rebuilding FZ-C should know that the recorded hash refers to the
    version before the guard.
- **Test edits.** The placeholder asserts became asserts against the frozen lists. The FZ-R pending form is rebuilt
  synthetically, and `test_fzx_fill_changed_only_placeholders` still checks the fill-only rule against git history.
  None of the edits weakens a check. `test_part2_fzx.py` and `test_registry_entries.py` pass (32 tests).

### 4. Amendment A1

- **Only non-fill change.**
  - Against main `db75052`, a textual diff of every registry YAML removes only placeholder and `pending` lines.
    The one exception is the H1 header comment, which A1 rewrites.
  - `registry.amendment_changes` for H1 and H1S lists only the declared paths, and the test enforces this.
  - The other entries are fill-only.
- **Counts and types.** The counts and types in the rationale are accurate: 318 silencing events, 23 LoF/CN
  (20 + 3), m from 366 to 27 and 339.
- **"Mostly lineage-restricted."** This is plausible, but I did not quantify it. M1 shows one clear exception.
- **Enforcement.**
  - `registry.validate` requires `headline: false` when `standing: exploratory`.
  - `test_h1s_is_exploratory_copy_of_h1` asserts exploratory and false.
  - The set-hash check stops H1 from being run on H1S units.
  - `controls.kill_tally` reads only the control entries, so H1S results cannot reach the kill rule.
- **Gaps.** These are notes, not leaks today:
  - No strongest-findings, claims or figure code exists yet that filters on `standing`.
  - H1S has `direction: H1`, so trials and ledger lines will carry H1. Any future summary by direction would merge
    H1 and H1S, and must key on the entry id or `standing`.
  - H1S inherits the H1 `fdr_rule` text verbatim: "BH ... pooled over every counted H1 test (all events)". The
    binding family is `family_m.value = 339`, an entry of its own. The wording should say "this entry".
  - Stale text remains in both entries: `expected_events "about 10-20"` and `implementation_note "...stays pending
    until then"`. The pairing is implemented: `h1.assoc_inputs` and the `pairs=` argument of `run_association`.

### 5. Controls and kill rule

| Control | Finding |
|---|---|
| PC7 | SMARCA2 is in no coarse unit |
| PC8 | STAG1 is in no coarse unit |
| PC18 | KIF18A is in no coarse unit |
| PC21 | CDK6 is unclustered. CCND3 is in Z0113 and CDK4/CCND1 in Z0037, with no sub-module containing any of them. |
| PC22 | BRD9 is in Z0186 (ncBAF) and ARID1A, ARID1B and DPF2 are in Z0109 (cBAF), with no sub-module. The map separates ncBAF and cBAF at the coarse level, so the sub-module contrast is impossible. |

All 15 complex-evaluable controls have a frozen unit containing their target, with complete-case counts of at least
12 / 12.

Kill-rule denominators:

| Direction | Evaluable controls |
|---|---|
| H3a | 2 / 2 |
| H3b | 10 / 11 |
| H1 | 3 / 5 (PC6, PC12, PC15) |
| H2 | 0 / 2 |

The registered rule handles all four:
- H2 with 0 evaluable controls is "n/a, labelled exploratory", as written in the H2 `kill_rule` and in
  `controls.kill_tally`. H2 loses its headline standing; the paper must say so.
- H3a, H3b and H1 apply as written. For H1 a single recovery is not enough: 1 of 3 is below 50%, so H1 is killed;
  2 of 3 is needed. PC15 is not lockbox-testable (n_L = 4), but the kill rule is A-based, so this does not matter.

### 6. Re-run

- **FZ-X.** I ran `scripts/39_freeze_fzx.py` in a fresh process, on a detached worktree at `408593e` with an isolated
  `data_cache/part2/fzx`. Every list and every registry YAML is byte-identical, and so is `h1_loss_table.parquet`
  (SHA-256 `89f6c359…`). FZX_manifest differs only in `created`, `parent_commit`, `runtime_s` and `timings_s`.
- **A1.** `scripts/40_amend_A1.py --timestamp 2026-10-28T22:54:50Z`, run on `84a0c43`, reproduces every file,
  manifest included, byte for byte.
- **Refusal.** `39` refuses to run once the freeze is amended.

## Notes (non-blocking)

- **Branch base.** The branch is based on `8f63b6e`, not main `db75052`. Main added `scripts/39_map_variability.py`,
  so after the merge two scripts carry the number 39. The merge has no file conflicts, and the FZ-C manifest is
  unchanged on main.
- **"Identical files" claim.** The claim that two runs give identical files holds for the lists, not for the
  manifest's timestamp fields. `created` is local time with no timezone.
- **Superseded hash.** The FZX_manifest `code` entry for `src/part2/fzx.py` was overwritten by A1 with the
  post-amendment hash. The FZ-X-time hash survives only in commit `408593e`.

## Not verified

- I did not re-implement the LoF functional check (Mann-Whitney on expression and proteomics) for the 893 extra
  genes. I used the producer's `h1_lof_check_extra.tsv` pass flags and checked only their scope and impact.
- I did not re-implement the HC3 full-rank check. With n ≥ 12 / 12 everywhere, no design failure was reported.
- I did not quantify the "lineage-restricted" share of the silencing events.
- I did not run the full test suite, only `test_part2_fzx.py` and `test_registry_entries.py`.
