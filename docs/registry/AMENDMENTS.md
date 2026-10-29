# Registry amendments

An amendment is a documented change to a committed, filled registry entry. It is not a freeze fill.

Each amendment has two records:
- a record here;
- an `amendments` item in every entry it touches, with fields `id`, `timestamp_utc`, `amends_commit`, `reason` and
  `changed_paths`.

`part2.registry.amendment_changes(committed_version, new_version)` must return `[]`: every changed field has to be
declared in `changed_paths`. `tests/test_part2_fzx.py` checks this against the amended commit. `fill_only_changes`
is unchanged and still governs freeze fills only.

## A1: H1 silencing events split into the exploratory entry H1S_silencing_GI

| | |
|---|---|
| Time (UTC) | 2026-10-28T22:54:50Z |
| Amends | FZ-X commit `408593e43c067e0a213e2dfadabb5c48b68bfbb2` |
| Decided by | PI, before any H1 run |
| Applied by | `scripts/42_amend_A1.py` (numbered `40_amend_A1.py` at the A1 commit; renamed with the review fixes, note c) |

**Decision.**
1. `H1_within_complex_GI` is the primary, headline panel. It keeps only the LoF-type and hom-del (CN-type) events of
   the FZ-X list: 23 events (20 LoF, 3 CN). m changes from 366 to 27, and the complex_set_hash is recomputed.
2. The new entry `H1S_silencing_GI` holds the 318 silencing-type events.
   - It uses the same estimand, score, classes, covariates, TOST and matched-set rules, and seeds. The
     matched-null seed is 7000 + the event index in `results/part2/fzx/h1_events.tsv`.
   - It has its own BH family, with m = 339.
   - Its standing is `exploratory` and it carries `headline: false`. It is never headlined; wherever it is reported
     it is labelled exploratory.
   - Its lockbox confirmations are run and reported. They never enter the kill rule or the main claims.
   - Its lifecycle `status` stays `open`/`closed`, as for every entry, because the lockbox read order needs it.

**Rationale (verbatim from the PI decision):**
- DESIGN v3.1 §H1 anticipated a panel of about 10–20 events. The registered silencing clause (TPM<1 where the lineage median is ≥3) yields 318 events, mostly lineage-restricted non-expression.
- Such events are confounded with lineage composition.
- A non-expressed gene is trivially non-dependent.
- At m=366 the intended LoF/CN panel would be diluted.
- The amendment was decided from event counts and types only. No dependency values, association statistics or lockbox outcomes were read.

**Confirmation by the implementing agent.** This is true. The decision used only the FZ-X event list
(`h1_events.tsv`: event counts and types). No H1, H1S or any other association has been run:
`results/part2/trials.jsonl` and `results/part2/lockbox_ledger.jsonl` do not exist. No lockbox outcome has been read.

At FZ-X, `scripts/41_freeze_fzx.py` (then numbered `39_freeze_fzx.py`) read only the missingness pattern (`notna`) of A's gene-effect matrix. It used
this for the complete-case evaluability counts of H3a, H3b, BRIDGE-L and the controls, and discarded the values.
The H1 enumeration uses no dependency data.

**Files.**
- `results/part2/fzx/h1_event_split_A1.tsv`
- `results/part2/fzx/units_h1_components_A1.txt` (H1, hash `4a5dcc9d02bc4b23`)
- `results/part2/fzx/units_h1s_components.txt` (H1S, hash `d6264047856ff8e9`)
- `results/part2/fzx/lockbox_testable.tsv`: the silencing events are now listed under H1S.
- `FZX_manifest.json`: the `amendments` record, plus the updated `files`, `summary.m` and `complex_set_hashes`.

The original FZ-X lists (`h1_events.tsv`, `h1_contexts.tsv`, `units_h1_components.txt`) are unchanged.

## A2: genetic-loss rule for H1 vs H1S, plus review notes on the H1/H1S entries

| | |
|---|---|
| Time (UTC) | 2026-10-28T23:28:31Z |
| Amends | A1 commit `84a0c43d773e56bee1c91fdcdff249c574ba0555` |
| Decided by | PI ruling on FZ-X review M1 (`results/part2/verify/fzx_review.md`), before any H1 run |
| Applied by | `scripts/43_amend_A2.py`, re-applied by `scripts/41_freeze_fzx.py` after the base freeze and A1 |

**Ruling (from the PI's message).** DESIGN §2.2 counts hom-del as LoF, and A1 intends genetic loss to be the primary
panel and silencing to be exploratory.
- An event is **genetic** (`H1_within_complex_GI`, primary) iff ≥ 50% of its A loss models carry a LoF call or a
  hom-del call. The two are combined: any model with either counts. Otherwise the event is **silencing**
  (`H1S_silencing_GI`).
- Genetic events keep a descriptive sub-label, LoF-majority or CN-majority, for reporting only.
- The rule is applied mechanically to all 341 FZ-X events.

**How it is applied.**
- The counts are the frozen per-gene loss-kind counts in `h1_contexts.tsv` (n_A_LoF + n_A_homdel over n_A). Each A
  loss model carries exactly one kind: hom-del first, then LoF that passed the functional check, then silencing.
  Counting every LoF call, including those of genes that failed the functional check, gives the same split
  (checked).
- The sub-label is CN-majority if hom-del calls make up ≥ 50% of the genetic loss models, else LoF-majority.
- The analysis specification of every event is unchanged: covariates and N/P components follow the frozen FZ-X
  rule-6 `type`. Each context item now also carries `loss_class`, `genetic_share` and, for genetic events,
  `genetic_sublabel`.

**Result.** Exactly one event moves: `Z0109|SMARCB1`, from H1S to H1.
- It has 17 A loss models: 5 LoF, 7 hom-del and 5 silencing. Its genetic share is 0.706 and it is CN-majority.
- Its FZ-X type is silencing, so it keeps the base covariate set. No Z0109 member lies within ±10 Mb of SMARCB1, so
  its N component would be the same under the CN-type rule.
- The highest genetic share among the remaining silencing events is BCR at 0.467.

| Entry | Events | Sub-labels | m | Hash | Lockbox-testable |
|---|---|---|---|---|---|
| H1 | 24 | 20 LoF-majority, 4 CN-majority | 28 | `686659438f7d781c` | 22 / 24 |
| H1S | 317 | – | 338 | `736ed74bc58bc04a` | 252 / 317 |

**Also part of A2, from the review notes:**
- H1S has its own direction label `H1S`. It counts as confirmatory for the lockbox read order only; it is never
  part of the kill rule.
- The H1S `fdr_rule` now says the BH family is this entry only.
- The stale `expected_events` ("about 10-20") and `implementation_note` ("stays pending") texts of both entries are
  replaced.

**Statement.** The ruling was decided from event counts and types only. No dependency values, association
statistics or lockbox outcomes were read. This is true: no association has been run, and neither
`results/part2/trials.jsonl` nor `results/part2/lockbox_ledger.jsonl` exists. The A2 script reads only the frozen
event and count tables.

**Applied with A2 (review M2 and note a; not a registry change).**
- `cscore.load_annotation` now derives paralogs from Ensembl gene IDs (DepMap Entrez to ENSG through the Ensembl
  coordinates file). Chromosomes now cover every Ensembl gene, so X-Y gametologs whose partner is not screened are
  found (RBMX to RBMY1A1, B, D, E, F and J).
- Re-running the freeze changes only `n_removed` / `n_used` in four rows of `evaluability_tests.tsv.gz`: LoY and
  sex_female on `Z0104` (H3b) and on `G0100` (BRIDGE-L). RBMX is now removed, so n_used goes from 12 to 11 and from
  8 to 7.
- n, evaluable status, every m, every unit set and hash, the control units and the H1 event list are unchanged.

## E1: erratum, sign typo in the H1 / H1S aligned-score formula

| | |
|---|---|
| Time (UTC) | 2026-10-29T00:41:42Z |
| Amends | A2 version of both entries, commit `fb264b8bd997ebf27c934aaed4f793a1433522fb` |
| Decided by | PI ruling R1, after the H1 discovery run (`d12c627`) |
| Entries | `H1_within_complex_GI`, `H1S_silencing_GI`; `changed_paths: [score_rule]` |

**Typo correction, semantics unchanged, identified after the H1 run, implementation predates it.**

- The `score_rule` read `y_N(m) = mean over j in N' of -sign(b_j) * effect_j(m)`. For an essential member
  (b_j < 0), more dependency is a more negative effect, so the literal formula gives kappa_N < 0 when dependency
  grows. That contradicts the entry's own stated meaning ("kappa_N > 0 when dependency magnitude grows in its
  baseline direction") and its class rules (residual sensitisation iff kappa_N > 0).
- Corrected text: `sign(b_j) * effect_j(m)`. Nothing else changes.
- The engine has implemented the stated semantics (`h1.ALIGN = +1`, `y_N = mean_j sign(b_j) * effect_j`) since
  commit `201d8bfa1814d5d7fa7f7981c39fc1bf33363b85` (2026-10-28T20:09:59Z, "Part 2 association engine"), tested on
  synthetic data from `040e0df` (`tests/test_part2_h1h2.py::test_alignment_sign_and_collapse_index`), i.e. before
  any real-data run. The H1 discovery run (`d12c627`) therefore used the registered semantics.
- The H1 trial line carries the pre-erratum registry blob (`ebfbaad4…`); the post-erratum blob differs only in
  `score_rule` and the `amendments` record. H1S has not been run before E1.

## E2: reporting rule when a non-neutral H1 class and TOST-neutral both hold (post hoc, conservative)

| | |
|---|---|
| Time (UTC) | 2026-10-29T01:07Z |
| Decided by | PI ruling R5, after the H1 run and its independent review (branch review-discovery, `bff6cad`) |
| Entries | `H1_within_complex_GI`, `H1S_silencing_GI`; **reporting layer only, no registry field changed** |

**Post hoc and conservative.** The registry gives no precedence between a non-neutral class (residual sensitisation,
collapse, paralog compensation) and the TOST-neutral criterion (90% CI inside ±0.1 Chronos). When both hold, the
reported label is "DETECTABLE, WITHIN EQUIVALENCE MARGIN (negligible)". The engine's class stays in its own column
(`class_complex_wide` / `class_paralog`); `label_final` / `label_final_paralog` carry the reported label
(`scripts/47_h1_label_final.py`). Such events are not headline findings. They still go to lockbox confirmation on
their registered statistic, but a confirmation cannot promote them. Engine class logic, q-values and trial lines are
unchanged.
- H1: Z0130|PTEN (collapse) and Z0158|TP53 (residual sensitisation) are relabelled.
- H1S (exploratory): Z0149|ERBB3, Z0194|PRR15L, Z0260|B3GNT5 and Z0273|CST6 are relabelled.
- Also in the reporting layer: Z0164|MAP2K4 (H1) is labelled NOT EVALUABLE, because its counted tests were not
  evaluable (ruling R2, applied without re-running H1).

## A3: lockbox permutation seeds; BENCH `within_line_rank` declared not run

| | |
|---|---|
| Time (UTC) | 2026-10-29T02:03:27Z |
| Amends | `21fe1be1cdfe9d9b5b68c02c1b5a46b2b5b66a69` (registry as of main `75bf5e6`, plus the dry-check stop note) |
| Decided by | PI rulings on `results/part2/lockbox_confirm/BLOCKED.md` |
| Applied by | `scripts/49_amend_A3.py` (`--check` re-derives the edit from the amended commit) |
| Statement | **before any lockbox read; no L values read** |

The lockbox ledger does not exist, no target list has been written, and no lockbox outcome has been loaded.

**(a) Lockbox permutation seeds.** The registry fixed the discovery permutation seeds but not the seed of the exact
PatientID-blocked Freedman-Lane permutation that a confirmation runs in L when fewer than 20 context (or loss) models
are in L. A3 adds `seeds.lockbox_permutation`:

| Entries | Seed |
|---|---|
| `H3a_lineage_atlas` | 6100 |
| `H3b_driver_family` | 6200 |
| `H1_within_complex_GI` | 6300 |
| `H1S_silencing_GI` | 6400 |
| `H2_submodule_divergence` | 6500 |
| `BRIDGE-L_confirmation_rate` | 6600 |
| `CTRL_*` | 6700 + registry index (the index of the registered control seed 4400 + index): CTRL_PC1 6700 … CTRL_PC19a 6718, CTRL_PC19b 6719 … CTRL_PC22 6722 |

- `confirm.confirm_entry` uses `seed + i`, where i is the hit's position in the target list's hit order, which is
  sorted by (context, unit). This is the same convention as the discovery runs (`perm_seed + row index`).
- Collision check: every value under `seeds` / `calibration.seed` in every entry (4000–4422, 4500, 4600, 4700, 4800,
  5100, 5200, 7000, 7100) and every seed named in DESIGN.md, RESEARCH_TREE.md and docs/registry/*.md (3000, 5000,
  7000). No base seed collides, so no block was shifted by +50.
- Offset ranges overlap other bases, as they already do for the discovery seeds (for example, H3a 6100 + up to 1,144
  hits). They are independent streams on different data.

**(b) BENCH `within_line_rank` is NOT RUN.** `cscore.score(method="within_line_rank")` ranks each line's effect among
*all* genes. The registered BENCH read is "all members of the complexes in the primary and secondary pair sets", so
the entry is internally inconsistent. The PI chose not to widen the lockbox read. `BENCH_complex_vs_gene` now carries
`multiverse.not_run: [within_line_rank]` with a reason. This is reported as a limitation. The other multiverse items
(eigengene, g* by |beta|, T1/T2 pairs) are unchanged.

`registry.amendment_changes(committed, amended) == []` holds for all 30 touched entries (29 seed entries plus BENCH),
and every declared path changed (`tests/test_part2_lockbox_driver.py`).
