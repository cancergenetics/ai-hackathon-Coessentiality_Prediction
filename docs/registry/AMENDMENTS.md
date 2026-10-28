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
| Applied by | `scripts/40_amend_A1.py --timestamp 2026-10-28T22:54:50Z` |

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

At FZ-X, `scripts/39_freeze_fzx.py` read only the missingness pattern (`notna`) of A's gene-effect matrix. It used
this for the complete-case evaluability counts of H3a, H3b, BRIDGE-L and the controls, and discarded the values.
The H1 enumeration uses no dependency data.

**Files.**
- `results/part2/fzx/h1_event_split_A1.tsv`
- `results/part2/fzx/units_h1_components_A1.txt` (H1, hash `4a5dcc9d02bc4b23`)
- `results/part2/fzx/units_h1s_components.txt` (H1S, hash `d6264047856ff8e9`)
- `results/part2/fzx/lockbox_testable.tsv`: the silencing events are now listed under H1S.
- `FZX_manifest.json`: the `amendments` record, plus the updated `files`, `summary.m` and `complex_set_hashes`.

The original FZ-X lists (`h1_events.tsv`, `h1_contexts.tsv`, `units_h1_components.txt`) are unchanged.
