# Analysis registry (FZ-R)

One file per registered analysis: `docs/registry/<id>.yaml`. An entry must be **committed** before its first run. `part2.assoc.run_association` and `part2.assoc.read_lockbox` refuse an entry if:
- the file is missing or fails validation;
- the file is untracked, or modified relative to HEAD;
- it still holds a freeze placeholder (see below);
- `status` is `closed` (refused by `run_association` only).

## Registered entries (FZ-R draft)

29 FZ-R entries plus `H1S_silencing_GI` (amendment A1), all `open`, validated by `tests/test_registry_entries.py`:
- `CTRL_PC1` … `CTRL_PC22` (with `CTRL_PC19a`, `CTRL_PC19b`): the 23 pre-registered controls of DESIGN §2.7, one per file, each with its kill-rule direction (§2.6).
- `H3a_lineage_atlas`, `H3b_driver_family`, `H1_within_complex_GI`, `H2_submodule_divergence`: the four directions in registry order.
- `BENCH_complex_vs_gene` (direction `benchmark`) and `BRIDGE-L_confirmation_rate` (direction `bridge`): their lockbox reads are allowed only after every confirmatory entry has a ledger line.

## FZ-X fill

At FZ-X every placeholder was filled by `scripts/39_freeze_fzx.py`, and every `pending` key was removed. The values
come from the frozen lists in `results/part2/fzx/`, described in its README:
- unit sets and complex-set hashes;
- m per family;
- control units;
- the H1 event list.

`frozen_inputs[2]` is `results/part2/construction/FZC_manifest.json` and `frozen_inputs[3]` is
`results/part2/fzx/FZX_manifest.json`.

`registry.fill_only_changes` against the FZ-R drafts is empty for all 29 entries. `tests/test_part2_fzx.py` re-checks
this from the commit history.

Placeholder-free entries can now be run, in registry order.

## Amendments

A documented change to a filled, committed entry is recorded in two places:
- in `docs/registry/AMENDMENTS.md`;
- as an `amendments` list in the entry, each item with `id`, `timestamp_utc`, `amends_commit` (40 hex),
  `reason` and `changed_paths`.

`registry.amendment_changes(committed, new)` must be `[]`: only declared paths may change. `fill_only_changes` is
not loosened.

Optional fields `standing` (`confirmatory` | `exploratory`) and `headline` (bool) mark entries that are never
headlined. An exploratory entry must carry `headline: false`. These fields are separate from the lifecycle
`status` (`open` / `closed`), which the lockbox read order relies on.

A1 split `H1_within_complex_GI` into two entries:
- `H1_within_complex_GI`: the LoF and hom-del events;
- `H1S_silencing_GI`: the silencing events, an exploratory entry. Its direction is `H1`, so it is read before
  benchmark and bridge.

## Freeze placeholders and `pending`

The entries are drafted before the complex map (FZ-C) and the context/event freeze (FZ-X) exist. A value that can only be known there is written as the placeholder `FZ-C` or `FZ-X`, and the field's dotted path (list items as `[i]`, e.g. `frozen_inputs[2]`, `family_m.value`, `control.complex_unit`, `context_family.contexts`) is listed in `pending`. Validation requires `pending` to list exactly the placeholder fields. `load_entry` (and so `run_association` and `read_lockbox`) refuses an entry that still holds a placeholder; `load_all` accepts it, so a pending entry still counts for the lockbox read order.

Filling at FZ-C / FZ-X: replace each placeholder, shorten `pending` (drop the key when empty), and check that nothing else changed with `part2.registry.fill_only_changes(old_entry, new_entry) == []` before committing. Placeholders allowed: `complex_set_hash` (hex of `set_hash(units)`), `context_family.contexts` (`FZ-X` only, and only with a `context_family.enumeration_rule`), and any value in the optional fields below.

The git blob SHA of the committed entry is written to every trial line (`results/part2/trials.jsonl`) and every lockbox-ledger line (`results/part2/lockbox_ledger.jsonl`). The design behind this is DESIGN.md v3.1 §0 (F4/F5), §2.5 and §2.6.

## Fields (all required)

| Field | Type | Meaning |
|---|---|---|
| `id` | string `[A-Za-z0-9_.-]` | Must equal the file name without `.yaml`. |
| `direction` | `H3a` \| `H3b` \| `H1` \| `H2` \| `controls` \| `benchmark` \| `bridge` | Registry order is §2.6. The lockbox may be read for `benchmark` and `bridge` only after every `controls`/`H3a`/`H3b`/`H1`/`H2` entry has a ledger line (round-3 B3). |
| `hypothesis` | string (block scalar allowed) | The claim, in words. |
| `sign` | `"+"` \| `"-"` \| `two-sided` \| `per-context` | `+` means **more dependent in context**, i.e. a lower (more negative) effect score, tested one-sided with `alternative="less"` on the regression coefficient. `-` means less dependent (`greater`). With `per-context`, every context carries its own `sign`. |
| `context_family` | mapping: `name`, `contexts` (list of `{id, definition[, sign]}`) | Context IDs are column names of the frozen context table (`data_cache/part2/contexts/context_table.parquet`; summary in `results/part2/contexts/context_summary.tsv`). The definition is restated in words. `run_association` refuses any context column not listed here. |
| `complex_set_hash` | 16–64 lowercase hex | `part2.registry.set_hash(units)`: the first 16 hex characters of the SHA-256 of the sorted, newline-joined IDs of the tested units (frozen complexes, sub-modules or genes). `run_association` refuses a score matrix whose column set has a different hash. |
| `score` | `S_C_mean_winsor3` \| `eigengene` \| `within_line_rank` \| `expression_masked` \| `pan_essential_referenced` \| `gene` \| `submodule_vs_sibling` \| `member_resolved` | Primary score (§2.4). The multiverse scores are reported per hit and never selected. |
| `covariates` | list of strings | Covariate names as built by the caller, e.g. `lineage_FE`, `library`, `ScreenType3DO`, `ModelType`, `log_burden`, `MSIScore`, `Ploidy`, `WGD`, `TP53_altered`. Follow the §2.2 covariate rule: no covariate may be a function of the context. |
| `fdr_family` | mapping: `method` (`benjamini_bogomolov` \| `bh` \| `stagewise_bb` \| `none`), `families` (`contexts` \| `complexes`, for BB), `q` | §2.5: Benjamini–Bogomolov per direction. Families are contexts, selected by Simes + BH at q; within each selected family, BH at q·R/M. Permutation p-values re-enter and BB is re-run. |
| `success_criterion` | string | E.g. discovery q < 0.05 in A, plus lockbox same sign with one-sided p < 0.05 (contexts with ≥ 5 L models). |
| `stopping_rule` | string | E.g. a single run, one lockbox read, the kill rule of §2.6. |
| `frozen_inputs` | list of manifest paths | Freeze manifests (JSON with a `files: {path: sha256}` map) that must be committed. Every `inputs=` path passed to `run_association` must be listed in one of them with its current SHA-256. `results/part2/lockbox/manifest.json` is always required. |
| `status` | `open` \| `closed` | Set to `closed` (in a commit) after the entry's lockbox read. A closed entry cannot hold placeholders. |

Field extensions (FZ-R):
- `complex_set_hash` may be the placeholder `FZ-C` or `FZ-X`.
- `context_family.contexts` may be `FZ-X` (with `context_family.enumeration_rule`). A context item may also carry `covariates` (a list that overrides the entry list for that context; the §2.2 covariate rule is applied per context), `role` (`primary`, the default, or `sensitivity`: registered so it can be run for discovered hits, never part of the family or of m), `removed_members`, `type`, `n_A`/`ref_A`/`n_L`/`ref_L` and `lockbox_testable` (FZ-LB counts from `results/part2/lockbox/context_counts.tsv`).
- `fdr_family.method` may be `stagewise_bb` (H2; needs `parent_directions`). `run_association` runs it only with a `parent_selection` (`assoc.parent_selection_from` the parent directions' discovery tables) and only when every parent direction has a discovery line in the committed `trials.jsonl`; it tests only the eligible (context, sub-module) pairs, BH within each (context, parent) family at q·R_x/m_x.

## Optional fields (type-checked when present)

| Field | Type | Meaning |
|---|---|---|
| `pending` | list of strings | Paths of the fields holding placeholders (see above). |
| `family_m` | mapping with `formula` (string) and `value` | The number of tests in the family, as a formula fixed now and a value filled at FZ-C/FZ-X. |
| `statistic` | string | Test statistic and evaluability rule. |
| `control` | mapping | Controls only: `design_id`, `level` (G, G+C or sub-module contrast), `targets`, `gene_units` (DepMap column names), `complex_unit` (frozen unit ID, or "not evaluable at complex level"), `complex_unit_rule`, `readout_C`, `removed_members`, `design_sign`, `sign_mapping`. |
| `kill_rule` | mapping | Controls: the direction it counts for (`H3a`, `H3b`, `H1`, `H2` or `none` for gene-level-only controls), `counted`, `level`, `recovered_if`, `evaluable_if`. Directions: the mapped `controls`. |
| `lockbox` | mapping | Confirmation rule, rare-context route, target-list path. |
| `multiverse` | mapping | Scores, covariate and context sensitivities, robustness checks; reported and never selected. |
| `calibration` | mapping | Negative controls (pseudo-contexts, matched random sets) and the pre-declared response. |
| `seeds` | mapping | Seeds used by the entry (permutations 4000+, bootstraps 5100+, matched sets 7000+). |
| `run_after` | list | Entries whose results must be committed first. |

Other descriptive keys (`score_rule`, `covariate_rule`, `fdr_rule`, `units`, `estimand`, `classes`, `matched_null`, `paralog_definition`, …) are free text.

## YAML subset

The pinned environment has no YAML library, so `part2.registry.load_yaml` accepts a strict subset:
- block mappings (`key: value`) and block sequences (`- item`), including a sequence of mappings whose continuation keys sit two spaces deeper than the dash;
- flow sequences of scalars (`[a, b]`);
- block scalars (`|`, `>`);
- quoted strings, ints, floats, `true`/`false`/`null`, and `#` comments.

These are refused: anchors and aliases, tags, flow mappings, nested flow collections, multi-document streams, tab indentation and duplicate keys. Quote any `+` or `-` sign: `sign: "+"`.

## Example (validated by `tests/test_part2_assoc.py`; not a registered entry)

```yaml
id: EXAMPLE_H3b_driver_family
direction: H3b
hypothesis: |
  Driver / genome-state contexts are associated with complex-level dependency
  (S_C) of the frozen complexes; tested against all frozen complexes.
sign: per-context
context_family:
  name: H3b driver / genome-state family (fixed at FZ-X)
  contexts:
    - id: hotspot:KRAS
      definition: "KRAS hotspot (hotspot matrix >= 1)"
      sign: two-sided
    - id: MSI_H
      definition: "MSIScore > 20"
      sign: two-sided
complex_set_hash: 0123456789abcdef
score: S_C_mean_winsor3
covariates: [lineage_FE, library, ScreenType3DO, ModelType, log_burden, MSIScore]
fdr_family:
  method: benjamini_bogomolov
  families: contexts
  q: 0.05
success_criterion: "discovery q < 0.05 in A; lockbox same sign, one-sided p < 0.05 (n_L >= 5)"
stopping_rule: "one discovery run; one lockbox read; kill rule of DESIGN §2.6"
frozen_inputs:
  - results/part2/contexts/checksums.json
  - results/part2/lockbox/manifest.json
status: open
```

## Lockbox target lists

`read_lockbox(analysis_id, targets_file)` needs a committed JSON file, conventionally `results/part2/lockbox_targets/<id>.json`:

```json
{"analysis_id": "<id>", "genes": ["SYMBOL (ENTREZ)", "..."], "units": {"<complex id>": ["SYMBOL (ENTREZ)", "..."]}, "contexts": ["..."]}
```

`genes` are the outcome columns returned for the lockbox models. Before any outcome is loaded, the read is appended to `results/part2/lockbox_ledger.jsonl` and committed. There is one read per entry and per source (`crispr`, `demeter2`).

An entry with no discovery is closed with a zero-target list (`"genes": [], "zero_targets": true`): the ledger line is committed and nothing is read. The confirmation step (`part2.confirm`: `write_target_list` → `confirm_entry` → `close_entry`) writes `results/part2/assoc/<id>/confirmation_<source>.json` and sets `status: closed` in a commit that changes nothing else. `benchmark` / `bridge` reads are refused until every confirmatory entry has a ledger line **and** is committed with `status: closed`.
