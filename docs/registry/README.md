# Analysis registry (FZ-R)

One file per registered analysis: `docs/registry/<id>.yaml`. The PI writes the entries; this folder holds only this README until then. An entry must be **committed** before its first run. `part2.assoc.run_association` and `part2.assoc.read_lockbox` refuse an entry if:
- the file is missing or fails validation;
- the file is untracked, or modified relative to HEAD;
- `status` is `closed` (refused by `run_association` only).

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
| `fdr_family` | mapping: `method` (`benjamini_bogomolov` \| `bh` \| `none`), `families` (`contexts` \| `complexes`, for BB), `q` | §2.5: Benjamini–Bogomolov per direction. Families are contexts, selected by Simes + BH at q; within each selected family, BH at q·R/M. Permutation p-values re-enter and BB is re-run. |
| `success_criterion` | string | E.g. discovery q < 0.05 in A, plus lockbox same sign with one-sided p < 0.05 (contexts with ≥ 5 L models). |
| `stopping_rule` | string | E.g. a single run, one lockbox read, the kill rule of §2.6. |
| `frozen_inputs` | list of manifest paths | Freeze manifests (JSON with a `files: {path: sha256}` map) that must be committed. Every `inputs=` path passed to `run_association` must be listed in one of them with its current SHA-256. `results/part2/lockbox/manifest.json` is always required. |
| `status` | `open` \| `closed` | Set to `closed` (in a commit) after the entry's lockbox read. |

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
