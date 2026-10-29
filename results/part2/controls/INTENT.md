# Controls run: pre-registered intent (single look)

This note is committed before the run. The controls are run exactly once (DESIGN v3.1 §2.7, kill rule §2.6;
`docs/registry/CTRL_*.yaml`). Any re-run would be a relook and must be disclosed here.

- Code and registry commit: `996a9a3797e3312b384203c221be2f135a09e7bc` (main). This note is committed on top of it,
  with no code, registry or test change.
- Planned start: 2026-10-28, about 23:56 UTC. The actual start and end times go in `runtime.json`.
- Working directory: the agent worktree. `datasets/`, `data_cache/` and `external/` are symlinks to the main
  repo; they are gitignored and never committed.
- Command (wrapped in `/usr/bin/time -l` to record wall time and peak RSS; stdout and stderr go to
  `results/part2/controls/run_log.txt`):

```
export PART1_DIR="<PART1_DIR>"; export PYTHONHASHSEED=0
<main repo>/.venv/bin/python scripts/34_run_controls.py \
    --map results/part2/fzx/zavg_engine_membership.tsv \
    --context-table data_cache/part2/contexts/context_table.parquet
```

- Inputs:
  - The map is listed in `results/part2/fzx/FZX_manifest.json` with sha256 `6a06223345989157…`.
  - The context table is listed in `results/part2/contexts/checksums.json` with sha256 `ba73a37195a2…`.
  - `run_association` checks both inputs against the committed manifests.
- Dependency data: A models only, through `assoc.load_A_effects` (1,007 A models). The lockbox is not read and no
  `read_lockbox` call is made.
- Excluded: no direction is run (H3a, H3b, H1, H1S, H2), and no hu.MAP, STRING or curated complexes are used.

## Dry check

The dry check did not touch dependency data: it used gene-effect column names and model IDs from the cache
metadata, a synthetic N(0,1) outcome, and patched writers, so nothing was written.
- All 23 `CTRL_*` entries load: committed, clean, no placeholders, no `pending`.
- `check_frozen` passes for both inputs under every entry's manifests.
- Every unit set matches its `complex_set_hash`, and every gene unit and context column exists.
- `evaluate_controls` runs end to end: 57 rows and 23 trial lines.
- The complex-level kill-rule structure is:

  | Direction | Counted controls | Evaluable (by context size and membership) |
  |---|---|---|
  | H3a | 2 | 2 |
  | H3b | 11 | 10 |
  | H1 | 5 | 3 |
  | H2 | 2 | 0 |

- Every entry has `fdr_family.method: none` and `n_perm = 0`, so there are no permutations.
- Expected runtime is under 2 minutes. The synthetic evaluation took 16 s.
- Tests that pass at this commit: `test_registry_entries`, `test_part2_lockbox_firewall`, `test_part2_firewall`,
  `test_part2_h1h2` and `test_part2_cscore`.

## Outputs, committed in one commit "Controls run (single look)"

- `results/part2/controls/controls_summary.tsv` and `kill_rule.json`
- `results/part2/trials.jsonl`, with 23 lines tagged `controls`
- `results/part2/assoc/CTRL_*/controls_summary.json`
- `run_log.txt` and `runtime.json`

The full per-test tables go to `data_cache/part2/assoc/CTRL_*/controls.parquet` (gitignored).
