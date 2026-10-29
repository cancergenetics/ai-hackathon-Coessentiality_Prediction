# H3a_lineage_atlas: A-side discovery run, pre-registered intent (single look)

Committed before the run. The discovery runs exactly once (DESIGN v3.1 §2.5–2.6; `docs/registry/H3a_lineage_atlas.yaml`,
blob `4caaa1f5d0b5dc229f33b7e5c620b32619d35ee2`). Any re-run is a relook and must be disclosed here.

- Code and registry commit: `beda09f862300249cd2961ccabf54413e5c9009f` (runner `scripts/44_run_direction.py`, added in
  that commit; engine, registry and tests unchanged since main `04ba423`). This note is committed on top of it with no
  other change.
- Planned start: 2026-10-29, about 00:25 UTC. Actual start and end go in `runtime.json`.
- Working directory: the agent worktree; `datasets/`, `data_cache/`, `external/` are symlinks to the main repo
  (never committed).
- Command (stdout and stderr to `results/part2/directions/H3a_lineage_atlas/run_log.txt`, wrapped in
  `/usr/bin/time -l` for wall time and peak RSS):

```
export PART1_DIR="/Users/metinyazar/Desktop/autoresearch/coessentiality_protein complexes"; export PYTHONHASHSEED=0
<main repo>/.venv/bin/python scripts/44_run_direction.py --id H3a_lineage_atlas
```

## What the run does (all engine calls; nothing selected)

1. Scores: `cscore.unit_scores` (S_C_mean_winsor3) for the 367 frozen coarse units (hash `53f9dbf46370d237`); no
   member removal for lineage contexts.
2. Covariates: `library`, `ScreenType3DO`, `ModelType` (lineage_FE excluded by the §2.2 rule); complete case.
3. Calibration (registry `calibration`, seed 4600): 100 pseudo-contexts, pseudo-context j matched to context j mod 23
   in registry order (seed 4600 + context index; random A sets of the same size, PatientID blocks, because every
   H3a context is lineage-type), each tested against all 367 units with the same HC3 machinery
   (`robust.pseudo_context_calibration`). Reported: pooled fraction of pseudo tests with p < 0.05 (expected 0.05)
   and BB rejections over the pseudo tests (expected ~0). Pre-declared response: if the fraction exceeds 0.075, every
   real test with raw p < 0.01 is also re-tested by permutation (`perm_screen_p=0.01`).
4. Discovery: one `assoc.run_association` call, 23 contexts × 367 units, HC3 + Satterthwaite, two-sided, evaluable
   iff ≥ 12 context and ≥ 12 reference models; BB (families = contexts, q = 0.05); screen = BB-adjusted q < 0.1
   (plus raw p < 0.01 if the response fires); exact Freedman–Lane permutation, 10,000 permutations, PatientID blocks,
   unrestricted across lineages (lineage contexts), seed 4000 + row; GPD tail only below 1e-4; BB re-run on the
   permutation p. Discovery = rejection after re-entry.
5. Not run (no engine implementation; reported-only negative controls that never enter discovery): random gene sets
   matched for size, essentiality and coherence, and dispersed never-expressed sets. Robustness (LOSO, leave-top-2-out,
   breadth) and multiverse scores are per-hit reports for a later step.

## Dry check (no dependency data)

`--dry-check` (output `dry_check.json`): real registry, manifests, contexts, covariates, map and A index; synthetic
N(-0.5, 1) outcome on the real gene columns (cache sidecar names only); trial ledger and writers disabled.
- Entry loads: committed, clean, `status: open`, no `pending`, no placeholder strings.
- `check_frozen` passes for the context table (`ba73a371…`), map (`6a062233…`) and `lof_homdel_shares.tsv`.
- A: 1,007 models, model-set hash `bd28f78d612bef9f` (= the FZ-LB A list).
- m = 8,441 evaluable tests in 23 families (registered 8,441). No context column missing.
- Synthetic pseudo-context calibration: 36,700 pseudo tests, fraction p < 0.05 = 0.046, 0 BB rejections.
- Runtime estimate: setup and calibration ~15 s; permutation 0.25 s per screened test (unrestricted design).
  Expected total 5–40 min depending on the number screened (all 8,441 would be ~35 min). Peak RSS ~3 GB.
- Synthetic check of the runner (`results/part2/directions/runner_synthetic_check.json`): H3a/H3b discoveries equal
  the reference pipeline's, no false positives, H1 class accuracy 1.0, second run refused.

## Excluded

No lockbox read, no `read_lockbox` call, no other direction run, no hu.MAP / STRING / curated complexes.

## Outputs (one commit "H3a_lineage_atlas discovery run (single look)")

`results/part2/directions/H3a_lineage_atlas/`: `summary.json`, `tests.tsv.gz`, `discoveries.tsv`, `calibration.json`,
`run_log.txt`, `runtime.json`; `results/part2/assoc/H3a_lineage_atlas/discovery_summary.json`; one line in
`results/part2/trials.jsonl`. The full parquet goes to `data_cache/` (gitignored); its SHA-256 is in the summary.
