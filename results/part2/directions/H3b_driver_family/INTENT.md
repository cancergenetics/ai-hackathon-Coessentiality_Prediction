# H3b_driver_family: A-side discovery run, pre-registered intent (single look)

Committed before the run. The discovery runs exactly once (DESIGN v3.1 §2.5–2.6; `docs/registry/H3b_driver_family.yaml`,
blob `9cf55f21099f108ade10211b3d54ec4e8d0fbaac`). Any re-run is a relook and must be disclosed here.

- Code and registry commit: runner `scripts/44_run_direction.py` as committed in `beda09f8623…` (unchanged since);
  HEAD at writing `b0dbd83b88ad99b36f59d7b94efe97ca1e099d67` (H3a outputs). Engine, registry and tests unchanged
  since main `04ba423`. This note is committed on top with no other change.
- Planned start: 2026-10-29, about 00:32 UTC. Actual start and end go in `runtime.json`.
- Command (stdout/stderr to `results/part2/directions/H3b_driver_family/run_log.txt`, wrapped in `/usr/bin/time -l`):

```
export PART1_DIR="<PART1_DIR>"; export PYTHONHASHSEED=0
<main repo>/.venv/bin/python scripts/44_run_direction.py --id H3b_driver_family
```

## What the run does (engine calls only)

1. Scores: `cscore.unit_scores` for the 367 frozen coarse units (hash `53f9dbf46370d237`) with each context's
   registered removals (`cscore.removal_for_context`; LoF ±10 Mb window from the frozen `lof_homdel_shares.tsv`,
   all < 0.5, so none applies).
2. Covariates: each context's registered list (`covariates.per_context`; §2.2 rule applied); complete case.
3. Calibration (seed 4700): 100 pseudo-contexts, pseudo-context j matched to primary context j mod 20 (registry order;
   seed 4700 + context index; lineage-stratified, PatientID-blocked relabelling), each tested against all 367 units
   with the same HC3 machinery and covariates. Reported: pooled fraction p < 0.05 (expected 0.05), BB rejections over
   the pseudo tests. Pre-declared response: if > 0.075, every real test with raw p < 0.01 is also permutation re-tested.
4. Discovery: one `assoc.run_association` call, the 20 primary contexts × 367 units (sensitivity-role contexts are not
   run), HC3 + Satterthwaite, two-sided, ≥ 12 / ≥ 12; BB (families = contexts, q = 0.05); screen BB q < 0.1;
   Freedman–Lane permutation, 10,000, PatientID blocks within OncotreeLineage, seed 4100 + row; GPD tail only below
   1e-4; BB re-run. Discovery = rejection after re-entry.
5. Not run (no engine implementation; reported-only negative controls outside discovery): matched random gene sets,
   dispersed never-expressed sets, the CN segment-matched null. Robustness and multiverse: later per-hit step.

## Dry check (no dependency data; `dry_check.json`)

Entry committed, clean, open, no `pending`/placeholders; `check_frozen` passes (context table `ba73a371…`, map
`6a062233…`, `lof_homdel_shares.tsv` `1eac6ce9…`); A = 1,007 models, hash `bd28f78d612bef9f`; m = 7,340 in 20 families
(registered 7,340); no context column missing. Synthetic calibration: fraction 0.041, 0 BB rejections.
Runtime estimate: ~25 s setup and calibration + ~1.0 s per screened test (within-lineage blocks); 5–40 min expected.
Peak RSS ~3 GB.

## Excluded

No lockbox read, no other direction, no hu.MAP / STRING / curated complexes. H3a outcomes (committed `b0dbd83`) do
not change anything in this run: code, registry, covariates and thresholds are as pre-registered.

## Outputs (one commit "H3b_driver_family discovery run (single look)")

As for H3a: `summary.json`, `tests.tsv.gz`, `discoveries.tsv`, `calibration.json`, `run_log.txt`, `runtime.json`;
`results/part2/assoc/H3b_driver_family/discovery_summary.json`; one `trials.jsonl` line.
