# H2_submodule_divergence: A-side discovery run, pre-registered intent (single look; EXPLORATORY)

Committed before the run. H2 is **exploratory**: its kill-rule controls PC21/PC22 are not evaluable (controls run,
`kill_rule.json`: 0 evaluable), so the rule is n/a (DESIGN §2.6). The discovery runs exactly once; any re-run is a
relook and must be disclosed here. Registry `docs/registry/H2_submodule_divergence.yaml`, blob
`484da0f11375db0cd17f7230782be3c082a75bd4`.

- Code: HEAD `06c8cd1` (H1S outputs; runner unchanged since `beda09f`, h1 fix `b3e8043` does not touch H2).
- Planned start: 2026-10-29, about 01:02 UTC. Command (stdout/stderr to `run_log.txt`, `/usr/bin/time -l`):

```
export PART1_DIR="<PART1_DIR>"; export PYTHONHASHSEED=0
<main repo>/.venv/bin/python scripts/44_run_direction.py --id H2_submodule_divergence
```

## What the run does

`assoc.parent_selection_from` the committed H3a (`b0dbd83`) and H3b (`5e1c1d1`) discovery tables (each read only if
its SHA-256 matches the committed summary; families = (context, discovered parent with ≥ 2 sub-modules), level
q·R_x/m_x) → `h2.h2_frames` (sub-module residualised on its sibling complement, A-fitted OLS, context removals) →
one `run_association` with `stagewise_bb` (q = 0.05), 43 contexts, the 109 registered sub-modules (hash
`141a70b7133ed8a4`; only eligible pairs tested), covariates per context + NNMD + CasActivity, ≥ 12 / ≥ 12; screen
q < 0.1 → Freedman–Lane 10,000 (within lineage; unrestricted for lineage contexts), seed 4300 + row; stage-wise re-run.

Calibration (random member splits, seed 4800; H3 pseudo-contexts through the stage-wise machinery) is not in the
engine at this commit; it is implemented and run afterwards under PI ruling R4 (D8).

## Dry check (`dry_check.json`; synthetic N(-0.5, 1) outcome, real committed parent tables)

Entry clean; frozen inputs pass; A hash `bd28f78d612bef9f`; 210 eligible (context, parent) families; m = 452
evaluable; ~0.9 s per screened test; dry runtime 20 s; expected total < 10 min.

## Outputs (commit "H2_submodule_divergence discovery run (single look)")

`summary.json`, `tests.tsv.gz`, `discoveries.tsv`, `run_log.txt`, `runtime.json`; assoc summary; one trials line.
