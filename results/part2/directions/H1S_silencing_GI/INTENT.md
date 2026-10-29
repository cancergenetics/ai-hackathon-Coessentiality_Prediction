# H1S_silencing_GI: A-side discovery run, pre-registered intent (single look; EXPLORATORY)

Committed before the run. H1S is exploratory (`standing: exploratory`, `headline: false`; amendments A1/A2) and is
never headlined. The discovery runs exactly once; any re-run is a relook and must be disclosed here.

- History: the first dry check crashed in `h1.classify` (`BLOCKED.md`). PI ruling R2 classed this an engine bug; fixed
  in `b3e8043` (NOT EVALUABLE axes skip the matched null; empty null → 'null unavailable'), with synthetic outputs
  identical before/after by SHA-256 (`../h1_fix_equivalence.json`). Erratum E1 (`c3c4422`) corrected the score_rule
  sign typo; semantics unchanged. No H1S outcome had been computed before either change.
- Code: HEAD `b3e8043fd6d37ef1ff0230cee58c26dd3d0f9374`; registry blob `57187c82b04cff9044842800a90376774af92d91`.
- Planned start: 2026-10-29, about 00:52 UTC. Command (stdout/stderr to `run_log.txt`, wrapped in `/usr/bin/time -l`):

```
export PART1_DIR="/Users/metinyazar/Desktop/autoresearch/coessentiality_protein complexes"; export PYTHONHASHSEED=0
<main repo>/.venv/bin/python scripts/44_run_direction.py --id H1S_silencing_GI
```

## What the run does

Identical machinery to H1 (`h1.assoc_inputs` → `run_association` with pairs, pooled BH q = 0.05 over this entry
only, screen BH q < 0.1, Freedman–Lane 10,000 within lineage, seed 4200 + row → `h1.classify`, K = 1000, seed 7000 +
FZ-X event index, coherence from the FZ-C A block sha256 `6cac0e53…`). 317 events, 338 registered components (hash
`736ed74bc58bc04a`), base covariates `[lineage_FE, library, ScreenType3DO, ModelType]`.

## Dry check (no dependency data; `dry_check.json`)

Entry clean, open, no placeholders; frozen inputs pass; A hash `bd28f78d612bef9f`. Synthetic N(-0.5, 1) outcome: 338
tests, 309 evaluable (the 28 leverage-1 designs of `degenerate_designs.tsv` plus label-driven complete-case shortfalls
are not evaluable and now classed NOT EVALUABLE), no crash. Runtime estimate ~2 min + ~0.9 s per screened test;
peak RSS ~4 GB.

## Excluded

No lockbox read, no hu.MAP / STRING / curated complexes. Nothing else changed after H3a/H3b/H1 outcomes.

## Outputs (commit "H1S_silencing_GI discovery run (single look)")

`summary.json`, `tests.tsv.gz`, `discoveries.tsv`, `classes.tsv`, `run_log.txt`, `runtime.json`; assoc summary; one
`trials.jsonl` line.
