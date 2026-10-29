# H1_within_complex_GI: A-side discovery run + D7 SMARCB1 sensitivity, pre-registered intent (single look)

Committed before the run. The discovery runs exactly once (DESIGN v3.1 §2.6; `docs/registry/H1_within_complex_GI.yaml`,
blob `ebfbaad407c552982c46586a3277602c76539803`, amendments A1/A2). Any re-run is a relook and must be disclosed here.

- Runner `scripts/44_run_direction.py` as committed in `beda09f8623…` (unchanged since). HEAD at writing: `5e1c1d1`
  (H3b outputs). Engine, registry and tests unchanged since main `04ba423`.
- Planned start: 2026-10-29, about 00:37 UTC; the sensitivity follows once the H1 outputs are committed.
- Commands (stdout/stderr to `run_log.txt` and `run_log_sensitivity.txt` in this folder, wrapped in `/usr/bin/time -l`):

```
export PART1_DIR="<PART1_DIR>"; export PYTHONHASHSEED=0
<main repo>/.venv/bin/python scripts/44_run_direction.py --id H1_within_complex_GI
# after the H1 outputs are committed:
<main repo>/.venv/bin/python scripts/44_run_direction.py --id H1_within_complex_GI --sensitivity SMARCB1_CN
```

## What the run does (engine calls only)

1. Events: the 24 genetic events of the frozen FZ-X list (`h1_events.tsv` rows whose A2 split entry is H1;
   `h1_event_split_A2.tsv`); loss columns from `data_cache/part2/fzx/h1_loss_table.parquet` (A rows only). Units: the
   28 registered components (hash `686659438f7d781c`); event covariates checked equal to the registry context items.
2. `h1.assoc_inputs` → `assoc.run_association` with `pairs` (each loss column only against its own components):
   HC3 + Satterthwaite, two-sided, ≥ 12 / ≥ 12, pooled BH q = 0.05; screen BH q < 0.1 → Freedman–Lane permutation,
   10,000, PatientID blocks within OncotreeLineage, seed 4200 + row, GPD tail only below 1e-4; BH re-run.
3. `h1.classify`: TOST (±0.1), P-vs-N contrast, member breadth, matched random-set null (K = 1000, seed 7000 + FZ-X
   event index; coherence = mean pairwise zavg_BEF on A from the FZ-C A block
   `data_cache/part2/blocks/zavg_BEF_bd28f78d612bef9f_S8c264ca513f8.npz`, sha256 `6cac0e53a9c0…`; essentiality = mean
   clipped effect over A, as in the synthetic pipeline). Class per event and axis to `classes.tsv`.
   Implementation note carried from the engine: `h1.ALIGN = +1` implements the registry's stated meaning of kappa_N
   (> 0 = dependency grows in its baseline direction), not the literal `-sign(b_j)` formula (open PI flag in h1.py).
4. Sensitivity (D7, labelled sensitivity, not a trial, no trial line): `Z0109|SMARCB1` re-classified with the
   registry's CN-type covariate set `[lineage_FE, library, ScreenType3DO, ModelType, Ploidy, WGD, TP53_altered]`
   (h1.COVARIATES["CN"]; none is a function of SMARCB1 loss) in place of the base set; same members, null and seeds;
   q / discovery flags from the primary run. Label check: n 960 → 732 complete cases, loss models 17 → 15; design
   full rank. Output `sensitivity_SMARCB1_CN.tsv` (primary row recomputed alongside for comparison).

## Dry check (no dependency data; `dry_check.json`)

Entry committed, clean, open, no `pending`/placeholders; `check_frozen` passes for the context table, map,
`lof_homdel_shares.tsv`, `h1_loss_table.parquet`, `h1_events.tsv`, `h1_event_split_A2.tsv`; A = 1,007 models, hash
`bd28f78d612bef9f`. Synthetic N(-0.5, 1) outcome with random pair scores: 24 events, 28 tests, 26 evaluable (2 below
12/12 after complete case by labels), end-to-end classify runs. Design check of all 24 events from labels: no
leverage-1 or rank-deficient design. Runtime estimate: ~20 s dry; permutation ~0.9 s per screened test; matched null
for 28 axes; total expected < 10 min. Peak RSS ~3.7 GB.

## Excluded

No lockbox read; no other direction; no hu.MAP / STRING / curated complexes. No change after the H3a/H3b outcomes.

## Outputs

Commit "H1_within_complex_GI discovery run (single look)": `summary.json`, `tests.tsv.gz`, `discoveries.tsv`,
`classes.tsv`, `run_log.txt`, `runtime.json`, `results/part2/assoc/H1_within_complex_GI/discovery_summary.json`, one
`trials.jsonl` line. Then commit "H1 SMARCB1 CN-covariate sensitivity (D7)": `sensitivity_SMARCB1_CN.tsv`,
`run_log_sensitivity.txt`.
