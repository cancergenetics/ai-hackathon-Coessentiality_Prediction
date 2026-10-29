# Screen-quality sensitivity of H3a / H3b: intent (pre-declared in D8)

Committed **before** anything below is computed on dependency data. The sensitivity was pre-declared in
RESEARCH_TREE.md D8 (commit `84dfa98`, spec fixed in `3f5e8fd`). It is a **sensitivity analysis, not a trial**:
- no line is written to `results/part2/trials.jsonl` (checked by SHA-256 before and after);
- no registered output, registry entry or threshold changes;
- the primary H3a / H3b discoveries stand.

Every output carries the label `sensitivity_screen_quality`.

## What changes

The only change is two added covariates, ScreenNNMD and CasActivity, taken per model from
`datasets/DepMAP/CRISPRConfounders.csv`. Growth pattern is excluded as collinear with lineage (D8). Each context's
covariate list is its registered list, after the §2.2 covariate-of-context rule, followed by `NNMD`, `CasActivity`:
- H3a: library, ScreenType3DO, ModelType, NNMD, CasActivity.
- H3b: each context's registered list (e.g. KRAS: lineage_FE, library, ScreenType3DO, ModelType, log_burden,
  MSIScore), then NNMD, CasActivity.

The engine already assembles these registry names (`covariates.assemble`; `NNMD` → `ScreenNNMD`), so the engine is
not changed.

**Keying and multiple screens (choice fixed here).**
- `CRISPRConfounders.csv` is keyed by **ModelID**, one row per model: 2,554 unique ModelIDs, all 1,007 A models
  present.
- `CRISPRGeneEffect.csv` is also one row per ModelID. DepMap has already combined a model's screens into one
  model-level row. The confounder file follows the same rule: library is multi-hot (171 A models carry 2–4 libraries),
  and ScreenNNMD / CasActivity hold one model-level value.
- So **we use the model-level row as published, joined on ModelID**, exactly as the registered covariates (library,
  ScreenType3DO) and the gene-effect matrix are. We do no screen-level aggregation of our own. No screen-level QC file
  is in `datasets/`, so the published model-level values cannot be re-derived here.
- This is the existing source table `covariates.source_table`, unchanged.

**Missing values: complete case.** This is the engine's design-builder rule. Covariates only, no outcomes:
- ScreenNNMD is complete on A.
- CasActivity is missing for 119 of 1,007 A models. That includes all 93 KY-only (Sanger) models, 18 of 145
  HumagneCD-only, 6 Avana-only and 2 others.
- In H3a, all 119 models drop out of every test, since the registered H3a covariates are complete.
- In H3b, a model drops out only if the context or a registered covariate did not already exclude it. That is 30–119
  models per context (`dropout_per_context.tsv`).
- Drop-outs per lineage (of n_A): CNS/Brain 19/99, Bowel 17/64, Head and Neck 16/58, Lung 8/93,
  Esophagus/Stomach 8/72, PNS 7/35, Bone 6/44, Pancreas 6/53, Breast 5/50, Ovary 5/46, Lymphoid 4/73, Prostate 4/14,
  Uterus 4/28, Skin 3/51, Myeloid 3/33, Kidney 2/26, Liver 1/19, Other 1/1; 0 for Biliary Tract, Bladder, Cervix,
  Pleura, Soft Tissue and the small lineages.

**Consequences for evaluability (registered ≥ 12 / ≥ 12 rule kept).**
- H3a: Prostate falls to 10 context models and becomes not evaluable. That gives 22 families and m = 8,074, against
  23 and 8,441. Prostate had 0 primary discoveries.
- H3b: LoF:KMT2D falls to 11 and becomes not evaluable. That gives 19 families and m = 6,973, against 20 and 7,340.
  The rare-context route is a confirmation rule and does not relax evaluability.

These counts come from `--dry-check`, which uses real covariates and contexts and a synthetic outcome. No dependency
value was read for them.

## Everything else as registered

Units, families, scores, removals, HC3 + Satterthwaite, BB at q = 0.05 (families = contexts), screen q < 0.1,
Freedman–Lane permutation (B = 10,000, PatientID blocks; unrestricted for lineage contexts, within lineage
otherwise), GPD tail only below 1e-4, BB re-run on the permutation p, and seeds:
- H3a: permutation 4000 + row, pseudo-contexts 4600 + context index.
- H3b: permutation 4100 + row, pseudo-contexts 4700 + context index.

The pseudo-context calibration (100 pseudo-contexts) is re-run with the sensitivity covariates. Its pre-declared
response (if the fraction exceeds 0.075, raw p < 0.01 is also permutation re-tested) applies mechanically.

The engine's run record is written to `<ID>/sensitivity_screen_quality/engine_record.jsonl`, not to the trial ledger.
This uses a `Project` whose `trials` path points there.

Code: `scripts/48_sensitivity_screen.py`, which reuses the runner `scripts/44_run_direction.py` and the engine
unchanged. Test: `tests/test_sensitivity_screen.py`, on synthetic data. It checks:
- the covariates are the registered ones plus the two added;
- dropout counts are correct;
- the trial ledger is untouched;
- the run is refused a second time;
- with constant added covariates the result equals the registered run.

## Commands (one run each, H3a then H3b)

Stdout and stderr go to `<ID>/sensitivity_screen_quality/run_log.txt`, wrapped in `/usr/bin/time -l`.

```
export PART1_DIR="/Users/metinyazar/Desktop/autoresearch/coessentiality_protein complexes"; export PYTHONHASHSEED=0
<main repo>/.venv/bin/python scripts/48_sensitivity_screen.py --id H3a_lineage_atlas
<main repo>/.venv/bin/python scripts/48_sensitivity_screen.py --id H3b_driver_family
```

## Reported (against the committed primary tables)

The primary tables are read from `data_cache/part2/assoc/<ID>/discovery.parquet`, with SHA-256 checked against the
committed `discovery_summary.json`.

- **Discoveries per lineage (H3a) and per context (H3b):** primary, sensitivity, surviving (in both), lost, new, and
  the surviving fraction.
- **Pseudo-context rate:** the fraction p < 0.05 and BB rejections, before (committed `calibration.json`) and after,
  with the change.
- **Lymphoid "scale effect":** the Spearman ρ between unit essentiality (mean S_C over all A models) and the
  lineage:Lymphoid β over its 367 evaluable units, before and after. The same ρ is given for every H3a lineage.
  Essentiality is the same unit property in both.
- **Complete-case dropout:** per lineage and per context.

Outputs go to `results/part2/directions/<ID>/sensitivity_screen_quality/`:
- `summary.json`, `tests.tsv.gz`, `discoveries.tsv`, `comparison_per_context.tsv`;
- `dropout_per_context.tsv`, `dropout_per_lineage.tsv`, `calibration.json`, `engine_record.jsonl`, `run_log.txt`;
- H3a only: `scale_effect_rho.tsv`.

The engine also writes `results/part2/assoc/<ID>/sensitivity_screen_quality_summary.json` and the full parquet in
`data_cache/` (gitignored).

## Excluded

No lockbox read, no other direction, no hu.MAP / STRING / curated complexes. H2, which is parented on the primary
H3a / H3b discoveries, is not re-run.
