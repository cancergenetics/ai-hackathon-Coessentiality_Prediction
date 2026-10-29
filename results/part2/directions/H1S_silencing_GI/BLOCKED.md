# H1S_silencing_GI: not run — engine crash found in the dry check (stop and report)

**Status.** No H1S discovery run, no trial line, no intent committed. H2_submodule_divergence (next in registry
order) is also not run; it is dry-checked and ready (see below). Nothing in the engine, registry or tests was changed.

## Evidence

`--dry-check` (real registry, contexts, covariates, map, A index and frozen event list; synthetic N(-0.5, 1) outcome,
no dependency value) aborts inside `h1.classify` (`dry_check_traceback.txt`):

```
h1.py:402 classify -> matched_null(E_A, fmap, e, "N", ...)
h1.py:347 matched_null: q025=float(np.quantile(nb, 0.025))
IndexError: index -1 is out of bounds for axis 0 with size 0
```

Cause: for 28 of the 317 H1S events the HC3 design (loss column + base covariates, complete case) has a
leverage-1 observation (`stats.HC3Design` raises "leverage 1 observation"). The label-only design check
(`scripts/verify/44_h1s_design_check.py`, output `degenerate_designs.tsv`) lists them, with the model at leverage 1.
Several are lineage-restricted silencing events (e.g. the B-cell genes CD19 / CD79A / CD79B / PAX5 of Z0319); the
same few models recur at leverage 1 (ACH-000133 in 6 events, ACH-001401, ACH-002040 and ACH-003711 in 2–3 each),
i.e. a model that is the only one of its level once lineage, library, ScreenType3DO and ModelType are combined.
`assoc._design` merges indicator levels with < 2 models but cannot remove a combination that isolates one model.

- `run_association` handles this: the test becomes `evaluable = False` (reason `design: …`) and leaves m.
- `h1.classify` does not: `fit_columns` returns NaN for every one of the K random sets, `nb` is empty, and
  `np.quantile` on an empty array raises. The whole run aborts **after** `run_association` has appended its trial
  line, so a real run would leave a trial line with no outputs and need a disclosed relook.
- The same design is used on real data (the design depends on labels only; real outcome NaNs can only add failures),
  so the crash is deterministic, not a dry-check artefact.
- H1 is unaffected: the same check finds 0 of 24 H1 events with a degenerate design, and the H1 run completed.

## Decision needed (PI)

Options (any of them is an engine change, so it needs a ruling and a re-review before the single H1S look):
1. `matched_null`: return `status = "null unavailable"` when no null fit succeeds (class → INCONCLUSIVE, which the
   registered matched-null rule already provides for), and/or `classify`: skip the axis when its counted test is not
   evaluable (class "not evaluable"), consistent with `run_association`.
2. Declare the 28 events not evaluable (labels only, before any outcome) by amendment.

## H2 readiness (dry check, synthetic outcome, real committed H3a/H3b parent tables)

Entry loads clean; A hash `bd28f78d612bef9f`; 210 eligible (context, parent) families, m = 452 sub-module tests;
permutation ~0.9 s per screened test; dry runtime 20 s. Not run because registry order puts H1S first.
