# R6 corrected refit: intent (NON-CONFIRMATORY relook)

This note is committed before any relook. It was written at 2026-10-29T10:15:59Z (UTC). The authorising PI note R6 is
quoted verbatim in `README.md`.

| | |
|---|---|
| Fix | `scripts/51_lockbox_confirm.l_covariates`. The L covariate design is assembled over the L rows read, using `covariates.per_context` and the conventions of `assoc._design` (FE levels with fewer than 2 models are merged into the reference), exactly as the A side assembled it over A rows. `confirm.confirm_entry` accepts the builder as a callable of the models read. |
| Unit test | `tests/test_part2_lockbox_driver.py::test_r6_l_covariates_assembled_on_L_rows` reproduces the leverage-1 failure and passes with the fix. |
| Access | `assoc.resume_lockbox_read(..., relook="R6")` reads the persisted reads only, never the loader, and writes no ledger line. It requires the committed R6 note and this note. |

## Commands (in order, from the worktree root)

```
export PART1_DIR="<PART1_DIR>"; export PYTHONHASHSEED=0
PY=<PROJECT_ROOT>/.venv/bin/python
# (3) reproduction check: every committed successful L fit, recomputed with the fixed builder; same seeds
caffeinate -i $PY scripts/51_lockbox_confirm.py --relook R6 --mode check --out <scratch>/r6_check.tsv
# (4) only if the check is identical: the failed fits
caffeinate -i $PY scripts/51_lockbox_confirm.py --relook R6 --mode refit \
    --out results/part2/lockbox_confirm/corrected_refit_noncomfirmatory.tsv
```

R6 (3) says: if any recomputed estimate or p-value differs from the committed one, STOP and report.

**Disclosed before the check (labels only, no outcome value).** The check can only pass for designs whose column span
is unchanged by the fix.
- I compared the two designs (all-model assembly subset to L, versus L-only assembly) on the context table and
  covariate sources for every (entry, context) of the confirmation tables.
- Same span: 2,524 hits (H3a 1,145, H2 101, BRIDGE-L 1,206, controls PC5/PC17/PC20, H1 5, H1S 58).
- Different span: 602 hits. 510 of these fitted successfully under the defect, but with one extra design column: a
  separate group for single-model lineages, because the all-model reference lineage is absent from the L subset.
  - H3b: all 122 successful hits;
  - H1: 17, including ARID1A P and TP53 N;
  - H1S: 220;
  - H2: 7;
  - BRIDGE-L GLS: 120;
  - controls: PC1–PC4, PC6–PC8, PC11, PC14, PC19a/b.
- For these 510 fits the reproduction check is therefore expected to differ. It will be run as R6 specifies. If it
  differs, the run stops there and the failed-fit refit (4) is not run without a further PI ruling.

## Failed fits (eligible hits on the lockbox route with no committed L estimate)

Of 53 failed fits, 50 are caused by the driver defect. The other 3 also fail with the fixed builder: CTRL_PC21 ×2
(heme_vs_breast, a single-model library level among 49 L models) and H1S Z0035|HNF1B|N (8 loss models). Those 3
remain NOT EVALUABLE: a genuine single-model covariate level in L.

| Entry | Context | Unit | n_context_L | Cause |
|---|---|---|---|---|
| BRIDGE-L_confirmation_rate | 9p21_loss | G0061, G0221 | 35 | driver defect (R6) |
| BRIDGE-L_confirmation_rate | WGD | G0005, G0024, G0040, G0076, G0100, G0272, G0315 | 173 | driver defect (R6) |
| CTRL_PC9 | 9p21_loss | PRMT5 (10419), WDR77 (79084), Z0228 | 35 | driver defect (R6) |
| CTRL_PC10 | FOCAD_del_or_MSI_H | HBS1L (10767), PELO (53918), Z0286 | 19 | driver defect (R6) |
| CTRL_PC12 | VPS4B_loss | VPS4A (27183), Z0115 | 85 | driver defect (R6) |
| CTRL_PC13 | female_or_LoY | DDX3X (1654), Z0137 | 212 | driver defect (R6) |
| CTRL_PC18 | WGD | KIF18A (81930) | 173 | driver defect (R6) |
| CTRL_PC21 | heme_vs_breast | CCND3 minus CCND1, CDK6 minus CDK4 | 33 | genuine (leverage 1 also with L-only covariates) |
| H1_within_complex_GI | H1loss:FOCAD | Z0286\|FOCAD\|N | 10 | driver defect (R6) |
| H1S_silencing_GI | H1loss:HNF1B | Z0035\|HNF1B\|N | 8 | genuine (leverage 1 also with L-only covariates) |
| H2_submodule_divergence | WGD | Z0026.2 | 160 | driver defect (R6) |
| H3b_driver_family | 9p21_loss | Z0069, Z0228 | 35 | driver defect (R6) |
| H3b_driver_family | WGD | Z0026, Z0035, Z0056, Z0065, Z0080, Z0081, Z0104, Z0141, Z0150, Z0163, Z0166, Z0187, Z0206, Z0235, Z0248, Z0257, Z0265, Z0276, Z0310, Z0319, Z0334, Z0341, Z0342, Z0354, Z0362, Z0365 | 173 | driver defect (R6) |

The output files are labelled `corrected_refit_noncomfirmatory` (the spelling of the R6 instruction). They never
overwrite the primary outputs.

## Addendum (R7): the scope is extended to all 560 affected fits

Written at 2026-10-29T10:41:42Z (UTC), before the corrected relook. PI ruling R7 is quoted verbatim in `README.md`. The reproduction
check (`68f0d3f`) found 510 changed fits, exactly the predicted set.

**Scope.** The corrected relook recomputes all 560 affected fits with the fixed L covariate design:
- the 510 committed fits whose design span the fix changes;
- the 50 lockbox-route eligible fits that failed through the defect.

The 3 genuine failures (CTRL_PC21 ×2, H1S Z0035|HNF1B|N) are refitted too, and are expected to stay NOT EVALUABLE.

**How it runs.** Through the R6 access path only: `assoc.resume_lockbox_read(..., relook="R6")`. The loader is never
called and no ledger line is written. Seeds, hit order, blocks and strata are those of the single read.

```
caffeinate -i $PY scripts/51_lockbox_confirm.py --relook R6 --mode corrected \
    --out results/part2/lockbox_confirm/corrected_refit_nonconfirmatory.tsv
$PY scripts/55_lockbox_r7_report.py      # primary relabel (R7) + primary-vs-corrected comparison; reads committed tables only
```

**Outputs** (separate files; primary outputs are never overwritten):
- `corrected_refit_nonconfirmatory.tsv`;
- `primary_relabelled_R7.tsv.gz` and `rates_primary_R7.tsv`;
- `comparison_primary_vs_corrected_R7.tsv`, `defect_sensitive_calls_R7.tsv` and `h1_events_primary_vs_corrected_R7.tsv`;
- `bridge_primary_vs_corrected_R7.json`.

**Defect-sensitive.** An eligible call whose confirmed status differs between the primary result (NOT EVALUABLE for
the failed fits) and the corrected relook. Such calls cannot support a claim.
