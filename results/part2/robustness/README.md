# results/part2/robustness

This folder holds the registered per-hit robustness checks (LOSO, leave-top-2-out, breadth) and the final
strongest-findings list. The rule is in DESIGN v3.1 §2.5 and the `success_criterion` of H3a, H3b and H1. It is
implemented in `src/part2/strongest.py` and run by `scripts/56_strongest_findings.py`. The intent note
(`INTENT.md`) and the candidate list (`candidates.tsv`) were committed in fd03b03, before any robustness value was
computed.

## Disclosure: computed after the lockbox read

**These checks were computed after the single lockbox read (5bc9444).** They were registered before it: the
strongest-findings rule, LOSO and leave-top-2-out are in DESIGN §2.4–2.5 and in each entry's `multiverse.robustness`.
They were not run at the discovery runs ("deferred for a later step"), which is a registry-to-code gap (LESSONS rule
15). Running them now does not change any discovery, confirmation, registry entry or ledger, and it is not a trial.

- **Data used:** the checks use A models only. No lockbox value, persisted read or loader was opened.
- **Columns copied from committed tables:** the L columns (`beta_L`, CI, `p_L`, `n_context_L`) come from
  `results/part2/lockbox_confirm/<ID>/confirmation.tsv`.
- **Remaining forking path:** the rule, the LOSO units, the candidate list and the ranking were fixed in `INTENT.md`
  before computation. What remains is a forking path in *reading* ambiguous registry text, not in choosing analyses.
  The two ambiguous points are the H1 LOSO unit and the meaning of "A statistic × tier", both described below.

## Rule (applied mechanically)

A finding is strongest if all of these hold:

- discovery q < 0.05;
- the PRIMARY lockbox call is "confirmed" or "rare-context route passed" (labelled in `route`);
- it is not defect-sensitive (R7);
- it is not a locus cluster;
- it is **LOSO-robust**: every evaluable fold has the same sign as the full A fit and p < 0.05 under the discovery
  test's alternative. A fold that leaves fewer than 2 context or 2 reference models is skipped and counted, not
  failed.

For H1, the call must also be non-neutral (eligible) and headline-eligible under E2.

The LOSO unit for each direction:

| Direction | LOSO unit |
|---|---|
| H3a | OncotreePrimaryDisease (registry text) |
| H3b | OncotreeLineage ("leave-one-lineage-out (LOSO)") |
| H1 | OncotreeLineage |

The H1 entry writes "LOSO" without expanding it. Every entry that expands the term reads it as leave-one-lineage-out,
and the H1 rare route used the same unit.

- **Leave-top-2-out and breadth:** reported for every hit, not used to select hits. The registry `use` field says
  "reported per hit, never selected", and the success criterion names LOSO only.
- **Ranking:** "A statistic × tier" has no registered numeric tier weights. Hits are ranked by tier class first (T1,
  then T2, then T3-neg/T3-nc), then by |t_A| in descending order.
- **Reproduction check:** the full A refit reproduced every committed discovery β (to 1e-6) and t (to 1e-4) for all
  575 hits.

## Results

| Direction | Candidates | LOSO-robust (strongest) | Units | Leave-top-2 robust among strongest | Manuscript "specific" subset surviving |
|---|---|---|---|---|---|
| H3a | 524 | **207** | 137 | 207 | 107 / 261 |
| H3b | 46 (41 lockbox, 5 rare route) | **45** | 41 | 45 | 28 / 29 |
| H1 | 2 (5 eligible calls computed) | **2** (ARID1A→ARID1B P; EP300→CREBBP P, rare route) | 2 | 2 | – |

The strongest findings by tier are:

| Direction | T1 | T2 | T3-neg |
|---|---|---|---|
| H3a | 124 | 56 | 27 |
| H3b | 27 | 9 | 9 |
| H1 | 1 | – | 1 |

The one H3b failure is MSI_H Z0334, which failed 1 of 30 lineage folds.

**H3a LOSO behaviour.**

- Of the 317 H3a failures, 308 fail exactly one of the 81–82 primary-disease folds.
- Failures cluster in lineages dominated by one primary disease. Examples:
  - Bowel: 0 of 29 robust. Colorectal adenocarcinoma is 60 of the 64 Bowel A models.
  - Head and Neck: 1 of 54 robust. HNSCC is 56 of 58.
  - Skin: 1 of 22 robust. Melanoma is 44 of 51.
  - Peripheral Nervous System: 5 of 44 robust. Neuroblastoma is 31 of 35.
- In these lineages, leaving out the dominant disease leaves only a handful of context models.
- When fewer than 2 context models remain, the fold is skipped instead of failed, as the engine rule requires. This is
  why Pancreas passes 17 of 17: 52 of its 53 models are one disease.
- The registered rule therefore tests "not carried by one primary disease". It is strict for lineages with a small
  non-dominant remainder and vacuous when the remainder has fewer than 2 models. This is a property of the registered
  rule, reported here, not altered.
- Lymphoid (92 of 119 robust) and Myeloid (21 of 37) have several sizeable diseases and pass more often.

Three caveats from the manuscript carry over to this list:

- **Lymphoid scale effect.** 55 of the 92 lymphoid strongest findings have positive β (lower lymphoid dependency on
  the unit). The manuscript's lymphoid scale-effect caveat (§4.5) applies to them.
- **Control readouts.** 2 rows are registered control readouts (`control` column), which the manuscript does not count
  as separate findings.
- **Control units in other contexts.** 21 rows involve control units read out in another context.

## Files

| File | Contents |
|---|---|
| `INTENT.md` | Rule, LOSO units, candidate list (committed before computation) |
| `candidates.tsv` | Every H3a/H3b discovery and every H1 eligible call, with the non-robustness criteria; `writer_specific_set` reproduces the manuscript's 261/29 |
| `robustness_per_hit.tsv` | Per hit: `loso_folds`, `loso_pass`, `loso_skipped`, `loso_robust`, `top2_beta`, `top2_p`, `top2_same_sign`, `top2_robust`, `breadth` (H3), `strongest`, `rank` |
| `strongest_findings.tsv` | The 254 strongest findings, ranked; A β/t, L β with 95% CI (committed), LOSO and leave-top-2 outcomes, members |
| `summary.json` | Counts, runtime (118 s wall, 2.5 GB peak RSS) |

## Reproduce

Set `PART1_DIR` and `PYTHONHASHSEED=0`, then run:

```
.venv/bin/python scripts/56_strongest_findings.py --candidates   # committed files only
.venv/bin/python scripts/56_strongest_findings.py --run          # needs INTENT.md + candidates.tsv committed
```

## Not covered

- H1S_silencing_GI (exploratory by registry standing) and H2_submodule_divergence (exploratory: its kill rule is n/a)
  are outside the strongest-findings rule.
- The other multiverse items (alternative scores, covariate sensitivity) remain not run, as before.
