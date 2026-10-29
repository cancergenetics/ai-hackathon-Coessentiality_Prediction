# Intent: registered per-hit robustness and the final strongest-findings list

Committed before any robustness value is computed. Code: `src/part2/strongest.py`, `scripts/56_strongest_findings.py`.

## Why

The registered strongest-findings rule (DESIGN v3.1 §2.5; `success_criterion` of H3a_lineage_atlas,
H3b_driver_family, H1_within_complex_GI) requires each finding to be LOSO-robust. The per-hit LOSO and leave-top-2-out
checks (`multiverse.robustness`) were deferred at the discovery runs and never computed, so the manuscript labels its
strongest set "provisional" (manuscript/main.md, limitation 15). This is a registry-to-code gap (LESSONS rule 15).

**Timing disclosure.** These checks run after the lockbox read (5bc9444). Their rule was registered before it. They
use A models only, and no lockbox value is read. The L summary statistics shown next to each hit are copied from the
committed confirmation tables.

## The rule, applied mechanically

A finding is **strongest** if it meets all of the following:

- discovery q < 0.05;
- the PRIMARY lockbox call is "confirmed" or "rare-context route passed" (labelled);
- it is not defect-sensitive (R7);
- it is not a locus cluster;
- it is **LOSO-robust**.

For H1, two further conditions apply:

- the call is a non-neutral class call (eligible);
- the call is headline-eligible under E2.

**LOSO-robust.** The full A fit is refitted with the discovery run's score, removals, covariates and alternative. Every
evaluable leave-one-group-out fold must have the same sign as the full fit and p < 0.05. A fold that leaves fewer than
2 context or 2 reference models, or has a singular design, is skipped and counted, not failed. A hit with 0 evaluable
folds is not robust. The engine functions are `robust.loso` and `robust.fit_one`, committed 201d8bf before FZ-X. The
design conventions (`assoc._design`) are applied on each fold's rows.

**LOSO unit, as registered:**

| Direction | LOSO unit | Registry source |
|---|---|---|
| H3a | leave-one-OncotreePrimaryDisease-out | H3a `multiverse.robustness` |
| H3b | leave-one-OncotreeLineage-out | H3b: "leave-one-lineage-out (LOSO)" |
| H1 | leave-one-OncotreeLineage-out | see below |

The H1 entry writes "LOSO" without expanding it. Every registry entry that does expand it uses
"leave-one-lineage-out (LOSO)": H3b, H2, all CTRL entries and the H1S A5 amendment. The H1 rare route was already run
with the lineage unit. Round-1 reviews once used LOSO to mean "leave-one-subunit-out", but no registry entry adopted
that meaning. This reading is disclosed.

**Reported per hit, not used in the rule:**

- leave-top-2-out: the 2 context models most extreme in the effect's direction are dropped, then same sign and
  p < 0.05 are checked (`robust.leave_top2_out`);
- breadth (`robust.breadth`) for H3. For H1, breadth comes from the engine's `classes.tsv`.

The registry's `use` field reads "reported per hit, never selected", and the success_criterion names LOSO only.

**Ranking ("A statistic × tier").** The registry gives no numeric tier weights. Findings are therefore ranked by tier
class first, then by |t_A| in descending order. The tier classes are T1, then T2, then T3 (T3-neg and T3-nc together),
which puts "complex" before "co-dependency module" (DESIGN §2.3).

**Scope.** H1S and H2 are exploratory (H1S by registry standing; H2 because its kill rule is n/a) and are not part of
the strongest-findings list.

## Candidate list

The list is in `results/part2/robustness/candidates.tsv`, computed from committed files only by
`scripts/56_strongest_findings.py --candidates`. Its SHA-256 is
`033d83b5a28a55eed21a51c9dfbf9f9bb9da4aeb04e68bd42a61165da8de0001`. Robustness is computed for every row with
`run = True`:

| Direction | Discoveries | Candidates (all non-LOSO criteria) | Rows computed |
|---|---|---|---|
| H3a | 1,145 | 524 | 524 |
| H3b | 161 | 46 (41 lockbox + 5 rare route, EWSR1 fusion) | 46 |
| H1 | 5 eligible calls | 2 (ARID1A P, EP300 P) | 5 (all eligible calls) |

**H3a candidates per context:**

| Context | Candidates |
|---|---|
| Lymphoid | 119 |
| CNS/Brain | 58 |
| Head and Neck | 54 |
| Peripheral Nervous System | 44 |
| melanoma | 40 |
| Myeloid | 37 |
| Bowel | 29 |
| Bone | 26 |
| Breast | 23 |
| Skin | 22 |
| Kidney | 18 |
| Pancreas | 17 |
| Lung | 16 |
| Esophagus/Stomach | 6 |
| Soft Tissue | 5 |
| Ovary/Fallopian Tube | 3 |
| Bladder/Urinary Tract | 2 |
| Liver | 2 |
| Uterus | 2 |
| Pleura | 1 |

**H3b candidates:**

| Context | Units |
|---|---|
| TP53_altered | Z0004 Z0010 Z0026 Z0041 Z0062 Z0066 Z0075 Z0115 Z0126 Z0140 Z0145 Z0146 Z0152 Z0157 Z0158 Z0166 Z0184 Z0185 Z0188 Z0222 Z0247 Z0253 Z0295 Z0303 Z0314 Z0337 Z0355 Z0365 |
| MSI_H | Z0157 Z0274 Z0286 Z0316 Z0334 Z0339 |
| fusion:EWSR1::FLI1\|ERG (rare route) | Z0016 Z0042 Z0059 Z0160 Z0343 |
| hotspot:BRAF | Z0018 Z0051 |
| hotspot:KRAS | Z0016 Z0051 |
| hotspot:NRAS | Z0132 |
| LoF:PTEN | Z0316 |
| LoF:RB1 | Z0069 |

**H1 eligible calls:**

| Call | Status |
|---|---|
| Z0109\|ARID1A\|P | confirmed; candidate |
| Z0186\|EP300\|P | rare route passed; candidate |
| Z0286\|FOCAD\|N | defect-sensitive; not a candidate |
| Z0130\|PTEN\|N | not confirmed; not a candidate |
| Z0158\|TP53\|N | confirmed, but E2 within margin and not headline-eligible; not a candidate |

**Reconciliation with the manuscript's "strongest specific" set (261 H3a and 29 H3b).** That set adds two D8/D11
reporting filters to the registered rule: "beats all matched sets" and "not a control readout". It is also limited to
lockbox confirmations. With those filters, `candidates.tsv` reproduces 261 and 29 exactly (column
`writer_specific_set`). The registered rule does not contain those filters, so the final list uses the full set of 524
and 46. The `writer_specific_set` column is carried through to the outputs, so the narrower set can still be read off.

## Checks before results are accepted

- `--run` refuses to start unless this note and `candidates.tsv` are committed and unchanged.
- `--run` recomputes the candidate list and requires it to be byte-identical to the committed file.
- For every hit, the full A refit must reproduce the committed discovery β (1e-6) and t (1e-4).

## Outputs

All in `results/part2/robustness/`:

- `robustness_per_hit.tsv`
- `strongest_findings.tsv`
- `summary.json`
- `README.md`

No discovery, confirmation, registry entry or ledger is changed, and this is not a trial.
