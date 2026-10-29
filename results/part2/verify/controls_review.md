# Controls run: independent review

**Verdict: ACCEPT-WITH-NOTES.** The process is clean (one look, with the intent note committed first). The
reviewer's own code reproduces all 15 tests in the six controls it recomputed, and the kill-rule tally applies
DESIGN §2.6 and the registry as written. None of the notes changes a result or a verdict.

- Reviewed: `4b5f800` "Controls run (single look)" on top of `e82b657` (intent) and `996a9a3` (code and registry).
- Script: `scripts/verify/controls_independent.py`. Output: `results/part2/verify/controls_independent.tsv`.
- The script does not import anything from `src/part2`. Its inputs are the registry text, the frozen map, the context
  table, `A_models.txt`, `CRISPRGeneEffect.csv` (A rows only), the raw DESIGN §2.5 covariate sources
  (`CRISPRConfounders.csv`, `Model.csv`) and the Ensembl-116 paralogue table.

## 1. Process

| Check | Result |
|---|---|
| Intent before results | Passes. `e82b657` (23:55:37 UTC) changes only `INTENT.md` on top of `996a9a3`. The run starts at 23:55:55, and `4b5f800` changes only `results/`, with no code or registry change. |
| Command | Passes. `runtime.json` records the INTENT command: `34_run_controls.py --map results/part2/fzx/zavg_engine_membership.tsv --context-table …/context_table.parquet`, with PYTHONHASHSEED=0 and `/usr/bin/time -l`. `run_log.txt` contains only stdout and the time output, so the command is recorded in `runtime.json`, abbreviated (`PART1_DIR=...`, `.venv/bin/python`). Minor. |
| trials.jsonl | Passes. It has 23 lines with 23 unique `analysis_id`, equal to the set of `CTRL_*.yaml`. All lines have tag `controls`, `head` = `e82b657`, `n_perm` = 0, and input sha256 equal to the registered manifests (`6a062233…`, `ba73a371…`). Every `registry_blob` equals the blob at `996a9a3`. |
| Lockbox firewall | Passes. `n_models` = 1007 on every line. `models_hash` `bd28f78d612bef9f` = sha256 of the sorted A list joined by newlines, which the reviewer recomputed. A ∩ L = ∅, and every `n` in the summary is ≤ 1007. `lockbox_ledger.jsonl` does not exist. |
| Single run | Passes. There are 23 cache parquets with mtimes 23:56:24–23:56:38, inside the one run window. Across all refs, only these two commits touch `trials.jsonl` and `controls/`. `runtime.json` records `n_runs` 1 and `relook` false. |

## 2. Independent recomputation

The reviewer recomputed six controls by following each registry entry:
- **Score.** The mean of the frozen `member` genes, each clipped to [−3, 3], after the registered removals.
- **PC6.** The paralog component is the Z0109 members that are Ensembl-116 paralogs of ARID1A, which is ARID1B
  alone. ARID1A is removed.
- **Model and test.** OLS on the registered covariates with an HC3 SE. The one-sided alternative comes from the
  registry. PC2 uses `greater`, following its sign mapping.
- **df.** Bell–McCaffrey/Satterthwaite df for the HC3 contrast.

| Control | Unit | n / n_ctx | β indep (prod) | SE indep (prod) | p indep (prod) |
|---|---|---|---|---|---|
| PC19a | Z0077 (13 members, none removed) | 834 / 82 | −0.041622 (−0.041622) | 0.02063 (0.02063) | 0.02307 (0.02307) |
| PC6 | Z0109 → ARID1B only | 895 / 37 | −0.37530 (−0.37530) | 0.06804 (0.06804) | 7.94e-7 (7.94e-7) |
| PC1 | Z0157 (8 members) | 1007 / 59 | −0.15227 (−0.15227) | 0.02459 (0.02459) | 1.68e-8 (1.68e-8) |
| PC2 | Z0166 (8 members; TP53 and its paralogs are not members) | 987 / 668 | 0.28510 (0.28510) | 0.01587 (0.01587) | 2.47e-53 (2.47e-53) |
| PC17 | Z0018 (21 members; no lineage FE) | 1007 / 50 | −0.40236 (−0.40236) | 0.03288 (0.03288) | 7.69e-18 (7.69e-18) |
| PC15 | Z0186 (5 members; SMARCB1, SS18 and SS18L1 are not members) | 899 / 17 | −0.21875 (−0.21875) | 0.05637 (0.05637) | 4.41e-4 (4.41e-4) |

- **Gene units.** All nine (WRN, PIK3CB, ARID1B, MDM2, PPM1D, USP7, SOX10, MITF, BRD9) also match.
- **Agreement.** Relative differences are ≤ 2e-8 for β and SE. The Satterthwaite df and the p-values agree to
  float32 precision, because the producer caches effects as f32.
- **One convention difference.** The producer merges indicator levels with fewer than 2 models in the analysed
  sample into the reference (`assoc._design`): here, 3 singleton OncotreeLineage levels.
  - With a full lineage FE instead, the singletons are absorbed. β then shifts by at most 1e-3 relative and p by at
    most 2% relative (PC19a: p = 0.02307 under both).
  - The `_merge1` columns reproduce the producer exactly. The shift is not material, but the convention is not
    written in DESIGN §2.5.
- **PC19a.** It is borderline (p = 0.023) but is independently confirmed below 0.05 under both conventions.

## 3. Kill-rule mapping

The reviewer tallied the kill rule independently from `kill_rule.{direction,counted}`, `control.complex_unit` and
`design_sign` in each registry entry and from the C rows of the summary. Evaluable means a frozen unit ID and
n_context ≥ 12. Recovered means the expected sign and one-sided p < 0.05. The tally reproduces `kill_rule.json`
exactly:

| Direction | Counted controls (DESIGN §2.6 = registry) | Evaluable | Recovered | Verdict |
|---|---|---|---|---|
| H3a | PC17, PC20 | 2 | 2 | passes (thin coverage, as the registry requires for 1–2 evaluable controls) |
| H3b | PC1, 2, 4, 5, 9, 10, 11, 13, 18, 19a, 19b | 10 (PC18 not evaluable) | 10 | passes |
| H1 | PC6, 7, 8, 12, 15 | 3 (PC7 and PC8 not evaluable) | 3 | passes |
| H2 | PC21, PC22 | 0 | 0 | n/a; exploratory, as registered in the H2 entry |

- **Gene-level controls.** PC3, PC14 and PC16 are not counted. This is registered in two places:
  - DESIGN §2.6 lists them as "gene-level only: PC3, PC14, PC16".
  - Each of their entries has `kill_rule.direction: none, counted: false` and the rule "gene-level only control
    (DESIGN 2.6): reported, never counted".

  PC16's MLLT1 gene unit (p = 0.063) therefore affects nothing.
- **Non-evaluable controls.** The "not evaluable" labels for PC7, PC8, PC18, PC21 and PC22 were pre-declared in the
  committed registry, and the map confirms each one:
  - SMARCA2, STAG1, KIF18A and CDK6 are unclustered.
  - CCND3 (Z0113) and CDK4/CCND1 (Z0037) are in different coarse clusters.
  - BRD9 (Z0186) and ARID1A/B and DPF2 (Z0109) are in different coarse clusters.
- **Alternatives.** The alternative in every C row agrees with the registry sign. PC5 and PC2 use `greater`.

## Notes (none blocking)

1. **The H1 "complex-level" recoveries are mostly single-gene.** Under the registered paralog-component readout,
   PC6 (ARID1B) and PC12 (VPS4A) each reduce to one gene. PC6's C row is identical to its G row. PC12's differs only
   through the ±3 clipping (one A model has VPS4A < −3). Only PC15 (5-member Z0186) is a multi-gene complex readout. This is the registered
   readout, but the paper should say that the H1 pass rests on 1 multi-gene recovery out of 3.
2. **The singleton-FE convention is code-level only.** Merging fewer than 2-model indicator levels into the
   reference is documented in `assoc._design` but not in DESIGN §2.5. It is immaterial here; DESIGN or the methods
   text should state it.
3. **A summary column is empty.** `covariates_dropped_by_rule` is empty in every row of `controls_summary.tsv`,
   while the gitignored parquets list the dropped levels (e.g. `lineage_FE:Fibroblast;…`). This is cosmetic.
4. **The command is recorded only in abbreviated form.** It appears only in `runtime.json`, with `PART1_DIR=...`
   and a relative `.venv` path, not echoed in `run_log.txt`. It is equivalent to the INTENT command.
