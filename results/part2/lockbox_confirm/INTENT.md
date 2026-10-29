# Lockbox confirmations: intent (single read)

This note is committed **before any lockbox read**. The ledger does not exist and no target list has been written. It
was written at 2026-10-29T06:46:29Z (UTC).

| | |
|---|---|
| Code commit | `7d3c23a` (R1, after the final GO in review-lockbox `209fd8d`) |
| Registry | The committed entries with amendments A3, A4 and A5 (`docs/registry/AMENDMENTS.md`). |
| Inputs | Committed discovery tables (data_cache parquet, hash-checked against the committed summaries); `classes.tsv`, `controls_summary.tsv`, the D8 calibration tables and `lockbox_testable.tsv` (FZ-X manifest hash), all committed. |
| Rehearsal | `REHEARSAL.md`: end to end in 2.7 h, peak 3.6 GB. |

## Command (one run; from the worktree root)

```
export PART1_DIR="<PART1_DIR>"; export PYTHONHASHSEED=0
nohup caffeinate -i <PROJECT_ROOT>/.venv/bin/python \
  scripts/51_lockbox_confirm.py --run > <scratch>/lockbox_run.log 2>&1 &
# wait by PID (kill -0); never pgrep -f
```

If the run fails mid-way, the only allowed continuation is the registered resume mode, which never calls the loader:

```
caffeinate -i <PROJECT_ROOT>/.venv/bin/python scripts/51_lockbox_confirm.py --run --resume
```

If resume cannot complete an entry, the run stops and the failure is reported. Nothing is improvised.

## Read order (registry; enforced by the driver and `assoc.read_lockbox`)

1. CTRL_PC1 … CTRL_PC22, by registry index. Lockbox permutation seeds are 6700 + index.
2. H3a (6100), H3b (6200), H1 (6300), H1S (6400), H2 (6500).
3. BENCH_complex_vs_gene: B 2000, seed 5100.
4. BRIDGE-L_confirmation_rate. Its GLS discovery on A runs first (seed 4500, expected counts 1,206 / 139, as in the
   rehearsal), then its read (seed 6600), then the estimate (B 2000, seed 5200).

## Hypotheses per entry

The full list is in `intent_hits.tsv.gz` (SHA-256 `bba0041ea86a046f756a4f69fe16f69d0594443eadc6c14a20091049df49cdb1`).
Each row has context, unit, beta_A, alternative, eligible and route. H1/H1S rows add the event, axis, engine class,
label_final and headline flag; BENCH rows add g\*, g\*(beta), breadth and tercile. The file was produced by the
driver's own spec builders on A only, with `read_lockbox` disabled.

| Entry | Hits | Eligible | Lockbox route | Rare-context route |
|---|---|---|---|---|
| 23 controls | 49 units | 49 | 44 | 5 (PC15 ×2, PC16 ×2, PC22; contexts SMARCB1_LoF_or_SS18SSX, KMT2A_r_vs_nonNPM1c) |
| H3a_lineage_atlas | 1,145 discoveries | 1,145 | 1,145 | 0 |
| H3b_driver_family | 161 discoveries | 161 | 150 | 11 (fusion EWSR1::FLI1\|ERG; KMT2D has no discovery) |
| H1_within_complex_GI | 24 event axes | 5 non-neutral calls | 4 | 1 (Z0186\|EP300\|P) |
| H1S_silencing_GI (exploratory) | 293 event axes | 24 non-neutral calls | 18 | 6 N axes, on the A CRISPR alignment (A5) |
| H2_submodule_divergence (exploratory, kill rule) | 109 discoveries | 109 | 109 | 0 |
| BENCH_complex_vs_gene | 1,295 primary pairs (H3a 1,145 + H3b 150) | – | – | – |
| BRIDGE-L_confirmation_rate | GLS discoveries from its registered pass | – | – | – |

**H1 eligible calls.** Each is confirmed iff it has the same sign and one-sided p < 0.05 in L, with ≥ 5/5 L models:

| Call | Engine class | label_final | Route |
|---|---|---|---|
| Z0109\|ARID1A\|P | PARALOG COMPENSATION | PARALOG COMPENSATION | lockbox |
| Z0186\|EP300\|P | PARALOG COMPENSATION | PARALOG COMPENSATION | rare-context route |
| Z0286\|FOCAD\|N | RESIDUAL SENSITISATION | RESIDUAL SENSITISATION | lockbox |
| Z0130\|PTEN\|N | COLLAPSE | DETECTABLE, WITHIN EQUIVALENCE MARGIN | lockbox (E2: never promoted) |
| Z0158\|TP53\|N | RESIDUAL SENSITISATION | DETECTABLE, WITHIN EQUIVALENCE MARGIN | lockbox (E2: never promoted) |

The remaining 19 H1 event axes are neutral, inconclusive or not evaluable. They get L effects with 95% and 90% CIs
and are not confirmed.

## Reporting (fixed now)

- **Confirmation rate per direction:** overall, non-locus, and split by "beats all matched sets" (H3a/H3b, D8 flags) or
  "beats all its random splits" (H2).
- **Rare-context hits:** labelled "not lockbox-testable", with pass, fail or NOT EVALUABLE (DEMETER2 not available).
- **H1/H1S per-event table:** `<ID>_events.tsv`.
- **BENCH:** the primary estimate with its CI and terciles; Spearman(Delta, breadth); the distinct-complex count;
  multiverse point estimates. `within_line_rank` and the A–L random-set correlation are NOT RUN.
- **BRIDGE-L:** the primary estimate with its CI, and the multiverse with bootstrap CIs.
- **Persisted reads:** kept under the retention rule (`README.md`).
