# results/part2/lockbox_confirm

Outputs of the lockbox confirmation step: the single read of lockbox L. Notes:
- `BLOCKED.md`: the dry-check stop;
- `GAPS.md`: the gaps and their A4 and A5 rulings;
- `REHEARSAL.md`: the dress rehearsal;
- `INTENT.md`: the intent note committed before the read.

Per entry, `<ID>/confirmation.tsv` and `<ID>/summary.json` are written by `scripts/51_lockbox_confirm.py`. The rest of
the read's record lives elsewhere:
- target lists: `results/part2/lockbox_targets/`;
- the ledger: `results/part2/lockbox_ledger.jsonl`;
- hashes of the persisted reads: `results/part2/lockbox_reads.jsonl`;
- the engine output: `results/part2/assoc/<ID>/confirmation_crispr.json`.

## Retention rule for persisted L reads (PI, R1)

Each lockbox read is persisted, for crash recovery only, to `data_cache/part2/lockbox_reads/<ID>__crispr.parquet`.
The directory is gitignored. Its SHA-256 is committed to `results/part2/lockbox_reads.jsonl` right after the read.
- The copies stay **only in the gitignored local `data_cache`** until the final reproduction check of the manuscript is
  done.
- They are **never copied, committed, or read by any other code**. Only `assoc.read_lockbox` (writing),
  `assoc.resume_lockbox_read` and `confirm.confirm_entry` (the registered `--resume` mode of the driver) may touch
  them. `tests/test_part2_lockbox_firewall.py` enforces this for every module and script.
- Deleting them later requires a **PI note in the research tree**. The committed hashes stay as provenance.

## PI ruling R6 (verbatim; the PI note that authorises the persisted-read access of the corrected refit)

"R6 (PI, 2026-10-29). The L-side covariate design was built over all models and then subset to L. That is an implementation defect, not an analysis choice; the registered analysis fits within L using the registered covariates. (1) Primary reporting: the 53 failed fits (H3b 28, H2 1, H1 FOCAD, H1S 1, controls 13 excluding PC21, BRIDGE-L GLS 9) are relabelled NOT EVALUABLE (driver defect), not 'not confirmed'. Primary rates use the evaluable denominator, and the all-in conservative rate is shown alongside. (2) A post hoc corrected refit is authorised as a disclosed, NON-CONFIRMATORY relook. It may read the persisted reads only through the driver's resume path (never the loader), only for the failed fits plus the reproduction check in (3), and only after the fix and an intent note are committed."

Access path: `scripts/51_lockbox_confirm.py --relook R6 --mode check|refit`. This calls
`assoc.resume_lockbox_read(..., relook="R6")`, which requires this note and `CORRECTED_REFIT_INTENT.md`, both
committed. It never calls the loader and writes no ledger line.
