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
