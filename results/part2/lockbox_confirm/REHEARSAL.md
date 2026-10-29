# Lockbox driver: dress rehearsal (review MUST-2c)

Date: 2026-10-29.

**No lockbox value was read.** The rehearsal used a dummy loader that returns random values: for the lockbox reads it
received only the 338 L model IDs and the gene keys.

| | |
|---|---|
| Code | `5082cc8` (the clone's HEAD at the start; the rehearsal then added its own ledger, target and output commits inside the clone) |
| Where | A throw-away `git clone` in the session scratchpad (not a worktree), marked `.lockbox_rehearsal`. Its `data_cache/part2` was a real directory, so the GLS tables and the persisted reads stayed inside the clone. The clone was deleted afterwards. |
| Command | `caffeinate -i .venv/bin/python scripts/verify/lockbox_rehearsal.py <out.json>` (with `PART1_DIR` and `PYTHONHASHSEED=0` set) |
| Inputs | Real A data: CRISPR and DEMETER2 on A through `assoc.load_A_effects`, plus the committed registry and frozen inputs. |
| Dummy loader | 30 calls, all `crispr`, all for the 338 L model IDs: 28 confirmatory entries, BENCH and BRIDGE-L. |

## Result

The whole driver ran end to end with no error. The preflight reported no gaps.

- All 30 entries have a ledger line, in registry order: 23 controls by index, then H3a, H3b, H1, H1S, H2, BENCH and
  BRIDGE-L.
- Every read was persisted, and its hash committed.
- All 30 entries were closed after their outputs were committed.

The estimates are meaningless (random L values) and are not reported.

| Stage | Wall time |
|---|---|
| Building the world (A load) | 15 s |
| Preflight, including BENCH A side (g\*, breadth, eigengene fits) and BRIDGE-L frozen inputs and coverage | 148 s |
| 23 controls | 37 s |
| H3a / H3b / H1 / H1S / H2 confirmations | 85 / 14 / 15 / 148 / 8 s |
| BENCH (read, primary bootstrap B 2000, multiverse point estimates) | 8,305 s (2.3 h) |
| BRIDGE-L (GLS discovery on A, about 10 min, then its read, confirmation and estimate with bootstrap CIs) | 849 s |
| **Total** | **9,734 s (2.7 h)**, peak RSS 3.64 GB |

**Disclosure.** The rehearsal ran the registered BRIDGE-L GLS-map discovery on real A data, with the dummy affecting
L only. The A-side outcome was printed in the log: 1,206 and 139 GLS discoveries.
- Its trial line and tables existed only in the deleted clone.
- The registered single pass will run in the real read, with identical, deterministic machinery and seed 4500.
- Seeing the count of an A-side, fully mechanical pass changes no choice.
