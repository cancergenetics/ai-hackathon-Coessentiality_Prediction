# Part 1 audit: no hu.MAP test information accessed during development

Auditor: independent agent. Script: `scripts/verify/audit_part1.py` (reproducible; writes `results/audit/part1_audit.json`).
Part 1 was only read (working tree and `.git`; history queries ran in the read-only clone `external/part1_locked`, HEAD c7f95ee, history identical to Part 1's 60 commits). The test lists were used for set membership only; no score, metric or prediction was computed on test pairs.

## Transition timeline (UTC)
| When | Event |
|---|---|
| 2026-10-23 11:13 / 11:56 | Part 1 `datasets/HuMAP/Train_positive_pairs.txt` / `Train_negative_pairs.txt` created (SHA-256 equal to those in Part 1 DATA.md). Folder holds only Train files. |
| 2026-10-27 12:16 | First Part 1 commit (M0). Dev-only design, Phases A-C, one internal-final-test (IFT) look at 17:36 (built from Train pairs). |
| 2026-10-28 12:55:43 | **Last Part 1 commit** c7f95ee (60 commits, branch `master` only). Newest non-Finder Part 1 file mtime 12:55:35. |
| 2026-10-28 13:55:21 / 13:55:59 | **Test files created** (birth time) in the new project; SHA-256 equal to the expected c4e3e68e... (pos) and c4a3e546... (neg). Gap to last commit: 59.6 min. |
| 2026-10-28 14:03 | Only later Part 1 modification: `.DS_Store` (Finder metadata). |

## Verdicts
**1. Timeline: PASS.** Test files were born 59.6 min after the last Part 1 commit; 0 commits (author or committer date) post-date their birth. 0 Part 1 files (excluding Finder `.DS_Store`) have mtime after the last commit; no file in `.git` changed after 12:55:43 (only the `.git` directory mtime moved, from transient lock files of my read-only git calls).

**2. Code/doc references: PASS (1 NOTE).** Searched all 60 commits and the working tree (excluding `.env/`, binaries skipped) for `Test_positive`, `Test_negative`, `Test_pairs`, `HuMAP/Test`, the file names and both SHA-256 prefixes. Hits: `tests/test_metrics.py` and `tests/test_baselines.py` function names `test_negative_control_*` (innocent); `REPORT.md` line 164 `--pairs Test_pairs.txt`, a generic placeholder in the instructions for the user to score the withheld set later (not a path to the real file). A broader pattern `test_pos|test_neg` matches 470 results JSON files and 9 code files: all are internal CV fold counts (`test_pos`/`test_neg` keys). No Part 1 file is byte-identical to either test file (checked by size and hash). Part 1 DATA.md, PROMPT.md and docs/AGENTS.md state the test set was withheld and never accessed.

**3. Content leakage: PASS (2 NOTES).**
- `data_cache/pairs.pkl` (634,779 rows) equals exactly Train pos plus Train neg; it contains 0 of the test pairs. The four `data_cache/layers/*.pairs.npz` gene-pair universes (596,597 pairs) equal the scorable Train pairs, and none of the test-only gene pairs is in them. `protein_assignment.csv` and `uniprot_map.csv` (1,151 proteins each) contain 0 test-only proteins.
- Text scan of 2,176 working-tree files and 762 historic blobs for test pairs on one line / list: no label, split, IFT, trial or doc file contains one. Hits only in (a) label-free score listings (`results/predictions/genome_wide_top_pairs.tsv.gz`, a genome-wide ranking of all gene pairs, expected to include test pairs, carries no labels; counts deliberately not reported: they would amount to precision/recall on test pairs), and (b) HGNC / `genes.tsv` lines where one gene lists two accessions (coincidence).
- Scoring code (`mvp/src/coessential/method.py`, `all_pairs.py`, `mvp/scripts/predict_all_pairs.py`, `scripts/predict_pairs.py`, `src/predict.py`): imports and code lines contain no label-file reads; `predict.py` uses only `data.load_matrix/build_matrix_cache/build_uniprot_map`, none of which references label files; run records state "no pair label is read". NOTE: `src/config.py` defines the Train file paths and `predict.py` imports `src.data`, but labels are never loaded on that path.

**4. Overlap descriptives: NOTE (information, no scores).**
| Quantity | Value |
|---|---|
| Test pos / neg (unordered, unique) | 3,312 / 597,426 (no duplicates, no self-pairs) |
| Pairs in both test pos and neg | 0 |
| Test pairs also in Part 1 Train (pos or neg) | 0 / 0 (pos), 0 / 0 (neg): pair-level disjoint |
| Prevalence test vs train | 0.551% vs 0.522% (train: 3,314 pos / 631,465 neg) |
| Test proteins (pos or neg) | 1,121; in Train proteins: 329 (322 in Train-positive proteins) |
| Test positives: C1 / C2 / C3 (both / one / no protein in Train positives) | 265 / 738 / 2,309 |
| Test proteins mapped to DepMap (HGNC, Entrez, `SYMBOL (ENTREZ)`) | 1,087 / 1,121 (29 not in DepMap, 3 not in HGNC, 2 ambiguous) |
| Scorable pairs (both mapped, distinct genes) | pos 3,159 / 3,312; neg 563,648 / 597,426 |
The mapping reproduces Part 1's `uniprot_map.csv` for all 1,151 Train proteins; new-project DepMap and HGNC files match the hashes in Part 1 DATA.md.

**5. Hyper-parameter selection: PASS.** All 845 logged trial starts have `universe = dev` (Train-derived cluster-grouped CV; 813 cluster, 12 random, 11 gene, 3 component, 6 unspecified). The only evaluation beyond dev was the IFT (847 pos / 40,786 neg, a Train-derived held-out set) with one pre-registered look (`results/ift_looks.jsonl`); the final method `zavg_BEF` has no fitted parameters. NOTE (D1): two round-1 reviewers computed label-free AP probes on all scorable Train pairs before the IFT existed; this compromises IFT independence, not the hu.MAP test.

**6. Other channels: NOTE.** Part 1 LITERATURE.md/DESIGN.md describe the hu.MAP 3.0 test construction from the public paper (random 50:50 split of Complex Portal complexes, negatives from different complexes) and a guessed test prevalence (about 0.5%; actual 0.551%). This is public knowledge and shaped the choice of complex/cluster-disjoint validation; it is not test data. Reviewer docs repeatedly forbid any use of test information.

## Overall
No evidence that test pairs, labels or test-derived information entered Part 1. Residual notes: public paper knowledge, the placeholder `Test_pairs.txt` command, genome-wide label-free score files that inevitably cover test pairs, and my own read-only access to Part 1 (one `git status` created and removed a lock file; index/refs unchanged).
