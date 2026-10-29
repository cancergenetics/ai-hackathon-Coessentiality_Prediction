# hu.MAP 3.0 test split (placeholder)

hu.MAP 3.0 test pair files (Fischer et al. 2025, *Mol Syst Biol*), <https://humap3.proteincomplexes.org/>. Read only by the single registered test evaluation (`scripts/20_score_test_pairs.sh test`, `scripts/21_evaluate_test.py`). The hu.MAP **train** pair files are read from the Part 1 method repository's `datasets/HuMAP/` (see `external/README.md`).

Put these files here. They are not distributed with this repository and are git-ignored; do not commit them.

| File | Version | SHA-256 of the file used | Licence (as recorded in DATA.md) |
|---|---|---|---|
| `Test_positive_pairs.txt` | hu.MAP 3.0 (user-stated) | `c4e3e68edcfc088516b3334dbfae248fb8f276c6475c13dc621abc8b4bba2d96` | CC BY 4.0 per hu.MAP 3.0 paper; licence of these files unverified |
| `Test_negative_pairs.txt` | as above | `c4a3e5468da0790690ee910af1822dcaaa301dea1afaad45816e34f5e40ae9de` | as above |

Check with `shasum -a 256 <file>` (or `scripts/check_data.sh`). Other releases will give somewhat different results.
