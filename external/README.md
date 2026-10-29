# external/ (placeholder — not distributed)

Everything here is git-ignored except this file.

## `external/part1_locked/` — the locked Part 1 method

`scripts/00_setup_part1.sh` creates it. Part 2 imports the method from this clone unchanged, and refuses to run on any
other code (`src/part2/part1.py`, `src/testeval/part1.py`).

- **Development route.** The original Part 1 repository at the locked commit
  `c7f95ee007576be9311c0a5f08c6ac63ca23a009`. Set `PART1_DIR` to its folder.
- **Public route.** The public release of the method, <https://github.com/cancergenetics/ai-hackathon-Coessentiality>,
  at commit `3c684c972f1e5c3ead722b4b3e2c508fa3cb5ddf`. Its git history is fresh, so `c7f95ee` is not in it. Every
  Part 1 path this project reads (`src/`, `mvp/`, `tests/`, `scripts/final_eval.py`, `scripts/predict_pairs.py` and
  three Part 1 result files) has the same git object id as at `c7f95ee`. The ids are listed as `LOCKED_OBJECTS` in
  `src/part2/paths.py` and `src/testeval/paths.py`, and the setup script and lock checks verify them.

## Part 1 outputs read in place (`PART1_DIR`)

Some inputs come from a Part 1 working copy, which `PART1_DIR` points to (default `../coessentiality_protein complexes`):

- `outputs/all_pairs/`: the genome-wide zavg_BEF scores and `genes.tsv`. Produced by
  `mvp/scripts/predict_all_pairs.py` of the method repository from DepMap 26Q3 `CRISPRGeneEffect.csv` and HGNC.
- `data_cache/splits/protein_assignment.csv`: produced by the Part 1 pipeline. Its SHA-256 must equal
  `results/splits/protein_assignment.sha256`.
- `datasets/HuMAP/Train_positive_pairs.txt`, `Train_negative_pairs.txt`: the hu.MAP 3.0 train pairs, with SHA-256
  pinned in `src/testeval/paths.py`.

For the public route, clone the method repository and follow its README to place its datasets and produce these
outputs. Then set `PART1_DIR` to that clone.
