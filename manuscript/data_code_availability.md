# Data and code availability

## Data

**No restricted data are shared, and no raw data are redistributed.** Every input is a public release that readers obtain from its source. The analysis is reproducible from these releases plus the committed code. <!-- DATA.md -->

| Input | Source | Version | Licence |
|---|---|---|---|
| CRISPRGeneEffect, ScreenGeneEffect, CRISPRConfounders, Model, ModelCondition | DepMap portal (depmap.org/portal/data_page) | 26Q3 (user-stated; files carry no version stamp) | CC BY 4.0 |
| OmicsSomaticMutationsMatrixDamaging / Hotspot, OmicsCNGeneWGS, OmicsExpressionTPMLogp1HumanProteinCodingGenes, OmicsProfiles, OmicsGlobalSignatures, OmicsFusionFiltered | DepMap portal | 26Q3 (user-stated) | CC BY 4.0 |
| DEMETER2 combined gene dependency | DepMap portal, "DEMETER2 Data v6" | v6 | CC BY 4.0 |
| CCLE proteomics (protein_quant_current_normalized) | Nusinow et al. 2020 | as published | publisher / DepMap terms (unverified) |
| BioPlex 293T and HCT116 networks | bioplex.hms.harvard.edu (BioPlex 3.0) | Dec 2019 | BioPlex terms (unverified) |
| Human paralogues and gene coordinates | Ensembl BioMart | Ensembl 116 | no restrictions |
| HGNC complete set | genenames.org | download date unknown | CC0 |
| hu.MAP 3.0 test positive and negative pairs | hu.MAP 3.0 (Fischer et al. 2025) | hu.MAP 3.0 test split | CC BY 4.0 per paper (file licence unverified) |

<!-- DATA.md table -->

- DATA.md records the SHA-256 checksum of every input file. `scripts/check_data.sh` verifies a local copy against these checksums before any analysis runs. The DepMap and HGNC checksums equal those of the Part 1 study. <!-- DATA.md -->
- Raw files live in a git-ignored `datasets/` folder, and derived caches in a git-ignored `data_cache/` folder. Neither is committed or distributed.
- The lockbox matrices persisted for crash recovery (`data_cache/part2/lockbox_reads/`) are subsets of the public DepMap CRISPRGeneEffect release. They are not distributed. They are kept locally for the final reproduction check under a retention rule, and their SHA-256 hashes are committed in `results/part2/lockbox_reads.jsonl`. <!-- results/part2/lockbox_confirm/README.md "Retention rule" -->
- CORUM, Complex Portal, STRING and hu.MAP predictions were not used. <!-- DATA.md "Pending external datasets" -->

**Derived results (committed, in the repository):**
- `results/test_eval/`: the single test look, its descriptive export and the independent re-derivation;
- `results/part2/lockbox/`: lockbox and construction-set membership (model identifiers only) and the split manifest;
- `results/part2/construction/`: complex maps, tiers, attachments, selection, variability and the bridge;
- `results/part2/contexts/`, `results/part2/fzx/`: context definitions, unit sets, H1 events and the lockbox-testable list;
- `results/part2/controls/`, `results/part2/directions/`: controls and A-side discoveries, calibrations and sensitivities;
- `results/part2/lockbox_confirm/`: lockbox confirmations, the R6/R7 records, the corrected relook, BENCH and BRIDGE-L;
- `results/part2/robustness/`: the registered per-hit robustness checks, the 254 strongest findings and their post hoc annotations;
- `results/part2/verify/`: independent reviews;
- `results/part2/trials.jsonl`, `results/part2/lockbox_ledger.jsonl`, `results/part2/lockbox_reads.jsonl`: trial and lockbox ledgers.

Manuscript tables and figures are in `manuscript/tables/` and `manuscript/figures/`.

## Code

All code is in the project repository:
- the Part 2 package `src/part2/`;
- the numbered pipeline scripts in `scripts/`;
- figure and table scripts (`scripts/fig_*.py`, `scripts/tab_*.py`);
- independent verification scripts (`scripts/verify/`);
- tests (`tests/`).

The repository also holds the frozen test protocol (`docs/TEST_PROTOCOL.md`), the analysis registry with all amendments (`docs/registry/`) and the design (`DESIGN.md`). The locked Part 1 method is used from its own repository at commit c7f95ee007576be9311c0a5f08c6ac63ca23a009, cloned read-only by `scripts/00_setup_part1.sh`. <!-- DATA.md "Part 1 (locked method)" -->

**To reproduce:**
1. Obtain the public releases above and place them under `datasets/` as listed in DATA.md.
2. Verify them with `scripts/check_data.sh`.
3. Set up the locked Part 1 clone with `scripts/00_setup_part1.sh`.
4. Run the numbered scripts in order, using the command recorded in each step's intent note.

Figures and tables regenerate from committed results alone (`manuscript/supplementary.md`, Supplementary Note 9), and a second run reproduces each figure byte for byte. Seeds are fixed (Methods, Reproducibility) and `PYTHONHASHSEED=0`.

The test look and the lockbox read were single-use by design, and the ledger refuses a second lockbox read. Re-running them therefore requires a fresh repository state. Any such re-run is non-confirmatory.

**Repository URL and archive identifier:** to be added at submission. No public archive has been created yet.
