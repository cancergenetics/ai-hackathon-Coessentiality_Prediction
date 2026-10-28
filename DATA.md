# DATA.md — datasets used

All raw files live in `datasets/` (gitignored; never committed). Derived caches in `data_cache/` (gitignored).
Checksums: SHA-256 via `shasum -a 256`, computed 2026-10-28 (re-check: `scripts/check_data.sh`).

| File | Source | Version | SHA-256 | License |
|---|---|---|---|---|
| DepMAP/CRISPRGeneEffect.csv | DepMap Portal (Broad), depmap.org/portal/data_page | 26Q3 (as stated by user; files carry no version stamp) | 1e0c492944a7297c0eda179081a1e81a1cfe41a97628facd7eddbdda1a7504ee | CC BY 4.0 (DepMap public data terms) |
| DepMAP/ScreenGeneEffect.csv | DepMap Portal | 26Q3 (user) | 08e15ab81d7450db70c48c7781f312f6cebc068da7077bea775569fb516a3f04 | CC BY 4.0 |
| DepMAP/CRISPRConfounders.csv | DepMap Portal | 26Q3 (user) | 0e098325d53403d7852f3a773477d20a508a9bf29f43e10fe28bf0f550b35779 | CC BY 4.0 |
| DepMAP/Model.csv | DepMap Portal | 26Q3 (user) | 65eb1895800e71a963a3ab56e5d08c8d006cb0d9e7aa68f7e34a3810a2a3cbc7 | CC BY 4.0 |
| DepMAP/ModelCondition.csv | DepMap Portal | 26Q3 (user) | a1edda78cd60718e5c99690ce8ab80315b65ddc9346643338d2c1c774fc51326 | CC BY 4.0 |
| HGNC/hgnc_complete_set.txt | HGNC, genenames.org | download date unknown (provided by user) | 0f2d3f035e0e6583c055a1d9e008f275376a4b66f72af380dfc2ac730272c247 | CC0 |
| huMAP/Test_positive_pairs.txt | hu.MAP 3.0 test split (Fischer et al. 2025), provided by user 2026-10-28 13:55 UTC | hu.MAP 3.0 (user-stated) | c4e3e68edcfc088516b3334dbfae248fb8f276c6475c13dc621abc8b4bba2d96 | CC BY 4.0 per hu.MAP 3.0 paper; licence of these files unverified |
| huMAP/Test_negative_pairs.txt | as above, 2026-10-28 13:55 UTC | as above | c4a3e5468da0790690ee910af1822dcaaa301dea1afaad45816e34f5e40ae9de | as above |

DepMap and HGNC checksums are identical to those recorded in Part 1 (`external/part1_locked/DATA.md`).

## Part 1 (locked method) — code and derived predictions
| Item | Source | Version | Checksum | Notes |
|---|---|---|---|---|
| Part 1 repository | local read-only folder `autoresearch/coessentiality_protein complexes` | git commit c7f95ee007576be9311c0a5f08c6ac63ca23a009 (2026-10-28 12:55:43 UTC) | git hash | cloned to `external/part1_locked/` (gitignored) by `scripts/00_setup_part1.sh` |
| Genome-wide pair scores (`outputs/all_pairs/zavg_BEF.npy`, `z_B/z_E/z_F.npy`, `genes.tsv`) | produced by Part 1 `mvp/scripts/predict_all_pairs.py` from commit 0287470 | 2026-10-28 12:44 UTC | see `external`-side `manifest.json`; recorded by scripts that load them | label-free; read in place from the Part 1 folder (not copied, not modified) |

## Pending external datasets
See DATA_REQUEST.md. Rows are added here when the user provides files.

## Notes / explicit assumptions
- hu.MAP pairs are UniProt accessions. Mapping as Part 1: UniProt → HGNC (Approved, `uniprot_ids`) → Entrez → DepMap column `SYMBOL (ENTREZ)`.
- Test files were created in this folder at 2026-10-28 13:55 UTC, one hour after the final Part 1 commit (12:55:43 UTC). Before the test protocol was frozen, the PI read the line counts (3,312 positive / 597,426 negative lines) and the first 3 lines of each file to check the format; no scores were computed (decision note D0 in RESEARCH_TREE.md).
