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
| DepMAP/OmicsSomaticMutationsMatrixDamaging.csv | DepMap Portal, provided by user 2026-10-28 | 26Q3 (user-stated) | 736b7a82b4d1f234dc597a409d104186c282f01b0caa56365c27ff987ae49e20 | CC BY 4.0 |
| DepMAP/OmicsSomaticMutationsMatrixHotspot.csv | DepMap Portal, user 2026-10-28 | 26Q3 (user) | 0cbdacdce2355d2dd600d3741a36ae093dbe77d7bee44bde97411b1e96973058 | CC BY 4.0 |
| DepMAP/OmicsCNGeneWGS.csv | DepMap Portal, user 2026-10-28 | 26Q3 (user); ploidy-relative linear CN ratio (diploid ≈ 1), not log2 (round-1 reviews) | 31cace66f202adaace306d975b5661c8b77ea82f07d2ed34045eec9c7958fcd6 | CC BY 4.0 |
| DepMAP/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv | DepMap Portal, user 2026-10-28 | 26Q3 (user) | 252d1fedc7014953b19ffe5d68151b4df973c8078b12c9e5ccd97e292548d177 | CC BY 4.0 |
| DepMAP/OmicsProfiles.csv | DepMap Portal, user 2026-10-28 | 26Q3 (user) | 4241e537f0b23139c0e525b4d78fb53ee65c5507771a33e2cd5debd78f7179ea | CC BY 4.0 |
| DepMAP/OmicsGlobalSignatures.csv | DepMap Portal, user 2026-10-28 | 26Q3 (user) | 51d242e99176e82238df3c90a62b848b40a459a4a4981a2b0a923ed46ad8fb93 | CC BY 4.0 |
| DepMAP/OmicsFusionFiltered.csv | DepMap Portal, user 2026-10-28 | 26Q3 (user) | 985ac6383a95c39f36856c583f32c7990c9f5464d982a2dd84a9dc07fcbc979b | CC BY 4.0 |
| DepMAP_RNAi/D2_combined_gene_dep_scores.csv | DepMap Portal "DEMETER2 Data v6" (McFarland et al. 2018), user 2026-10-28 | DEMETER2 v6 (user-stated) | 6bfe7b24ef191debccd502b7dfe61c4b0763caba3e02c6c14e521c8226bb98a7 | CC BY 4.0 |
| Proteomics/protein_quant_current_normalized.xlsx | CCLE proteomics, Nusinow et al. 2020 (Cell), user 2026-10-28 | as published (xlsx, not csv as requested) | a0ec29ad609d309d884852c268e3382ea6b38bca171af2cf59e4f38a518d21ff | per publisher / DepMap terms (unverified) |
| BioPlex/BioPlex_293T_Network_10K_Dec_2019.tsv | BioPlex 3.0 (Huttlin et al. 2021, Cell), bioplex.hms.harvard.edu, user 2026-10-28 | Dec 2019 release | 5a98def087fc77c7bf07aac79e2553fd64ea1a27a7d0c4617870476f3f120833 | BioPlex terms (free for research; unverified) |
| BioPlex/BioPlex_HCT116_Network_5.5K_Dec_2019.tsv | as above | Dec 2019 release | b37be6de6fe15e948151ab9480cfcb88a6795780cf7a4059b386b07f3c1e760d | as above |
| Paralogs/ensembl116_human_paralogues.tsv | Ensembl BioMart human paralogues, user 2026-10-28 | Ensembl 116 | 003adc6bb6e0f1ad3500aee537a4df99b6d63fafc32ef024835c76f4c9570700 | Ensembl: no restrictions |

DepMap and HGNC checksums are identical to those recorded in Part 1 (`external/part1_locked/DATA.md`).

## Part 1 (locked method) — code and derived predictions
| Item | Source | Version | Checksum | Notes |
|---|---|---|---|---|
| Part 1 repository | local read-only folder `autoresearch/coessentiality_protein complexes` | git commit c7f95ee007576be9311c0a5f08c6ac63ca23a009 (2026-10-28 12:55:43 UTC) | git hash | cloned to `external/part1_locked/` (gitignored) by `scripts/00_setup_part1.sh` |
| Genome-wide pair scores (`outputs/all_pairs/zavg_BEF.npy`, `z_B/z_E/z_F.npy`, `genes.tsv`) | produced by Part 1 `mvp/scripts/predict_all_pairs.py` from commit 0287470 | 2026-10-28 12:44 UTC | see `external`-side `manifest.json`; recorded by scripts that load them | label-free; read in place from the Part 1 folder (not copied, not modified) |

## Pending external datasets
All requested items 1–7 and 10–13 provided by 2026-10-28 16:20 UTC. Optional items 8–9 (CORUM, Complex Portal) not provided: post hoc interpretation uses HGNC gene_group + literature only.

## Notes / explicit assumptions
- hu.MAP pairs are UniProt accessions. Mapping as Part 1: UniProt → HGNC (Approved, `uniprot_ids`) → Entrez → DepMap column `SYMBOL (ENTREZ)`.
- Test files were created in this folder at 2026-10-28 13:55 UTC, one hour after the final Part 1 commit (12:55:43 UTC). Before the test protocol was frozen, the PI read the line counts (3,312 positive / 597,426 negative lines) and the first 3 lines of each file to check the format; no scores were computed (decision note D0 in RESEARCH_TREE.md).
