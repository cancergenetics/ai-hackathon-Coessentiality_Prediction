# Part 2 context table (DESIGN v3.1 §2.2)

These files were built by `scripts/32_build_contexts.py`, using `src/part2/contexts.py`, at code commit 10b5f90.
- **Inputs:** DepMap Omics (damaging and hotspot mutation matrices, OmicsCNGeneWGS, expression, global signatures, fusions, OmicsProfiles), Model.csv, Ensembl 116 coordinates, HGNC locations, and CCLE proteomics (used only for the LoF functional check).
- **No dependency value is read.** U (1,345 models) comes from the ModelID column of the CRISPR file only.
- **Runtime:** 88 s from raw CSVs and 40 s from the parquet cache (laptop, single process).
- **Reproducibility:** two consecutive builds gave byte-identical outputs (`checksums.json`).

## Committed here

| File | Content |
|---|---|
| `context_summary.tsv` | One row per context column: kind, n_U, n non-NaN, n positive and definition. |
| `lof_functional_check.tsv` | Pre-registered LoF functional check. It covers every gene with ≥ 8 MSS LoF lines in U, plus the named LoF genes. For each gene it compares expression (log1p TPM) and protein (CCLE proteomics) in MSS LoF lines vs MSS non-LoF lines, using one-sided Mann–Whitney. A gene passes if expression is lower (p < 0.05, median difference < 0) **or** protein is lower (p < 0.05, ≥ 3 LoF lines with protein). Result: 118 of 225 pass. |
| `segments_ge8.tsv` | Hom-del and high-amp segments with ≥ 8 events in U. Adjacent genes in Ensembl order are merged when their event vectors have Jaccard ≥ 0.9. A model carries a segment event if ≥ 50% of the segment's observed genes have the event. Segments are named by HGNC cytoband and first–last gene. The loss and amp segment tables (16k and 14k mostly single-gene segments) are only in the cache. |
| `build_log.json` | Default-profile counts, the CN-scale check, the LoY threshold, the TP53 split, fusion counts, the segments of the named CN genes and the functional-check summary. |
| `checksums.json` | SHA-256 of every output, including the gitignored parquet files under `data_cache/part2/contexts/`. It also holds content hashes, input SHA-256s, the universe hash and the code commit. It is the freeze manifest that `run_association` checks inputs against. |

## Decisions and checks worth knowing

- **Default profiles.** For each Omics matrix, the row with `IsDefaultEntryForModel == Yes` is used; it is unique per model. Coverage in U:
  - DNA (mutations, MSIScore): 1,345 models.
  - WGS copy number: 1,208 models. The other 137 are WES-only, so they have no CN.
  - RNA (expression, fusions): 1,276 models.
- **Ploidy and WGD.** These come from the default DNA profile. 299 WGS lines have NaN Ploidy and WGD, so WGD is non-NaN in 1,032 models. The primary analysis keeps these as NaN, per the complete-case rule. `Ploidy_anyprofile` and `WGD_anyprofile` borrow the values from another DNA profile of the same model; they are a sensitivity check only.
  - In those lines, hom-del uses r < 0.25, and partial loss is NaN.
  - So `VPS4B_loss` is non-NaN in only 909 models.
- **CN scale check.** The median r of housekeeping genes (73 RPL/RPS genes) is 1.036 in WGD lines and 1.041 in non-WGD lines, so r is ploidy-relative, as the design assumes.
- **Sex chromosomes.** In male lines, the median r is 0.57 for X and 0.30 for Y.
  - For male X (non-PAR) and Y genes, partial-loss and any-loss events are set to NaN. Hom-del (a < 0.5) is kept.
  - In non-male lines, chrY CN events and chrY LoF are NaN.
  - Male X-linked LoF is damaging ≥ 1.
- **LoY.** The score is the sum of log1p TPM over 43 male-specific chrY protein-coding genes. The threshold is the Gaussian-KDE antimode of the male distribution, 16.54. The modes are at 1.1 and 37.5; the female 99th percentile is 5.0, and no female line lies above the threshold.
  - 290 of 694 male lines have LoY.
  - The call agrees with "DDX3Y < 1" in 97.4% of male lines.
- **TP53 status.** 887 lines are altered and 428 are WT. 15 lines with a single damaging hit (no loss, hom-del or hotspot) are NaN as ambiguous, and 15 more are unknown.
- **Fusions.** Gene1 and Gene2 are parsed from 'SYMBOL (ENSG…)'. The other formats present, 'SYMBOL (.)' and viral IDs, are tolerated. Counts in U match the expected minimums exactly: EWSR1::FLI1|ERG 19, SS18::SSX 5, KMT2A-r 13, PAX3|7::FOXO1 7, BCR::ABL1 12.
- **Named LoF genes failing the functional check.** These are dropped from H1, per DESIGN §2.2:
  - CREBBP: n = 12; neither expression nor protein is lower.
  - Genes with no or almost no MSS LoF lines: CUL3 (0), VPS4B (0), and the MMR genes MLH1, MSH2, PMS2, PMS1, MLH3, MSH6, EPCAM.
  - KMT2D, KEAP1 and EP300 pass on protein only.
- **Hom-del list.** Common germline deletion polymorphisms (LCE3B/C, UGT2B17, SIGLEC14) appear as frequent hom-dels. They fail the functional check, because these genes are not expressed.
