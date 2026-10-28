# DATA_REQUEST.md — external datasets requested from the user

Per PROMPT.md, the team does not download datasets itself. The files below are needed for Part 2
(cancer-context complex dependencies). Please place them at the listed paths (folder `datasets/` is gitignored).

## Needed for analysis (evidence-generating)
All DepMap files should come from the **same release as CRISPRGeneEffect.csv (26Q3)**, DepMap portal → Downloads → "DepMap Public 26Q3".

| # | File | Put at | Why |
|---|---|---|---|
| 1 | `OmicsSomaticMutationsMatrixDamaging.csv` | `datasets/DepMAP/` | Damaging (LoF) mutation status per model × gene: mutation contexts (e.g. ARID1A, SMARCA4, KEAP1, STK11). |
| 2 | `OmicsSomaticMutationsMatrixHotspot.csv` | `datasets/DepMAP/` | Hotspot (GoF) mutation status: KRAS/BRAF/PIK3CA/NRAS contexts. |
| 3 | `OmicsCNGeneWGS.csv` *or* `OmicsCNGene.csv` (gene-level relative copy number; whichever the release has) | `datasets/DepMAP/` | Copy-number contexts (deletion / partial loss, CYCLOPS-type collateral vulnerability, amplification). |
| 4 | `OmicsExpressionProteinCodingGenesTPMLogp1.csv` (name may be `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv`) | `datasets/DepMAP/` | Expression of subunits: controls for non-expressed genes, paralog-expression contexts, LoF-by-silencing. |
| 5 | `OmicsProfiles.csv` (or `OmicsDefaultModelProfiles.csv`) | `datasets/DepMAP/` | Maps omics profile IDs to ModelIDs if the matrices are keyed by profile. |
| 6 | `OmicsSignatures.csv` (if present; MSI score, ploidy, aneuploidy, LoH fraction) | `datasets/DepMAP/` | Genome-instability contexts (MSI, whole-genome doubling). Optional. |
| 7 | BioPlex 3.0 interaction networks: `BioPlex_293T_Network_10K_Dec_2019.tsv` and `BioPlex_HCT116_Network_5.5K_Dec_2019.tsv` (bioplex.hms.harvard.edu → Downloads) | `datasets/BioPlex/` | Experimental AP-MS PPIs to integrate with co-essentiality for complex construction (allowed by PROMPT). |

## Requested for post hoc interpretation only (never used to build or validate complexes)
| # | File | Put at | Why |
|---|---|---|---|
| 8 | CORUM 5.x human complexes (`corum_humanComplexes.txt` or `allComplexes.txt`, mips.helmholtz-muenchen.de/corum → Download) | `datasets/posthoc/` | Post hoc naming/annotation of discovered complexes and novelty assessment. |
| 9 | Complex Portal human complexes (`9606.tsv`, ebi.ac.uk/complexportal → Download) | `datasets/posthoc/` | Same, second source. Optional. |

Notes
- If 26Q3 files are not available, any DepMap release is acceptable *provided* you tell us which; we will restrict to ModelIDs present in CRISPRGeneEffect.csv and record the mismatch in DATA.md.
- Items 1–5 and 7 are the priority. Items 6, 8, 9 are optional.

## Round 2 requests (2026-10-28, after LITERATURE.md §8)
User confirmed items 1–7 will be provided; CORUM/Complex Portal (8, 9) skipped.

| # | File | Put at | Why |
|---|---|---|---|
| 10 | DepMap 26Q3 `OmicsFusionFiltered.csv` | `datasets/DepMAP/` | Fusion contexts with known complex biology: SS18::SSX → ncBAF/BRD9, EWSR1::FLI1, KMT2A-r → menin/DOT1L/ENL, PAX3/7::FOXO1. |
| 11 | DepMap RNAi DEMETER2 `D2_combined_gene_dep_scores.csv` (McFarland 2018; DepMap portal "DEMETER2 Data v6") | `datasets/DepMAP_RNAi/` | Orthogonal perturbation modality to replicate complex–context associations. |
| 12 | Ensembl BioMart human paralogues (Ensembl gene ID, paralogue gene ID, %identity; any release ≥110) | `datasets/Paralogs/` | Sequence-based paralog pairs for paralog-switching within complexes (not curated complexes). |
| 13 | CCLE proteomics `protein_quant_current_normalized.csv` (Nusinow 2020) | `datasets/Proteomics/` | Optional: subunit co-abundance / collateral loss validation. |
