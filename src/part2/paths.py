"""File locations for Part 2 (one place; no hard-coded paths elsewhere in the package)."""
from __future__ import annotations

import os
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[2]
DATASETS = PROJECT / "datasets"
DEPMAP = DATASETS / "DepMAP"
GENE_EFFECT_CSV = DEPMAP / "CRISPRGeneEffect.csv"
MODEL_CSV = DEPMAP / "Model.csv"
CONFOUNDERS_CSV = DEPMAP / "CRISPRConfounders.csv"
HGNC_FILE = DATASETS / "HGNC" / "hgnc_complete_set.txt"
GENE_COORDS = DATASETS / "Paralogs" / "ensembl116_gene_coordinates.tsv"   # Ensembl 116 BioMart (chr, start, end)
PARALOGS = DATASETS / "Paralogs" / "ensembl116_human_paralogues.tsv"      # Ensembl 116 paralogues (H1 definition)

PART1 = PROJECT / "external" / "part1_locked"              # locked clone (read-only), commit LOCKED_COMMIT
PART1_MVP_SRC = PART1 / "mvp" / "src"                      # package `coessential`
LOCKED_COMMIT = "c7f95ee007576be9311c0a5f08c6ac63ca23a009"
# Public release of the locked method, https://github.com/cancergenetics/ai-hackathon-Coessentiality. It has a fresh git
# history, so LOCKED_COMMIT is not in it. A clone is accepted only at PUBLIC_MIRROR_COMMIT and only if every Part 1 path
# read by this project has the same git object id as at LOCKED_COMMIT (checked 2026-10-29: all identical).
PUBLIC_MIRROR_URL = "https://github.com/cancergenetics/ai-hackathon-Coessentiality"
PUBLIC_MIRROR_COMMIT = "3c684c972f1e5c3ead722b4b3e2c508fa3cb5ddf"
LOCKED_OBJECTS = {"src": "0cdda7c11991e30b3bbf55b81c8368ba132b75cd",
                  "mvp": "20656be9532ae7f6add07df1dda5d79b3d8392d1",
                  "tests": "83b853f40cc5d6f734e7c6acf43d8862a14df532",
                  "scripts/final_eval.py": "748a70f1a6a31be9d69eb2dafa141962ac10815b",
                  "scripts/predict_pairs.py": "cfed9113c6b905bee44123bf4c9a79e1c40d0d41",
                  "results/ift/look01_2026-10-27.json": "f9062f2ba4f1edce048e7b81659e263f26ef626f",
                  "results/revision/ift_posthoc_spec_k20_tinf_std.json": "5beb584c9d291220ccad2c8c03e96eb212ead389",
                  "results/splits/protein_assignment.sha256": "0bb7670f41bc7ddd474f0d68c97f9027b7fe7d95"}
# Part 1 genome-wide outputs (read in place, never written): genes.tsv, zavg_BEF.npy, manifest.json
PART1_ORIG = Path(os.environ.get("PART1_DIR", PROJECT.parent / "coessentiality_protein complexes"))
PART1_ALL_PAIRS = PART1_ORIG / "outputs" / "all_pairs"
GENES_TSV = PART1_ALL_PAIRS / "genes.tsv"
ZAVG_NPY = PART1_ALL_PAIRS / "zavg_BEF.npy"

CACHE = PROJECT / "data_cache" / "part2"                   # gitignored
TOPK_CACHE = CACHE / "topk"
RESULTS = PROJECT / "results" / "part2"

# ---- Part 2 context inputs (DESIGN v3.1 §2.2; F3 "context inputs": DepMap Omics, Model.csv, Ensembl, HGNC) ----
OMICS_PROFILES = DEPMAP / "OmicsProfiles.csv"
OMICS_DAMAGING = DEPMAP / "OmicsSomaticMutationsMatrixDamaging.csv"
OMICS_HOTSPOT = DEPMAP / "OmicsSomaticMutationsMatrixHotspot.csv"
OMICS_CN = DEPMAP / "OmicsCNGeneWGS.csv"
OMICS_EXPR = DEPMAP / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
OMICS_SIGNATURES = DEPMAP / "OmicsGlobalSignatures.csv"
OMICS_FUSIONS = DEPMAP / "OmicsFusionFiltered.csv"
# CCLE proteomics: used ONLY by the pre-registered LoF functional check of contexts.py (DESIGN §2.2), and later as a
# validation input (tiering on A models). Never a construction/selection input. See tests/test_part2_firewall.py.
CCLE_PROTEOMICS_XLSX = DATASETS / "Proteomics" / "protein_quant_current_normalized.xlsx"
# Lockbox outcome-only input: DEMETER2 (read for L models only through assoc.read_lockbox).
DEMETER2_CSV = DATASETS / "DepMAP_RNAi" / "D2_combined_gene_dep_scores.csv"

CONTEXT_CACHE = CACHE / "contexts"                         # context_table.parquet + genome-wide matrices (gitignored)
CONTEXT_RESULTS = RESULTS / "contexts"                     # committed summaries + checksums
LOCKBOX_DIR = RESULTS / "lockbox"                          # FZ-LB: lockbox_models.txt, A_models.txt, manifest.json
LOCKBOX_MODELS = LOCKBOX_DIR / "lockbox_models.txt"
A_MODELS = LOCKBOX_DIR / "A_models.txt"
LOCKBOX_MANIFEST = LOCKBOX_DIR / "manifest.json"
LOCKBOX_LEDGER = RESULTS / "lockbox_ledger.jsonl"          # git-tracked; one line per lockbox read
TRIALS = RESULTS / "trials.jsonl"                          # every run_association call
REGISTRY_DIR = PROJECT / "docs" / "registry"
