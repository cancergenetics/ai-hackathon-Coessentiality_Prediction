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

PART1 = PROJECT / "external" / "part1_locked"              # locked clone (read-only), commit LOCKED_COMMIT
PART1_MVP_SRC = PART1 / "mvp" / "src"                      # package `coessential`
LOCKED_COMMIT = "c7f95ee007576be9311c0a5f08c6ac63ca23a009"
# Part 1 genome-wide outputs (read in place, never written): genes.tsv, zavg_BEF.npy, manifest.json
PART1_ORIG = Path(os.environ.get("PART1_DIR", PROJECT.parent / "coessentiality_protein complexes"))
PART1_ALL_PAIRS = PART1_ORIG / "outputs" / "all_pairs"
GENES_TSV = PART1_ALL_PAIRS / "genes.tsv"
ZAVG_NPY = PART1_ALL_PAIRS / "zavg_BEF.npy"

CACHE = PROJECT / "data_cache" / "part2"                   # gitignored
TOPK_CACHE = CACHE / "topk"
RESULTS = PROJECT / "results" / "part2"
