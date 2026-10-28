# Part 2 construction infrastructure (status 2026-10-28)

This folder covers infrastructure only. No lockbox has been frozen, no construction set A exists, and no context labels, hu.MAP files, BioPlex, proteomics or HGNC `gene_group` were read. The specification is DESIGN.md v3.1 §0, §2.1 and §2.3, plus docs/reviews/round2_algorithms.md (its defaults are adopted, and every parameter is configurable).

## What exists

| Module (`src/part2/`, import as `part2`) | Content |
|---|---|
| `paths.py`, `part1.py` | File locations. Read-only import of the locked Part 1 code: `coessential` (CoessentialityModel, `Config()`) and `src.baselines` (GLS). `assert_locked()` checks that the clone is at c7f95ee and clean. |
| `data.py` | Loads CRISPRGeneEffect exactly as Part 1 does (pandas → float32), with an .npy cache keyed by the CSV's SHA-256. Also loads Model.csv and CRISPRConfounders. `universe_models()` returns U = 1,361 − 16 engineered/drug-adapted = **1,345** models. `GeneEffect.subset(ids)` returns rows in CSV order, so a subset is a function of the *set* of IDs. `model_set_hash`, `blocked_subsample` (PatientID blocks) and `blocked_split` (PatientID blocks, stratified e.g. by lineage). Gene keys come from Part 1's genes.tsv. |
| `scores.py` | Lean locked recompute. `score_block(model_ids, method, S)` returns the **S×S float32 block** (per round-2 MUST 2) plus, for zavg, the per-gene top-50 lists over all protein genes. `scores_topk(...)` gives exact top-K within S. Results are cached in `data_cache/part2/blocks/<method>_<model-hash>_S<S-hash>.npz`. The JSON sidecar records model IDs and hash, runtime, stage timings, peak RSS, code commit and the **null moments recomputed on the subset**. Null moments use the locked procedure: same Config, n_null = 200k, seed 0, seeded ARPACK v0. S genes that are not background genes in a run get NaN for that run only. `gls` is Part 1's dev-selected `gls_signed` (exact pairwise-complete per missingness pattern, whitening 'drop'). |
| `genefilter.py` | S = effect < −0.5 in ≥1% of models, non-NaN ≥95%, not Part 1 low-coverage. An optional +0.5 arm exists (round-2 SHOULD 8) and is off by default. Also `selective_mask` (5–90% of models). |
| `graph.py` | Mutual top-k ∧ score ≥ z0; weight min(s, z0+3) − z0 + 1. `matched_graph` is the GLS bridge: same S and k, matched edge count, attainability flag, quantile-matched weights. Also a degree-preserving double-edge-swap rewire that carries each edge's per-run attributes. |
| `consensus.py` | Stage 1: Leiden-CPM per run × SEEDS (0–4) over γ = 0.05·√2^j, j = 0…10. Stage 2: union of co-clustered pairs, frequency over all runs, CPM γ₂ = 0.5. Persistence: Jaccard ≥ 0.75 to the next γ, support ≥ 0.6, ≥3 adjacent γ. Coarse-first, one cluster per gene. Size 3–60 and cohesion: all-pairs mean ≥ z0 for coarse clusters, ≥ z0+1 for sub-modules. Violators are re-clustered top-down on the same frequencies; sub-modules are re-clustered inside each coarse cluster. Rule counts are returned. `ensure_hashseed()` re-execs scripts with PYTHONHASHSEED=0. |
| `selection.py` | Split-half statistic: R splits; each half is one run × 5 seeds; half-1 clusters of size 5–30 recovered at Jaccard ≥ 0.5; selective-gene coverage averaged over both directions; locus clusters excluded. It also reports the locus-included statistic and a rewired baseline. Tie-break goes to higher z0, then lower k. Also stability by size bin vs a rewired baseline. |
| `locus.py` | Locus flag = ≥50% of members within 10 Mb (Ensembl 116 coordinates, Entrez join, primary chromosomes, midpoints; 99.9% of genes positioned) **or** on one arm (HGNC). The HGNC band is reported for interpretation. |
| `toy.py`, `pipeline.py` | Toy simulator (see below) and the glue shared by toy and DepMap runs. |

Scripts:
- `scripts/30_toy_construction.py` writes `toy_results.json` (21 min).
- `scripts/31_score_blocks.py --set all|U|sub --method zavg_BEF|gls` caches a block and appends to `timings.jsonl`.

## Tests (`tests/test_part2_*.py`; 29 pass in 69 s; `.venv/bin/python -m pytest -q tests/test_part2_*.py`)

- **Bit-exact:** the all-model (1,361) zavg block over S (6,163 genes) equals Part 1's stored `zavg_BEF.npy` on 200 random rows × all S columns (Δ = 0). The all-gene top-50 lists are identical. Recomputed null moments equal Part 1's manifest (|Δ| ≤ 1e-12; in fact equal).
- **GLS:** compared with Part 1's `src/predict.py` `BaselineEngine.compute('gls')` on 4,999 random S pairs. The stored float32 equals float32(Part 1) for all pairs (max |Δ| 2.3e-7 = float32 rounding). On toy data with NaN patterns it matches `baselines.gls_pairs`.
- **Data:** U = 1,345; gene keys match genes.tsv and `protein_genes()`; subset is order-invariant; no patient straddles a subsample or split; |S| on all models = 6,163 (matches the round-2 review).
- **Graph:** mutual-kNN, thresholds, matched counts, rewire preserves degree, locus arm and window rules, coordinate join.
- **Consensus:** hand-computed co-clustering; planted cliques and nesting recovered; null graph gives 0 clusters; size and cohesion rules; **determinism in fresh subprocesses** (PYTHONHASHSEED 0, 0 and 12345 give identical maps).
- **Toy:** planted ARI ≥ 0.8; null gives 0 clusters; arm-block clusters are all locus-flagged and score 0 on split-half; a dyad map scores 0.
- **Firewall:** grep finds no huMAP / test_eval / BioPlex / proteomics / gene_group / curated-complex references in src/part2 or scripts/3x.

## Toy results (`toy_results.json`)

Setup: 300 lines × 3,000 genes, Config(k=10, d=64), B = 8 bootstraps of 80%, 5 seeds, z0 = 5, k = 10.

**Planted scenario (14 complexes, sizes 3–40):**
- ARI 0.875.
- Best-match F1 = 1.0 for 8 of 14 complexes and ≥ 0.75 for all complexes of size ≥ 5.
- Both planted sub-modules are recovered (F1 = 1).
- **Misses:** complexes of size 3 or 4 (3 of 4 missed).
- The nested 40-gene parent is not recovered as one coarse cluster: its cross-sub-module mean z is 4.6. Its two halves are found as separate coarse clusters.

**Split-half selection:**
- The criterion picks (6, 20). Its ARI of 0.878 is the best in the grid; the grid spans 0.870–0.878.
- The rewired baseline scores 0 for every configuration.

**Null scenario:** 0 clusters, and the split-half statistic is 0.

**Arm-block scenario:**
- All 6 clusters are locus-flagged.
- The split-half statistic is 0 with the locus exclusion. Without the exclusion it is 0.11–0.19, so the flag is what stops arm blocks being rewarded.

## Timings (serial, laptop, load average 2.4–4.2, toy job running on 1 core)

| Run | Models | \|S\| | Runtime | Peak RSS | ARPACK |
|---|---|---|---|---|---|
| zavg_BEF, all models | 1,361 | 6,163 | 109 s | 3.09 GB | 77 s |
| zavg_BEF, U | 1,345 | 6,132 | 116 s | 3.18 GB | 83 s |
| zavg_BEF, U 80% PatientID-blocked (seed 1000) | 1,076 | 6,132 | 120 s | 3.05 GB | 88 s |
| GLS, U | 1,345 | 6,132 | 22 s | 2.06 GB | – |

- Each block is 150 MB.
- Indicative clustering time on the U block (1 run × 5 seeds): stage 1 over the 11-value γ grid takes 21 s; persistence and nesting take 8 s.
- Estimate for the frozen map: 20 bootstraps × 5 seeds ≈ 7 min of stage 1 per (z0, k).
- The real-scale smoke run is timing only. Its cluster list was not saved, because U includes the future lockbox.

## Open issues for the design (not decided here)

1. **"Support" as specified drops small complexes.** With support defined as "a stage-1 run holds a ≥ 0.75-Jaccard match", 3–4-gene complexes are consistently co-clustered at coarse γ but sit inside larger, varying stage-1 clusters. They reach support ≥ 0.6 only at the finest γ, so they fail the 3-adjacent-γ rule. Two alternatives fix the toy misses:
   - `support_mode="containment"` (≥ 75% of members in one run cluster). Cost: 1 false size-3 cluster in the null and 1 non-locus cluster in the arm-block scenario.
   - Extending the grid to 13 steps (to 3.2). This produced no false clusters.
2. **Rewired stability is high for small clusters.** The rewired chance baseline has stability 0.84–0.90 for size 3–4 clusters, because per-run edge presence is carried with each edge. Small-cluster stability is therefore not informative.
3. **The split-half criterion is nearly flat.** On the toy it varies little (0.137–0.167). On real half-data the cohesion rule rejects many candidates, because half-sample z-scores are lower; for example, 205 rejections in a toy half.
4. **Greedy one-per-gene placement drops many candidates.** On the U smoke run, 381 of 567 persistent candidates were dropped by trimming.
5. **Locus flags are frequent.** 30% of clusters in the U smoke run are locus-flagged (arm OR 10 Mb).
6. **"Reproduces Part 1 scores on all U models" cannot hold literally.** Part 1 was computed on all 1,361 models, so the reproduction test runs on all 1,361 models. U runs use the same code, with null moments recorded in the sidecar.
7. **Not implemented:** the attachment rule (DESIGN v3.1), Infomap/MCL, BioPlex tiers and the full bootstrap set (waiting for FZ-LB).
