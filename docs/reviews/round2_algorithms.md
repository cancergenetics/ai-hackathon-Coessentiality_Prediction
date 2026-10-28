# Round 2 review: algorithms and compute (role: algorithms)

Scope: DESIGN v2 §2.3, §2.4 (±10 Mb rule), §2.7 (class-C controls), §5. No hu.MAP test file, curated complex or context association was touched.

**What I ran (label-free; scratch scripts outside the repo):**
- Filter and graph on the existing all-model `zavg_BEF` matrix, restricted to the v2 gene filter.
- Leiden-CPM timing and cohesion over a γ grid.
- BioPlex file hygiene.

**Round-1 status.** My MUSTs 1–3 (dependency filter, persistence instead of stability maximisation, bootstrap plan) are accepted and correctly summarised. SHOULDs 4–9 are accepted. The items below are what an implementer would otherwise have to guess, plus four new defects.

## MUST

1. **The locked predictor will refuse to run on A or a subsample.**
   - Problem: `predict.py` asserts that the null moments equal the full-data manifest, and the null pairs are tied to the full gene set. On A or an 80% subsample the moments differ, so the run aborts or, if the assertion is bypassed, uses full-data moments silently.
   - Change: state the rule.
     - Use the same 200k null pairs (SEEDS[0]).
     - Recompute the moments on each run's own data, with `check_null=False`.
     - Log the moments per run.
     - Fix the gene universe to A, so a gene that is low-coverage in a subsample gets NaN scores in that run only.
     - Seed ARPACK `v0` per run; PYTHONHASHSEED does not control it.
     - `reproduce.sh` compares edge lists and clusters, not matrix hashes. BLAS threading is not bitwise-stable across machines.
   - Why: z0 ∈ {4,5,6} must mean the same quantity in every run.

2. **Store the S×S score block per run, not top-50 lists.**
   - Problem: v2 says neighbours within the filter set are "exact for k ≤ 20" from top-50 lists. Measured on the full matrix, the filter (≥ 1% effect < −0.5, NaN ≤ 5%, not low-coverage) gives |S| = 6,163.
     - 26% of S genes have < 20 in-set neighbours in their top-50.
     - 4% have < 10 in-set neighbours in their top-50.
     - Edges lost by truncation:
       - k=20, z0=4: 191 of 37,928.
       - k=20, z0=5: 29 of 34,638.
       - k=10, z0=4: 12 of 17,245.
   - Change: the lean run writes the S×S float32 block (≈ 150 MB; 20 runs ≈ 3 GB, no n×n). Then:
     - top-k, mutual edges and the cohesion rule are exact;
     - the whole (z0, k) grid costs nothing extra;
     - GLS uses the same block format.
   - Why: this removes the "exact" claim and the need for a separate pair-score store.

3. **The clustering pipeline has two unspecified γ scales and undefined persistence terms.** Fix all of the following in FZ-C:
   - **Stage 1, per run and seed:** Leiden-CPM on edge weights in [1,4]. Use 5 seeds (v2 says both "SEEDS 0–4" and "Leiden seeds 0–19"; use SEEDS).
   - **γ grid:** one geometric grid, γ_j = 0.05·2^(j/2), j = 0…10 (0.05…1.6).
     - On the real filtered graph this spans median size 24 down to 7.
     - Coarse γ gives clusters of > 60 (18 at 0.05, 3 at 0.1), which are re-clustered.
     - Three adjacent √2 steps span a 2× range, which is a sensible "persistence".
   - **Stage 2:** the pair set is the union of pairs co-clustered in at least one run. Frequencies use denominator = all runs, and a pair absent from a run (gene NaN there) counts as 0. Run Leiden-CPM on the frequency graph at the same γ index, but with γ_2 = 0.5 fixed (cluster = average co-clustering ≥ 0.5). Then the γ-dependence lives in stage 1 only.
   - **Persistence:** cluster c at γ_j recurs at γ_{j+1} if the best-match Jaccard is ≥ 0.75. "Support ≥ 0.6" means the fraction of stage-1 runs (at that γ) in which a ≥ 0.75-Jaccard match to c exists.
   - **Nesting:** child ⊂ parent if |child ∩ parent| / |child| ≥ 0.8. A gene may sit in one cluster per level.
   - **Cohesion:** mean of the A-score over **all** within-cluster pairs, computed from the full-A S×S block. Measured pass rates for pair-mean ≥ z0+1 at k=20, z0=5:
     - 23% at γ=0.05;
     - 38% at γ=0.1;
     - 78% at γ=0.4;
     - 100% at γ=1.6.
     - So the rule silently removes coarse clusters. Apply it after persistence and report how many it removes.
   - Leiden runtime is 0.2–0.3 s on this graph (≈ 35k edges). The full grid is 6 combos × 11 γ × 100 runs ≈ 30 min single-core; it parallelises across processes.

4. **Define the split-half selection statistic.** "Recovery of selective-gene coverage" is undefined, and the obvious version favours the loosest graph (z0=4, k=20: 6,026 of 6,038 nodes in one giant component).
   - Change:
     - Use R = 4 PatientID-blocked, lineage-stratified splits. Each half is one full run with 5 seeds, not a bootstrap consensus: 8 runs, about 15–20 min.
     - Run the same persistence step on each half.
     - A half-1 cluster of size 5–30 counts as recovered if some half-2 cluster has Jaccard ≥ 0.5.
     - Criterion = the fraction of selective genes (5–90% of models below −0.5, computed on A) lying in recovered clusters, averaged over both directions and the R splits.
     - Report this next to the same statistic on degree-rewired graphs.
     - Ties go to the smaller graph (higher z0, then lower k).
   - Why: size 5–30 stops blobs and trios from counting.

5. **The GLS bridge "matched by edge count" is under-determined.**
   - Problem: the GLS scale is not z, so the weight cap z0+3 is meaningless for it. The mutual top-k cap may bind before the threshold does, so the counts may not be matchable. "Same rule applied to GLS" (re-selecting (z0,k)) breaks the match.
   - Change:
     - Use the same gene set S and the same k as the selected zavg graph.
     - Choose the GLS threshold that reproduces the zavg edge count (±1%); if the mutual-top-k cap makes that impossible, report the largest attainable count.
     - Map GLS edge weights to [1,4] by quantile-matching to the zavg weight distribution.
     - Run the bootstrap and persistence pipeline unchanged.
     - Tuning GLS's own (z0, k) is a secondary variant.
     - Compare the bridge metrics at matched number of clustered genes and matched cluster-size distribution. Otherwise T1/T2 fractions and context-effect correlations differ by coverage, not score quality.
   - Cost: 24 s per GLS run.

6. **The "within 10 Mb" and "±10 Mb" rules cannot be implemented with the files provided.**
   - HGNC has only cytoband (`location`, e.g. 19q13.43). The Ensembl paralog file has no coordinates. The CN file's columns are in Ensembl-ID order, not genomic order.
   - Change: add one Ensembl 116 gene-coordinates table (chr, start; BioMart) to DATA_REQUEST. Until then, define the locus flag by arm plus the CN-segment rule (|r| ≥ 0.9 of CN across A), and say that ±10 Mb in §2.4 is approximated by same-arm.
   - Why: this is the confound protection for D5/CN contexts.

## SHOULD

7. **BioPlex T1 is now implementable. Specify the following.**
   - Data: the files are pre-thresholded at pInt ≥ 0.75, with no self-loops or duplicates.
     - 293T: 118,162 edges, 8,995 baits.
     - HCT116: 70,966 edges, 4,866 baits.
     - Union: 167,374 edges; 21,754 shared.
   - Gene IDs: 24% of UniProt IDs are isoforms, so map by Entrez ID (GeneA/GeneB) to the CRISPR column's Entrez, never by symbol.
   - Per cluster: restrict to members that are BioPlex nodes (no edges can mean "not tested" or "not detected"). Require ≥ 3 such members and ≥ 2 observed edges.
   - Null: draw random sets from the BioPlex-node ∩ S universe, with each member replaced by a gene in the same BioPlex-degree bin (log-spaced) and the same essentiality class, 10⁴ draws.
   - Multiple testing: p < 0.05 uncorrected will assign T1 to about 5% of null clusters. Use BH q < 0.1 over clusters and report the expected false-T1 count.
   - Bait K-fold applies only to the integrated map (BioPlex edges used in construction). Assign each gene a fold and hold out every edge touching a fold gene. T1 for the label-free map needs no holdout.

8. **The negative-only filter drops growth-suppressor members.** TP53, CDKN1A, USP28, TBC1D7 and CDKN2A have effect < −0.5 in < 1% of models, so they are outside S. PC3 class C is then "not evaluable", and the TSC complex falls to a dyad (TSC1/TSC2 are in S, but TBC1D7 is not). Either add genes with effect > +0.5 in ≥ 1% of A (743 genes on the full matrix) or declare PC3-C non-evaluable now.

9. **Chance baseline for stability.**
   - Rewire the union graph once (double-edge swaps) and carry each edge's per-run presence vector with it. The rewired copy keeps degree and edge-noise structure, so stability is not trivially zero.
   - Soft-overlap null: use the closed form. For gene g, the number of members with score ≥ z0 is hypergeometric(|S|, deg_g, |C|), and the membership test is a tail at ⌈|C|/2⌉. This avoids permutation over millions of gene×cluster tests.
   - Infomap/MCL: compare to the persistence level with the nearest cluster count; report ARI and best-match F1.

10. **Update §5 and the seeds.**
    - Total serial compute:
      - 20 A-bootstraps and 1 full-A run: ≈ 46 min.
      - 8 split-half runs: ≈ 18 min.
      - GLS, 21 runs: ≈ 8 min.
      - Clustering grid: ≈ 30 min.
      - Total: ≈ 1 h 45 min. v2's "30 runs ≈ 65 min" omits the halves and GLS.
    - Do not run ARPACK jobs in parallel (round 1: 90 s → 300–530 s per job under contention). Leiden and null draws can use 8 processes.

## COULD
- `OmicsCNGeneWGS` can define the locus segment data-driven, which also handles the 10 Mb problem where gene coordinates are missing.
- Keep the all-model `zavg_BEF.npy` out of Part 2 (it includes L). The scripts check the manifest.

**Verdict.** Items 1–6 are needed before FZ-C; none changes the design's scope. With them, the construction is correct, runs in about 2 h serial on the laptop, and has no remaining guess points that I can see.
