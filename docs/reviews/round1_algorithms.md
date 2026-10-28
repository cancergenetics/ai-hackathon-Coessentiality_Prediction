# Round 1 review: algorithms and compute (role: algorithms)

Scope: DESIGN.md v1 §2.2, §2.6, §2.8, §4. Nothing was run on hu.MAP test pairs, curated complexes or context labels.

**What I ran (all label-free, scratch scripts outside the repo; `zavg_BEF.npy` read-only):**
1. Locked `CoessentialityModel` (read-only import) on an 80% model subsample. A lean variant keeps only per-gene top-50 lists and never writes an n×n matrix. On full data it reproduces the stored matrix bit-exactly (200 random rows: max |Δz| = 0, top-50 sets identical).
2. Parallel timing, 6 jobs (5 bootstraps + full) and 3 jobs.
3. Mutual-kNN graphs from the top-50 lists of full data and 9 bootstraps. Stats: gene essentiality class (fraction of models with effect < −0.5, from CRISPRGeneEffect only) and chromosome arm (HGNC `location`). Leiden-CPM stability vs the 9 bootstraps, with weights w = z − z0 + 1.

## Findings that drive the changes

**Compute (design is feasible, but not as written)**
- One lean bootstrap takes 132 s and 2.0 GB, not 6.6 GB.
  - ARPACK (d = 512) is 90 s (68%).
  - The other stages: init 4 s, F/B row stats 17 s, block scan with top-50 19 s.
- Process parallelism does not help here.
  - 3 jobs with one BLAS thread: 372 s wall = 124 s/job.
  - 6 jobs: 650 s wall = 108 s/job. ARPACK went from 90 s to 300–530 s per job (Accelerate/AMX or memory-bound).
  - The machine was shared with other processes (load average about 10), so I cannot separate the two causes.
- Budget serial time: B = 20 → 45 min, B = 30 → 65 min.
- The mutual-kNN graph is fully determined by per-gene top-K lists. Store 18,429 × 50 (index, z) per run (about 7 MB).

**The graph is mostly proximity artefact among non-essential genes**
- 11,947 of 18,429 genes (65%) are never essential (effect < −0.5 in < 1% of models).
- Mutual-kNN edges (full data):

| z0, k | edges | same-arm edges | both never-essential | either never-essential |
|---|---|---|---|---|
| 5, 20 | 80,301 | 48% | 36% | 63% |
| 6, 20 | 67,555 | 51% | 32% | 59% |

- A random pair is on the same arm 3% of the time.
- Same-arm share rises with z0. For never–never edges it is 70% at z0 = 5 and 81% at z0 = 6. For edges between genes essential in ≥ 5% of models it is 7%.
- Leiden-CPM at z0 = 5, k = 20, γ = 0.5 gives 1,522 clusters of ≥ 3 genes, median size 4. Of these, 46% have a majority of genes on one arm and 76% are majority never-essential.

## Changes

**MUST**

1. **Dependency-bearing gene filter and a proximity flag (§2.2 gene filter).**
   - Problem: "SD below the 5th percentile" removes almost nothing. Never-essential genes fill neighbour slots and form locus-linked clusters that are useless for dependency testing and confound CN contexts (D5).
   - Change:
     - Build the graph only on genes with effect < −0.5 in ≥ 1% of models (≈ 6.5k genes; threshold fixed now, from the matrix only).
     - Take neighbours within this set. Filtering top-50 lists is exact for k ≤ 20, since ≥ 76% of slots for these genes are in-set.
     - Flag any cluster with ≥ 50% of members on one arm (HGNC location is allowed under F3) as "locus-linked". Report it separately and add CN as a covariate if it is tested.
   - Evidence: on the filtered set, stability (mean best-match Jaccard) rises from 0.24–0.60 to 0.45–0.75, and selective-gene coverage stays 0.79–0.99. Same-arm edges fall to 25–29%, so the flag is still needed.

2. **The stability-maximising rule is degenerate (§2.2 parameter selection).**
   - Problem: stability rises monotonically with γ and z0 and falls with k. At z0 = 3, k = 20 on the full graph, it goes from 0.24 to 0.49 for γ from 0.1 to 2. "Maximise stability subject to median size 3–15" therefore selects the smallest clusters (median 4–5), and 4 is inside the constraint. Locus-linked clusters are *more* stable than others (0.47 vs 0.43).
   - Change: do not pick one (z0, k, γ).
     - Use a HiDeF-style persistence hierarchy. Cluster the consensus at a γ grid, keep clusters that recur (best-match Jaccard ≥ 0.75, support ≥ 0.6 of bootstrap × seed runs) over ≥ 3 adjacent γ values, and nest by containment.
     - Choose (z0, k) by split-half recovery of selective-gene coverage, not raw stability.
     - Report stability by size bin and class; give a chance baseline from degree-preserving rewires.
     - Add minimum size 3 and cohesion (mean internal z ≥ z0 + 1) as hard rules.
     - Validate the rule on the toy (item 11): it must recover the planted γ and must not choose dyads.

3. **Bootstrap plan (§2.2, §4).** Run B = 20 to 30 subsamples of 80% of models with the locked class through the lean pipeline above (2 GB, 132 s each, serial, about 45–65 min). Store top-50 lists only. Use B = 10 for parameter selection and the full B once for freezing. Draw subsamples by patient, not model: Model.csv has duplicate patients, so lines from one patient would otherwise straddle subsamples and inflate stability.

**SHOULD**

4. **Consensus without n×n storage.**
   - Compute co-clustering frequency only on the union of bootstrap edges (a few 10⁵ pairs), with the denominator = all B × seeds.
   - Run Leiden-CPM on these frequencies directly with γ = 0.5 (a cluster must average ≥ 50% co-clustering). This replaces "threshold at 0.5, then Leiden".
   - Run it per γ for the persistence step. Cost: 600 Leiden runs (20 bootstraps × 5 seeds × 6 γ) at ≈ 1–1.5 s each, about 10–15 min.

5. **Nesting by construction.** Leiden-CPM at two γ values does not nest. Re-cluster inside each coarse cluster at a finer γ, top-down, on the same consensus frequencies. This is the hu.MAP stage-2 logic.

6. **Edge weights and overlap.**
   - Weights z − z0 + 1 reach about 28 (99.9999th percentile). A few paralog edges dominate CPM. Cap at z0 + 3 or use ranks. Define γ on this scale.
   - The soft-membership rule "mean zavg ≥ 99.9th percentile" is too lax: 18k genes × thousands of clusters gives tens of thousands of false memberships, and one paralog partner can drive a mean. Replace it with "≥ 50% of members have z ≥ z0 to g" plus a size-matched null at FDR 0.01 across all gene × cluster tests.

7. **BioPlex (not yet in `datasets/`).**
   - Keep co-essentiality-only complexes (C1) as primary for D1–D6. Treat the BioPlex-integrated set (C2) as a sensitivity analysis and for annotation. This avoids coupling findings to 293T and HCT116 detectability.
   - Held-out evaluation:
     - Hold out baits, not just cell lines. 293T → HCT116-only edges leak through triadic closure, and HCT116-only edges are partly 293T false negatives. Use K-fold over bait proteins within each line, plus the cross-line check.
     - The null must match BioPlex degree, essentiality class, expression and cluster size, otherwise detectability inflates enrichment.
   - Variants are chosen on stability and split-half recovery; held-out BioPlex enrichment is reported once, for the final choice.
   - Fix z1 now.

8. **Cross-fit construction (§2.4).** Clusters built on all models make the validation half non-independent. Both halves informed membership, so replication rates inflate for lineage-driven modules. Build complexes on half A and test them in half B (and swap), using lineage-stratified halves with the same pipeline (2 extra half-size runs per split). The headline "replicated" status uses cross-fit complexes. The permutation p-value on the full-data complexes remains valid, because clusters depend on X, not on labels.

9. **Alternatives, within the environment.**
   - Run Infomap (igraph, Boyle precedent) and a sparse-matrix MCL with an inflation sweep as sensitivity checks. Report ARI/overlapping NMI to the consensus, not "better".
   - ClusterONE and HiDeF are not installed. Ask the PI only if wanted. The persistence step in item 2 is about 80 lines in numpy.

## Toy simulation design (§2.8)

Run the real `CoessentialityModel` on ≈ 3,000 genes × 300 lines (seconds). Plant:
- complexes of size 3–40, nested sub-modules and one shared subunit;
- pan-essential cores;
- 60% never-essential genes with arm-block correlated noise (CN-like) and a few hubs;
- lineage low-rank structure and library batch.

Checks:
- Recovery (ARI, best-match F1) by size and loading.
- The item-2 rule recovers the planted γ, avoids dyads and does not reward arm blocks (arm blocks only → no non-locus clusters; no planted complexes → zero persistent clusters).
- Planted context effects (prevalence 5–30%, 8–200 models): calibrated permutation p (QQ) and power.
- A lineage-driven planted module shows replication inflation without cross-fit (item 8).
- LOSO flags a one-subunit effect.
