# FZ-C pre-freeze algorithm decisions

Logged **before** the zavg / GLS maps were built and committed, and before any physical-interaction or
protein-abundance file was opened. Every decision is label-free: it was settled on the toy simulation (planted truth,
`scripts/35_toy_decisions.py` → `toy_decisions.json`) and, where noted, on label-free accounting on A's own score
blocks. None of them was revisited after the validation step (tiers / bridge).

The spec is DESIGN.md v3.1 §2.3 together with docs/reviews/round2_algorithms.md. The open issues are those listed in
results/part2/infra/README.md.

## D1 Persistence support: containment (infra issue 1)

- **Rule.** A stage-2 cluster c at γ_j counts as supported in a stage-1 run when one run cluster holds ≥ 75% of c's
  members (`support_mode="containment"`). The round-2 wording asked for a ≥ 0.75-Jaccard match in that run. The other
  persistence terms are unchanged: Jaccard ≥ 0.75 recurrence at γ_{j+1}, support ≥ 0.6, and ≥ 3 adjacent γ.
- **Why.** Small complexes of 3–4 genes are co-clustered at coarse γ, but in each run they sit inside a larger,
  varying stage-1 cluster. A Jaccard-based support therefore fails them at every γ except the finest, and the
  3-adjacent-γ rule then removes them.
- **Toy evidence** (3 toy seeds, 20 subsamples, z0 = 5, k = 10): see the table below.
  - Containment raises recovery of the planted complexes (best-match F1 ≥ 0.8 at any level), mainly the 3–4-gene ones.
  - Its price is false clusters, counted three ways: clusters on the null scenario; clusters with < 50% of members
    from any one planted complex; and non-locus clusters in the arm-block scenario.
- **Label-free check on A** (6 subsamples, z0 = 5, k = 20; exploratory, not saved). Level-1 gene coverage was 2,960
  genes with Jaccard support and 4,180 with containment, from 231 and 439 persistent candidates respectively.

## D2 Placement: persistence-first, containment = nesting (infra issue 4)

- **Rule** (`placement="nested"`). Candidates are visited most persistent first: adjacent-γ count, then support,
  then size, then coarser γ. Each candidate is compared with the clusters already accepted:
  - **Duplicate** (Jaccard ≥ 0.75 with an accepted cluster): dropped.
  - **Nested child** (≥ 80% of it inside one accepted level-1 cluster): kept as a sub-module candidate of that cluster.
  - **Parent** (it holds ≥ 80% of one or more accepted clusters): it becomes the level-1 cluster, with members equal
    to the candidate ∪ those clusters, and they become its sub-modules. This applies only if the merged parent
    passes size ≤ 60 and the level-1 cohesion rule (mean ≥ z0); otherwise the candidate is dropped.
  - **Partial overlap:**
    - If more than 50% of the candidate lies in other accepted clusters, the candidate is dropped and the overlap is
      resolved to the more persistent cluster.
    - Otherwise the overlapping genes are trimmed from it.
- **Unchanged from the spec:**
  - Size and cohesion violators are re-clustered on the same frequencies.
  - Level 2 (sub-modules) is still the top-down re-clustering inside each level-1 cluster at finer γ (stage-2
    γ₂ = 0.5). Its candidates are pooled with the nested children found above and placed by the same rule (one
    partition per parent, sub-module cohesion mean ≥ z0 + 1, no third level).
- **Toy evidence:** see the table below.
  - On every toy scenario the greedy and nested placements recover the same planted complexes and sub-modules. On
    the toy, greedy's dropped candidates were nested children, and greedy's level-2 re-clustering found them again.
  - The extra hierarchical scenario (`hier`) shows the same: three cohesive 24-gene parents with 3 × 8 sub-modules,
    and two pairs of 10-gene complexes sharing 3 subunits.
- **Label-free accounting on A** (6 subsamples, z0 = 5, k = 20, containment support). This is what infra issue 4
  measured. Under greedy placement 201 of 439 candidates were "overlap-dropped". Under the nested rule the same
  candidates are accounted for as follows:

  | Outcome | Candidates |
  |---|---|
  | Nested children (giving 92 sub-modules) | 65 |
  | Children absorbed into a more inclusive parent | 45 |
  | Trimmed | 28 |
  | Overlap conflicts | 8 |
  | Parent-rule failures | 17 |
  | Duplicates | 5 |

  Level-1 gene coverage is the same under both rules: 4,180 (greedy) and 4,192 (nested).
- **Conclusion.** "Greedy loses most candidates" was mostly an accounting artefact: the dropped candidates were
  nested children. The nested rule is adopted because it keeps those children as sub-modules directly and resolves
  conflicts by persistence rather than by γ order. It does not change recovery on the toy.

## D3 Attachment rule for genes outside S (DESIGN §2.3; infra issue 7)

- **Rule.** Gene g (a protein gene outside S that has a score in A) is attached to cluster C (any level) if both hold:
  - x = #members with score ≥ z0 to g is ≥ ⌈|C|/2⌉;
  - the size-matched null gives BH q < 0.01 over all gene × cluster tests.
- **Null.** Closed form (round-2 SHOULD 9): x ~ Hypergeometric(N_g, K_g, |C|), where N_g is the number of S genes
  with a finite score to g and K_g the number of those with score ≥ z0. Pairs with x = 0 have p = 1 and stay in the
  BH denominator.
- **Scores.** From the full-A zavg_BEF run, which also stores the (non-S) × S rectangle (`scores.zavg_block(rect=True)`;
  tested equal to the full block).
- **Status.** Attached genes are not members. They may define H1 events.
- **Locus guard.** The arm-block toy attaches never-essential genes to arm-block clusters (CN artefact, table below).
  An attachment is therefore flagged `locus_attach` when any of these holds, and such attachments are not used to
  define H1 events:
  - the cluster is a locus cluster;
  - the gene sits on the cluster's dominant arm;
  - the gene lies within 10 Mb of a member.
- **Toy evidence:**
  - Before clustering, one member of every planted complex of size ≥ 5 is held out of S. These held-out members must
    attach to their own complex.
  - Never-, pan- and selective-essential genes must not attach.
  - Random size-matched member sets on the null scenario must collect no attachments.
  - Results are in the table below.

## D4 GLS comparator scale (round-2 MUST 5 detail)

- **Graph.** The GLS graph of each run (each subsample or half) uses the same S and k as the selected zavg graph.
  It keeps the top mutual-top-k GLS pairs up to the edge count of the zavg graph built on the *same model set*, with
  weights quantile-matched to that zavg graph's weights.
- **Cohesion.** The cohesion rule ("mean ≥ z0", "≥ z0 + 1") is applied to the GLS block rank-mapped onto the zavg
  block of the same model set (`construct.rank_map`). The rule therefore removes the same quantile in both maps.
- **Selection.** GLS's own split-half grid (each grid point edge-matched to zavg's graph at that (z0, k)) is reported.
  The primary GLS map uses zavg's selected (z0, k), so the bridge compares maps at matched edge count.

## D5 Bridge details fixed before validation (DESIGN §2.3, round-3 B1)

- **Clusters and binning:**
  - Unit: level-1 clusters (sub-modules are nested and would double-count pairs). All level-1 clusters are used,
    including locus clusters.
  - Size bins: {3–4, 5–8, 9–16, 17–60} by the cluster's member count.
  - Precision within a bin is pair-pooled: BioPlex edges / pairs, over pairs with both proteins BioPlex-observed.
    BioPlex-observed means an Entrez node of the 293T ∪ HCT116 union.
  - Bin weights: the share of pooled zavg + GLS level-1 clusters (after coverage matching) in each bin. A bin
    without pairs in a map is dropped and the remaining weights are renormalised.
- **Matched gene coverage.** The map that clusters more genes at level 1 drops its least stable level-1 clusters
  (ascending stability, then cohesion, then descending id) until its coverage is ≤ the other map's.
- **Bootstrap and decision:**
  - CI: protein-node bootstrap, B = 2000, `default_rng(5000)`. Nodes are the BioPlex-observed proteins in either
    map's within-cluster pairs, resampled multinomially. Clusters are fixed, and each pair is weighted by the
    product of its two nodes' multiplicities. The paired difference is zavg − GLS.
  - Decision: "better" iff the 95% percentile CI excludes 0.
- **Sensitivities and secondary endpoints:**
  - Sensitivities: no coverage matching; locus clusters excluded.
  - Secondary endpoints:
    - protein co-abundance with the same standardisation (pairwise-complete Pearson over A models with ≥ 20
      measured models per protein; bootstrap seed 5001);
    - stability by size bin;
    - pair-level AP / AUROC against BioPlex over all S × S pairs with both genes observed, for zavg_BEF, GLS and
      Pearson on A.
- **Tiers** (DESIGN §2.3):
  - **T1:**
    - Null: each BioPlex-observed member is replaced by a random BioPlex-observed S gene of the same log2-degree
      bin (degree in the union network) and A-essentiality class (pan ≥ 90% / selective 5–90% / other). Sparse
      strata are merged to ≥ 10 genes. 10⁴ draws, seed 4100 + level.
    - Testable when ≥ 3 members are observed. T1 needs ≥ 2 within-cluster edges and BH q < 0.05 within the level.
  - **T2:**
    - Statistic: mean pairwise protein correlation among ≥ 3 measured members (CCLE proteomics, A models only,
      mapped by gene symbol).
    - Null: random measured S genes of the same count, 10⁴ draws, seed 4200 + level. BH q < 0.05 within the level.
    - Assigned only to clusters that are not T1.
  - **T3:** T3-neg if testable by either source, T3-nc otherwise.

## D6 Split-half selection (DESIGN §2.3; infra issue 3)

- **Statistic.** Unchanged from round-2 MUST 4: R = 4 PatientID-blocked, lineage-stratified 50/50 splits of A
  (seeds 2000–2003). Each half is one run × 5 Leiden seeds through the same pipeline (D1, D2). Half-1 clusters of
  size 5–30 that are not locus clusters count as recovered when a half-2 cluster has Jaccard ≥ 0.5. The statistic is
  the selective-gene coverage of recovered clusters, averaged over both directions and the 4 splits.
- **Cohesion.** Applied on each half's own block.
- **Ties.** Broken by higher z0, then lower k.
- **Baseline.** The degree-rewired baseline (seed 4000 + 2·split) is reported next to the statistic.

## Toy results (filled from `toy_decisions.json`)


**Setup** (`toy_decisions.json`; runtime 2,529 s):
- 3 toy seeds × 4 scenarios.
- 20 row subsamples of 80% (seeds 1000–1019) × Leiden SEEDS 0–4.
- Config(k = 10, d = 64), z0 = 5, k = 10.
- Planted truth: 14 complexes per seed (sizes 3–40) plus, in `hier`, 3 parents with 3 × 8 sub-modules and 2
  overlapping pairs. That makes 35 truth sets per seed, 105 over the 3 seeds.

**Reading:**

- **D1 (containment).**
  - Recovery of planted complexes rises from 67/105 to 85–86/105. For 3–4-gene complexes it rises from 6/30 to
    24–25/30.
  - Null scenario: 0 clusters in all 3 seeds. With only 8 subsamples, one 3-gene null cluster had appeared (seed 0;
    infra README). With the 20 subsamples used for FZ-C it does not.
  - Arm-block scenario: 1 non-locus cluster of 25. All the rest are locus-flagged, and the locus flag keeps them out
    of claims and out of the split-half statistic.
  - The 6 "false clusters" (< 50% of members from one planted complex) are not noise. Inspection of seed 0 shows
    coarse level-1 merges of two small planted complexes, bridged by the simulated hub genes (each hub loads on
    three complexes). The constituent complexes are recovered as sub-modules inside them. Jaccard support avoids
    these merges only by losing the small complexes altogether.
- **D2 (nested).**
  - Recovery equals greedy's (86 vs 85/105). Sub-modules: 39/39 for both. Overlapping complexes: 12/12 for both.
  - Parents: only 2 of the 9 planted `hier` parents reach level 1 at F1 ≥ 0.8. The others fail the level-1 cohesion
    rule (cross-sub-module mean z < z0), which is the specified behaviour, so their sub-modules appear as level-1
    clusters.
- **D3 (attachment).**
  - Held-out members: all 81 attach to their own complex (10 + 17 per seed). There are 3 extra wrong-cluster pairs
    (seeds 0–1, planted).
  - Spurious attachments in planted / hier: 6 never-essential genes over all seeds.
  - Null scenario: random size-matched sets collect 0 attachments in 75k tests.
  - Arm-block: 92–116 never-essential genes per seed attach to arm-block clusters. These are all CN-arm artefacts on
    locus-flagged clusters, hence the locus guard in D3.

| Variant | Planted complexes recovered (F1 >= 0.8, any level; planted + hier) | of which size 3-4 | Sub-modules recovered | hier parents at level 1 | hier overlapping complexes | False clusters (< 50% from one complex) | Null-scenario clusters | Arm-block non-locus clusters |
|---|---|---|---|---|---|---|---|---|
| V0_jaccard_greedy | 67/105 | 6/30 | 39/39 | 1/9 | 12/12 | 0 | 0 | 0 of 19 |
| V1_containment_greedy | 85/105 | 24/30 | 39/39 | 2/9 | 12/12 | 6 | 0 | 1 of 25 |
| V2_containment_nested | 86/105 | 25/30 | 39/39 | 2/9 | 12/12 | 6 | 0 | 1 of 25 |
| V3_jaccard_nested | 67/105 | 6/30 | 39/39 | 1/9 | 12/12 | 0 | 0 | 0 of 19 |

Attachment (D3, V2 map):

| Toy seed | Scenario | Held-out members | Attached to own complex | Wrong-cluster pairs | Other genes attached (by class) | Tests |
|---|---|---|---|---|---|---|
| 0 | planted | 10 | 10 | 1 | {'never': 1} | 30702 |
| 0 | hier | 17 | 17 | 0 | {'never': 2} | 61642 |
| 0 | null, 14 random size-matched sets | 0 | - | - | 0 attached | 25144 |
| 0 | armblock | 0 | 0 | 0 | {'never': 92} | 12572 |
| 1 | planted | 10 | 10 | 2 | {'never': 1} | 30668 |
| 1 | hier | 17 | 17 | 0 | {'never': 2} | 57952 |
| 1 | null, 14 random size-matched sets | 0 | - | - | 0 attached | 25116 |
| 1 | armblock | 0 | 0 | 0 | {'never': 114} | 14352 |
| 2 | planted | 10 | 10 | 0 | {} | 28896 |
| 2 | hier | 17 | 17 | 0 | {} | 58016 |
| 2 | null, 14 random size-matched sets | 0 | - | - | 0 attached | 25144 |
| 2 | armblock | 0 | 0 | 0 | {'never': 116} | 17960 |

