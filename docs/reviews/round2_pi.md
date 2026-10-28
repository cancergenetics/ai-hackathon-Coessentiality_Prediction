# Design review round 2 — external PI (role: pi)

No hu.MAP test file opened; no complex–context association computed.

**Feasibility check** (label-free, 30 s; script `h1_feas.py` in my scratchpad):
- **What I counted:** biallelic-LoF/hom-del events under the v2 §2.2 definitions, in 1,145 non-MSI models. The X-linked male rule was omitted.
- **Result:**

| Genes with total LoF in | Genes | Of which pass the graph filter |
|---|---|---|
| ≥16 models | 95 | 26 |
| ≥20 models (≥5 in L) | 80 | 17 |

- TP53, CDKN2A, MTAP and RB1 are never dependency-bearing, so they cannot be cluster members.

## Round-1 MUSTs
- **Resolved:** all of them (pi, statistics, ai_failure, biology, data, algorithms). D1 is disclosed in protocol §0/§5.4.
- **Gap:** §9 cites `docs/reviews/round1_response.md`, which does not exist.

## Verdict
v2 is coherent, feasible and much tighter. Three problems remain before execution:
1. the bridge cannot yet fail cleanly;
2. H1 is a small panel, not a systematic map;
3. PROMPT's mutation and CN contexts have no registered direction.

## MUST

**M1. Make the bridge decidable.**
- **Problem:**
  - §2.3 gives two incompatible GLS graph rules: "matched by edge count" vs "same rule applied".
  - Five endpoints are listed, with no primary and no decision rule.
  - Endpoint (v), the A-vs-L effect correlation, is dominated by lineage. It is also confounded by the maps' different cluster sizes.
- **Change:**
  - **GLS graph:** GLS uses zavg's k, with its threshold matched by edge count. GLS's own split-half selection is secondary.
  - **Primary endpoint:** within-cluster BioPlex precision.
    - Count only pairs where both genes are observed in BioPlex.
    - Compare at matched gene coverage.
    - Take the CI from the paired 20 bootstrap maps.
    - "Better" means the CI of the difference excludes 0.
  - **Secondary endpoints:**
    - proteomic co-abundance;
    - the H3 lockbox confirmation rate in matched size bins.
  - **Fail paragraph:** pre-write it.
  - **Disclose:** physical endpoints partly re-measure the Complex Portal-trained gain.
- **Why:** this is the claim that unifies the paper.

**M2. Register a driver/genome-state family.**
- **Problem:** H3 covers lineage only, and H1 covers only in-graph subunits. KRAS/BRAF/PIK3CA hotspots, TP53, MSI, WGD, 9p21, MYC amplification, LoY and fusions therefore appear only as controls, yet PROMPT names them.
- **Change:** split H3 into two families, each with its own BH:
  - **H3a:** lineage/subtype;
  - **H3b:** about 25 driver/genome-state contexts, fixed at FZ-X, tested against all frozen complexes.
- **Why:** this is the core requested analysis. It is well powered and lockbox-testable.

**M3. Re-scope H1 and make its classification well posed.**
- **Problem:**
  - Only about 26 events are eligible, and 17 are lockbox-testable; roughly 10–20 will land in clusters.
  - The best-known events are already controls (PC6–8, PC15).
  - "Neutral" has no equivalence margin.
  - The selection rule's lockbox requirement excludes about 35% of H1.
- **Change:**
  - **Unit:** the pair (lost subunit g, complex C).
  - **Primary statistic:** the score of C without g.
  - **Secondary:** a member-resolved Δ profile, tested stage-wise. The paralog arm becomes an annotation ("is the most-shifted member g's paralog?"), not a separate family.
  - **Classes:**
    - buffered: Δ < 0, confirmed;
    - collapsing: Δ > 0, confirmed;
    - equivalent: TOST 90% CI within ±0.1 Chronos;
    - otherwise inconclusive.
  - **Confirmation:** events with fewer than 5 L models use the LOSO + DEMETER2 + leave-one-model-out route and are listed as a separate class.
  - **Scope:** move complex×complex to COULD, and frame H1 as a mechanistic panel.

## SHOULD

**S1. Tiers.**
- Problem: per-cluster p < 0.05 over about 1,000 clusters gives about 50 false T1 clusters. T3 also mixes negative clusters with uncovered ones.
- Change: use BH q < 0.1, and split T3 into T3-neg (≥50% of members observed) and T3-nc (not covered).
- Why: T3-neg answers "not complexes".

**S2. H2 score.**
- Problem: sub-module − parent is mechanically coupled, because the parent contains the sub-module.
- Change: use sub-module − (parent \ sub-module).

**S3. Kill rule.**
- Problem: H2 has no complex-level controls and H3 has one.
- Change: map each control to a direction. With fewer than 3 evaluable controls, the rule is "n/a" and the direction is labelled exploratory.

**S4. Selection.**
- Change:
  - Replace "effect × tier" with a lexicographic rank: tier first, then lockbox effect.
  - For n_L < 5, report effects from A with a winner's-curse caveat.

**S5. Cut.**
- Cut:
  - H4;
  - the BioPlex-integrated map, which blurs the role of BioPlex as validation;
  - Infomap (keep MCL).
- Keep the multiverse and LOSO.

**S6. Figure plan and stage gates.**
- Figures:
  - Fig 1: Part 1;
  - Fig 2: map, tiers and bridge;
  - Fig 3: H3a/b;
  - Fig 4: H1;
  - Fig 5: H2.
- Gates:
  - G1: the look;
  - G2: FZ-C, then the bridge;
  - G3: controls;
  - G4: the directions, in registry order.
- Pre-write the null versions of Figs 2 and 4.

## COULD
- **C1.** Stratify the lockbox on lineage × the ~20 pre-listed subunit-LoF genes (label-free).
- **C2.** Report a genome-wide pair-level AP against BioPlex for zavg vs GLS vs Pearson; it supports Part 1 orthogonally.

## No change needed
F1–F6, the lockbox/registry, the context definitions, HC3 + Freedman–Lane, the families, the controls and the compute plan.
