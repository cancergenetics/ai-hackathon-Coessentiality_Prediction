# DESIGN.md — v1 (2026-10-28, pre-review)

Fixed by PROMPT.md: (Part 1) single, locked evaluation of the Part 1 final method `zavg_BEF` and the Part 1 baselines (Pearson, Spearman, GLS) on the withheld hu.MAP test pairs, AUPRC primary; (Part 2) data-driven complexes from the locked co-complex scores (optionally + BioPlex) and complex-level, cancer-context-specific dependencies; no STRING / hu.MAP predictions / curated complexes in construction or validation; submission-ready manuscript. Ideas from LITERATURE.md §7 are priors, not commitments.

## 0. Firewalls (apply to every agent)
- F1 Part 1 is frozen at commit c7f95ee (external/part1_locked). Nothing learned in this project (test results or Part 2) may change it. Code bugs found later are reported, not fixed in the scored run.
- F2 hu.MAP test pairs are used only by the Part 1 test look (§1) and its verifier. Part 2 never reads them.
- F3 Construction / validation inputs for Part 2: DepMap (CRISPRGeneEffect, ScreenGeneEffect, confounders, Model, ModelCondition, Omics files requested), DEMETER2 RNAi, BioPlex 3.0, Ensembl paralogs, CCLE proteomics, HGNC *non-curated* fields (symbol, location, locus type). Post hoc only (never select parameters or claim validation): HGNC `gene_group`, literature knowledge used for naming, and any curated complex resource.
- F4 Positive controls (known biology) are pre-registered (§2.7) and reported, never used to tune construction or test parameters.
- F5 Every headline number is re-derived by an agent that did not produce it.

## 1. Part 1 — independent test evaluation
Pre-registered in docs/TEST_PROTOCOL.md (frozen by commit before any test score is computed; commit hash recorded in the results). Summary: locked predictor `scripts/predict_pairs.py` from external/part1_locked; methods zavg_BEF (primary), Part 1's dev-selected Pearson / Spearman / GLS variants, B (Part 1's pre-declared secondary). Primary metric AP at the test file's natural prevalence over all test pairs (intention-to-score: unmapped pairs get Part 1's minimum−1 rule); secondary AUROC, P/R/F1 at decision rules fixed from development data only, calibration (isotonic/Platt fitted on train pairs), macro per-cluster AP; paired cluster-bootstrap CIs (B=2000) for each AP and ΔAP vs baselines. Transition documentation and leakage audit: docs/part1_audit.md (PASS). Single look; independent re-derivation.

## 2. Part 2 — context-specific complex dependencies

### 2.1 Data and cell-line set
- Models: CRISPRGeneEffect (1,361 models). Exclude engineered / drug-adapted models (16) from context analyses; keep organoids but carry ModelType as covariate. Model metadata from Model.csv; screen covariates from CRISPRConfounders (library indicators, NNMD, doubling time, CasActivity, 3D screen).
- Contexts (C-types):
  - L: lineage (OncotreeLineage, ≥15 models) and subtype (OncotreePrimaryDisease / OncotreeCode with ≥8 models) one-vs-rest, and within-lineage subtype contrasts.
  - M: damaging (LoF) mutations and hotspot mutations of genes mutated in ≥8 models (≥5 within a lineage for within-lineage tests).
  - CN: deep deletion / partial loss / amplification of genes (gene-level relative CN; thresholds fixed a priori: deep del log2(rel CN) < −1.2 ⇒ ≈homozygous; partial loss −1.2…−0.4; amp > 0.8 — confirmed against the file's scale before any association is computed).
  - FU: recurrent fusions (≥5 models).
  - SIG: MSI status, ploidy / WGD from OmicsSignatures (if present).
  - EXP: low expression (TPM log1p < 1) of a complex subunit or paralog — used for paralog switching (D4) and as a guard (non-expressed subunits flagged).
- Seeds: global list SEEDS = [0, 1, 2, 3, 4]; bootstrap seeds 1000–1199; split-half seeds 2000–2009. PYTHONHASHSEED=0 for all graph code.

### 2.2 Complex construction (label-free, frozen before any context test)
- Input: genome-wide zavg_BEF (Part 1 outputs/all_pairs; recomputable with the locked predictor).
- Graph G(z0, k): edge (i,j) iff zavg ≥ z0 AND j in i's top-k AND i in j's top-k (mutual kNN; removes outlier-driven hubs observed in the profile: non-essential genes with degree ≥100). Weight w = zavg − z0 + 1.
- Gene filter (label-free): exclude low-coverage genes (Part 1 flag) and genes whose effect SD across lines is below the 5th percentile of all genes (no information); keep pan-essential genes (they form the core machines).
- Clustering: Leiden CPM over γ grid; consensus over 20 seeds and over 20 cell-line bootstrap recomputations of zavg_BEF (resample 80% of models without replacement, locked predictor configuration unchanged; ~150 s each) → co-clustering frequency matrix → final clusters = Leiden on consensus graph (edges with co-cluster freq ≥ 0.5). Multiscale: two γ levels define complexes (coarse) and sub-modules (fine) with nesting enforced.
- Size: 3–60 members; larger clusters re-clustered at finer γ.
- Overlap (secondary): soft membership — gene g added to cluster C if mean zavg(g, C) ≥ the 99.9th percentile of a size-matched random-set null; used only for moonlighting annotation.
- Parameter selection (z0 ∈ {3,4,5,6}, k ∈ {10,20,50}, γ grid) by label-free criteria only: (i) stability = mean best-match Jaccard of clusters across cell-line bootstraps; (ii) split-half replication (clusters from lineage-stratified half A recovered in half B); (iii) coverage of genes with selective dependency. Rule: maximise stability subject to median cluster size 3–15 and ≥1,000 genes covered; ties → simpler (higher z0).
- BioPlex integration (variant, compared not assumed): G' = G ∪ {BioPlex edges with zavg ≥ z1} (BioPlex edges supported by weaker co-essentiality). Cross-validated by cell line of the AP-MS: integrate 293T edges, evaluate enrichment of held-out HCT116-only edges within clusters (and vice versa). The integrated variant is adopted only if it raises stability and held-out BioPlex enrichment without reducing gene coverage of selective genes; otherwise BioPlex is used for annotation/validation only.
- Construction quality (reported, not optimised): within-cluster BioPlex edge enrichment vs degree- and size-matched random sets; within-cluster CCLE proteomic co-abundance vs random; post hoc agreement with HGNC gene groups (labelled post hoc).
- Freeze: complexes C1 (+ sub-modules) committed with checksums before §2.3 starts.

### 2.3 Complex-level dependency score
- Per model m and complex C: S_C(m) = mean over members of CRISPRGeneEffect (Chronos scale; primary). Sensitivity: eigengene of standardised member profiles (sign-aligned), and singscore-like within-line rank score. Members flagged non-expressed in m are excluded from that model's mean (once expression arrives; sensitivity analysis without this rule).
- Leave-one-subunit-out (LOSO): every reported complex association must hold (same sign, p<0.05) after dropping any single member; else reported as "subunit-driven".

### 2.4 Association testing (shared machinery)
- Model: S_C ~ context + lineage fixed effects (for non-lineage contexts) + library indicators + ScreenType3DO; sensitivity models add NNMD, doubling time, CasActivity, ModelType, GrowthPattern. For lineage contexts: one-vs-rest with library covariates (lineage is the context).
- Inference: OLS t-statistic; p-value by permutation of context labels within lineage strata (≥10,000 permutations or until p resolved, fixed seeds); BH FDR within each context family; effect size = Δ mean Chronos (context vs rest) and standardised d.
- Discovery / validation: models split into two halves stratified by lineage (and context) — 10 split seeds; a finding is "replicated" if FDR<0.1 in discovery and same-sign p<0.05 in validation in the canonical split (seed 2000) and the median over 10 splits passes. Primary inferential claim = full-data FDR<0.05 AND split-half replication AND LOSO.
- Orthogonal replication: DEMETER2 RNAi (overlapping models, same model), and ScreenGeneEffect screens from a different library for the same model (technical).
- Complex vs gene-level: the same pipeline on single genes provides (i) the benchmark for whether complex-level aggregation gives more power (number of replicated associations, replication rate, held-out effect size) and (ii) the comparison "complex signal not explained by its single best subunit" (partial test: complex score residualised on best subunit).
- Negative controls: size- and variance-matched random gene sets (null distribution of association counts); label-permuted contexts; olfactory-receptor / non-expressed gene sets (should show nothing).

### 2.5 Research directions (to be pursued in order of priority; each with pre-registered positive controls)
- D1 Lineage/subtype complex-dependency atlas (rank 1, data present).
- D2 Sub-module divergence: (sub-module score − parent score) ~ context; which part of a complex carries the context dependency; sub-modules built on half A and tested on half B.
- D3 Genotype-conditioned complex dependencies: complex-level disruption ("any core subunit with LoF/deep deletion") × complex dependency, including complex×complex (restricted a priori to complex pairs connected in zavg_BEF top edges or BioPlex), lineage FE, MSI/mutation burden covariates.
- D4 Paralog switching within complexes: in-complex paralog pairs (Ensembl); dep(P1) ~ status(P2) (expression / LoF / CN), complex-level score stable while subunit identity switches.
- D5 Copy-number dosage (CYCLOPS in CRISPR): dep(g) ~ CN(g) for complex subunits vs segment-matched non-complex genes; partner-subunit dependency ~ CN of lost subunit; RNAi replication.
- D6 (exploratory) Complex coherence rewiring across lineages with variance-matched null. D7/D8 only if time.

### 2.6 "Method" optimisation and runtime (auto-research)
- The Part 2 method = graph → clusters → complex score → association. Variants (graph z0/k, BioPlex integration, consensus vs single run, score type) are compared only on label-free criteria (§2.2) and on split-half replication rate of context associations (a held-out, data-internal criterion). Runtime and peak memory are logged for every stage.

### 2.7 Pre-registered positive controls (reported as recovered / not recovered; never tuned on)
- Lineage: MYB/ CBFB–RUNX1 / Mediator kinase or MLL–menin in myeloid / lymphoid; SOX10 / MITF modules in melanoma (Skin); HNF1B/ HNF4A-type TF modules in GI lineages; OXPHOS in AML.
- Mutation / fusion / CN: ARID1A-LoF → ARID1B; SMARCA4-LoF → SMARCA2; KEAP1-LoF → NFE2L2-axis loss of dependency; MSI → WRN; MTAP deletion → PRMT5/WDR77 (and MAT2A); STAG2-LoF → STAG1; SS18::SSX → BRD9/ncBAF; KMT2A-r → MEN1/KMT2A/DOT1L module; EWSR1::FLI1 → lineage modules; PSMC2/SF3B1 partial loss → increased dependency (CYCLOPS).
- Negative: olfactory-receptor modules (no context effect); random gene sets.

### 2.8 Verification
- Toy simulations with planted complexes (known membership, noise, hubs) → clustering recovers them (ARI); planted context effects of known size in complex scores → pipeline recovers effect and controls FDR under null (calibration of permutation p-values, QQ).
- Independent re-derivation of each headline finding from raw files by a different agent.

## 3. Deliverables
- results/test_eval (Part 1 test), results/part2/*, manuscript/ (main text, figures, tables, supplementary data/text, cover letter, data/code availability), REPORT.md.

## 4. Compute budget
Laptop (8 cores, 24 GB). zavg_BEF genome-wide recompute ~150 s, 6.6 GB peak → 20 bootstraps ≈ 50 min. Leiden ~1 s per run. Permutation tests: vectorised, ≤10 min per context family.
