# DESIGN.md — v3.1 (2026-10-28, after review round 3 — final)

History: v1 → v2 after round 1 (6 reviewers), v2 → v3 after round 2 (5 reviewers), v3.1 after round 3 (final; docs/reviews/round3_final.md: B1 bridge CI by protein bootstrap with size standardisation, B2 genotype-only lockbox strata, B3 lockbox read order); responses in §9 and docs/reviews/round{1,2}_response.md. Fixed by PROMPT.md: (Part 1) a single locked evaluation of `zavg_BEF` + Part 1 baselines on the withheld hu.MAP test pairs, AUPRC primary; (Part 2) data-driven complexes from the locked co-complex scores and complex-level, cancer-context-specific dependencies; no STRING / hu.MAP / curated complexes in construction or validation; submission-ready manuscript.

**Paper thesis (pre-registered):** a label-free co-complex predictor, validated on withheld complexes, yields a curation-free complex map that agrees better with orthogonal physical evidence than the GLS-module standard; complex-level analysis of cancer contexts reveals driver- and lineage-linked complex dependencies and how loss of a subunit rewires dependency on the rest of its complex (paralog compensation, residual sensitisation, collapse) and on sub-modules.

## 0. Firewalls and process guardrails
- F1 Part 1 frozen at c7f95ee (external/part1_locked). Nothing learned here changes it.
- F2 Part 2 code must not read `datasets/huMAP/`, `results/test_eval/` or hu.MAP train labels (grep test). Part 2 parameters may not cite test-set metrics.
- F3 Construction/selection inputs: CRISPRGeneEffect on set A (+ locked predictor), Model, Confounders, HGNC location/symbols, Ensembl coordinates. Validation inputs (never used to select): BioPlex 3.0, CCLE proteomics (A models only for tiering), DEMETER2, lockbox L. Context inputs: DepMap Omics, Model.csv, Ensembl 116 paralogs. Post hoc only: HGNC `gene_group`, literature naming. Disclosure: zavg_BEF was selected on hu.MAP train labels (Complex Portal-derived), so agreement with physical/curated complexes partly re-measures that training.
- F4 Freeze points (commit + SHA-256 manifest; hashes reported to the user in chat; scripts verify):
  - FZ-P1 test protocol (marker file) → look runs immediately after.
  - FZ-LB lockbox model list (needs the context table's *model labels* for stratification, not outcomes).
  - FZ-C complex map on A + tiers + member/attachment rules.
  - FZ-X context table + family enumeration (m per family) + lockbox-testable list + H1 event list.
  - FZ-R registry entries (`docs/registry/*.yaml`, machine-readable) committed before first run; a lockbox target list committed before each lockbox read.
- F5 All associations via `run_association(analysis_id, ...)` (refuses unregistered IDs/unfrozen inputs; appends to `results/part2/trials.jsonl`: commit, input hashes, n tests, summary counts — no top hits). Exploration uses an A-only effect matrix. Lockbox outcome data (CRISPR, DEMETER2 of L models) are only read via `read_lockbox(analysis_id)` with a ledger; grep test that nothing else touches L outcomes. One lockbox read per closed registry entry. Order: all registered confirmations (controls, H3a, H3b, H1, H2) first; the complex-vs-gene benchmark and the bridge's lockbox-based secondary endpoint only after every confirmation is closed.
- F6 Every headline number re-derived by an agent that did not produce it, from raw files. References and biological claims checked by an agent that did not write them.

## 1. Part 1 — independent test evaluation
Fully specified in docs/TEST_PROTOCOL.md: single look after freeze; primary AP over all test pairs at natural prevalence (ITT); the locked Part 1 "clearly beats" verdict computed on S1 (mapped-only, Part 1's evaluation universe) governs the outcome reading, ITT quoted alongside and any disagreement reported; AUROC, P/R/F1 at train-fitted thresholds; isotonic/Platt calibration fitted on train; cluster-bootstrap CIs (B=2000); strata C1/C2/C3. All prior contact with test files disclosed (D0, D1). Outcome paragraphs pre-written both ways. Part 2 uses zavg_BEF regardless.

## 2. Part 2

### 2.1 Model sets
- U = CRISPRGeneEffect models minus engineered/drug-adapted (Model.csv EngineeredModel).
- Lockbox L (FZ-LB): 25% of U, seed 3000, PatientID-blocked, multilabel-stratified by OncotreeLineage × a genotype list fixed now (MSI-H; LoF of ARID1A, SMARCA4, STAG2, SMARCB1, PTEN, RB1, KMT2D, CREBBP, EP300; NRF2-active; KMT2A-r; 9p21/MTAP hom-del; VPS4B loss; KRAS, NRAS, BRAF, PIK3CA hotspots), so each has ≥5 L models where n_U permits. Stratification uses model labels only (no outcomes); H1 events enumerated later at FZ-X inherit whatever L counts result. A = U \ L.
- Seeds: clustering SEEDS=[0,1,2,3,4]; bootstraps 1000–1019; split-half 2000–2003; permutation 4000+; PYTHONHASHSEED=0.

### 2.2 Context table (FZ-X; complete-case per test, n reported, no imputation)
- MSI-H: MSIScore > 20. Log mutation burden: damaging-matrix row sum.
- LoF: damaging == 2 (biallelic) OR damaging ≥1 + copy loss (abs CN < ploidy − 0.75) OR homozygous deletion; X-linked genes in males: damaging ≥1. Primary excludes MSI-H models from non-MMR LoF contexts. Pre-registered functional check: LoF lines show lower expression (truncating/deletion) or lower protein (CCLE proteomics, where measured) of the gene than non-LoF lines; genes failing the check are reported and dropped from H1.
- Hotspot: hotspot matrix ≥1, separate from LoF.
- CN: r = OmicsCNGeneWGS (ploidy-relative linear, confirmed: housekeeping median 1.04 in WGD and non-WGD), clipped at 10; a = r × Ploidy. Hom-del: a < 0.5 (r < 0.25 if Ploidy NaN) AND expression log1p TPM < 1 where measured; partial loss: 0.5 ≤ a ≤ Ploidy − 0.75; amp: r ≥ 2; high: r ≥ 3. Adjacent genes with Jaccard ≥ 0.9 context vectors merged into segments named by locus (e.g. "9p21 loss").
- Fusions: unordered partner pairs from a pre-registered list (EWSR1::FLI1/ERG; SS18::SSX1/2/4; KMT2A::{MLLT3, AFF1, MLLT1, MLLT10, AFDN, ELL}; PAX3/PAX7::FOXO1; BCR::ABL1), matched on gene symbols parsed from "SYMBOL (ENSG…)"; self-fusions/read-throughs excluded; counts asserted against expected minimums (fails loudly if 0).
- Other: sex; loss of Y (male line whose summed chrY gene expression falls below a threshold fixed at FZ-X from the bimodal male distribution, female lines as reference); WGD; TP53 status; MYC/MYCN amp; NRF2-active (KEAP1/CUL3 LoF or NFE2L2 hotspot); lineage (OncotreeLineage) and melanoma (OncotreePrimaryDisease, any site).
- Sizes: confirmatory ≥12 context models in A; 8–11 exploratory; <8 descriptive.
- Covariate rule: no covariate that is a function of the context (e.g. MSIScore/burden are not covariates for MSI-H; ploidy not for WGD).

### 2.3 Complex construction (on A; label-free; FZ-C)
- Scores: zavg_BEF recomputed on A with the locked configuration (null moments recomputed per run with the locked procedure; on all U models it reproduces Part 1 scores to ≤1e-12 — test); stored as the S×S block over the filtered set S. Bridge comparator: GLS (Part 1 dev-selected variant) on A, same S; Pearson secondary.
- Gene set S: effect < −0.5 in ≥1% of A models, non-NaN in ≥95% of A models, not Part 1 low-coverage (≈6.2k genes). Attachment (for context genes outside S, e.g. LoF subunits): gene g is attached to complex C if ≥50% of C's members have score ≥ z0 to g and a size-matched null gives FDR < 0.01 over all gene×complex tests; attached genes can define H1 events but are not members.
- Graph: score ≥ z0 AND mutual top-k within S; weight = min(score, z0+3) − z0 + 1. GLS graph: same S and k, threshold set to match the zavg edge count, weights quantile-matched.
- Clustering: consensus co-clustering frequencies over 20 PatientID-blocked 80% bootstraps of A × 5 Leiden seeds on the union of edges; Leiden-CPM on frequencies over the γ grid {0.05·√2^i, up to 1.6}; persistence (best-match Jaccard ≥ 0.75, support ≥ 0.6, ≥3 adjacent γ); coarse level = persistent clusters at the coarsest persistent γ; sub-modules = top-down re-clustering inside each coarse cluster at γ=0.5 on the same frequencies; size 3–60; cohesion = all-pairs mean score ≥ z0 + 1 (applied to sub-modules; coarse clusters must pass at mean ≥ z0). Exact specification follows docs/reviews/round2_algorithms.md.
- Selection of (z0, k) ∈ {4,5,6} × {10,20}: split-half recovery (R=4 PatientID-blocked halves of A; clusters of size 5–30 recovered at Jaccard ≥ 0.5) weighted by coverage of selective genes; ties → smaller graph. Validated on toy data (no dyads, no arm-block reward, null → nothing).
- Locus flag: ≥50% of members within 10 Mb (Ensembl coordinates) or on one arm → "locus cluster", excluded from biological claims.
- Tiers (fixed at FZ-C): T1 within-cluster BioPlex edges (Entrez-mapped, 293T ∪ HCT116; ≥2 edges) above a node-restricted, degree-preserving null, BH over clusters q<0.05; T2 (not T1) proteomic co-abundance (A models only) above a size-matched null, BH q<0.05; T3-neg (tested, not supported) and T3-nc (insufficient BioPlex/proteomic coverage). "Complex" for T1/T2; "co-dependency module" otherwise.
- Bridge (primary methods endpoint, pre-registered): size-standardised within-cluster BioPlex edge precision (fraction of within-cluster pairs that are BioPlex edges, among pairs with both proteins BioPlex-observed; computed within cluster-size bins {3–4, 5–8, 9–16, 17–60} and averaged with the bin weights of the pooled zavg+GLS cluster-size distribution) for the frozen zavg vs GLS maps at matched gene coverage; CI for the paired difference by protein-node bootstrap (B=2000, seed 5000, clusters fixed, pairs re-weighted by node multiplicity); claim "better" iff the 95% CI excludes 0. Map-to-map variability over the 20 bootstrap maps reported descriptively. Secondary: proteomic co-abundance; stability by size bin; genome-wide pair-level AP against BioPlex for zavg vs GLS vs Pearson; H3 lockbox confirmation rate in matched size bins. Fail paragraph pre-written.
- Sensitivity only (COULD): Infomap/MCL agreement (ARI), BioPlex-integrated map.

### 2.4 Complex scores
- Primary: S_C(m) = mean of fixed members' CRISPRGeneEffect, winsorised ±3; context gene(s), their paralogs within C, and members within ±10 Mb of a CN-context gene removed for that test.
- Multiverse (reported per hit, never selected): eigengene, within-line rank, expression-masked, pan-essential-referenced (screen-quality control).
- Robustness: LOSO and leave-top-2-out (same sign, p<0.05); breadth.

### 2.5 Association testing
- Model: S_C ~ context + lineage FE (non-lineage contexts) + library multi-hot + ScreenType3DO + ModelType; mutation contexts + log burden + MSIScore; CN contexts + ploidy + WGD + TP53 status; covariate rule of §2.2 applies. Sensitivity: NNMD, doubling time, CasActivity, GrowthPattern.
- Statistic: HC3 t, Satterthwaite df; calibration by PatientID-blocked within-lineage Freedman–Lane permutation on random subsets and negative controls; hits (q<0.1) re-tested by exact Freedman–Lane permutation (GPD tail <1e-4), and BH re-run with permutation p-values.
- FDR: Benjamini–Bogomolov per direction: families = contexts; family selection by Simes; within-family BH at q·R/M. Sub-modules stage-wise within passing complexes.
- Discovery on A. Confirmation in L: same sign, one-sided p<0.05 (contexts with ≥5 L models). Rarer: LOSO + DEMETER2 (A models) same sign + leave-one-model-out → "not lockbox-testable". Ranking by A statistic × tier; L effects with CIs reported for all hits, confirmed or not.
- Paper's strongest findings (rule fixed): discovery q<0.05 AND lockbox-confirmed (or the rare-context route, labelled) AND LOSO-robust AND not a locus cluster; ranked by A statistic × tier. All directions reported incl. nulls, with trial counts.
- Complex vs gene benchmark: per (complex, registered context), L z of the complex score minus L z of its best member chosen in A; PatientID-blocked bootstrap over L; expectation stratified by breadth in A (complex wins when breadth is high). Secondary: A–L correlation vs matched random sets.
- Negative controls: pseudo-contexts matched for size and lineage composition; random gene sets matched for size, essentiality and coherence; dispersed never-expressed sets.

### 2.6 Directions (registry order; each a registry entry)
- H3a Lineage/subtype atlas (foundation; one figure; OXPHOS/mito confounders folded in as supplement). Family per lineage.
- H3b Driver/genome-state family (~25 contexts fixed at FZ-X): KRAS, NRAS, BRAF, PIK3CA hotspots; TP53 status; MSI-H; WGD; 9p21 loss; MYC amp; MYCN amp; LoY; sex; NRF2-active; ARID1A, SMARCA4, STAG2, PTEN, RB1, KMT2D, CREBBP LoF; registered fusions with n ≥ 12; tested against all frozen complexes.
- H1 Within-complex genetic interactions (headline; a panel, expected 10–20 events enumerated at FZ-X): event = (complex C, lost gene g ∈ C or attached), loss = LoF / hom-del / silencing. Estimand β = effect of loss(g) on S_{C∖g}, member-resolved profiles reported. Classes: paralog compensation (β < 0 driven by paralogs of g; reported separately), residual sensitisation (non-paralog members more dependent), collapse (members' dependency shrinks toward 0 relative to baseline sign), neutral (90% CI within ±0.1 Chronos, TOST), inconclusive. A class is called only if β lies outside the 2.5/97.5% quantiles of the same context's effects on random sets matched on essentiality, coherence and size. Complex×complex: COULD.
- H2 Sub-module divergence (headline): sub-module score residualised on its sibling complement (rest of the parent) ~ lineage/subtype/H3b contexts; NNMD and CasActivity primary covariates; controls CDK6/CCND3 vs CDK4/CCND1 (haematopoietic vs breast) and ncBAF vs cBAF (SMARCB1-null/SS18::SSX), mapped to sub-modules after FZ-C.
- Cut: H4 CN dosage (remaining idea), v1 D6/D7.
- Kill rule: a direction is dropped (reported) if <50% of its evaluable complex-level controls are recovered. Control→direction map: H3a: PC17, PC20; H3b: PC1, PC2, PC4, PC5, PC9, PC10, PC11, PC13, PC18, PC19a/b; H1: PC6, PC7, PC8, PC12, PC15; H2: PC21, PC22; gene-level only: PC3, PC14, PC16.

### 2.7 Pre-registered controls (one script, run once after FZ-C + FZ-X; HC3 one-sided test at α=0.05 in A for direction checks, lockbox confirmation reported; fixed denominator; "not evaluable at complex level" pre-declared when the gene is unclustered). Sign + = more dependent in context.
| ID | Context | Target | Expected |
|---|---|---|---|
| PC1 | MSI-H | WRN (G; C: WRN cluster) | + |
| PC2 | TP53 WT vs altered | MDM2, PPM1D, USP7 (C: MDM2 cluster) | + |
| PC3 | TP53 WT | TP53 / CDKN1A / USP28 (G) | − (growth gain) |
| PC4 | KRAS or NRAS hotspot | SHOC2 (C: SHOC2 cluster) | + |
| PC5 | BRAF hotspot vs RAS hotspot | SHOC2 | − |
| PC6 | ARID1A LoF (MSS) | ARID1B | + |
| PC7 | SMARCA4 LoF (MSS) | SMARCA2 | + |
| PC8 | STAG2 LoF | STAG1 | + |
| PC9 | 9p21 / MTAP hom-del | PRMT5, WDR77 (C: PRMT5 cluster) | + |
| PC10 | FOCAD deletion or MSI-H | PELO / HBS1L (C: PELO cluster) | + |
| PC11 | NRF2-active | NFE2L2 | + |
| PC12 | VPS4B partial/hom loss | VPS4A | + |
| PC13 | Female or LoY male | DDX3X | + |
| PC14 | RB1 LoF | CDK4/CDK6/CCND1 (G) | − |
| PC15 | SMARCB1 LoF or SS18::SSX | BRD9 (C: ncBAF cluster) | + |
| PC16 | KMT2A-r (canonical) | MEN1 / MLLT1 (G) | + |
| PC17 | Melanoma (any site) | SOX10 / MITF cluster | + (expression-explained) |
| PC18 | WGD | KIF18A | + |
| PC19a | PTEN LoF | PIK3CB | + |
| PC19b | PIK3CA hotspot | PIK3CA | + |
| PC20 | Myeloid/lymphoid lineage | MYB cluster | + |
| PC21 | Haematopoietic vs breast | CDK6/CCND3 vs CDK4/CCND1 sub-module contrast | + (heme: CDK6/CCND3) |
| PC22 | SMARCB1 LoF or SS18::SSX | ncBAF vs cBAF sub-module contrast | + (ncBAF) |

## 3. Verification
Toy simulations (planted complexes, sub-modules, arm blocks, lineage structure, planted context effects) for construction, H1 classifier, FDR control and HC3/permutation calibration; independent re-derivation of FZ-C and every headline; clean-clone `scripts/verify/reproduce.sh`.

## 4. Figures and stage gates
Fig 1 Part 1 test look; Fig 2 complex map, tiers, bridge; Fig 3 H3a/b; Fig 4 H1; Fig 5 H2. Gates: G1 test look → G2 FZ-C + bridge → G3 controls → G4 directions in registry order (H3a, H3b, H1, H2). Null versions of Figs 2 and 4 pre-written.

## 5. Compute
Lean locked recompute ~132 s / 2 GB per run (+ S×S block ~150 MB); 20 bootstraps + 8 half runs + GLS ≈ 1 h 45 min serial. Association screening analytic; permutations only for hits.

## 9. Responses
Round 1 (all MUSTs accepted; summary in v2 history, detail docs/reviews/round1_response.md). Round 2 accepted: bridge endpoint = BioPlex precision at matched coverage with paired bootstrap CI (pi M1); H3b driver family (pi M2); H1 re-scoped as panel with TOST neutral class, matched-set null, paralog split, member-resolved (pi M3, stats B-M3, biology 2); T1 BH, T3 split, H2 vs sibling complement (pi S, biology 4); controls mapped to directions, figure plan, gates (pi S6, biology 5); lockbox stratified by rare contexts (stats B-M1, biology 5); Benjamini–Bogomolov (stats); benchmark estimand (stats B-M4); winner's-curse ranking on A (stats S1); covariate-of-context rule (biology 1); PC5/PC11/PC17/PC19 fixes, LoF functional check, fusion assertions (biology); A-only effect matrix + lockbox ledger + machine-readable registry (ai_failure M-C1); subset null moments + reproduction test (ai_failure M-C2, algorithms 1); S×S storage, clustering specification, split-half statistic, GLS matching (algorithms 2–4); coordinates requested from user (algorithms 5). Cut per pi: H4, BioPlex-integrated map and Infomap as primary work (COULD only). Not adopted: user escrow of test files (procedural guarantee stated; hashes reported to the user).
