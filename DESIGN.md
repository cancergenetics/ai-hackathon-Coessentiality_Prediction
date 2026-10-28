# DESIGN.md — v2 (2026-10-28, after review round 1)

Changes from v1 are driven by docs/reviews/round1_*.md (pi, statistics, ai_failure, biology, data_priorwork, algorithms); the response table is §9. Fixed by PROMPT.md: (Part 1) a single locked evaluation of `zavg_BEF` + Part 1 baselines on the withheld hu.MAP test pairs, AUPRC primary; (Part 2) data-driven complexes from the locked co-complex scores (optionally + BioPlex) and complex-level, cancer-context-specific dependencies; no STRING / hu.MAP / curated complexes in construction or validation; submission-ready manuscript.

**Paper thesis (pre-registered, PI M2):** a label-free co-complex predictor, validated on withheld complexes, yields a curation-free complex map that is better than the GLS-module standard on orthogonal physical evidence; analysing dependencies at the complex level shows how cancers rewire dependency *within* complexes (paralog buffering vs complex collapse; sub-module specialisation) beyond published lineage atlases.

## 0. Firewalls and process guardrails
- F1 Part 1 frozen at c7f95ee (external/part1_locked). Nothing learned here changes it.
- F2 hu.MAP test pairs/results: used only by the Part 1 look and its verifier. Part 2 code must not read `datasets/huMAP/` or `results/test_eval/` (grep test enforces). hu.MAP *train* labels are also excluded from Part 2 construction and parameter choice. Part 2 parameters may not cite test-set metrics.
- F3 Construction/selection inputs: DepMap CRISPRGeneEffect (+ locked predictor), Model/Confounders. Validation inputs (never used to select): BioPlex 3.0, CCLE proteomics, DEMETER2, lockbox models. Context inputs: DepMap Omics (mutation, CN, expression, signatures, fusions), Model.csv, Ensembl 116 paralogs. Post hoc only: HGNC `gene_group`, literature naming. Disclose that zavg_BEF was selected on hu.MAP train labels (Complex Portal-derived), so agreement with curated complexes is partly built in.
- F4 Freeze points (each = commit with SHA-256 manifest; scripts check them):
  - FZ-P1: TEST_PROTOCOL frozen (marker file; see protocol).
  - FZ-LB: lockbox model list (before any Part 2 computation that uses models).
  - FZ-C: complex map (construction set only) + tiers + member-inclusion rule.
  - FZ-X: context table (all context definitions, thresholds after the CN-scale check), family enumeration with m per family.
  - FZ-R: each registry entry (docs/ANALYSIS_REGISTRY.md) committed before its first run.
- F5 All associations go through one wrapper `run_association(analysis_id, ...)` that refuses unregistered IDs or unfrozen inputs and appends to `results/part2/trials.jsonl` (commit, input hashes, n tests, summary counts; no top hits). Lockbox reads go through `read_lockbox(analysis_id)` with its own ledger; one read per registered hypothesis.
- F6 Every headline number is re-derived by an agent that did not produce it, from raw files, importing only the frozen inputs. References and biological claims in the manuscript are checked by an agent that did not write them.

## 1. Part 1 — independent test evaluation
Fully specified in docs/TEST_PROTOCOL.md (frozen before any test score; freeze marker + look-intent commit + git-tracked look ledger; hashes reported to the user in chat). Primary: AP at the natural prevalence of all test pairs (intention-to-score); zavg_BEF vs Part 1's dev-selected Pearson/Spearman/GLS on identical pairs; Part 1's locked "clearly beats" rule reported exactly as Part 1 defined it, plus CI-excludes-0; AUROC, P/R/F1 at train-fitted thresholds, isotonic/Platt calibration fitted on train, cluster-bootstrap CIs (B=2000); strata C1/C2/C3 (proteins seen in train). All prior contact with test files disclosed (D0, D1). Pre-written outcome paragraphs for pass and fail. Part 2 uses zavg_BEF regardless of the outcome; the bridge experiment (§2.3) carries "does it matter downstream".

## 2. Part 2 — complexes and context-specific dependencies

### 2.1 Model sets
- Universe U: CRISPRGeneEffect models (1,361) minus engineered/drug-adapted (16).
- Lockbox L (FZ-LB): 25% of U, sampled stratified by OncotreeLineage, blocked by PatientID (all models of a patient on one side), seed 3000. Lockbox models are excluded from construction, exploration, controls and direction triage.
- Construction/discovery set A = U \ L (~1,010 models).
- Seeds: SEEDS=[0,1,2,3,4] (clustering), bootstrap seeds 1000–1019, Leiden seeds 0–19, permutation seed 4000+, PYTHONHASHSEED=0.

### 2.2 Context table (FZ-X; complete-case per test; n reported per test; nothing imputed)
- Hypermutation: MSI-H = MSIScore > 20 (OmicsGlobalSignatures); log mutation burden from the damaging matrix row sums.
- LoF (tumour suppressors, complex subunits): damaging == 2 (biallelic) OR damaging ≥ 1 with copy loss (abs CN < ploidy − 0.75) OR homozygous deletion; X-linked genes in males: damaging ≥ 1. Primary analysis excludes MSI-H models from non-MMR LoF contexts; sensitivity includes them with burden + MSIScore covariates. Sensitivity "any damaging".
- Hotspot (oncogenes): hotspot matrix ≥ 1 (165 genes in the file), separate from LoF.
- CN: OmicsCNGeneWGS is a ploidy-relative linear ratio r (scale check logged before freezing thresholds; clip r at 10); absolute a = r × Ploidy where Ploidy exists. Homozygous deletion: a < 0.5 (or r < 0.25 when Ploidy is NaN) AND expression log1p TPM < 1 when expression exists; partial loss: 0.5 ≤ a ≤ Ploidy − 0.75; amplification: r ≥ 2; high-level r ≥ 3. Neighbouring genes with near-identical context vectors (Jaccard ≥ 0.9) are merged into segments ("9p21 loss"); attribution to a single gene only via conditional models or expression-null lines. "Has-CN" indicator as sensitivity covariate.
- Expression: low expression = log1p TPM < 1.
- Fusions: pre-registered unordered partner pairs only: EWSR1::FLI1/ERG; SS18::SSX1/2/4; KMT2A::{MLLT3, AFF1, MLLT1, MLLT10, AFDN, ELL}; PAX3/PAX7::FOXO1; BCR::ABL1. n<8 → descriptive.
- Other: sex (Model.csv) and loss of Y (male with summed chrY gene expression below a threshold set from female lines' distribution); WGD (OmicsGlobalSignatures); TP53 status (LoF or hotspot); MYC/MYCN amplification; lineage (OncotreeLineage ≥15 models in A) and subtype (OncotreePrimaryDisease ≥15 in A).
- Minimum sizes: confirmatory tests need ≥12 context models in A; 8–11 exploratory (except pre-registered controls); <8 descriptive.

### 2.3 Complex construction (on A only; label-free; FZ-C)
- Scores: recompute zavg_BEF on A with the locked predictor configuration (no Part 1 setting changed). Bridge comparator: GLS (Part 1's dev-selected variant) genome-wide on A; Pearson secondary.
- Gene filter (fixed now, from CRISPRGeneEffect on A only): graph genes = dependency-bearing genes, effect < −0.5 in ≥1% of A models (≈6.5k genes), non-NaN in ≥95% of A models, not Part 1 low-coverage. Neighbours taken within this set (exact for k ≤ 20 from per-gene top-50 lists). Never-essential genes (65% of genes; 70–81% of their strong edges are same-arm proximity artefacts, algorithms review) are excluded from construction.
- Compute: bootstraps via the lean locked pipeline (top-50 lists per gene, 2 GB, ~132 s each, serial; B=10 for (z0,k) choice, B=20 for the frozen map); subsamples = 80% of A drawn by PatientID.
- Graph: edge iff score ≥ z0 AND mutual top-k; weight = min(score, z0+3) − z0 + 1 (capped so paralog edges do not dominate CPM). For GLS/Pearson the threshold is matched by edge count to the zavg graph.
- Clustering (primary): consensus co-clustering frequencies on the union of bootstrap edges (denominator = all bootstrap × seed runs, 20 × 5 seeds); Leiden-CPM on frequencies over a γ grid; HiDeF-style persistence: keep clusters recurring (best-match Jaccard ≥ 0.75, support ≥ 0.6) over ≥3 adjacent γ values; nesting by top-down re-clustering inside each coarse cluster (complexes → sub-modules). Hard rules: size 3–60 (larger re-clustered), cohesion mean internal score ≥ z0+1.
- Selection (label-free only): (z0, k) ∈ {4,5,6} × {10,20} chosen by split-half recovery of selective-gene coverage (clusters built on PatientID-blocked halves of A recovered in the other half), not raw stability (raw stability is degenerate: it favours dyads). Stability reported by size bin with a degree-preserving-rewire chance baseline. Rule validated on toy data (must recover planted structure, not pick dyads, not reward arm blocks). Same rule applied to GLS.
- Sensitivity algorithms: Infomap and sparse MCL (inflation sweep) — ARI/overlapping NMI to the consensus reported, not "better". Soft overlap: g joins C if ≥50% of members have score ≥ z0 to g, size-matched null, FDR 0.01 over all gene×cluster tests; annotation only.
- Locus artefact flag: clusters with ≥50% members on one chromosome arm or within 10 Mb (HGNC location) are flagged "locus clusters" and excluded from biological claims.
- Evidence tiers (fixed at FZ-C, before any context test): T1 within-cluster BioPlex (293T ∪ HCT116) edge density above a degree-preserving null (p<0.05); T2 (not T1) CCLE proteomic co-abundance (mean within-cluster correlation) above a size-matched null; T3 co-essentiality only. Text uses "complex" for T1/T2 and "co-dependency module" for T3.
- Bridge experiment (primary methodological claim of Part 2): zavg vs GLS maps compared on (i) BioPlex within-cluster enrichment, (ii) proteomic co-abundance, (iii) bootstrap stability, (iv) fraction of genes in T1/T2 clusters, (v) downstream: correlation of context effect estimates between A and L for the same registered context families (threshold-free), and DEMETER2 concordance. Reported, never used to select.
- BioPlex-integrated map (G ∪ BioPlex edges with score ≥ z1, z1 = 2 fixed now) only as a sensitivity analysis; never used for tiering of itself. BioPlex enrichment nulls match BioPlex degree, essentiality class, expression and cluster size; held-out evaluation is by bait-protein K-fold, not by cell line.

### 2.4 Complex-level dependency score
- Primary (frozen): S_C(m) = mean over fixed members of CRISPRGeneEffect, winsorised at ±3; members excluded per the inclusion rule; context gene(s) and members within ±10 Mb of a CN-context gene are removed from the score for that test. Coarse level primary; sub-modules for D2.
- Sensitivity (multiverse, reported per hit, never selected): eigengene, within-line rank score, expression-masked members, score relative to the line's median over a disjoint pan-essential reference set (screen-quality control for core machines).
- Robustness per hit: leave-one-subunit-out (LOSO) and leave-top-2-out same sign with p<0.05; breadth (fraction of members with the same sign).

### 2.5 Association testing
- Model: S_C ~ context + lineage fixed effects (for non-lineage contexts) + library multi-hot + ScreenType3DO + ModelType; mutation contexts + log burden + MSIScore; CN contexts + ploidy + WGD (+ TP53 status). Sensitivity covariates: NNMD, doubling time, CasActivity, GrowthPattern.
- Statistic: HC3 robust t; screening p from t with Satterthwaite df; calibration checked by within-lineage, PatientID-blocked Freedman–Lane permutation on a random subset (QQ) and on negative controls. Hits (FDR<0.1) confirmed by exact Freedman–Lane permutation (GPD tail for p<1e-4).
- Families (m enumerated at FZ-X): per direction; BH within family; hierarchical BH over families (family-level Simes then within-family) for the cross-direction claim. Sub-modules tested stage-wise only inside complexes that pass.
- Discovery in A; confirmation in L: a finding is confirmed if lockbox effect has the same sign with one-sided p<0.05 (contexts with ≥5 lockbox models); rarer contexts: LOSO + DEMETER2 same sign + leave-one-model-out, labelled "not lockbox-testable". Effect sizes reported from L (cross-fitted, avoids winner's curse) with CIs.
- Selection rule for the paper's strongest findings (fixed now): discovery FDR<0.05 AND lockbox-confirmed AND LOSO-robust AND not a locus cluster; ranked by lockbox effect size × tier (T1>T2>T3). Every registered direction is reported, including nulls, with trial counts.
- Complex vs gene benchmark (pre-registered expectation: complex level wins for coherent shifts, loses for single-subunit effects): same gene universe, contexts, covariates; primary metric = correlation of A vs L effect estimates and lockbox confirmation rate at matched numbers of tests; nulls = random gene sets matched on size and mean within-set zavg.
- Negative controls: pseudo-contexts matched on size and lineage composition; coherence-matched random gene sets; dispersed never-expressed gene sets.

### 2.6 Directions (priority order; each a registry entry before running)
- H1 (headline) Within-complex genetic interactions (merges v1 D3+D4): when subunit g of complex C is lost (biallelic LoF / homozygous deletion / silencing), do the remaining members become more essential (buffering: paralog takeover, e.g. ARID1A→ARID1B) or less essential (collapse: complex non-functional)? Classify complex×subunit-loss events as buffered / collapsing / neutral. Paralog arm: in-complex Ensembl paralog pairs (both members in C or one in C and paralog anywhere), dep(P1) ~ loss/low-expression of P2. Secondary: complex×complex (restricted a priori to complex pairs connected by top zavg or BioPlex edges).
- H2 (headline) Sub-module divergence (cross-fitted): (sub-module score − parent score) ~ lineage/subtype/genotype; sub-modules from A, tested in L.
- H3 (foundation, one figure) Lineage/subtype complex-dependency atlas, with mitochondrial/OXPHOS confounder analysis (doubling time, medium, growth pattern) folded in as supplement.
- H4 (COULD) CN dosage: trans partner effects only (partner dependency ~ partial loss of a subunit), DEMETER2 primary, segment-matched null.
- Cut: v1 D6 coherence rewiring, D7 BioPlex cell-line rewiring.
- Kill rule: a direction is dropped (reported as such) if <50% of its evaluable pre-registered complex-level controls are recovered after the freezes.

### 2.7 Pre-registered controls (evaluated once, by a separate script, after FZ-C and FZ-X; fixed denominator; failures reported)
Class G (gene-level pipeline checks) and class C (complex-level: "the frozen cluster containing gene X"; "not evaluable at complex level" is a pre-declared outcome if X is unclustered). Sign: + = more dependent (more negative Chronos) in context.
| ID | Context | Target | Expected |
|---|---|---|---|
| PC1 | MSI-H | WRN (G) | + |
| PC2 | TP53 WT vs TP53 LoF/hotspot | MDM2, PPM1D, USP7 (G; C: cluster of MDM2) | + |
| PC3 | TP53 WT | TP53 / CDKN1A / USP28 knockout | − (growth gain) |
| PC4 | KRAS or NRAS hotspot | SHOC2 (C: SHOC2–MRAS–PPP1CB cluster) | + |
| PC5 | BRAF V600 | SHOC2 | no effect / − |
| PC6 | ARID1A LoF (MSS) | ARID1B | + |
| PC7 | SMARCA4 LoF (MSS) | SMARCA2 | + |
| PC8 | STAG2 LoF | STAG1 | + |
| PC9 | MTAP homozygous deletion (9p21 segment) | PRMT5, WDR77 (C: PRMT5 cluster) | + |
| PC10 | FOCAD deletion or MSI-H | PELO / HBS1L (C: cluster of PELO) | + |
| PC11 | KEAP1 LoF | NFE2L2 | + (gain of dependency) |
| PC12 | VPS4B partial/hom loss | VPS4A | + |
| PC13 | Female or loss of Y (male) | DDX3X | + |
| PC14 | RB1 LoF | CDK4/CDK6/CCND1 | − (loss of dependency) |
| PC15 | SMARCB1 LoF or SS18::SSX | BRD9 (C: ncBAF cluster) | + |
| PC16 | KMT2A-r (canonical) | MEN1 / MLLT1 (C: cluster of MEN1) | + |
| PC17 | Skin lineage (melanoma) | SOX10 / MITF cluster | + (expression-explained) |
| PC18 | WGD | KIF18A | + |
| PC19 | PTEN LoF vs PIK3CA hotspot | PIK3CB (PTEN) vs PIK3CA (PIK3CA) | + |

## 3. Verification
- Toy simulations with planted complexes (hubs, noise) → clustering ARI; planted context effects → effect recovery and FDR control; HC3/permutation type-I under unequal variance and small n1 (reuse statistics review simulation).
- Independent re-derivation of every headline result and of FZ-C from raw files; clean-clone `scripts/verify/reproduce.sh` regenerates committed outputs and compares hashes.

## 4. Deliverables
results/test_eval, results/part2/*, docs/ANALYSIS_REGISTRY.md, manuscript/ (main text, figures, tables, supplementary data/text, cover letter, data/code availability), REPORT.md.

## 5. Compute
Laptop (8 cores, 24 GB). Lean locked zavg recompute (top-50 lists) ~132 s / 2 GB → 30 runs ≈ 65 min serial (parallel runs slow ARPACK; not relied on). Leiden ~1 s/run. Association screening analytic; permutation only for hits.

## 9. Round-1 response (summary; detail in docs/reviews/round1_response.md)
Accepted: lockbox + construction on A (ai_failure M7, biology M6, data M3); bridge vs GLS (pi M2); scope cut + kill rule (pi M3); tiers T1–T3, BioPlex as validation (pi M4); hypermutation/biallelic LoF (pi M5, stats M2, biology M2); HC3 + Freedman–Lane + family enumeration (stats M3, M4); selection on label-free criteria only + multiverse (stats M5, ai_failure M6); self-association removal (stats M6); CN scale (biology M1, data M1); fusions list (biology M3, data M4); corrected controls (biology M4, ai_failure M9); locus flags (biology M5); registry + ledger (ai_failure M8); NaN/library inclusion rule (data M3); PatientID blocking (stats S3). Partially accepted: escrow of test files by the user replaced by marker + intent commit + hashes reported to user in chat (procedural guarantee stated honestly). Algorithms review: dependency-bearing gene filter, persistence hierarchy instead of single-parameter stability maximisation, (z0,k) by split-half recovery, weight cap, top-down nesting, consensus on edge union, Infomap/MCL sensitivity, bait-held-out BioPlex nulls (all accepted). Deferred: OmicsSomaticMutations long format (not requested; damaging==2 used for biallelic; limitation).
