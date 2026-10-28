# LITERATURE.md: Part 2 (complexes and cancer-context complex dependencies) and Part 1 test evaluation

Compiled 2026-10-28 by the literature agent. Scope is fixed by PROMPT.md. Part 1's review, `external/part1_locked/LITERATURE.md`, already covers co-essentiality network construction, whitening/GLS, DepMap confounders, hu.MAP 3.0 split construction and pair-evaluation pitfalls; it is cited here, not repeated.

**How I verified things**
- Every reference in §9 has had its bibliographic record (title, authors, journal, year, DOI, and PMID where one exists) checked against Crossref or Europe PMC on 2026-10-28.
- Content tags:
  - **[FT]**: the statement was checked in the full text (PMC/Europe PMC).
  - **[ABS]**: checked against the abstract only.
  - **[BIB]**: only the record was verified. Statements about that paper's content are background knowledge; treat them as such.
  - **[unverified]**: a claim I could not check.
- No datasets were downloaded. Data facts quoted from `datasets/` (line counts by lineage etc.) were computed from Model.csv restricted to the 1,361 CRISPRGeneEffect models.

---

## 0. What Part 1 established (compressed)

**The score**
- `zavg_BEF` is label-free. It is the equal-weight mean of z-scored first-order correlation (F), second-order neighbourhood cosine (B) and a Laplacian-eigenmap cosine (E).
- All three are computed on gene-standardised, top-20-PC-removed, ridge-whitened 26Q3 CRISPRGeneEffect profiles.

**Performance**

| set | zavg_BEF | GLS |
|---|---|---|
| dev HT-AP | 0.470 | 0.312 |
| internal final test HT-AP | 0.345 | 0.242 |

- The gain is concentrated in large, dense complexes. On small complexes, macro per-cluster AP is about equal to GLS (0.365 vs 0.369).

**Background from Part 1** (see Part 1 §1–4)
- The co-essentiality → co-complex literature: Pan 2018, Boyle 2018, Kim 2019, Wainberg 2021, Gheorghe & Hart 2022, Rahman 2021 (FLEX), Hassan 2023, Rohde 2025 (BaCoN).
- The mitochondrial dominance of naive networks.
- Library, copy-number and proximity artefacts (Lazar 2024).

**Consequences for Part 2**
1. The pair scores are strongest for large, well-connected complexes. Small complexes and peripheral subunits will be noisier. Any clustering must tolerate that, and complex-size-stratified reporting is needed.
2. The score was built on *whitened, PC-filtered* profiles. That deliberately removes the dominant line-to-line covariance, which includes much lineage structure. **Complex definitions are therefore "pan-context". Context-specific dependency must be quantified on the unwhitened CRISPRGeneEffect scale (Chronos units), not on the whitened profiles used for the pair score.**
3. Because the score is label-free, Part 2 can construct complexes without touching curated labels.

---

## 1. From pairwise scores to complexes

### 1.1 Algorithms

**Markov clustering (MCL)** (Enright 2002 [BIB])
- Flow simulation on a weighted graph, with inflation *I* controlling granularity. It gives a hard partition with no overlaps.
- In benchmark studies on PPI graphs it was the most robust of the methods tested to noise and missing edges:
  - Brohée & van Helden 2006 [BIB]: MCL most robust to edge addition and removal.
  - Vlasblom & Wodak 2009 [BIB]: MCL comparable to or better than affinity propagation on weighted graphs.
- Kim 2019 used MCL on their coessentiality network [FT]:
  - They swept *I* and chose it by GO-BP log-likelihood.
  - Clusters larger than 50 genes, which they note are dominated by OXPHOS and mito-ribosome, were excluded from the selection set.
  - That selection step uses curated annotation and must not be copied here (§5).

**ClusterONE** (Nepusz 2012 [BIB])
- Greedy growth from seeds, maximising a cohesiveness score (in-weight vs boundary weight), followed by merging of overlapping clusters. It allows **overlapping** membership.
- The original paper also introduced the maximum matching ratio (MMR) for comparing predicted complexes to a reference.
- Wainberg 2021 [FT]:
  - Ran ClusterONE on the 17,634² GLS matrix with edge weight 1 − FDR.
  - Used default settings except density *d* ∈ {0.2, 0.5, 0.9}, merged into one list of 5,218–5,229 modules of 4–741 genes.
  - Overlaps were retained deliberately, for example the two MTOR modules corresponding to mTORC1 and mTORC2.

**hu.MAP pipeline** (Drew 2017, 2021; Fischer 2025 [FT])
- Two stages: threshold the pair-confidence network, apply ClusterONE to find dense regions, then apply MCL inside each region.
- Parameter sweep:
  - confidence threshold 0.03–0.95
  - ClusterONE overlap 0.6/0.7 and density 0.1–0.4
  - MCL inflation 1–15
- After clustering, proteins without edges above the threshold are removed and clusters larger than 100 are discarded.
- Parameter sets were selected by k-clique weighted precision/recall against **training gold-standard complexes**. Six sets were chosen manually, from "Extremely high" to "Medium" confidence, and their union is released as a confidence-tiered, nested, overlapping complex map.
- **Lesson:** the two-stage design plus a union of resolutions is a sensible way to obtain nested, overlapping complexes. But its parameter selection is supervised by curated complexes, which we must replace with label-free criteria (§1.3).

**Modularity / Louvain** (Blondel 2008 [BIB]) and **Leiden** (Traag 2019 [BIB])
- Leiden guarantees connected communities; Louvain can return badly connected ones.
- A resolution parameter γ (Reichardt & Bornholdt 2006 [BIB]) or the constant Potts model (CPM; Traag 2011 [BIB]) sets the scale.
- Modularity has a resolution limit (Fortunato & Barthélemy 2007 [BIB]): small complexes merge in large graphs. CPM is "resolution-limit-free" in the sense of Traag 2011, so **CPM with a γ sweep is preferable to modularity** for a 18k-node graph containing complexes of 2–100 members.
- Pan 2018 used Louvain per complex to compare sub-module assignments across networks [FT]. Boyle 2018 used Infomap (igraph `cluster_infomap`) on 1−|r|-weighted co-functional edges, giving 2,857 communities [FT].

**Overlapping / local methods**
- MCODE (Bader & Hogue 2003 [BIB]): seeds and k-cores, unweighted-graph oriented. It also introduced the overlap score ω = |A∩B|²/(|A||B|).
- CFinder / clique percolation (Palla 2005; Adamcsek 2006 [BIB]): k-clique communities. Thresholding a dense weighted graph makes it very sensitive to the threshold.
- Link communities (Ahn 2010 [BIB]): cluster edges, giving multiscale overlap.
- OSLOM (Lancichinetti 2011 [BIB]): significance-based communities, overlapping, with a "homeless" option.
- OCG (Becker 2012 [BIB]): overlapping clusters for multifunctional proteins.
- In the DREAM disease-module challenge (Choobdar 2019 [BIB]), top methods came from diverse families, and no method was uniformly best.

**Hierarchical / multiscale**
- WGCNA-style average linkage on (1 − similarity) plus dynamic tree cut (Langfelder 2008a,b [BIB]).
- HiDeF (Zheng 2021 [BIB]) finds communities that *persist* across a resolution sweep and assembles them into a hierarchy. It underlies the MuSIC and U2OS multimodal maps:
  - Qin 2021 [BIB]: BioPlex + images.
  - Schaffer 2025 [BIB]: multimodal U2OS map.
- Consensus and multiresolution clustering:
  - Lancichinetti & Fortunato 2012 [BIB]: consensus over repeated stochastic runs.
  - Jeub 2018 [BIB]: multiresolution consensus.
  - Monti 2003 [BIB]: resampling consensus.
- shinyDepMap (Shimada 2021 [FT]) used t-SNE + DBSCAN at three density values (d = 65, 100, 141 → 338/608/879 clusters) and linked them into a hierarchy where nested.

**Prior co-essentiality "module" definitions, summarised**

| paper | data | module definition | how granularity was chosen |
|---|---|---|---|
| Pan 2018 [FT] | Avana CRISPR + RNAi | hierarchical clustering (Euclidean, complete linkage) within known complexes; Louvain per complex | per-complex inspection |
| Boyle 2018 [FT] | 436 Avana screens, OR-PC corrected | Infomap communities; community score = PC1 of members | none (Infomap default) |
| Kim 2019 [FT] | BAGEL BF, 276 lines | MCL | GO-BP LLS (curated) |
| Wainberg 2021 [FT] | GLS, 485 lines | ClusterONE, d = 0.2/0.5/0.9 union, overlapping | multiple d, union |
| Shimada 2021 [FT] | CRISPR + shRNA | t-SNE + DBSCAN, 3 densities, hierarchy | multiple densities |
| Chen 2023 [FT] | 1,086 lines, 22Q2 | Pearson matrix cut at r = 0.34; 145 modules of ≤15 genes | threshold tuned on **CORUM** recovery |
| Pan 2022 "Webster" [BIB] | 2,921 high-variance genes | graph-regularised dictionary learning (220 overlapping "functions") | model selection |
| Camacho & Fallahi-Sichani 2026 [ABS] | DepMap | **curated** epigenetic complex annotations | n/a (complexes given) |

### 1.2 Practical implications for converting zavg_BEF into complexes

**Graph construction**
- zavg_BEF is a z-scale score that can be negative. MCL, ClusterONE and CPM need non-negative weights.
- Options:
  - A sparse mutual-kNN or top-q graph, e.g. q = 0.25% to match the B/E construction in Part 1, weighted by a monotone transform such as max(0, z − z₀) or a logistic.
  - Or a global threshold set from a label-free null: the genome-wide random-pair z distribution used for Part 1's z-moments.
- The threshold *is* a hyperparameter. Treat it like resolution and sweep it.

**Recommended primary algorithm**
- Leiden-CPM over a γ grid, with persistence (HiDeF-style) or consensus across 50–100 seeds and across bootstrap resamples of *cell lines*: recompute zavg_BEF on a resample, which costs about 137 s per run per Part 1.
- This gives a nested hierarchy: complex → sub-module.
- Pair it with ClusterONE (overlapping) as a sensitivity analysis, since this is the Wainberg precedent.
- MCL inside persistent regions, as in hu.MAP, is a third option.

**Overlap**
- Hard partitions misassign moonlighting subunits, for example shared SWI/SNF subunits across cBAF, PBAF and ncBAF, or mTORC1/2.
- Options: ClusterONE overlap; post hoc soft membership, where a gene is assigned to every cluster whose mean within-cluster score exceeds a null quantile; or link communities.

**Size handling**
- Cap cluster size, as hu.MAP did with ≤100 and its Complex Portal gold standard with ≤30.
- Recursively re-cluster giants: the OXPHOS/mito-ribosome "mega-module" is expected (Kim 2019; Rahman 2021).

### 1.3 Evaluating clusters without circularity

**Not allowed**
- Selecting resolution or threshold by agreement with CORUM, Complex Portal, hu.MAP, STRING, GO or HGNC gene groups. This would make them construction inputs.
- Chen 2023 tuned on CORUM, Kim 2019 on GO, and hu.MAP on Complex Portal training complexes. All three are fine for their goals but disallowed here.
- Using the withheld hu.MAP test pairs for anything in Part 2.

**Label-free internal criteria, which can select parameters**
- Stability across cell-line bootstraps and seeds: adjusted Rand or overlapping NMI (Lancichinetti 2009 [BIB]); per-cluster Jaccard best-match recovery.
- Persistence across γ (HiDeF).
- CPM quality versus degree-preserving rewired null graphs.
- Split-half replication: build clusters on lines A and recover them on lines B, stratified by lineage.
- Cross-library replication via ScreenGeneEffect (Avana vs KY vs Humagne screens of the same models). This checks technical reproducibility only, not biological independence.

**Independent external evidence, which can validate**
- **BioPlex 3.0 AP-MS edges** as an orthogonal physical assay (§2), and OpenCell if provided. Report the enrichment of within-cluster BioPlex edges against size- and degree-matched random gene sets.
- Caveat: BioPlex 3.0 is an *input feature* of hu.MAP 3.0 (Fischer 2025 [FT]: "includes … ~10k new AP-MS experiments" from BioPlex3). So BioPlex agreement and hu.MAP agreement are not independent of each other. That does not create circularity with our labels as long as hu.MAP is used only post hoc.

**Post hoc only, clearly labelled**
- CORUM 5 (Tsitsiridis 2023), Complex Portal (Meldal 2022), hu.MAP 3.0, STRING (Szklarczyk 2025).
- Metrics: k-clique weighted precision/recall (Drew 2017/2021), MMR (Nepusz 2012), Brohée's Sn/PPV/accuracy, and the ω overlap score.
- **HGNC `gene_group` is literature-curated and contains complex-defined groups.** Examples from our file: "Mediator complex" (26 genes), "BAF complex subunits" (13), "Mitochondrial complex IV" subunits, spliceosomal complexes. It must be treated as post hoc only, like CORUM. HGNC symbols, locations and locus types remain fine for annotation and confounder control (chromosome/arm).

---

## 2. Integrating co-essentiality with experimental PPIs

**Resources**
- BioPlex 3.0 (Huttlin 2021 [ABS/FT]):
  - 293T: 10,128 baits, 118,162 interactions among 14,586 proteins.
  - HCT116: 5,522 IPs (about 70k interactions [unverified exact count]).
  - Shared interactions concentrate in core complexes and essential proteins. Cell-specific interactions *link* complexes, "rewiring" subnetworks. Direct (structurally contacting) pairs replicate far more often than indirect pairs.
  - The paper already correlated DepMap fitness profiles across BioPlex edges [FT]:
    - Positively co-essential edges were enriched for edges shared across both cell lines, and the positive subnetwork contains numerous complexes.
    - Negatively correlated edges matched background and form no complexes.
  - This is the closest prior integration, but it was pan-cancer only and not lineage-resolved.
- Earlier versions: BioPlex 1/2 (Huttlin 2015, 2017 [BIB]).
- OpenCell (Cho 2022 [ABS]): about 1,300 endogenously tagged HEK293T proteins [unverified count] with IP-MS and localisation. Endogenous expression reduces over-expression artefacts.
- Hein 2015 [BIB]: HeLa, stoichiometry-aware IP-MS.
- Go 2021 [BIB]: BioID, a proximity rather than physical map.
- HuRI (Luck 2020 [BIB]): binary Y2H, largely orthogonal to complexes.

**Integration methods in prior work**
- Pan 2018 [FT]: compared co-essentiality sub-modules with BioPlex 2.0, PDB and hu.MAP networks (Louvain per complex, Fisher tests).
- Lord 2020 [FT]: used PPI adjacency to restrict the genetic-interaction search space, which reduced the multiple-testing burden. Robust SL pairs were enriched for physically interacting pairs.
- Pacini 2024 [ABS]: combined dependency–marker associations with PPI networks for target prioritisation.
- Network fusion (SNF, Wang 2014 [BIB]); MuSIC (Qin 2021) for images + AP-MS.
- Supervised integration such as SPOC/Predictomes (Schmid & Walter 2025, Part 1) uses curated labels and is therefore off-limits for construction here.

**Options for us, all label-free**
- (a) **Evidence union with tiers:** a complex is "co-essential + physical" if within-cluster BioPlex density exceeds a null.
- (b) **Edge reweighting:** w = f(zavg) + λ·1[BioPlex edge]. λ must be fixed a priori or chosen by a label-free stability criterion, never by CORUM.
- (c) **Validation only:** keep BioPlex out of construction and use it only as an independent test. This is the cleanest option, and preferable if BioPlex is the main validation source.
- (b) and (c) cannot both be applied to the same edges. If BioPlex is used for construction, validation must come from split-half/cross-library replication and post hoc curated resources.

**Caveats**
- *Cell-type specificity.* 293T is an adrenal/neuronal-like embryonic kidney line, and HCT116 is MSI colorectal. Neither represents most DepMap lineages. Absence of an edge is weak evidence.
- *Bait bias and coverage.* HCT116 has about half the baits of 293T. Bait over-expression (BioPlex uses ORF over-expression) inflates sticky, abundant and chaperone interactions. CompPASS scoring (Sowa 2009 [BIB]) and CRAPome contaminant frequencies (Mellacheruvu 2013 [BIB]) mitigate this but do not eliminate it.
- *Degree bias.* Hubs dominate enrichment statistics. Use degree-preserving nulls.
- *Detection requires expression in the AP cell line.* Lineage-restricted subunits, such as neural or immune complexes, are missing from 293T and HCT116.
- *Indirect co-complex pairs.* BioPlex edges are bait–prey, not pairwise co-membership, so complex-level comparison should use within-cluster edge density rather than pair-by-pair labels.

---

## 3. Context-specific dependency: contexts, examples and statistics

### 3.1 Contexts available in our data

- **Lineage, subtype and fusion.**
  - Model.csv provides `OncotreeLineage`, `OncotreePrimaryDisease`, `OncotreeSubtype`, `OncotreeCode` and `ModelSubtypeFeatures`. The latter is sparse, for example EWS::FLI 14, BCR::ABL1 9, HER2+ 13, MSI 10, MYC_Amplified 21.
  - Line counts in the 1,361 screened models:

| group | lines |
|---|---|
| CNS/Brain | 133 |
| Lung | 126 |
| Lymphoid | 97 |
| Esophagus/Stomach | 96 |
| Bowel | 85 |
| Skin | 75 |
| Pancreas | 71 |
| Breast | 67 |
| Myeloid | 45 |
| AML (code) | 23 |
| NBL | 41 |
| ES (Ewing) | 28 |
| SCLC | 25 |
| PCM (myeloma) | 19 |
| CCRCC | 16 |
| ARMS | 8 |
| SYNS (synovial sarcoma) | 6 |
| MRT (rhabdoid) | 5 |

- **MSI from Model.csv is badly under-ascertained** (10 lines flagged, whereas DepMap contains dozens of MSI lines). Use `OmicsSignatures` (DATA_REQUEST #6).
- **Culture covariates:** `GrowthPattern` (938 adherent, 169 suspension, ...), `OnboardedMedia` (150 formulations) and `ModelCondition` media fields. **CRISPRConfounders** provides screen QC (NNMD, ROC-AUC, FPR, essential/non-essential median depletion and MAD), CasActivity, ScreenDoublingTime, library indicators (Avana, KY, Humagne, TKOv3, Brunello) and 2D/3D screen type.
- **Not yet present:** mutation, copy number, expression, MSI score and fusions (DATA_REQUEST #1–6).

### 3.2 Established complex-relevant context dependencies (positive controls)

**SWI/SNF paralogs and sub-complexes**
- ARID1A-mutant → ARID1B (Helming 2014 [BIB]).
- SMARCA4-deficient → SMARCA2: Hoffman 2014 [BIB]; Wilson 2014 [BIB] on residual SMARCA2 complexes.
- SMARCB1-loss (rhabdoid) and SS18-SSX (synovial sarcoma) → ncBAF/BRD9 (Michel 2018; Brien 2018; Wang 2019 [BIB]).
- SMARCB1-mutant → DCAF5. DCAF5 degrades incompletely assembled SWI/SNF complexes, and the work was done with DepMap screens of 14 SMARCB1-mutant lines (Radko-Juettner 2024 [ABS]). This is a clean "complex quality-control" example.
- POU2F3-SCLC and myeloma → mSWI/SNF ATPases (He 2024; Duplaquet 2024 [ABS/BIB]), an example of a lineage TF complex driving complex dependency.

**Metabolic collateral loss**
- MTAP deletion (9p21, with CDKN2A) → PRMT5 and the MAT2A/PRMT5/RIOK1 axis (Kryukov 2016; Mavrakis 2016; Marjon 2016 [BIB]).
- PRMT5 and its obligate partner WDR77 should co-move as a unit [unverified as a DepMap observation].

**MSI**
- MSI → WRN (Behan 2019 [ABS]; Chan 2019; Lieb 2019; Kategaya 2019 [BIB]), with mechanism via TA-repeat expansions (van Wietmarschen 2020 [BIB]) and clinical-stage inhibitors (Picco 2024 [BIB]).
- **Complex–complex SL:** a destabilised SKI complex (in MSI-H and 9p21-deleted contexts) → dependence on the PELO–HBS1L ribosome-rescue complex (Prindle 2025 [ABS]). This is the best published template for "context-dependent dependency between two complexes".

**CYCLOPS / partial copy loss**
- Nijhawan 2012 [ABS] (RNAi) found 56 genes whose suppression selectively kills lines with *partial* copy loss. They are enriched for **spliceosome, proteasome and ribosome** components.
- Mechanism for PSMC2: normal cells keep excess PSMC2 in a PSMC1/PSMD2/PSMD5 reservoir complex that is absent after partial loss.
- POLR2A (17p, co-deleted with TP53) is a related example (Liu 2015 [BIB]).
- VPS4A dependency with 18q/16q loss of VPS4B (Neggers 2020 [BIB]) is a paralog-CYCLOPS hybrid.
- For complex members, collateral loss is often visible in protein rather than mRNA (Ryan 2017 [ABS]). Aneuploid cancers show extensive protein dosage compensation, strongest for complex subunits (Schukken & Sheltzer 2022 [BIB]).
- Amplification-related gain of sensitivity (ARGOS; Rendo 2025 [ABS]) is the amplification analogue.

**Paralog buffering**
- Paralogs are less often essential. Pairs from whole-genome duplication **and pairs that function in protein complexes** are significantly more likely to be synthetic lethal (De Kegel & Ryan 2019 [ABS]).
- Shared PPIs and conservation predict paralog SL (De Kegel 2021 [ABS]).
- Combinatorial screens: Dede 2020 (enCas12a), Thompson 2021, Parrish 2021, Ito 2021 (DUSP4/6) [BIB].
- Cohesin STAG1/STAG2 (van der Lelij 2017 [BIB]); EP300/CREBBP (Ogiwara 2016 [BIB]).
- Incomplete paralog compensation, e.g. TRA2A (Lee 2025 [BIB]).
- In yeast, paralogs that heteromerise can impart *fragility* through dependency (Diss 2017 [BIB]). In human, there is evidence that paralog dependency indirectly shapes robustness (Dandage & Landry 2019 [BIB]).

**KEAP1/NRF2**
- KEAP1 loss → NRF2 pathway activation and glutaminolysis dependence (Romero 2017 [BIB]) and NRF2-dependent druggable vulnerabilities (Bar-Peled 2017 [BIB]).
- NFE2L2 dependency in KEAP1-mutant lung lines is a widely used DepMap positive control [unverified for 26Q3]. The CUL3–KEAP1 E3 ligase is the relevant complex.

**VHL/HIF**
- HIF-2α dependence of VHL-null ccRCC (Chen 2016 [BIB]).
- HIF-independent CDK4/6 synthetic lethality with VHL loss (Nicholson 2019 [BIB]).
- Expect weak EPAS1 signal in 2D CRISPR screens [unverified]. Only 16 CCRCC lines are available.

**Lineage TF / co-activator complexes**
- Neuroblastoma core regulatory circuitry (Durbin 2018 [BIB]).
- SAGA KAT module in MYCN-amplified neuroblastoma (Malone 2024 [ABS]).
- Mediator kinase:
  - AML (Pelish 2015 [BIB]).
  - Fusion-positive alveolar rhabdomyosarcoma (Zhang 2026 [BIB]); only 8 ARMS lines are available.
- MLL / menin / ENL / KAT6A in MLL-rearranged leukaemia (Krivtsov 2019; Wan 2017; Erb 2017; Yan 2022 [BIB]).
- PRC2 × MLL-MEN1 in EZH2-mutant lymphoma (Chen 2023 [ABS]).
- Set1C/COMPASS (CXXC1) in melanoma (Camacho 2026 [ABS]).

**Mitochondria / OXPHOS by lineage**
- AML sensitivity to OXPHOS and mitochondrial translation inhibition (Molina 2018; Škrtić 2011 [BIB]).
- OXPHOS essentiality depends on medium and carbon source (Arroyo 2016; Rossiter 2021 [BIB]).
- The OXPHOS signal in DepMap tracks screen duration and doubling time (Rahman 2021, Part 1).

**Other complex-level discoveries from DepMap**
- XPR1–KIDINS220 phosphate-export complex in ovarian cancer (Bondeson 2022 [BIB]).
- SKA complex in chromosomally unstable cancers (Uiterkamp 2026, bioRxiv [ABS-title only, unverified]).
- Aneuploidy → spindle assembly checkpoint/KIF18A (Cohen-Sharir 2021 [BIB]).
- Pan-cancer proteasome subunits PSMB5/6 as context biomarkers (Kim JH 2026, *Molecules* [ABS-title only]; low evidential weight).

### 3.3 Statistical practice in biomarker/context discovery

**Data and models**
- **Gene-effect estimation:** CERES (Meyers 2017 [BIB]) then Chronos (Dempster 2021 [BIB]), which corrects copy-number/cutting toxicity and screen quality. The DepMap pipeline additionally regresses CN post hoc and aligns arm medians (Part 1 §4).
- Copy-number artefact origin: Aguirre 2016; Munoz 2016 [BIB].
- Proximity/arm-truncation artefacts: Lazar 2024 [BIB].
- **Linear models with tissue and MSI covariates.** Lord 2020 [FT] used multiple regression with *tissue type, MSI status and driver status* as covariates, explicitly to avoid confounding such as SOX10 dependency appearing BRAF-associated because melanoma is BRAF-mutant.
  - A driver was tested only if altered in ≥5 lines; BH FDR.
  - Associations had to replicate across **independent datasets and non-overlapping cell-line panels**.
  - PPI adjacency was used to restrict the passenger-gene search.
- Behan 2019 and Pacini 2024 used ANOVA-style tissue-stratified marker tests within a prioritisation framework [ABS; details unverified].
- Tsherniak 2017 [BIB] introduced DepMap-scale predictive modelling of dependencies.
- Dempster 2019 (*Nat Commun* [BIB]) showed Broad/Sanger concordance of biomarker associations. Dempster 2019 (bioRxiv 720243 [BIB]) describes Achilles screen-quality metrics.
- **Module-level lineage tests:** Wainberg 2021 [FT] ran per-gene GLS with 20 cancer-type indicators (an ANOVA under GLS) and aggregated gene p-values per module with the **Cauchy combination test / ACAT** (Liu & Xie 2020; Liu 2019 [BIB]). 442 modules were selectively essential in 16 cancer types.
  - Caveat: because GLS whitens with the line×line covariance, which encodes lineage similarity, it down-weights lineage-shared variance. It is a conservative choice for lineage tests.
- Kim 2019 [FT]: cluster essentiality = mean BF per line; cases are lines with median cluster BF > 5; then differential expression, CN and mutation tests.
- Boyle 2018 [FT]: community score = PC1, compared by lineage and genotype.
- Chen 2023 [FT]: module dependency vs cancer type (FDR < 0.1, 137 cancer-type-specific modules), then drug sensitivity (GDSC).

**Differential co-essentiality** (Kim 2022 [FT])
- dPCC = PCC(all lines) − PCC(all lines except context).
- Significance comes from 1,000 random removals of the same number of lines, giving an empirical p-value at P ≤ 0.001 (null distributions are non-Gaussian).
- They found 30 features with gain and 133 with loss of interaction, for example BRAF-mutant loss of BRAF–RAF1.
- They note that small contexts make within-context networks meaningless, citing 7 CML lines as an example.

**Multiple testing and permutation**
- BH (Benjamini & Hochberg 1995 [BIB]).
- Stage-wise complex → subunit testing (stageR, Van den Berge 2017 [BIB]) suits "complex first, then which subunit".
- Permutation with nuisance covariates: Freedman–Lane or within-lineage label shuffles (Winkler 2014 [BIB]).
- limma-style empirical-Bayes moderation (Ritchie 2015 [BIB]) stabilises variance for many complexes × contexts.

**Robustness, power and confounders**
- Robustness across datasets and libraries: Pacini 2021; Vinceti 2024 (Part 1).
- **Power with few mutant lines** (my calculation, not from a paper): two-sample z approximation, n₁ altered vs about 1,000 others, two-sided α = 10⁻⁴:

| altered lines | Cohen's d | power |
|---|---|---|
| 10 | 1.0 | ≈ 0.23 |
| 10 | 1.5 | ≈ 0.80 |
| 25 | 1.0 | ≈ 0.86 |

  - Averaging *k* subunits whose noise correlates at ρ_ε reduces noise variance by a factor (1+(k−1)ρ_ε)/k. That is the statistical case for complex-level scores, provided the subunits share the context effect.
- **Confounder subtleties**
  - Screen-quality metrics such as NNMD and essential-median depletion are *functions of essential-gene depletion*. Adjusting for them partly adjusts away genuine variation in core-essential complexes (proteasome, ribosome, spliceosome). Use them as a sensitivity analysis, or use within-line rank normalisation, rather than as unconditional covariates [reasoning; no published source found].
  - ScreenDoublingTime correlates with OXPHOS/mito-ribosome dependency (Rahman 2021).
  - Suspension vs adherent growth and medium are near-collinear with lymphoid/myeloid lineage.
  - Library (Avana/KY/Humagne) composition differs by lineage because Sanger KY screens over-represent some tissues [unverified magnitude].

---

## 4. Complex-level dependency quantification and rewiring

**Score options** (applied to unwhitened Chronos gene effects)
1. **Mean or trimmed-mean gene effect** of members. It stays on the Chronos scale (−1 is the median common essential), so it is interpretable. It is dominated by the most essential subunits. Kim 2019 used the mean BF [FT].
2. **First PC / eigengene** of standardised member profiles, sign-aligned to the mean (Boyle 2018 [FT]; Langfelder & Horvath 2007 [BIB]; PLAGE, Tomfohr 2005 [BIB]). It is dominated by high-variance subunits and loses the absolute scale.
3. **Rank-based single-sample scores:** singscore (Foroutan 2018 [BIB]), GSVA (Hänzelmann 2013 [BIB]) and ssGSEA (Barbie 2009 [BIB]). Within-line ranking removes per-screen offsets and quality scaling. Note that these were designed for expression, and gene-set size and membership overlap matter.
4. **P-value aggregation** across members (Wainberg's ACAT) tests "any member" rather than "the complex as a unit".
5. **Multivariate complex-level analysis:** Camacho 2026 [ABS] analysed curated chromatin complexes multivariately across lineages.

**Recommendations**
- Report (1) as primary for interpretability and (2) or (3) as sensitivity analyses.
- Require leave-one-subunit-out stability, so that a "complex" effect is not one subunit.
- Weight members by membership confidence.
- Exclude or flag members that are non-expressed (once expression arrives) or common-essential with near-zero variance.

**Within-complex divergence**
- A complex is not always a dependency unit:
  - modular sub-complexes: Mediator kinase vs core, SAGA KAT vs core, cBAF/PBAF/ncBAF;
  - paralog alternatives: ARID1A/B, SMARCA2/4, STAG1/2;
  - quality-control relationships: DCAF5.
- Pan 2018 [FT] showed that co-essentiality resolves BAF into core BAF, PBAF and a novel module (SMARCD1, BRD9, GLTSCR1), which was later defined as ncBAF.
- The test is: subunit effect minus complex score (residual) versus context, with a stage-wise FDR.

**Differential co-essentiality / rewiring**
- Prior work:
  - Gene-pair level: Kim 2022 dPCC [FT]; FIREWORKS (Amici 2021 [ABS]) for context subsets.
  - Module level by context: Chen 2023.
  - Billmann 2026 *Cell* [ABS]: a HAP1 genetic-interaction map complements DepMap co-essentiality and identifies determinants of cell-line-specific dependencies.
  - General methods: Ideker & Krogan 2012 [BIB]; DiffCoEx (Tesson 2010 [BIB]); DINGO (Ha 2015 [BIB]); direct differential-network estimation (Zhao 2014 [BIB]).
- **Key pitfall:** correlation depends on signal variance. A complex that is uniformly essential within a lineage, or uniformly dispensable, shows low within-lineage correlation even if physically intact. "Loss of co-essentiality" must therefore be tested conditional on subunit variance. Use variance-matched permutation, or compare regression slopes and covariances rather than correlations.
- The Fisher-z standard error is 1/√(n−3): about 0.19 at n = 30 lines. Lineage-level pair rewiring is underpowered except in the largest lineages. Complex-level coherence (mean pairwise r over k members) helps but is not independent across pairs.

---

## 5. Circularity: what must not be used, and how prior work leaked

**Why these resources must stay out of construction and validation**
- STRING aggregates curated complex and pathway databases, text mining, co-expression, and experiments including AP-MS (Szklarczyk 2023, 2025 [BIB/ABS]).
- hu.MAP 3.0 is trained on Complex Portal and built from AP-MS (including BioPlex), CF-MS and proximity data (Fischer 2025 [FT]).
- CORUM and Complex Portal are literature-curated. HGNC gene groups include curated complexes (§1.3).
- Our Part 1 labels (hu.MAP train/test pairs) derive from Complex Portal.
- Using any of these to choose clustering parameters, define complexes, weight edges or judge "validity" makes the evaluation partly a re-discovery of the curation. The withheld hu.MAP test must stay untouched by Part 2.

**Prior leakage patterns**

| study | what was tuned or evaluated against curated data |
|---|---|
| Chen 2023 [FT] | network cutoff tuned on CORUM recovery |
| Kim 2019 [FT] | MCL inflation chosen by GO-BP enrichment |
| Gheorghe & Hart 2022 (Part 1) | threshold set at hu.MAP 2.0 precision; evaluated on Reactome |
| Wainberg 2021 | GLS benchmarked against CORUM, hu.MAP and STRING (legitimate for benchmarking, but those modules cannot then be "validated" by the same sources) |
| Pacini 2024; Lord 2020 | PPI networks used for prioritisation (fine, and disclosed) |
| hu.MAP | parameters chosen on training complexes, evaluated on held-out complexes (a proper split, but still curated) |

**Post hoc interpretation practice**
1. Freeze complexes, scores and context results (commit hash) before loading CORUM, Complex Portal, hu.MAP, STRING or HGNC gene groups.
2. Use them only to *name* complexes, to measure novelty (fraction of a cluster's pairs absent from all curated resources) and for descriptive concordance.
3. Never feed post hoc results back into thresholds.
4. Report concordance as description, not validation, and keep a "discovery evidence" column (co-essentiality, BioPlex, replication) separate from an "annotation" column.

---

## 6. Brief: final held-out test evaluation of a ranking score on dependent pairs (Part 1 hu.MAP test)

**Metric**
- AP (average precision) is the primary metric, as locked in Part 1. Boyd 2013 [BIB]: AP is a sound estimator. Interpolated areas must not be used.
- PR is more informative than ROC under imbalance (Saito 2015 [BIB]; Davis & Goadrich 2006 [BIB]).
- The test files have 3,312 positives and 597,426 negatives, a prevalence of about 0.55%. This is the natural prevalence *of the hu.MAP construction* (negatives are cross-complex pairs among complex members), not genome-wide.
- Report AP at this prevalence as primary. Optionally also report precision re-weighted to other prevalences via π-adjustment (precision = TPR·π / (TPR·π + FPR·(1−π))), labelled as a derived quantity.
- If Part 1's locked metric is HT-weighted AP, apply exactly that definition; do not switch.

**Strata (pre-registered)**
- C1/C2/C3 protein overlap with Part 1 training proteins (Park & Marcotte 2012, Part 1).
- Complex-size bins, essentiality strata, and the largest component excluded.
- zavg_BEF is label-free, but its design was *selected* on train labels, so overlap strata still matter.

**Uncertainty**
- Pairs share proteins, so i.i.d. bootstrap CIs are too narrow.
- Use a **cluster bootstrap over connected components of the positive graph** (Field & Welsh 2007 [BIB]), which matches the complex-level hu.MAP split.
- As a sensitivity analysis, use a **two-way / pigeonhole bootstrap over proteins** for dyadic data (Owen 2007 [BIB]; Davezies 2021 [BIB]; multiway clustering, Cameron 2011 [BIB]).
- Compute method differences with *paired* resampling: the same resample for zavg_BEF, GLS, Pearson and Spearman. Use 2,000 or more replicates with fixed seeds.

**Calibration**
- zavg_BEF and the baselines are not probabilities.
- Calibrate with a monotone map such as isotonic regression (Zadrozny & Elkan 2002 [BIB]) or Platt-type logistic (Niculescu-Mizil & Caruana 2005 [BIB]). The map must be **fitted on Part 1 development data only** and then applied to test scores.
- Report reliability diagrams, the Brier score, calibration intercept and slope (Van Calster 2019 [BIB]), and the prevalence mismatch between dev and test.
- AP and AUROC are invariant to monotone calibration.

**Thresholded metrics**
- Precision, recall and F1 need a threshold chosen **on development data only**: the dev F1-optimal threshold, or the threshold giving a pre-stated dev precision. The test set's own F1-optimal threshold may be shown only as a labelled oracle.

**Baselines and protocol**
- Pearson, Spearman and GLS must be scored on exactly the same pairs with the same unmappable-pair rule.
- Commit the protocol before scoring. Score once.

---

## 7. Candidate research directions (ranked by interpretability × robustness × novelty)

Common elements apply to every direction:
- complexes are frozen from zavg_BEF clustering (§1);
- complex scores are on raw Chronos (§4);
- the model is score ~ context + lineage (or within-lineage) + library + sensitivity covariates;
- inference uses within-stratum permutation, BH plus stage-wise FDR, and **discovery/validation on lineage-stratified, non-overlapping halves of the cell lines** (Lord 2020 style);
- negative controls are size- and variance-matched random gene sets, the same complexes with shuffled context labels, and olfactory-receptor or non-expressed gene sets;
- positive controls are pre-registered before any result is seen.

**Data legend**
- **[P]** runs on data present now (CRISPRGeneEffect, ScreenGeneEffect, CRISPRConfounders, Model, ModelCondition, HGNC).
- **[R]** needs requested data (DATA_REQUEST.md).
- **[X]** needs extra data (§8).

### D1. Lineage/subtype complex-dependency atlas with explicit confounder control [P] — rank 1

**Hypothesis**
- Data-driven complexes show lineage- and subtype-selective dependency that is not explained by culture or screen confounders.
- Selective complexes are enriched for transcriptional co-activator and chromatin complexes, while core-machine complexes vary mainly with quality or doubling time.

**Data:** CRISPRGeneEffect, Model (Oncotree lineage, code and subtype), CRISPRConfounders, ModelCondition (medium).

**Approach**
- Complex score (mean, eigengene, singscore) ~ one-vs-rest lineage indicators.
- Sensitivity models add library, ScreenDoublingTime and GrowthPattern, with NNMD only as a sensitivity covariate (§3.3).
- Run nested tests within lineage at OncotreeCode level, e.g. AML vs other myeloid, SCLC vs NSCLC.
- Replicate (i) across split halves of lines and (ii) across libraries via ScreenGeneEffect for multiply-screened models.

**Positive controls**
- Mediator kinase, MLL/menin, ENL/KAT6A in AML/myeloid.
- ncBAF/BRD9 in synovial sarcoma (6 lines; descriptive) and rhabdoid.
- SAGA KAT module in neuroblastoma.
- mSWI/SNF in myeloma.
- Lineage TF complexes, e.g. SOX10/melanoma modules (Wainberg 2021 found these).
- OXPHOS in AML.

**Confounders:** suspension/medium collinear with heme lineages; library-by-lineage imbalance; doubling time.

**Novelty vs prior work**
- Wainberg 2021: 485 lines, ClusterONE modules, GLS+ACAT, cancer-type level.
- Chen 2023: 145 small modules, CORUM-tuned cutoff.
- Camacho 2026: curated chromatin complexes only.
- What is new here: about 2.8× more lines, curation-free complexes from a better co-complex predictor, subtype (OncotreeCode) resolution, explicit confounder decomposition and replication.
- Novelty is moderate, but this is the foundation for D2–D5.

**Risk:** low. Main risk is "rediscovering lineage TFs"; mitigate by reporting complexes rather than single genes and by stratifying by complex class.

### D2. Sub-module divergence: which part of a complex carries the context dependency [P] — rank 2

**Hypothesis:** for modular complexes, context selectivity resides in specific sub-modules or subunits while the core stays pan-essential. Examples: Mediator kinase module vs core, SAGA KAT vs core, cBAF/PBAF/ncBAF, Integrator, CCR4-NOT, COMPASS family, PRC1 variants.

**Approach**
- Use the multiscale hierarchy (§1.2) to define parent complexes and child sub-modules.
- Test (child score − parent score), or subunit residuals, against lineage and subtype, with a stage-wise FDR (complex, then sub-module, then subunit).
- Quantify "selectivity heterogeneity" per complex as the variance of member context effects.

**Positive controls**
- SAGA KAT module in MYCN-amplified neuroblastoma (Malone 2024).
- Mediator kinase in AML and ARMS (Pelish 2015; Zhang 2026).
- ncBAF (BRD9/GLTSCR1) vs cBAF in synovial sarcoma and rhabdoid (Michel 2018; Brien 2018).
- Set1C/CXXC1 vs other COMPASS in melanoma (Camacho 2026).

**Confounders**
- Sub-modules built from the same data that defines selectivity. Mitigate by building the hierarchy on half A and testing on half B.
- Differences in essentiality floor between subunits (Chronos saturation near −1 to −2).

**Novelty and data**
- Pan 2018 showed BAF sub-modules qualitatively. A systematic, hierarchy-wide, context-resolved sub-module analysis has not been published as far as I could find (searches: Europe PMC "complex-level", "sub-module", "differential co-essentiality").
- Highly interpretable, and present data suffice.

**Risk:** moderate. Small sub-modules (2–3 genes) are where zavg_BEF is weakest (Part 1 macro-AP ≈ GLS), so BioPlex support would help.

### D3. Genotype-conditioned complex dependencies, including complex-to-complex synthetic lethality [R] — rank 3

**Hypothesis:** damaging mutation or deletion of a subunit of complex A creates dependency on (i) its paralog in A, (ii) a sister sub-complex, or (iii) a functionally coupled complex B. Examples: SKI → PELO–HBS1L; SMARCB1 → DCAF5/ncBAF; KEAP1 → NRF2 axis; MTAP → PRMT5–WDR77/MAT2A.

**Data:** mutation (damaging/hotspot), CN deletions, OmicsSignatures MSI, Model lineage. Optional: fusions [X].

**Approach**
- Collapse genotypes to the *complex* level: "complex A disrupted" if any core subunit has a damaging alteration. This aggregates rare alterations and raises power, but report the subunit level too.
- Test all complex-A-disrupted × complex-B-score pairs with lineage fixed effects, MSI and library covariates, within-lineage permutation and split-half replication.
- The hypothesis space is restricted a priori to A×B pairs connected in zavg_BEF or BioPlex, analogous to Lord 2020's PPI restriction.

**Positive controls:** ARID1A→ARID1B; SMARCA4→SMARCA2; SMARCB1→BRD9/DCAF5; MTAP-del→PRMT5/WDR77/MAT2A; MSI→WRN; MSI/9p21→PELO/HBS1L; KEAP1→NFE2L2; STAG2→STAG1; CREBBP→EP300.

**Confounders**
- Lineage–genotype collinearity (e.g. SMARCA4 in NSCLC).
- Hypermutation in MSI lines inflates "damaging" calls, so require recurrent or truncating alterations and adjust for MSI and mutation burden.
- CN correction for deletions.

**Novelty:** gene-level SL screens are extensive (Behan, Pacini, Lord, Kim 2022, De Kegel). A systematic **complex×complex** genetic-interaction map with curation-free complexes is new. Prindle 2025 is a single-case precedent.

**Risk:** moderate, and it depends on data arrival.

### D4. Paralog switching within complexes [P partial / R / X] — rank 4

**Hypothesis:** for complexes containing paralog alternatives, dependency switches between paralogs across lineages or genotypes. The switch is predicted by the *other* paralog's loss, low expression or low dependency, and it leaves the complex score constant while subunit identity changes.

**Data**
- Present data allow a dependency-only version: anti-correlated paralog dependencies within a complex, with paralogs approximated from HGNC gene families. HGNC families are curated gene families rather than complexes, but are still literature-derived; prefer sequence-based paralogs.
- Full version: expression and mutation [R], plus Ensembl paralogue pairs or De Kegel 2021 Table S1 [X].

**Approach**
- For in-complex paralog pairs (both genes in the same frozen complex), test the interaction dep(P1) ~ status(P2) × lineage.
- Run a "switch test": negative within-complex residual correlation between paralogs alongside a positive complex-level correlation.

**Positive controls:** ARID1A/B, SMARCA2/4, STAG1/2, EP300/CREBBP, VPS4A/B (complex: ESCRT-III/VPS4), DUSP4/6 (not a complex; negative-class control).

**Confounders:** paralog sgRNA cross-reactivity (shared guides); same-arm proximity for tandem paralogs; expression-driven non-dependency (a non-expressed gene is trivially dispensable).

**Novelty:** paralog SL is well studied gene-pair-wise (De Kegel 2019, 2021; Dede 2020; Thompson 2021; Parrish 2021). Placing it inside data-driven complexes and across lineage contexts, and asking which complexes tolerate subunit substitution, is new.

**Risk:** moderate. Needs paralog annotation; De Kegel 2019 already reported that in-complex paralogs are SL-enriched.

### D5. Dosage and copy-number vulnerabilities at complex level: CYCLOPS revisited in CRISPR [R] — rank 5

**Hypothesis:** partial copy loss of a subunit increases dependency on that subunit or on its complex partners, and this is strongest for stoichiometric "machine" complexes (proteasome, spliceosome, ribosome, RNA Pol II), as in CYCLOPS. Amplification may produce the reverse (ARGOS).

**Data:** gene-level relative CN [R]. Ideally absolute CN and ploidy plus arm-level calls [X], and expression [R]. Proteomics [X] would test the collateral-loss mechanism.

**Approach**
- dep(g) ~ CN(g) + lineage, compared against a **segment-matched null**: genes in the same CN segment that are not in the complex, which controls residual cutting and correction artefacts and proximity bias.
- At complex level: dependency of partner subunits ~ CN of the lost subunit.

**Positive controls:** PSMC2, SF3B1 and other CYCLOPS genes (Nijhawan 2012); POLR2A with 17p loss (Liu 2015); VPS4A with 18q loss (Neggers 2020).

**Confounders**
- Chronos and post hoc CN correction can absorb or create CN–dependency slopes.
- Arm-level truncation (Lazar 2024).
- Ploidy and whole-genome doubling.

**Novelty:** CYCLOPS was RNAi-based (2012). A CRISPR-era, complex-aware re-analysis with segment-matched nulls has not been done to my knowledge [unverified absence].

**Risk:** moderate-high, because the CN correction is in the pipeline itself. DEMETER2 RNAi [X] gives an orthogonal check, since partial knockdown is closer to partial loss.

### D6. Context-specific rewiring of complex coherence (differential co-essentiality at complex level) [P] — rank 6

**Hypothesis**
- Some complexes are coherent dependency units in some lineages but fragment in others, which would indicate cell-type-specific composition.
- Some gain new members: lineage-specific partners, such as FURIN/CPD with IGF1R/MET (Kim 2022).

**Approach**
- Within-lineage coherence = mean pairwise covariance or slope among members, on Chronos with lineage-mean centring.
- Contrast with a variance-matched permutation null (§4 pitfall) and Kim 2022's random-removal bootstrap.
- Restrict to lineages with ≥40 lines: CNS, Lung, Lymphoid, Esophagus/Stomach, Bowel, Head & Neck, Skin, Pancreas, Breast, Ovary, Bone, Soft tissue, PNS, Myeloid.
- Complexes come from pan-cancer clustering, so selection is label-free and label permutation is valid.

**Positive controls:** BRAF–RAF1 loss of interaction in BRAF-mutant lines (needs mutation; Kim 2022); KRAS/NRAS–RAF1; lineage-restricted complexes expected to fragment where not expressed.

**Confounders:** variance differences; lineage size; screen quality differing by lineage.

**Novelty:** gene-pair dPCC exists (Kim 2022; FIREWORKS), so this extends it to complex-level coherence with variance control.

**Risk:** high, because power is low (SE of r ≈ 0.19 at n = 30) and variance confounding is strong. Ranked lower for robustness.

### D7. Cell-type-specific interactome rewiring vs lineage-specific co-essentiality (BioPlex 293T vs HCT116) [R] — rank 7

**Hypothesis**
- HCT116-specific BioPlex interactions are more co-essential in colorectal (Bowel, 85 lines; MSI-high subset) than in other lineages, and 293T-specific edges are more co-essential in kidney or neural lineages.
- Shared edges are pan-lineage co-essential.

**Approach**
- Restrict to baits screened in both lines and to preys expressed in both [R expression].
- Compare lineage-specific co-essentiality (raw-profile partial correlation within lineage) for shared, 293T-only and HCT116-only edges.
- Use degree-preserving nulls and logistic or ordinal regression on edge class.

**Positive controls:** shared core complexes are co-essential everywhere (Huttlin 2021 pan-cancer result).

**Confounders:** 293T lineage mismatch; detection false negatives (much of the apparent "cell specificity" is technical, per Huttlin 2021); bait over-expression.

**Novelty:** Huttlin 2021 did only a pan-cancer co-essentiality comparison. The lineage-matched test is new.

**Risk:** high. The effect may be small, and the design depends on BioPlex arriving.

### D8. Mitochondrial complexes: lineage vs culture determinants [P] — rank 8

**Hypothesis:** dependence on respiratory-chain complexes I–V, the mito-ribosome and mitochondrial import varies by lineage beyond doubling time, assay length and medium. Heme lineages remain OXPHOS-dependent after adjustment.

**Approach**
- Complex scores for each OXPHOS module from the hierarchy, modelled on lineage + ScreenDoublingTime + medium (ModelCondition/Model formulation) + growth pattern + library.
- Use variance decomposition (partial R²), and matched comparisons within identical media.

**Positive controls:** AML sensitivity (Molina 2018; Škrtić 2011); Rahman 2021's doubling-time association.

**Confounders:** non-identifiability (lineage ≈ medium ≈ growth pattern).

**Novelty:** low-moderate. It is useful mainly as a confounder-characterisation chapter that protects D1 from mitochondrial dominance.

**Risk:** low technical, high interpretive.

---

## 8. External datasets that would substantially help (not in DATA_REQUEST.md; do not download—PI to request)

| # | Dataset (exact file, as named on source) | Why | Directions |
|---|---|---|---|
| E1 | DepMap 26Q3 `OmicsFusionFiltered.csv` | Fusion contexts: SS18::SSX (synovial sarcoma → ncBAF), EWSR1::FLI1, KMT2A-r (MLL/menin/ENL), PAX3/7::FOXO1 (ARMS → Mediator kinase). Model.csv `ModelSubtypeFeatures` is sparse. | D1–D3 |
| E2 | DepMap 26Q3 `OmicsAbsoluteCNGene.csv` + `OmicsArmLevelCNA.csv` (names per recent releases; verify in the 26Q3 file list) | Homozygous vs partial loss, ploidy/WGD, arm-level loss (VPS4B 18q, 17p/POLR2A) | D3, D5 |
| E3 | DepMap RNAi DEMETER2 v6 `D2_combined_gene_dep_scores.csv` (McFarland 2018) | Orthogonal perturbation for replication of complex-context findings; partial knockdown suits CYCLOPS | D1, D3, D5 |
| E4 | DepMap 26Q3 `CRISPRInferredCommonEssentials.csv`, `AchillesCommonEssentialControls.csv`, `AchillesNonessentialControls.csv` | Label-free classes for screen-quality-aware modelling and negative controls (not complex-curated) | all |
| E5 | Ensembl BioMart human paralogue export (gene ID, paralogue ID, % identity, Ensembl release ≥110), *or* De Kegel et al. 2021 Table S1 (paralog pairs with features) | Sequence-based paralog definitions, avoiding curated HGNC families | D4 |
| E6 | Combinatorial paralog-screen supplementary tables: Dede 2020, Thompson 2021, Parrish 2021, Ito 2021 | Independent post hoc validation of paralog switching | D4 |
| E7 | Cancer cell-line proteomics: CCLE (Nusinow 2020, `protein_quant_current_normalized.csv`) and/or Sanger ProCan-DepMapSanger protein matrix (Gonçalves 2022; exact file name to confirm on cellmodelpassports/PRIDE) | Subunit co-abundance and collateral loss (Ryan 2017; Schukken 2022); protein-level context for D5 | D3, D5 |
| E8 | OpenCell interaction table (`opencell-protein-interactions` download, opencell.czbiohub.org) | Second AP-MS source with endogenous tagging (less over-expression bias) | §2, D7 |
| E9 | PRISM Repurposing secondary/extended viability matrix (Corsello 2020; DepMap portal "PRISM Repurposing Public") | Post hoc pharmacological support (BRD9, menin, PRMT5/MAT2A, WRN, KAT2A/B degraders) | D1–D3, post hoc |

---

## 9. References (all bibliographically verified 2026-10-28 via Crossref/Europe PMC unless marked)

**Clustering and complex detection**
1. Enright AJ, van Dongen S, Ouzounis CA. An efficient algorithm for large-scale detection of protein families. *Nucleic Acids Res* 2002;30:1575–1584. doi:10.1093/nar/30.7.1575
2. Nepusz T, Yu H, Paccanaro A. Detecting overlapping protein complexes in protein-protein interaction networks. *Nat Methods* 2012;9:471–472. doi:10.1038/nmeth.1938
3. Blondel VD, Guillaume J-L, Lambiotte R, Lefebvre E. Fast unfolding of communities in large networks. *J Stat Mech* 2008:P10008. doi:10.1088/1742-5468/2008/10/P10008
4. Traag VA, Waltman L, van Eck NJ. From Louvain to Leiden: guaranteeing well-connected communities. *Sci Rep* 2019;9:5233. doi:10.1038/s41598-019-41695-z
5. Traag VA, Van Dooren P, Nesterov Y. Narrow scope for resolution-limit-free community detection. *Phys Rev E* 2011;84:016114. doi:10.1103/PhysRevE.84.016114
6. Fortunato S, Barthélemy M. Resolution limit in community detection. *PNAS* 2007;104:36–41. doi:10.1073/pnas.0605965104
7. Reichardt J, Bornholdt S. Statistical mechanics of community detection. *Phys Rev E* 2006;74:016110. doi:10.1103/PhysRevE.74.016110
8. Bader GD, Hogue CWV. An automated method for finding molecular complexes in large protein interaction networks. *BMC Bioinformatics* 2003;4:2. doi:10.1186/1471-2105-4-2
9. Palla G, Derényi I, Farkas I, Vicsek T. Uncovering the overlapping community structure of complex networks in nature and society. *Nature* 2005;435:814–818. doi:10.1038/nature03607
10. Adamcsek B, et al., Vicsek T. CFinder: locating cliques and overlapping modules in biological networks. *Bioinformatics* 2006;22:1021–1023. doi:10.1093/bioinformatics/btl039
11. Brohée S, van Helden J. Evaluation of clustering algorithms for protein-protein interaction networks. *BMC Bioinformatics* 2006;7:488. doi:10.1186/1471-2105-7-488
12. Vlasblom J, Wodak SJ. Markov clustering versus affinity propagation for the partitioning of protein interaction graphs. *BMC Bioinformatics* 2009;10:99. doi:10.1186/1471-2105-10-99
13. Choobdar S, et al., Marbach D. Assessment of network module identification across complex diseases. *Nat Methods* 2019;16:843–852. doi:10.1038/s41592-019-0509-5
14. Monti S, Tamayo P, Mesirov J, Golub T. Consensus clustering. *Mach Learn* 2003;52:91–118. doi:10.1023/A:1023949509487
15. Lancichinetti A, Fortunato S. Consensus clustering in complex networks. *Sci Rep* 2012;2:336. doi:10.1038/srep00336
16. Jeub LGS, Sporns O, Fortunato S. Multiresolution consensus clustering in networks. *Sci Rep* 2018;8:3259. doi:10.1038/s41598-018-21352-7
17. Langfelder P, Horvath S. WGCNA: an R package for weighted correlation network analysis. *BMC Bioinformatics* 2008;9:559. doi:10.1186/1471-2105-9-559
18. Langfelder P, Zhang B, Horvath S. Defining clusters from a hierarchical cluster tree: the Dynamic Tree Cut package for R. *Bioinformatics* 2008;24:719–720. doi:10.1093/bioinformatics/btm563
19. Ahn Y-Y, Bagrow JP, Lehmann S. Link communities reveal multiscale complexity in networks. *Nature* 2010;466:761–764. doi:10.1038/nature09182
20. Lancichinetti A, Radicchi F, Ramasco JJ, Fortunato S. Finding statistically significant communities in networks. *PLoS ONE* 2011;6:e18961. doi:10.1371/journal.pone.0018961
21. Lancichinetti A, Fortunato S, Kertész J. Detecting the overlapping and hierarchical community structure in complex networks. *New J Phys* 2009;11:033015. doi:10.1088/1367-2630/11/3/033015
22. Becker E, et al., Brun C. Multifunctional proteins revealed by overlapping clustering in protein interaction network. *Bioinformatics* 2012;28:84–90. doi:10.1093/bioinformatics/btr621
23. Zheng F, et al., Ideker T. HiDeF: identifying persistent structures in multiscale 'omics data. *Genome Biol* 2021;22:21. doi:10.1186/s13059-020-02228-4
24. Qin Y, et al., Ideker T. A multi-scale map of cell structure fusing protein images and interactions. *Nature* 2021;600:536–542. doi:10.1038/s41586-021-04115-9
25. Schaffer LV, et al., Ideker T. Multimodal cell maps as a foundation for structural and functional genomics. *Nature* 2025;642:222–231. doi:10.1038/s41586-025-08878-3
26. Zheng F, et al., Ideker T. Interpretation of cancer mutations using a multiscale map of protein systems. *Science* 2021;374:eabf3067. doi:10.1126/science.abf3067
27. Wang B, et al., Goldenberg A. Similarity network fusion for aggregating data types on a genomic scale. *Nat Methods* 2014;11:333–337. doi:10.1038/nmeth.2810
28. Meinshausen N, Bühlmann P. Stability selection. *J R Stat Soc B* 2010;72:417–473. doi:10.1111/j.1467-9868.2010.00740.x

**Complex maps, PPI resources and curated databases**
29. Drew K, et al., Marcotte EM. Integration of over 9,000 mass spectrometry experiments builds a global map of human protein complexes. *Mol Syst Biol* 2017;13:932. doi:10.15252/msb.20167490
30. Drew K, Wallingford JB, Marcotte EM. hu.MAP 2.0. *Mol Syst Biol* 2021;17:e10016. doi:10.15252/msb.202010016
31. Fischer SN, et al., Drew K. hu.MAP3.0: atlas of human protein complexes by integration of >25,000 proteomic experiments. *Mol Syst Biol* 2025;21:911–943. doi:10.1038/s44320-025-00121-5 [FT]
32. Huttlin EL, et al., Gygi SP. Dual proteome-scale networks reveal cell-specific remodeling of the human interactome. *Cell* 2021;184:3022–3040. doi:10.1016/j.cell.2021.04.011. PMID 33961781 [FT]
33. Huttlin EL, et al., Harper JW. Architecture of the human interactome defines protein communities and disease networks. *Nature* 2017;545:505–509. doi:10.1038/nature22366
34. Huttlin EL, et al., Gygi SP. The BioPlex network. *Cell* 2015;162:425–440. doi:10.1016/j.cell.2015.06.043
35. Cho NH, et al., Leonetti MD. OpenCell. *Science* 2022;375:eabi6983. doi:10.1126/science.abi6983. PMID 35271311
36. Hein MY, et al., Mann M. A human interactome in three quantitative dimensions. *Cell* 2015;163:712–723. doi:10.1016/j.cell.2015.09.053
37. Go CD, et al., Gingras A-C. A proximity-dependent biotinylation map of a human cell. *Nature* 2021;595:120–124. doi:10.1038/s41586-021-03592-2
38. Luck K, et al., Calderwood MA. A reference map of the human binary protein interactome. *Nature* 2020;580:402–408. doi:10.1038/s41586-020-2188-x
39. Sowa ME, Bennett EJ, Gygi SP, Harper JW. Defining the human deubiquitinating enzyme interaction landscape. *Cell* 2009;138:389–403. doi:10.1016/j.cell.2009.04.042
40. Mellacheruvu D, et al., Nesvizhskii AI. The CRAPome. *Nat Methods* 2013;10:730–736. doi:10.1038/nmeth.2557
41. Szklarczyk D, et al., von Mering C. The STRING database in 2023. *Nucleic Acids Res* 2023;51:D638–D646. doi:10.1093/nar/gkac1000
42. Szklarczyk D, et al., von Mering C. The STRING database in 2025. *Nucleic Acids Res* 2025;53:D730–D737. doi:10.1093/nar/gkae1113. PMID 39558183
43. Tsitsiridis G, et al., Ruepp A. CORUM: the comprehensive resource of mammalian protein complexes–2022. *Nucleic Acids Res* 2023;51:D539–D545. doi:10.1093/nar/gkac1015
44. Meldal BHM, et al., Orchard S. Complex Portal 2022. *Nucleic Acids Res* 2022;50:D578–D586. doi:10.1093/nar/gkab991

**Co-essentiality networks, modules and rewiring**
45. Pan J, et al., Kadoch C. Interrogation of mammalian protein complex structure, function, and membership using genome-scale fitness screens. *Cell Syst* 2018;6:555–568. doi:10.1016/j.cels.2018.04.011 [FT]
46. Boyle EA, Pritchard JK, Greenleaf WJ. High-resolution mapping of cancer cell networks using co-functional interactions. *Mol Syst Biol* 2018;14:e8594. doi:10.15252/msb.20188594 [FT]
47. Kim E, et al., Hart T. A network of human functional gene interactions from knockout fitness screens in cancer cells. *Life Sci Alliance* 2019;2:e201800278. doi:10.26508/lsa.201800278 [FT]
48. Wainberg M, et al., Bassik MC. A genome-wide atlas of co-essential modules assigns function to uncharacterized genes. *Nat Genet* 2021;53:638–649. doi:10.1038/s41588-021-00840-z [FT]
49. Shimada K, Bachman JA, Muhlich JL, Mitchison TJ. shinyDepMap. *eLife* 2021;10:e57116. doi:10.7554/eLife.57116 [FT]
50. Chen X, et al., Lu C. Context-defined cancer co-dependency mapping identifies a functional interplay between PRC2 and MLL-MEN1 complex in lymphoma. *Nat Commun* 2023;14:4259. doi:10.1038/s41467-023-39990-5. PMID 37460547 [FT]
51. Pan J, et al., Hahn WC. Sparse dictionary learning recovers pleiotropy from human cell fitness screens. *Cell Syst* 2022;13:286–303. doi:10.1016/j.cels.2021.12.005. PMID 35085500
52. Kim E, et al., Hart T. Dynamic rewiring of biological activity across genotype and lineage revealed by context-dependent functional interactions. *Genome Biol* 2022;23:140. doi:10.1186/s13059-022-02712-z. PMID 35768873 [FT]
53. Amici DR, et al., Mendillo ML. FIREWORKS: a bottom-up approach to integrative coessentiality network analysis. *Life Sci Alliance* 2021;4:e202000882. doi:10.26508/lsa.202000882. PMID 33328249 [ABS]
54. Billmann M, et al., Boone C. Global genetic interaction network of a human cell maps conserved principles and informs functional interpretation of gene co-essentiality profiles. *Cell* 2026;189:3736–3759. doi:10.1016/j.cell.2026.03.044. PMID 42049019 [ABS]
55. Camacho LQ, Fallahi-Sichani M. Systematic multivariate analysis of chromatin complex dependencies reveals Set1C/COMPASS as a melanoma-enriched epigenetic vulnerability. *PLoS Comput Biol* 2026;22:e1014018. doi:10.1371/journal.pcbi.1014018. PMID 42658887 [ABS]
56. Ideker T, Krogan NJ. Differential network biology. *Mol Syst Biol* 2012;8:565. doi:10.1038/msb.2011.99
57. Tesson BM, Breitling R, Jansen RC. DiffCoEx. *BMC Bioinformatics* 2010;11:497. doi:10.1186/1471-2105-11-497
58. Ha MJ, Baladandayuthapani V, Do K-A. DINGO: differential network analysis in genomics. *Bioinformatics* 2015;31:3413–3420. doi:10.1093/bioinformatics/btv406
59. Zhao SD, Cai TT, Li H. Direct estimation of differential networks. *Biometrika* 2014;101:253–268. doi:10.1093/biomet/asu009

**DepMap methods and statistics**
60. Tsherniak A, et al., Hahn WC. Defining a cancer dependency map. *Cell* 2017;170:564–576. doi:10.1016/j.cell.2017.06.010. PMID 28753430
61. Meyers RM, et al., Tsherniak A. Computational correction of copy number effect improves specificity of CRISPR–Cas9 essentiality screens. *Nat Genet* 2017;49:1779–1784. doi:10.1038/ng.3984
62. Dempster JM, et al., Iorio F. Agreement between two large pan-cancer CRISPR-Cas9 gene dependency data sets. *Nat Commun* 2019;10:5817. doi:10.1038/s41467-019-13805-y
63. Dempster JM, et al. Extracting biological insights from the Project Achilles genome-scale CRISPR screens in cancer cell lines. *bioRxiv* 2019. doi:10.1101/720243
64. Dempster JM, et al., McFarland JM. Chronos. *Genome Biol* 2021;22:343. doi:10.1186/s13059-021-02540-7. PMID 34930405
65. Behan FM, et al., Garnett MJ. Prioritization of cancer therapeutic targets using CRISPR–Cas9 screens. *Nature* 2019;568:511–516. doi:10.1038/s41586-019-1103-9. PMID 30971826 [ABS]
66. Pacini C, et al., Garnett MJ. A comprehensive clinically informed map of dependencies in cancer cells and framework for target prioritization. *Cancer Cell* 2024;42:301–316. doi:10.1016/j.ccell.2023.12.016. PMID 38215750 [ABS]
67. Pacini C, et al., Iorio F. Integrated cross-study datasets of genetic dependencies in cancer. *Nat Commun* 2021;12:1661. doi:10.1038/s41467-021-21898-7
68. Arafeh R, et al., Vazquez F. The present and future of the Cancer Dependency Map. *Nat Rev Cancer* 2025;25:59–73. doi:10.1038/s41568-024-00763-x
69. Lord CJ, Quinn N, Ryan CJ. Integrative analysis of large-scale loss-of-function screens identifies robust cancer-associated genetic interactions. *eLife* 2020;9:e58925. doi:10.7554/eLife.58925. PMID 32463358 [FT]
70. Krill-Burger JM, et al., Tsherniak A. Partial gene suppression improves identification of cancer vulnerabilities when CRISPR-Cas9 knockout is pan-lethal. *Genome Biol* 2023;24:192. doi:10.1186/s13059-023-03020-w
71. Vinceti A, et al., Iorio F. CoRe: a robustly benchmarked R package for identifying core-fitness genes. *BMC Genomics* 2021;22:828. doi:10.1186/s12864-021-08129-5
72. Hart T, et al., Moffat J. High-resolution CRISPR screens reveal fitness genes and genotype-specific cancer liabilities. *Cell* 2015;163:1515–1526. doi:10.1016/j.cell.2015.11.015. PMID 26627737
73. Hart T, et al., Moffat J. Evaluation and design of genome-wide CRISPR/SpCas9 knockout screens. *G3* 2017;7:2719–2727. doi:10.1534/g3.117.041277
74. Aguirre AJ, et al., Hahn WC. Genomic copy number dictates a gene-independent cell response to CRISPR/Cas9 targeting. *Cancer Discov* 2016;6:914–929. doi:10.1158/2159-8290.CD-16-0154
75. Munoz DM, et al., Schlabach MR. CRISPR screens provide a comprehensive assessment of cancer vulnerabilities but generate false-positive hits for highly amplified genomic regions. *Cancer Discov* 2016;6:900–913. doi:10.1158/2159-8290.CD-16-0178
76. Lazar NH, et al., Haque IS. High-resolution genome-wide mapping of chromosome-arm-scale truncations induced by CRISPR–Cas9 editing. *Nat Genet* 2024;56:1482–1493. doi:10.1038/s41588-024-01758-y. PMID 38811841
77. McFarland JM, et al., Tsherniak A. Improved estimation of cancer dependencies from large-scale RNAi screens. *Nat Commun* 2018;9:4610. doi:10.1038/s41467-018-06916-5
78. Ghandi M, et al., Sellers WR. Next-generation characterization of the Cancer Cell Line Encyclopedia. *Nature* 2019;569:503–508. doi:10.1038/s41586-019-1186-3
79. Corsello SM, et al., Golub TR. Discovering the anticancer potential of non-oncology drugs by systematic viability profiling. *Nat Cancer* 2020;1:235–248. doi:10.1038/s43018-019-0018-6
80. McDonald ER 3rd, et al., Sellers WR. Project DRIVE. *Cell* 2017;170:577–592. doi:10.1016/j.cell.2017.07.005
81. Dharia NV, et al., Stegmaier K. A first-generation pediatric cancer dependency map. *Nat Genet* 2021;53:529–538. doi:10.1038/s41588-021-00819-w
82. Shi X, et al., Dezso Z. Building a translational cancer dependency map for The Cancer Genome Atlas. *Nat Cancer* 2024;5:1176–1194. doi:10.1038/s43018-024-00789-y
83. Benjamini Y, Hochberg Y. Controlling the false discovery rate. *J R Stat Soc B* 1995;57:289–300. doi:10.1111/j.2517-6161.1995.tb02031.x
84. Van den Berge K, Soneson C, Robinson MD, Clement L. stageR. *Genome Biol* 2017;18:151. doi:10.1186/s13059-017-1277-0
85. Winkler AM, Ridgway GR, Webster MA, Smith SM, Nichols TE. Permutation inference for the general linear model. *NeuroImage* 2014;92:381–397. doi:10.1016/j.neuroimage.2014.01.060
86. Ritchie ME, et al., Smyth GK. limma powers differential expression analyses. *Nucleic Acids Res* 2015;43:e47. doi:10.1093/nar/gkv007
87. Liu Y, Xie J. Cauchy combination test. *J Am Stat Assoc* 2020;115:393–402. doi:10.1080/01621459.2018.1554485
88. Liu Y, et al., Lin X. ACAT: a fast and powerful p value combination method. *Am J Hum Genet* 2019;104:410–421. doi:10.1016/j.ajhg.2019.01.002

**Module and complex scoring**
89. Langfelder P, Horvath S. Eigengene networks for studying the relationships between co-expression modules. *BMC Syst Biol* 2007;1:54. doi:10.1186/1752-0509-1-54
90. Tomfohr J, Lu J, Kepler TB. Pathway level analysis of gene expression using singular value decomposition. *BMC Bioinformatics* 2005;6:225. doi:10.1186/1471-2105-6-225
91. Foroutan M, et al., Davis MJ. Single sample scoring of molecular phenotypes. *BMC Bioinformatics* 2018;19:404. doi:10.1186/s12859-018-2435-4
92. Hänzelmann S, Castelo R, Guinney J. GSVA. *BMC Bioinformatics* 2013;14:7. doi:10.1186/1471-2105-14-7
93. Barbie DA, et al., Hahn WC. Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1. *Nature* 2009;462:108–112. doi:10.1038/nature08460

**Context biology (positive controls)**
94. Helming KC, et al., Roberts CWM. ARID1B is a specific vulnerability in ARID1A-mutant cancers. *Nat Med* 2014;20:251–254. doi:10.1038/nm.3480
95. Hoffman GR, et al., Jagani Z. Functional epigenetics approach identifies BRM/SMARCA2 as a critical synthetic lethal target in BRG1-deficient cancers. *PNAS* 2014;111:3128–3133. doi:10.1073/pnas.1316793111
96. Wilson BG, et al., Roberts CWM. Residual complexes containing SMARCA2 (BRM) underlie the oncogenic drive of SMARCA4 (BRG1) mutation. *Mol Cell Biol* 2014;34:1136–1144. doi:10.1128/MCB.01372-13
97. Michel BC, et al., Kadoch C. A non-canonical SWI/SNF complex is a synthetic lethal target in cancers driven by BAF complex perturbation. *Nat Cell Biol* 2018;20:1410–1420. doi:10.1038/s41556-018-0221-1
98. Brien GL, et al., Armstrong SA. Targeted degradation of BRD9 reverses oncogenic gene expression in synovial sarcoma. *eLife* 2018;7:e41305. doi:10.7554/eLife.41305
99. Wang X, et al., Roberts CWM. BRD9 defines a SWI/SNF sub-complex and constitutes a specific vulnerability in malignant rhabdoid tumors. *Nat Commun* 2019;10:1881. doi:10.1038/s41467-019-09891-7
100. Radko-Juettner S, et al., Roberts CWM. Targeting DCAF5 suppresses SMARCB1-mutant cancer by stabilizing SWI/SNF. *Nature* 2024;628:442–449. doi:10.1038/s41586-024-07250-1. PMID 38538798 [ABS]
101. He T, et al., Chinnaiyan AM. Targeting the mSWI/SNF complex in POU2F-POU2AF transcription factor-driven malignancies. *Cancer Cell* 2024;42:1336–1351. doi:10.1016/j.ccell.2024.06.006 [ABS]
102. Duplaquet L, et al., Oser MG. Mammalian SWI/SNF complex activity regulates POU2F3 and constitutes a targetable dependency in small cell lung cancer. *Cancer Cell* 2024;42:1352–1369. doi:10.1016/j.ccell.2024.06.012
103. Kryukov GV, et al., Garraway LA. MTAP deletion confers enhanced dependency on the PRMT5 arginine methyltransferase in cancer cells. *Science* 2016;351:1214–1218. doi:10.1126/science.aad5214
104. Mavrakis KJ, et al., Sellers WR. Disordered methionine metabolism in MTAP/CDKN2A-deleted cancers leads to dependence on PRMT5. *Science* 2016;351:1208–1213. doi:10.1126/science.aad5944
105. Marjon K, et al., Marks KM. MTAP deletions in cancer create vulnerability to targeting of the MAT2A/PRMT5/RIOK1 axis. *Cell Rep* 2016;15:574–587. doi:10.1016/j.celrep.2016.03.043
106. Chan EM, et al., Bass AJ. WRN helicase is a synthetic lethal target in microsatellite unstable cancers. *Nature* 2019;568:551–556. doi:10.1038/s41586-019-1102-x
107. Lieb S, et al., Wöhrle S. Werner syndrome helicase is a selective vulnerability of microsatellite instability-high tumor cells. *eLife* 2019;8:e43333. doi:10.7554/eLife.43333
108. Kategaya L, Perumal SK, Hager JH, Belmont LD. Werner syndrome helicase is required for the survival of cancer cells with microsatellite instability. *iScience* 2019;13:488–497. doi:10.1016/j.isci.2019.02.006
109. van Wietmarschen N, et al., Nussenzweig A. Repeat expansions confer WRN dependence in microsatellite-unstable cancers. *Nature* 2020;586:292–298. doi:10.1038/s41586-020-2769-8
110. Picco G, et al., Garnett MJ. Novel WRN helicase inhibitors selectively target microsatellite-unstable cancer cells. *Cancer Discov* 2024;14:1457–1475. doi:10.1158/2159-8290.CD-24-0052
111. Prindle V, et al., Plotnik JP. Synthetic lethality of mRNA quality control complexes in cancer. *Nature* 2025;638:1095–1103. doi:10.1038/s41586-024-08398-6. PMID 39910291 [ABS]
112. Nijhawan D, et al., Hahn WC. Cancer vulnerabilities unveiled by genomic loss. *Cell* 2012;150:842–854. doi:10.1016/j.cell.2012.07.023. PMID 22901813 [ABS]
113. Liu Y, et al., Lu X. TP53 loss creates therapeutic vulnerability in colorectal cancer. *Nature* 2015;520:697–701. doi:10.1038/nature14418. PMID 25901683
114. Neggers JE, et al., Aguirre AJ. Synthetic lethal interaction between the ESCRT paralog enzymes VPS4A and VPS4B in cancers harboring loss of chromosome 18q or 16q. *Cell Rep* 2020;33:108493. doi:10.1016/j.celrep.2020.108493
115. Ryan CJ, Kennedy S, Bajrami I, Matallanas D, Lord CJ. A compendium of co-regulated protein complexes in breast cancer reveals collateral loss events. *Cell Syst* 2017;5:399–409. doi:10.1016/j.cels.2017.09.011. PMID 29032073 [ABS]
116. Schukken KM, Sheltzer JM. Extensive protein dosage compensation in aneuploid human cancers. *Genome Res* 2022;32:1254–1270. doi:10.1101/gr.276378.121
117. Rendo V, et al., Foijer F. A compendium of Amplification-Related Gain Of Sensitivity genes in human cancer. *Nat Commun* 2025;16:1077. doi:10.1038/s41467-025-56301-2 [ABS]
118. Cohen-Sharir Y, et al., Ben-David U. Aneuploidy renders cancer cells vulnerable to mitotic checkpoint inhibition. *Nature* 2021;590:486–491. doi:10.1038/s41586-020-03114-6
119. De Kegel B, Ryan CJ. Paralog buffering contributes to the variable essentiality of genes in cancer cell lines. *PLoS Genet* 2019;15:e1008466. doi:10.1371/journal.pgen.1008466 [ABS]
120. De Kegel B, Quinn N, Thompson NA, Adams DJ, Ryan CJ. Comprehensive prediction of robust synthetic lethality between paralog pairs in cancer cell lines. *Cell Syst* 2021;12:1144–1159. doi:10.1016/j.cels.2021.08.006 [ABS]
121. Dede M, McLaughlin M, Kim E, Hart T. Multiplex enCas12a screens detect functional buffering among paralogs. *Genome Biol* 2020;21:262. doi:10.1186/s13059-020-02173-2
122. Thompson NA, et al., Adams DJ. Combinatorial CRISPR screen identifies fitness effects of gene paralogues. *Nat Commun* 2021;12:1302. doi:10.1038/s41467-021-21478-9
123. Parrish PCR, et al., Berger AH. Discovery of synthetic lethal and tumor suppressor paralog pairs in the human genome. *Cell Rep* 2021;36:109597. doi:10.1016/j.celrep.2021.109597
124. Ito T, et al., Sellers WR. Paralog knockout profiling identifies DUSP4 and DUSP6 as a digenic dependence in MAPK pathway-driven cancers. *Nat Genet* 2021;53:1664–1672. doi:10.1038/s41588-021-00967-z
125. Lee AR, Tangiyan A, Singh I, Choi PS. Incomplete paralog compensation generates selective dependency on TRA2A in cancer. *PLoS Genet* 2025;21:e1011685. doi:10.1371/journal.pgen.1011685
126. Dandage R, Landry CR. Paralog dependency indirectly affects the robustness of human cells. *Mol Syst Biol* 2019;15:e8871. doi:10.15252/msb.20198871
127. Diss G, et al., Landry CR. Gene duplication can impart fragility, not robustness, in the yeast protein interaction network. *Science* 2017;355:630–634. doi:10.1126/science.aai7685
128. van der Lelij P, et al., Petronczki M. Synthetic lethality between the cohesin subunits STAG1 and STAG2 in diverse cancer contexts. *eLife* 2017;6:e26980. doi:10.7554/eLife.26980
129. Ogiwara H, et al., Kohno T. Targeting p300 addiction in CBP-deficient cancers causes synthetic lethality by apoptotic cell death due to abrogation of MYC expression. *Cancer Discov* 2016;6:430–445. doi:10.1158/2159-8290.CD-15-0754
130. Romero R, et al., Papagiannakopoulos T. Keap1 loss promotes Kras-driven lung cancer and results in dependence on glutaminolysis. *Nat Med* 2017;23:1362–1368. doi:10.1038/nm.4407
131. Bar-Peled L, et al., Cravatt BF. Chemical proteomics identifies druggable vulnerabilities in a genetically defined cancer. *Cell* 2017;171:696–709. doi:10.1016/j.cell.2017.08.051. PMID 28965760
132. Chen W, et al., Brugarolas J. Targeting renal cell carcinoma with a HIF-2 antagonist. *Nature* 2016;539:112–117. doi:10.1038/nature19796
133. Nicholson HE, et al., Kaelin WG. HIF-independent synthetic lethality between CDK4/6 inhibition and VHL loss across species. *Sci Signal* 2019;12:eaay0482. doi:10.1126/scisignal.aay0482
134. Durbin AD, et al., Stegmaier K. Selective gene dependencies in MYCN-amplified neuroblastoma include the core transcriptional regulatory circuitry. *Nat Genet* 2018;50:1240–1246. doi:10.1038/s41588-018-0191-z
135. Malone CF, et al., Stegmaier K. The KAT module of the SAGA complex maintains the oncogenic gene expression program in MYCN-amplified neuroblastoma. *Sci Adv* 2024;10:eadm9449. doi:10.1126/sciadv.adm9449. PMID 38820154 [ABS]
136. Pelish HE, et al., Shair MD. Mediator kinase inhibition further activates super-enhancer-associated genes in AML. *Nature* 2015;526:273–276. doi:10.1038/nature14904
137. Zhang S, et al., Stegmaier K. CDK8 inhibition releases the muscle differentiation block in fusion-driven alveolar rhabdomyosarcoma. *Cancer Discov* 2026 (online first). doi:10.1158/2159-8290.CD-25-1171
138. Krivtsov AV, et al., Armstrong SA. A menin-MLL inhibitor induces specific chromatin changes and eradicates disease in models of MLL-rearranged leukemia. *Cancer Cell* 2019;36:660–673. doi:10.1016/j.ccell.2019.11.001
139. Wan L, et al., Shi X. ENL links histone acetylation to oncogenic gene expression in acute myeloid leukaemia. *Nature* 2017;543:265–269. doi:10.1038/nature21687
140. Erb MA, et al., Bradner JE. Transcription control by the ENL YEATS domain in acute leukaemia. *Nature* 2017;543:270–274. doi:10.1038/nature21688
141. Yan F, et al., Blanco MA. KAT6A and ENL form an epigenetic transcriptional control module to drive critical leukemogenic gene-expression programs. *Cancer Discov* 2022;12:792–811. doi:10.1158/2159-8290.CD-20-1459
142. Wang T, et al., Sabatini DM. Gene essentiality profiling reveals gene networks and synthetic lethal interactions with oncogenic Ras. *Cell* 2017;168:890–903. doi:10.1016/j.cell.2017.01.013
143. Molina JR, et al., Marszalek JR. An inhibitor of oxidative phosphorylation exploits cancer vulnerability. *Nat Med* 2018;24:1036–1046. doi:10.1038/s41591-018-0052-4
144. Škrtić M, et al., Schimmer AD. Inhibition of mitochondrial translation as a therapeutic strategy for human acute myeloid leukemia. *Cancer Cell* 2011;20:674–688. doi:10.1016/j.ccr.2011.10.015
145. Arroyo JD, et al., Mootha VK. A genome-wide CRISPR death screen identifies genes essential for oxidative phosphorylation. *Cell Metab* 2016;24:875–885. doi:10.1016/j.cmet.2016.08.017
146. Rossiter NJ, et al., Cantor JR. CRISPR screens in physiologic medium reveal conditionally essential genes in human cells. *Cell Metab* 2021;33:1248–1263. doi:10.1016/j.cmet.2021.02.005
147. Bondeson DP, et al., Golub TR. Phosphate dysregulation via the XPR1–KIDINS220 protein complex is a therapeutic vulnerability in ovarian cancer. *Nat Cancer* 2022;3:681–695. doi:10.1038/s43018-022-00360-7

**Proteomics (context and extra data)**
148. Nusinow DP, et al., Gygi SP. Quantitative proteomics of the Cancer Cell Line Encyclopedia. *Cell* 2020;180:387–402. doi:10.1016/j.cell.2019.12.023
149. Gonçalves E, et al., Reddel RR. Pan-cancer proteomic map of 949 human cell lines. *Cancer Cell* 2022;40:835–849. doi:10.1016/j.ccell.2022.06.010
150. Sousa A, et al., Beltrao P. Multi-omics characterization of interaction-mediated control of human protein abundance levels. *Mol Cell Proteomics* 2019;18:S114–S125. doi:10.1074/mcp.RA118.001280

**Evaluation statistics (test set)**
151. Saito T, Rehmsmeier M. The precision-recall plot is more informative than the ROC plot when evaluating binary classifiers on imbalanced datasets. *PLoS ONE* 2015;10:e0118432. doi:10.1371/journal.pone.0118432
152. Davis J, Goadrich M. The relationship between precision-recall and ROC curves. *ICML* 2006:233–240. doi:10.1145/1143844.1143874
153. Boyd K, Eng KH, Page CD. Area under the precision-recall curve: point estimates and confidence intervals. *ECML PKDD*, LNCS 8190, 2013:451–466. doi:10.1007/978-3-642-40994-3_29
154. Field CA, Welsh AH. Bootstrapping clustered data. *J R Stat Soc B* 2007;69:369–390. doi:10.1111/j.1467-9868.2007.00593.x
155. Owen AB. The pigeonhole bootstrap. *Ann Appl Stat* 2007;1:386–411. doi:10.1214/07-AOAS122
156. Davezies L, D'Haultfœuille X, Guyonvarch Y. Empirical process results for exchangeable arrays. *Ann Stat* 2021;49:845–862. doi:10.1214/20-AOS1981
157. Cameron AC, Gelbach JB, Miller DL. Robust inference with multiway clustering. *J Bus Econ Stat* 2011;29:238–249. doi:10.1198/jbes.2010.07136
158. Zadrozny B, Elkan C. Transforming classifier scores into accurate multiclass probability estimates. *KDD* 2002:694–699. doi:10.1145/775047.775151
159. Niculescu-Mizil A, Caruana R. Predicting good probabilities with supervised learning. *ICML* 2005:625–632. doi:10.1145/1102351.1102430
160. Van Calster B, et al., Steyerberg EW. Calibration: the Achilles heel of predictive analytics. *BMC Med* 2019;17:230. doi:10.1186/s12916-019-1466-7
161. Park Y, Marcotte EM. Flaws in evaluation schemes for pair-input computational predictions. *Nat Methods* 2012;9:1134–1136. doi:10.1038/nmeth.2259

**Cited via Part 1 (verified there)**
- Rahman 2021 (doi:10.15252/msb.202010013), Gheorghe & Hart 2022 (doi:10.1186/s12859-022-05078-y), Hassan 2023 (doi:10.15252/msb.202311657), Rohde 2025 (doi:10.1038/s44320-025-00103-7), Vinceti 2024 (doi:10.1186/s13059-024-03336-1), Schmid & Walter 2025 (doi:10.1016/j.molcel.2025.01.034).
- Their records were re-checked in Crossref today.

**Unverified or weak items**
- Uiterkamp et al. 2026, SKA complex in CIN cancers. bioRxiv; title and DOI seen only in a Europe PMC listing.
- Kim JH et al. 2026, *Molecules*, proteasome PSMB5/6 (doi:10.3390/molecules31111954; listing only).
- Neither is load-bearing.
- Exact BioPlex HCT116 interaction count and OpenCell target count: [unverified].
- Behan 2019 and Pacini 2024 statistical details: [unverified beyond abstract].
