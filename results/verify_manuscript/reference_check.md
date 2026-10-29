# Reference check: manuscript/main.md, supplementary.md, cover_letter.md

Independent reference checker, 2026-10-29. Scope: every bracketed citation in `manuscript/main.md` (64 reference-list entries; 41 citation groups), checked against `LITERATURE.md` §9 and §10. `manuscript/supplementary.md` and `manuscript/cover_letter.md` contain **no** bracketed citations; they were checked for statements that need a citation and for STRING / hu.MAP / CORUM use.

Method:
- Metadata (authors, year, title, journal, volume, issue, pages): Crossref REST API record for every DOI, cross-checked against the Europe PMC record (PMID, volume, pages) where indexed. The three no-DOI items (Holm 1979, Kish 1965, Bell & McCaffrey 2002) were checked by web search against JSTOR / Statistics Canada / secondary records.
- Retraction: Crossref `updated-by` / `filter=updates:<DOI>` notices and the Europe PMC publication type for every DOI. **No retraction, withdrawal or expression of concern was found for any cited item.** Two items have non-retraction updates: Wang 2019 has a published correction (10.1038/s41467-019-12524-8); De Kegel & Ryan 2019 has only a self "new_version" record.
- Claim support: Europe PMC abstracts; Europe PMC full text for Prindle 2025, Fischer 2025, Chen 2023 (open access); PMC author manuscript for Wainberg 2021. For classical statistics and methods papers the claimed use is the paper's named method, checked against title and abstract.

Status key: OK, FIX-METADATA (the reference-list entry is wrong or truncated), CLAIM-MISMATCH (the paper does not support the sentence as written or placed), NOT-FOUND, RETRACTED.

## Summary counts

| Status | n |
|---|---|
| OK | 53 |
| FIX-METADATA | 9 |
| CLAIM-MISMATCH | 2 (both minor) |
| NOT-FOUND | 0 |
| RETRACTED | 0 |

All 13 §10 methods references (M1–M13), including the 6 marked [VERIFY], exist and their metadata as given are correct. For M10, an Ensembl paper that matches release 116 more closely now exists; see below.

**STRING / hu.MAP / CORUM as evidence: none found.** CORUM, Complex Portal and STRING are cited only in the Introduction as background. In Methods they are stated as not used. hu.MAP 3.0 appears only as the pre-registered Part 1 held-out test set and its data source. Neither use counts as Part 2 evidence. No change is needed.

## Table

| Ref key | Status | Fix | Source URL |
|---|---|---|---|
| Aguirre 2016 [74] | OK | – | https://doi.org/10.1158/2159-8290.CD-16-0154 |
| Arafeh 2025 [68] | OK | – | https://doi.org/10.1038/s41568-024-00763-x |
| Bell & McCaffrey 2002 [M7] | OK | Remove [VERIFY]. Surv Methodol 2002;28(2):169–181 confirmed; no DOI. | https://www150.statcan.gc.ca/n1/en/catalogue/12-001-X20020029058 |
| Benjamini & Bogomolov 2014 [M1] | OK | – | https://doi.org/10.1111/rssb.12028 |
| Benjamini & Hochberg 1995 [83] | FIX-METADATA | The title is truncated. Full title: "Controlling the false discovery rate: a practical and powerful approach to multiple testing". | https://doi.org/10.1111/j.2517-6161.1995.tb02031.x |
| Blondel 2008 [3] | OK | – | https://doi.org/10.1088/1742-5468/2008/10/P10008 |
| Boyd 2013 [153] | FIX-METADATA | Add the series volume: "ECML PKDD 2013, LNCS 8190:451–466". | https://doi.org/10.1007/978-3-642-40994-3_29 |
| Boyle 2018 [46] | OK | – | https://doi.org/10.15252/msb.20188594 |
| Chan 2019 [106] | OK | – | https://doi.org/10.1038/s41586-019-1102-x |
| Chen 2023 [50] | OK | – | https://doi.org/10.1038/s41467-023-39990-5 |
| Davis & Goadrich 2006 [152] | OK | – | https://doi.org/10.1145/1143844.1143874 |
| De Kegel & Ryan 2019 [119] | OK | – | https://doi.org/10.1371/journal.pgen.1008466 |
| Dempster 2021 [64] | FIX-METADATA | The title is truncated to "Chronos". Full title: "Chronos: a cell population dynamics model of CRISPR experiments that improves inference of gene fitness effects". | https://doi.org/10.1186/s13059-021-02540-7 |
| Dyer 2025 [M10] | OK (update recommended) | The metadata are correct: NAR 53(D1):D948–D957, 91 authors, first author Dyer SC. Ensembl 116 is described in Ensembl 2026; see edit E10. | https://doi.org/10.1093/nar/gkae1071 |
| Field & Welsh 2007 [154] | OK | – | https://doi.org/10.1111/j.1467-9868.2007.00593.x |
| Fischer 2025 [31] | OK | – (21(7):911–943 confirmed) | https://doi.org/10.1038/s44320-025-00121-5 |
| Fortunato & Barthélemy 2007 [6] | OK | – | https://doi.org/10.1073/pnas.0605965104 |
| Ghandi 2019 [78] | OK | – | https://doi.org/10.1038/s41586-019-1186-3 |
| Helming 2014 [94] | OK | – | https://doi.org/10.1038/nm.3480 |
| Holm 1979 [M3] | OK | Remove [VERIFY]. Scand J Stat 6(2):65–70 confirmed; JSTOR 4615733; no DOI. | https://www.jstor.org/stable/4615733 |
| Hubert & Arabie 1985 [M13] | OK | – | https://doi.org/10.1007/BF01908075 |
| Huttlin 2021 [32] | OK | – | https://doi.org/10.1016/j.cell.2021.04.011 |
| Jeub 2018 [16] | OK | – | https://doi.org/10.1038/s41598-018-21352-7 |
| Kim 2019 [47] | OK | – | https://doi.org/10.26508/lsa.201800278 |
| Kim 2022 [52] | OK | – | https://doi.org/10.1186/s13059-022-02712-z |
| Kish 1965 [M12] | OK | Remove [VERIFY]. Kish L. Survey Sampling. New York: Wiley; 1965. A book, with no DOI. The (Σw)²/Σw² effective size is conventionally attributed to it. | https://search.r-project.org/CRAN/refmans/svyweight/html/eff_n.html (secondary record citing Kish 1965) |
| Knijnenburg 2009 [M8] | OK | – | https://doi.org/10.1093/bioinformatics/btp211 |
| Kryukov 2016 [103] | OK | – | https://doi.org/10.1126/science.aad5214 |
| Lancichinetti & Fortunato 2012 [15] | OK | – | https://doi.org/10.1038/srep00336 |
| Lazar 2024 [76] | OK | – | https://doi.org/10.1038/s41588-024-01758-y |
| Lieb 2019 [107] | OK | – | https://doi.org/10.7554/eLife.43333 |
| Long & Ervin 2000 [M5] | OK | Remove [VERIFY]. Am Stat 54(3):217–224 and the DOI are confirmed. | https://doi.org/10.1080/00031305.2000.10474549 |
| MacKinnon & White 1985 [M4] | OK | – | https://doi.org/10.1016/0304-4076(85)90158-7 |
| Mavrakis 2016 [104] | OK | – | https://doi.org/10.1126/science.aad5944 |
| McFarland 2018 [77] | FIX-METADATA | The title is truncated. Full title: "Improved estimation of cancer dependencies from large-scale RNAi screens using model-based normalization and data integration". | https://doi.org/10.1038/s41467-018-06916-5 |
| Meldal 2022 [44] | FIX-METADATA | The title is truncated. Full title: "Complex Portal 2022: new curation frontiers". | https://doi.org/10.1093/nar/gkab991 |
| Meyers 2017 [61] | FIX-METADATA | The title is truncated. Full title ends "…essentiality screens in cancer cells". | https://doi.org/10.1038/ng.3984 |
| Michel 2018 [97] | OK | – | https://doi.org/10.1038/s41556-018-0221-1 |
| Monti 2003 [14] | FIX-METADATA | The title is truncated. Full title: "Consensus clustering: a resampling-based method for class discovery and visualization of gene expression microarray data". | https://doi.org/10.1023/A:1023949509487 |
| Munoz 2016 [75] | OK | – | https://doi.org/10.1158/2159-8290.CD-16-0178 |
| Niculescu-Mizil & Caruana 2005 [159] | OK | – | https://doi.org/10.1145/1102351.1102430 |
| Nusinow 2020 [148] | OK | – | https://doi.org/10.1016/j.cell.2019.12.023 |
| Ogiwara 2016 [129] | OK | – | https://doi.org/10.1158/2159-8290.CD-15-0754 |
| Owen 2007 [155] | CLAIM-MISMATCH (placement) | The citation is attached to "Monte-Carlo seeds 1–4 show seed dependence". Owen 2007 is about resampling rows and columns of crossed (two-way) data, not about seed dependence. Move it to the protein (node) bootstrap clause (edit E3). | https://doi.org/10.1214/07-AOAS122 |
| Pan 2018 [45] | OK | – | https://doi.org/10.1016/j.cels.2018.04.011 |
| Park & Marcotte 2012 [161] | OK | – | https://doi.org/10.1038/nmeth.2259 |
| Prindle 2025 [111] | OK | – | https://doi.org/10.1038/s41586-024-08398-6 |
| Saito & Rehmsmeier 2015 [151] | OK | – | https://doi.org/10.1371/journal.pone.0118432 |
| Satterthwaite 1946 [M6] | OK | – (2(6):110–114) | https://doi.org/10.2307/3002019 |
| Seal 2023 [M11] | OK | Remove [VERIFY]. NAR 51(D1):D1003–D1009 and the DOI are confirmed. | https://doi.org/10.1093/nar/gkac888 |
| Sechidis 2011 [M9] | OK | Remove [VERIFY]. Pages 145–158 and the DOI are confirmed. The DOI prefix 978-3-642-23808-6 is ECML PKDD 2011 Part III (LNCS 6913). | https://doi.org/10.1007/978-3-642-23808-6_10 |
| Simes 1986 [M2] | OK | – | https://doi.org/10.1093/biomet/73.3.751 |
| Szklarczyk 2023 [41] | FIX-METADATA | The title is truncated. Full title: "The STRING database in 2023: protein–protein association networks and functional enrichment analyses for any sequenced genome of interest". | https://doi.org/10.1093/nar/gkac1000 |
| Traag 2011 [5] | OK | – | https://doi.org/10.1103/PhysRevE.84.016114 |
| Traag 2019 [4] | OK | – | https://doi.org/10.1038/s41598-019-41695-z |
| Tsherniak 2017 [60] | CLAIM-MISMATCH (minor, Methods table only) | Tsherniak 2017 is an RNAi (Achilles/DEMETER) paper. In the Introduction it is fine as the origin of DepMap. In the Methods data table it is the only primary citation for the CRISPRGeneEffect row; cite the CRISPR/Chronos papers there (edit E4). | https://doi.org/10.1016/j.cell.2017.06.010 |
| Tsitsiridis 2023 [43] | OK | – (cited as a resource only; see edit E7) | https://doi.org/10.1093/nar/gkac1015 |
| Van Calster 2019 [160] | OK | – | https://doi.org/10.1186/s12916-019-1466-7 |
| Van den Berge 2017 [84] | FIX-METADATA | The title is truncated to "stageR". Full title: "stageR: a general stage-wise method for controlling the gene-level false discovery rate in differential expression and differential transcript usage". | https://doi.org/10.1186/s13059-017-1277-0 |
| Wainberg 2021 [48] | OK | – | https://doi.org/10.1038/s41588-021-00840-z ; https://pmc.ncbi.nlm.nih.gov/articles/PMC8763319/ |
| Wang 2019 [99] | OK (note) | Not retracted. An author correction exists (10.1038/s41467-019-12524-8); it does not affect the cited claim. | https://doi.org/10.1038/s41467-019-09891-7 |
| Winkler 2014 [85] | OK | – | https://doi.org/10.1016/j.neuroimage.2014.01.060 |
| Zadrozny & Elkan 2002 [158] | OK | – | https://doi.org/10.1145/775047.775151 |
| Zheng 2021 [23] | OK | – | https://doi.org/10.1186/s13059-020-02228-4 |

## Claim check (manuscript sentence → what the paper says)

Introduction (main.md l.32–36)
- "Genome-scale CRISPR loss-of-function screens across large cancer cell-line panels, such as the Cancer Dependency Map (DepMap), measure the fitness effect of knocking out each gene…" [Tsherniak 2017; Meyers 2017; Dempster 2021; Arafeh 2025]. Supported as a group. Tsherniak: an RNAi-based dependency map of 501 lines. Meyers: CERES copy-number correction of CRISPR screens. Dempster: the Chronos model of CRISPR fitness effects. Arafeh: a review of DepMap.
- "…co-essentiality has been used to recover complexes, pathways and the functions of uncharacterised genes" [Pan 2018; Boyle 2018; Kim 2019; Wainberg 2021]. Supported. Wainberg: the atlas "recapitulates diverse pathways and protein complexes" and assigns 108 uncharacterised genes. The others build co-essential or co-functional networks recovering complexes.
- "Generalised least squares (GLS) on whitened profiles, which removes the dominant line-to-line covariance, is a widely used standard…" [Wainberg 2021]. Supported. The paper says GLS accounts for the non-independence of cell lines.
- "pair-input predictions are vulnerable to optimistic evaluation when test pairs share proteins with training pairs" [Park & Marcotte 2012]. Supported: shared proteins between training and test pairs inflate the performance estimates.
- "downstream complex-level analyses … usually start from curated complex resources (CORUM, Complex Portal) or functional networks (STRING)" [Tsitsiridis 2023; Meldal 2022; Szklarczyk 2023]. The three are correct resource papers, cited only to name the resources. They do not themselves show that downstream analyses "usually" start from them (edit E7). This is background, not evidence, so it is allowed.
- "evaluated once on the hu.MAP 3.0 test pairs" [Fischer 2025]. Supported. hu.MAP 3.0 evaluates on a "leave out test set" of gold-standard co-complex interactions.
- "physical (BioPlex 3.0 [Huttlin 2021])" and "proteomic (CCLE proteomics [Nusinow 2020])". Supported: 293T and HCT116 interactomes; quantitative proteomics of CCLE lines.

Results
- l.58 "Kish effective number of clusters 21.3" [Kish 1965]. Supported: the (Σw)²/Σw² effective size from Kish's design-effect framework.
- l.94 "the signature of copy-number and CRISPR-cutting artefacts" [Aguirre 2016; Munoz 2016; Lazar 2024]. Supported. Aguirre: copy number causes a gene-independent response to Cas9 cutting. Munoz: false-positive hits in amplified regions. Lazar: "proximity bias" from chromosome-arm truncations.
- l.106 "adjusted Rand index (ARI)" [Hubert & Arabie 1985]. Supported. This is the source of the corrected-for-chance ARI.
- l.120 "WRN complex in MSI-H lines … PRMT5 complex with 9p21 loss … BRD9-containing complex … with SMARCB1 loss or SS18::SSX fusion … PELO–HBS1L complex with FOCAD deletion or MSI-H …, each matching published dependencies" [Chan 2019; Lieb 2019; Kryukov 2016; Mavrakis 2016; Michel 2018; Wang 2019; Prindle 2025]. Each part is supported:
  - Chan and Lieb: WRN is selectively required in MSI-H cells.
  - Kryukov and Mavrakis: MTAP deletion confers PRMT5 dependence.
  - Michel: ncBAF/BRD9 is synthetic lethal in synovial sarcoma and malignant rhabdoid tumours.
  - Wang: BRD9 is a vulnerability in SMARCB1-mutant rhabdoid tumours.
  - Prindle: in 9p21.3-deleted and MSI-H tumours, a destabilised SKI complex creates dependence on PELO–HBS1L. The full text names FOCAD as the 9p21.3 gene.
- l.174 "ARID1A → ARID1B … reproduces a published dependency" [Helming 2014]. Supported (the paper's title).
- l.175 "the reciprocal of the published CREBBP-deficient → EP300 dependency" [Ogiwara 2016]. Supported: CBP-deficient cancers depend on p300.
- l.175 "paralogues that act within complexes are known to be enriched for synthetic lethality" [De Kegel & Ryan 2019]. Supported. Paralog pairs functioning in protein complexes are significantly more likely to be synthetic lethal.
- l.176 "the published dependency of SKI-complex-destabilised cancers on PELO–HBS1L" [Prindle 2025]. Supported (as above).
- l.232 "FURIN and CPD are reported lineage-specific co-dependency partners of IGF1R" [Kim 2022]. Supported. The paper shows lineage-specific function through FURIN and CPD maturation of IGF1R and MET.

Discussion
- l.281 "Earlier module-level lineage analyses did not calibrate against such matched sets" [Wainberg 2021; Chen 2023]. The premise is supported. Wainberg found "444 modules … differentially essential … in cancers derived from 16 distinct tissue types". Chen found "137 cancer type-specific modules". Neither abstract or method text describes a matched random-gene-set calibration; for Wainberg this was confirmed from the author manuscript summary.
- l.287 Helming, Ogiwara, Prindle, De Kegel & Ryan: as above. "paralog buffering acts within complexes" matches De Kegel & Ryan's complex-membership finding.

Methods
- l.320–326 data table:
  - Tsherniak/Arafeh for the DepMap portal: see the CLAIM-MISMATCH row above.
  - Ghandi 2019 for DepMap Omics: supported (CCLE multi-omic characterisation).
  - McFarland 2018 for DEMETER2: supported.
  - Huttlin 2021 for BioPlex: supported.
  - Dyer 2025 for Ensembl: correct but superseded (edit E10).
  - Seal 2023 for HGNC: supported.
- l.335 "signed GLS t-statistics" [Wainberg 2021]. Supported.
- l.343 average precision [Saito & Rehmsmeier 2015; Davis & Goadrich 2006; Boyd 2013]. Supported: PR curves are informative under imbalance; PR–ROC relationship; AUCPR point estimates and CIs.
- l.343 "components > 40 proteins split by Louvain" [Blondel 2008]. Supported.
- l.345 isotonic and Platt calibration [Zadrozny & Elkan 2002; Niculescu-Mizil & Caruana 2005]. Supported: isotonic calibration of scores; comparison of Platt and isotonic calibration. Calibration measures [Van Calster 2019] are supported: calibration-in-the-large, slope and calibration curves. Note: the manuscript's calibration-in-the-large is a ratio Σp/Σy, while Van Calster frames it as mean predicted vs observed. It is worth stating "expressed as a ratio".
- l.347 paired cluster bootstrap [Field & Welsh 2007]. Supported. Node bootstrap and seed dependence [Owen 2007]: misplaced (edit E3).
- l.361 "iterative stratification" [Sechidis 2011]. Supported. This is the source of the iterative stratification algorithm for multi-label data.
- l.375 Leiden-CPM, resolution limit, consensus and persistence [Traag 2019; Traag 2011; Fortunato & Barthélemy 2007; Monti 2003; Lancichinetti & Fortunato 2012; Jeub 2018; Zheng 2021]. Each is supported:
  - Traag 2019: Leiden.
  - Traag 2011: CPM is resolution-limit-free.
  - Fortunato & Barthélemy: modularity resolution limit.
  - Monti: resampling consensus clustering.
  - Lancichinetti & Fortunato: consensus clustering of networks.
  - Jeub: multiresolution consensus.
  - Zheng: HiDeF persistence across resolutions.
- l.397 "HC3 heteroskedasticity-consistent standard errors" [MacKinnon & White 1985; Long & Ervin 2000]. Supported. MacKinnon & White introduced the HC2/HC3 small-sample estimators. Long & Ervin recommend HC3 for small samples.
- l.397 "Satterthwaite (Bell–McCaffrey) degrees of freedom" [Satterthwaite 1946; Bell & McCaffrey 2002]. Supported. Bell & McCaffrey pair bias-reduced linearisation with a Satterthwaite df approximation.
- l.399 "PatientID-blocked Freedman–Lane permutation" [Winkler 2014]. Supported. The paper covers Freedman–Lane and exchangeability blocks.
- l.399 "generalised Pareto tail approximation" [Knijnenburg 2009]. Supported. The paper fits a GPD to the permutation tail for small p.
- l.401:
  - BB selective FDR [Benjamini & Bogomolov 2014]: supported. The paper controls the average within-family FDR over selected families.
  - Simes combination [Simes 1986]: supported.
  - BH [Benjamini & Hochberg 1995]: supported.
  - Stage-wise testing [Van den Berge 2017]: supported (stage-wise screening and confirmation).
  - Holm for three ΔAP comparisons [Holm 1979]: supported.

Supplementary and cover letter: no citations. Supplementary Notes 2 and 4 name methods without references; see M-list below.

## Statements needing a citation but having none

- **N1 (main.md Methods "Post hoc power aid", l.415; supplementary.md Note 4, l.96):** "the empirical-Bayes posterior of the A estimate (Tweedie's formula …)" has no citation. Add Efron B. Tweedie's formula and selection bias. *J Am Stat Assoc* 2011;106(496):1602–1614. doi:10.1198/jasa.2011.tm11181 (verified in Crossref).
- **N2 (main.md §5.1, §5.4 and Methods H1: "neutral, by equivalence testing … 90% CI inside ±0.1"):** this is a two-one-sided-tests equivalence criterion and has no citation. Add Schuirmann DJ. A comparison of the two one-sided tests procedure and the power approach for assessing the equivalence of average bioavailability. *J Pharmacokinet Biopharm* 1987;15(6):657–680. doi:10.1007/BF01068419 (verified in Crossref).
- **N3 (main.md Discussion, l.268):** "the method was selected on hu.MAP training labels derived from Complex Portal" has no citation. Fischer 2025 states that its training set is "derived from manually curated protein complexes in the Complex Portal". Add [Fischer 2025].
- **N4 (main.md Discussion "What survives", l.285):** "several of these match known biology: nucleotide and haem synthesis in lymphoid lines, intestinal lineage transcription, and the cyclin E–CDK2 axis after RB1 loss…" has no citation for any of the three. The writer should either cite primary sources for each (to be verified before use) or soften to "are consistent with familiar pathways (post hoc reading)". I have not proposed specific papers because none were verified here.
- **N5 (main.md Introduction l.34):** "downstream complex-level analyses of cancer dependencies usually start from curated complex resources … or functional networks". The cited items are the resources, not the analyses. Add one or two examples of such analyses, or change "usually" to "often" (edit E7).
- **Optional (not required):**
  - Freedman & Lane 1983, the original Freedman–Lane procedure (currently covered by Winkler 2014).
  - Platt 1999 for Platt scaling (covered by Niculescu-Mizil & Caruana).
  - A reference for equal-mass ECE.
  - Software citations (statsmodels, leidenalg/igraph, NumPy).
  - Dempster 2021 at §4.5 "compressed Chronos dynamic range".
  - Seal 2023 at §2.3 "post hoc annotations from HGNC gene groups".
- **Supplementary:** Notes 2 and 4 name HC3, Satterthwaite df, Freedman–Lane, Benjamini–Bogomolov and Tweedie's formula without references. A single pointer such as "references for these methods are given in the main-text Methods" is enough.
- **Cover letter:** no citation is needed. A factual note outside citation scope: "All inputs are public DepMap, BioPlex, Ensembl, HGNC and hu.MAP releases" omits CCLE proteomics (Nusinow 2020), which the study uses for T2 tiers and co-abundance.

## Proposed edits for the writer (exact)

E1. main.md reference list: complete the truncated titles.
- `Benjamini Y, Hochberg Y. Controlling the false discovery rate. *J R Stat Soc B*` → `Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J R Stat Soc B*`
- `Dempster JM, et al. Chronos. *Genome Biol*` → `Dempster JM, et al. Chronos: a cell population dynamics model of CRISPR experiments that improves inference of gene fitness effects. *Genome Biol*`
- `McFarland JM, et al. Improved estimation of cancer dependencies from large-scale RNAi screens. *Nat Commun*` → `McFarland JM, et al. Improved estimation of cancer dependencies from large-scale RNAi screens using model-based normalization and data integration. *Nat Commun*`
- `Meldal BHM, et al. Complex Portal 2022. *Nucleic Acids Res*` → `Meldal BHM, et al. Complex Portal 2022: new curation frontiers. *Nucleic Acids Res*`
- `CRISPR–Cas9 essentiality screens. *Nat Genet* 2017` → `CRISPR–Cas9 essentiality screens in cancer cells. *Nat Genet* 2017`
- `Monti S, et al. Consensus clustering. *Mach Learn*` → `Monti S, et al. Consensus clustering: a resampling-based method for class discovery and visualization of gene expression microarray data. *Mach Learn*`
- `Szklarczyk D, et al. The STRING database in 2023. *Nucleic Acids Res*` → `Szklarczyk D, et al. The STRING database in 2023: protein–protein association networks and functional enrichment analyses for any sequenced genome of interest. *Nucleic Acids Res*`
- `Van den Berge K, et al. stageR. *Genome Biol*` → `Van den Berge K, et al. stageR: a general stage-wise method for controlling the gene-level false discovery rate in differential expression and differential transcript usage. *Genome Biol*`
- `*ECML PKDD* 2013:451–466. [153]` → `*ECML PKDD* 2013, LNCS 8190:451–466. [153]`

E2. main.md reference list ordering:
- Move the `Dyer SC, et al. Ensembl 2025…` line above `Field CA, Welsh AH…`.
- Move the `Kish L. *Survey Sampling*…` line below `Kim E, et al. Dynamic rewiring…`.

E3. main.md Methods "Uncertainty", replace:
`Protein (node) and stratified pair bootstraps are secondary, the latter anti-conservative; Monte-Carlo seeds 1–4 show seed dependence [Owen 2007].`
→
`Protein (node) [Owen 2007] and stratified pair bootstraps are secondary, the latter anti-conservative; Monte-Carlo seeds 1–4 show seed dependence.`

E4. main.md Methods data table, CRISPRGeneEffect row: `DepMap portal [Tsherniak 2017; Arafeh 2025]` → `DepMap portal [Meyers 2017; Dempster 2021 (Chronos); Arafeh 2025]`. Keep Tsherniak 2017 in the Introduction.

E5. main.md Methods "Post hoc power aid": `(Tweedie's formula on the z-scores` → `(Tweedie's formula [Efron 2011] on the z-scores`.
- Add to the reference list: `Efron B. Tweedie's formula and selection bias. *J Am Stat Assoc* 2011;106:1602–1614. doi:10.1198/jasa.2011.tm11181`
- Add it to LITERATURE.md §10 as M14.

E6. main.md §5.1, first sentence: `or nothing changes beyond ±0.1 Chronos units (neutral, by equivalence testing)` → `or nothing changes beyond ±0.1 Chronos units (neutral, by equivalence testing with two one-sided tests [Schuirmann 1987])`.
- Add to the reference list: `Schuirmann DJ. A comparison of the two one-sided tests procedure and the power approach for assessing the equivalence of average bioavailability. *J Pharmacokinet Biopharm* 1987;15:657–680. doi:10.1007/BF01068419`
- Add it to LITERATURE.md §10 as M15.

E7. main.md Introduction: `usually start from curated complex resources (CORUM, Complex Portal) or functional networks (STRING) [Tsitsiridis 2023; Meldal 2022; Szklarczyk 2023]` → `often start from curated complex resources (CORUM [Tsitsiridis 2023], Complex Portal [Meldal 2022]) or functional networks (STRING [Szklarczyk 2023])`. Optionally add verified examples of such analyses.

E8. main.md Discussion: `Second, the method was selected on hu.MAP training labels derived from Complex Portal, so` → `Second, the method was selected on hu.MAP training labels derived from Complex Portal [Fischer 2025], so`.

E9. main.md Discussion "What survives": either add verified citations for the three "known biology" readings, or replace `Read post hoc from the gene symbols, several of these match known biology:` with `Read post hoc from the gene symbols, several of these are consistent with familiar pathways:`.

E10. main.md Methods data table, Ensembl row, and reference list (recommended): replace Dyer 2025 with `Yates AD, et al. Ensembl 2026. *Nucleic Acids Res* 2026;54:D1053–D1060. doi:10.1093/nar/gkaf1239` (verified in Crossref). This paper describes release 116 as the last release on the current platform. If Dyer 2025 is kept, remove `[VERIFY]`.

E11. main.md reference list, remove the `VERIFY` tags: `[M7, VERIFY]`, `[M10, VERIFY]`, `[M3, VERIFY]`, `[M12, VERIFY]`, `[M5, VERIFY]`, `[M11, VERIFY]`, `[M9, VERIFY]` → `[M7]`, `[M10]` (or the E10 replacement), `[M3]`, `[M12]`, `[M5]`, `[M11]`, `[M9]`. Replace the preamble sentence `Methods references M1–M13 are in LITERATURE.md §10, added at drafting and not yet verified; those marked [VERIFY] there need checking before submission.` with `Methods references M1–M13 are in LITERATURE.md §10 and were verified against Crossref on 2026-10-29 (results/verify_manuscript/reference_check.md).`

E12. LITERATURE.md (PI or writer; not edited by this checker):
- Mark §10 as verified 2026-10-29 and drop the [VERIFY] tags.
- Complete the truncated titles in §9 #64 (Chronos), #84 (stageR); also #30 (hu.MAP 2.0), #35 (OpenCell) and #92 (GSVA) have short titles, but they are not cited in the manuscript.
- Add M14 (Efron 2011) and M15 (Schuirmann 1987) if E5 and E6 are adopted.

E13. Optional (journal style): add DOIs to the §9 entries in the manuscript reference list. They are already in LITERATURE.md §9.

E14. supplementary.md, Note 2 (association engine paragraph) and Note 4 (power aid): add `References for these methods are given in the main-text Methods.`

E15. cover_letter.md: `All inputs are public DepMap, BioPlex, Ensembl, HGNC and hu.MAP releases.` → `All inputs are public DepMap (including CCLE proteomics), BioPlex, Ensembl, HGNC and hu.MAP releases.`
