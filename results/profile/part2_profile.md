# Part 2 data profile (descriptive; no labels)

Script `scripts/05_profile_part2.py` (seed 0); total 44.7 s, peak RSS 2.71 GB. No huMAP files or curated complex DBs read.

## 1. Models

CRISPRGeneEffect: **1361 models x 18434 genes** (NaN 3.2%); all join to Model.csv (1193 Model.csv rows have no screen). 30 lineages, 87 primary diseases (31 with >=10), 182 subtypes (36 with >=10). ScreenGeneEffect has 1611 screens (rows `SC-id`, 1363 models; not linkable to ModelCondition).

| Lineage (n>=15) | n | Lineage (n>=15) | n |
|---|---|---|---|
| CNS/Brain | 133 | Lung | 126 |
| Lymphoid | 97 | Esophagus/Stomach | 96 |
| Bowel | 85 | Head and Neck | 77 |
| Skin | 75 | Pancreas | 71 |
| Breast | 67 | Ovary/Fallopian Tube | 62 |
| Bone | 59 | Soft Tissue | 50 |
| Peripheral Nervous System | 47 | Myeloid | 45 |
| Uterus | 38 | Biliary Tract | 36 |
| Kidney | 35 | Bladder/Urinary Tract | 34 |
| Liver | 25 | Pleura | 21 |
| Prostate | 19 | Cervix | 18 |
| Eye | 15 |  |  |


Lineages >=15 models: 23 (1331 models); >=30: 18 (1233 models). Primary diseases >=30: 15; >=15: 24; subtypes >=15: 29.

| Field | Distribution (CRISPR models) |
|---|---|
| Sex | {'Male': 740, 'Female': 555, 'Unknown': 66} |
| PrimaryOrMetastasis | {'Primary': 793, 'Metastatic': 423, 'NA': 119, 'Unknown': 13, 'Recurrent': 13} |
| AgeCategory | {'Adult': 993, 'Unknown': 194, 'Pediatric': 171, 'Fetus': 3} |
| PediatricModelType | {'False': 1116, 'True': 245} |
| GrowthPattern | {'Adherent': 938, 'Suspension': 169, 'Dome': 108, 'Mixed': 65, 'Unknown': 41, 'Spheroid': 40} |
| ModelType | {'Cell Line': 1253, 'Organoid': 108} |
| EngineeredModel (flag to exclude) | {'NA': 1345, 'drug adapted': 8, 'engineered clone': 5, 'genetic knockout': 1, 'genetic knock down': 1, 'genetic knock-in': 1} |
| Age | n=1126, median 54.0, range 0.0-90.0 |


**Molecular annotation in Model.csv.** No mutation/CN/expression columns. Only free text: ModelSubtypeFeatures n=223 (99 distinct), PatientSubtypeFeatures n=82; any feature: 258 models (19%). Keyword hits (feature fields):

| Keyword | n models | top lineages |
|---|---|---|
| KRAS | 4 | Myeloid:3, Lymphoid:1 |
| NRAS | 11 | Myeloid:8, Lymphoid:3 |
| PTEN | 3 | Myeloid:2, Bone:1 |
| TP53 | 58 | Myeloid:26, Lymphoid:14, Bowel:9, Pancreas:4 |
| MYC (not MYCN) | 33 | Peripheral Nervous System:18, CNS/Brain:8, Lymphoid:5, Myeloid:1 |
| ERBB2/HER2 | 21 | Breast:21 |
| ER/PR (ER+, PR+) | 14 | Breast:14 |
| TNBC | 31 | Breast:31 |
| MSI | 10 | Bowel:10 |
| MGMT | 22 | CNS/Brain:22 |
| KMT2A (MLL) | 15 | Myeloid:11, Lymphoid:4 |
| BCR-ABL1 | 12 | Myeloid:12 |
| EWSR1/EWS fusion | 17 | Bone:17 |
| '::' fusion (any) | 58 | Myeloid:24, Bone:18, Lymphoid:15, Soft Tissue:1 |
| amplified/amplification | 21 | Peripheral Nervous System:18, CNS/Brain:3 |
| deletion/del/loss | 4 | Myeloid:2, Bone:2 |
| mutation/mutant (text) | 11 | Myeloid:6, Lymphoid:3, Peripheral Nervous System:1, CNS/Brain:1 |


Subtype-name proxies: {'IDH-mutant': 19, 'IDH-wildtype': 85, 'BCR-ABL1+': 10, 'ALK positive/negative': 8, 'NUT': 3, 'HER2/ER/TNBC in name': 0, 'MSI/MSS in name': 0, 'p53/MYCN in name': 0, 'Ewing': 28, 'Melanoma (BRAF/NRAS-driven proxy only)': 76}. Conclusion: only a few binary contexts reach >=10 models (MSI, HER2+/TNBC/ER, MYC-amp, BCR::ABL1, EWS fusions, IDH-mutant/GBM names, MGMT), mostly lineage-nested; KRAS/BRAF/TP53 status etc. cannot be derived without an external mutation file (candidate DATA_REQUEST).

## 2. Confounders and ModelCondition

CRISPRConfounders: 2554 rows x 22 cols; covers all CRISPR models (0% missing for quality cols among CRISPR models; CasActivity 11%, DoublingTime 14% missing). Libraries: {'LibraryAvana': 1048, 'LibraryHumagneCD': 221, 'LibraryKY': 316, 'LibraryTKOv3': 4, 'LibraryBrunello': 5}; screen type {'ScreenType2DS': 1205, 'ScreenType3DO': 94}. Lineage(>=15) x library Cramer's V = 0.31 (p=4.1e-33).

| Covariate | eta^2 vs lineage |
|---|---|
| ScreenNNMD | 0.036 |
| ScreenROCAUC | 0.037 |
| ScreenFPR | 0.037 |
| ScreenMedianEssentialDepletion | 0.065 |
| ScreenMedianNonessentialDepletion | 0.053 |
| ScreenMADEssentials | 0.087 |
| ScreenMADNonessentials | 0.087 |
| CasActivity | 0.092 |
| ScreenDoublingTime | 0.034 |
| LibraryAvana | 0.139 |
| LibraryHumagneCD | 0.125 |
| LibraryKY | 0.122 |
| ScreenType3DO | 0.195 |
| pan-essential mean effect (derived) | 0.059 |


Effect-matrix PCA: var explained PC1-5 = [0.0302, 0.0188, 0.0177, 0.0139, 0.0112]; lineage eta^2 of PC1-5 = [0.3, 0.14, 0.3, 0.21, 0.06]; PC1-3 Spearman with NNMD [-0.05, -0.11, 0.18], with KY library [0.05, -0.02, 0.04]. Pan-essential strength per model vs NNMD rho=0.38. ModelCondition: 3385 rows, Broad/Sanger/other = {'BROAD': 1367, 'SANGER': 623, 'OTHER': 18}; models by source set {'BROAD': 743, 'BROAD+SANGER': 546, 'SANGER': 61, 'BROAD+OTHER+SANGER': 9, 'BROAD+OTHER': 2}; Sanger~library KY {'BROAD': {'n': 743, 'frac_KY': 0.0}, 'BROAD+SANGER': {'n': 546, 'frac_KY': 0.46}, 'BROAD+OTHER+SANGER': {'n': 9, 'frac_KY': 0.67}, 'SANGER': {'n': 61, 'frac_KY': 1.0}, 'BROAD+OTHER': {'n': 2, 'frac_KY': 0.0}}. Useful only as site/medium batch labels.

## 3. Genome-wide score matrix (zavg_BEF, n=18429 genes, 908 low-coverage)

| Quantile | 0.0001 | 0.001 | 0.01 | 0.05 | 0.1 | 0.25 | 0.5 | 0.75 | 0.9 | 0.95 | 0.99 | 0.999 | 0.9999 | 0.99999 | 0.999999 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| z | -1.92 | -1.55 | -1.16 | -0.84 | -0.68 | -0.40 | -0.09 | 0.25 | 0.72 | 1.19 | 2.25 | 4.91 | 13.85 | 21.40 | 27.97 |


| z | pairs | genes deg>=1 | >=3 | >=10 | median deg | max deg | edge frac low-cov | edge frac both pan-ess | edge frac both non-pan |
|---|---|---|---|---|---|---|---|---|---|
| >=3 | 667,677 | 18429 | 18429 | 18429 | 69.0 | 252 | 0.1122 | 0.0264 | 0.0992 |
| >=4 | 278,986 | 18428 | 18373 | 15822 | 26.0 | 124 | 0.0991 | 0.038 | 0.1532 |
| >=5 | 163,445 | 17257 | 14100 | 10072 | 15.0 | 100 | 0.0823 | 0.0428 | 0.1902 |
| >=6 | 114,513 | 13610 | 10917 | 8238 | 14.0 | 88 | 0.0679 | 0.0436 | 0.2113 |
| >=8 | 67,727 | 9883 | 8423 | 5788 | 12.0 | 62 | 0.0512 | 0.0419 | 0.235 |


Baseline low-coverage-any pair fraction 0.0961. Top hubs z>=5: [['ACP7', 100], ['SIRT2', 100], ['LRFN3', 98], ['WDR87', 94], ['ANKRD36B', 92], ['RASGRP4', 91], ['ZNF569', 89], ['SPAG5', 89]]. Degree histogram at z>=5: {'0': 1172, '1': 1654, '2-3': 2570, '4-9': 2961, '10-24': 4122, '25-99': 5948, '>=100': 2}.


Gene classes by frac(effect<-0.5): {'pan(>=90%)': 1014, 'selective(5-90%)': 2423, 'weak(1-5%)': 3179, 'never(<1%)': 11813}. Edges by class pair:

| z | pan(>=90%)-pan(>=90%) | pan(>=90%)-selective(5-90%) | pan(>=90%)-weak(1-5%) | pan(>=90%)-never(<1%) | selective(5-90%)-selective(5-90%) | selective(5-90%)-weak(1-5%) | selective(5-90%)-never(<1%) | weak(1-5%)-weak(1-5%) | weak(1-5%)-never(<1%) | never(<1%)-never(<1%) |
|---|---|---|---|---|---|---|---|---|---|---|
| 3 | 17605 | 24182 | 7722 | 18586 | 39500 | 38767 | 67838 | 29514 | 142410 | 281553 |
| 4 | 10610 | 14058 | 2925 | 4367 | 27033 | 21987 | 27564 | 16213 | 57393 | 96836 |
| 5 | 7001 | 9436 | 1508 | 1415 | 20204 | 14938 | 15076 | 11108 | 32968 | 49791 |
| 6 | 4995 | 6820 | 982 | 613 | 15974 | 11184 | 9722 | 8401 | 22918 | 32904 |
| 8 | 2837 | 4160 | 560 | 234 | 10667 | 6979 | 5130 | 5415 | 13378 | 18367 |


Median degree by class at z>=5: {'pan(>=90%)': 25.0, 'selective(5-90%)': 33.0, 'weak(1-5%)': 21.0, 'never(<1%)': 6.0}; hubs (deg>=50 at z>=5) by class: {'pan(>=90%)': 58, 'selective(5-90%)': 364, 'weak(1-5%)': 234, 'never(<1%)': 302}; Spearman(deg, gene SD)=0.381, (deg, frac<-0.5)=0.43.

| z | edges | nodes | components | largest 3 sizes | components >=3 | frac nodes in giant |
|---|---|---|---|---|---|---|
| 3 | 667677 | 18429 | 1 | [18429] | 1 | 1.0 |
| 4 | 278986 | 18428 | 1 | [18428] | 1 | 1.0 |
| 5 | 163445 | 17257 | 54 | [17144, 4, 3] | 7 | 0.993 |
| 6 | 114513 | 13610 | 197 | [13150, 9, 7] | 41 | 0.966 |


Leiden on z>=5 graph (17257 nodes, 163445 edges; seed 0, 2 iterations):

| Method | res | s | clusters | single | size>=3 | >=10 | >=50 | median(>=3) | largest | pan-ess-dominated (>=80%) |
|---|---|---|---|---|---|---|---|---|---|---|
| RBConfiguration | 0.25 | 0.94 | 71 | 0 | 24 | 18 | 18 | 459.0 | [4549, 2525] | 0.0 |
| RBConfiguration | 0.5 | 0.97 | 93 | 0 | 46 | 40 | 39 | 305.5 | [1722, 1308] | 0.0 |
| RBConfiguration | 1.0 | 1.01 | 117 | 0 | 70 | 64 | 63 | 191.0 | [1055, 871] | 0.0 |
| RBConfiguration | 2.0 | 1.01 | 149 | 0 | 102 | 96 | 95 | 140.0 | [945, 452] | 0.0 |
| RBConfiguration | 5.0 | 1.05 | 201 | 0 | 154 | 148 | 137 | 105.5 | [687, 259] | 0.0 |
| CPM | 0.05 | 1.19 | 2457 | 1180 | 946 | 363 | 77 | 7.0 | [162, 122] | 0.005 |
| CPM | 0.1 | 1.39 | 3159 | 1488 | 1260 | 384 | 52 | 5.0 | [144, 110] | 0.006 |
| CPM | 0.25 | 1.48 | 4383 | 1854 | 1458 | 377 | 22 | 4.0 | [100, 85] | 0.009 |
| CPM | 0.5 | 1.71 | 6063 | 3032 | 980 | 367 | 3 | 6.0 | [54, 52] | 0.023 |


## 4. Essentiality landscape

| Quantity | Genes |
|---|---|
| pan-essential (frac<-0.5 >= 90%) | 1014 |
| pan-essential (mean<-0.5) | 1712 |
| effect<-1 in >=15 models | 2404 |
| non-pan, effect<-0.5 in >=15 models | 5267 |
| selective: skew<-1, frac>=2%, non-pan | 680 |
| bimodal (BC>0.555), frac 5-95%, non-pan | 45 |
| lineage eta^2 >=0.2 / >=0.3 | 59 / 10 |
| frac<-0.5 in 0 (never < -0.5) | 2143 |
| frac<-0.5 in (0,1%) | 9673 |
| frac<-0.5 in [1,5%) | 3180 |
| frac<-0.5 in [5,25%) | 1282 |
| frac<-0.5 in [25,50%) | 484 |
| frac<-0.5 in [50,90%) | 658 |
| frac<-0.5 in >=90% (pan-essential) | 1014 |


## 5. Timings (s) and peak RSS (GB, cumulative max)

| Step | s | peak RSS GB |
|---|---|---|
| load_depmap_matrix | 3.49 | 0.72 |
| models | 0.14 | 0.73 |
| confounders | 0.1 | 0.73 |
| essentiality | 10.94 | 2.05 |
| score_matrix_scan | 18.45 | 2.62 |
| graph_components | 0.21 | 2.69 |
| leiden_z5 | 11.39 | 2.71 |
