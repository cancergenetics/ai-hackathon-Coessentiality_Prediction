# Round 1 review: biology (contexts, controls, directions, confounders)

**What I ran.** Two label-free scripts: `/private/tmp/claude-502/-Users-metinyazar--claude-spares-spare-gqTg9n/5a84c62a-d8a2-48ad-b26c-5dd43a72db34/scratchpad/ctx_feas.py` and `ctx_feas2.py`. They compute context sizes and file scales only. They read no dependencies, complexes or hu.MAP files, and test no association. Findings:
- `OmicsCNGeneWGS` is **linear** relative CN: median 1.02, 99th percentile 2.33, maximum 252. MTAP rel < 0.25 occurs in 222 lines, whose MTAP log1p TPM median is 0.1 (4.9 otherwise).
- MSIScore is bimodal. 79 lines have MSIScore > 20.
- Damaging calls fall disproportionately in MSI lines (MSI = 6% of lines): ARID1A 51/146, KMT2D 52/126, CREBBP 27/67, EP300 30/86. TP53 has 889 damaging calls, so the matrix includes missense variants.
- 870 fusion names occur in ≥5 models. The most frequent are artefacts or read-throughs (TFAM--TFAM 323, ACIN1--ACIN1 272, KRT7--KRT8 140). Names are unordered (FOXO1--PAX3, ABL1--BCR, AFF1--KMT2A).
- Real fusion contexts: EWSR1::FLI1 16, SS18::SSX 5, KMT2A with a canonical partner ≈13, BCR::ABL1 13, PAX3::FOXO1 6.
- 9p21: 48 lines have a FOCAD deletion, 47 of them also MTAP-deleted. 217/222 MTAP deletions are co-deleted with CDKN2A. 16 lines are MTAP-null by expression without a deletion.
- 256/700 male lines have lost Y expression.
- DEMETER2 shares 549 lines with CRISPR via `CCLEName`.

## MUST

**1. CN scale and ploidy.**
- Problem: §2.1 uses log2 thresholds on a linear file. 740/1,045 lines are WGD, so if the values are ploidy-relative, rel 0.25 can mean one remaining copy.
- Change: compute absolute CN = rel × `Ploidy` (OmicsGlobalSignatures). Then define:
  - homozygous deletion: abs < 0.5, falling back to rel < 0.25, plus expression < 1 log1p TPM;
  - partial loss: 0.5 ≤ abs ≤ ploidy − 0.75;
  - amplification: rel ≥ 2; high-level: rel ≥ 3.

  Confirm the normalisation in the release notes and record it in DATA.md.
- Why: collateral-null and CYCLOPS-dosage biology differ.

**2. LoF = biallelic loss, with hypermutation control.**
- Change:
  - LoF = damaging==2, or damaging + copy loss, or homozygous deletion. For X-linked genes (STAG2, KDM6A, ATRX), any damaging call counts in males.
  - Keep hotspot contexts separate from LoF contexts.
  - Exclude MSI-H lines from non-MMR LoF contexts in the primary analysis; include them in a sensitivity analysis with a mutation-burden covariate.
  - Request `OmicsSomaticMutations.csv` (long format), which gives truncating vs missense and VAF.
- Why: about 35% of ARID1A/CREBBP/EP300 calls are likely MSI repeat passengers, and paralog synthetic lethality needs protein loss.

**3. Fusions.**
- Change: replace "≥5 models" with a pre-registered list, matched on the unordered Gene1/Gene2 pair, with self-fusions and read-throughs excluded:
  - EWSR1::FLI1/ERG;
  - SS18::SSX;
  - KMT2A with a canonical partner (MLLT3, AFF1, MLLT1, MLLT10, AFDN, ELL);
  - PAX3/7::FOXO1;
  - BCR::ABL1.

  Contexts with n < 8 are descriptive.
- Why: the current rule yields about 870 mostly artefactual contexts.

**4. Correct the §2.7 controls** (see the audit table). One has the wrong direction (KEAP1), one uses the wrong context gene (PELO), and four are ill-posed for CRISPR knockouts.

**5. Locus artefacts and regional co-deletion.**
- Change:
  - Flag clusters with ≥50% of members within 10 Mb (HGNC location) as locus clusters; they are not biological claims.
  - Drop members within ±10 Mb of a CN-context gene from the complex score.
  - Use segment-matched non-complex sets as the null for every CN context.
  - Name regional contexts by segment ("9p21 loss"). Attribute an effect to MTAP, CDKN2A or FOCAD only through conditional models or expression-null lines.
- Why: residual cutting toxicity and proximity bias (Lazar 2024) create co-essential "complexes" and CN associations.

**6. Construction absorbs context variance.**
- Problem: lineage modules (SOX10/MITF, haematopoietic TFs) cluster *because* lines of that lineage co-depend, so lineage p-values computed on all lines are optimistic. The same applies to D2 splits.
- Change: cross-fit the D1/D2 headline claims: build on lineage-stratified half A, test on half B, and swap.

## SHOULD

**7. Expression.**
- Keep fixed membership as the primary score, and make per-model exclusion of non-expressed members a sensitivity analysis. This reverses §2.3, because per-line composition aliases with paralog switching.
- For every hit, report the effect within lines that express all members, and label the hit "expression-explained" or "expression-independent".

**8. Proliferation and Chronos scaling.**
- Chronos fixes each line's common-essential median at −1, so effects on core machines (ribosome, proteasome, spliceosome, Pol II) are relative and partly zero-sum. Label them as such and require DEMETER2 concordance.
- Use doubling time as a covariate for OXPHOS and mitochondrial-ribosome complexes.
- Use TP53 status as a covariate for CN contexts, because TP53-WT lines have stronger cutting toxicity.

**9. Add high-n, clean controls.**

| Context | Expected result |
|---|---|
| TP53-WT (≈470) | MDM2/MDM4/PPM1D/USP7 dependency, and a growth gain on TP53/CDKN1A/USP28/TP53BP1 knockout |
| KRAS/NRAS hotspot (209/64) | SHOC2–MRAS–PPP1CB dependency |
| BRAF-V600 (95) | No SHOC2 dependency |
| RB1-null (≈70) | Loss of CDK4/6–CCND dependency |
| PTEN-null vs PIK3CA-mutant | PIK3CB vs PIK3CA |
| Female or loss of Y | DDX3X (DDX3Y paralog); KDM6A/UTY |
| VPS4B rel < 0.5 (113) | VPS4A |
| SMARCB1-null (≈15–20) | BRD9 and DCAF5 |
| WGD | KIF18A |

**10. Negative controls.** Olfactory-receptor genes are tandem arrays, so they control for artefacts, not for biology. Add dispersed, never-expressed gene sets as the biological null.

**11. Add contexts.**
- Sex and loss of Y; WGD, ploidy and CIN; TP53 status; MYC/MYCN amplification.
- Co-mutation covariates: KRAS/STK11/KEAP1 in LUAD, ARID1A/PIK3CA in endometrial cancer.

**12. D3 collapsing.**
- Problem: "any subunit disrupted" mixes opposite consequences (ARID1A loss → ARID1B; SMARCA4 loss → SMARCA2).
- Change: collapse by paralog position or sub-module, and require biallelic loss.

## COULD

**13. Priority.** Highest novelty × interpretability first:
1. D4: paralog switching, including sex/LOY, PI3K, SWI/SNF, cohesin and VPS4.
2. D3: restricted to biallelic loss in MSS lines. The SKI→PELO template is the most novel class.
3. D2, cross-fitted.
4. D5, with RNAi primary. In CRISPR, test only *trans* partner effects, because a cis knockout of a pan-essential gene is lethal at any copy number.
5. D1 as the foundation.
6. Deprioritise D6–D8.

**14. Presentation.**
- Give each hit a mechanism class: lineage addiction, paralog buffering, dosage/collateral, complex quality control, or complex–complex coupling.
- Give each hit a novelty tier.
- For each finding, show:
  - member graph (zavg vs BioPlex edges);
  - complex score by context;
  - member heatmap with LOSO;
  - CRISPR vs RNAi;
  - expression-stratified effect.

## Audit of claimed controls

| Claim | Verdict |
|---|---|
| ARID1A → ARID1B | Real (Helming 2014). Use biallelic loss in MSS lines. |
| SMARCA4 → SMARCA2 | Real. Truncating/biallelic only; lineage-nested. |
| KEAP1-LoF → NFE2L2 "loss of dependency" | **Wrong direction.** KEAP1/CUL3 loss or an NFE2L2 hotspot → **gain** of NFE2L2 dependency. Confounded with STK11/KRAS. |
| MSI → WRN | Real; the strongest control. |
| MTAP → PRMT5/WDR77, MAT2A | Real but modest in CRISPR (PRMT5 near pan-essential). 98% co-deleted with CDKN2A; nested in CNS and pancreas. Add RIOK1/CLNS1A. |
| STAG2 → STAG1 | Real. Zygosity depends on sex. |
| SS18::SSX → BRD9 | Real; n=5, so descriptive. Pool with SMARCB1-null lines. |
| KMT2A-r → MEN1/DOT1L | Real for canonical partners. Add MLLT1. |
| MLL–menin and Mediator kinase as myeloid one-vs-rest | Ill-posed. Menin dependency needs KMT2A-r/NPM1c. CDK8/19, MED12/12L and MED13/13L redundancy blunts single knockouts. |
| MYB, CBFB–RUNX1, SOX10/MITF | Real; expected to be expression-explained. |
| "HNF1B/HNF4A in GI" | Imprecise. HNF4A/CDX2 are GI; HNF1B is kidney. |
| EWSR1::FLI1 → "lineage modules" | Undefined. FLI1 itself is technical (sgRNAs cut the fusion). |
| PSMC2/SF3B1 partial loss | Real (Nijhawan 2012; Paolella 2017), but an RNAi phenotype. Expect near-null in CRISPR. |
| OXPHOS in AML | Confounded (suspension, medium, doubling time); not a control. |
| 9p21 → PELO–HBS1L | Goes through **FOCAD**, not MTAP. Use FOCAD deletion, SKI-subunit LoF or MSI. |
| CREBBP → EP300 | Real but modest. 40% of CREBBP calls are in MSI lines. |
