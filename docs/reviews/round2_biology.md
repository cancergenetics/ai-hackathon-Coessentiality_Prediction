# Round 2 review: biology

**What I ran.** Two label-free scratchpad scripts. Neither opens hu.MAP or computes any association between contexts and complexes.
- `r2/ctx_v2.py` (31 s) computes the v2 context sizes in U (1,345 models). n_A and n_L are the median and 5th percentile over 200 stratified, PatientID-blocked 75/25 splits (seeds 0–199).
- `r2/h1_pool.py` (25 s) computes the graph-filter marginals and genome-wide LoF counts.

**Round-1 status.**
- M1 (CN scale): resolved. The housekeeping median r is 1.04 in both WGD and non-WGD lines, so CN is ploidy-relative.
- M2 (LoF): resolved, apart from the long-format mutation data.
- M3: resolved; see item 6.
- M4: mostly resolved.
- M5: partly resolved; see item 9.
- M6: resolved by the lockbox.

## Controls (n_A median/p5; n_L median/p5)

| ID | n_A | n_L | Verdict |
|---|---|---|---|
| PC1 MSI-H→WRN | 59/53 | 20/14 | Correct; but see item 1 |
| PC2/3 TP53 | 337 WT / 666 mut | 111 / 221 | Correct. Exclude the 15 ambiguous single-hit lines; unsequenced lines are NA |
| PC4 RAS→SHOC2 | 203 | 66 | Correct. MRAS fails the graph filter (0.2%) |
| PC5 BRAF V600 | 66 (any hotspot) | 22 | **Not constructible** (item 3) |
| PC6 ARID1A→ARID1B | 37/32 | 12/7 | Correct |
| PC7 SMARCA4→SMARCA2 | 30/26 | 10/6 | Correct |
| PC8 STAG2→STAG1 | 16/13 | 5/2 | Correct. |
| PC9 MTAP→PRMT5/WDR77 | 148 | 49 | Correct. Add MAT2A (G) |
| PC10 FOCAD-del or MSI→PELO | 88 (FOCAD 32) | 30 | Correct. FOCAD-del lies within MTAP-del (47/48), so also test FOCAD-del vs MTAP-del with FOCAD intact |
| PC11 KEAP1→NFE2L2 | 14/10 | 4/1 | Sign correct; n is thin. Use NRF2-active (KEAP1 or CUL3 LoF, or NFE2L2 hotspot): 31 |
| PC12 VPS4B→VPS4A | 255 | 85 | Correct. Homozygous deletion n = 0, so this is partial loss of the 18q segment |
| PC13 female or LOY→DDX3X | 597 | 198 | Correct. DDX3Y expression < 1 is cleaner (650) |
| PC14 RB1→CDK4/6/CCND1 | 50 | 16 | Correct (−) |
| PC15 SMARCB1 or SS18::SSX→BRD9 | 16/12 | 5/2 | Correct |
| PC16 KMT2A-r→MEN1/MLLT1 | 10/7 | 3/0 | Correct. Drop menin-dependent NPM1c lines (frameshift calls) from the reference group |
| PC17 skin→SOX10/MITF | 51 | 17 | Use "melanoma, any site" (52). Skin has 8 non-melanomas, misses 8 uveal |
| PC18 WGD→KIF18A | 549 | 181 | Correct. Add TP53 status as a covariate |
| PC19 PTEN / PIK3CA | 83 / 85 | 27 / 30 | Split into PC19a and PC19b; drop the 7 lines that have both |

## MUST

**1. Collinear covariates.**
- Problem: §2.5 adds log burden and MSIScore to every mutation context. For PC1 and the MSI arm of PC10, the context *is* MSIScore > 20, so the covariate absorbs the effect and the control fails spuriously.
- Change: no covariate may be a function of the context. Hypermutation contexts get neither burden nor MSIScore.

**2. The H1 classes conflate mechanisms.**
- Problem: "more essential = buffering, less = collapse" is wrong in three ways.
  1. Tumour-suppressor complexes (KEAP1–CUL3, TSC1/2, NF2–Hippo, STK11) have positive baseline Chronos scores. When they collapse, members move toward 0, and the sign rule reads that as "more essential".
  2. Non-paralog members becoming more essential is residual or dosage sensitisation, not buffering.
  3. The mean S_C dilutes a single paralog's takeover; ARID1B is one of about 15 BAF members.
- Change:
  - Classify each member relative to its baseline magnitude: collapse means |effect| shrinks toward 0.
  - Use four classes:
    - paralog compensation, requiring an in-complex Ensembl paralog of g;
    - residual sensitisation, requiring breadth across non-paralog members;
    - collapse;
    - neutral.
  - The H1 statistic is member-resolved: the paralog contrast plus non-paralog breadth, not S_C. Use it as the complex-level readout for PC6–8, PC12, PC13 and PC15.

**3. PC5.**
- Problem: the hotspot matrix is gene-level, so V600 cannot be resolved. Class-3 BRAF mutants are RAS- and SHOC2-dependent, and "no effect / −" is not a testable outcome.
- Change: redefine as BRAF hotspot vs RAS hotspot within MAPK-mutant lines, expected −; or drop it.

## SHOULD

**4. Functional-LoF check.**
- Problem: biallelic loss of a non-redundant essential subunit is not viable.
- Change: pre-register that an H1 event is evaluable only if g's own dependency is attenuated in its LoF lines.

**5. The H1 pool is small.**
- 100 genes have MSS biallelic LoF in ≥16 U models. Only 27 of these pass the graph filter, and some are not subunits (HLA-A/B/C, CDH1, APC). That leaves about 15–25 events, many with n_L ≤ 5.
- Keep the paralog-expression arm.
- Stratify the lockbox by lineage × MSI × membership of the pre-registered H1 contexts. These are genotypes, not outcomes, so this does not peek at results.

**6. Fusion parsing.**
- Problem: Gene1 and Gene2 are stored as "SYMBOL (ENSG…)". Naive matching gives n = 0 silently.
- Change: assert these counts in U: EWSR1::FLI1/ERG 19, SS18::SSX 5, KMT2A-r 13, PAX3/7::FOXO1 7, BCR::ABL1 12.

**7. H2 is scale-confounded.**
- Problem: sub-module minus parent is a part–whole difference. Its spread tracks each line's dynamic range (NNMD, Cas activity), and that range varies by lineage.
- Change:
  - Compare against the sibling complement, residualised on the complement across A.
  - Make NNMD and CasActivity primary covariates.
  - Pre-register two H2 controls:
    - haematopoietic CDK6/CCND3 vs CDK4/CCND1;
    - ncBAF vs cBAF in SMARCB1-null or SS18::SSX lines.

**8. Map controls to directions** so the kill rule is defined:
- H1: PC6–8, PC12, PC13, PC15;
- H2: the controls in item 7;
- H3: PC17;
- H4: PC9, PC12.

  As written, H2 cannot be killed.

**9. Segment-matched null.** Apply it to every CN context (PC9, PC10, PC12), not only to H4.

## COULD

**10. Collapse is exploratory.** The collapse class has no established positive control, so label collapse findings exploratory unless the lockbox confirms them.

H1/H2 are well posed after items 2 and 7; nothing else needed.
