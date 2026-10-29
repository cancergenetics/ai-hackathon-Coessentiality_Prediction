# Supplementary Table S9: the 254 registered strongest findings, with post hoc annotations

Registered rule (DESIGN §2.5; results/part2/robustness/README.md): discovery q < 0.05; primary lockbox call confirmed, or rare-context route passed; not defect-sensitive (R7); not a locus cluster; robust to leave-one-out by primary disease (H3a) or lineage (H3b, H1); for H1 also headline-eligible (E2). Ranked by tier class, then |t_A|. Computed on A after the lockbox read (disclosed).

Counts: 207 H3a, 45 H3b, 2 H1; tiers T1 152, T2 65, T3-neg 37.

Annotation columns are post hoc and descriptive (none is a filter or a test): beats_all_matched_sets, control_recovery, loso_partially_testable, unit_essentiality_A and the scale_effect_* columns (scripts/57_annotate_strongest.py). The two reporting labels below are also post hoc (D16):
- lead_example (37 rows): H3a, T1, no control annotation, |residual z| >= 2.
- possible_lineage_scale_artefact (56 rows): H3a in Lymphoid or Myeloid, scale-effect-concordant, |residual z| < 2 (beta in line with the lineage-wide essentiality trend; median unit essentiality -0.63 against -0.30 over all 367 units).

Revision round 1 (POST HOC, descriptive; results/part2/revision, R-2 and R-3). The D16 labels above used a linear trend and applied the artefact label only in Lymphoid and Myeloid. The revised labels use a natural-spline trend per context and one rule in every lineage:
- rev_scale_artefact_possible (95 rows): sign as the trend predicts and |spline residual z| < 2, in any lineage.
- rev_lead (37 rows): H3a, T1, no control readout, |spline residual z| >= 2 and not rev_scale_artefact_possible. 26 of them lie beyond the trend in the trend's own direction (rev_concordant_nl); the rule cannot flag these, so a residual beyond the trend does not exclude a scale artefact.
- rev_lineage_less_reliable: the context is one of the three where never-expressed gene sets reach p < 0.05 most often (melanoma, Skin, Lymphoid; 0.27-0.33); 9 rev_lead rows. The flag does not mark the other lineages as reliable: the never-expressed rate of the other lead contexts is 0.01-0.22 against a nominal 0.05, and only Myeloid (0.01) is close to calibrated, so read the numeric rate (rev_never_expressed_frac). Revision round 2: the four leads named in the main text (Ovary/Fallopian Tube Z0035, 0.22; Bone Z0303, 0.19; Head and Neck Z0051, 0.19; Myeloid Z0254, 0.01) all lie beyond the trend in its own direction (rev_concordant_nl), so the trend rule cannot exclude a scale artefact for them; their main independent support is the DEMETER2 check (rev_D2_p_dir 0.0003-0.009).
- rev_D2_*: the same complex score and design on DEMETER2 (RNAi) for A models with RNAi data; not affected by Chronos dynamic range; no FDR.

Revised lead rows (rev_lead; full table: tableS9_strongest_findings.tsv). Members are listed in full in the TSV; beta in Chronos units, negative = more dependent in the context. 'D16 lead' marks the earlier label.

| Rank | Context | Unit | Members (first 8) | A β | L β [95% CI] | n L | Spline residual z | In trend direction | Never-expressed rate | DEMETER2 β (one-sided p; n) | D16 lead |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | Lymphoid | Z0016 | ATRN, BCAR1, CRK, CRKL, DOCK5, ELMO2, FERMT2, ILK… | +0.182 | +0.173 [+0.133, +0.213] | 24 | +2.67 | no | 0.27 (less reliable) | +0.145 (3.2e-10; 20) | yes |
| 2 | Lymphoid | Z0190 | E4F1, PSMB5, PSMB6, PSMB7, THAP1 | +0.484 | +0.565 [+0.403, +0.727] | 24 | +5.92 | yes | 0.27 (less reliable) | +0.062 (0.14; 16) | yes |
| 4 | Myeloid | Z0016 | ATRN, BCAR1, CRK, CRKL, DOCK5, ELMO2, FERMT2, ILK… | +0.138 | +0.169 [+0.093, +0.244] | 11 | +2.49 | no | 0.01 | +0.167 (1.8e-08; 12) | yes |
| 6 | Lymphoid | Z0196 | PPP1CC, PPP1R11, PPP1R2, PPP1R7, SUGT1 | +0.179 | +0.171 [+0.116, +0.225] | 24 | +2.28 | yes | 0.27 (less reliable) | +0.043 (0.0087; 16) | yes |
| 10 | Lymphoid | Z0363 | ACTR1A, DYNC1LI1, DYNC1LI2, DYNLRB1, NDE1 | +0.199 | +0.238 [+0.153, +0.323] | 24 | +2.43 | yes | 0.27 (less reliable) | +0.075 (0.026; 16) | yes |
| 13 | Myeloid | Z0001 | ABCB7, CIAO1, CIAO2B, CIAO3, CIAPIN1, FDX2, FXN, GFER… | +0.145 | +0.151 [+0.051, +0.251] | 11 | +2.13 | yes | 0.01 | +0.040 (0.055; 12) | yes |
| 15 | CNS/Brain | Z0254 | AGFG1, AP2A1, AP2B1, AP2M1, AP2S1, C10orf88, CCDC32, DNM2… | +0.138 | +0.091 [+0.022, +0.160] | 33 | +2.78 | yes | 0.17 | +0.078 (0.0016; 33) | yes |
| 17 | Lymphoid | Z0014 | BTAF1, CCNC, DRAP1, KAT2A, KDM5C, MED12, MED13, MED13L… | -0.200 | -0.174 [-0.251, -0.098] | 24 | -2.68 | yes | 0.27 (less reliable) | -0.082 (1.9e-05; 20) | yes |
| 19 | Head and Neck | Z0051 | CALM1, CCNJ, CHCHD2, EGFR, FKBP9, GAB1, GRB2, LRRC57… | -0.133 | -0.107 [-0.155, -0.059] | 19 | -2.79 | yes | 0.19 | -0.043 (0.00027; 18) | yes |
| 33 | melanoma | Z0071 | DENND4C, EHBP1L1, EXOC1, EXOC2, EXOC3, EXOC4, EXOC5, EXOC7… | +0.110 | +0.065 [+0.008, +0.123] | 19 | +2.28 | no | 0.33 (less reliable) | +0.027 (0.031; 17) | yes |
| 37 | melanoma | Z0126 | BCAS2, CCNL1, CDC5L, CTNNBL1, PPP1R8, PQBP1, PRPF19, SNRPA… | +0.101 | +0.066 [+0.000, +0.132] | 19 | +2.28 | no | 0.33 (less reliable) | +0.069 (0.002; 22) | yes |
| 41 | Soft Tissue | Z0080 | CCT4, MSTO1, PDCD5, PDCL3, PFDN1, PFDN4, PFDN5, TBCA… | -0.156 | -0.091 [-0.202, +0.019] | 12 | -3.79 | yes | 0.12 | +0.016 (0.77; 6) | yes |
| 42 | Kidney | Z0033 | ARMC7, CENATAC, DDX59, DHX9, FAM204A, PDCD7, PPIH, RBM48… | +0.127 | +0.070 [-0.018, +0.157] | 9 | +2.17 | yes | 0.19 | +0.065 (0.01; 11) | yes |
| 46 | Pleura | Z0067 | B3GALT6, B3GAT3, B4GALT7, EXT1, EXT2, EXTL3, FGFR1, FRS2… | -0.265 | -0.182 [-0.408, +0.044] | 5 | -4.90 | yes | 0.13 | -0.044 (0.13; 6) | yes |
| 49 | Myeloid | Z0254 | AGFG1, AP2A1, AP2B1, AP2M1, AP2S1, C10orf88, CCDC32, DNM2… | -0.189 | -0.239 [-0.380, -0.098] | 11 | -2.91 | yes | 0.01 | -0.120 (0.00079; 12) | yes |
| 51 | Bone | Z0023 | ADRM1, PSMB2, PSMC1, PSMC2, PSMC3, PSMC4, PSMC5, PSMC6… | +0.093 | +0.084 [+0.021, +0.147] | 15 | +2.26 | no | 0.19 | +0.032 (0.031; 12) | yes |
| 53 | Bone | Z0303 | SNAPC1, SNAPC2, SNAPC3, SNAPC4, SNAPC5 | -0.210 | -0.251 [-0.370, -0.132] | 15 | -3.70 | yes | 0.19 | -0.095 (0.0089; 11) | yes |
| 55 | CNS/Brain | Z0081 | CHMP5, CHMP6, HGS, PTPN23, SNF8, STAMBP, TFG, TSG101… | +0.106 | +0.103 [+0.034, +0.172] | 33 | +2.85 | no | 0.17 | +0.032 (0.033; 26) | yes |
| 56 | Pancreas | Z0283 | AUP1, DERL2, FAF2, SEC62, SEL1L, SLC35B1, SYVN1, TMEM41B… | +0.104 | +0.101 [+0.041, +0.162] | 18 | +2.25 | yes | 0.09 | -0.003 (0.55; 12) | yes |
| 58 | Breast | Z0009 | ACO2, BOLA3, DLD, DLST, FDX1, FH, GCSH, GLRX5… | -0.117 | -0.128 [-0.219, -0.037] | 17 | -2.22 | yes | 0.10 | -0.025 (0.046; 35) | yes |
| 66 | Bone | Z0001 | ABCB7, CIAO1, CIAO2B, CIAO3, CIAPIN1, FDX2, FXN, GFER… | -0.120 | -0.158 [-0.249, -0.068] | 15 | -2.00 | yes | 0.19 | +0.017 (0.77; 12) | yes |
| 67 | Pancreas | Z0059 | ARHGEF7, CD151, DLG5, EFR3A, GIT2, ITGA3, ITGB1, MARK2… | -0.114 | -0.080 [-0.149, -0.010] | 18 | -2.45 | yes | 0.09 | -0.052 (0.0046; 24) | yes |
| 76 | Lung | Z0174 | ANKRD49, APC, AXIN1, CSNK1A1, CTNNBIP1, FKBPL | +0.086 | +0.077 [+0.032, +0.122] | 31 | +2.97 | yes | 0.12 | +0.024 (0.018; 77) | yes |
| 80 | CNS/Brain | Z0146 | DDX21, DNAJC2, DRG1, HSPA14, NOLC1, NUCLEOLIN, PINX1, ZC3H15… | -0.078 | -0.063 [-0.106, -0.020] | 33 | -2.33 | no | 0.17 | -0.041 (0.073; 33) | no |
| 81 | melanoma | Z0190 | E4F1, PSMB5, PSMB6, PSMB7, THAP1 | -0.218 | -0.274 [-0.382, -0.167] | 19 | -4.03 | yes | 0.33 (less reliable) | -0.037 (0.096; 18) | yes |
| 83 | Peripheral Nervous System | Z0190 | E4F1, PSMB5, PSMB6, PSMB7, THAP1 | -0.171 | -0.184 [-0.340, -0.028] | 12 | -2.75 | yes | 0.11 | not evaluable | yes |
| 85 | Ovary/Fallopian Tube | Z0035 | AFDN, ARHGAP29, CDC42, CDH1, CDH2, CRB3, CTNNA1, CTNND1… | -0.140 | -0.202 [-0.298, -0.107] | 16 | -4.36 | yes | 0.22 | -0.065 (0.00041; 25) | yes |
| 97 | Lung | Z0001 | ABCB7, CIAO1, CIAO2B, CIAO3, CIAPIN1, FDX2, FXN, GFER… | +0.062 | +0.061 [-0.002, +0.124] | 31 | +2.08 | yes | 0.12 | -0.024 (1; 77) | yes |
| 100 | Breast | Z0253 | ADAT2, ADAT3, DCPS, DNLZ, DOHH, DOLK, FDXR, HEMK2… | +0.089 | +0.124 [+0.060, +0.187] | 17 | +2.04 | no | 0.10 | +0.022 (0.033; 35) | no |
| 105 | melanoma | Z0097 | BRF1, CRCP, POLR1C, POLR3A, POLR3B, POLR3C, POLR3D, POLR3E… | +0.094 | +0.088 [+0.018, +0.158] | 19 | +2.29 | no | 0.33 (less reliable) | +0.048 (0.00063; 18) | yes |
| 111 | Soft Tissue | Z0003 | BCS1L, COA5, COA6, COA7, COX10, COX11, COX15, COX16… | +0.081 | +0.117 [+0.051, +0.182] | 12 | +2.25 | no | 0.12 | +0.021 (0.27; 6) | yes |
| 114 | Myeloid | Z0040 | ATP6AP1, ATP6AP2, ATP6V0B, ATP6V0C, ATP6V0D1, ATP6V1A, ATP6V1B2, ATP6V1C1… | -0.112 | -0.119 [-0.210, -0.028] | 11 | -2.43 | no | 0.01 | -0.009 (0.37; 11) | yes |
| 124 | Lung | Z0150 | ACD, CTC1, DCLRE1B, PNKP, POT1, STN1, TEN1, TINF2 | +0.070 | +0.072 [-0.005, +0.150] | 31 | +2.32 | yes | 0.12 | +0.009 (0.14; 77) | yes |
| 127 | Kidney | Z0323 | ARL8B, SNAPIN, SPG21, VPS16, VPS18, VPS33A, VPS39, VPS41 | -0.143 | -0.086 [-0.193, +0.021] | 9 | -2.43 | yes | 0.19 | +0.002 (0.55; 14) | yes |
| 129 | Lung | Z0087 | ABI1, ACTR2, ACTR3, ARPC2, ARPC3, ARPC4, BRK1, CYFIP1… | +0.062 | +0.105 [+0.051, +0.159] | 31 | +2.04 | yes | 0.12 | +0.012 (0.16; 74) | yes |
| 132 | Kidney | Z0362 | RNASEH2A, RNASEH2B, RNASEH2C | -0.120 | -0.242 [-0.465, -0.020] | 9 | -2.04 | yes | 0.19 | not evaluable | yes |
| 133 | Soft Tissue | Z0336 | DDX1, PYROXD1, RTCB, RTRAF, TRMT6, TSLIG1, TSLIG2, TSLIG3B | -0.086 | -0.082 [-0.163, -0.001] | 12 | -2.11 | yes | 0.12 | +0.039 (0.96; 8) | yes |

D16 leads not kept by the revised rule: Bone Z0130 (rank 64); Myeloid Z0278 (rank 86).
