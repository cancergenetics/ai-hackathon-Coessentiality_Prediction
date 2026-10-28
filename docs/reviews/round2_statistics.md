# Round 2 review: statistics

Reviewer role: biostatistics. Objects: DESIGN v2 §2 and TEST_PROTOCOL DRAFT v1.1.

## What I ran
No hu.MAP test file was opened, and no context–dependency association was computed.
- **Scripts** (in the scratchpad): `r2_lockbox_power.py` (15 s) and `r2_lof_counts.py` (7 s).
- **Context sizes** in U (1,345 models), using the §2.2 definitions (approximate).
- **1,000 simulated lockboxes** (25% per lineage, PatientID-blocked), giving the distribution of n_L for each context.
- **Power:** discovery in A at α = 1e-4, and confirmation in L at one-sided 0.05.
- **Welch check:** t(5) noise with variance ratio 3.

## A. Protocol: not ready to freeze until two one-line fixes are made

My round-1 items M1, S-P1 and C-P1 are resolved, and so are the AI-failure MUST-1 and MUST-2.

**P-M1 (MUST). The freeze gate will refuse the protocol.**
- *Problem:* `protocol.py` rejects any text matching `re.search("draft", text, re.I)`. After the header edit, **line 540** still matches ("…a draft header and changed hashes are refused…"); I simulated this. The response's "verified" is wrong.
- *Change:* reword line 540, for example to "a pre-freeze header".

**P-M2 (MUST). The pair set behind the §7.3 readings is ambiguous.**
- *Problem:* the CI-excludes-0 statistic is computed on both ITT and S1, but §7.3 does not say which set selects the reading. That is a forking path at write-up.
- *Change:* add one sentence. "The reading is selected by the S1 locked rule and the S1 CI statistic. The quoted AP and ΔAP are ITT, then S1. If ITT and S1 disagree, S1 governs and the disagreement is reported."

**P-S1 (SHOULD). The verifier tolerance is too strict.**
- *Problem:* "CI endpoints equal" can fail spuriously, because float differences of about 1e-12 can arise.
- *Change:* tolerance ≤ 1e-9.

With P-M1 and P-M2 made, the protocol is **ready to freeze**. No other change is needed.

## B. DESIGN v2 §2

### Round-1 status
**Resolved:**
- M2 (burden/MSIScore, ploidy/WGD, MSS-primary)
- M3 (HC3, Freedman–Lane, GPD)
- M5, M6
- S3, S6, S8
- the other reviewers' statistics-relevant MUSTs

**Still open:**
- M4: the family structure is still loose (B-M2).
- S4: lost its best-subunit and hit-mapping arms (B-M4).
- §9 cites `round1_response.md`, which does not exist.

### Lockbox power

n_L is the median (5–95% over draws). Confirmation power is P(confirm | discovered).

| context | n_U | n_L | P(n_L ≥ 5) | confirmation power, d = 0.8 / 1.2 |
|---|---|---|---|---|
| MSI-H | 79 | 20 (15–26) | 1.00 | 0.96 / 1.00 |
| ARID1A LoF | 50 | 12 (8–18) | 1.00 | 0.86 / 0.99 |
| SMARCA4 LoF | 40 | 10 (6–14) | 0.98 | 0.77 / 0.95 |
| STAG2 LoF | 21 | 5 (2–8) | 0.62 | 0.39 / 0.56 |
| KEAP1 LoF | 18 | 4 (2–7) | 0.49 | 0.30 / 0.44 |
| SMARCB1 LoF | 18 | 4 (2–8) | 0.46 | 0.28 / 0.41 |
| KMT2A-r | 13 | 3 (1–6) | 0.22 | 0.13 / 0.20 |

- **Large contexts:** TP53, RAS, BRAF, PIK3CA, PTEN, RB1, MTAP, WGD and sex all have n_L ≥ 16 and power ≥ 0.92.
- **Welch at n1_L = 5, d = 1, variance ratio 3:** power 0.34, type-I error 0.044.
- **H1 LoF gene counts:** 102 genes have MSS LoF in ≥ 16 U models (about ≥ 12 in A). Of these, 74 are homozygous-deletion driven (segments). Only 24 reach that size by biallelic mutation alone.

### MUST

**B-M1. Lockbox-testability is a coin flip at the confirmatory boundary.**
- *Problem:* a context with n_A = 12 has n_L ≈ 4. Contexts with n_U of 16–21 are confirmatory in A, but they are testable in L with probability of only 0.2–0.6. This affects PC8, PC11, PC15 and PC16.
- *Change:*
  - Before drawing FZ-LB, stratify by lineage × a pre-registered list of rare contexts (n_U < 30), so that n_L = round(0.25 n).
  - Freeze the list of lockbox-testable contexts at FZ-X.
  - Label the remaining contexts "A-only"; they are ineligible for the strongest findings.
  - For lineage-nested contexts (EWS, KMT2A-r, SMARCB1), compute the minimum detectable effect from the within-lineage comparators.

**B-M2. Specify the hierarchical FDR.**
- *Problem:* the families, the within-family level, the "cross-direction claim" and the re-entry of permutation p-values are all undefined.
- *Change:*
  - **Budget:** one error budget per direction, at q = 0.05; no cross-direction claim.
  - **Families:** one family per context. The family p-value is the Simes p over complexes. Apply BH over contexts. Within the R selected contexts, apply BH at q·R/M (Benjamini–Bogomolov).
  - **Sub-modules:** only inside parents selected in the same context, at the analogous level.
  - **Permutation re-entry:** substitute p_perm for the screening candidates (BH < 0.1), then re-run BH on all p. Discovery means q = 0.05.

**B-M3. H1 classification statistics.**
- *Problem:* "neutral" means "not significant", which confuses lack of power with evidence of no effect. Three further problems:
  - A global shift in dependency in the LoF lines (growth rate, or screen quality: pan-essential strength vs NNMD, ρ = 0.38) produces spurious "collapse".
  - One paralog can drive "buffering".
  - Losses of essential subunits are observed only where they are tolerated (survivorship).
- *Change:*
  - **Estimand:** β = the effect of loss(g) on S_{C∖g}.
  - **Four classes:**
    - buffered or collapsing: directional BH;
    - neutral: the 90% CI lies within ±δ (TOST), with δ fixed now, for example 0.25 residual SD;
    - inconclusive: everything else.
  - **Calls:** a class is called only if β lies beyond the 2.5/97.5% quantiles of the same context's effects on random sets matched on essentiality, coherence and size.
  - **Paralog split:** report β separately for the paralogs of g and for the non-paralog members. A "complex-wide" call needs the non-paralog component.
  - **Enumeration:** enumerate the (complex, g) events and their m at FZ-X. Expect tens: most tumour suppressors fail the dependency-bearing filter. Pre-declare a membership rule for lost genes.

**B-M4. The complex-vs-gene benchmark needs one primary estimand.**
- *Problem:*
  - It names two primaries.
  - Its expectation is two-sided.
  - There is no gene-to-complex mapping.
  - Averaging members inflates A–L correlation trivially (Spearman–Brown).
- *Change:* the primary estimand is, per (complex, registered context), the L z of the complex score minus the L z of its best member chosen in A.
  - Summarise with a PatientID-blocked bootstrap over L models.
  - Stratify the expectation by breadth measured in A.
  - Secondary: the A–L correlation in residual-SD units against matched random sets, and gene hits mapped to unique complexes.

### SHOULD

1. **Winner's curse.** Rank by the A statistic × tier, and report L effects and CIs for *all* discovered hits, confirmed or not. Ranking by L and then reporting L re-introduces the curse.
2. **The L test.** Use exact within-lineage Freedman–Lane permutation when n1_L < 20. Use Bell–McCaffrey df for HC3.
3. **The kill rule.**
   - "Recovered" means the same sign with one-sided p < 0.05 in A.
   - Each control is mapped to its directions now.
4. **State the lockbox's limits.** It guards against noise and forking paths, not against confounding, which replicates in L.

No further change is needed beyond these items.
