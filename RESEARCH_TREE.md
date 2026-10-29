# RESEARCH_TREE.md

Legend: [ ] planned · [~] in progress · [x] done · [-] abandoned/negative · [!] key finding

## Status summary (2026-10-28 ~22:30, M3: Part 1 test evaluation complete; FZ-C maps frozen; map variability and FZ-X in progress)
Part 1 is complete. One confirmatory look: zavg_BEF ITT AP 0.333 [0.256, 0.410], against GLS 0.212, Spearman 0.204 and Pearson 0.197 (ΔAP vs GLS +0.121 [+0.079, +0.159]). An independent agent reproduced it. For Part 2, the lockbox (FZ-LB) and the complex maps (FZ-C, 27fa9b2) are frozen on construction set A. zavg map: 367 complexes covering 3,642 genes. GLS map: 379 complexes covering 2,862 genes. The pre-registered bridge endpoint is **null** (Δ BioPlex precision −0.025 [−0.075, 0.031]). The deficit comes from locus-flagged clusters; excluding them gives +0.003. An independent review accepted FZ-C with notes: byte-identical rebuild, bridge reproduced exactly, no firewall breach found.

The association engine is merged and tested on synthetic data. An H3a FDR failure was traced to a non-registered GPD-tail rule plus 1,000 permutations; it is fixed to the registry rule (10,000 permutations, GPD only below 1e-4). BB error rate is 0.046 over 1,000 simulations; pooled FDP 0.068 is reported and not claimed (D6). Next: map-to-map variability (needed before Fig 2 / G2); FZ-X (registry fills, H1 events, lockbox-testable list); controls; H3a → H3b → H1 → H2 on A; lockbox confirmations; BENCH/BRIDGE-L.

## Previous status (M2)
Design review closed after 3 rounds (6 + 5 + 1 reviewers; docs/reviews/). TEST_PROTOCOL v1.2 + harness judged ready to freeze by the round-3 reviewer (34 tests; dry run reproduces Part 1 IFT look 1 to 1.5e-8). Freeze + single look scheduled as soon as the working tree is clean (construction-infra agent still writing files). Audit PASS; all prior test-file contact disclosed (D0, D1).

DESIGN v3.1: lockbox L (25%, PatientID-blocked, lineage × genotype strata) and construction set A; complexes from locked zavg_BEF on dependency-bearing genes with consensus/persistence clustering; GLS bridge (size-standardised BioPlex precision, protein bootstrap CI) as the main methods claim; tiers T1/T2/T3; directions H3a lineage atlas, H3b ~25 driver/genome-state contexts, H1 within-complex genetic interactions (paralog compensation / sensitisation / collapse / neutral), H2 sub-module divergence; 22 pre-registered controls mapped to directions with kill rule; registry + trial ledger + lockbox ledger. Next: G1 test look → FZ-LB → FZ-C (+bridge) → FZ-X → controls → directions.

## Decision notes
- D0 (M0, 2026-10-28): Test files arrived 13:55 UTC, after final Part 1 commit (12:55 UTC). PI read line counts (3,312 pos / 597,426 neg lines) and first 3 lines (format check) before protocol freeze; no scoring. Test evaluation will run only after docs/TEST_PROTOCOL.md is committed (freeze commit hash recorded with results). Part 1 code is used from a locked clone; any Part 1 code bug found later is reported, not fixed in the scored run.

- D1 (M1, 2026-10-28): The Part 1 audit agent, in an intermediate run, counted how many test pairs (per file) appear in Part 1's label-free top-pair table (zavg ≥ 3). This is equivalent to precision/recall at the z ≥ 3 rule on test pairs. Counts were deleted from all outputs; the PI never saw them; the z ≥ 3 rule pre-dates this (Part 1 manifest) and no protocol decision was informed by them. Disclosed in TEST_PROTOCOL and the paper.
- D2 (M1): User will provide DATA_REQUEST items 1–7 and 10–13 (DepMap omics incl. signatures, BioPlex, fusions, DEMETER2, Ensembl paralogs, CCLE proteomics); no CORUM/Complex Portal → post hoc interpretation via HGNC gene_group + literature only.

- D3 (M2, 2026-10-28): Locked verdict computed on S1 (mapped-only, Part 1's evaluation universe); ITT AP is the primary *metric* and is quoted alongside; disagreement reported. No user escrow of test files (procedural guarantee + hashes reported to user).
- D4 (M2): Scope cut per PI review: H4 CN dosage, coherence rewiring, BioPlex-cell-line rewiring, BioPlex-integrated map are remaining ideas, not planned work. BioPlex is validation/tiering only (not construction) so it can serve as orthogonal evidence.
- D5 (M3, 2026-10-28): FZ-C accepted with notes (results/part2/verify/fzc_review.md). Disclosed: contexts.py read proteomics for the LoF functional check (18:57) before the FZ-C commit. DESIGN permits this and it does not feed construction. Before G2: compute map-to-map variability. "zavg more stable" stays descriptive (no rewired baseline). The bridge null is reported with the pre-written fail paragraph.
- D6 (M3): H3a/H3b error control is claimed as Benjamini–Bogomolov selective FDR (average over selected context families), not pooled FDP. The single-seed `false_pos <= 1` assertion was removed because BB does not guarantee it (holds in 69% of seeds); multi-sim FDR tests stay (results/part2/verify/h3a_fdr_diagnosis.md).
- D7 (M3, 2026-10-28): H1 amendments, decided from event counts and labels only, before any H1 outcome (docs/registry/AMENDMENTS.md). A1: the 318 silencing events (TPM<1 rule; lineage-confounded, trivially non-dependent, would dilute the panel DESIGN expected at 10–20) become the exploratory family H1S, never headlined. A2: an event is genetic loss if ≥50% of its loss models carry LoF or hom-del (DESIGN §2.2 counts hom-del as LoF), which moves SMARCB1 to H1. Pre-declared before the H1 run: a CN-covariate sensitivity for SMARCB1, which keeps the base covariate set in the primary analysis.

## Tree
- [x] M0 Setup: PROMPT.md, DATA.md, DATA_REQUEST.md, git, locked Part 1 clone, venv
- [x] M1 Literature → LITERATURE.md (161 refs verified; 8 ranked directions D1–D8)
- [x] M1 Part 2 data profile (results/profile/part2_profile.md): 1,361 models, 23 lineages ≥15; Model.csv has almost no molecular annotation; z≥5 graph 163k edges with non-essential outlier hubs → mutual-kNN filter needed
- [x] M2 Design reviews
  - [x] R1: 6 reviewers (docs/reviews/round1_*) → DESIGN v2, TEST_PROTOCOL v1.1
  - [x] R2: 5 reviewers → DESIGN v3, TEST_PROTOCOL v1.2
  - [x] R3: final reviewer → DESIGN v3.1 (3 blocking fixes); protocol ready to freeze
- [x] P1 Test-set evaluation
  - [x] Part 1 audit → docs/part1_audit.md: PASS (timeline, references, content, hyper-parameters); test set pair-disjoint from train; 1,121 test proteins (329 shared with train)
  - [x] TEST_PROTOCOL v1.2 frozen 039a843 (marker sha 524c4d68…); intent ab674b7; started 7afcf06 (before labels); results 0861e59 (test_results.json sha e7617d79…); runtime 33 min eval + 3.3 min scoring
  - [!] Single look, confirmatory: ITT AP zavg_BEF 0.333 [0.256,0.410] vs GLS 0.212, Spearman 0.204, Pearson 0.197 (B 0.310; post hoc PC-removed Pearson 0.239). ΔAP vs GLS +0.121 [+0.079,+0.159] (+57%), vs Spearman +0.129 [+0.023,+0.216], vs Pearson +0.136 [+0.039,+0.220]. AUROC 0.789 vs 0.730/0.756/0.749. Locked "clearly beats" verdict PASS on S1 (AP 0.349) and ITT; B fails vs GLS (as on IFT). Macro per-cluster AP ≈ GLS (0.297 vs 0.292): gain concentrated in large complexes, as in Part 1. Rule b P/R/F1 0.553/0.328/0.411.
  - [x] Independent re-derivation REPRODUCED (c69d828 numbers committed before comparison; 8d009f2): verdict-relevant numbers within 2.6e-15; 36 secondary calibration fields differ ≤4.1e-7 (ECE bin remainder rule; optimiser tolerance) — explained; 1% raw re-scoring agrees ≤3.5e-6
- [ ] P2 Context-specific complex dependencies
  - [x] Construction infrastructure (29 tests; bit-exact zavg reproduction; toy ARI 0.875, null → 0 clusters)
  - [x] FZ-LB lockbox (bd865f8): L=338, A=1,007; contexts table (88 s); stats machinery (HC3 = statsmodels; FL type-I 0.05; BB-FDR 0.0975 @q=0.1)
  - [x] Registry draft (18c17bb): 29 entries (23 controls, H3a, H3b, H1, H2, BENCH, BRIDGE-L), pending FZ-C/FZ-X fills
  - [x] FZ-C (27fa9b2; tiers/bridge 9498d05; merged 262e7b6): selected z0=6, k=10 (split-half 0.343). zavg 367 complexes / 3,642 genes / 132 sub-modules, 28% locus-flagged. Tiers: 159 T1 / 74 T2 / 124 T3-neg / 10 T3-nc. Runtime: score runs 70 min, selection 69 min, map build ~4 min.
  - [-] Bridge (primary methods endpoint): ΔBioPlex precision zavg − GLS −0.025 [−0.075, 0.031], null. Non-locus +0.003. Protein co-abundance favours GLS (−0.061 [−0.101, −0.018]). Pair-AP vs BioPlex: zavg 0.0345, Pearson 0.0223, GLS 0.0215.
  - [x] FZ-C independent review (e0d17c0): ACCEPT-WITH-NOTES; byte-identical rebuild (444 s)
  - [x] Map-to-map variability (134a31c; PI re-ran analyse stage, identical): median best-match Jaccard of frozen complexes in bootstrap maps, zavg 0.69 (73% ≥0.5; 3–4 genes 0.50, 17–60 genes 0.82) vs GLS 0.67 (71%); random size-matched sets 0.05–0.08 (96% of zavg complexes above random P95); pairwise ARI zavg 0.754 / GLS 0.747 (covered genes). Runtime 325 s + 32 s. Not independent of the frozen map (same subsamples); random baseline is a floor.
  - [x] Association engine (48d1076; synthetic). H3a FDR fix 02f62ee (BB rate 0.046, pooled 0.068 over 1,000 simulations)
  - [x] FZ-X (408593e; A1 84a0c43; A2 + fixes fb264b8; merged bba129c; manifest sha f23a5b51…). Lineage atlas (H3a): 23 lineages, m 8,441. Driver/genome-state contexts (H3b): 20 contexts, m 7,340; KMT2D and EWSR1 fusion go by the rare-context route; CREBBP, KMT2A-r, BCR::ABL1, PAX3/7::FOXO1 and SS18::SSX are dropped (<12 A models). H1: 24 genetic-loss events (m 28, 22 lockbox-testable). H1S: 317 silencing events (m 338, exploratory). BRIDGE-L: 8,717 + 7,577. Controls: PC7/8/18/21/22 not evaluable, so H2 has no controls and is exploratory. Freeze 23 s, deterministic.
  - [x] FZ-X independent review (f6e7ffd; delta ec3513e): ACCEPT. Independent recomputation matches exactly; no outcomes read; no trials or lockbox ledger.
  - [x] Controls, single look (intent e82b657 → results 4b5f800, 28 s). All 15 evaluable complex-level controls recovered (p 1e-88 … 0.023 for PC19a). Kill rule: H3a 2/2, H3b 10/10, H1 3/3 pass; H2 has none evaluable, so exploratory. PC16 is gene-level only and mixed (MLLT1 p=0.063; not counted, as registered).
  - [x] Controls independent review (fd17c32): ACCEPT-WITH-NOTES. Own HC3/Satterthwaite code reproduces 15 tests exactly. Notes: lineage levels with <2 models are merged into the reference (engine convention not stated in DESIGN §2.5; at most 2% change in p; disclosed in Methods). H1 control coverage is thin: PC6 and PC12 are single paralog genes, only PC15 is multi-gene.
  - [ ] A-side discovery, registry order: H3a → H3b → H1 (+ SMARCB1 CN sensitivity) → H1S → H2
  - [ ] Lockbox confirmations → BENCH / BRIDGE-L
- [ ] M-final Manuscript + peer review + REPORT.md
