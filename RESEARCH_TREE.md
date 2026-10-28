# RESEARCH_TREE.md

Legend: [ ] planned · [~] in progress · [x] done · [-] abandoned/negative · [!] key finding

## Status summary (2026-10-28 ~18:30, M2 done: design frozen for execution)
Design review closed after 3 rounds (6 + 5 + 1 reviewers; docs/reviews/). TEST_PROTOCOL v1.2 + harness judged ready to freeze by the round-3 reviewer (34 tests; dry run reproduces Part 1 IFT look 1 to 1.5e-8). Freeze + single look scheduled as soon as the working tree is clean (construction-infra agent still writing files). Audit PASS; all prior test-file contact disclosed (D0, D1).

DESIGN v3.1: lockbox L (25%, PatientID-blocked, lineage × genotype strata) and construction set A; complexes from locked zavg_BEF on dependency-bearing genes with consensus/persistence clustering; GLS bridge (size-standardised BioPlex precision, protein bootstrap CI) as the main methods claim; tiers T1/T2/T3; directions H3a lineage atlas, H3b ~25 driver/genome-state contexts, H1 within-complex genetic interactions (paralog compensation / sensitisation / collapse / neutral), H2 sub-module divergence; 22 pre-registered controls mapped to directions with kill rule; registry + trial ledger + lockbox ledger. Next: G1 test look → FZ-LB → FZ-C (+bridge) → FZ-X → controls → directions.

## Decision notes
- D0 (M0, 2026-10-28): Test files arrived 13:55 UTC, after final Part 1 commit (12:55 UTC). PI read line counts (3,312 pos / 597,426 neg lines) and first 3 lines (format check) before protocol freeze; no scoring. Test evaluation will run only after docs/TEST_PROTOCOL.md is committed (freeze commit hash recorded with results). Part 1 code is used from a locked clone; any Part 1 code bug found later is reported, not fixed in the scored run.

- D1 (M1, 2026-10-28): The Part 1 audit agent, in an intermediate run, counted how many test pairs (per file) appear in Part 1's label-free top-pair table (zavg ≥ 3). This is equivalent to precision/recall at the z ≥ 3 rule on test pairs. Counts were deleted from all outputs; the PI never saw them; the z ≥ 3 rule pre-dates this (Part 1 manifest) and no protocol decision was informed by them. Disclosed in TEST_PROTOCOL and the paper.
- D2 (M1): User will provide DATA_REQUEST items 1–7 and 10–13 (DepMap omics incl. signatures, BioPlex, fusions, DEMETER2, Ensembl paralogs, CCLE proteomics); no CORUM/Complex Portal → post hoc interpretation via HGNC gene_group + literature only.

- D3 (M2, 2026-10-28): Locked verdict computed on S1 (mapped-only, Part 1's evaluation universe); ITT AP is the primary *metric* and is quoted alongside; disagreement reported. No user escrow of test files (procedural guarantee + hashes reported to user).
- D4 (M2): Scope cut per PI review: H4 CN dosage, coherence rewiring, BioPlex-cell-line rewiring, BioPlex-integrated map are remaining ideas, not planned work. BioPlex is validation/tiering only (not construction) so it can serve as orthogonal evidence.

## Tree
- [x] M0 Setup: PROMPT.md, DATA.md, DATA_REQUEST.md, git, locked Part 1 clone, venv
- [x] M1 Literature → LITERATURE.md (161 refs verified; 8 ranked directions D1–D8)
- [x] M1 Part 2 data profile (results/profile/part2_profile.md): 1,361 models, 23 lineages ≥15; Model.csv has almost no molecular annotation; z≥5 graph 163k edges with non-essential outlier hubs → mutual-kNN filter needed
- [x] M2 Design reviews
  - [x] R1: 6 reviewers (docs/reviews/round1_*) → DESIGN v2, TEST_PROTOCOL v1.1
  - [x] R2: 5 reviewers → DESIGN v3, TEST_PROTOCOL v1.2
  - [x] R3: final reviewer → DESIGN v3.1 (3 blocking fixes); protocol ready to freeze
- [ ] P1 Test-set evaluation
  - [x] Part 1 audit → docs/part1_audit.md: PASS (timeline, references, content, hyper-parameters); test set pair-disjoint from train; 1,121 test proteins (329 shared with train)
  - [~] TEST_PROTOCOL v1.2 + harness ready; freeze + look pending clean tree
  - [ ] Score test pairs: zavg_BEF + Pearson/Spearman/GLS; metrics + bootstrap CIs + calibration
  - [ ] Independent re-derivation
- [ ] P2 Context-specific complex dependencies
  - [~] Construction infrastructure (lean locked recompute, consensus/persistence clustering, toy sims)
- [ ] M-final Manuscript + peer review + REPORT.md
