# RESEARCH_TREE.md

Legend: [ ] planned · [~] in progress · [x] done · [-] abandoned/negative · [!] key finding

## Status summary (2026-10-28 ~17:00, M2: design review round 2)
Part 1 test look is prepared but not run: docs/TEST_PROTOCOL.md DRAFT v1.1 + harness (31 tests; dry run reproduces Part 1 IFT look 1 to 1e-8) revised after 3 round-1 reviews (freeze marker, two-step look with git-tracked ledger, full disclosure of prior test-file contact incl. D1). Freeze after round-2 statistics + AI-failure reviews. Audit PASS.

Part 2: all requested data arrived (DATA.md, 19 checksummed files). DESIGN v2 after 6 round-1 reviews: lockbox (25% models) + construction on set A; construction on dependency-bearing genes with persistence hierarchy; GLS bridge experiment as main methods claim; evidence tiers (BioPlex/proteomics); headline H1 within-complex genetic interactions (buffered vs collapsing), H2 sub-module divergence, H3 lineage atlas; corrected control table; registry + trial ledger. Round 2 (5 reviewers) running; construction infrastructure being built in parallel.

## Decision notes
- D0 (M0, 2026-10-28): Test files arrived 13:55 UTC, after final Part 1 commit (12:55 UTC). PI read line counts (3,312 pos / 597,426 neg lines) and first 3 lines (format check) before protocol freeze; no scoring. Test evaluation will run only after docs/TEST_PROTOCOL.md is committed (freeze commit hash recorded with results). Part 1 code is used from a locked clone; any Part 1 code bug found later is reported, not fixed in the scored run.

- D1 (M1, 2026-10-28): The Part 1 audit agent, in an intermediate run, counted how many test pairs (per file) appear in Part 1's label-free top-pair table (zavg ≥ 3). This is equivalent to precision/recall at the z ≥ 3 rule on test pairs. Counts were deleted from all outputs; the PI never saw them; the z ≥ 3 rule pre-dates this (Part 1 manifest) and no protocol decision was informed by them. Disclosed in TEST_PROTOCOL and the paper.
- D2 (M1): User will provide DATA_REQUEST items 1–7 and 10–13 (DepMap omics incl. signatures, BioPlex, fusions, DEMETER2, Ensembl paralogs, CCLE proteomics); no CORUM/Complex Portal → post hoc interpretation via HGNC gene_group + literature only.

## Tree
- [x] M0 Setup: PROMPT.md, DATA.md, DATA_REQUEST.md, git, locked Part 1 clone, venv
- [x] M1 Literature → LITERATURE.md (161 refs verified; 8 ranked directions D1–D8)
- [x] M1 Part 2 data profile (results/profile/part2_profile.md): 1,361 models, 23 lineages ≥15; Model.csv has almost no molecular annotation; z≥5 graph 163k edges with non-essential outlier hubs → mutual-kNN filter needed
- [~] M2 Design reviews
  - [x] R1: 6 reviewers (docs/reviews/round1_*) → DESIGN v2, TEST_PROTOCOL v1.1
  - [~] R2: 5 reviewers (pi, statistics, ai_failure, biology, algorithms)
- [ ] P1 Test-set evaluation
  - [x] Part 1 audit → docs/part1_audit.md: PASS (timeline, references, content, hyper-parameters); test set pair-disjoint from train; 1,121 test proteins (329 shared with train)
  - [~] TEST_PROTOCOL v1.1 draft + harness (dry run reproduces Part 1 IFT exactly); freeze after R2
  - [ ] Score test pairs: zavg_BEF + Pearson/Spearman/GLS; metrics + bootstrap CIs + calibration
  - [ ] Independent re-derivation
- [ ] P2 Context-specific complex dependencies
  - [~] Construction infrastructure (lean locked recompute, consensus/persistence clustering, toy sims)
- [ ] M-final Manuscript + peer review + REPORT.md
