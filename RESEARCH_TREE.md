# RESEARCH_TREE.md

Legend: [ ] planned · [~] in progress · [x] done · [-] abandoned/negative · [!] key finding

## Status summary (2026-10-28, M0: setup)
Project = Part 1 locked-method evaluation on the withheld hu.MAP test set + Part 2 cancer-context complex dependencies + submission package. Locked Part 1 method: zavg_BEF (label-free; Part 1 commit c7f95ee, cloned read-only to external/part1_locked). Genome-wide pair scores exist in the Part 1 folder (outputs/all_pairs, 18,429 protein genes).

Setup done: git repo, PROMPT.md verbatim, DATA.md with checksums (DepMap/HGNC identical to Part 1), venv (Part 1 env + leidenalg/igraph). External data for Part 2 (DepMap mutations/CN/expression, BioPlex, post hoc CORUM) requested from the user in DATA_REQUEST.md. Next: literature review, test-set protocol freeze, Part 1 audit.

## Decision notes
- D0 (M0, 2026-10-28): Test files arrived 13:55 UTC, after final Part 1 commit (12:55 UTC). PI read line counts (3,312 pos / 597,426 neg lines) and first 3 lines (format check) before protocol freeze; no scoring. Test evaluation will run only after docs/TEST_PROTOCOL.md is committed (freeze commit hash recorded with results). Part 1 code is used from a locked clone; any Part 1 code bug found later is reported, not fixed in the scored run.

- D1 (M1, 2026-10-28): The Part 1 audit agent, in an intermediate run, counted how many test pairs (per file) appear in Part 1's label-free top-pair table (zavg ≥ 3). This is equivalent to precision/recall at the z ≥ 3 rule on test pairs. Counts were deleted from all outputs; the PI never saw them; the z ≥ 3 rule pre-dates this (Part 1 manifest) and no protocol decision was informed by them. Disclosed in TEST_PROTOCOL and the paper.
- D2 (M1): User will provide DATA_REQUEST items 1–7 and 10–13 (DepMap omics incl. signatures, BioPlex, fusions, DEMETER2, Ensembl paralogs, CCLE proteomics); no CORUM/Complex Portal → post hoc interpretation via HGNC gene_group + literature only.

## Tree
- [x] M0 Setup: PROMPT.md, DATA.md, DATA_REQUEST.md, git, locked Part 1 clone, venv
- [x] M1 Literature → LITERATURE.md (161 refs verified; 8 ranked directions D1–D8)
- [x] M1 Part 2 data profile (results/profile/part2_profile.md): 1,361 models, 23 lineages ≥15; Model.csv has almost no molecular annotation; z≥5 graph 163k edges with non-essential outlier hubs → mutual-kNN filter needed
- [~] M2 DESIGN.md v1 drafted; review rounds (≤3)
- [ ] P1 Test-set evaluation
  - [x] Part 1 audit → docs/part1_audit.md: PASS (timeline, references, content, hyper-parameters); test set pair-disjoint from train; 1,121 test proteins (329 shared with train)
  - [ ] TEST_PROTOCOL.md frozen (commit) before any score
  - [ ] Score test pairs: zavg_BEF + Pearson/Spearman/GLS; metrics + bootstrap CIs + calibration
  - [ ] Independent re-derivation
- [ ] P2 Context-specific complex dependencies
- [ ] M-final Manuscript + peer review + REPORT.md
