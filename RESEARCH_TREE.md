# RESEARCH_TREE.md

Legend: [ ] planned · [~] in progress · [x] done · [-] abandoned/negative · [!] key finding

## Status summary (2026-10-28, M0: setup)
Project = Part 1 locked-method evaluation on the withheld hu.MAP test set + Part 2 cancer-context complex dependencies + submission package. Locked Part 1 method: zavg_BEF (label-free; Part 1 commit c7f95ee, cloned read-only to external/part1_locked). Genome-wide pair scores exist in the Part 1 folder (outputs/all_pairs, 18,429 protein genes).

Setup done: git repo, PROMPT.md verbatim, DATA.md with checksums (DepMap/HGNC identical to Part 1), venv (Part 1 env + leidenalg/igraph). External data for Part 2 (DepMap mutations/CN/expression, BioPlex, post hoc CORUM) requested from the user in DATA_REQUEST.md. Next: literature review, test-set protocol freeze, Part 1 audit.

## Decision notes
- D0 (M0, 2026-10-28): Test files arrived 13:55 UTC, after final Part 1 commit (12:55 UTC). PI read line counts (3,312 pos / 597,426 neg lines) and first 3 lines (format check) before protocol freeze; no scoring. Test evaluation will run only after docs/TEST_PROTOCOL.md is committed (freeze commit hash recorded with results). Part 1 code is used from a locked clone; any Part 1 code bug found later is reported, not fixed in the scored run.

## Tree
- [x] M0 Setup: PROMPT.md, DATA.md, DATA_REQUEST.md, git, locked Part 1 clone, venv
- [ ] M1 Literature → LITERATURE.md
- [ ] M2 DESIGN.md + review rounds (≤3)
- [ ] P1 Test-set evaluation
  - [ ] Part 1 audit (no test access) → docs/part1_audit.md
  - [ ] TEST_PROTOCOL.md frozen (commit) before any score
  - [ ] Score test pairs: zavg_BEF + Pearson/Spearman/GLS; metrics + bootstrap CIs + calibration
  - [ ] Independent re-derivation
- [ ] P2 Context-specific complex dependencies
- [ ] M-final Manuscript + peer review + REPORT.md
