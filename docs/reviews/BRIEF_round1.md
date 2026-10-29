# Design review round 1 — common brief

Project root: <PROJECT_ROOT>. Read, in order: PROMPT.md (fixed goal; do not change scope), DESIGN.md (v1, the object of review), docs/TEST_PROTOCOL.md (Part 1 test-look protocol, DRAFT, freezes after this round), docs/part1_audit.md, LITERATURE.md (esp. §1.2–1.3, §3.3, §4, §7), results/profile/part2_profile.md, DATA_REQUEST.md (files the user is providing: DepMap omics mutations/CN/expression/signatures/fusions, BioPlex 3.0, DEMETER2 RNAi, Ensembl paralogs, CCLE proteomics; not CORUM). Part 1 (locked, read-only): external/part1_locked (REPORT.md, DESIGN.md).

Rules:
- Do NOT open datasets/huMAP/Test_*.txt and compute no metric, score or count involving hu.MAP test pairs. Compute no association between complexes and contexts (that would be results-peeking). You may run small label-free/feasibility computations (e.g., timing, size of files, power simulation on synthetic data) if they sharpen a recommendation; say what you ran.
- Propose concrete changes that stay within the project's scope (numbered, each: problem → change → why), or state that none are needed. Prioritise: MUST (blocking), SHOULD, COULD. Max ~1,200 words.
- Write your review to docs/reviews/round1_<role>.md. Do not edit any other file. Do not commit.
- Reply with ≤200 words: your top 5 changes, each one line.
