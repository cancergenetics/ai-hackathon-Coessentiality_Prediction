# LESSONS.md — process rules (rule + evidence)

Inherited from Part 1 (external/part1_locked/LESSONS.md), still applicable:
1. Reviewers and non-evaluating agents get no test-label access; say "compute no metric on test pairs" in every brief until the test look. (Part 1 D1 contamination.)
2. Fix PYTHONHASHSEED and seeds for any graph/community code; regenerate in a subprocess to test determinism. (Part 1 Louvain tie-break bug.)
3. Run a power simulation before freezing any statistical gate. (Part 1 D2.)
4. Verifiers re-implement from raw data and sweep neighbouring configs. (Part 1 a2f8622/460fd36.)
5. Pre-register the exact look command before a held-out look; one look. (Part 1 IFT look #1 clean.)
6. Update the tree in the same commit as each result. (Part 1 tree went stale.)
7. torch and LightGBM in one pytest process can deadlock on OpenMP; separate processes. (Part 1 final suite.)

Project 2:
8. Audit briefs must forbid counting test pairs inside any score-derived list (top-k tables are metrics in disguise). (D1: audit agent computed test hits in the z≥3 table.)
9. Request all external data in one batch right after the literature review's dataset section, not piecemeal. (Two request rounds were needed on day 1.)
