# Round 3 (final): statistics + AI-agent failure modes

**Ran.** The harness pytest: 34 passed. I opened no test file, ran no look, and computed no metric or association.

## 1. TEST_PROTOCOL v1.2 + harness: **ready to freeze**

All round-2 MUSTs are resolved:
- **M-B1/P-M1.** "draft" occurs only in the header (TEST_PROTOCOL.md:3). A test applies the real header edit.
- **M-B2.** Escrow is dropped consistently (§8.4:393, DESIGN §9).
- **M-B3.** `started` is committed before `test_files()` (21_evaluate_test.py:253–257). The scan covers `--all` and `--reflog` (protocol.py:162). A real-git test covers reset and a deleted branch.
- **P-M2** at §7.3:313–315; **P-S1** at :383; **S-B4** at paths.py:15–16; **S-B5** at §8.2:352–354.

These need fixing in the pre-freeze edit, but none blocks:
- **Version label.** TEST_PROTOCOL.md:5 and :336, and 22_freeze_protocol.py:7, say "FROZEN v1.1" / "v1.1". They should say v1.2. The gate accepts either.
- **SHOULD: a failed `started` commit leaves the record on disk.** For example, `index.lock` held by a Part 2 agent makes protocol.py:173 fail. The way out is then reverting the ledger (an "erasure"), and §8.3:371 would make the look non-confirmatory although no label was read. Fix: restore `looks.jsonl` from HEAD on failure, and state that this spends nothing.
- **SHOULD: untracked Part 2 files.** Files under `src/`, `scripts/` and `tests/` (present now) make `require_frozen` refuse. Commit Part 2 work-in-progress before the freeze, and pause Part 2 agents during steps 2–4.
- **COULD.** Check the score-TSV SHA-256s against the committed run records, or have Verifier A do it.

## 2. DESIGN v3: three blocking items

**B1. The bridge inference is ill-defined. Fix before FZ-C (DESIGN.md:48).**
- Problems with "CI over the 20 bootstrap maps":
  - 20 replicates cannot give 2.5/97.5% percentiles.
  - 80% subsamples without replacement understate variance unless rescaled.
  - The spread measures per-run instability, not the uncertainty of the final consensus maps' precision.
- Problem with "matched gene coverage" alone: within-cluster precision falls with cluster size, so a different size distribution can decide "better" by itself. Algorithms item 5 required matched sizes.
- Change:
  - Use a paired protein-node bootstrap over the BioPlex-observed genes of the two final maps (B = 2,000, fixed seed).
  - Size-standardise precision: per size-bin values, weighted by the pooled size distribution.
  - Report the 20-run spread as stability only.

**B2. The lockbox stratification order is circular. Fix before FZ-LB (DESIGN.md:27 vs :15).**
- Problem: L is stratified on "the H1 event genes listed at FZ-X". Those genes need complexes, and the complexes are built on A = U∖L.
- Change: stratify on a genotype-only list fixed now from Omics, without CRISPR outcomes. For example: genes with MSS LoF or homozygous deletion in ≥16 U models.

**B3. A full-L read can come before the confirmations. Fix before the first lockbox read (DESIGN.md:17, :48, :62).**
- Problem: the complex-vs-gene benchmark and the bridge's lockbox-confirmation secondary read L for all complexes × contexts, for both maps. "One read per closed registry entry" lets them run first, which compromises every later "lockbox-confirmed" claim.
- Change: allow these reads only after all H3a/H3b/H1/H2 confirmatory reads are committed (ai_failure M-C1, FZ-L).

All other round-2 MUSTs are resolved.

## Non-blocking
1. DESIGN:46 and :52 cite Ensembl coordinates, but none are in DATA.md; `locus.py` uses the HGNC arm/band. State this fallback before FZ-C.
2. DESIGN:44 puts the γ grid on the frequencies, while algorithms item 3 puts it in stage 1 with γ₂ = 0.5. Write out the two-stage rule at FZ-C.
3. H3a and H2 each have only 2 controls. Adopt PI S3: fewer than 3 evaluable controls means "n/a", and the direction is exploratory.
4. Seed ARPACK `v0` per run, and make `reproduce.sh` compare edge lists, not matrix hashes.
5. Also report the H1 neutral margin (±0.1 Chronos) in residual-SD units.
