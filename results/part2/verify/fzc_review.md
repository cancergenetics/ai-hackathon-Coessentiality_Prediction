# FZ-C independent review

This review covers:
- the FZ-C map commit `27fa9b2`;
- the post-freeze tiers and bridge commit `9498d05` on branch `worktree-agent-abed910aca3ac8712`;
- the pre-freeze code commit `f0f940e`.

The reviewer did not produce this work. No producer file was changed. The reviewer's scripts are in `scripts/verify/fzc_*.py` and their outputs are in `results/part2/verify/fzc/`.

## Verdict: **ACCEPT-WITH-NOTES**

- **Maps.** Both maps rebuild byte-identically in a fresh process.
- **Bridge.** The primary endpoint, its CI, the sensitivities and the protein co-abundance secondary all reproduce exactly with an independent implementation.
- **Tiers.** Tier counts reproduce to within Monte-Carlo noise.
- **Firewalls.** I found no firewall breach in the code paths or artefacts.
- **Notes.** None of the notes below blocks FZ-X. Notes N1–N3 must be handled in the Fig 2 write-up.

## 1. Firewall and order

| Check | Result |
|---|---|
| Commit order | `f0f940e` (code, 19:58) → `27fa9b2` (maps and `MANIFEST.sha256`, 21:41:06) → tiers files written 21:41:25–21:42:06 → `9498d05` (21:44). `tiers/summary.json` and `bridge.json` record `map_commits = 27fa9b2` for both maps. |
| Guard in script 37 | `validation.require_frozen` must succeed before `load_bioplex` or `read_proteomics` run (`test_validation_requires_frozen_maps` passes). `require_frozen` needs the manifests at `HEAD`. No commit before `27fa9b2` contains them, so script 37 could not have read BioPlex earlier. |
| Construction code | Scripts 31, 34, 35 and 36, plus consensus, construct, graph, selection, attach, runsets, scores, locus and genefilter, never name or import BioPlex, proteomics or `validation`. I checked this by grep and by the import graph; `tests/test_part2_firewall.py`, `test_part2_lockbox_firewall.py` and `test_part2_validation.py` pass (14 tests). `src/testeval` (hu.MAP) is not imported by any Part 2 file. No STRING, CORUM or Complex Portal reference exists anywhere in `src/part2` or `scripts/3*`. |
| `contexts.py` / script 32 | These read CCLE proteomics at 18:57 (parquet cache mtime) for the LoF functional check. That use is pre-registered as a context input (DESIGN F3, §2.2). The output is per-gene own-protein levels for LoF calls. It is not co-abundance and does not feed construction (no construction module imports `contexts`). The literal rule "no proteomics read before FZ-C" is therefore not met, but DESIGN permits this use. **Harmless.** |
| Score blocks are A-only | All 59 `timings.jsonl` runs (29 zavg, 29 GLS, 1 Pearson) share `S_hash` 8c264ca513f8. Every block's `model_ids` is a subset of `A_models.txt` (1,007 models), with 0 overlap with `lockbox_models.txt` (338 models; A ∩ L = ∅). The run sets have 1,007, 806 (80% subsamples) and 502–505 (halves) models. Script 31 only accepts `runsets.run_names()`, which are all derived from manifest-verified A. |
| Other blocks in the cache | Five infra blocks (all, U and U-sub; 1,076–1,361 models, including L) exist. They have different S hashes, so the S-keyed loaders cannot have picked them up. They were created at 17:00–17:37, before FZ-LB was drawn (19:00). They are the scorer reproduction tests, which DESIGN sanctions. |
| What I could not verify | That no ad hoc process (REPL, deleted script) read BioPlex before 21:41. BioPlex atime is 16:19:46, which matches the DATA.md SHA-256 at download. APFS relatime-style atime cannot show reads after that. This is an unavoidable limit of trusting the process. |
| Restricted data | No `datasets/`, `data_cache/` or `external/` path is tracked by git. Only derived TSV and JSON files are committed. |

## 2. Determinism

Script: `scripts/verify/fzc_rebuild_maps.py`. It runs a fresh process with `PYTHONHASHSEED=0`, calls the committed `36_build_maps.stage_build` with the committed `selection.json`, writes only to the scratchpad, and uses 3 processes.

| File | zavg | GLS |
|---|---|---|
| clusters.tsv | identical (SHA-256) | identical |
| membership.tsv | identical | identical |
| attachments.tsv | identical | – |
| parameters.json | identical apart from `code_commit` | identical apart from `code_commit` |

- Counts: zavg 367 level-1 clusters, 3,642 genes, 132 sub-modules, 104 locus clusters. GLS 379 clusters, 2,862 genes, 108 sub-modules.
- Runtime: 444 s wall (zavg 208 s, GLS 235 s), peak RSS 2.9 GB, load average 5–11.
- Not re-run: the 58 score runs and the split-half selection (about 2.3 h of compute). The selection itself was re-checked from `split_half_per_split.tsv`:
  - mean statistic for zavg at (6, 10) is 0.3428 (SD 0.005);
  - (4, 10) gives 0.3407 and (5, 10) gives 0.3381;
  - all rewired baselines are 0.
- Reproduction note: running from a worktree needs `PART1_DIR` set explicitly. Otherwise `paths.PART1_ORIG` resolves to `.claude/worktrees/coessentiality_protein complexes`.

## 3. Independent tiers and bridge

Script: `scripts/verify/fzc_independent_bridge.py`.
- It does not import `validation.py` or script 37.
- It has its own BioPlex parser, network, strata, nulls, BH correction, D5 matching, bin weights, standardised precision and bootstrap.
- It reuses the producer's code only for S/essentiality and for the cached proteomics parse. Correlations are computed with pandas pairwise-complete correlation.

**Inputs.** Every input count reproduces:
- BioPlex union: 167,374 edges and 14,543 nodes;
- S genes observed in BioPlex: 5,267, with 35,565 edges among them;
- proteomics: 239 A models and 4,687 measured S proteins.

**Bridge.**

| | Producer | Reviewer |
|---|---|---|
| Primary: zavg / GLS / Δ | 0.1570 / 0.1822 / −0.0252 | identical to 1e-16 |
| Primary 95% CI (multinomial, seed 5000) | [−0.0745, 0.0308] | [−0.0745, 0.0308] (identical) |
| Primary 95% CI (different RNG stream) | – | [−0.0740, 0.0302] |
| Sensitivity: no matching | −0.020 [−0.054, 0.017] | identical |
| Sensitivity: locus excluded | +0.003 [−0.060, 0.074] | identical |
| Secondary: protein co-abundance | −0.061 [−0.101, −0.018] | identical; alternative stream [−0.102, −0.019] |

Coverage matching gives 252 zavg vs 379 GLS clusters, 2,862 genes each. Bin weights are 0.306 / 0.269 / 0.306 / 0.119. The decision is "no difference claimed", and it is confirmed.

**Tiers.** These use different seeds (7100/7200 + level) and a simpler rule for merging strata.

| | Producer (T1 / T2 / T3-neg / T3-nc) | Reviewer |
|---|---|---|
| zavg level 1 | 159 / 74 / 124 / 10 | 160 / 75 / 122 / 10 |
| zavg level 2 | 59 / 28 / 42 / 3 | 59 / 26 / 44 / 3 |
| GLS level 1 | 159 / 81 / 131 / 8 | 160 / 80 / 131 / 8 |

- Concordance: 494/499 zavg clusters and 485/487 GLS clusters get the same tier.
- Observed edge counts and the numbers of observed and measured members are identical for every cluster.
- All 7 discordant clusters have q values within about 0.01 of 0.05. This is Monte-Carlo noise.

I did not recompute the pair-level AP: it needs the score blocks only and its code is simple. Note that it covers S × S pairs (13.87 M), while DESIGN says "genome-wide". It should be labelled as S-restricted.

## 4. Spot checks

Script: `scripts/verify/fzc_spotcheck.py`, seed 6100. Output: `fzc/spotcheck.tsv`.

**Locus flag.** I recomputed the flag independently for every cluster: Ensembl 116 midpoints by Entrez ID, a 10-Mb window, and the HGNC arm. It agrees for 499/499 zavg and 487/487 GLS clusters, and window and arm shares are identical.

**Five random T1 complexes.** All are biologically coherent:

| Cluster | Content |
|---|---|
| Z0021 | MCM2–7, GINS1–4, CDC7/DBF4, ORC1/6, TICRR, MTBP, DONSON (replication origin licensing and firing). CDK8 is an odd member. |
| Z0061 | PRC2 (EZH2, EED, SUZ12, MTF2) + PRC1 (BMI1, RING1, PCGF1/3, PHC2, BCOR) + linker H1. Two complexes merged into one level-1 cluster (cf. producer note 7). |
| Z0109 | Canonical BAF (SMARCA4, SMARCB1, SMARCC1/2, SMARCD2, SMARCE1, ARID1A/B, DPF2, SS18); 35 edges. |
| Z0276 | ALG/DPM N-glycosylation pathway. Only 2 edges, q = 0.028: a pathway module that barely passes T1. |
| Z0289 | BORC/BLOC-1 (BORCS5/7, BLOC1S2, KXD1) + MTMR6/9. |

**Five random locus-flagged clusters.** All are correctly flagged as obvious copy-number/arm blocks:

| Cluster | Arm share | Window share |
|---|---|---|
| Z0005 | 17q 100% | 74% |
| Z0017.1 | 19q 100% | 100% |
| Z0199 | 10q | 100% |
| Z0225 (LCE cluster) | 1q | 100% |
| Z0233 | 9q | 100% |

## 5. Open items

**Coverage matching drops 115 zavg clusters (780 genes).** My probe found no bias against zavg:

| | Dropped zavg clusters | Kept zavg clusters |
|---|---|---|
| Standardised precision | 0.151 | 0.157 |
| Fraction T1 | 0.36 | 0.47 |
| Median stability | 0.41 | 0.71 |
| Locus share | 20% | 32% |

- The 115 dropped clusters are mostly small: 44 / 36 / 34 / 1 across the four size bins.
- Dropping random zavg clusters to matched coverage (50 draws) gives Δ = −0.020 (2.5–97.5%: −0.038 to −0.004).
- With no matching, Δ = −0.020.
- So matching by stability does not drive the null. The one direction of bias it adds is a slightly higher locus share among the kept zavg clusters (stable copy-number blocks are kept).
- **Documented limitation; nothing needs re-running.**

**Map-to-map variability over the 20 subsample maps (not done).**
- DESIGN §2.3 asks for it "descriptively" as part of the bridge report (Fig 2). It does not enter any decision, and FZ-X (context table, families, lockbox-testable list, H1 events) does not depend on it.
- **It can wait until after FZ-X, but it must be done before Fig 2 / gate G2 is closed.** Otherwise the deviation must be stated explicitly in the manuscript.

**No rewired stability baseline (not done).**
- Infra issue 2 shows a rewired chance stability of 0.84–0.90 for 3–4-gene clusters. That is *higher* than the observed values (zavg 0.49, GLS 0.36).
- Absolute stability is therefore not interpretable. The secondary claim "zavg is more stable in every size bin" has no chance reference and depends on how each graph was built.
- **Documented limitation; not needed before FZ-X.** The README and Fig 2 should present stability as relative and descriptive only, never as evidence for zavg.

## Notes for the write-up

These do not block FZ-X.

- **N1. Locus clusters drive the primary-endpoint gap.**
  - 104 of 367 zavg level-1 clusters are locus clusters (28%), against 84 of 379 for GLS (22%).
  - After matching, the locus share is 32% of the kept zavg clusters.
  - Excluding locus clusters gives Δ = +0.003.
  - Among non-locus level-1 clusters, zavg has 155 T1 out of 263 (59%) and GLS has 151 out of 295 (51%).
  - The pre-registered primary result is still a null, and it stands. The interpretation should say that zavg's extra copy-number/arm artefact load costs it, not its physical-complex content.
- **N2. T2 is not locus-aware.**
  - 14 zavg locus clusters are T2, because copy-number dosage makes neighbouring proteins co-abundant.
  - Locus clusters are already excluded from claims, but T2 counts that include them should not be reported as "complex support".
- **N3. Reporting fixes.**
  - Pair-level AP is S × S, not genome-wide.
  - The "stability favours zavg" statement needs the caveat about the missing chance baseline.
  - Some T1 clusters merge complexes (Z0061: PRC1 + PRC2) or are pathway modules (Z0276: N-glycosylation, 2 edges). "Complex" in the T1 label should be read as "physically supported module".
- **N4. Reproduction convenience.** Document `PART1_DIR` for runs from a worktree or clean clone. `paths.PART1_ORIG` defaults to `PROJECT.parent`.

## Reproduce

```
export PART1_DIR=".../coessentiality_protein complexes"; export PYTHONHASHSEED=0
.venv/bin/python scripts/verify/fzc_rebuild_maps.py --out <scratch>/rebuild --n-proc 3   # ~7.5 min
.venv/bin/python scripts/verify/fzc_independent_bridge.py --out <scratch>/indep           # ~45 s
.venv/bin/python scripts/verify/fzc_spotcheck.py --out <scratch>/spot                     # ~10 s
```
