# ai-hackathon-Coessentiality_Prediction — held-out test of a co-essentiality co-complex score, and cancer-context complex dependencies in DepMap (AI hackathon project)

## Disclaimer

(Human-written text starts)

This is an exercise in assessing the bioinformatics capabilities of agentic AI. This exercise was performed by the Ryan lab at University College Dublin, collaborating as part of a hackathon run by the Beltrao Lab at ETH Zurich, where we used Claude to run bioinformatics projects. All results in this project were generated entirely by AI with minimal human intervention. This is not intended to be viewed as a traditional publication and bioinformatics tool.

This project is the second part of [ai-hackathon-Coessentiality](https://github.com/cancergenetics/ai-hackathon-Coessentiality), which developed a method to predict co-complex membership from gene co-essentiality. Here, the locked method is first evaluated on the hu.MAP 3.0 test pairs that were withheld during its development; it is then used to build protein complexes and to study complex dependencies associated with cancer contexts (lineage, driver mutations, copy number and other genome states). It was executed by Claude Opus 5.5 with delegation to Opus 5.5/Sonnet 5.5 subagents, over about 26 hours (28–29 October 2026) on a personal 13-inch MacBook Pro equipped with an Apple M2 chip and 24 GB of unified memory.

Currently none of the results have been verified by a human. We plan to verify the agent’s methodology and benchmarking, and we will update this section accordingly.

The rest of this repository is entirely AI generated. It includes the code, the pre-registered test protocol and analysis registry, all committed results, and an AI-written report and manuscript in the style of a journal article, with two rounds of AI peer review.

(Human-written text ends)

## What is here

**New readers:** start with [`docs/explainer/`](docs/explainer/README.md). It explains what the method does (with a worked example on real data: the SNAPC complex in bone cancers), how the agents arrived at it, and what was learned. `docs/explainer/explainer.html` is the same content as one interactive page; download it and open it in a browser.

| file / directory | contents |
|---|---|
| `PROMPT.md` | the project brief, verbatim |
| `REPORT.md` | final report: approach, test-set results, Part 2 findings, peer review and revision, limitations, follow-ups, runtimes |
| `RESEARCH_TREE.md` | research log: status summary, decision notes D0–D23, tree of what was tried |
| `DESIGN.md` | Part 2 design (v3.1, after three review rounds) |
| `LITERATURE.md` | literature review |
| `DATA.md` | data sources, versions, SHA-256 checksums and licences |
| `LESSONS.md` | process lessons for running agentic research projects |
| `manuscript/` | AI-written manuscript (`main.md`, `supplementary.md`), figures, tables, cover letter, data/code availability, peer-review reports and responses (`review/`), clean submission copy (`submission/`) |
| `docs/` | plain-language explainer (`explainer/`), frozen test protocol (`TEST_PROTOCOL.md`), Part 1 leakage audit, analysis registry and amendments (`registry/`), design reviews, retrospectives |
| `src/` | `testeval/` (test-set evaluation harness, built on the locked Part 1 code) and `part2/` (complex construction, contexts, association engine, statistics, lockbox, revision analyses) |
| `scripts/` | numbered pipeline scripts (`00_`–`62_`), figure and table scripts (`fig_*.py`, `tab_*.py`), manuscript checks, independent verification scripts (`verify/`) |
| `tests/` | unit, synthetic end-to-end and firewall tests (test-label and lockbox access, forbidden resources) |
| `results/` | committed results: the single test evaluation (`test_eval/`), its dry run, Part 2 construction, contexts, directions, lockbox confirmation, robustness, revision analyses, independent reviews (`part2/verify/`), manuscript number checks |
| `datasets/` | placeholders only: where to put the DepMap, DEMETER2, HGNC, hu.MAP, CCLE proteomics, BioPlex and Ensembl files (git-ignored) |
| `external/` | placeholder only: where `scripts/00_setup_part1.sh` puts the locked Part 1 method (git-ignored) |

## Results in brief

**Part 1 method on the withheld hu.MAP 3.0 test pairs** (one pre-registered look, `docs/TEST_PROTOCOL.md`; intention-to-treat set of 600,738 pairs, 0.55% positives; 95% paired cluster-bootstrap CIs, B = 2,000):

| method | AP (primary) | AUROC | precision | recall | F1 |
|---|---|---|---|---|---|
| **zavg_BEF (locked Part 1 method)** | **0.333** [0.256, 0.410] | 0.789 | 0.553 | 0.328 | 0.411 |
| GLS (Wainberg et al. 2021) | 0.212 [0.156, 0.273] | 0.730 | 0.567 | 0.205 | 0.301 |
| Spearman | 0.204 [0.106, 0.352] | 0.756 | 0.306 | 0.272 | 0.288 |
| Pearson | 0.197 [0.103, 0.338] | 0.749 | 0.273 | 0.282 | 0.277 |

ΔAP zavg_BEF − GLS = +0.121 [+0.079, +0.159]; the pre-registered verdict passed. The gain comes from pair-rich test clusters: macro per-cluster AP is no better than GLS's. Scoring all six methods took 179.7 s (zavg_BEF 112.6 s).

**Part 2 (cancer-context complex dependencies).** The locked score built a map of 367 complexes without curated complexes, STRING or hu.MAP. Neither bridge from pair accuracy to complex-level value held: the primary BioPlex bridge is null and the secondary held-out-line bridge favours GLS. Context associations of complexes (lineage, drivers, genome states) are real but mostly not complex-specific: coherent co-dependency modules of non-member genes associate as often. A registered rule leaves 254 strongest findings, held-out-line confirmed, each with post hoc caveats (lineage-scale artefact, RNAi concordance). Details, numbers and limitations: `REPORT.md` and `manuscript/`.

## Data

No data files are included in this repository. To reproduce, download the files listed in `DATA.md` and place them under `datasets/` as described in [`datasets/README.md`](datasets/README.md). All data, caches (`data_cache/`) and the Part 1 clone (`external/`) are git-ignored. `scripts/check_data.sh` verifies the checksums.

## Reproduce

Python 3.12.

```bash
pip install -r requirements.txt
bash scripts/00_setup_part1.sh            # locked Part 1 method -> external/part1_locked (see external/README.md)
PYTHONHASHSEED=0 python -m pytest tests -q
```

`scripts/00_setup_part1.sh` clones the Part 1 method. Without a local Part 1 repository it clones the public release, <https://github.com/cancergenetics/ai-hackathon-Coessentiality>, at commit `3c684c9`. That commit is verified to be identical to the locked commit `c7f95ee` for every Part 1 file used here. Part 2 also reads three things from a Part 1 working copy pointed to by `PART1_DIR`: the genome-wide scores, the cluster assignment and the hu.MAP train pairs (`external/README.md`).

The pipeline order, commands, seeds and runtimes are in `REPORT.md` (§11) and `RESEARCH_TREE.md`. Two steps were single-use by design: the test look (`scripts/21_evaluate_test.py`) and the lockbox read (`scripts/51_lockbox_confirm.py`). Their committed ledgers refuse a second run, so re-running either needs a fresh repository state and would be non-confirmatory. Figures and tables regenerate from committed results alone (`scripts/fig_*.py`, `scripts/tab_*.py`).

## Paths in committed records

Machine-specific absolute paths in committed logs, run records and notes are written as placeholders: `<PROJECT_ROOT>` (this repository), `<PART1_DIR>` (the Part 1 working copy, see `external/README.md`) and `<SCRATCH>` (an agent's temporary folder, not distributed). The exceptions are `results/test_eval/fitted_rules.json` and `results/test_eval/test_results.json`. Their SHA-256 is pinned by the frozen test protocol and the test-look ledger, so they keep the original path of the Part 1 clone.

## Data sources and attribution

The results are derived from:
- DepMap data (Broad Institute; CC BY 4.0): DepMap Public 26Q3 (Dempster et al. 2021, *Genome Biol* 22:343), DEMETER2 v6 (McFarland et al. 2018) and CCLE proteomics (Nusinow et al. 2020, *Cell*);
- hu.MAP 3.0 (Fischer et al. 2025, *Mol Syst Biol*);
- BioPlex 3.0 (Huttlin et al. 2021, *Cell*);
- HGNC;
- Ensembl 116.

Credit these sources when reusing the results. Committed results contain summary statistics and gene- or complex-level estimates only, not the source data.
