# Which protein complexes do particular cancers depend on?

**What it does.** The method takes two inputs:
- the DepMap CRISPR knockout data: gene effects for 1,345 cancer cell lines;
- a locked score of how likely two genes' proteins share a complex, from the Part 1 project.

Three steps follow:
1. It builds a map of 367 protein complexes using three quarters of the cell lines.
2. It tests every complex against cancer contexts defined only from omics data: lineage, driver mutations and genome states such as whole-genome doubling.
3. It tests each finding once on 338 held-out cell lines (findings in contexts too rare for that test use a pre-registered rare-context route instead).

The output is a list of 254 registered "strongest findings", each labelled with its caveats.

**How it differs from the standard approach.** Most published work does one of two things. It tests genes one at a time and groups the hits into complexes afterwards. Or it tests complexes taken from curated catalogues such as CORUM, often tuning its clustering on those same catalogues. This method differs in four ways:
- It builds its complexes without any curated complex list.
- It pre-registers every test in a machine-readable registry.
- It keeps a quarter of the cell lines sealed until a single confirmation read.
- It compares every complex with matched gene sets: other complexes' fragments of the same size, essentiality and co-dependency coherence, plus never-expressed genes.

This design measures how often a finding replicates, and whether it belongs to the complex at all, rather than assuming either.

![Pipeline](figures/pipeline.svg)

**The six steps.** Each step's reasons are in [METHOD.md](METHOD.md#step1).

1. **Seal a lockbox.** Split the cell lines into a construction set (1,007) and a lockbox (338), keeping lines from the same patient together and balancing lineages and key mutations. Everything is fitted on the construction set; the lockbox is read once.
2. **Build the complex map.** Link genes whose pair score is high (z ≥ 6) and who are among each other's 10 best partners. Cluster the resulting graph many times on resampled cell lines (Leiden algorithm) and keep clusters that recur: 367 complexes. Flag clusters that sit on one chromosome stretch.
3. **Score complexes in contexts.** A complex's score in a cell line is the mean knockout effect of its members. Contexts come from omics only, and a member that defines the context (e.g. the mutated gene) is removed.
4. **Test associations.** For each complex and context, regress the complex score on the context with covariates such as screen quality. Get robust standard errors (HC3), permutation p-values, and a false-discovery rate controlled within each context family (Benjamini–Bogomolov).
5. **Calibrate.** Run the same tests on matched gene sets, never-expressed gene sets and random pseudo-contexts, to see how often non-complexes "find" the same contexts.
6. **Confirm once.** Test each discovery once in the lockbox, then apply a pre-registered rule to keep findings that are significant, confirmed, robust and not chromosome-locus artefacts.

**Results.**

| Question | Answer | Strength |
|---|---|---|
| Is the input pair score better than the standard co-essentiality scores? | Yes: AP 0.333 against 0.212 (GLS), 0.204 (Spearman) and 0.197 (Pearson) on the held-out hu.MAP 3.0 test pairs | One pre-registered look; the gain is in large complexes |
| Is the complex map better supported by physical interactions than a GLS-built map? | No: −0.025 [−0.075, +0.031] | Null |
| Do its context findings replicate more often than the GLS map's? | No: 0.514 against 0.578 | Negative |
| Do context findings replicate in held-out lines? | 52.5% of 1,145 lineage findings, about what statistical power predicts | Moderate |
| Are complexes more context-associated than matched gene modules? | No: 28.1% against 28.1% of lineage tests at p < 0.05 | The central negative result |

**When not to use it.** Do not read one finding as "complex X is a drug target in cancer Y" without a check outside DepMap's CRISPR data. Never-expressed genes "associate" with some lineages in up to a third of tests, so lineage findings can be effect-scale artefacts. The lockbox has been read, so new confirmatory claims need new held-out data. See [Limitations](METHOD.md#limitations).

## The three documents

| Document | What it answers |
|---|---|
| [METHOD.md](METHOD.md) | What the method is and does: steps with reasons, a worked example (SNAPC complex in bone cancers) traced through real data, pseudocode, code locations, results, evidence strength, limitations |
| [JOURNEY.md](JOURNEY.md) | How the AI agents arrived at it: eight turning points on a timeline, the research process, dead ends |
| [INSIGHTS.md](INSIGHTS.md) | What was learned: scientific, methodological and engineering insights, ranked |
| [REVIEW_NOTES.md](REVIEW_NOTES.md) | Discrepancies between documents and code, open questions and uncertain choices, for the human reviewer |

An interactive version of the first three is in [`explainer.html`](explainer.html): a single self-contained file to download and open in a browser, with a clickable pipeline, a step-by-step worked example and an expandable timeline.

**Provenance.** The project was run by AI agents: a principal-investigator agent, Claude Opus 5.5, with 55 subagents. Human input was minimal, and no result has yet been verified by a human. Every number comes from a committed result file or from running the code, and sources are cited inline. The method builds on the locked pair score of [ai-hackathon-Coessentiality](https://github.com/cancergenetics/ai-hackathon-Coessentiality) (Part 1).

## Reproducing the figures and the page

```bash
PART1_DIR="<path to the Part 1 working copy>" PYTHONHASHSEED=0 .venv/bin/python docs/explainer/worked_example.py
```

```bash
.venv/bin/python docs/explainer/make_figures.py && .venv/bin/python docs/explainer/make_journey_figure.py
```

```bash
PANDOC=pandoc .venv/bin/python docs/explainer/build_html.py
```

The worked example needs the DepMap files and the construction caches (`data_cache/`), neither of which is distributed; see `datasets/README.md`. It reads construction-set lines only, takes about 25 s and peaks at 3.7 GB of memory. Its per-line table (`data/worked_example_scores_A.tsv`, raw DepMap gene effects of the five SNAPC genes) and the permutation draws (`.npy`) are not distributed, because they hold row-level DepMap data. `make_figures.py` needs them, so run `worked_example.py` first.
