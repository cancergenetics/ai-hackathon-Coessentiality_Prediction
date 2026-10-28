# Part 2 complex construction on A (FZ-C)

This folder covers DESIGN v3.1 §2.3: the zavg_BEF and GLS complex maps on the construction set A (1,007 models),
their tiers, and the bridge endpoint.

| Item | Value |
|---|---|
| FZ-C map commit | `27fa9b275a83fd76a9a2f35107fc0b4f3aaa9ac4` |
| Pre-freeze decisions commit | `f0f940e` |
| Freeze manifest | `FZC_manifest.json` (SHA-256 of the map files, selection, DECISIONS.md and construction code) |
| Per-map manifests | `*/MANIFEST.sha256` |

The validation inputs (BioPlex, CCLE proteomics) were first read *after* the map commit, by
`scripts/37_validate_maps.py`. That script refuses to run unless both map manifests are committed and unchanged
(`validation.require_frozen`). `tests/test_part2_firewall.py` checks that no construction file names these inputs.
Other guarantees:
- Only A's CRISPRGeneEffect was used. `scripts/31` accepts only A run sets (`src/part2/runsets.py`, A read through
  the manifest-verified `assoc.a_models`).
- No hu.MAP file, context label, curated complex resource or lockbox outcome was read.

## Pipeline and files

| Step | Script | Output |
|---|---|---|
| Pre-freeze algorithm decisions D1–D6, toy-validated | `35_toy_decisions.py` | `DECISIONS.md`, `toy_decisions.{json,md,log}` |
| Score blocks: A, 20 PatientID-blocked 80% subsamples (seeds 1000–1019), 4 lineage-stratified split halves (seeds 2000–2003); zavg_BEF (locked config, null moments recomputed per run) and GLS (Part 1 dev variant); Pearson on A | `34_run_A_blocks.sh` → `31_score_blocks.py` | `runs/timings.jsonl`, `runs/*.log` (blocks are in `data_cache/part2/blocks`) |
| Split-half selection of (z0, k) and the consensus/persistence maps | `36_build_maps.py` | `selection/`, `zavg_map/`, `gls_map/`, `build_maps.log` |
| Tiers and bridge (after the FZ-C commit) | `37_validate_maps.py` | `tiers/`, `bridge/`, `validate_maps.log` |
| Descriptive summary (post hoc HGNC gene_group column, labelled) | `38_describe_maps.py` | `summary/` |

Each map directory holds:
- `clusters.tsv`: id, level, parent, size, γ span, support, stability, cohesion, locus flag, members as symbols and
  Entrez IDs.
- `membership.tsv`.
- `parameters.json`: per-run model-set hashes, edge counts, null moments, rule counts.
- `runtime.json`.
- `MANIFEST.sha256`.
- `attachments.tsv` (zavg only).

## Selection

The rule is DESIGN §2.3 with D6: the split-half statistic over R = 4 splits.

| (z0, k) | zavg statistic | GLS statistic (edge-matched to zavg) | Rewired baseline |
|---|---|---|---|
| (4, 10) | 0.341 | **0.315** | 0 |
| (4, 20) | 0.204 | 0.264 | 0 |
| (5, 10) | 0.338 | 0.297 | 0 |
| (5, 20) | 0.216 | 0.253 | 0 |
| (6, 10) | **0.343** | 0.263 | 0 |
| (6, 20) | 0.256 | 0.240 | 0 |

- **Selected: z0 = 6, k = 10** for zavg (bold in the zavg column). The primary GLS map uses the same S and k, with
  per-run GLS edge counts matched to the zavg graph of the same model set (all 20 runs attainable, 13,966–14,402
  edges per run).
- GLS's own split-half preference would be (4, 10) (bold in the GLS column). It is reported only.
- **Caveat (infra issue 3).** For k = 10 the criterion is flat: 0.338–0.343, with split SD 0.005–0.02. The choice
  among z0 values is not well determined. k = 10 is clearly preferred over k = 20.

## Maps (level 1 = complexes, level 2 = sub-modules; S = 6,022 genes)

| | zavg_BEF (FZ-C map) | GLS (bridge comparator) |
|---|---|---|
| Level-1 clusters | 367 (3,642 genes) | 379 (2,862 genes) |
| Sub-modules (parents with any) | 132 (74 parents) | 108 (65 parents) |
| Level-1 size: median (IQR, max) | 9 (5–13, 29) | 5 (4–9, 37) |
| Level-1 size bins 3–4 / 5–8 / 9–16 / 17–60 | 85 / 85 / 149 / 48 | 152 / 121 / 78 / 28 |
| Coverage of pan-essential S genes (≥ 90% of A < −0.5) | 594/965 (62%) | 410/965 (42%) |
| Coverage of selective S genes (5–90%) | 1,722/2,269 (76%) | 1,494/2,269 (66%) |
| Coverage of other S genes | 1,326/2,788 (48%) | 958/2,788 (34%) |
| Locus-flagged, level 1 / level 2 | 28.3% / 23.5% | 22.2% / 25.0% |

**Attachment (zavg, D3):**
- 2,563 gene–cluster pairs, covering 2,000 non-S genes.
- 11,499 genes were tested, 5.7 M tests, BH q < 0.01, so the expected number of false pairs is about 26.
- Not locus-guarded: 712 pairs, 532 genes. Only these may define H1 events.

**Rule counts.** Of the 519 persistent level-1 candidates in the zavg map:
- 62 were nested children and 66 were absorbed into parents;
- 9 were duplicates and 9 were overlap conflicts;
- 1 failed cohesion.

So nothing is "lost" silently.

## Tiers (`tiers/`; BH within level)

| | zavg level 1 | zavg level 2 | GLS level 1 | GLS level 2 |
|---|---|---|---|---|
| T1 (BioPlex) | 159 (155 non-locus) | 59 | 159 | 42 |
| T2 (protein co-abundance) | 74 | 28 | 81 | 24 |
| T3-neg | 124 | 42 | 131 | 38 |
| T3-nc | 10 | 3 | 8 | 4 |
| Expected false T1 (5% of T1) | 8.0 | 3.0 | 8.0 | 2.1 |

Coverage of the validation inputs:
- BioPlex 293T ∪ HCT116 union: 167,374 edges and 14,543 nodes. 5,267 of the 6,022 S genes are observed, with 35,565
  edges among them.
- Proteomics: 239 A models; 4,687 S proteins measured in ≥ 20 of those models.
- Among the 263 non-locus zavg complexes: 155 are T1, 60 T2, 47 T3-neg and 1 T3-nc.

## Bridge (`bridge/bridge.json`)

### Primary endpoint (pre-registered)

The endpoint is size-standardised within-cluster BioPlex precision of the level-1 clusters at matched gene coverage
(D5). The zavg map dropped its 115 least-stable clusters to match GLS's 2,862 genes.

| | zavg | GLS | Difference (zavg − GLS) |
|---|---|---|---|
| Precision | 0.157 | 0.182 | **−0.025** |
| 95% CI (protein-node bootstrap, B = 2000, seed 5000) | – | – | **[−0.075, 0.031]** |

The CI includes 0, so **no difference is claimed**. The pre-registered "zavg map agrees better with physical
evidence than the GLS-module standard" is **not supported**.

Per size bin:

| Size bin | zavg | GLS |
|---|---|---|
| 3–4 | 0.237 | 0.228 |
| 5–8 | 0.086 | 0.156 |
| 9–16 | 0.163 | 0.194 |
| 17–60 | 0.099 | 0.095 |

Sensitivities (difference and 95% CI):

| Sensitivity | Difference | 95% CI |
|---|---|---|
| No coverage matching | −0.020 | [−0.054, 0.017] |
| Locus clusters excluded | +0.003 | [−0.060, 0.074] |

### Secondary endpoints

**Protein co-abundance** (same standardisation, seed 5001):

| | zavg | GLS | Difference | 95% CI |
|---|---|---|---|---|
| Mean within-cluster r | 0.211 | 0.272 | −0.061 | [−0.101, −0.018] |

On this secondary endpoint GLS clusters are more co-abundant.

**Stability** (mean best-match Jaccard against the stage-1 runs, level 1). zavg is higher in every size bin:

| Size bin | zavg | GLS |
|---|---|---|
| 3–4 | 0.49 | 0.36 |
| 5–8 | 0.57 | 0.48 |
| 9–16 | 0.65 | 0.61 |
| 17–60 | 0.77 | 0.73 |

**Pair-level ranking against BioPlex** over all 13.87 M S × S pairs with both genes observed (prevalence 0.26%):

| Method | AP | AP / prevalence | AUROC |
|---|---|---|---|
| zavg_BEF | 0.0345 | 13.5× | 0.606 |
| GLS | 0.0215 | 8.4× | 0.579 |
| Pearson | 0.0223 | 8.7× | 0.595 |

So the locked score ranks physical pairs better than GLS and Pearson pair by pair. After clustering at matched edge
count and gene coverage, the advantage is not detectable in within-cluster precision.

## Largest T1 non-locus zavg complexes

Members are gene symbols. The group label in brackets is the post hoc HGNC gene_group (most frequent among
members). It is an annotation only and was not used anywhere upstream.

| Cluster | Size | Members (selection) | POST HOC HGNC gene_group (verbatim) |
|---|---|---|---|
| Z0001 | 29 | ABCB7, CIAO1, CIAO2B, CIAO3, CIAPIN1, FDX2, FXN, GLRX3, HSCB, LYRM4, MMS19, NFS1, NUBP1, NUBP2, … | Cytosolic iron-sulfur assembly components |
| Z0002 | 27 | AARS2, EARS2, FARS2, HARS2, LARS2, MARS2, NARS2, PARS2, RARS2, SARS2, TARS2, VARS2, YARS2, POLG, POLG2, POLRMT, TWNK, TFB2M, … | Aminoacyl tRNA synthetases, Class II |
| Z0003 | 27 | COX10/11/15/16/19, COX4I1, COX5A/5B, SCO1/2, CYC1, CYCS, UQCRB, UQCRC1/2, UQCRFS1, BCS1L, … | Mitochondrial respiratory chain complex assembly factors |
| Z0004 | 27 | PLK4, SASS6, STIL, CPAP, CEP135, CEP152, CEP295, POC1A/B, POC5, RTTN, PCNT, … | WD repeat domain containing |
| Z0008 | 26 | BUB1, BUB1B, BUB3, MAD1L1, MAD2L1, TTK, CDC20, MIS12, DSN1, NSL1, PMF1, KNL1, ZWINT, … | MIS12 kinetochore complex |
| Z0009 | 25 | OGDH, DLST, DLD, SDHA–D, LIAS, LIPT1/2, GLRX5, ISCA1/2, NFU1, BOLA3, … | Mitochondrial iron-sulfur assembly components |
| Z0010 | 25 | TRAPPC1/3/4/5/8/11/13, SEC23B, SEC24A, SAR1A, RAB1A/B, TMED2/9/10, SCFD1, … | Trafficking protein particle complex subunits |
| Z0012 | 24 | NDUFA/B/C/S/V subunits, TIMMDC1 | NADH:ubiquinone oxidoreductase supernumerary subunits |
| Z0011 | 24 | AURKB, INCENP, BIRC5, CDCA8, NCAPD2/G/H, SMC2/4, KIF20A/23, PRC1, RACGAP1, ECT2, … | Condensin I subunits |
| Z0014 | 23 | MED12/13/13L/15/16/19/23/24, CCNC, TADA1/2B/3, SUPT7L, SUPT20H, TAF5L/6L, KAT2A, SGF29 | Mediator complex |

Full table: `summary/zavg_top_T1_complexes.tsv`.

## Runtimes (serial unless stated; laptop, load average 5–20 from concurrent agents)

| Step | Time | Peak RSS |
|---|---|---|
| zavg_BEF on A (with the non-S × S attachment rectangle) | 159 s | 3.99 GB |
| zavg_BEF subsamples (20) | 124–174 s each, 3,018 s total | ≤ 3.1 GB |
| zavg_BEF halves (8) | 119–129 s each | ≤ 2.9 GB |
| GLS (29 runs) | 1–3 s each (A has only 4 missingness patterns) | ≤ 1.75 GB |
| Pearson on A | 2 s | – |
| All 58 score runs | 4,218 s | – |
| Split-half selection (48 jobs, 3 processes) | 4,115 s wall | ≤ 2.5 GB per worker |
| Map build: zavg / GLS (stage 1 on 3 processes) | 250 s / 242 s | – |
| Tiers + bridge + pair AP | 52 s | 2.8 GB |
| Toy decision validation | 2,529 s | – |

## Issues and limitations

1. **The bridge's primary endpoint is null.** The secondary protein co-abundance endpoint favours GLS. Fig 2's
   pre-written null version applies.
2. **Selection is flat across z0 at k = 10** (above). The statistic was applied as registered.
3. **Locus-flagged clusters are frequent** (28% at level 1). The flag excludes them from biological claims. With
   locus clusters excluded the bridge difference is +0.003 (CI spans 0).
4. **Coverage matching removes 115 of 367 zavg clusters**, because the zavg map clusters 27% more genes than GLS at
   the same edge count. Without matching the result is unchanged.
5. **Not done:**
   - map-to-map variability over the 20 single-subsample maps: now done (section below);
   - a rewired chance baseline for map stability (infra issue 2: uninformative for small clusters);
   - Infomap/MCL (COULD);
   - the H3 lockbox confirmation-rate secondary endpoint (after the confirmations, per F5).
6. **Reproduction.** Determinism of the full maps was not re-run in a fresh process. The code path is tested for
   determinism in `tests/test_part2_construction.py`, and the maps should be re-derived by an independent agent (F6).
7. **The toy hub-gene merges** (DECISIONS D1) suggest that some large zavg level-1 clusters may merge two complexes
   via shared or hub genes. Two candidates are Z0011 and Z0014. Reading them as CPC + condensin I + cytokinesis and as Mediator kinase
   module + SAGA is post hoc interpretation from gene symbols, not validated here. The sub-modules resolve some of these.

## Map-to-map variability over the 20 bootstrap subsamples (`variability/`, DESIGN §2.3, descriptive)

Script `scripts/39_map_variability.py` (added after the FZ-C freeze; frozen maps and producer code untouched; A model
set only; no validation data, external complex resource or lockbox model is read).

- **Bootstrap map b** = the frozen consensus + persistence + nesting pipeline (`construct.build_map`, same
  `ClusterConfig`, z0 = 6, k = 10, 5 Leiden seeds, cohesion on the full-A block) applied to the single graph of
  PatientID-blocked subsample b (seeds 1000–1019, cached blocks reused). So a bootstrap map has no across-subsample
  consensus; it shows how much one 80% subsample varies. GLS: per-subsample graph edge-matched to the zavg graph of the
  same subsample, as in the frozen build.
- **Pipeline check** (`variability/verify_frozen_rebuild.json`): rebuilding the zavg map from all 20 graphs with this
  code gives the 367 frozen level-1 clusters exactly (identical member sets).
- Level 1 only. Complexes are counted in size bins 3–4 / 5–8 / 9–16 / 17–60.
- **(a)** For each frozen level-1 complex, best-match Jaccard against each bootstrap level-1 partition (uncovered genes
  unassigned). Per-complex values: `*_per_complex_jaccard.tsv`.
- **Chance reference.** For every frozen complex, 100 random gene sets of the same size drawn uniformly from S (6,022
  genes; a second pool, the genes in the frozen map's complexes, gives the same numbers) scored the same way
  (seed 7000). "Stable" can be read against its 95th percentile. This baseline is dispersed random sets; it does not
  represent a set of co-clustered but unrelated genes (hubs), so it is a floor, not a null for the clustering.
- **(b, c)** Pairwise ARI/NMI over the 190 pairs of bootstrap maps on genes covered by both maps (`ari_both`). The
  sensitivity counts genes covered by either map, uncovered genes as singletons (`*_either_singletons`).

Best-match Jaccard, pooled over complex x bootstrap map (median, IQR, fraction ≥ 0.5) against random size-matched
sets from S (median, fraction ≥ 0.5, 95th percentile):

**zavg_BEF** (frozen level-1 complexes: 367, genes 3,642)

| Size bin | Complexes | Median | IQR | Fraction ≥ 0.5 | Random median | Random ≥ 0.5 | Random P95 | Complexes with median > random P95 |
|---|---|---|---|---|---|---|---|---|
| 3-4 | 85 | 0.50 | 0.16–0.75 | 0.52 | 0.08 | 0.000 | 0.17 | 0.88 |
| 5-8 | 85 | 0.60 | 0.33–0.75 | 0.64 | 0.08 | 0.000 | 0.12 | 0.96 |
| 9-16 | 149 | 0.75 | 0.56–0.88 | 0.83 | 0.06 | 0.000 | 0.08 | 1.00 |
| 17-60 | 48 | 0.82 | 0.68–0.91 | 0.93 | 0.05 | 0.000 | 0.07 | 1.00 |
| all | 367 | 0.69 | 0.44–0.84 | 0.73 | 0.06 | 0.000 | 0.12 | 0.96 |

**GLS** (frozen level-1 complexes: 379, genes 2,862)

| Size bin | Complexes | Median | IQR | Fraction ≥ 0.5 | Random median | Random ≥ 0.5 | Random P95 | Complexes with median > random P95 |
|---|---|---|---|---|---|---|---|---|
| 3-4 | 152 | 0.60 | 0.27–0.80 | 0.65 | 0.14 | 0.000 | 0.20 | 0.97 |
| 5-8 | 121 | 0.67 | 0.46–0.80 | 0.75 | 0.10 | 0.000 | 0.14 | 0.98 |
| 9-16 | 78 | 0.70 | 0.50–0.83 | 0.78 | 0.07 | 0.000 | 0.09 | 1.00 |
| 17-60 | 28 | 0.73 | 0.48–0.86 | 0.74 | 0.04 | 0.000 | 0.06 | 1.00 |
| all | 379 | 0.67 | 0.43–0.80 | 0.71 | 0.10 | 0.000 | 0.17 | 0.98 |

Pairwise agreement of the 20 bootstrap maps (median over 190 pairs, IQR in brackets) and agreement of each bootstrap map with the frozen map:

| | ARI (both covered) | NMI (both covered) | ARI (either, singletons) | NMI (either, singletons) | ARI vs frozen | NMI vs frozen | Level-1 clusters per bootstrap map | Genes covered |
|---|---|---|---|---|---|---|---|---|
| zavg_BEF | 0.754 (0.745–0.764) | 0.935 (0.933–0.937) | 0.626 (0.618–0.635) | 0.899 (0.898–0.901) | 0.843 (0.839–0.849) | 0.964 (0.962–0.965) | 348 (344–352) | 3912 (3877–3928) |
| GLS | 0.747 (0.735–0.763) | 0.939 (0.937–0.940) | 0.529 (0.518–0.543) | 0.899 (0.897–0.900) | 0.868 (0.854–0.873) | 0.970 (0.968–0.971) | 538 (531–558) | 3476 (3440–3504) |

Reading (descriptive only; nothing here enters a decision):
- A single 80% subsample recovers the typical frozen complex (zavg pooled median Jaccard 0.69, 73% of complex-by-map
  values ≥ 0.5), far above the size-matched random reference (median 0.06, none ≥ 0.5). Recovery rises with size and is
  weakest for 3–4-gene complexes, where the spread is widest.
- The frozen maps are partly built from these subsamples, so agreement with the frozen map is not independent of them.
  The pairwise ARI between bootstrap maps is the cleaner measure of subsample-to-subsample variability.
- The single-graph GLS maps have many more level-1 clusters than the frozen GLS map (median 538 against 379): the
  across-subsample consensus is what prunes them. For zavg the numbers are close (348 against 367).
- The README stability table above (mean best-match Jaccard against stage-1 runs) is a different quantity.

Runtime: bootstrap maps, 20 subsamples x 2 methods on 4 processes, 325 s wall (62–68 s per subsample,
no score runs needed: the cached blocks were reused); peak RSS 1.97 GB per worker.
Analysis 32 s, peak 0.63 GB. Frozen-map rebuild check 150 s,
peak 1.6 GB. Reproduce: `PYTHONHASHSEED=0 .venv/bin/python scripts/39_map_variability.py --stage all --n-proc 4`
(`--stage verify` for the rebuild check). Files: `variability/variability.json`, `*_per_complex_jaccard.tsv`,
`*_pairwise_ari_nmi.tsv`, `*_bootstrap_vs_frozen.tsv`, `*_bootstrap_level1_labels.tsv.gz`, `maps_runtime.json`.
