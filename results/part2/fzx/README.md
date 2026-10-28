# FZ-X freeze (DESIGN v3.1 §0 F4)

Built by `scripts/41_freeze_fzx.py`, which was numbered `39_freeze_fzx.py` at the FZ-X commit and renamed after
merging main's `39_map_variability.py` (review note c). Helpers are in `src/part2/fzx.py` and H1 rules in
`src/part2/h1.py`.

After the base freeze, the script re-applies the registered amendments in order: `42_amend_A1.py`, then
`43_amend_A2.py`. A single run therefore rebuilds the amended freeze. Two fresh-process runs give byte-identical
lists, registry entries and H1 loss table. In the manifest, only `created`, `runtime_s` and `timings_s` differ.
A run takes about 30 s. The freeze manifest is `FZX_manifest.json`. Its `files` map lists the SHA-256 of every
frozen list below, of the gitignored H1 loss table `data_cache/part2/fzx/h1_loss_table.parquet`, and of the FZ-C map
and tier files it was built from. `run_association` checks its inputs against this map.

**Inputs, all verified before use:**
- the FZ-C map files, checked against `results/part2/construction/FZC_manifest.json` (map commit 27fa9b2);
- the context table, LoF matrix and hom-del matrix, checked against `results/part2/contexts/checksums.json`;
- A, checked against the FZ-LB manifest;
- Ensembl 116 paralogs and coordinates, and HGNC.

**What was not read:**
- no dependency value;
- no association statistic;
- no lockbox outcome.

The complete-case part of the registered evaluability rule needs the missingness pattern of A's gene-effect
matrix, so the script reads it as `notna()` and discards the values at once. L appears only as label counts, with
L = U \ A taken from the context-table index.

## Files

| File | Content |
|---|---|
| `zavg_engine_membership.tsv`, `gls_engine_membership.tsv` | The FZ-C maps in the engine format read by `cmap.FrozenMap.from_membership` (unit_id, level, parent, gene, role, locus, tier). Unit IDs are the frozen IDs. Genes are DepMap gene-effect column keys, mapped by Entrez through Part 1 `genes.tsv`; 8 genes have a DepMap symbol that differs from the FZ-C (HGNC) symbol. Role `attached` covers only zavg level-1 attachments that are not locus-flagged (485 pairs; FZ-C DECISIONS D3). |
| `units_zavg_coarse.txt` / `units_gls_coarse.txt` | The tested coarse units (367 / 379). Their `registry.set_hash` is the `complex_set_hash` of H3a, H3b and BENCH (zavg) and of BRIDGE-L (GLS). |
| `h2_submodules.tsv` | The H2 unit set: 109 sub-modules of the 51 coarse complexes that have ≥ 2 sub-modules. |
| `evaluability_tests.tsv.gz` | One row per (entry, context, unit) for H3a, H3b and BRIDGE-L: members after removals, complete-case n / n_context / n_ref, evaluable, reason. m is the number of evaluable rows. |
| `h3a_contexts.tsv`, `h3b_contexts.tsv`, `bridgeL_gls_contexts.tsv` | The final context lists: FZ-LB n_A / ref_A / n_L / ref_L, complete-case n, evaluable tests, and whether the context stays in the family. |
| `lof_homdel_shares.tsv` | For each LoF context, the share of positive A lines that are hom-del calls. All are < 0.5, so no ±10 Mb window removal applies. |
| `control_units.tsv` | For each CTRL entry: the frozen complex unit (or "not evaluable at complex level", with the reason), the readout members, complete-case n, family m and complex-set hash. |
| `h1_events.tsv`, `h1_contexts.tsv`, `units_h1_components.txt` | The H1 event list from the registered enumeration rule (rules 1–8), one context per lost gene, and the counted components. |
| `h1_candidates.tsv` | Every linked gene with at least one loss model in A, its loss-kind counts and its outcome. |
| `h1_lof_check_extra.tsv` | The pre-registered LoF functional check (`contexts.lof_functional_check`) for linked genes that have MSS LoF lines but are not in `results/part2/contexts/lof_functional_check.tsv`: 893 genes, 130 pass. |
| `lockbox_testable.tsv` | Every registered (entry, context) and every H1 event, with n_L, ref_L, lockbox-testable (≥ 5 and ≥ 5) and the confirmation route. |

## Results

### Unit sets

| Map | Coarse complexes | Sub-modules | complex_set_hash |
|---|---|---|---|
| zavg | 367 | 132 | `53f9dbf46370d237` |
| GLS | 379 | – | `14b6c53a01146a61` |

The 109 H2 sub-modules have hash `141a70b7133ed8a4`.

### H3a and H3b

**H3a.** M = 23 families (22 lineages plus melanoma). m = 8,441 = 23 × 367, so every test is evaluable.

**H3b.** M = 20 families and m = 7,340 = 20 × 367; every test is evaluable.
- No primary context falls below 12 A context models after complete-case. The smallest are:

  | Context | A context models after complete-case |
  |---|---|
  | LoF:KMT2D | 12 |
  | fusion:EWSR1::FLI1\|ERG | 15 |
  | LoF:STAG2 | 16 |
  | MYCN_amp | 23 (28 before complete-case) |

- Five listed contexts were already excluded at FZ-R because n_A < 12: LoF CREBBP, KMT2A-r, BCR::ABL1, PAX3/7::FOXO1 and SS18::SSX.
- 18 of the 20 contexts are lockbox-testable. LoF:KMT2D and the EWSR1 fusion (n_L = 4 each) take the rare-context route.

**BRIDGE-L (GLS map).** m = 8,717 for H3a and 7,577 for H3b. Three tests are not evaluable because fewer than 2 members remain after removal.

**H2.** m is set by the committed H3a and H3b discovery runs, as registered (stage-wise).

### Controls

**Evaluable at complex level (15):** PC1, PC2, PC4, PC5, PC6, PC9, PC10, PC11, PC12, PC13, PC15, PC17, PC19a, PC19b and PC20.

**Not evaluable, by the registered rule (5):**

| Control | Reason |
|---|---|
| PC7 | SMARCA2 is unclustered |
| PC8 | STAG1 is unclustered |
| PC18 | KIF18A is unclustered |
| PC21 | CDK6 and CCND3 are in no sub-module |
| PC22 | BRD9 and BICRA are in no sub-module |

This gives the following kill-rule denominators:

| Direction | Evaluable controls |
|---|---|
| H3a | 2 of 2 |
| H3b | 10 of 11 |
| H1 | 3 of 5 (PC6, PC12, PC15) |
| H2 | 0 of 2, so H2 is exploratory under the registered kill rule |

### H1

The registered rule gives **341 events**, on 337 genes, with 366 counted components.

| Event type | Events |
|---|---|
| LoF | 20 |
| CN | 3 |
| silencing | 318 |

Further breakdowns:
- **Link:** 273 events are by membership and 68 by a non-locus attachment.
- **Locus clusters:** 124 events are in locus clusters, of which 123 are silencing-type. They are flagged and carry no claims.
- **Lockbox-testable:** 274 events.
- **Paralog axis (|P| ≥ 1):** 25 events.

This is far more than the "about 10–20" anticipated in DESIGN §2.6 and the registry. Almost all of the excess comes
from rule 1(c), silencing: log1p TPM < 1 where the lineage median is ≥ 3. The 23 LoF and CN events are close to the
anticipated panel. Examples are ARID1A (P = ARID1B), SMARCA4, STAG2, KMT2D (P = KMT2C), EP300 (P = CREBBP), PTEN, RB1,
NF1, STK11, SMAD4, FOCAD hom-del and MAP2K4 hom-del (P = MAP2K7). The rule was applied as registered, without
judgement.

### Amendment A1 (2026-10-28T22:54:50Z; PI decision, before any H1 run)

A1 splits the FZ-X event list into two entries:

| Entry | Events | m | complex_set_hash | Standing |
|---|---|---|---|---|
| `H1_within_complex_GI` | 23 (20 LoF, 3 CN; 4 with a paralog axis) | 27 | `4a5dcc9d02bc4b23` | headline |
| `H1S_silencing_GI` | 318 silencing | 339 | `d6264047856ff8e9` | exploratory, never headlined |

H1S has its own BH family. Its lockbox confirmations are reported, but they never enter the kill rule or the main
claims.

The decision used event counts and types only. The rationale and a confirmation are in
`docs/registry/AMENDMENTS.md`.

`scripts/42_amend_A1.py` (numbered `40_amend_A1.py` at the A1 commit) applies the amendment on top of the base
freeze. It writes:
- `h1_event_split_A1.tsv`;
- `units_h1_components_A1.txt` and `units_h1s_components.txt`;
- a re-assigned `lockbox_testable.tsv`;
- an updated `FZX_manifest.json` (`amendments`, `files`, `summary.m`, `complex_set_hashes`).

The original FZ-X lists are unchanged.

### Amendment A2 (2026-10-28T23:28:31Z; PI ruling on review M1, before any H1 run)

**Rule.** An event is genetic, and goes to H1, iff ≥ 50% of its A loss models carry a LoF or hom-del call. Otherwise
it is silencing and goes to H1S. The rule was applied mechanically to all 341 events.

**Moved:** only `Z0109|SMARCB1`, from H1S to H1. Its genetic share is 0.706 (5 LoF, 7 hom-del, 5 silencing), and it
is CN-majority. It keeps its FZ-X covariate set and components.

| Entry | Events | m | complex_set_hash | Lockbox-testable |
|---|---|---|---|---|
| `H1_within_complex_GI` | 24 (20 LoF-majority, 4 CN-majority) | 28 | `686659438f7d781c` | 22 / 24 |
| `H1S_silencing_GI` | 317 | 338 | `736ed74bc58bc04a` | 252 / 317 |

**Other changes in A2:**
- H1S now has its own direction label `H1S` and an fdr_rule that names this entry only.
- The stale texts are replaced.
- New files: `h1_event_split_A2.tsv`, `units_h1_components_A2.txt` and `units_h1s_components_A2.txt`.
- The A1 files are kept as the record of A1.

The details are in `docs/registry/AMENDMENTS.md`.

### Review M2 and note a: annotation from Ensembl gene IDs

Paralogs are now mapped by Ensembl gene ID, using DepMap Entrez → ENSG from the coordinates file. Chromosomes cover
every Ensembl gene, so RBMX is recognised as an X–Y gametolog of RBMY1*. Renamed DepMap genes such as
TSLIG3A/B/C and TIPE2 keep their paralogs.

Effect on the freeze: four rows of `evaluability_tests.tsv.gz` change. These are LoY and sex_female on `Z0104` (H3b)
and on `G0100` (BRIDGE-L), where RBMX is now removed. n_removed goes from 0 to 1, and n_used from 12 to 11 and from 8
to 7. Covariates are unaffected; they do not depend on the annotation.

Unchanged:
- n and evaluable status;
- every m (H3a 8,441; H3b 7,340; BRIDGE-L 8,717 + 7,577);
- every unit set and hash;
- the control units;
- the H1 event list and the H1 loss table.

### Script 31

`scripts/31_score_blocks.py` differs from the SHA-256 recorded in the FZ-C manifest. The reason is the A-only guard
that FZ-X added (`assert_A_only`, a no-op for A run sets). This is a code-only change and no map is affected; the
independent FZ-C review confirmed a byte-identical rebuild. The note is also recorded in `FZX_manifest.json`
(`notes.script31`).

### Lockbox-testable

| Entry | Testable / total |
|---|---|
| H3a | 23 / 23 |
| H3b | 18 / 20 |
| H2 | 41 / 43 |
| BENCH | 41 / 41 |
| BRIDGE-L | 41 / 43 |
| H1 (after A2) | 22 / 24 events |
| H1S (after A2) | 252 / 317 events |

**Controls:** 20 of 23 are testable. PC15 and PC22 (n_L = 4) and PC16 (n_L = 3) are not.

## Registry fill

Every placeholder in `docs/registry/*.yaml` was filled from these lists, and the `pending` keys were removed.
`registry.fill_only_changes` against the committed FZ-R drafts returns no changes. `tests/test_part2_fzx.py`
re-checks this from git history.
