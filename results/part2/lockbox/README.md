# FZ-LB lockbox (DESIGN v3.1 §2.1)

Built by `scripts/33_make_lockbox.py`, using `src/part2/lockbox.py`, in 0.6 s.
- **Inputs:** model labels only. These are the context table (its SHA-256 is checked against `results/part2/contexts/checksums.json`) and Model.csv PatientID. No outcome is read.
- **Determinism:** `--verify` re-draws the lists in a fresh process and compares them byte for byte with the committed lists and manifest hashes. It is run by `tests/test_part2_lockbox.py` with a different PYTHONHASHSEED.

| File | Content |
|---|---|
| `lockbox_models.txt` | L: 338 models (25.1% of U), sorted ModelIDs. |
| `A_models.txt` | A = U \ L: 1,007 models, sorted. |
| `strata.tsv` | n_U, n_A, n_L and the target n_L for every stratification label (lineages and §2.1 genotypes). |
| `context_counts.tsv` | Positives and reference counts in U/A/L for every binary context, the §2.2 size class in A, and whether n_L ≥ 5. |
| `manifest.json` | SHA-256 of the files above, seed, algorithm, input hashes and a balance summary. |

## Algorithm

The draw uses iterative stratification (Sechidis et al. 2011) on PatientID blocks, with fold ratios L/A = 0.25/0.75 and `numpy.random.default_rng(3000)`.

1. **Labels.** One label per OncotreeLineage, plus the 18 genotype labels of §2.1: MSI-H; LoF of ARID1A, SMARCA4, STAG2, SMARCB1, PTEN, RB1, KMT2D, CREBBP and EP300; NRF2-active; KMT2A-r; 9p21/MTAP hom-del; VPS4B loss; and KRAS, NRAS, BRAF and PIK3CA hotspots. Labels use the primary context definitions; NaN counts as label-absent.
2. **Targets.** Each label's L target is round-half-up(0.25 n_U). The exception is a genotype with 17 ≤ n_U < 20, whose target is 5. This rule implements "≥5 L models where n_U permits": n_U ≥ 17 is the smallest size that also leaves ≥ 12 models in A.
3. **Assignment.** The draw repeatedly takes the rarest label still carried by unassigned blocks. It visits those blocks in seeded order and puts each one into the fold with the largest remaining demand for that label. Ties go to the larger remaining total demand, then to a seeded coin.

## Balance

- **Genotypes:** every genotype label hits its target exactly (maximum |n_L − target| = 0).
- **Lineages:** maximum |n_L − target| = 1 (Soft Tissue: 12 vs 13).
- **Patients:** no patient straddles A and L.
- **Genotypes with n_L < 5 (n_U ≤ 16, so not lockbox-testable by size):**
  - LoF SMARCB1 (16 → 12/4)
  - LoF KMT2D (16 → 12/4)
  - LoF CREBBP (12 → 9/3)
  - LoF EP300 (14 → 10/4)
  - KMT2A-r (13 → 10/3)
- **Fusions outside the §2.1 list.** These were not stratified and inherit their split:
  - SS18::SSX: 5/0 (A/L)
  - PAX3|7::FOXO1: 7/0
  - BCR::ABL1: 11/1
  - EWSR1::FLI1|ERG: 15/4

  So PC15/PC22 (SMARCB1 LoF or SS18::SSX) have n_L = 4, and no registered fusion is lockbox-testable.
