# Revision analyses R-1 to R-4: intent (POST HOC, DESCRIPTIVE)

Committed before any of the values below is computed. Requested by the round-1 referee (major issues M1–M3,
`manuscript/review/peer_review_round1.md`).

**Standing of every analysis here.** All four analyses are **post hoc and descriptive**. None is registered, none is a
trial, none is a hypothesis test with error control, and none may change a registry entry, discovery, confirmation,
primary result, ledger or `trials.jsonl` line. No `run_association` call is made (it writes the trials ledger); the
discovery machinery is re-implemented read-only from the same engine functions (`assoc._design`, `stats.HC3Design`,
`stats.benjamini_bogomolov`, `assoc.permutation_p`).

**Data.**
- Construction set A only (1,007 models, `assoc.load_A_effects`, manifest-verified). DEMETER2 is read through the same
  A-only loader (`load_A_effects(source="demeter2")`).
- No lockbox L value is read. No persisted read and no L row is opened. R-4 uses only committed L **sample sizes**
  (`n_L`, `n_context_L` from `lockbox_confirm/<ID>/confirmation.tsv`, read by column) and the committed primary call
  labels (`primary_relabelled_R7.tsv.gz`, `primary_call` column only) for the observed rate.
- No hu.MAP, STRING, CORUM, Complex Portal or other curated complex resource.
- Inputs: the frozen zavg map (`results/part2/fzx/zavg_engine_membership.tsv`), the GLS map
  (`gls_engine_membership.tsv`), the FZ-C zavg_BEF A block (S × S, 6,022 genes), the committed discovery tables
  (data_cache, SHA-checked against their committed summaries), `calibration_D8_tests.tsv.gz`,
  `robustness/strongest_findings.tsv` (non-L columns), `robustness/unit_essentiality_A.tsv`, H3a `tests.tsv.gz`.

**Code.** `src/part2/revision.py` (functions, unit-tested in `tests/test_part2_revision.py`) and one script per
analysis: `scripts/58_rev_matched_sets.py` (R-1), `scripts/59_rev_scale_effect.py` (R-2),
`scripts/60_rev_demeter2.py` (R-3), `scripts/61_rev_bridge_power.py` (R-4). Outputs in `results/part2/revision/`
(README.md + TSV/JSON). Every script records wall time and peak RSS. `PYTHONHASHSEED=0`; all random streams are
`numpy.default_rng` with the fixed seeds below.

---

## R-1 Matched-set reference (referee M1)

### R-1(a) Matched sets drawn from genes outside every frozen complex (H3a, H3b)

- **Non-member pool N.** S genes (keys of the FZ-C A block) that are A-effect columns and appear **nowhere** in the
  zavg membership table (any unit, any level, member or attached). Locus-cluster genes are excluded too.
- **Matching (unchanged acceptance rule).** For each of the 367 coarse units, as in `calib.matched_random_sets`:
  size k = the unit's scored members; essentiality (mean clipped A effect) within ±0.1 of the unit; coherence (mean
  pairwise zavg_BEF on A) within ±1.0 z; R = 10 distinct sets; up to 20,000 attempts; if fewer than R are accepted,
  both tolerances are doubled once, else the unit is "unavailable". The only change is the **candidate generator**,
  because the pool has no units to draw from. Three generators are pre-declared:
  - **N-nbr (main).** A seed gene uniform from N; the candidate is the seed plus k − 1 genes drawn uniformly from the
    seed's top-max(2k, 10) neighbours within N by zavg_BEF on A. This is the non-member analogue of "a piece of another
    coherent module". Stream `default_rng([seed, unit index, 1])`.
  - **N-unif.** A uniform k-subset of N, same acceptance (coherence is expected to be rarely met; reported as matching
    quality). Stream `[seed, unit index, 2]`.
  - **N-ess** (referee M1.2(b)). A uniform k-subset of N, matched on size and essentiality only (coherence not
    matched). Stream `[seed, unit index, 3]`.
  Seeds: H3a 4600, H3b 4700 (the entries' calibration seeds).
- **Scoring and testing (unchanged).** `calib.score_sets` with each context's removals (primary score, ≥ 2 members),
  `calib.test_columns` with the context's covariates and registered alternative (≥ 12 / ≥ 12).
- **Matching quality reported** per generator: status counts (ok / tolerances doubled / unavailable), and the median
  and 95th percentile of |Δ essentiality| and |Δ coherence| between each accepted set and its unit.
- **Reported:** pooled and per-context fraction p < 0.05 for (i) real units (committed discovery table, evaluable),
  (ii) complex-member matched sets (committed `calibration_D8_tests.tsv.gz`) and (iii) non-member sets per generator;
  BB rejections over each set of random tests (families = contexts, entry q). Also on the **common-unit subset** (units
  that have non-member sets of that generator), so the three rates compare the same units.

### R-1(b) Selection-aware reference for "beats all matched sets" (H3a, H3b)

- **Random atlases.** The complex-member matched sets are regenerated with `calib.matched_random_sets` and the entry
  seed. Their t must reproduce the committed `calibration_D8_tests.tsv.gz` (|Δt| < 1e-6); otherwise the analysis stops.
  Replicate atlas r (r = 0..9) replaces each unit by its r-th matched set (units without sets are absent).
- **Same discovery selection.** Each replicate atlas goes through the discovery pipeline as in `run_association`:
  analytic HC3 / Satterthwaite tests → BB (families = contexts, entry q) → Freedman–Lane permutation re-test (10,000,
  PatientID blocks within lineage strata, unrestricted for lineage contexts, GPD tail) of every test with BB-adjusted
  q < 0.1 (the discovery runs' screen; their calibration response did not fire, so no raw-p screen) → BB re-run →
  discovery. Permutation streams `default_rng([perm_seed, r + 1, row])` with the entry's permutation seed (H3a 4000,
  H3b 4100). No ledger line is written.
- **Cost rule.** If, after replicate 0 of H3a, the projected wall time for both directions exceeds 2 h, the documented
  approximation is used instead: units passing the same raw-p screen (p < 0.05) are compared by |t|, and the README
  says so.
- **Reported (conditional on discovery):**
  1. **Symmetric flag reference** (referee M1.3). For each discovered random set (context c, unit u, replicate r): does
     its |t| exceed every other member of its group (the other 9 matched sets and the real unit, in c)? The rate among
     random discoveries is set beside the real rate (real discoveries with sets whose |t| exceeds all 10 sets).
  2. **Effect-size distributions.** |t| and |β| of real discoveries against random-atlas discoveries, per context and
     pooled. Headline: the fraction of real discoveries beating random discoveries at matched selection, defined as
     the mean over real discoveries of P(|t_real| > |t_random|) against the random discoveries of the **same context**
     (ties count ½); 0.5 means no difference. The same with |β|, and against random discoveries of the **same unit's**
     sets in the same context where any exist.
  3. The number of random-atlas discoveries per replicate against the real 1,145 (H3a) / 161 (H3b).

## R-2 Lineage scale effect, non-linear trend (referee M2)

- **Trend.** For every H3a context (22 lineages + melanoma), OLS of β on unit essentiality (`unit_essentiality_A.tsv`)
  over **all** evaluable tested units of that context (`tests.tsv.gz`), with a **natural cubic spline basis with
  4 df** (plus intercept): boundary knots at the range, interior knots at the 25/50/75% quantiles of that context's
  unit essentiality (the `ns(x, df = 4)` convention). Residual SD with ddof = 5. The fit is fixed; no LOESS.
- **Per hit (all 207 H3a strongest findings):** nonlinear fitted β, residual β, residual z (residual / residual SD),
  concordant_nl = sign(β) = sign(fitted β). The linear z (D16) is shown beside it.
- **Uniform artefact rule (any context):** `scale_artefact_possible` = concordant_nl AND |nonlinear residual z| < 2.
- **Corrected lead set:** H3a, tier T1, not a control readout (`control_recovery` does not start with "readout"),
  |nonlinear residual z| ≥ 2, not `scale_artefact_possible`. (The last condition is implied by the third; it is kept
  so the rule reads as requested.) Leads are listed by |nonlinear residual z| (descending) with their registered rank.
- **Reported:** per context, the number of strongest findings, |nonlinear residual z| ≥ 2, `scale_artefact_possible`,
  corrected leads, D16 leads (linear), and the transitions of the D16 labels (37 leads, 56 possible artefacts).
- **Never-expressed reliability column.** Per context, the p < 0.05 rate of the 100 dispersed never-expressed sets
  (committed `calibration_D8_tests.tsv.gz`, H3a; expected 0.05) next to the corrected lead count. Descriptive; it
  discounts nothing mechanically.

## R-3 DEMETER2 orthogonal check (referee M2(e))

- **Hits.** Every corrected lead from R-2 (headline); the other H3a strongest findings are reported in the same table
  as a supplement.
- **Data.** DEMETER2 combined gene scores on A models that have RNAi data (A-only loader).
- **Score.** The same complex score as the discovery and the rare-context route: `cscore.unit_scores` on the DEMETER2
  matrix with the context's removals (clip ±3, primary score, ≥ 2 members).
- **Test.** The same HC3 / Satterthwaite design (`assoc._design`) with the context's registered covariates as
  available on those models (covariates constant or collinear on the subset are dropped by the design rules). Fitted
  when ≥ 3 context and ≥ 3 reference models remain; hits with fewer than 12 context models are flagged `low_n`.
- **Reported per hit:** n, n_context, β_D2, t, same sign as β_A, one-sided p in the direction of β_A (p_dir) and
  two-sided p. **Rates:** same-sign rate and p_dir < 0.05 rate among leads (and among all 207). Descriptive; no FDR.

## R-4 BRIDGE-L power comparison (referee M3.3)

- **Hits.** As the R7 primary-evaluable BRIDGE-L (−0.065 [−0.107, −0.025]): eligible, lockbox-testable H3a + H3b
  discoveries of the zavg map, and the BRIDGE-L GLS-map discoveries, excluding NOT EVALUABLE calls. Observed
  confirmed = `primary_call == "confirmed"`.
- **Power (the lockbox results review's power aid, reused).** Per hit, the power of the one-sided α = 0.05 L test:
  - effect: empirical-Bayes posterior from Tweedie's formula on the z of **all evaluable A tests of the same map and
    family** (zavg H3a / H3b discovery tables; GLS `BRIDGE-L_confirmation_rate/discovery_H3a|H3b.parquet`), log-density
    by Efron's Lindsey method (degree-7 Poisson GLM, 120 bins) — the same `tweedie` function;
  - versions: (c) EB posterior mean, (d) EB posterior mean + variance (**headline**), and (a) naive for reference;
  - L SE: **predicted from A only**, A SE × sqrt((1/n_ctx,L + 1/n_ref,L) / (1/n_ctx,A + 1/n_ref,A)) with committed L
    counts (the observed L SE is **not** used);
  - L df: A Satterthwaite df × (n_ctx,L − 1)/(n_ctx,A − 1), clipped to ≥ 2 (no L df is read). Sensitivity: df = n_ctx,L − 1.
- **Reported per map:** mean expected rate (unweighted) and expected rate under the bridge's size-bin
  standardisation (`bench.SIZE_BINS`, pooled bin weights, as `bench.bridge_rate`), observed rates (same weighting),
  the expected difference zavg − GLS with its Poisson-binomial SD, against the observed −0.065; plus median A |z|,
  EB shrinkage and n_ctx,L per map. **Question answered:** does the power-expected difference account for the observed
  BRIDGE-L difference?
