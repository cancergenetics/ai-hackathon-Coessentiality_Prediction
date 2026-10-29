# Intent: referee round 2, post hoc computations A1–A3 (committed before any value is computed)

Referee report: `manuscript/review/peer_review_round2.md` (minor revision). Code: `scripts/62_rev_round2.py`.
Outputs: `results/part2/revision/r5_*`. Summary added to `results/part2/revision/README.md` ("R-5").

**Standing.** All three computations are POST HOC and DESCRIPTIVE. They read committed result tables only: no
dataset, no A or L dependency matrix, no lockbox value, no persisted read, no test label, no download. They change no
registered result, label, call, rank or ledger, and they are not trials. No new statistic is fitted: each is a filter,
a count or a re-sort of committed columns.

## A1. H1 neutrality at a ±0.05 margin (referee N1 / round 1 M6)

- Input: `results/part2/directions/H1_within_complex_GI/classes.tsv` (committed with the H1 run).
- Rows: the events whose final complex-wide label is NEUTRAL (`label_final == "NEUTRAL"`; expected 11, all with
  `N_tost == True` at the registered ±0.1).
- Rule: an event "stays neutral at ±0.05" if its committed 90% CI lies strictly inside the margin:
  `kappa_N_ci90_lo > -0.05` and `kappa_N_ci90_hi < 0.05` (the same two-one-sided-tests reading as the registered
  ±0.1 rule, with the margin halved). The count with non-strict inequalities is reported too, to show whether any
  event sits on the boundary.
- Also reported per event: the smallest symmetric margin that would contain its 90% CI, max(|lo|, |hi|).
- The referee reports 0 of 11; this is a verification. The registered label stays NEUTRAL at ±0.1 whatever the
  count.

## A2. Strongest-findings counts under the alternative registry readings (referee round 1 M7.2; round 2 Section 3)

The two readings the PI accepted after the robustness run (`results/part2/robustness/README.md`) were:
(i) the H1 leave-one-out (LOSO) unit is the lineage; (ii) "A statistic × tier" means tier class first, then |t_A|.
Logic: `src/part2/strongest.py` (`apply_rule`, `TIER_CLASS`, `LOSO_UNIT`) as used by
`scripts/56_strongest_findings.py`, applied to committed tables (`candidates.tsv`, `robustness_per_hit.tsv`,
`strongest_findings.tsv`). Nothing is refitted.

(i) **H1 LOSO unit.** Alternative reading: leave-one-OncotreePrimaryDisease-out, the unit H3a uses. (The other
reading once used in early reviews, leave-one-subunit-out, was never adopted by any registry entry and is undefined for
a single-paralogue P axis; it is noted, not computed.)
- The LOSO unit for H1 enters only the H1 rows. H3a (207) and H3b (45) are unaffected by construction; the script
  checks this from the committed per-hit table (their `loso_unit` column).
- Among the 5 H1 eligible calls, only the candidates (all non-LOSO criteria met) can change status. The script counts
  them from `candidates.tsv` (expected 2: ARID1A P and EP300 P) and records why the other three are not candidates.
- The per-fold H1 fits under the disease unit are not in committed outputs. Computing them would need an A-side refit,
  which this round does not make. The result is therefore reported as a bound: total under the alternative reading
  in [254 − n_H1_candidates, 254], with the H3a and H3b counts unchanged; whether either H1 call changes is stated as
  not determined. The committed lineage-unit folds (folds tested, passed, skipped) of the two calls are reported.

(ii) **Ranking.** The 254 findings are fixed by the rule; ranking cannot add or remove any. Two alternative readings
of "A statistic × tier" are applied to the committed `strongest_findings.tsv`:
- R-alt1: |t_A| alone (tier ignored);
- R-alt2: |t_A| × tier weight, with weights T1 = 3, T2 = 2, T3 (T3-neg / T3-nc) = 1. These weights are not
  registered; they are one illustrative product reading, chosen here before computing.
Ties are broken as in `apply_rule` (context, unit). Reported for each: the number of rows whose rank changes, the
Spearman correlation with the registered rank, the overlap of the top 10 and top 30, and the new ranks of the two
H1 calls and of the text's highlighted findings (Z0035 Ovary, Z0303 Bone, Z0051 Head and Neck, Z0254 Myeloid,
Z0069 RB1 loss).

## A3. Range of the R-1(b) "beats all" rate across the 10 random atlases (referee N5)

- Input: `results/part2/revision/r1b_random_discoveries.tsv.gz` (committed).
- Per direction and atlas (`rep`), over random-atlas discoveries (`discovery == True`): the fraction with
  `exceeds_all_others == True`. Reported: per-atlas values, min, max and the pooled rate.
- Check: the pooled rate must reproduce `r1_summary.json` (`b.summary.random_symmetric_flag` /
  `n_random_discoveries_total`: H3a 3,649 / 11,087; H3b 667 / 1,381), and the per-atlas discovery counts must match
  `b.replicates[*].n_discovery`. The script stops if either check fails.
- This is a spread across atlases, not a confidence interval for the real rate (0.464 H3a; 0.669 H3b).

## Environment

`PYTHONHASHSEED=0`; no random number is drawn. Run: `.venv/bin/python scripts/62_rev_round2.py`.
