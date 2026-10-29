# Supplementary Table S6: lockbox confirmation rates per direction

| direction | subset | hits (lockbox route) | NOT EVALUABLE (primary) | confirmed (primary) | rate, primary evaluable | rate, primary all-in | confirmed (corrected) | rate, corrected (non-confirmatory) | rate, corrected all-in | defect-sensitive calls | expected (post hoc power aid) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| H3a_lineage_atlas | overall | 1145 | 0 | 601 | 0.525 | 0.525 | 601 | 0.525 | 0.525 | 0 | 0.47–0.52 |
| H3a_lineage_atlas | non_locus | 942 | 0 | 524 | 0.556 | 0.556 | 524 | 0.556 | 0.556 | 0 |  |
| H3a_lineage_atlas | beats_all_matched_sets | 502 | 0 | 287 | 0.572 | 0.572 | 287 | 0.572 | 0.572 | 0 | 0.52–0.53 |
| H3a_lineage_atlas | not_beats_all | 581 | 0 | 284 | 0.489 | 0.489 | 284 | 0.489 | 0.489 | 0 | 0.41–0.43 |
| H3b_driver_family | overall | 150 | 28 | 60 | 0.492 | 0.4 | 69 | 0.46 | 0.46 | 30 | 0.43–0.44 |
| H3b_driver_family | non_locus | 110 | 17 | 42 | 0.452 | 0.382 | 46 | 0.418 | 0.418 | 19 |  |
| H3b_driver_family | beats_all_matched_sets | 94 | 22 | 42 | 0.583 | 0.447 | 48 | 0.511 | 0.511 | 22 | 0.50–0.51 |
| H3b_driver_family | not_beats_all | 50 | 6 | 15 | 0.341 | 0.3 | 18 | 0.36 | 0.36 | 8 | 0.29–0.30 |
| H1_within_complex_GI | overall | 4 | 1 | 2 | 0.667 | 0.5 | 3 | 0.75 | 0.75 | 1 |  |
| H1_within_complex_GI | non_locus | 4 | 1 | 2 | 0.667 | 0.5 | 3 | 0.75 | 0.75 | 1 |  |
| H1S_silencing_GI | overall | 18 | 1 | 9 | 0.529 | 0.5 | 9 | 0.529 | 0.5 | 0 |  |
| H1S_silencing_GI | non_locus | 14 | 1 | 6 | 0.462 | 0.429 | 6 | 0.462 | 0.429 | 0 |  |
| H2_submodule_divergence | overall | 109 | 1 | 76 | 0.704 | 0.697 | 77 | 0.706 | 0.706 | 1 |  |
| H2_submodule_divergence | non_locus | 97 | 1 | 67 | 0.698 | 0.691 | 68 | 0.701 | 0.701 | 1 |  |
| H2_submodule_divergence | beats_all_matched_sets | 40 | 1 | 26 | 0.667 | 0.65 | 27 | 0.675 | 0.675 | 1 |  |
| H2_submodule_divergence | not_beats_all | 69 | 0 | 50 | 0.725 | 0.725 | 50 | 0.725 | 0.725 | 0 |  |
| controls (23) | overall | 44 | 13 | 29 | 0.935 | 0.659 | 39 | 0.929 | 0.886 | 11 |  |
| controls (23) | non_locus | 43 | 12 | 29 | 0.935 | 0.674 | 38 | 0.927 | 0.884 | 10 |  |
| BRIDGE-L_confirmation_rate (GLS hits) | overall | 1335 | 9 | 764 | 0.576 | 0.572 | 769 | 0.576 | 0.576 | 10 |  |

- Primary = the single read as executed (5bc9444) under ruling R7. Fits that failed because of the driver defect (the L covariate design was assembled over all models and then subset to L) are NOT EVALUABLE; 'evaluable' excludes them from the denominator, 'all-in' counts them as not confirmed. Three further fits failed for a genuine reason (CTRL_PC21 ×2, H1S HNF1B).
- Corrected = the disclosed, NON-CONFIRMATORY relook of all 560 affected fits (510 whose design changed, 50 that had failed) through the registered resume path, after the fix and an intent commit (R6, R7). A call that differs between primary and corrected is defect-sensitive and carries no claim.
- H3b is always reported as evaluable / all-in / corrected together (D14): 26 of its 28 NOT EVALUABLE fits are in whole-genome doubling, which is effectively not evaluable in the primary.
- Expected (post hoc, descriptive; not registered, not a test): mean power of the one-sided L test over hits, with the A effect shrunk by an empirical-Bayes (Tweedie) posterior; the range spans the posterior mean / mean-plus-variance variants and the two SE choices (independent review §4). The gap between 'beats all matched sets' and the rest is about what this expectation predicts and is not evidence of complex-specificity.
- For H2 the subset 'beats_all_matched_sets' means: beats all 20 of its own random member splits.
- Rates exclude rare-context hits (H3b 11, H1 1, H1S 6), which are reported separately (Table S8). The 'no sets' subset (H3a 30/62, H3b 3/6) is in results/part2/lockbox_confirm/rates.tsv.
- Sources: results/part2/lockbox_confirm/comparison_primary_vs_corrected_R7.tsv, rates.tsv; results/part2/verify/lockbox_results_review/tables_review.json (power, power_by_beats_all); script scripts/tab_lockbox.py.
