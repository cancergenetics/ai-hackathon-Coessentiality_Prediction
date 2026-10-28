# Round 1 review: AI-agent failure modes (role: ai_failure)

Scope: DESIGN.md v1, docs/TEST_PROTOCOL.md (draft), harness (`src/testeval`, `scripts/19–21`), docs/part1_audit.md, RESEARCH_TREE.md.
What I ran: code reading; `git log`/`git ls-files`; `testeval.protocol.require_frozen()` and `git_state()` (git commands only). No test file was opened. No score, metric or context association was computed.

## A. Part 1 test look: ways around the single look

**MUST-1. The draft already passes the freeze gate.**
- Problem: `require_frozen()` checks only that TEST_PROTOCOL.md is tracked and that `src/`, `scripts/` and the protocol are clean. The draft was committed in 7a842f2 with the header "draft for PI freeze". Running the gate today returns no error, so `20_score_test_pairs.sh test` and `21 --mode test` would run against the unreviewed draft.
- Change: the freeze needs an explicit marker that a draft cannot carry. Either an annotated tag `test-protocol-frozen-v1.0` or a `docs/TEST_PROTOCOL.FROZEN` file holding the protocol SHA-256 and the fitted_rules SHA-256. `require_frozen()` must then check all of these:
  - the tag or marker exists, and HEAD descends from it;
  - the protocol's SHA-256 equals the marker;
  - the protocol text contains no "DRAFT";
  - `tests/`, `results/test_eval/fitted_rules.json` and the `.venv` lock (pip freeze hash) are unchanged since the tag.
  Add a test that fails when the protocol says DRAFT.
- Why: without this, any agent that runs the documented command early spends the only look.

**MUST-2. The look ledger can be deleted, and the labels are readable outside the harness.**
- Problem: `looks.jsonl` is untracked until after the look. Deleting it (or `git clean`, or `git checkout`) resets the guard. The labels are plain files, and genome-wide scores for all pairs already exist (`outputs/all_pairs`, Part 1's `genome_wide_top_pairs.tsv.gz`). Any agent can therefore compute test AP with no harness at all. `test_only_entry_points_touch_test_files` greps only the repo. It also skips `scripts/verify/audit_part1.py`, which reads the test files, and it cannot see scratchpad scripts or `python -c` calls. `--mode dryrun` would also evaluate any file copied into `data_cache/test_eval/dryrun/`.
- Change:
  1. **Escrow.** Ask the user to move the two Test files out of the project. The PI sends the user the freeze-tag hash, and the user restores the files only after that. The look then asserts the files' SHA-256 equals the DATA.md values before it reads them.
  2. **Commit the intent first.** Commit the `started` record (with the tag hash) *before* labels are read: a "look intent" commit. The script then refuses to run if any commit in `git log -- results/test_eval/test_results.json results/test_eval/looks.jsonl` already exists. Git history is harder to erase silently than an untracked file.
  3. **Notarise outside the agent's reach.** Report the protocol hash, intent-commit hash and results hash to the user in chat. They then exist somewhere the agents cannot rewrite.
  4. **Close the partial-look gap.** Mark a relook `confirmatory:false` if `test_results.json` exists, not only when a `done` record exists. Otherwise a look killed after the results are dumped but before `done` is written could be re-run as confirmatory.
- Why: the single-look guarantee is currently procedural. These steps make a bypass either impossible or visible.

**MUST-3. Undisclosed prior contact with the test files, and a contaminated operating point.**
- Problem 1: RESEARCH_TREE D1 says an audit agent counted test pairs in Part 1's z ≥ 3 table, which is equivalent to precision and recall at rule (a). The tree says this is "Disclosed in TEST_PROTOCOL". It is not: §0 does not mention it.
- Problem 2: §0 says "not computed: overlap counts, protein lists, mapping rates". But docs/part1_audit.md (commit 7d1fbdd, before the protocol draft) reports:
  - the overlap with train (0), so S2 equals the ITT set and R4 is moot;
  - the C1/C2/C3 counts;
  - the scorable pairs *by label* (pos 3,159/3,312, neg 563,648/597,426). That is the unmapped rate per class, which bears directly on the R5 ITT-versus-mapped decision.
- Problem 3: §0 and §9.5 still call the audit "pending".
- Change:
  - Rewrite §0 to list everything any agent computed on the test files, with commit references.
  - Declare rule (a) (z ≥ 3) **non-blind**: report it, but label it as compromised.
  - State that R4/R5/S2 were drafted after the label-stratified mapping counts were known. The primary choices stay as they are, and the reason is given.
  - Update the audit status.
- Why: a protocol that misstates what has been seen breaks the "explicitly document the transition" requirement. A reviewer who finds the gap would discount the whole look.

**SHOULD-4. The verifier is not independent enough.**
- Problem: the verifier re-derives AP, AUROC and ΔAP from the harness's own score TSVs. It does not recompute the verdict, and nothing stops it from seeing the harness numbers first. That invites anchoring and copying.
- Change: the verifier receives only the frozen protocol, the label files and the score TSVs. It must independently implement R1–R5, the floor rule, the test clusters, the seed-0 cluster bootstrap and the §7 verdict. It reports these before it sees `test_results.json`. Pre-register the tolerances: point metrics 1e-9, CI endpoints exact (same RNG stream). As a separate check, a second agent re-scores a 1% sample of test pairs from the raw CRISPRGeneEffect file through the locked predictor.

**SHOULD-5. Remove the remaining discretion at freeze.** §9 leaves options open ("keep or revert", "recommended"). The frozen text must contain only the chosen options; the alternatives go into a changelog. The reporting language for the post hoc reference is fixed in advance as well. If `spec_k20` beats zavg_BEF on test, the result is reported as "a post hoc reference, not a candidate"; it is never written up as "a better method exists".

## B. Part 2: results-peeking, forking paths, cherry-picking

**MUST-6. §2.6 selects methods on association outcomes.**
- Problem: the pipeline variants (z0, k, γ, BioPlex, consensus, score type) are compared on the "split-half replication rate of context associations". The same models then supply the full-data FDR claims. This is selection on the outcome. It also contradicts the §2.2 freeze ("complexes frozen before §2.3"). With roughly 4×3×|γ|×2×2×3 variants, the replication rate is maximised by search.
- Change:
  - Construction and score type are chosen **only** on the label-free criteria of §2.2.
  - Freeze construction with the hash of the frozen complexes **before any context table is built**.
  - Variant comparison on association replication becomes a reported sensitivity analysis and never a selection rule.
- Why: otherwise the discoveries are tuned on themselves.

**MUST-7. Add a model-level lockbox; the split halves are not held out.**
- Problem: agents will see full-data results (the primary claim) and the 10 split seeds. The "validation half" is therefore seen while hypotheses are iterated. PROMPT's "keep going until you run out of ideas" maximises forking paths.
- Change:
  - Set aside a lineage-stratified lockbox of about 25% of models (seed fixed now), and commit the hash of its model list.
  - Exclude the lockbox from complex construction: recompute zavg_BEF on the other models with the locked configuration. This is about 150 s and changes no Part 1 setting.
  - Exclude the lockbox from all exploration, controls and direction triage.
  - Read the lockbox once per registered hypothesis, through a gated script with a look ledger (reuse the MUST-1/2 machinery).
  - Headline claims require lockbox replication (same sign, pre-set α); DEMETER2 stays an orthogonal check.

**MUST-8. Add an analysis registry and an automatic trial ledger; this is how to stop cherry-picking among directions.**
- Problem: D1–D8 plus "D7/D8 only if time" leaves the PI free to report the directions that worked. Context resets (the tree is the memory) cause analyses to be silently re-run, which hides multiplicity.
- Change:
  1. Create `docs/ANALYSIS_REGISTRY.md` with one entry per analysis ID, committed **before** its first run. Each entry gives:
     - the hypothesis and its direction;
     - the context family and its definition;
     - the complex set (by hash), the test, the covariates and the FDR family;
     - the success criterion and the stopping rule.
  2. Route every association through one wrapper, `run_association(analysis_id, …)`. It refuses an unregistered ID or an uncommitted context table. It appends to `results/part2/trials.jsonl` the git commit, the input hashes, the number of tests and the summary counts (not the top hits).
  3. Fix the cross-direction rule now. Write down the rule for choosing the "strongest findings" for the paper, for example lockbox-replicated plus a minimum effect, ranked by lockbox effect. Correct across directions with BH over all registered families, or with a hierarchical correction (family, then hypothesis).
  4. Report every registered direction, including null and abandoned ones, in a supplementary table with its trial counts.

**MUST-9. Make the positive controls testable and outcome-blind; fix vague or wrong ones.**
- Problem:
  - Several controls are disjunctive or vague, which lets any hit count as a success: "MYB / CBFB–RUNX1 / Mediator kinase **or** MLL–menin"; "HNF1B/HNF4A-type"; "EWSR1::FLI1 → lineage modules".
  - "KEAP1-LoF → NFE2L2-axis loss of dependency" has the sign wrong or ambiguous. The expected effect is *increased* NFE2L2 dependency, and LITERATURE.md itself marks it "[unverified for 26Q3]".
  - The olfactory-receptor negative control cannot fail, because the SD filter and the expression filter remove those genes.
- Change: give each control one row with:
  - the context definition;
  - the target, defined as "the frozen complex containing gene X", or gene-level if X is unclustered (and then "not evaluable at complex level" is a pre-declared outcome);
  - the sign;
  - the statistic and the threshold.
  Evaluate the controls **once**, by a separate script, after the construction and context freezes. Report the hit rate over a fixed denominator. If a control fails, that is reported; any later pipeline change is a logged deviation, with results reported before and after the change. Replace the olfactory-receptor control with pseudo-contexts matched for size and lineage composition and with gene sets matched for coherence (mean within-set zavg). Those can fail.

**SHOULD-10. Freeze the inputs, not just the code.**
- Fix these before any association, as hashed files with a README:
  - the context table, including the CN thresholds after the scale check, with the check's output logged;
  - the definition of "selective dependency" used in the §2.2 coverage criterion;
  - the expression rule.
- Analyses run before the expression data arrive are labelled preliminary and re-run, never mixed with later results.

**SHOULD-11. No flow from test results into Part 2.** F2 bars test *pairs* from Part 2 but not test *results*. Add to F2: Part 2 parameters (z0 and the rest) may not cite test-set metrics. The construction scripts must not read `results/test_eval/`; enforce this with a grep test.

**SHOULD-12. Verify references, biology and reproducibility.**
- An agent that did not write LITERATURE.md checks every manuscript citation against its DOI or PMID, and every biological claim in §2.7 against the cited source. Unverifiable claims are dropped.
- A clean-clone `run_all` regenerates all results; output hashes are compared to the committed ones (`scripts/verify/reproduce.sh`).

**COULD-13. Commit the hash of each post hoc interpretation run** (HGNC gene_group, literature naming) after the construction freeze tag. The script refuses to run before that tag, so post hoc labels cannot steer construction informally.

**COULD-14. Disclose the circularity.** zavg_BEF was selected on hu.MAP train labels, which derive from Complex Portal. Cluster agreement with curated complexes is therefore partly built in. State this wherever post hoc comparisons appear.
