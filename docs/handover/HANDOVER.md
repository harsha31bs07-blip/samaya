# Handover: MILP work (from the cloud session, 2026-09-26)

Read `CLAUDE.md` first; its rules hold without exception. This note says where the MILP work
stands and how to continue it.

## Working rules agreed with the project owner

- **Stay on the fundamentals; never tune to the test set.**
  - Fixed instance sets, chosen before any run, and every result reported, bad ones included.
  - No growing or cherry-picking instances until the others look slow.
  - Parameters are chosen from the literature or from principle, and any tuning is stated.
- **Never weaken tests or verifier tolerances.** A new feature needs a differential test
  against `tests/reference_milp.hpp` and a planted bug that the test catches; say both in the
  commit message.
- **Before every commit:** the debug, release and asan presets, zero warnings, with the exact
  test output quoted.
- **Benchmarks:**
  - Never edit a script while a benchmark is running it.
  - Never touch `results/*` branches; the owner pushes laptop results there.
  - Don't build or test on a machine while it runs a timing benchmark.
- **No PR or merge without asking the owner.**
- **Kill any wait loops you leave behind.**

Branch: `claude/kind-hamilton-jr2b6u`. The goal is to be better than HiGHS, the leading
open-source solver, by honest means.

## State on 27 September 2026, evening (IST)

- **Commits since the handover note below:**
  - 9343abf: continuous-column propagation (MILP).
  - 457101d: interior point for convex QP and `--lp-method barrier`, on our own sparse LDLᵀ.
  - b5dab94: conflict-pool skip.
  - 670bab2: the teammate's WP1/WP2 (PDLP, GPU), fixed and integrated.
  - f600682: MIQP.
  - 8297b1e: a time limit in barrier or PDLP is final; fallbacks get only the remaining time.
  - 687448b: PDLP restarts on the KKT error with a primal weight; `bench/gpu_lp.test`.
- **Plans:**
  - `docs/SIH_48H_PLAN.md`: the SIH submission work, with progress.
  - `docs/SIH_GAP_PLAN.md`: the requirement gaps.
  - `docs/IMPROVEMENT_PLAN.md`: MILP items, with what was tried and dropped.
- **Results:**
  - `docs/results/qp.md`: Maros–Mészáros, samaya 106/138 against HiGHS 104.
  - `docs/results/robustness.md`.
  - `docs/results/large_lp.md`: GPU against CPU on large LPs; filled from the night runs of
    27–28 September.
- **GPU on the laptop:** RTX 3050 6 GB. `sudo apt install nvidia-cuda-toolkit` gives CUDA 12.0
  with its own GCC 12 host wrapper; then `cmake --preset cuda` (sm_86). WSL2 restarts clear
  `/tmp`: keep logs in the home directory.
- **Open items:**
  - A crossover for barrier and PDLP. The degenerate LP-like QPs fail strict verification without
    it.
  - The LISWET family (slow interior-point steps).
  - GPU-side evaluation of PDLP's residuals (the host evaluation limits the GPU on small LPs).
  - The generated refinery LP (467,376 x 1,001,520, 4.1M nonzeros): GPU PDLP on the RTX 3050
    stopped at its 600 s limit on 27 September (17:48 IST, a build of the work-in-progress
    tree), unconverged. The night runs retry it at 1000 s with 687448b.
  - The MILP items of IMPROVEMENT_PLAN (presolve aggregator, clique table).

## What is in the branch

Recent MILP work, each item with tests and planted-bug checks (details in `docs/MILP_PLAN.md`):

- **Feasibility Jump** (`src/mip/feasibility_jump.cpp`): LP-free local search before the root
  cuts. The root still runs the objective pump while FJ's point is the incumbent.
- **Dive fix:** dives stop on columns already within the tolerance of an integral bound.
- **Cuts in the tree:** a pool of root cuts, plus fresh c-MIR/cover cuts every 10th depth
  derived with the root bounds. A node keeps its cuts only if they close 1% of its gap
  (they're removed exactly otherwise), and kept cuts may grow the LP by at most 25% of its rows.
  Bound A/B: 23 instances better, 9 worse. Largest gains neos17, enlight_hard, binkar10_1;
  largest loss neos-911970.
- **MIP start:**
  - `--mip-start FILE` in the CLI and `Params::mip_start` in the API.
  - Re-planning variants: `cases/mrpl.py generate --update SEED`.
  - Measured gain on our cases is small; see `cases/README.md`.
- **Verification write-up:** `docs/results/comparison.md`, last section.
- **Bound A/B tools:** `bench/ab_bound.sh` and `bench/ab_report.py`.
- **Conflict analysis** (`src/mip/conflicts.cpp`; dual proofs, Witzig, Berthold and Heinz 2017):
  - An infeasible node LP gives a proof from its Farkas ray: (A'y)'x must lie within
    [min y's, max y's] over the row bounds. A node LP cut off by the incumbent gives one from
    its duals: d'x <= cutoff - offset - min y's, with d = c - A'y.
  - Both are valid for any multipliers, since they use only the rows and their global bounds.
    Proofs sparser than 15% of the columns + 10 are kept (at most 1000, the least recently used
    replaced) and propagated at the start of every node.
  - Round-off entries of the ray (up to 1e-12 of the largest, in the scaled space) are zeroed
    before mapping, or a one-sided row would make the proof's range infinite.
  - Tests: `mip_conflicts_never_exclude_an_optimal_solution` (98 random models with their rows
    scaled by 1e-3..1e3: 556 proofs, 3746 prunes, 0 violations) and
    `mip_farkas_rays_prove_their_nodes_infeasible` (60 lot-sizing models: 265 infeasible node
    LPs, every ray accepted by `verify_infeasibility` on its node's bounds and kept as a proof).
  - The first test's families never reach the Farkas code (0 infeasible LPs): in
    multi_knapsack propagation prunes those nodes, and in fixed_charge and random equality
    models strong branching proves infeasible children before they get an LP. The lot-sizing
    test therefore runs without strong branching.
  - Planted bugs, all caught: the row scaling of the ray dropped (338 of 362 rays rejected by
    the verifier; such a bug can only lose proofs, never make a wrong one), a Farkas proof put
    on the wrong side (181 violations), a sign flip in the cut-off proof (violations and wrong
    optima in both tests).
  - Strong-branching children proven infeasible or cut off give proofs too. With an integral
    objective, cut-off proofs are stated against the integer below the cutoff (the limit
    `effective_bound` prunes by), and a cut-off proof is kept only if it excludes the LP it came
    from. Bound A/B against conflict analysis alone: neutral (solved 8 -> 9, within the +-1 noise
    of two runs of one binary; bound better on 11, worse on 15; mean distance 22.25% both).
  - Bound A/B, 60 s (both builds at the same time): solved 7 -> 9 (neos17, nu25-pr12); bound
    better on 11 instances, worse on 9; mean distance from the optimum 22.44% -> 22.29%, median
    8.34% -> 7.65%. Largest gains neos17 (solved), timtab1 36.5% -> 34.1%, neos-911970 8.4% ->
    7.1%; largest losses enlight_hard 16.2% -> 18.9%, reblock115 1.7% -> 3.5%.

## Next tasks, in the order agreed

The detailed, ordered plan after the 27 September research (targets, tests, gates) is
`docs/IMPROVEMENT_PLAN.md`; the list below is the earlier summary.

1. **Conflict analysis, remaining limit:** a round-off coefficient on a column with an infinite
   bound loses a Farkas proof (the proof cannot be relaxed validly there). Strong-branching
   proofs are in (A/B neutral); their pool pass costs every node, so watch it on large models.
2. **No first solution on 9 feasible MIPLIB instances** (cost266-UUE, csched007/008,
   ic97_potential, neos-3024952-loue, neos-5107597-kakapo, the two rococo, fhnw-binpack4-48),
   where HiGHS, SCIP and CBC mostly find one in 60 s. Diagnosis on rococoC10-001000: Feasibility
   Jump and the pump fail, RENS is infeasible, and the tree restarts every failed plunge at the
   best bound near the root. Tried (worktree branch `claude/next-wip`, commit "depth first until
   the first incumbent"): depth first with backtracking until the first incumbent. It finds
   solutions on rococoC10/B10, cost266-UUE and kakapo (A/B: solutions 52 -> 54, better 13 /
   worse 10, solved 12 both) but traps enlight_hard in a deep infeasible region (no incumbent in
   60 s; best bound finds the optimum), so it is not merged. The hybrid that alternates after
   each failed plunge between backtracking and restarting at the best bound (commit "alternate
   backtracking and best-bound restarts" on that branch) fixes the trap: enlight_hard finds its
   optimum, rococoC10-001000 finds 16490. Bound A/B against f0cbe85: solved 12 both, better 14 /
   worse 8, mean distance 20.00% -> 19.86%, geomean 47.42 -> 47.29 s; the reference tests catch
   a planted node-dropping bug. Not merged yet: its best-bound pick scans the whole stack
   (O(open nodes) per pick; ic97_potential lost 15% of its node rate at 60 s). Before merging,
   index the stack by bound with lazy deletion (a multimap of live positions, tombstones popped
   from the back), then rerun the A/B and the presets.
3. **Presolve gap (measured):** HiGHS's presolve removes far more on about a dozen instances
   (columns after presolve, samaya vs HiGHS: app1-1 2477 vs 1059, milo-v12-6-r2-40-1 2688 vs
   1698, b1c1s1 3872 vs 2718, beasleyC3 2500 vs 1704, nu25-pr12 5834 vs 4263). Turning HiGHS's
   rules off one at a time attributes most of it to the aggregator (substituting columns out
   through equations with limited fill; samaya substitutes only singleton columns). On
   beasleyC3, nu25-pr12 and n5-3 the reductions come from HiGHS's basic rules, not attributed
   yet.
4. **Clique table beyond cuts:** the conflict graph now feeds only the clique separator; HiGHS
   also uses it in propagation and probing. Implied-bound cuts are still missing.
5. **Symmetry** (fhnw-binpack4-4, graph20-20-1rand): orbital fixing on a simple detector.
   Only if time allows.
6. **Structure cuts for the MRPL models**, e.g. (l,S) inequalities. The cases already solve in
   seconds, so this is low priority.
7. **Final results:**
   - Netlib, the cases and MIPLIB at 60 s from the laptop are in `docs/results/comparison.md`.
     Still to add: MIPLIB at 600 s.
   - The first laptop 600 s run switched samaya builds mid-run (the release binary was rebuilt
     at 12:26 UTC while it ran), so samaya is rerun alone on one build.
   - The laptop's WSL2 clock ran about 7% fast, which gave samaya about 6% less time than the
     other solvers (they stop on the real time). Check the clock before a laptop run (see
     `docs/results/comparison.md`, laptop section) and restart WSL if it drifts.
   - SCIP's "gaplimit" status counts as optimal (`bench/harness.py`). Keep numbers from different
     machines apart.

## How to measure

- **Build and test:** see `CLAUDE.md`.
- **Screening (MIPLIB, 62 instances at 60 s):** `bench/run_miplib.sh 4 60 bench/miplib_small.test
  BIN OUT`. Its CSV has no bound, so it measures only incumbents, which vary by several points
  run to run.
- **Judging cuts, conflicts or search changes:** `bench/ab_bound.sh`, both builds at the same
  time. Report better/worse counts and the mean and median bound distance, and name the
  largest losses too.
- **Comparing with the other solvers:** `bench/compare.sh SET` (header lists the sets); HiGHS
  gets the same thread count.
