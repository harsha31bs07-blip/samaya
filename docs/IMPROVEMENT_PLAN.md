# Improvement plan: closing the gap to HiGHS and SCIP (from 2026-09-27)

This plan turns the research of 26–27 September into ordered work items. That research compared
HiGHS's and SCIP's source with ours, analysed the structure of the 62 MIPLIB instances and ran
HiGHS with its presolve rules switched off one at a time. The rules in `CLAUDE.md` and
`docs/handover/HANDOVER.md` hold for every item: a differential test, a planted bug, the three
presets, a bound A/B on the fixed set, and no tuning to the test set.

## 1. Where we lose (measured, laptop, 600 s, equal real time)

samaya solves 17 of 62, HiGHS 32 and SCIP 25. HiGHS or SCIP solve these 18 and samaya doesn't
(solve time in seconds; TL = time limit):

| Instance | HiGHS | SCIP | Likely missing piece (hypothesis, to confirm in step 0) |
|---|---|---|---|
| neos-3381206-awhea | 6 | 2 | presolve (1900 integer columns with implied bounds) |
| neos-911970 | 9 | 10 | presolve and propagation on continuous columns |
| n5-3 | 34 | 41 | presolve (HiGHS's basic rules, not yet attributed) |
| mc11 | 35 | 14 | flow cuts, presolve aggregator |
| beasleyC3 | 39 | 26 | presolve aggregator (2500 vs 1704 columns), flow cuts |
| timtab1 | 63 | 68 | presolve, cuts on general integers |
| neos-3627168-kasai | 100 | 220 | unknown |
| rococoC10-001000 | 172 | TL | first solution late (heuristics, node selection) |
| csched008 | 185 | TL | clique table in propagation and probing |
| gmu-35-40 | 200 | TL | clique table |
| neos-4738912-atrato | 220 | TL | unknown |
| pg5_34 | 229 | TL | cuts (bound stalls), LP speed |
| tr12-30 | 236 | TL | lot-sizing cuts, presolve |
| gen-ip054 | 596 | TL | LP speed (pure integer, tree size) |
| fastxgemm-n2r6s0t2 | 585 | TL | symmetry, big-M structure |
| ran14x18-disj-8 | 636 | TL | disjunctive big-M, cuts |
| graph20-20-1rand | TL | 8 | symmetry, cliques |
| mcsched | TL | 315 | clique table |

The first six are solved by both in about a minute. They are the best targets: the gap there is
a missing technique, not the last few percent of engineering.

## 2. Step 0: attribute before building (half a day, mostly machine time)

Measured on presolve only so far. Before each large item, switch the matching component off in
HiGHS and SCIP on the 18 instances above, at 600 s, and record which instances stop solving.

- **HiGHS:** `presolve=off`; `mip_detect_symmetry=false`; `mip_heuristic_effort=0`; the
  `mip_heuristic_run_*` switches for single heuristics.
- **SCIP:** `presolving/emphasis off`, `separating/emphasis off`, `heuristics/emphasis off`,
  `misc/usesymmetry=0`, `conflict/enable=false`.
- **How to run it:** a small script in `bench/` next to `harness.py`, 4 at a time, never
  alongside a timing benchmark.
- **What we get:** a table instance × component → solved/unsolved, so each item below starts from
  the instances it should move, stated in advance. An item whose instances don't move in the A/B
  isn't merged, however good the idea is.

## 3. Work items, in order

Effort is in working hours, including tests, presets (about 1.5 h of machine time) and a 30 min
A/B.

### A. First solutions (in progress; 9 feasible instances have none at 60 s)

1. **Best-estimate node selection until the first incumbent: tried, not merged.** Best
   estimate from pseudocosts, with a best-bound pick every 10th node (SCIP's estimate selector,
   HiGHS's every-10th rule), on an indexed queue with lazy deletion.
   - **A/B against f0cbe85, 60 s:** solved 12 → 12, solutions 54 → 54, bound better on 7
     instances and worse on 14, mean distance 19.96% → 20.00%.
   - **Largest gain:** rococoC10-001000, 19.8% → 16.6%.
   - **Largest losses:** enlight_hard, mc11 and b1c1s1.
   - The code is on the local branch `claude/estimate-wip` (commit 5d2e00b).
2. **Randomized rounding** (HiGHS `randomizedRounding`) at every node until the first incumbent:
   **tried, not merged.** It rounds by locks, else floor(x + U(0.1, 0.9)), propagates after
   every fixing and solves the LP over the continuous columns.
   - **Tests:** the reference tests catch a planted bug that leaks its fixings into the node.
   - **A/B against f0cbe85, 60 s:** solved 12 → 12, solutions 54 → 52 (kakapo and rococoC10
     lost), bound better on 8 and worse on 12.
   - **Why it fails:** at every node it costs too much; node counts fall, e.g. csched007
     14,645 → 1,032.
   - Local branch `claude/rr-wip`. A root-only variant is untested.
3. **RENS in the tree without an incumbent: tried, not merged.** RENS reruns on the node LP,
   on RINS's schedule, while there is no incumbent.
   - **Tests:** a test that it runs, with a planted bug (trigger inverted: 0 runs) caught.
   - **A/B, 60 s, on top of the conflict skip:** solved 12 → 12, solutions 53 → 54 (reblock115,
     which the same base build had solved in the previous A/B, so noise), bound better on 10 and
     worse on 7, mean 17.72% → 17.76%.
   - **Why it fails:** its targets don't move (rococoC10 19.70% → 19.78%, csched007
     15.59% → 15.62%).
   - Local branch `claude/tree-rens2-wip`.
4. **Later:** shifting heuristic and zi-round. Central rounding needs an interior-point
   solution, which we don't have, so it is out.

### B. Propagation on continuous columns (6 h)

`BranchAndBound::propagate` (`src/mip/branch_and_bound.cpp:297`) tightens only integer columns.
The research counted many bounds on continuous columns that the rows imply but that are never
derived: ic97_potential 2,674, binkar10_1 49,000, b1c1s1 14,800, cost266-UUE 1.36 million.
Those bounds feed back into the integer columns and into the LP.

- **Change:** tighten continuous columns too, but only when the change is at least 5% of the
  domain, or of max(|bound|, 1) when the domain is infinite (SCIP's `boundstreps`, 0.05). This
  stops endless tiny steps. Keep the row budget.
- **Tolerances:** a new continuous bound is relaxed by the feasibility tolerance times its scale.
  Never set it tighter than the row implies; the LP must stay a relaxation.
- **Conflict analysis** reads only global row bounds, so it is unaffected. Check that
  `effective_bound` and reduced-cost fixing still see integer columns only.
- **Test:** `mip_propagation_never_excludes_an_optimal_solution` on the random families plus a
  continuous-heavy family (fixed-charge network). Every node's bounds must contain the reference
  optimum whenever the node's ancestors do.
- **Planted bug:** drop the relaxation, or use the wrong side for a < 0; the test must fail.
- **Targets:** neos-911970, beasleyC3, mc11, tr12-30, ic97_potential.

### C. Presolve aggregator (the largest measured gap; 2–3 days)

HiGHS keeps far fewer columns after presolve:

| Instance | samaya | HiGHS |
|---|---|---|
| app1-1 | 2477 | 1059 |
| milo-v12-6-r2-40-1 | 2688 | 1698 |
| b1c1s1 | 3872 | 2718 |
| beasleyC3 | 2500 | 1704 |
| nu25-pr12 | 5834 | 4263 |

Switching rules off in HiGHS puts most of the difference on substitution through equations.
samaya substitutes only free singleton columns (`src/presolve/presolve.cpp`). Build it in three
parts, each committed and measured on its own:

1. **Implied-free columns.** A column whose bounds are implied by the other rows
   (activity-bound argument) can be treated as free in one of its equations. Detection only,
   feeding the existing free-singleton substitution. That substitution gets more columns at once.
2. **Doubleton equations** `a x + b y = c`. Substitute x = (c − b y) / a into x's other rows,
   when:
   - x is continuous; or
   - both are integer and a divides b and c exactly (then y's integrality keeps x integral).

   y takes the bounds implied by x's bounds, rounded when y is integer.
   - **Postsolve:** x from y. For the LP, the row dual is restored from x's reduced cost, which
     must be zero if x was strictly between its bounds. Otherwise the dual is computed so that
     x's reduced cost has the right sign.
   - **Testing** follows the existing presolve tests: random LPs through presolve, postsolve and
     `verify_lp_optimality`, and random MILPs against `reference_milp`.
3. **General substitution with a fill limit.** Substitute an implied-free column out of an
   equation with up to about 10 nonzeros when the fill-in is at most about 10. That is HiGHS's
   default limit; state it as taken from HiGHS, not tuned.
   - **Numerics:** reject pivots below 1e-3 of the row's largest coefficient, since a small
     pivot amplifies errors in postsolve.

- **Measure:** columns and rows after presolve on all 62 instances against HiGHS (the script
  from the research, moved into `bench/presolve_compare.py`), then the bound A/B.
- **Planted bugs:** a wrong sign in the postsolve of x, and a skipped divisibility check. The
  MILP differential test must report a wrong optimum or a verifier failure.
- **Risk:** postsolve bugs give wrong answers. The verifier catches them, but only after the
  solve, so the differential tests must be thorough before any benchmark.

### D. Clique table in propagation and probing (1.5 days)

The conflict graph (`src/mip/cuts.cpp`) now feeds only the clique separator. The research found
large cliques on csched007/008 (64/58), gmu-35-40/50 (34/64), graph20-20-1rand (38), mad (20)
and cvs16r128-89 (16). Probing also derives 8,000–40,000 implications on csched, gmu, mcsched,
cvs16r128 and fhnw-binpack.

1. **Keep a clique table:**
   - sources: set-packing rows, rows whose two largest coefficients exceed the slack, and
     probing implications x = 1 ⇒ y = 0;
   - merge cliques that are subsets of larger ones;
   - store each clique once, with a literal → cliques index.
2. **Propagate:** a literal fixed to 1 fixes every other literal of its cliques to 0. Hook it
   into `propagate` next to the rows.
3. **Probing:** it already applies bounds implied by both sides (`src/mip/probing.cpp`), but
   throws away the one-sided implications. Keep x = 1 ⇒ y = 0 as clique edges and
   x = 1 ⇒ y ≤ u as implied bounds. Probing is capped at 400 columns, so raise the cap for
   columns in large cliques only if the time spent stays under its present share.
4. **Implied-bound cuts** from the implications x = 1 ⇒ y ≤ u (x binary, y continuous):
   y ≤ u + (U − u)(1 − x).

- **Test:** propagation never excludes the reference optimum on set-packing and
  assignment-heavy families.
- **Planted bug:** fix the literal itself to 0.
- **Targets:** csched008, gmu-35-40, mcsched, graph20-20-1rand.

### E. Branching score with inferences (4 h)

SCIP's default hybrid score adds inference history (how many bounds a branching on j propagated)
to the pseudocost product, with a small weight. It helps where pseudocosts are uninformative,
like pure feasibility models and early in the tree. samaya doesn't count inferences yet. Count the bound changes that `propagate` makes after each
branching, keep their average per column and direction, and add it with SCIP's default weight
(1e-4 × the average scaled score). No new test is needed beyond the reference tests. It is a change to the
search order, and the A/B decides.

### F. Zero-half cuts (1 day)

Zero-half cuts are {0, 1/2}-Chvátal–Gomory cuts: a mod-2 aggregation of integer rows, following
Koster, Zymolka and Kutschka (2009), with Gaussian elimination mod 2 on the rows with small
slack.
- **Targets:** enlight_hard (parity structure), timtab1, and the general-integer models.
- **Test:** cut validity against the reference optimum on random integer families, as for the
  other separators.
- **Planted bug:** a wrong rounding of the right-hand side.

### G. Symmetry (3–4 days; only after A to F)

- **Detection:** build the column–row colored graph, refine colors to a stable partition, and
  find generators by individualisation–refinement, written from scratch (no nauty or bliss).
  Restrict it to permutations of binary columns that map rows to rows.
- **Use:** orbital fixing in the tree. This means fixing a column to 0 when a column of its orbit
  that is lexicographically first under the stabilizer is 0; Ostrowski, Linderoth, Rossi and
  Smriglio (2011).
- **Targets:** fhnw-binpack4-4, graph20-20-1rand, fastxgemm-n2r6s0t2.
- **Risk:** the largest and easiest to get wrong. The test is the MILP differential on families
  generated with planted symmetry (identical machines, bins).

### H. Not planned now

- **Dual-simplex speed:** it matters on gen-ip054 and pg5_34, but needs a profile first.
- **Cut pool ageing:** small.
- **MRPL structure cuts:** the cases already solve in seconds.

## 4. Order and gates

| Order | Item | Effort | Gate to merge |
|---|---|---|---|
| done | A1 estimate selection | not merged (better 7 / worse 14) | solutions ≥ 54, solved ≥ 12, better ≥ worse |
| done | A2 randomized rounding | not merged (solutions 54 → 52) | solutions up or first-solution time down, no bound loss |
| 1 | Step 0 attribution runs | 0.5 day (machine) | the table in section 2 |
| 2 | B continuous propagation | 6 h | bound A/B better ≥ worse, node rate within 10% |
| 3 | C1–C3 presolve aggregator | 2–3 days | fewer columns, A/B, all presolve tests and differential tests |
| 4 | D clique table | 1.5 days | targets move, no wrong answer |
| done | A3 RENS in the tree | not merged (neutral, targets unmoved) | solutions up |
| 6 | E inference branching | 4 h | A/B |
| 7 | F zero-half cuts | 1 day | A/B, enlight_hard |
| 8 | G symmetry | 3–4 days | targets move |

- **Results:** after C and again after D, rerun MIPLIB at 60 s and 600 s (`bench/compare.sh`),
  with the clock checked, and write them into `docs/results/comparison.md`.
- **Realistic aim:** at 600 s, the six instances both HiGHS and SCIP solve within a minute, plus
  csched008 and rococoC10. That would take samaya from 17 to about 23–25 of 62, level with SCIP.
  Matching HiGHS's 32 needs G and LP speed as well.
