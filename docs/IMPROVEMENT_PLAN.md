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

### Step 0 results (7–8 October 2026, laptop, alone on the machine)

`bench/attribute.py`: each baseline on the instances it solves, with its defaults and with one
component switched off, 600 s, one thread, 4 runs at a time; 118 runs. Seconds, or TL (not
solved in 600 s). Full table: `bench/results/attribution/summary.md` (not committed, like all
of `bench/results/`).

| Component off | HiGHS instances lost (of 13 it solves) | SCIP instances lost (of 9) |
|---|---|---|
| Presolve | **7**: mc11, neos-3627168-kasai, rococoC10-001000, csched008, gmu-35-40, neos-4738912-atrato, tr12-30 (and beasleyC3 38 → 506 s) | 4: mc11, beasleyC3, neos-3627168-kasai, graph20-20-1rand |
| Cuts (SCIP `separating`) | – | **7**: neos-3381206-awhea, neos-911970, n5-3, mc11, beasleyC3, timtab1, neos-3627168-kasai |
| Heuristics | 3: gmu-35-40, neos-4738912-atrato, pg5_34 | 0 (slower on six) |
| Symmetry | 1: neos-4738912-atrato | 2: neos-3627168-kasai, graph20-20-1rand |
| Conflict analysis | – | 0 |

HiGHS's defaults lose gen-ip054, fastxgemm-n2r6s0t2 and ran14x18-disj-8 on this machine (they
were near the limit in the earlier run), so they attribute nothing.

What it says for the items below:
- **C and D first.** HiGHS's presolve is its largest single gain. Its rules are more than the
  aggregator (probing, clique merging, dominated columns), so after C, attribute again with the
  HiGHS presolve rules switched off one by one (`presolve_rule_off` bit mask), as the research
  of 26–27 September started to.
- **Cuts matter more than the plan assumed.** SCIP without separation loses 7 of 9. samaya has
  Gomory, c-MIR, covers, flow covers and cliques, so the next step is to measure which of SCIP's
  separators carry those instances (`separating/<name>/freq = -1` one at a time) before choosing
  between F (zero-half) and stronger aggregation in c-MIR.
- **G (symmetry)** is confirmed for three instances; it stays last.
- **Conflict analysis** carries no instance here; no further work on it.

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

#### Design (13 October 2026, revised after review)

Presolve now works on a mutable matrix and postsolves from column snapshots (270568e), so
substitutions can change rows. The first draft of this design was reviewed by Codex (astra,
high); the review found three cases where it would give wrong answers (approximate
divisibility, integer bound transfer, dropping small coefficients), and this version fixes them.

**Equivalence is proved, not approximated.** Every reduction must leave a model with the same
feasible set and optimum, not a relaxation. Where floating-point arithmetic cannot settle a
case (a ratio that is "nearly" an integer, an integer bound that sits within rounding of an
integer, a residual coefficient that is "nearly" zero), the reduction is **skipped**, never
taken with a tolerance. Each skip is counted in `PresolveStats` so the measurements show what
the caution costs.

**Mode.** MIP mode (`PresolveOptions::mip`: the root, B&B restarts and the sub-MIP heuristics)
restores x only. C1 and C3 are valid in both modes (their dual postsolve below is a one-line
rule); C2 starts in MIP mode and gets its dual postsolve as a separate step. The verifier
re-checks returned solutions, but it cannot see a better solution or a feasible region that a
wrong reduction removed (root infeasibility, a restart concluding optimal, a sub-MIP finding
nothing), so the tests compare status and optimum against the reference solvers, not only the
returned vector.

**Implied-free test (C1, C3).** For column j in row i with coefficient a: the activity bounds of
the row's *other* columns, kept as finite sums plus counts of infinite contributions (never
by subtracting j's term from a total that may be infinite), and the row bounds [rl, ru] give
j's implied range [(rl - U)/a, (ru - L)/a], swapped when a < 0. Any infinite or NaN part, or
an empty range, means "not implied free". The proof needs the implied range inside
[l_j, u_j] *shrunk* by an outward error allowance for the sums and the division (relative 1e-12
of the terms' magnitudes), not "inside up to the feasibility tolerance". It is rechecked right
before each substitution and never stored, since bounds keep changing. Redundant bounds can
still be attained; what the proof gives is that they need no multiplier in the restored duals.

**C1. Implied-free column singletons.** `process_singleton_cols` today requires
l_j = -inf and u_j = +inf. It now also takes a singleton column that is implied free in its
row. Integer columns stay excluded in MIP mode. The existing kSubstituteCol / kSlackCol
postsolve is valid as it is: the restored value satisfies the row, so it lies in the implied
range, and its reduced cost is zero.
- Tests: random-family differential (presolve on/off, LP and MILP; status, optimum and
  verification compared) with bounds made redundant on purpose; one and two infinite
  contributions in the row; a column tightened after the reduction. Assert that the rule fires.
- Planted bug: the implied range computed without the a < 0 swap.

**C2. Doubleton equations** `a x + b y = c` (row with exactly two active entries, equality).
Eliminate x: x = (c - b y) / a.
- **Which one is x:** prefer a continuous column; then the one with fewer other entries; then
  the larger |coefficient|. Refuse when |a| < kSubstitutionPivotRatio (0.01, the singleton
  substitution's value) · max(|a|, |b|); count rejects to decide later whether to relax it.
- **Integer x:** only when y is integer and a, b, c are exactly integral doubles of magnitude
  at most 2^31 with b and c divisible by a in integer arithmetic (b % a == 0, c % a == 0).
  Anything else is skipped (or the roles swapped when y is continuous).
- **Bound transfer:** y's new bounds are its bounds intersected with
  [(c - a u_x)/b, (c - a l_x)/b] (swapped when a/b < 0). For integer y the transferred
  endpoints are not rounded with `set_bounds`' integrality tolerance (an allowance in y is
  multiplied by |b/a| in x). Each candidate integer endpoint Y is tested by computing
  x(Y) = (c - b Y)/a with a rounding-error bound e = 4 · 2^-52 · (|c| + |b Y|) / |a|: Y is
  accepted when x(Y) is inside [l_x, u_x] by more than e, rejected (move one integer inward)
  when it is outside by more than e, and in between the substitution is skipped. Integer endpoints that
  cross prove infeasibility only when the comparison is unambiguous; otherwise skip.
- **Rows:** for every other active row r with coefficient a_rx on x: a_ry += -a_rx · b / a and
  both row bounds shift by -a_rx · c / a. A result that is exactly zero removes the entry; a
  nonzero residual is kept, and if it is numerically unacceptable (below 1e-9 of the row's
  largest coefficient) the substitution is rejected, never dropped: dropping can make a feasible
  model infeasible (εy + z = 1 with ε = 2^-41 and y up to 2^42). Updates are staged, then
  applied, since the helpers mutate the lists being traversed; all results must be finite
  before anything is committed.
- **Objective:** c_y += -c_x · b / a; offset += c_x · c / a.
- **Record** (saved before any mutation): row, x, y, a, b, c; x's cost and column snapshot;
  x's and y's bounds before the transfer; y's cost and column snapshot right after the
  substitution (for the LP dual postsolve). **Primal postsolve:** x = (c - b y) / a.
- **LP dual postsolve (second step):** with q = c_x - Σ_(r≠e) a_rx λ_r over x's snapshot and
  d'_y the reduced cost of y right after the substitution: if d'_y has a valid sign at y's
  bounds before the transfer, λ_e = q / a; otherwise the bound multiplier moves to x:
  λ_e = q / a + d'_y / b, d_y = 0, d_x = -(a/b) d'_y. Example: min y, x - y = 0,
  1 <= x <= 3, 0 <= y <= 10: after elimination y = 1 with d'_y = 1, and the correct duals are
  λ_e = -1, d_x = 1, d_y = 0.
- **Order:** after the singleton and duplicate-row passes in each round.
- Tests: the review's counterexamples as fixed cases (near-divisible ratio, the infeasible
  integer transfer x + 1000y = 999.9995, the ε-cancellation); both signs of a/b; integer pairs
  divisible and not; a survivor tightened, fixed or eliminated later; aggregation before and
  after forcing and duplicate rows; random MILPs and LPs with planted doubletons against
  `reference_milp` / `reference_lp`. For the LP step: survivors at transferred lower and upper
  bounds, coincident old and new bounds, fixed survivors, maximization.
- Planted bugs: the postsolve sign of b y; the divisibility check skipped; the multiplier
  transfer skipped (LP step).

**C3. Implied-free substitution with a fill limit.** Column x implied free in an equality row e;
row_r -= (a_rx / a_ex) · row_e for each other row r of x, with the same pivot test (against the
largest entry of e), exact-zero removal, rejection of unacceptable residuals, staged updates,
bound shifts and objective transfer as C2; x's bounds are dropped (implied). Integer x under the
same exact-integer divisibility rule for every coefficient of e and its right-hand side, and
every other column of e integer. Limits, samaya's own choices for a first version (HiGHS limits
net fill at 10 and treats length-two rows and columns specially; it has no length caps):
e has at most 10 entries, x at most 10 other entries, net nonzero growth at most 10; each limit
counts its rejects. Primal postsolve: x = (b_e - Σ a_ek x_k) / a_ex. **Dual postsolve, both
modes:** λ_e = q / a_ex with q = c_x - Σ_(r≠e) a_rx λ_r over x's snapshot; every other
reduced cost is unchanged (x's reduced cost is zero, which an implied-free column allows even
when it lands on a bound).

**Known related risk.** Duplicate-row detection treats rows as parallel within a relative
1e-12 (`kParallelTol`); y + z = 0 and y + (1 + ε)z = 1 would be merged although they are jointly
feasible. Aggregation can create such nearly parallel rows. Make the parallel test exact
(identical normalized coefficients) before C2 lands, and measure what it costs.

**Measure** after each of C1, C2 and C3, before the next:
- columns and rows after presolve on all 62 instances against HiGHS (the script from the
  research, moved into `bench/presolve_compare.py`), and the skip counters;
- then the bound A/B and MIPLIB at 60 s; no wrong answer anywhere.
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

## 5. Local branches (status on 2026-10-07)

None of these is pushed. They are kept until the owner says otherwise; `results/*` branches are
left alone.

| Branch | Status |
|---|---|
| `claude/activity-wip`, `claude/cont-prop-wip` | Landed: same tree as 9343abf (propagation on continuous columns) |
| `claude/conflicts-sb-wip` | Landed: same tree as 28f1ce5 |
| `claude/cutset-wip` | Landed: same tree as 0bff207 |
| `claude/fallback-time` | Landed: same tree as 8297b1e |
| `claude/integration`, `claude/miqp-wip` | Landed: f600682 (MIQP) |
| `claude/pdlp-restarts` | Landed: same tree as 687448b |
| `claude/qp-next` | Landed: same tree as 457101d |
| `claude/qp-ipm` | Earlier draft of the QP interior point; superseded by the reworked 457101d |
| `claude/wp1-wp2-merge` | Earlier attempt at the WP1/WP2 merge; superseded by 670bab2 |
| `claude/conflict-speed-wip` | Draft of the conflict skip, superseded by b5dab94 |
| `claude/tree-rens2-wip` | A3, RENS in the tree: tried, not merged (section 3A) |
| `claude/estimate-wip` | A1, best-estimate node selection: tried, not merged (section 3A) |
| `claude/rr-wip` | A2, randomized rounding: tried, not merged (section 3A) |
| `claude/next-wip` | The alternating first-solution hybrid: not merged; see `docs/handover/HANDOVER.md` (d20a9ba) |
| `claude/tree-rens-wip` | Fully merged; can be deleted |
