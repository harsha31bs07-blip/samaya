# samaya against the open-source solvers

samaya compared with HiGHS 1.15 (highspy), SCIP 10.0 (PySCIPOpt), CBC 2.10.11 and GLPK 5.0.
Everything was run with `bench/compare.sh` on the same machine, one instance at a time per shard,
every solver single-threaded with the same relative MIP gap (1e-4).
- Unsolved instances count at the time limit in the shifted geometric mean (shift 10 s).
- The test sets were fixed before any solver was run and are reported in full.

Machine: the 4-core cloud container (Intel Xeon, 2.8 GHz), 4 instances at a time. The laptop
runs (Dell G15) are added as they arrive; numbers from the two machines are never mixed.

## Netlib LP (93 instances, 300 s)

| Solver | Solved | Shifted geomean |
|---|---|---|
| HiGHS | 93/93 | 0.17 s |
| CBC (Clp) | 93/93 | 0.19 s |
| **samaya** | **93/93** | **0.32 s** |
| SCIP (SoPlex) | 93/93 | 0.41 s |
| GLPK | 92/93 | 0.76 s |

GLPK's one miss (e226) was a harness bug: an objective constant was lost when converting the
model to LP format for GLPK. It is fixed; after the fix all five agree on e226.

## MRPL case studies and generated MILPs (16 instances, 300 s)

These are `cases/` (planning LP, crude scheduling MILP, utility MILP; three sizes each) and the
seven generated refinery-scheduling, knapsack and facility models.

| Solver | Solved | Shifted geomean |
|---|---|---|
| SCIP | 16/16 | 3.18 s |
| CBC | 16/16 | 3.68 s |
| **samaya** | **16/16** | **3.71 s** |
| HiGHS | 16/16 | 4.63 s |
| GLPK | 12/16 | 19.6 s |

No solver contradicts another on any instance.

## MIPLIB 2017, 62 small "easy" instances, 60 s

Two same-machine runs: before and after the overnight work of 25–26 September.

| Solver | Solved | Shifted geomean |
|---|---|---|
| SCIP | 18/62 | 39.91 s |
| HiGHS | 16/62 | 42.57 s |
| **samaya**, build cda171e | **7/62** | **50.99 s** |
| CBC | 7/62 | 55.26 s |
| samaya, build c9a2a31 (before the night) | 5/62 | 54.57 s |

The HiGHS, SCIP and CBC numbers above are from the second run; in the first they were 16, 18
and 7 with geomeans within 1%.

Solve times of every instance samaya solves (seconds; "–" = not solved in 60 s):

| Instance | samaya | HiGHS | SCIP | CBC |
|---|---|---|---|---|
| markshare_4_0 | **14.4** | – | – | – |
| mas76 | **25.5** | – | – | 35.0 |
| pk1 | **43.0** | – | – | 45.9 |
| neos859080 (infeasible) | **0.1** | 1.8 | 0.8 | – |
| exp-1-500-5-5 | 3.5 | 4.2 | **3.0** | – |
| app1-1 | 13.5 | 25.5 | 8.5 | **6.3** |
| sp150x300d | 0.8 | **0.1** | 0.4 | – |

- samaya solves three instances that neither HiGHS nor SCIP solves in 60 s (markshare_4_0,
  mas76, pk1), and proves neos859080 infeasible faster than any of them.
- HiGHS or SCIP solve 14 instances that samaya does not:

| Cause | Instances | What the logs show |
|---|---|---|
| Big-M network gaps | p200x1188c, mc11, beasleyC3, n5-3 | Root bound still far below the optimum (p200x1188c 6,020 vs 15,078; mc11 1,156 vs 11,689) |
| No solution for general-integer / symmetric models | neos-3381206-awhea, enlight_hard, fhnw-binpack4-4, graph20-20-1rand | No incumbent: the pump ignores interior general integers, no symmetry handling |
| Slow bound progress | binkar10_1, mik-250-20-75-4, nu25-pr12, pg, neos17, neos-911970 | Good solution, bound closing slowly (no cuts in the tree) |

**Fixed overnight (all in docs/MILP_PLAN.md):**
- Aggregated c-MIR with variable bounds: exp-1-500-5-5 in 3.5 s.
- Flow covers: sp150x300d in 0.8 s.
- Reduced-cost fixing and a root restart: mas76 from 53.6 s to 25.5 s.
- A cap on strong-branching effort.
- Three simplex speed-ups: a warm-started solve factorized three times, now once.

After this run, the child-selection change (90f8df0) raised samaya's screening to 48 feasible
(from 46) at the same 7 solved. The next same-machine comparison will include it.

## Laptop runs (Dell G15, 26 September)

The same sets on the owner's laptop, never mixed with the cloud numbers above.
- **Machine:** Intel Core i5-13450HX (10 cores: 6 performance, 4 efficiency; 16 threads), 12 GB
  of memory in WSL2.
- **Solvers:** HiGHS 1.15.1, SCIP 10.0, CBC 2.10.11, GLPK 5.0. Every solver single-threaded with
  the 1e-4 gap, 4 instances at a time, as above. SCIP's "gaplimit" status (stopped at the
  requested gap) counts as solved.
- **Builds:** netlib and the cases use 51be19c (tree cuts and the MIP start); MIPLIB at 60 s uses
  0ffd224 (Feasibility Jump, before the tree cuts). The LP code is the same in both.
- **Load:** netlib and the cases ran while the 600 s MIPLIB run was going, so up to 8 solver
  processes shared the 10 cores. MIPLIB at 60 s ran alone. Every solver in a set ran under the
  same load.
- **Clock:** the WSL2 kernel clock on this laptop ran about 7% fast (90 s on it took 84 s of real
  time). The harness times every solver on that clock, so all laptop times read about 7% high
  and stay comparable with each other. But samaya times its limit on the same clock, while HiGHS,
  SCIP, CBC and GLPK use the real time: at the 60 s limit they stopped at a median of 63.5 s
  (CBC 64.7 s) as measured, samaya at 60.0 s. samaya thus had about 6% less time on every
  unsolved instance. One solve falls in that margin: HiGHS on neos-3627168-kasai at 63.4 s.
- No two solvers disagree on any instance, and every optimal samaya result passed the verifier.

### Netlib LP (93 instances, 300 s)

| Solver | Solved | Shifted geomean |
|---|---|---|
| HiGHS | 93/93 | 0.19 s |
| CBC (Clp) | 93/93 | 0.20 s |
| **samaya** | **93/93** | **0.36 s** |
| GLPK | 93/93 | 0.42 s |
| SCIP (SoPlex) | 93/93 | 0.47 s |

The order matches the cloud run, except that GLPK (93/93 after the e226 harness fix) moves ahead
of SCIP.

### MRPL case studies and generated MILPs (16 instances, 300 s)

| Solver | Solved | Shifted geomean |
|---|---|---|
| **samaya** | **16/16** | **3.43 s** |
| SCIP | 16/16 | 3.59 s |
| CBC | 16/16 | 4.02 s |
| HiGHS | 16/16 | 5.18 s |
| GLPK | 11/16 | 20.8 s |

samaya has the lowest shifted geomean here, about 4.5% below SCIP's; on the cloud machine it was
third. The two instances with the largest differences, which roughly cancel against SCIP:

| Instance | samaya | HiGHS | SCIP | CBC | GLPK |
|---|---|---|---|---|---|
| knapsack_80x5 | **17.2** | 284.0 | 55.0 | 30.9 | – |
| refsched_8c_6p_3u_26t | 122.6 | **47.0** | 51.3 | 132.2 | – |

(Seconds; "–" = not solved in 300 s.) The largest refinery-scheduling model is the one case
where both HiGHS and SCIP are more than twice as fast as samaya (elsewhere only on models solved
in under a second). On the nine MRPL cases every solver but GLPK finishes within 5 s.

### MIPLIB 2017, 62 small "easy" instances, 60 s

Two runs on this laptop. The final one (26-27 September, build 528b3ab: conflict analysis, the
cut-set flow covers) gave samaya the same real time as the others (`SAMAYA_CLOCK_FACTOR=1.052`,
measured against Windows' clock just before; see the clock note above). The earlier one (build
0ffd224, before tonight's work) ran with the plain limit, so samaya had about 6% less time.

| Solver | Solved | Shifted geomean | With a solution |
|---|---|---|---|
| SCIP | 19/62 | 38.10 s | 56 |
| HiGHS | 18/62 | 40.97 s | 57 |
| **samaya** (528b3ab) | **13/62** | **46.69 s** | 54 |
| CBC | 8/62 | 54.06 s | 55 |
| samaya (0ffd224, earlier run) | 7/62 | 50.33 s | 53 |

The baselines of the earlier run: SCIP 19, HiGHS 18, CBC 8, the same counts.

A rerun on the final build f0cbe85 (27 September, 06:40-07:38 IST) found the laptop's clock
steady (factor 1.000; the baselines stopped at a median of 60.00 s, samaya at 59.99 s): samaya
12/62 (46.85 s), HiGHS 18 (40.20 s), SCIP 19 (37.43 s), CBC 8 (53.64 s). The one difference from
the 528b3ab run is enlight_hard, solved there in 29 s and not here; its bound varies widely from
run to run (16% to 46% from the optimum after 60 s in tonight's A/Bs).

Solve times of every instance samaya solves (seconds, final run; "–" = not solved in 60 s):

| Instance | samaya | HiGHS | SCIP | CBC |
|---|---|---|---|---|
| markshare_4_0 | **15.1** | – | – | – |
| neos5 | **42.0** | – | – | – |
| mas76 | **19.1** | – | – | 33.5 |
| pk1 | **29.7** | – | – | 38.3 |
| exp-1-500-5-5 | **1.1** | 3.4 | 2.5 | – |
| neos859080 (infeasible) | **0.4** | 1.4 | 0.6 | – |
| p200x1188c | 3.0 | **0.6** | 4.1 | – |
| app1-1 | 13.5 | 18.9 | 6.4 | **4.7** |
| neos17 | 14.1 | 7.0 | **6.1** | – |
| enlight_hard | 29.1 | 11.6 | **0.0** | – |
| binkar10_1 | 49.4 | **21.0** | 26.5 | – |
| nu25-pr12 | 54.5 | 5.6 | **4.2** | 28.4 |
| sp150x300d | 0.6 | **0.1** | 0.3 | – |

- samaya solves four instances that neither HiGHS nor SCIP solves in 60 s: markshare_4_0,
  neos5, mas76 and pk1.
- HiGHS or SCIP solve 11 that samaya does not: beasleyC3, mc11, n5-3 (big-M networks, much
  closer since the cut-set fix), graph20-20-1rand and fhnw-binpack4-4 (symmetry),
  mik-250-20-75-4, neos-911970, pg, timtab1, neos-3381206-awhea and neos-3627168-kasai.
- No two solvers disagree on any instance, and every optimal samaya result passed the verifier.
  samaya's infeasibility claim on neos859080 is not certificate-checked (the verifier has no
  certificate for MILP infeasibility); HiGHS, SCIP and the MIPLIB reference agree with it.

Since this run: rotating aggregation starts and clique cuts (f0cbe85; bound A/B in
docs/MILP_PLAN.md).

### MIPLIB 2017, 62 instances, 600 s

HiGHS and SCIP from the full 600 s run (26 September, 15:30-21:44 IST; they stop on the real-time
clock, so they had 600 real seconds). samaya from a rerun alone on the final build f0cbe85
(27 September, 04:14-06:38 IST): the samaya results of the full run are not used because its
binary was rebuilt partway through (two builds in one run). samaya's limit was converted with
the clock factor measured before the run (1.071); over the run the monotonic clock ran 4.7-6.3%
fast (monotonic seconds per real second, from each shard's file times), so samaya had about
1-2% more real time than the baselines.

| Solver | Solved | Shifted geomean |
|---|---|---|
| HiGHS | 32/62 | 181.3 s |
| SCIP | 25/62 | 180.1 s |
| **samaya** | **17/62** | **293.2 s** |

Every instance samaya solves (seconds; "–" = not solved in 600 s):

| Instance | samaya | HiGHS | SCIP |
|---|---|---|---|
| mas74 | **323** | – | – |
| neos5 | **45** | 91 | – |
| markshare_4_0 | **16** | 123 | 110 |
| mas76 | **19** | 106 | 62 |
| pk1 | **33** | 209 | 120 |
| exp-1-500-5-5 | **1** | 3 | 2 |
| neos859080 (infeasible) | **0** | 1 | 1 |
| binkar10_1 | 40 | **29** | 43 |
| app1-1 | 13 | 18 | **6** |
| p200x1188c | 3 | **1** | 4 |
| sp150x300d | 1 | **0** | **0** |
| neos17 | 29 | 7 | **6** |
| nu25-pr12 | 57 | 5 | **3** |
| mik-250-20-75-4 | 96 | **14** | 35 |
| qap10 | 144 | 102 | **96** |
| pg | 212 | **6** | 21 |
| enlight_hard | 320 | 11 | **0** |

- samaya solves mas74, which neither HiGHS nor SCIP solves in 600 s, and is 4-8x faster than
  both on markshare_4_0, mas76 and pk1.
- HiGHS or SCIP solve 19 that samaya does not: beasleyC3, mc11, n5-3 (big-M networks),
  graph20-20-1rand, fhnw-binpack4-4 (symmetry), csched008, rococoC10-001000 (no first solution),
  gmu-35-40, gen-ip054, mcsched, fastxgemm-n2r6s0t2, neos-3381206-awhea, neos-3627168-kasai,
  neos-4738912-atrato, neos-911970, pg5_34, ran14x18-disj-8, timtab1 and tr12-30.
- No wrong answer: every optimal samaya result passed the verifier and matches the MIPLIB value.

MIPLIB with 6 threads (miplib600-t6): not run; at about 10 hours on this laptop it did not fit.

## Verified results: what samaya checks, and what the others do

Every answer samaya reports is checked by `src/verify/`, an independent checker that shares no
code with the solvers and recomputes everything from the original model (before presolve) in
`long double`. It is on by default, and the JSON output carries the outcome (`verified`,
`max_primal_violation`, `max_dual_violation`).

| Claim | Evidence checked |
|---|---|
| LP optimal | Primal feasibility, and dual feasibility of the reduced costs d = c − Aᵀy recomputed from the returned duals (with the objective sense) |
| LP infeasible | A Farkas certificate: the returned row multipliers must prove that no point satisfies the rows and bounds |
| LP unbounded | A ray: it improves the objective and never violates a finite bound or row |
| MILP solution | Bounds, rows and integrality of the returned point on the original model; the objective is recomputed |

If an LP answer fails the check, samaya re-solves with tighter tolerances and then without
scaling. An answer that still fails is reported as `numerical_error`, never as a result, and
the checker's tolerances are never loosened to make it pass.

**What this does not cover:** the MILP bound. "Optimal" for a MILP means the search closed the
gap to 1e-4, and nothing re-checks that the pruned subtrees held nothing better.

**What the other solvers do:**
- **HiGHS** reports its own primal violations (number, sum and maximum) and corrects the model
  status when a claimed optimum breaks its tolerances. That check runs inside the solver.
- **SCIP 10** added an exact rational mode that writes VIPR certificates. An independent tool
  checks these, including the MILP bound, which makes it the strongest guarantee of the four.
  It is optional, and 7–10× slower than SCIP's default settings, per the SCIP 10 report.
- **CBC and GLPK** check their own solutions.

Mature solvers still ship wrong answers now and then. HiGHS examples:
- [ERGO-Code/HiGHS#3273](https://github.com/ERGO-Code/HiGHS/issues/3273): HiGHS 1.15.1 cut off
  the optimal solution with presolve off (faulty conflict analysis). Opened September 2026,
  fixed.
- [ERGO-Code/HiGHS#1710](https://github.com/ERGO-Code/HiGHS/issues/1710): HiGHS 1.7.0 reported
  an infeasible MILP as optimal with presolve on.

samaya's verifier would not have caught the first: it is a bound error. It would have caught
the second: the returned point fails the row check.

**The honest claim:** every samaya result is checked by default, at the cost of one pass over the
matrix, against certificates for everything except the MILP bound. The benchmark harness then
cross-checks every final answer against the other solvers' and against the MIPLIB reference
values. There has been no disagreement in the runs reported here.
