# Numerical robustness

SIH26119 asks for "a clear demonstration of numerical robustness by solving challenging
large-scale optimization problems involving degeneracy, weak LP relaxations or ill-conditioned
constraint matrices". This page collects that evidence.

Every answer below passed samaya's independent verifier (`src/verify/`):
- It recomputes row activities and reduced costs from the original model, in extended
  precision.
- It checks bounds, rows, integrality and the optimality conditions, or the Farkas certificate
  of an infeasible model.
- Its tolerances are 1e-6 relative and are never loosened.

A solve whose answer fails the verifier is re-solved with tighter settings, and if it still
fails it is reported as `numerical_error`, never as optimal. **No wrong answer has been reported
on any instance in any run** (Netlib, the infeasible Netlib set, MIPLIB 2017, the MRPL cases,
Maros–Mészáros).

Machine: Dell G15 laptop (Intel Core i5-13450HX), WSL2. Times are from the laptop comparison
run in [comparison.md](comparison.md); the residuals and bounds on this page are from runs of
build 9343abf on 27 September 2026.

## 1. Ill-conditioned matrices

The 12 Netlib models with the widest coefficient ranges (largest |a_ij| over the smallest
nonzero |a_ij|), default settings. Residuals are the verifier's largest relative violations.

| Model | Coefficient range | samaya objective | Primal residual | Dual residual |
|---|---|---:|---:|---:|
| pilotnov | 3.0e12 | -4497.2761882 | 6.1e-12 | 6.1e-16 |
| pilot.ja | 3.0e12 | -6113.1364656 | 2.6e-11 | 8.8e-14 |
| pilot87 | 1.0e9 | 301.71036566 | 6.9e-11 | 5.0e-08 |
| tuff | 1.0e9 | 0.29214776509 | 7.5e-16 | 2.9e-18 |
| pilot4 | 7.5e8 | -2581.1392589 | 3.7e-12 | 2.7e-12 |
| pilot.we | 5.0e8 | -2720107.5328 | 1.4e-12 | 2.0e-14 |
| perold | 4.0e8 | -9380.7552782 | 8.1e-12 | 3.2e-12 |
| maros | 2.0e8 | -58063.743701 | 6.8e-13 | 3.0e-14 |
| pilot | 1.0e8 | -557.48970616 | 9.0e-12 | 8.6e-09 |
| cycle | 9.0e7 | -5.2263930249 | 1.5e-16 | 1.5e-14 |
| wood1p | 3.3e7 | 1.4429024116 | 2.6e-14 | 2.4e-16 |
| agg3 | 2.0e7 | 10312115.935 | 2.2e-16 | 8.9e-13 |

All 12 are verified on the first attempt; no re-solve was needed on any of the 93 Netlib models.
On pilot and pilot87 samaya's optimum differs from the published value by 4e-8 and 6e-8
relative (−557.4897062 vs −557.4897293, 301.7103657 vs 301.7103473), well inside the 1e-6
tolerance. HiGHS and SCIP report the published values; GLPK reports the same values as samaya.

## 2. Degenerate LPs

The Netlib models known for degeneracy (cycling risk in the simplex), with all five solvers.
Every solver reports the same optimum to the digits shown (pilot and pilot87 within 1e-7
relative, as above).

| Model | Optimum | samaya | HiGHS | SCIP | CBC | GLPK |
|---|---:|---:|---:|---:|---:|---:|
| degen2 | -1435.178 | 0.02 s | 0.02 s | 0.06 s | 0.03 s | 0.03 s |
| degen3 | -987.294 | 0.41 s | 0.27 s | 1.29 s | 0.24 s | 0.41 s |
| cycle | -5.2263930 | 0.08 s | 0.44 s | 0.16 s | 0.04 s | 0.09 s |
| greenbea | -72555248.13 | 0.73 s | 0.21 s | 3.88 s | 0.72 s | 0.60 s |
| greenbeb | -4302260.261 | 1.73 s | 0.63 s | 2.34 s | 0.63 s | 0.72 s |
| pilot | -557.48973 | 2.01 s | 1.63 s | 4.47 s | 1.32 s | 2.96 s |
| pilot87 | 301.71035 | 7.93 s | 4.46 s | 7.78 s | 4.07 s | 6.38 s |
| perold | -9380.755278 | 0.05 s | 0.06 s | 0.33 s | 0.10 s | 0.07 s |
| d2q06c | 122784.2108 | 1.87 s | 1.29 s | 5.17 s | 2.03 s | 2.09 s |

The dual simplex handles degeneracy with cost perturbation and cost shifting, a bound-flipping
("long step") ratio test and Harris's two-pass ratio test with dual steepest-edge pricing.

## 3. Infeasibility, proven with certificates

On the 29 infeasible Netlib models, samaya proves 28 infeasible. Each of the 28 comes with a
Farkas certificate y: the verifier checks that y'Ax − y'r excludes zero over all column bounds
and row bounds.

The 29th, cplex2, ends in `numerical_error`: the certificate the simplex found did not pass the
verifier, so samaya does not claim infeasibility. HiGHS reports infeasible. Details:
[netlib.md](netlib.md).

For MILP, samaya's conflict analysis learns from infeasible node LPs through their Farkas rays.
In the test suite every ray (265 of 265 on the lot-sizing family) is accepted by the verifier
before it is used.

## 4. Weak LP relaxations

MIPLIB 2017 instances whose LP relaxation is far from the integer optimum. The root bound is
the LP relaxation before and after samaya's cutting planes: Gomory, c-MIR with path
aggregation, flow covers including cut-set covers on network rows, knapsack covers and
cliques. The last column is the result at 60 s.

| Instance | Optimum | Root LP bound | After root cuts | Gap closed at the root | samaya at 60 s |
|---|---:|---:|---:|---:|---|
| p200x1188c | 15078 | 5678.6 | 12828.7 | 76% | optimal, 5.7 s |
| exp-1-500-5-5 | 65887 | 28427.0 | 65289.7 | 98% | optimal, 2.2 s |
| sp150x300d | 69 | 4.9 | 66.4 | 96% | optimal, 1.2 s |
| mc11 | 11689 | 608.8 | 10697.0 | 91% | time limit |
| beasleyC3 | 754 | 40.4 | 639.5 | 84% | time limit |
| timtab1 | 764772 | 28694 | 329034 | 41% | time limit |
| markshare_4_0 | 1 | 0 | 0 | 0% (by construction) | optimal, 31.5 s |
| neos-911970 | 54.76 | 23.26 | 23.26 | 0% | time limit |

- markshare_4_0 is designed so that no cut helps (its LP bound is 0). samaya proves it optimal
  by search in 31.5 s; neither HiGHS nor SCIP does so within 60 s
  ([comparison.md](comparison.md)).
- Where samaya stops at the time limit, the solution it reports still passes the verifier; only
  optimality is not proven.

## 5. Convex QP (Maros–Mészáros)

A screening run on the 138 Maros–Mészáros convex QPs, with the interior-point method and 60 s
per model:
- **105 of 138 solved and verified** by `verify_qp_optimality`, which recomputes
  d = c + Qx − A'y.
- **0 wrong answers.**
- VALUES is refused as not convex. Its Q has a negative eigenvalue (−1.27e-5, computed
  independently with numpy); samaya's check factorizes Q + εI, and a nonpositive pivot means an
  indefinite Q. HiGHS returns a point for VALUES.
- **14** end in `numerical_error`: the interior point converged, but neither its point nor the
  polished point passed the strict verifier. These are mostly the degenerate, LP-like
  Q-versions of Netlib models (QSHIP*, QSCFXM*, QSHARE1B); a crossover would finish them.
- **14** reach the iteration limit (the LISWET family and YAO), and **3** the time limit.

Infeasible QPs get a verified Farkas certificate from the simplex on their constraints. The
same interior point solves LPs (`--lp-method barrier`): on Netlib, 71 of 93 answers are
verified from the interior point alone and the dual simplex completes the other 22, all 93
verified. The timed run against HiGHS's QP solver goes in [qp.md](qp.md).

## How to reproduce

```sh
cmake --preset release && cmake --build --preset release
build/release/apps/cli/samaya --json bench/instances/netlib/pilot87.mps      # section 1
bench/compare.sh netlib 4                                                    # section 2
bench/fetch_instances.sh netlib-infeas                                       # section 3
build/release/apps/cli/samaya --time-limit 60 bench/instances/miplib_small/mc11.mps   # section 4
```
