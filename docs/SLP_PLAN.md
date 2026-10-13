# SLP: blending qualities by successive linear programming (FINALE_PLAN P2.3)

Refinery planning tools (Aspen PIMS and similar) model blend qualities with bilinear terms:
the quality of a pool times the flow out of it. They solve these models by successive linear
programming ("distributive recursion" in PIMS, SLP in general). samaya solves LP, MILP and
convex QP, but not these. This plan adds them: quadratic terms in constraint rows, read from MPS
in the CPLEX/Gurobi QCMATRIX convention, solved by a trust-region penalty SLP on the dual
simplex, with every LP verified and the final point checked by an independent KKT verifier.

The first draft was reviewed by Codex (astra, high) on 13 October; this version takes its
corrections: what the final check certifies (first-order stationarity, not local optimality),
an independent KKT verifier, complete initialization, a precise merit function and penalty
rule, a stopping rule that does not mistake a collapsed trust region for convergence, and the
Haverly test facts.

## 1. What the answer means

Bilinear models are nonconvex. SLP finds a **feasible, first-order stationary point (KKT
point)** from its start; that is what PIMS reports too. Neither SLP nor its KKT check certifies
a local optimum (at u = v = 0, min z with z = uv is a KKT point although z = -t² is feasible
nearby), let alone a global one. So:

- New `Status::kFirstOrderStationary` ("first_order_stationary"), appended to the enum (and the
  C API's status list) without renumbering. The report says "feasible, first-order stationary;
  local optimality not certified". Never "optimal".
- The row duals reported are the **local KKT multipliers** at that point: the marginal values
  planners read, valid for that local solution under the usual regularity, possibly not unique.
  `dual_bound` stays unknown.

## 2. Model and input

- **Model** (`include/samaya/model.hpp`): `struct RowQuadraticTerm { Index row; Index col1;
  Index col2; double value; };`, `std::vector<RowQuadraticTerm> row_quadratic;` with
  col1 <= col2, one entry per pair; row i's activity is (A x)_i + Σ value · x_col1 · x_col2.
  `ProblemClass::kQCP` when there are row terms. Before dispatch: row terms together with
  integer columns, or with a quadratic objective (Q), return `Status::kNotImplemented` with a
  message (the MILP and QP code would otherwise ignore the rows).
- **MPS** (`src/io/mps_reader.cpp`): `QCMATRIX <row>` sections in the **CPLEX/Gurobi
  convention**: the full symmetric matrix Q of x'Qx (no factor 1/2; MOSEK's convention differs
  and is not supported), both (i, j) and (j, i) listed with equal values. Stored as
  value = Q_ii on the diagonal and Q_ij + Q_ji off it. Its own parser (the objective's
  `parse_quadratic` keeps one triangle and does not check symmetry). Errors with a message: an
  unknown row or column, a missing or unequal mirror entry.
- **Start:** `Params::mip_start` becomes the starting point for QCPs too (NaN = unknown); its
  comment and the CLI help say so. No new field.
- **Generators:** `cases/` gets the Haverly instances and a pooling/blending case written as
  MPS with QCMATRIX.

## 3. Verification

- **`verify_primal`** includes the quadratic terms in each row's activity, summed in extended
  precision, with the existing fixed tolerances.
- **New `verify_qcp_kkt(model, x, y)`** in `src/verify/`: recomputes the activities and the
  Jacobian J(x) **from the original model** (not from the SLP's linearization, so a wrong
  linearizer cannot pass with its own LP), forms d = ∇f - J(x)'y, and checks: x feasible;
  row duals' signs and complementarity against the true activities; reduced costs' signs
  against the **original** bounds only. Same tolerance constants as `verify_lp_optimality`.
- A QCP answer is **verified** when x passes `verify_primal` and (x, y) passes
  `verify_qcp_kkt`, with y from the last LP whose own status is a verified optimum (the
  `verified` flag alone is not enough: the helper passes limit statuses and runs nothing when
  `Params::verify` is off). If the final point is feasible but not KKT, the status is a limit
  status with the message "stalled before stationarity", never stationary.

## 4. The algorithm (`src/slp/slp.{hpp,cpp}`)

Penalty SLP with a trust region (Zhang, Kim and Lasdon 1985, Management Science 31(10), which
treats refinery models; l1 merit as in Fletcher's Sl1QP). The update rules below are ours; the
paper's convergence result is not claimed for them.

1. **Initialization.** Start from `mip_start` where given; else each column in a bilinear term
   at the midpoint of its finite bounds (or its finite bound, or 0). Then a **projection LP**:
   minimize Σ |x_j - x0_j| (weighted by the scales below) over all linear rows and bounds,
   with the nonlinear rows left out. Its answer is the first iterate. If that LP is infeasible,
   the model is infeasible (verified certificate on the linear part); that is the only
   infeasible verdict SLP gives.
2. **Scales:** s_j = max(1, |x_0j|) unless the column has a finite range, then that range
   (capped at 1e6 · max(1, |x_0j|)). Trust bounds |x_j - x_k,j| <= δ · s_j on every column in a
   bilinear term (diagonal terms included), intersected with the original bounds; fixed
   columns get none. δ starts at 0.1.
3. **Linearize** at x_k: v · u · w → v · (u_k w + w_k u - u_k w_k); v · u² → v · (2 u_k u - u_k²).
4. **LP_k:** linear rows; linearized rows with elastics p_i, n_i >= 0 at cost μ · w_i, with
   row weights w_i = 1 / max(1, ||gradient of row i at x_0||_inf) fixed at the start; the
   trust bounds. Objective σ c'x with σ = +1 for minimization, -1 for maximization, so the
   elastics are always penalized. Dual simplex, warm-started from LP_(k-1)'s statuses (the
   factorization is rebuilt, since coefficients change). Remaining time passed down.
5. **Merit and ratio:** V(x) = Σ w_i dist(h_i(x), [rl_i, ru_i]) over the nonlinear rows,
   Φ(x) = σ c'x + μ V(x). Predicted reduction pred = -σ c's + μ (V(x_k) - V_lin(s)), computed
   from the step (no objective offsets subtracted). If pred is not positive or negligible, go
   to the stationarity check instead of dividing. ρ = (Φ(x_k) - Φ(x_k + s)) / pred; accept when
   ρ > 0.1; δ × 2 when ρ > 0.75 and the step reached the trust region, δ / 4 when ρ < 0.25 or
   on rejection (keeping x_k). Whenever μ changes, Φ and pred are recomputed with the new μ.
6. **Penalty steering:** when LP_k leaves elastics > 0, solve a **feasibility LP** (elastics
   only, same center and radius). If it can reduce V appreciably, raise μ by 10 until LP_k
   achieves at least 10% of that reduction (cap 1e8, an operational limit). If even the
   feasibility LP cannot reduce V, report restoration failure (limit status, message "no
   feasible point found from this start"), never infeasible.
7. **Second-order correction** after a rejected step that is nearly feasible: one LP at the
   trial point with the trial point's Jacobian, the hard rows and bounds, small radius,
   minimizing the nonlinear residual; accepted only on the true merit.
8. **Stopping:** max_j |s_j| / s_j small (1e-7), V(x) within the verifier's tolerance, and the
   independent KKT check passing → `kFirstOrderStationary`. A radius below 1e-10 or the
   iteration limit (default 200) without passing it → `kIterationLimit` with the message
   "stalled before stationarity". The time limit as everywhere else.
9. **Empty pools** (zero throughput leaves a quality undetermined) are a known SLP difficulty;
   the generator's cases avoid them at the start, and the Studio report says when a pool's
   flow is zero.

Logging: level 1 a summary line, level 2 one line per iteration (k, objective, V, δ, μ, ρ,
simplex iterations), as `Params::log_level` documents. Results: x and the original objective
and activities recomputed from the model, elastics removed; `slp_iterations` and
`max_nonlinear_violation` (relative, as the verifier measures it) in `Result` and the JSON.

**Plumbing:** `solve_lp_verified` (anonymous namespace in `src/core/solver.cpp`) moves to a
small internal header so `src/slp/` can call it; `LpResult` gains the final basis statuses
from the simplex that produced the answer, for the warm starts; the dispatch hunk in
`solver.cpp` stays small.

## 5. Tests

- **Haverly 1–3** in the p-formulation, coordinates (a, b, u, v, z1, z2, p): a, b into the pool,
  u, v from the pool to products 1 and 2, z1, z2 direct component C, p the pool's sulfur:

      min 6a + c_B b + 10(z1 + z2) - 9(u + z1) - 15(v + z2)
      a + b = u + v;  3a + b = p (u + v);
      p u + 2 z1 <= 2.5 (u + z1);  p v + 2 z2 <= 1.5 (v + z2);
      u + z1 <= D1;  v + z2 <= 200;  flows >= 0;  1 <= p <= 3

  | Instance | D1 | c_B | Global min | A global solution | Trap start | Local point | Local value |
  |---|---|---|---|---|---|---|---|
  | Haverly 1 | 100 | 16 | -400 | (0,100,0,100,0,100,1) | (49,0,49,0,49,0,3) | (50,0,50,0,50,0,3) | -100 |
  | Haverly 2 | 600 | 16 | -600 | (300,0,300,0,300,0,3) | (0,99,0,99,0,99,1) | (0,100,0,100,0,100,1) | -400 |
  | Haverly 3 | 100 | 13 | -750 | (50,150,0,200,0,0,1.5) | (74.25,24.75,99,0,0,0,2.5) | (75,25,100,0,0,0,2.5) | -125 |

  The global values are published (GAMS pooling library, SAS pooling example); the trap starts
  and local points are the reviewer's derivations, to be confirmed by the tests themselves.
  Tests: from starts near a global solution, the global value; from the trap starts, the
  stationary local point and value, reported as stationary (never optimal); `verify_qcp_kkt`
  passes at every listed point with the right multipliers and fails when one multiplier is
  perturbed.
- **The z = uv example:** reported stationary, not optimal.
- **KKT verifier:** analytical fixtures with hand-computed multipliers; a wrong Jacobian entry
  in the candidate y is caught.
- **MPS QCMATRIX:** diagonal 2 and 5 with symmetric off-diagonal 3 give 2u² + 6uv + 5v², which
  is 89 at (2, 3); L, G and E rows; several QCMATRIX blocks; a diagonal-only row; unknown
  column; unequal and missing mirror entries refused.
- **Initialization:** u, v in [0, 10] with u + v = 1 and a bilinear row (the midpoint start
  (5, 5) is infeasible for the linear row; the projection LP must fix it).
- **Planted bugs:** the linearization's constant term with the wrong sign; one column of a
  term left without trust bounds (caught by a direct check on LP_k's bounds, not by Haverly's
  objective); σ dropped for a maximization (elastics rewarded).

## 5a. What the implementation changed (13 October)

Found while testing Haverly 1–3 from the default start:
- **Empty pools at the start.** Flows at their lower bound 0 leave the pool quality without a
  first-order effect, and the first LP stops at objective 0, a valid but useless KKT point.
  Columns in quadratic terms with one finite bound now start one unit inside it.
- **Scales follow the iterate.** A column without a usable finite range is scaled by
  max(initial scale, |x_k|), not by its value at the start: flows grown from 1 to 100 were
  still limited to steps of 0.2 (Haverly 3 crawled for 200 iterations).
- **Second-order correction is needed, not optional.** Without it, steps that satisfy the
  linearized rows leave a curvature residual the next LP again predicts to remove; the ratio
  stayed near 0.63, the radius never grew, and Haverly 3 stalled at -718 of -750. With it:
  -750 in 12 iterations. The planted "correction never taken" bug fails that test.
- **The projection pulls every column the caller gave a value for,** not only the columns in
  quadratic terms; otherwise a linear column lands on an arbitrary vertex and the loop starts
  with a violation the start did not have.

Open: columns whose linearized gradient is zero (u, v near u = v = 0 in z = uv) can jump to
arbitrary trust-region corners. A small L1 proximal term (|x - x_k| through two more columns per
such column) would keep them at the center. Haverly does not need it; add it when a case does.

Results now (Windows release): Haverly 1, 2, 3 from near-global starts: -400, -600, -750 (the
global optima); from the trap starts: -100, -400, -125 (stationary, as derived); from the
default start: -400 (global), -400 (stationary local), -750 (global). All verified by
`verify_qcp_kkt`.

## 6. Who does what

| Part | Who |
|---|---|
| Model fields, ProblemClass, Status (+ C API), dispatch and rejections | Claude |
| `verify_primal` with quadratic rows, `verify_qcp_kkt` | Claude |
| The SLP loop (sections 4.1–4.9), the plumbing | Claude |
| MPS QCMATRIX parser and its tests (section 2 and 5) | Codex (sol/high) |
| Haverly 1–3 and trap fixtures as MPS, the pooling/blending generator and demo case | Codex (sol/medium) |
| JSON, CLI and Studio display of SLP results | Codex (sol/medium) |
| Review of the finished SLP branch | Codex (astra/high) |

## 7. Gate

Haverly 1–3: global values from near-global starts, the stationary local points from the trap
starts; every LP in the loop a verified optimum; the final point passing `verify_primal` and
`verify_qcp_kkt`; planted bugs caught; the four presets.
