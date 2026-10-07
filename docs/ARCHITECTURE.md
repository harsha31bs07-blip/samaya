# samaya: architecture and work packages

This document is for anyone joining the solver work. The first half explains how the code fits
together today and the rules every change follows. The second half splits the remaining work
into self-contained **work packages** (WP) that can be built in parallel without stepping on each
other. [PLAN.md](../PLAN.md) holds the overall 16-week plan; this document is the engineering
view of the same plan, kept in sync with the code.

If you take a work package, read sections 1–5, then your package in section 6.

---

## 1. Where we are

| Area | State | Evidence |
|---|---|---|
| Model, MPS/QPS reader (free and fixed format), C/C++ API, CLI | done | unit tests |
| Scaling (geometric mean + equilibration, powers of two) | done | unit tests |
| Sparse LU (Markowitz, threshold pivoting, Forrest–Tomlin updates) | done | tests against a dense LU |
| Dual simplex (dual steepest edge, bound flipping, Harris, perturbation) + primal cleanup | done | about 2,400 random LPs versus a reference simplex |
| Independent verifier: optimality, Farkas certificates, unbounded rays | done | every solve runs it |
| LP presolve + postsolve (primal and dual) | done | planted-reduction tests |
| MILP branch-and-bound, propagation, reliability branching, heuristics | done | 1,900 random MILPs versus a reference MILP solver |
| Root cuts: Gomory mixed-integer, c-MIR, knapsack covers | done | cuts checked against known optima |
| **Netlib LP** | **93/93 optimal, all objectives match HiGHS**; 28/29 infeasible set proven | [docs/results/netlib.md](results/netlib.md) |
| Generated MILPs (refinery scheduling, knapsack, facility location) | 7/7 match HiGHS | [docs/results/generated-mip.md](results/generated-mip.md) |
| MIPLIB 2017, 62 small "easy" instances, 10 min each | 8/62 with root cuts (7 without), 0 wrong answers | [docs/results/miplib.md](results/miplib.md) |
| **Barrier (interior point), PDLP, GPU, QP** | **not started** | these are the work packages below |

The core is dependency-free C++20. **No external solver or linear-algebra library is linked**
(no BLAS, no SuiteSparse, no cuSPARSE for factorizations). HiGHS is used only by the benchmark
harness, as a comparison baseline.

---

## 2. Build, test, benchmark

```sh
cmake --preset release && cmake --build --preset release && ctest --preset release
cmake --preset debug   && cmake --build --preset debug      # warnings are errors
cmake --preset asan    && cmake --build --preset asan       # AddressSanitizer + UBSan
cmake --preset cuda                                         # reserved for WP2 (SAMAYA_CUDA=ON)
```

CI (`.github/workflows/ci.yml`) builds `{gcc, clang} x {debug, release, asan}` and runs all
tests on every pull request. All six jobs must be green before merging.

```sh
build/release/apps/cli/samaya --json model.mps                 # solve, verify, JSON summary
build/release/apps/cli/samaya --lp-method pdlp model.mps       # method switch (WP1 wires it)
bench/fetch_instances.sh netlib netlib-infeas                  # public instances
bench/fetch_instances.sh miplib-list bench/miplib_small.test
bench/generate_lps.py --set lp; bench/generate_lps.py --set mip
bench/harness.py bench/instances/netlib --baseline highspy --time-limit 300
bench/report.py bench/results/<file>.csv --title "..."         # Markdown table
```

`pip install highspy` gives the comparison baseline. Instances are downloaded and never committed.

---

## 3. How a solve flows

```
 MPS file ──> io/mps_reader ──> Model ──> Solver::solve (core/solver.cpp)
                                              │
                    validate, crossed bounds ─┤
                                              ├── LP ───> presolve ──> solve_lp_verified ──> postsolve
                                              │                          │  (lp/lp_solver: scale,
                                              │                          │   dual simplex, unscaled
                                              │                          │   cleanup, unscale)
                                              │                          └── verifier; retry with tighter
                                              │                              tolerances, then unscaled
                                              ├── MILP ─> presolve (MIP mode) ──> BranchAndBound
                                              │              (root cuts, node LPs via dual simplex)
                                              │           ──> postsolve x ──> verify_primal
                                              └── QP ───> not_implemented   (WP4)
                                                            │
                                              Result: status, objective, x, y, d, certificates,
                                                      verified flag, violations, iterations, nodes
```

**Rule: nothing is reported without the verifier.** An LP optimum must pass
`verify_lp_optimality`, an infeasibility claim `verify_infeasibility` (a Farkas certificate), an
unbounded claim `verify_unbounded_ray`, and a MILP solution `verify_primal` (bounds, rows,
integrality). An outcome that fails becomes `numerical_error`. New solvers plug in *before* the
verifier, so they inherit this guarantee for free.

### Directory map

| Path | Contents |
|---|---|
| `include/samaya/` | Public API: `model.hpp`, `solver.hpp` (`Result`), `params.hpp`, `status.hpp`, `io.hpp`, `verify.hpp`, `sparse_matrix.hpp`, `types.hpp` |
| `include/samaya_c.h` | Stable C API |
| `src/core/` | `Solver` dispatch (`solver.cpp`), model validation and statistics, logging (`Logger`, `Timer`), C API |
| `src/io/` | MPS/QPS reader |
| `src/linalg/` | `SparseMatrix` (CSC), `Scaling`, `BasisFactor` (LU with FT updates) |
| `src/lp/` | `Simplex` (dual + primal), `lp_solver` (scaling, cleanup, unscaling, duals and rays) |
| `src/presolve/` | `Presolve` (reductions + postsolve stack) |
| `src/mip/` | `BranchAndBound`, `cuts` |
| `src/verify/` | Independent checks in `long double` |
| `apps/cli/` | `samaya` command-line tool |
| `tests/` | Unit and differential tests, reference solvers, generators |
| `bench/` | Harness, instance generator, download script, report writer |
| `docs/` | This file, benchmark results |

---

## 4. Contracts you must respect

### 4.1 Model

```cpp
struct Model {                       // include/samaya/model.hpp
  ObjSense sense;  double obj_offset;
  std::vector<double> obj, col_lower, col_upper;  std::vector<VarType> col_type;
  std::vector<double> row_lower, row_upper;       // row_lower <= A x <= row_upper
  SparseMatrix A;                    // m x n, CSC (col_start / row_index / values)
  SparseMatrix Q;                    // lower triangle of Q; empty for LP/MILP
};
```

Infinite bounds are `±kInf` (`types.hpp`); values beyond ±1e30 in files map to infinity. `Index`
is `int32_t` and `NnzIndex` is `int64_t`. An equality row has `row_lower == row_upper`.

### 4.2 Sign conventions (the most common source of bugs)

- **Reduced costs:** `d = c - Aᵀy`, always computed from the original model.
- **Minimization optimality:** a column strictly between its bounds has `d = 0`; at its lower
  bound `d ≥ 0`; at its upper bound `d ≤ 0`.
- **Row duals** behave like reduced costs of the row activity `r = a_i x`: at the row's lower
  bound `y_i ≥ 0`, at its upper bound `y_i ≤ 0`, strictly inside `y_i = 0`.
- **Maximization:** the same conditions hold for `sense · d` and `sense · y` (`sense = -1`). The
  verifier implements exactly this; test your solver against it.
- **Farkas certificate:** a row vector `y` for which `yᵀAx − yᵀr` cannot be zero for any `x` and
  `r` within their bounds.
- **Unbounded ray:** a direction that improves the objective and keeps every bound and row
  feasible.

### 4.3 Internal LP form (what the simplex sees)

```cpp
struct LpProblem {                  // src/lp/simplex.hpp
  Index m, n;  SparseMatrix A, At;  // scaled A and its transpose (row-wise)
  std::vector<double> cost, lower, upper;   // length n + m
};  // min costᵀv  s.t.  [A  -I] v = 0,  lower <= v <= upper;  v = (x, row activities)
```

`lp_solver.cpp` builds it from a `Model`: it negates the objective for maximization, applies
`Scaling` (R A C with power-of-two factors, exact to unscale), solves, runs an **unscaled
cleanup** from the optimal basis, and maps `x`, `y`, `d`, rays and certificates back. New LP
methods should reuse the scaling and the result mapping instead of writing their own.

### 4.4 Result

`Result` (`include/samaya/solver.hpp`) has `status`, `objective`, `dual_bound`, `col_value`,
`row_activity`, `row_dual`, `col_dual`, `infeasibility_certificate`, `unbounded_ray`,
`verified`, `max_primal_violation`, `max_dual_violation`, `simplex_iterations`,
`barrier_iterations` and `nodes`. **Fill `barrier_iterations` from WP1/WP3** (PDLP iterations
count there too).

### 4.5 Parameters already reserved for you

`Params::lp_method` (`kAuto`, `kDualSimplex`, `kPrimalSimplex`, `kBarrier`, `kPdlp`,
`kConcurrent`), `Params::use_gpu` and `Params::threads`, plus the CLI flags `--lp-method` and
`--gpu`, already exist. `solve_lp_model` in `src/core/solver.cpp` currently logs "not available
yet" for barrier and PDLP. **That `if` is your integration point.**

---

## 5. Engineering rules

1. **Correctness before speed.** Every new algorithm gets a differential test against an
   independent reference: `tests/reference_lp.hpp` (dense two-phase simplex in `long double`) and
   `tests/reference_milp.hpp`. Compare statuses and objectives on hundreds of random models from
   `tests/lp_generators.hpp` (feasible, infeasible, unbounded, degenerate and boxed families).
2. **Prove your tests can fail.** Before trusting a green suite, plant a bug (flip a sign,
   drop a term) and confirm a test catches it. Every existing module was checked this way; a
   test that cannot fail is not evidence.
3. **Fallbacks must not hide bugs.** When a solver falls back (retry, cold start, other
   method), add a test that bypasses the fallback or counts how often it fires.
4. **Warning-free** on GCC and Clang in debug (`-Werror`), and clean under ASan/UBSan.
5. **Tests use the built-in framework** (`tests/test_framework.hpp`: `TEST`, `CHECK`,
   `CHECK_EQ`, `CHECK_NEAR`, `CHECK_THROWS`, `REQUIRE`), because GoogleTest cannot be downloaded
   in every environment. Add files to `tests/CMakeLists.txt`; tests may include `src/` headers.
6. **Add sources to the `add_library(samaya ...)` list in `CMakeLists.txt`.** CUDA sources go
   behind `if(SAMAYA_CUDA)`, with a CPU fallback always built.
7. **Style:** match the surrounding code. That means 2-space indent, `snake_case` functions,
   `trailing_underscore_` members, `kConstants`, a comment on every non-obvious numerical choice
   (tolerances, thresholds) and named constants instead of magic numbers.
8. **Benchmarks are reported honestly.** Unsolved instances count at the time limit (shifted
   geometric mean, shift 10 s). Wrong answers are listed, not hidden. Run benchmarks on an idle
   machine and never rebuild the binary under test while a benchmark runs.
9. **Workflow:** one branch per work package, small commits, a PR into `main` when a milestone's
   tests pass on all presets. The PR description states what works, the test results and the
   known gaps.

---

## 6. Work packages

Each package lists its goal, where it plugs in, the algorithm, its milestones and the definition
of done. They touch different directories, so they can run in parallel. **WP1 → WP2** and
**WP3 → WP4** are the only hard dependencies.

### WP1: PDLP on the CPU (`src/lp/pdlp.{hpp,cpp}`)

**Why:** PDLP is the algorithm that makes the GPU story real (WP2). Building it on the CPU first
gives a correct reference to test the GPU kernels against.

**Interface:**

```cpp
struct PdlpOptions { double tol = 1e-8; long long max_iterations = -1; double time_limit = kInf;
                     bool use_gpu = false; };
struct PdlpResult  { SimplexStatus status; std::vector<double> x, y; long long iterations;
                     double primal_residual, dual_residual, gap; };
PdlpResult solve_pdlp(const LpProblem& lp, const PdlpOptions&, const Logger&);
```

Work on the scaled `LpProblem` so that `lp_solver.cpp` can reuse its unscaling. Map the logical
bounds to row bounds of the saddle-point problem.

**Algorithm** (Applegate et al., "Practical large-scale linear programming using primal-dual
hybrid gradient", NeurIPS 2021; the cuPDLP papers of Lu and Yang, 2023):
1. Preconditioning: Ruiz equilibration (about 10 passes), then Pock–Chambolle (α = 1), applied
   on top of our scaling.
2. PDHG iteration with **adaptive step size** and **primal weight** updates.
3. **Adaptive restarts** based on the normalized duality gap, restarting to the average or the
   current iterate.
4. Termination by relative primal residual, dual residual and gap ≤ `tol`. Detect
   infeasibility from the difference of iterates (this gives a ray / Farkas candidate that the
   verifier then checks).
5. Optional feasibility polishing.

**Milestones** (step-by-step prompts in [GETTING_STARTED.md](GETTING_STARTED.md)):
- **M1:** plain PDHG solves tiny LPs.
- **M2:** restarts, adaptive steps and primal weight added; matches the reference simplex to
  1e-6 relative on the random LP families.
- **M3:** infeasibility and unboundedness detection, with rays that pass the verifier.
- **M4:** `--lp-method pdlp` wired into `solve_lp_model`, results verified. PDLP output is lower
  accuracy, so a PDLP point that fails the strict verifier is polished by the simplex from a
  basis guessed from it (see WP3 M5 for the full crossover); the verifier is never loosened
  silently.
- **M5:** Netlib at `tol = 1e-4` and `1e-8`, reporting iterations and time per instance.

**Done when:** the random families pass, Netlib runs at `tol = 1e-4` with a solved-count report,
the verifier is wired in, and all presets are green.

**Pitfalls:** free variables and infinite bounds in the projection; the step size must respect
‖A‖; do not compute ‖A‖ exactly (use power iteration); test degenerate and infeasible families,
not just feasible ones.

### WP2: PDLP on the GPU (`src/gpu/`, CUDA 12, `SAMAYA_CUDA=ON`)

**Why:** this is the "GPU-accelerated" requirement of the problem statement. PLAN.md §5 limits
the GPU to places where the gain is measurable, and PDLP is that place.

**Scope:**
- Device-resident CSR copies of A and Aᵀ.
- Our own SpMV kernel (a CSR-vector kernel with a warp per row; cuSPARSE is allowed only as a
  benchmark reference).
- Fused vector kernels: projection onto bounds, axpby, and the reductions for norms and gaps.
- The restart logic stays on the host, with scalar reductions copied back once per check.
- The same `solve_pdlp` interface, with `PdlpOptions::use_gpu`.

**Milestones:**
- **M1:** SpMV and vector kernels, unit-tested against the CPU on random matrices.
- **M2:** the full PDLP loop on the device, with iterates matching the CPU version to 1e-10 on
  small LPs over 100 iterations.
- **M3:** benchmark CPU against GPU on large LPs from `bench/generate_lps.py --scale 10+` and
  Mittelmann LP instances. Report time, iterations and speedup, including host–device transfer.

**Done when:** at least 5 large LPs show a GPU speedup at equal tolerance (the PLAN.md Phase 3
exit criterion); the CPU fallback is used automatically when no GPU is present; and the CUDA
preset is in CI if a GPU runner is available (otherwise documented as a manual step).

### WP3: Barrier (interior point) method (`src/lp/ipm.{hpp,cpp}`, `src/linalg/cholesky*`, `amd*`)

**Why:** robust on large, degenerate LPs where the simplex stalls, and it is the base for QP
(WP4).

**Algorithm:**
- Homogeneous self-dual embedding (HSD) with Mehrotra predictor–corrector and Gondzio
  centrality correctors.
- Normal equations `A D Aᵀ` with a **supernodal Cholesky** that we write ourselves, with
  regularization for tiny pivots.
- Dense-column handling (split or augmented system).
- Infeasibility and unboundedness detected from τ/κ.

**New linear algebra (ours, no libraries):**
- **AMD ordering.** Approximate minimum degree (Amestoy, Davis, Duff) is a well-documented
  algorithm; write it and test the fill against a naive minimum degree.
- **Symbolic factorization** (elimination tree, column counts, supernodes).
- **Numeric supernodal Cholesky** with dense kernels on the supernodes.

**Milestones:**
- **M1:** dense Cholesky IPM on small LPs, matching the reference.
- **M2:** sparse Cholesky + AMD, tested against the dense version on random SPD matrices.
- **M3:** HSD with predictor–corrector solves Netlib to 1e-8.
- **M4:** `--lp-method barrier` wired in.
- **M5:** **crossover** (primal push, dual push, then simplex cleanup through `Simplex::solve(statuses)`,
  which already accepts a starting basis), giving a basic solution that passes the strict
  verifier. The same crossover serves PDLP (WP1).
  - **Status (2026-10-07):** built without the push phases. `crossover_basis`
    (`src/lp/crossover.cpp`) ranks the n + m variables by distance to their bounds against their
    reduced cost; the m best are basic, and `Simplex::solve(statuses)` repairs and finishes it.
    `--crossover` reports that basic solution; a barrier or PDLP point that fails verification
    goes to the simplex from the same basis either way.
  - **Measured:** on random LPs the simplex needs about 18x fewer iterations than from the slack
    basis. On Netlib, barrier + crossover gives 93/93 verified optima, 73 of them as crossover
    basic solutions. The other 20 are models where the barrier stops at its iteration limit;
    they go to the dual simplex from the slack basis, since a crossover from an unconverged
    iterate was tried and made most of them slower (dfl001 failed).
  - **Still open:** the push phases, if the guess proves too weak on large LPs (W2), and the 20
    barrier non-convergences (an IPM item, not a crossover one).

**Done when:** Netlib 93/93 by barrier plus crossover, and the time comparison with the dual
simplex is reported.

### WP4: Convex QP (`src/qp/`)

**Depends on WP3** (the factorization and the IPM skeleton).

**Scope:**
- The `Model::Q` lower triangle is already read from QPS (`QUADOBJ` and `QMATRIX`).
- A QP IPM on the regularized augmented system `[-(Q + X⁻¹Z + ρI)  Aᵀ; A  δI]` with an LDLᵀ
  factorization (quasi-definite, so the ordering is static).
- Convexity check: a shifted Cholesky of Q, reporting non-convex models cleanly.
- Extend the verifier: the dual condition becomes `Qx + c − Aᵀy = d`.

**Benchmark:** Maros–Mészáros (add a `maros` set to `bench/fetch_instances.sh`).

**Done when:** at least 80% of Maros–Mészáros is solved and verified (the PLAN.md Phase 3 exit
criterion).

### WP5: MRPL case studies and demo (`cases/`, `bench/`)

**Why:** judges care most about problems that look like MRPL's.

**Scope:**
- **Parametric generators** in the style of `bench/generate_lps.py refinery_scheduling`, using
  public figures: MRPL capacity of about 15 MMTPA; typical crude slate and assays (sulfur, API);
  units (CDU, VDU, FCC, hydrocracker, reformer, DHDS); product specifications (BS-VI diesel
  sulfur ≤ 10 ppm, octane, RVP); inventory and tankage; cargo lots; unit on/off; startup costs.
- **Families:**
  - multi-period planning (LP);
  - crude scheduling with cargo arrival windows (MILP);
  - product blending with pooling linearized (LP);
  - a unit-commitment style utility system (MILP).
- **Deliverables:**
  - the models;
  - `bench/results/cases.csv` comparing us with HiGHS;
  - a short narrative per case (decision, objective, what the solver found);
  - one end-to-end demo script (`cases/demo.sh`) that generates, solves, verifies and prints a
    readable plan.

**Done when:** each family has three sizes, every size agrees with HiGHS within the gap, and the
demo script runs in under 2 minutes.

### WP6: Simplex speed (`src/lp/`, `src/linalg/`)

For the simplex owner.

**Scope:**
- Hyper-sparse FTRAN and BTRAN (symbolic reach and sparse right-hand sides).
- A crash basis (triangular).
- Row-wise pricing with a partial list.
- LU update tuning.
- Presolve additions: doubleton equations, implied free columns, dominated columns via implied
  bounds.

**Target:** the Netlib shifted geomean within 1.5× of HiGHS (it is about 2× now); `dfl001`
and `maros-r7` are the slow outliers.

**Done when:** the Netlib and generated LP sets are faster, with no correctness regressions.

### Ownership and conflict map

| Package | Creates | Edits (coordinate first) |
|---|---|---|
| WP1 | `src/lp/pdlp.*`, `tests/test_pdlp.cpp` | `src/core/solver.cpp` (method dispatch), `CMakeLists.txt` |
| WP2 | `src/gpu/*`, `tests/test_gpu.cpp` | `CMakeLists.txt` (CUDA block) |
| WP3 | `src/lp/ipm.*`, `src/lp/crossover.*`, `src/linalg/{cholesky,amd}.*` | `src/core/solver.cpp`, `include/samaya/solver.hpp` |
| WP4 | `src/qp/*` | `src/verify/verify.cpp`, `src/core/solver.cpp` |
| WP5 | `cases/*`, `bench/generate_*.py` | `bench/harness.py` |
| WP6 | — | `src/lp/simplex.*`, `src/linalg/basis_factor.*`, `src/presolve/*` |
| MILP (ongoing) | `src/mip/*` | `src/core/solver.cpp` (`solve_mip_model`) |

`src/core/solver.cpp` is shared by everyone. Keep edits there to small dispatch hunks and rebase
often.

---

## 7. Definition of done (every package)

- [ ] Differential tests against a reference, on randomized families, including infeasible and
      unbounded models where they apply.
- [ ] A planted-bug check shows the tests catch errors (say so in the PR).
- [ ] Every reported result passes the verifier; any tolerance relaxation is explicit and logged.
- [ ] GCC and Clang, debug/release/asan are all green in CI.
- [ ] Benchmark numbers in `docs/results/`, generated by `bench/report.py`, including
      failures.
- [ ] README status table and this file updated.

## 8. References

- Applegate, Díaz, Hinder, Lu, Lubin, O'Donoghue, Schudy. *Practical Large-Scale Linear
  Programming using Primal-Dual Hybrid Gradient.* NeurIPS 2021.
- Lu, Yang. *cuPDLP.jl: A GPU Implementation of Restarted Primal-Dual Hybrid Gradient for Linear
  Programming in Julia.* 2023.
- Mehrotra. *On the Implementation of a Primal-Dual Interior Point Method.* SIAM J. Optim., 1992.
- Gondzio. *Multiple centrality corrections in a primal-dual method for linear programming.*
  1996.
- Xu, Hung, Ye. *A simplified homogeneous and self-dual linear programming algorithm.* 1996.
- Amestoy, Davis, Duff. *An Approximate Minimum Degree Ordering Algorithm.* SIAM J. Matrix
  Anal. Appl., 1996.
- Davis. *Direct Methods for Sparse Linear Systems.* SIAM, 2006 (elimination trees, symbolic
  and supernodal Cholesky).
- Koberstein. *The dual simplex method, techniques for a fast and stable implementation.* PhD
  thesis, 2005 (our simplex follows it).
- Achterberg. *Constraint Integer Programming.* PhD thesis, 2007 (branching, cuts, heuristics).
- Maros, Mészáros. *A repository of convex quadratic programming problems.* 1999.
