# Plan: closing the gaps to the SIH26119 problem statement (from 2026-09-27)

This plan compares SIH26119 ("Indigenous GPU-Accelerated Optimization Solver (Sovereign
Alternative to Xpress / CPLEX)", MRPL) against what samaya can show today, and orders the work
that closes the difference. The MILP search improvements are in `docs/IMPROVEMENT_PLAN.md`; the
work packages are specified in `docs/ARCHITECTURE.md` §6.

The rules of `CLAUDE.md` hold for every step:
- no external solver or linear-algebra code;
- every result goes through `src/verify/`;
- differential tests plus a planted bug;
- the three presets before each commit;
- fixed benchmark sets, chosen before running.

## 1. What the problem statement asks, and where we stand

| Requirement (quoted or paraphrased from SIH26119) | Evidence today | Gap |
|---|---|---|
| LP, MILP **and QP** "as the initial focus" | LP and MILP solved and verified; QPS files are read (`Model::Q`) but QP returns `kNotImplemented` | **QP solver missing** |
| Modular, extensible to MIQP, NLP, MINLP | `ProblemClass` has kQP and kMIQP; one B&B engine | MIQP follows from QP (section 4, step 4) |
| Revised simplex **and interior-point** methods | Dual simplex with our own LU; PDLP (WP1, teammate, nearly done) | **No barrier (interior point)** |
| Branch-and-cut, cuts, presolve, heuristics, node selection | All present, measured item by item | Keep improving (IMPROVEMENT_PLAN) |
| Sparse techniques, multi-core, GPU where it helps | Sparse LU, parallel tree search; GPU PDLP (WP2, teammate, nearly done) | Measured GPU speedups still to publish |
| Benchmarks: MIPLIB, Netlib, **Mittelmann**, QPLIB; compared with an established solver | Netlib 93/93, MIPLIB 62, the cases; against HiGHS, SCIP, CBC, GLPK | **No Mittelmann large LP set, no QP set** |
| Industrial scale: "thousands to millions of variables" | Largest instances solved: MIPLIB tens of thousands of columns | **No million-variable run** |
| "Clear demonstration of numerical robustness": degeneracy, weak relaxations, ill-conditioning | The pieces exist (degenerate Netlib models, 28/29 infeasibility certificates, verifier, re-solve at tighter tolerances), but nowhere as one demonstration | **Robustness report missing** |
| Refinery scheduling, blending, planning, logistics, power dispatch | Three MRPL case families plus generators; faster than HiGHS on them | Blending with QP (pooling-free) once QP exists |

## 2. Timeline this plan is built around

- **Idea screening:** September–October; nomination deadline about 30 September.
- **Results:** October.
- **Grand finale:** December, a 36-hour software sprint.

Two consequences:
- The **idea submission** must state each requirement's status honestly, as built and measured
  or as in progress with a date.
- The **finale** needs everything in sections 4–6 working and demonstrable.

## 3. Phase 0: submission-ready evidence (by 29 September)

1. **Requirement-to-evidence table** in `docs/sih/SIH_PPT_CONTENT.md`: the table of section 1,
   with numbers and links. Slide 4 gets a "Robustness" row.
2. **Robustness report** `docs/results/robustness.md`, built from runs we already have plus one
   fixed set:
   - **Degenerate LPs:** the Netlib models known for degeneracy and cycling risk (degen2,
     degen3, cycle, greenbea, greenbeb, pilot, pilot87, perold, d2q06c). For each, report
     samaya's result against HiGHS, the degenerate-pivot share, and the verifier's residuals.
   - **Infeasibility:** 28 of 29 infeasible Netlib models proven, each with a Farkas certificate
     checked by `verify_infeasibility`.
   - **Weak LP relaxations:** MIPLIB instances with root gaps above 50% (p200x1188c, mc11,
     markshare, neos-911970), showing the root bound before and after our cuts and the time to
     optimality.
   - **Ill-conditioning:** the models on which the verifier triggered the tighter-tolerance
     re-solve, with the condition estimate the LU already computes.
   - **Wrong answers:** zero over every run so far, with the counts.
3. **Today's MILP gains** go into the results and the slides once committed: continuous
   propagation, the conflict skip, and RENS in the tree if its A/B passes.
4. **WP1/WP2 (teammate):** one slide line with the measured GPU speedup and the instance sizes,
   taken from the teammate's run.

Effort: about 4 h of writing plus reruns already scheduled. No new solver code.

## 4. Phase 1: interior point and convex QP (about two weeks, October)

This closes the two largest gaps with one piece of linear algebra.

**Design decision: one augmented-system interior-point method for LP and QP.**
`ARCHITECTURE.md` WP3 plans normal equations (A D Aᵀ with Cholesky) for LP and a separate
augmented system for QP (WP4). The plan is to build only the **regularized quasi-definite
augmented system** and use it for both:

  [ -(Q + Θ⁻¹ + ρI)   Aᵀ ]
  [        A          δI ]

- **Quasi-definite matrices** (Vanderbei 1995) have an LDLᵀ factorization for any symmetric
  permutation. The ordering is therefore fixed once (AMD) and reused every iteration, with no
  pivoting.
- **Primal and dual regularization** (ρ, δ ≈ 1e-8 relative; Friedlander and Orban 2012) keep it
  quasi-definite with rank-deficient A or free variables.
- **No dense-column problem:** normal equations make A D Aᵀ dense when one column of A is dense;
  the augmented system never forms it.
- **LP is the case Q = 0.** One code path serves `--lp-method barrier` and QP.

Clarabel (an interior-point method) and OSQP solve their systems the same way: a quasi-definite
KKT matrix, AMD ordering and QDLDL. That is why the change from WP3's text is justified. `ARCHITECTURE.md` WP3/WP4 is updated with it.

### Steps

**1. Sparse LDLᵀ with AMD ordering** (`src/linalg/ldl.{hpp,cpp}`, `src/linalg/amd.{hpp,cpp}`),
3 days.
- **AMD:** approximate minimum degree on the pattern of K + Kᵀ (Amestoy, Davis and Duff 1996),
  with quotient graph and element absorption.
- **Symbolic:** elimination tree and column counts, done once per pattern.
- **Numeric:** up-looking LDLᵀ (as in QDLDL: row by row along the elimination tree), with
  static regularization. A pivot of the wrong sign or too small is replaced by ±δ, counted and
  reported.
- **Tests:**
  - random quasi-definite matrices against a dense LDLᵀ in `tests/dense_reference.hpp`
    (residual ‖Kx − b‖ ≤ 1e-10 ‖b‖);
  - AMD fill checked against natural order and a naive minimum degree on grid and random
    patterns.
- **Planted bugs:** a skipped update in the numeric phase, and a wrong elimination-tree parent.
  The residual test must fail.

**2. Interior-point method** (`src/lp/ipm.{hpp,cpp}`), 4 days.
- **Method:** Mehrotra predictor–corrector on the bound-constrained form: all our row and column
  bounds, ranged rows as slack columns, free variables kept free thanks to the regularization.
- **Step:** a fraction-to-boundary step of 0.995.
- **Stopping:** relative primal and dual residuals and gap ≤ 1e-8.
- **Infeasibility:** detected from diverging iterates (the homogeneous embedding is a later
  refinement if needed). Every "infeasible" answer must come with a ray that
  `verify_infeasibility` accepts, as the simplex's answers do.
- **Scaling:** Ruiz equilibration (the same as WP1's), then the existing `Scaling` on top.
- **Tests:**
  - differential against `tests/reference_lp.hpp` on the random LP families (feasible,
    infeasible, unbounded, degenerate, boxed);
  - Netlib at 1e-8.
- **Planted bug:** a sign error in the corrector right-hand side (the tests must see wrong optima
  or non-convergence).

**3. Convex QP** (`src/qp/qp.{hpp,cpp}`, dispatch hunk in `src/core/solver.cpp`), 3 days.
- **The method:** the same interior-point method with Q in the (1,1) block. `Model::Q` is the
  lower triangle, and the matrix assembly reads it symmetrically.
- **Convexity check:** an LDLᵀ of Q + εI. A negative pivot beyond the tolerance means
  non-convex, reported cleanly as `kNotConvex`, never solved wrongly.
- **Verifier:** new `verify_qp_optimality` in `src/verify/`: the reduced cost is d = Qx + c − Aᵀy,
  with the same sign conventions as LP (ARCHITECTURE §4.2). Tolerances are named constants,
  documented, and never loosened.
- **Reference for tests:** `tests/reference_qp.hpp`, a small dense active-set QP solver (exact
  on tiny problems), on random families:
  - convex with Q ⪰ 0 of low rank;
  - strictly convex;
  - with equality and ranged rows;
  - infeasible;
  - non-convex, which must be refused.
- **Planted bugs:** Q applied without its upper triangle (the lower-triangle convention), and
  a wrong sign of Qx in the dual residual. The verifier and the reference comparison must catch
  both.

**4. MIQP**, 2 days. Branch-and-bound over QP relaxations, reusing `BranchAndBound`'s tree
without LP-specific parts:
- no cuts at first;
- no dual proofs, since conflict analysis assumes LP duals;
- the relaxation is warm-started by the interior-point method from the parent's point.

It stays small and is honestly labelled as the first step of the "modular extension" the
problem statement asks for.

**5. Crossover for LP**, 3 days, shared with PDLP (WP3 M5).
- **Method:** a primal and dual push to a basis, then a simplex cleanup through
  `Simplex::solve(statuses)`.
- **Why:** it gives interior-point and PDLP answers a basic, strictly verified solution. This is
  what lets the teammate's GPU PDLP report strict-tolerance optima.

**Benchmarks:**
- **Maros–Mészáros convex QP set:** 138 models, a `maros` set in
  `bench/fetch_instances.sh`, compared against HiGHS's QP solver in the harness.
  - **Target** (WP4's exit criterion): at least 80% solved and verified.
- **Netlib by barrier plus crossover:** 93/93 target, with times against the dual simplex.
- **QPLIB:** the convex continuous subset, as the problem statement names it.

**Risks:**
- **LDLᵀ numerics on ill-conditioned Maros–Mészáros models:** the regularization plus iterative
  refinement (two steps) is the standard remedy, and the verifier guards the answers.
- **Time:** if the interior-point method runs late, the QP story still needs step 3, so steps 1–3
  come before 4–5.

## 5. Phase 2: scale and GPU (teammate's WP1/WP2 plus us, late October)

1. **Merge WP1 and WP2** when the teammate is done:
   - run the three presets and the CUDA preset;
   - check that `--lp-method pdlp` answers go through the verifier, polished by the simplex or
     the crossover (step 5 above).
2. **Mittelmann large LP set:** a fixed subset chosen now from Mittelmann's "LPopt" and "feasible
   LP" benchmarks, of sizes 10⁵ to 10⁷ nonzeros (for example qap15, nug08-3rd, rail4284,
   dbic1, ns1688926, stormG2_1000, L1_sixm250obs, zib03).
   - **Solvers:** dual simplex, barrier, PDLP on the CPU and PDLP on the GPU, against HiGHS.
   - **Report:** time, iterations and verified status.
   - **Hardware:** the laptop's RTX 3060 and, if available, a data-centre GPU.
3. **Million-variable demonstration:** `bench/generate_lps.py refinery_scheduling --scale` at
   10⁶+ columns, and a multi-period refinery planning model from the case generator at full
   year scale. Solve with GPU PDLP plus crossover, then verify. This answers "millions of
   variables" directly.
4. **Parallel speedup table:** the tree search on 1, 4 and 8 threads over the MIPLIB set.
   The code exists; the measurement doesn't.

## 6. Phase 3: MILP competitiveness and the finale (November)

- **From `docs/IMPROVEMENT_PLAN.md`:**
  - the presolve aggregator;
  - the clique table in propagation and probing;
  - zero-half cuts;
  - symmetry, if time allows;
  - WP6 simplex speed (hyper-sparse FTRAN/BTRAN), since the profile shows the LP is where the
    time goes.
- **Final results:** netlib, cases, MIPLIB 60 s and 600 s, Maros–Mészáros and Mittelmann, all on
  the same machine with the clock checked, written into `docs/results/`.
- **Finale kit:**
  - `cases/demo.sh` extended with a QP blending case and a GPU LP case;
  - a 5-minute scripted demo;
  - the robustness report as the backup slide;
  - a list of features safe to build during the 36-hour sprint (small, testable, visible: e.g.
    MIQP cuts, a solution-pool option, a sensitivity report for the planning LP).

## 7. Order and gates

| # | Item | Effort | Gate |
|---|---|---|---|
| 0 | Evidence table, robustness report, slides | 4 h | every number traceable to a run in `docs/results/` |
| 1 | Sparse LDLᵀ + AMD | 3 days | residual tests on random quasi-definite matrices; planted bugs caught |
| 2 | Interior point for LP | 4 days | random families match the reference; Netlib solved at 1e-8 and verified |
| 3 | Convex QP + `verify_qp_optimality` | 3 days | reference QP tests; non-convex refused; ≥ 80% of Maros–Mészáros |
| 4 | MIQP (B&B over QP) | 2 days | random MIQP match a brute-force reference |
| 5 | Crossover (LP) | 3 days | barrier and PDLP points become verified basic optima on Netlib |
| 6 | Merge WP1/WP2; Mittelmann set; million-variable run | 3 days (with teammate) | verified results, speedups reported with transfer time |
| 7 | MILP items from IMPROVEMENT_PLAN, WP6 simplex speed | rest of November | A/B gates as in that plan |
| 8 | Final results and finale kit | 3 days | the same machine, the clock checked, everything reproducible |

- **Critical path:** 1 → 2 → 3 (QP by about 12 October if started on 28 September), then 5 → 6
  (large LPs with GPU PDLP).
- **Items 4 and 7** can move to fill gaps.
- **Owners:** 1–5, 7 and 8 are ours; 6 is shared with the teammate (WP1/WP2).
