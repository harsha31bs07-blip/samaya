# Plan after the idea submission: everything left out at the deadline (from 2026-10-07)

The SIH 2026 idea submission (PS 26119) went in on 29 September. This plan collects every item
that was dropped, postponed or only half measured to meet that deadline, and orders them up to
the grand finale (December, a 36-hour sprint).

Sources:
- `docs/IMPROVEMENT_PLAN.md`: the MILP items, with what was tried and not merged.
- `docs/SIH_GAP_PLAN.md`: the gaps to the problem statement, phases 2 and 3.
- `docs/handover/HANDOVER.md`: the open items of 27–28 September.
- `docs/WINDOWS_PLAN.md` and `docs/STUDIO_REDESIGN_PLAN.md`: what was not built.
- The judge-style review of 28 September: MILP strength, scale, the GPU claim, synthetic cases,
  and fit with how refineries plan (recursion/SLP).

The rules in `CLAUDE.md` hold for every item:
- a differential test against the reference solvers;
- a planted bug that the test catches, mentioned in the commit message;
- the debug, release and asan presets with zero warnings (and the CUDA preset for GPU code);
- every reported number through `src/verify/`, with tolerances never loosened.

Benchmarks run on one machine with the clock checked (`bench/clock_factor.py`), and no script is
edited while a benchmark is running it.

## 1. Where things stand (measured, laptop)

| Area | Today | Gap |
|---|---|---|
| LP (Netlib, 93) | 93/93, geomean 0.36 s | HiGHS 0.19 s; infeasibility proofs 0.10 s vs 0.008 s; `cplex2` unproven |
| Refinery-type cases (16) | 16/16, geomean 3.43 s (fastest of the five solvers on the laptop) | All generated, all small (largest 3,328 × 4,732); `refsched_8c_6p_3u_26t` 122.6 s vs HiGHS 47.0 s |
| MIPLIB 2017 (62 easy) | 60 s: 12–13; 600 s: 17 | HiGHS 18 / 32, SCIP 19 / 25 |
| Convex QP (Maros–Mészáros) | 106/138, HiGHS 104 | LISWET stalls; large CVXQP*_L; degenerate LP-like QPs need a crossover |
| Large LP and GPU | PDLP on the CPU and GPU works, with a CPU fallback | No crossover; the 467,376 × 1,001,520 refinery LP was unconverged at 600 s on the RTX 3050; no Mittelmann results; `docs/results/large_lp.md`, cited by `bench/gpu_lp.test`, was never written |
| Parallel tree search | Works | Only a 60 s screening on an old build (1 thread: 5 solved, 6 threads: 6); no speedup table |
| Refinery fit | LP/MILP/QP from MPS | No recursion (SLP) for blend qualities and pooling; no sensitivity ranging; Studio lacks tank capacity lines |
| Windows | Release zip, `setup.bat`, Studio, GPU option | No one-file user installer; unsigned binaries; no `.mps.gz` |

## 2. Workstreams

Effort is in working days, including tests, the presets and the A/B run.

### W0. Housekeeping (0.5 day, first)

1. The empty benchmark files in the repo root (`cases.txt`, `netlib.txt`, `miplib600-t6.txt`)
   and the two MIPLIB summaries (`miplib60.txt`, `miplib600.txt`):
   - move the summaries into `bench/results/` with the build and date they came from;
   - delete the empty files.
2. The local `claude/*-wip` branches (estimate, rr, tree-rens, tree-rens2, conflict-speed and
   others):
   - list each with its A/B verdict in `docs/IMPROVEMENT_PLAN.md`, which mostly exists already;
   - keep the branches and touch no `results/*` branch.
3. Slide and video fixes for the finale deck:
   - "Express / CEPLEX" on slide 1 becomes "Xpress / CPLEX";
   - the YouTube link points to the Prabhat-voice cut; re-upload the own-voice cut if wanted.

### W1. Crossover for barrier and PDLP (3 days; unlocks W2 and part of W5)

Barrier and PDLP end at an interior or approximate point; the verifier then rejects some answers
(degenerate LP-like QPs) and there is no basic solution to warm-start from or to report as
"optimal basis".

- **Method:** primal and dual push phases to a basis (Megiddo; as in HiGHS and Clp), then a few
  simplex iterations to clean up with our own LU.
- **Gate:**
  - Netlib by barrier plus crossover: 93/93, verified basic optima, times next to the dual
    simplex;
  - PDLP plus crossover on the same set at 1e-4, then polished to verified optima.
- **Tests:**
  - random LP families: crossover from the reference optimum perturbed into the interior must
    return a basis that `verify_lp_optimality` accepts;
  - planted bugs: a wrong sign in a dual push step, a skipped bound flip.
- **Also:** re-run Maros–Mészáros; the LP-like QPs that fail strict verification today should
  pass. Update `docs/results/qp.md`.

### W2. Large LPs and the GPU (4 days; answers "millions of variables")

1. **GPU-side residual evaluation for PDLP** (handover open item). The host evaluation now
   copies the iterate back every check, which limits the GPU. Compute the primal and dual
   residuals and the gap with GPU kernels, and copy only scalars.
   - **Test:** CPU and GPU residuals agree to 1e-12 relative on random LPs.
2. **The fixed Mittelmann subset** (`bench/mittelmann_lp.test`, fixed on 27 September, before
   any run) and the large refinery LP (`bench/gpu_lp.test`).
   - **Runs:** dual simplex, barrier, PDLP on the CPU, PDLP on the GPU, against HiGHS.
   - **Report:** time, iterations, verified status, and the host–device transfer time stated
     separately.
3. **The 467,376 × 1,001,520 refinery LP:**
   - re-run with the residuals on the GPU and the 687448b restarts, at 1,800 s;
   - then crossover (W1) and the verifier.
   - If it still does not converge, say so in the results file.
4. **Write `docs/results/large_lp.md`.** Every number comes from these runs.
   - **State plainly where the GPU wins and where it doesn't.** On the small models the CPU is
     faster today.

**Gate:** at least one verified solve above 10⁶ columns, by whichever method is fastest. The
slides then quote that run, and nothing else, for scale.

### W3. MILP competitiveness (IMPROVEMENT_PLAN items C–G; about 2.5 weeks)

The order, the targets and the merge gates are in `docs/IMPROVEMENT_PLAN.md`, sections 2–4. In
short:

| # | Item | Effort | Instances it should move |
|---|---|---|---|
| 0 | Step 0: attribution runs (switch components off in HiGHS and SCIP on the 18 lost instances) | 0.5 day, machine | the table that every later item cites |
| C | Presolve aggregator: implied-free columns, doubleton equations, general substitution with a fill limit | 2–3 days | neos-3381206-awhea, n5-3, beasleyC3, app1-1 |
| D | Clique table in propagation and probing; implied-bound cuts | 1.5 days | csched008, gmu-35-40, mcsched, graph20-20-1rand |
| E | Inference history in the branching score | 0.5 day | A/B only |
| F | Zero-half cuts | 1 day | enlight_hard, timtab1 |
| G | Symmetry: our own graph refinement and orbital fixing | 3–4 days | fhnw-binpack4-4, graph20-20-1rand, fastxgemm |
| A4 | Shifting and zi-round heuristics | 1 day | the instances with no solution at 60 s |

- **Item B** (propagation on continuous columns) is done (9343abf).
- **Not merged** (A1 estimate selection, A2 randomized rounding, A3 RENS in the tree) stays
  documented and is not retried without a new idea.
- **After C and after D:** re-run MIPLIB at 60 s and 600 s and the cases, and update
  `docs/results/comparison.md`.
- **Aim at 600 s:** from 17 to 23–25 of 62, level with SCIP. Matching HiGHS (32) needs G and W4
  as well.

### W4. Simplex speed (WP6; 1.5 weeks, profile first)

The LP is where MILP time goes (gen-ip054, pg5_34), and Netlib runs at about 2× HiGHS.

1. **Profile** Netlib and the slowest MIPLIB node LPs. Record where time goes (FTRAN, BTRAN,
   pricing, ratio test, refactorization) before changing anything.
2. **Hyper-sparse FTRAN and BTRAN:**
   - a symbolic pass (depth-first search on the column graph) when the right-hand side is
     sparse;
   - hyper-sparse CHUZC for pricing.
3. **Faster infeasibility proofs:** take the dual ray as soon as phase 1 stalls with a
   certificate candidate. The target is the 0.10 s vs 0.008 s gap, and `cplex2` (infeasible by
   less than the tolerance) stays honestly "unproven".
4. **Barrier ordering (found 8 October):** `minimum_degree_order` (`src/linalg/ordering.cpp`)
   updates exact external degrees, which takes about 15 s on a 400 x 400 transportation LP
   (160,000 columns) before the first barrier iteration, and the barrier's setup does not check
   the time limit. Replace the degree update with AMD's approximate degrees (Amestoy, Davis and
   Duff 1996), and check the limit between setup steps. Gate: the same fill on the random
   quasi-definite tests, the ordering under 1 s on that LP.
5. **Gate:**
   - Netlib geomean from 0.36 s to 0.25 s or better, with no lost instance and every answer
     verified;
   - the node rate on the MIPLIB set up, and the bound A/B no worse.

### W5. Refinery fit: what MRPL's planners actually run (2 weeks)

Refinery planning tools (Aspen PIMS and similar) solve blend-quality and pooling models by
successive linear programming ("distributive recursion"). samaya has LP, MILP and convex QP, but
no recursion loop.

1. **SLP / recursion driver** on top of the LP engine (`src/slp/`, a new work package):
   - **Loop:**
     - linearize the bilinear quality terms (flow × quality) at the current point;
     - solve the LP;
     - update the qualities;
     - stop when the qualities and the objective stop changing.
   - **Safeguards:** a step bound (trust region) and a stop on oscillation.
   - **Warm start:** each LP warm-starts from the previous basis, which needs W1 only for
     barrier starts.
   - **Verification:** every LP in the loop is verified, and the final point is checked against
     the nonlinear quality constraints, with the residual reported.
   - **Tests:** the Haverly pooling problems (known optima, including the local optimum trap);
     a generated multi-crude blending case.
   - **Planted bug:** a wrong-sign linearization term.
2. **QP blending case:** the "blending with QP" row of `SIH_GAP_PLAN.md` §1, a pooling-free
   blend with quadratic quality penalties, added to `cases/mrpl.py` and `cases/demo.sh`.
3. **Sensitivity report (ranging)** for the planning LP:
   - objective-coefficient and right-hand-side ranges from the optimal basis;
   - in the CLI (`--ranging`), the JSON, Studio's report and the Excel export.

   Planners read these ("how much can Arab Light's price move before the slate changes").
   - **Test:** finite differences on random LPs.
4. **Studio leftovers** from `STUDIO_REDESIGN_PLAN.md`:
   - tank capacity lines on the stock chart;
   - utilization against capacity in the plan.

   Both need column bounds next to the solution: add them to the JSON result.
5. **A realistic, public benchmark:**
   - scale `cases/mrpl.py` to full-year planning and a 4-week, 20-cargo schedule;
   - add a published refinery scheduling instance if one with a licence for reuse is found.

   The README keeps calling the cases generated.

### W6. Parallel speedup (2 days, mostly machine time)

- MIPLIB 62 at 600 s with 1, 4 and 8 threads, on the current build, alone on the machine.
- Write the table in `docs/results/comparison.md` (solved, shifted geomean, speedup on the
  instances all three solve).
- The 60 s screening of September (old build: 5 → 6 solved) is the baseline to beat.
- If the speedup is poor, profile lock contention in `src/mip/parallel.cpp` before tuning.

### W7. Windows distribution (3 days)

1. **One-file user installer** (`samaya-setup.exe`, WINDOWS_PLAN W6):
   - per-user install with Start-menu and desktop shortcuts;
   - an optional "Open .mps with samaya Studio" association;
   - the WebView2 bootstrapper for old Windows 10;
   - uninstall from Settings.
2. **Code signing** (WINDOWS_PLAN decision 4): the owner decides on a certificate before a public
   release, since unsigned binaries trigger SmartScreen.
3. **`.mps.gz` input** with our own inflate (MIPLIB ships gzip). It needs its own tests against
   zlib-compressed fixtures, generated in the test, not shipped.
4. **GPU on Windows:** measure the WDDM launch overhead and timeout detection on the RTX 3050,
   and record it in `docs/results/windows.md`.
5. **Release v0.2.0** with W1–W5, built by `tools/windows/package.ps1`. The asset is replaced
   only after the owner's OK.

### W8. Finale preparation (from 1 December)

1. **Freeze features** a week before the finale; only fixes after that.
2. **Final results:**
   - Netlib, cases, MIPLIB 60 s and 600 s, Maros–Mészáros, Mittelmann and the large LP;
   - one machine, the clock checked, written into `docs/results/`, every number traceable to a
     run.
3. **Finale kit:**
   - `cases/demo.sh` with the QP blending case, an SLP pooling case and a GPU LP case;
   - a 5-minute scripted demo in Studio;
   - the robustness report as the backup slide;
   - the deck and the video updated with the new numbers.
4. **"Defend the internals" sheets:** one page per component (sparse LU, dual simplex, barrier,
   PDLP, branch-and-cut, verifier, SLP). For each: what it does, why that design, and the
   questions a judge will ask. Every team member presents each one once before the finale.
5. **List of features safe to build in the 36-hour sprint:** small, testable and visible. For
   example: MIQP cuts, a solution pool, extra Studio views, or a ranging export.

## 3. Calendar (IST)

| Weeks | Dates | Work |
|---|---|---|
| 1 | 7–13 Oct | W0; W3 step 0 (machine time, overnight); W1 crossover starts |
| 2 | 14–20 Oct | W1 done; W2 GPU residuals; W2 runs overnight |
| 3–4 | 21 Oct – 3 Nov | W2 report; W3 C presolve aggregator, then D clique table; MIPLIB re-run |
| 5–6 | 4–17 Nov | W5 SLP, QP blending case, ranging, Studio leftovers; W3 E and F |
| 7–8 | 18–30 Nov | W4 simplex speed; W3 G symmetry if the gates allow; W6 parallel runs; W7 installer, v0.2.0 |
| 9+ | from 1 Dec | W8 freeze, final results, finale kit |

**Critical path:** W1 → W2, which answers scale and the GPU claim. W3 C and D are the biggest
MILP gains.

- If the October results are delayed, the order stays the same; everything here improves the
  solver either way.
- If the team is not shortlisted, W3, W4 and W5 remain the useful core.

## 4. Gates summary

| Item | Merge gate |
|---|---|
| W1 crossover | Netlib 93/93 by barrier + crossover, verified; random-family tests; planted bugs caught |
| W2 large LP | At least one verified solve above 10⁶ columns; `large_lp.md` with transfer times |
| W3 C–G | The A/B gates of IMPROVEMENT_PLAN §4; no wrong answer anywhere |
| W4 simplex | Netlib ≤ 0.25 s geomean, nothing lost, node rate up |
| W5 SLP | Haverly optima reproduced; final point checked against the nonlinear constraints |
| W5 ranging | Finite-difference test on random LPs |
| W6 parallel | Table on the current build, alone on the machine |
| W7 release | Clean install on a fresh Windows user; the owner's OK before replacing the asset |

## 5. Decisions for the owner

1. Whether to re-upload the own-voice video, and update the YouTube link in the README and the
   slides.
2. A code-signing certificate (cost) before a public v0.2.0.
3. Whether SLP (W5.1) comes before simplex speed (W4). This plan puts SLP first because it is
   the refinery-fit question judges will ask; W4 helps MILP more.
4. Whether to drop G (symmetry) if C and D land late.
