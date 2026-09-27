# Fast plan: closing the SIH26119 gaps (from Sun 27 Sep 2026, 12:24 IST)

This compresses `docs/SIH_GAP_PLAN.md` for the idea nomination deadline (about 30 Sep).
- **Target:** done by **Mon 28 Sep, 18:00 IST** (about 30 h).
- **Hard limit:** **Tue 29 Sep, 12:24 IST** (48 h); the hours after the target are buffer, not
  planned work.

WP1 (PDLP) and WP2 (GPU PDLP) are merged by the teammate this evening.

The rules of `CLAUDE.md` hold throughout:
- tests with planted bugs;
- the three presets before every commit;
- the verifier on every result;
- fixed benchmark sets;
- nothing pushed without the owner's OK.

## Scope

**Must:**
1. **Convex QP** (interior point on a regularized quasi-definite system, our own sparse LDLᵀ),
   with `verify_qp_optimality`, tests against a reference QP solver, clean refusal of
   non-convex models, and a Maros–Mészáros run against HiGHS.
2. **Interior point for LP** (`--lp-method barrier`): the same code with Q = 0, verified on
   Netlib.
3. **Robustness report:** degenerate, infeasible, weak-relaxation and ill-conditioned models,
   as the problem statement asks.
4. **WP1/WP2 merged and verified;** a Mittelmann large-LP subset and a million-variable
   generated refinery LP with the GPU speedup measured.
5. **Slide content** with a requirement-to-evidence table, every number traced to a run.

**Out:**
- MIQP;
- crossover (PDLP and barrier points are polished by the simplex, as WP1 M4 specifies);
- AMD's approximate degrees (plain minimum degree first);
- further MILP work.

## How the time is saved

- **One code path for LP and QP:** the barrier is the QP solver with Q = 0, so there is no
  separate LP interior-point phase.
- **Coding overlaps timing runs.** Writing code never disturbs a benchmark; only building and
  testing do, so those wait for the gaps between runs.
- **The heavy benchmarks run overnight and do not need the QP code:** the robustness data, Netlib
  by PDLP, Mittelmann and the million-variable LP. That frees the day for QP.
- **One set of presets per commit, all three in parallel** (about 55 min, bounded by asan).
  Commits are batched so that presets run twice, not five times.

## Schedule (IST)

| # | Target | Latest | Work |
|---|---|---|---|
| 1 | Sun 16:00 | Sun 18:00 | **MILP work in flight.** A/Bs for the conflict skip and RENS in the tree (after the presets now running); commit continuous propagation plus whichever pass. **While they run, write:** `src/linalg/ldl.{hpp,cpp}` (elimination tree, column counts, up-looking LDLᵀ with static regularization, the QDLDL scheme) and `src/linalg/ordering.{hpp,cpp}` (minimum degree on the quotient graph, exact external degrees). |
| 2 | Sun 20:00 | Sun 23:00 | **WP1/WP2 merge:** debug, release, asan and CUDA presets; Netlib with `--lp-method pdlp` at 1e-4 and 1e-8, every answer through the verifier; fix integration faults. **In between builds:** `tests/test_ldl.cpp` (random quasi-definite matrices against a dense LDLᵀ, residual ≤ 1e-10; fill against natural order), with planted bugs (a skipped update; a wrong elimination-tree parent). |
| 3 | Mon 02:00 | Mon 06:00 | **Interior point and QP together** (`src/lp/ipm.*`, `src/qp/qp.*`, dispatch hunks in `src/core/solver.cpp`), detailed below. |
| 4 | Mon 02:00–08:00 | same | **Overnight timing runs, no builds:** (a) the robustness data: Netlib's degenerate models (degen2, degen3, cycle, greenbea, greenbeb, pilot, pilot87, perold, d2q06c) on samaya and HiGHS with pivot statistics; (b) Netlib by PDLP on the CPU and GPU; (c) the Mittelmann subset, fixed now (qap15, nug08-3rd, rail4284, dbic1, ns1688926, stormG2_1000, L1_sixm250obs, zib03), 1,000 s each, with dual simplex, CPU PDLP, GPU PDLP and HiGHS; (d) a generated refinery-scheduling LP at 10⁶+ columns by GPU PDLP, verified. The clock is checked before and after. |
| 5 | Mon 11:00 | Mon 16:00 | **Presets for the QP/IPM commit** (all three in parallel), fixes, commit. Netlib by barrier (fast: minutes). |
| 6 | Mon 14:00 | Mon 20:00 | **Maros–Mészáros:** 138 convex QPs, samaya against HiGHS's QP solver, 300 s, 4 at a time (most finish in seconds). |
| 7 | Mon 18:00 | Tue 12:24 | **Results and submission material:** `docs/results/{qp,robustness,large_lp}.md` (every instance, failures included); `SIH_PPT_CONTENT.md` (the evidence table, QP and interior-point rows, GPU speedup with transfer time, robustness row, today's MILP numbers); `HANDOVER.md`. The branch goes to the owner to review and push. |

## Step 3 in detail: the interior-point method and QP (one code path)

**Linear system:** `[-(Q + Θ⁻¹ + ρI)  Aᵀ; A  δI]`, with ρ and δ at 1e-8 relative.
- It is quasi-definite (Vanderbei 1995), so the ordering is computed once and reused.
- Two steps of iterative refinement per solve.

**Algorithm:** Mehrotra predictor–corrector on our bound form (row and column bounds; ranged
rows as slack columns; free variables stay free thanks to the regularization).
- Ruiz equilibration.
- Fraction to boundary 0.995.
- Stopping when the relative primal and dual residuals and gap are ≤ 1e-8.

**QP handling:**
- `Model::Q` holds the lower triangle and is applied symmetrically.
- **Convexity:** an LDLᵀ of Q + εI. A negative pivot means `kNotConvex`, reported cleanly, never
  solved.

**Infeasibility:** detected from diverging iterates, and reported only with a ray that
`verify_infeasibility` accepts.

**Verifier:** `verify_qp_optimality` in `src/verify/`: d = Qx + c − Aᵀy, with the LP sign
conventions (ARCHITECTURE §4.2) and named tolerances.
- **LP by barrier:** accepted only if the strict verifier accepts it, else polished by the dual
  simplex from a basis guessed from the point (as WP1 M4). Never loosened.

**Tests:**
- the random LP families against `tests/reference_lp.hpp`;
- random QPs against a new `tests/reference_qp.hpp`, a dense active-set solver for tiny
  problems. Families: strictly convex, low-rank PSD, equalities and ranges, infeasible, and
  non-convex (must be refused).

**Planted bugs:**
- a corrector sign error;
- Q used as its lower triangle only;
- a wrong sign of Qx in the dual residual.

## Checkpoints and fallbacks

| Checkpoint | If it slips |
|---|---|
| Sun 23:00: LDLᵀ passes its residual tests | Natural order plus reverse Cuthill–McKee instead of minimum degree: more fill, still correct |
| Mon 06:00: IPM matches the references on LP and QP families | QP only, on small and medium models; no barrier-for-LP claim |
| Mon 16:00: presets green | Nothing uncommitted is claimed; the slides call QP "in progress", with the test results so far |
| Mon 20:00: benchmarks finished | Report what finished, the rest marked unfinished; never extrapolate |
| Tue 12:24: hard limit | Whatever is green and measured goes into the slides; nothing else |

## Needs from the owner

- **Push OK before the merge (about 17:00 today):** `claude/kind-hamilton-jr2b6u` has 9 local
  commits, plus today's. Alternatively, name the branch the teammate merges into, and it is
  merged locally.
- **The laptop stays on,** plugged in and in performance mode, from tonight to Monday evening.
- **If the usage limit interrupts:** work resumes from the last commit or saved WIP branch. Each
  row of the schedule ends in one.

## Progress (updated Sun 27 Sep, 13:25 IST)

| Item | State |
|---|---|
| Continuous-column propagation (MILP) | committed 9343abf: bound A/B mean gap 19.99% → 17.78%, better 21 / worse 8 |
| Conflict-pool skip, RENS in the tree (MILP) | presets green; A/Bs this afternoon, then committed if their gates pass |
| Sparse LDLᵀ + minimum-degree ordering | done, in the QP commit |
| Interior point for LP and convex QP | done: 450/450 random QPs match the reference; Maros–Mészáros screen 105/138 verified, 0 wrong; Netlib by barrier 71/93 verified from the interior point alone (93/93 with the simplex); presets running |
| Robustness report | `docs/results/robustness.md` written |
| Slide content | evidence table added to `docs/sih/SIH_PPT_CONTENT.md`; timed numbers pending |
| Mittelmann subset | fixed in `bench/mittelmann_lp.test` (13 instances), downloading |
| WP1/WP2 merge | waiting for the teammate (this evening) |

Tried and not kept today: best-estimate node selection (A/B worse 14 / better 7), randomized
rounding at every node (solutions 54 → 52), incremental row activities (2–9% slower: the profile
shows the LP, not propagation, is where the time goes), Gondzio correctors (no gain on LISWET,
lost STADAT1).
