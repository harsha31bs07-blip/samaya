# MILP implementation plan (our part)

The teammate owns WP1 and WP2 (PDLP on the CPU, then on the GPU). This file covers the rest of
our work: making the MILP solver competitive (Track A), then the MRPL case studies (Track B).
It also says where each benchmark runs: in the cloud session or on the Dell G15.

## Where we start

The baseline is `main` at ba0046e. Full numbers are in [results/miplib.md](results/miplib.md).

| Measure (62 small "easy" MIPLIB 2017 instances) | Now |
|---|---|
| Solved in 10 min (first cut version) | 8/62 |
| Solved in 60 s | 3/62 |
| Feasible solution found in 60 s | 37/62 |
| Median gap to the known optimum at 60 s | 2.8% |
| Wrong answers | 0 |

**Why instances fail today:**
- **No feasible solution.** On 25/62 instances no feasible solution is found in 60 s, for
  example b1c1s1, cost266-UUE, glass4, lotsize, timtab1 and the rococo instances.
  - Cause: the only heuristics are simple rounding and round-and-solve.
- **Weak bounds.** Many instances have a good solution but a bound that doesn't move: mas74,
  p200x1188c, gen-ip*, pg5_34, neos5.
  - Cause: no flow covers, no probing or clique information, and cuts only at the root.
- **One core only.** All searches are single-threaded; the laptop has 16 hardware threads.

## Targets (after the deadline)

| Measure | Now | Target |
|---|---|---|
| Solved in 10 min | 8/62 | **≥ 20/62** |
| Feasible solution in 60 s | 37/62 | **≥ 52/62** |
| Median gap at 60 s | 2.8% | **≤ 1%** |
| Wrong answers | 0 | **0**, always |

All targets are for a single thread. The parallel search (A6) is reported separately.

## 48-hour schedule (deadline: 2026-09-27)

The full plan below is more than 48 hours of work, so it is cut to what shows best to the
judges. Code freezes at **T+36 h**; the last 12 hours are final benchmarks, docs and the demo.

| Hours | Work | Benchmark (ask: laptop or cloud) |
|---|---|---|
| 0–2 | A1: exact time limit | A0 10-min baseline on the laptop, in parallel |
| 2–12 | A3: diving, feasibility pump, RENS (sub-MIP), shift-and-propagate | 60 s screening |
| 12–20 | A2 (lite): probing on binaries, coefficient tightening | 60 s screening |
| 20–28 | Track B: two MRPL case families and `cases/demo.sh` | HiGHS cross-check |
| 28–36 | A4 (lite): flow covers, clique cuts from probing | 60 s screening |
| 36–40 | Final 10-min run, samaya and HiGHS, same machine | laptop |
| 40–48 | Results docs, merge the teammate's PDLP/GPU work, demo rehearsal, buffer | – |

**Dropped for the deadline:** A5 (conflict analysis), cuts in the tree, zero-half cuts, root
restarts, barrier and QP. (A6, the parallel search, was dropped at first and then done because
the other steps finished early.) If a step runs over, the next one is cut,
never the tests or the verifier.

**Targets for 48 hours (single thread):**
- solved in 10 minutes: at least 12/62;
- a feasible solution within 60 s: at least 48/62;
- median gap at 60 s: at most 1.5%;
- no wrong answers.

The longer-term targets below still stand for after the deadline.

## Progress (updated as work lands)

60 s screening of the 62 instances, single-threaded, four at a time on the 4-core cloud machine.
"Median gap" counts an instance without a solution as 100%; that is stricter than the 2.8%
above, which ignores such instances.

| Build | Feasible | Median gap | Solved | Wrong |
|---|---|---|---|---|
| `main` (tuned cuts) | 37 | 46.9% | 3 | 0 |
| + A1 time limit, A3 heuristics | 43 | 12.8% | 3 | 0 |
| + A2 probing, coefficient tightening, nearly-integral node fix | 45 | 11.4% | 3 | 0 |
| + parallel search (run with 1 thread), pseudocost fix, dive backoff | 46 | 11.6% | 3 | 0 |
| + RINS time budget | 46 | 10.2% | 4 | 0 |
| + LP speed-ups (warm-start refactorization, propagation, optimality check) | 48 | 12.0% | 5 | 0 |
| + aggregated c-MIR with variable bounds, strong-branching cap | 47 | 9.3% | 6 | 0 |
| + simple generalized flow covers | 46 | 12.8% | 7 | 0 |
| + reduced-cost fixing at every node | 46 | 10.0% | 7 | 0 |
| + restart after the root (20% of integers fixed) | 46 | 10.0% | 7 | 0 |
| + plunge into the child with the smaller pseudocost estimate | 48 | 9.7% | 7 | 0 |
| + Feasibility Jump, dive step limit (same-day rerun of the previous build: 46, 12.3%, 7) | 51 | 10.2% | 7 | 0 |
| + cuts in the tree (pool + fresh c-MIR/covers), MIP start | 51 | see below | 7 | 0 |

**Done:**
- **A1:** runs end within the time limit (the node release time is reserved).
- **A3:** feasibility pump, four diving rules, RENS and RINS.
- **A2 lite:** probing on binaries, coefficient tightening.
- **Track B:** `cases/` (three MRPL families, the demo, 9/9 agree with HiGHS).
- **A6:** parallel tree search, `--threads N`.
- **A correctness fix:** a node whose LP point was integral only within the tolerance could be
  pruned (found by the MRPL crude case).
- **Feasibility Jump** (Luteberget and Sartor 2023, also in HiGHS 1.13+): LP-free weighted local
  search before the root cuts. New solutions on neos-3381206-awhea (474; optimum 453),
  graph20-20-1rand, neos-4738912-atrato, rococoC10-001000 and reblock115. Its point ignores the
  objective, so the root still runs the objective pump while it is the incumbent.
- **A dive fix it exposed:** a dive could repeat a step on a column already within the tolerance
  of its integral bound, without an LP iteration, until the time limit (neos-2657525-crna). The
  column is now snapped, and every step counts against the dive budget.
- **Cuts in the tree:**
  - Root cuts that end up non-binding go to a pool; nodes add the pool cuts their LP point
    violates.
  - Every 10th depth separates fresh c-MIR and cover cuts with the root bounds (so they are valid
    everywhere; Gomory cuts stay at the root, their tableau rows are only locally valid).
  - A node keeps its cuts only if they close 1% of its gap. They are the LP's last rows, so
    removing them restores the LP exactly.
  - Kept cuts may grow the LP by at most 25% of its rows.
  - The metric for cuts is the proven bound, not the incumbent: the screening's median primal
    gap moved by up to 5 points between two runs of the same binary.
  - Bound A/B, both builds at the same time, 60 s: the bound improved on 23 instances and got
    worse on 9. The mean distance of the bound from the optimum went from 22.4% to 22.0% (the
    median from 8.1% to 8.3%). Solved 7 and 7, same shifted geomean.
  - Largest gains: neos17 21.0% -> 8.5%, enlight_hard 16.2% -> 10.8%, binkar10_1 0.32% ->
    0.03%. Largest loss: neos-911970 5.5% -> 8.4% (the kept cuts slow its node LPs).
- **MIP start (`--mip-start`, `Params::mip_start`):** see `cases/README.md`, re-planning.
- **Conflict analysis** (dual proofs; Witzig, Berthold and Heinz 2017), `src/mip/conflicts.cpp`:
  - An infeasible node LP gives a proof from its Farkas ray, a node LP cut off by the incumbent one
    from its duals. Both use only the rows and their global bounds, so they hold at every node;
    sparse ones are kept (at most 1000) and propagated at every node.
  - Every Farkas ray is checked in the tests with `verify_infeasibility` on its node's bounds: a
    wrong mapping of the ray to the original rows only loses proofs (any multipliers give a valid
    one), so nothing else would catch it.
  - Bound A/B, 60 s: solved 7 -> 9 (neos17, nu25-pr12); bound better on 11, worse on 9; mean
    distance from the optimum 22.44% -> 22.29%, median 8.34% -> 7.65%. Largest loss enlight_hard
    16.2% -> 18.9%.

**Next, from the per-instance analysis** ([results/comparison.md](results/comparison.md)):
- neos-3381206-awhea now has a solution (Feasibility Jump) but its bound is stuck at 416 with no
  cuts found (optimum 453). HiGHS takes 8 s and SCIP 3 s.
- Path cuts on big-M networks: p200x1188c and mc11 keep large root gaps. Flow covers help
  sp150x300d but lower the root bound there.
- Symmetry handling (fhnw-binpack4-4, graph20-20-1rand): not before the deadline.

**Tried and dropped:**
- Single-row variable-upper-bound substitution in c-MIR (A4): no MIPLIB root bound moved. It
  pays off only with multi-row aggregation, which is now in.
- Sorted bound-flip groups in the ratio test: no measurable gain.
- Best-estimate node selection (every 5th pick): neutral (14 better, 16 worse).
- RENS retry with only at-bound fixings: neos-3381206-awhea's neighbourhood is still infeasible.
- Parameter tuning on the 62 screening instances: not done on purpose. It would overfit the
  set the results are reported on.
- Build flags: LTO +0.5%, `-march=native` +2.5%, profile-guided optimization about +4% on mas76.
  All within noise; the time is in sparse, memory-bound code, so the build stays simple.

Same-machine comparison with HiGHS, SCIP, CBC and GLPK:
[results/comparison.md](results/comparison.md).

## Rules for every milestone

- **Tests.** Each milestone adds tests to `tests/test_mip.cpp`, or to a new
  `tests/test_mip_*.cpp`:
  - a differential check against `tests/reference_milp.hpp` on the random families;
  - a debug-solution check: no reduction, cut or heuristic may cut off the known optimum;
  - a planted-bug check, mentioned in the commit message.
- **Build.** The debug, release and asan presets must be green with zero warnings before every
  push.
- **Screening.** A 60 s run of the 62 instances, the new build against the previous one on the
  same machine. A change is kept only if:
  - there are no wrong answers;
  - solved and feasible counts do not drop;
  - the shifted geomean does not get worse by more than 3% (the noise level).
- **One PR per milestone.**
- **Full runs.** A 10-minute run at A0, after A3 and at the end (A7). Results go into
  `docs/results/miplib.md`, with the machine named.
- **Same machine.** Numbers from the laptop and from the cloud are never compared with each
  other. Every comparison reruns both builds on the same machine.

## Track A: MILP

Ordered by expected gain per unit of work.

### A0. Baseline with the tuned cuts (benchmark only)
- 10-minute run of the current `main` (it has not been run with the tuned cuts).
- If possible, run HiGHS single-threaded on the same machine so we know the gap to close.
- **Decides:** the reference numbers for everything below.

### A1. Small fixes in the search (`src/mip/branch_and_bound.cpp`)
- **Time limit:**
  - check the clock inside strong branching and the cut loop;
  - give the node LP an iteration or time budget so a run stops within 1 s of the limit
    (it can overshoot by up to 8 s today).
- **Root restart:** if the root cuts and propagation fix more than 20% of the integers,
  presolve again and restart once.
- **Unit tests** for the time-limit contract.

### A2. MIP presolve (`src/presolve/`, MIP mode)
- **Probing on binaries:** fix each to 0 and 1, then propagate.
  - A side that is infeasible fixes the column.
  - Collect implications and implied bounds.
  - Detect equal and complementary binaries and merge them.
  - A work limit keeps the cost down.
- **Coefficient tightening** on ≤ rows with binaries (a_j := min(a_j, rhs − minact) style).
- **Clique detection:** set-packing rows plus probing implications give a clique table, shared
  with A4 and propagation.
- **Implied integers:** a continuous column that is forced to be integral is marked as such.
- **Dual reductions for integer columns** and dominated-column fixing.
- **Tests:**
  - postsolve checks with the verifier;
  - debug-solution checks;
  - planted reductions on random MILPs (as in `test_presolve.cpp`).
- **Expected to help:** p200x1188c, gen-ip*, mas*, neos5, set-partitioning-style instances.

### A3. Primal heuristics (`src/mip/heuristics.{hpp,cpp}`, a new file)
- **Diving:** fractional, coefficient (lock-based), guided (toward the incumbent) and
  pseudocost diving.
  - Diving runs from the root and periodically in the tree, using the node LP with a node
    budget.
- **Feasibility pump** at the root, with perturbation and restarts.
  - This targets the instances where we find no solution at all.
- **Sub-MIP heuristics:**
  - RENS at the root: fix integers that are already integral, round the domains of the rest;
  - RINS in the tree: fix where the incumbent and the LP solution agree.
  - Both are solved with our own `BranchAndBound` on the reduced model, with a node and time
    limit. This needs A1's limits to be exact.
- **Shift-and-propagate** as a cheap LP-free start heuristic.
- **Scheduling:** each heuristic tracks its success rate and the time it has spent; the budget
  is about 10% of the solve time.
- **10-minute benchmark checkpoint.**

### A4. More cuts, and cuts in the tree (`src/mip/cuts.cpp`)
- **Flow covers** on variable-upper-bound structures (x ≤ u·y). These are the p200x1188c,
  fixed-charge and network-design instances.
- **Clique cuts** from the A2 clique table, and **implied-bound cuts** from probing.
- **Zero-half cuts** (simple version) for pure binary instances.
- **Cut pool:**
  - cuts are stored with their age;
  - they are separated again at nodes up to depth 5 or so, and every k nodes deeper;
  - local cuts stay valid only in their subtree.
- **Tests:** the debug-solution check on every separator, and the planted-bug check.

### A5. Tree search (`src/mip/branch_and_bound.cpp`)
- **Conflict analysis:** when a node is infeasible after propagation, learn a clause and
  propagate it.
- **Node LP speed:**
  - reuse the factorization when bounds change only by a few columns;
  - skip strong branching deep in the tree when the pseudocosts are reliable.
- **Symmetry** is out of scope unless A0–A4 leave time. It needs orbit detection, which is a
  large piece of work.

### A6. Parallel tree search (`src/mip/parallel*`, std::thread only)
- **Structure:**
  - each worker owns a copy of the scaled LP and a Simplex;
  - the open-node queue, incumbent and pseudocosts are shared, behind a mutex, in batches;
  - the root runs once, and workers start from its basis.
- **Deterministic mode** (fixed synchronization points) for reproducible tests; an
  opportunistic mode for speed.
- **Option:** `--threads N`, default 1, so single-thread results stay comparable.
- **Target:** at least 2.5× speed-up at 8 threads on the solved instances.
  - This is measured on the laptop, since the cloud machine has only 4 cores.

### A7. Final numbers and documentation
- **10-minute runs**, single-threaded and 8 threads, next to HiGHS on the same machine.
- **Update:**
  - `docs/results/miplib.md`;
  - the status table in `docs/ARCHITECTURE.md`;
  - a short "how the MILP solver works" section for the judges.

## Track B: MRPL case studies (WP5)

This starts after A3, when the solver finds solutions reliably, or earlier in parallel if
there is time.
- **Model families** (defined in ARCHITECTURE.md WP5):
  - multi-period planning (LP);
  - crude scheduling with cargo windows (MILP);
  - blending (LP);
  - a utility-system unit commitment (MILP).
- **Sizes:** three per family. Each size is cross-checked against HiGHS.
- **Demo:** `cases/demo.sh` generates, solves, verifies and prints a readable plan in under
  2 minutes.

Barrier (WP3) and QP (WP4) are not in this plan. We decide after A4 whether there is time for
them. They are large, and PDLP (the teammate's WP1) already gives us a second LP method.

## Where benchmarks run: cloud or laptop

**Before every benchmark I'll ask you "local or cloud?"**, with these estimates:

| Run | Cloud (4 cores, 4 at a time) | Dell G15 (5 at a time) |
|---|---|---|
| 60 s screening, one build | ≈ 17 min | ≈ 13 min, and faster per core |
| 60 s screening, A/B builds | ≈ 35 min | ≈ 26 min |
| 10 min full run, one build | ≈ 2.6 h | ≈ 2.1 h |
| HiGHS baseline, 10 min | ≈ 2.6 h | ≈ 2.1 h |
| Parallel search (A6) scaling | not possible (4 cores) | **laptop only** |

**When the laptop is worth it:**
- **Always for full 10-minute runs and A6.** While the laptop benchmarks, I keep coding in the
  cloud. In the cloud, a benchmark and my builds share 4 cores, which makes the timings noisy.
- **The cloud is fine for 60 s screening** between small steps, so you don't have to be at the
  laptop.
- **Memory:** 16 GB on the laptop, and one run can hold up to about 2 GB of open nodes. Keep it
  to 5 runs at a time.

### Running a benchmark on the laptop

On Linux or WSL2, from a clone of the repository:

```sh
git fetch origin && git checkout <branch I name> && git pull
cmake --preset release && cmake --build --preset release
bench/fetch_instances.sh miplib-list bench/miplib_small.test     # first time only
bench/run_miplib.sh 5 600 bench/miplib_small.test build/release/apps/cli/samaya \
  bench/results/<run-name>
```

Send the results back on a branch, and I'll fetch and write up the report:

```sh
git checkout -b results/<run-name>
git add -f bench/results/<run-name>/*.csv && git commit -m "Benchmark results: <run-name>"
git push -u origin results/<run-name>
```

**Keep the laptop in its performance power mode, plugged in, and otherwise idle during the
run.**

For the HiGHS baseline, first run `pip install highspy`. Then add `--baseline highspy` to a
`bench/harness.py` call; I'll give you the exact command when we get there.
