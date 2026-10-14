# Simplex speed (FINALE_PLAN P4 / W4), from the large-LP run of 13 October

The dual simplex is the bottleneck of the headline run: after GPU PDLP reaches its 1e-4 point
in about 13 s, the crossover's dual simplex makes 86,912 iterations in about 59 minutes on the
467,376-row refinery LP, about **41 ms per iteration**, and does not finish within the hour
(`docs/results/large_lp.md`). On Netlib the same simplex is about 2x slower than HiGHS.

An independent code reading (Codex, astra, 13 October; no profiling was possible that night)
found **no quadratic step but no hyper-sparsity anywhere**: every iteration runs one BTRAN and
two FTRANs over the whole factor order, and several passes over all m basic rows. Its findings,
ranked by expected gain per day of work, are the plan below. The gains are its estimates, not
measurements; step 0 replaces them with measurements.

## 0. Measure first (0.5 day)

Counters and timers in the simplex, printed at log level 2 and returned in `SimplexStats`:
iterations per phase (dual, primal cleanup), time in FTRAN, BTRAN, CHUZR, the pivot row, the
ratio test, the primal and DSE updates, the Forrest-Tomlin update and rebuilds (with their
count and reason); densities of rho, alpha_col, tau and the flip right-hand side; factor and
update nonzeros; ratio-test candidate and group counts; rank repairs at a warm start.
Then a 10-minute run of the crossover on the large LP and the Netlib set. **Every later step
starts from that table.**

Build the release CLI, then run `bash bench/simplex_profile.sh MODEL --lp-method dual`
(add `--no-presolve` to measure the original matrix). `SAMAYA_BIN=/path/to/samaya` overrides
the default `build/release/apps/cli/samaya`. For crossover use the usual barrier/PDLP arguments
and `--crossover`. The driver extracts every `simplex profile begin` / `simplex profile end`
block, including any unscaled cleanup. Level 1 is unchanged; C++ callers can collect without
printing via `SimplexOptions::collect_stats`. Profiles are per `Simplex::solve`; the original
work counters retain their lifetime semantics. Warm-start iterations overlap the phase counts:
the current algorithm has no distinct repair iteration phase. Dual phase 1 includes the auxiliary
feasibility search. Density counts use exact nonzero comparisons, include failed pivot attempts,
and describe rho/alpha/tau after solving and the flip RHS before solving. Factor snapshots include
stored L/R eta entries and U's diagonal, both before and after each factorization attempt (including
rank-repair retries); a new factor's pre-snapshot is empty, while reuse retains the old updates.
Ratio counts cover dual ratio-test attempts.
The `rest` row includes setup, recomputation outside rebuilds, and measurement overhead; rebuild
timing includes factorization and primal/dual recomputation. Confirm the overhead gate on the idle
Linux Netlib run by comparing timing collection against the same binary without collection/logging,
and against the uninstrumented release baseline.
Local Release measurements exceeded the overhead gate when timing every iteration, so component
timers default to every 32nd iteration attempt (including retries). Counters and densities remain
exact, rebuild timers are unsampled, and the table labels the extrapolated component estimates.
Their common scale is attempts/samples, capped to the wall-time budget excluding rebuilds; raw
sampled seconds are also returned in `SimplexStats::sampled_seconds`. C++ callers can set
`SimplexOptions::timing_sample_interval = 1` for a full trace, or a larger stride to lower overhead.

## 1. Indexed vectors through the iteration (2–4 days; est. 1.3–2.5x here, 1.0–1.3x Netlib)

Today: CHUZR scans all m basics (`simplex.cpp:312`), the BTRAN right-hand side and the entering
column are zeroed over m (`:479`, `:80`), the primal update writes every basic value (`:534`),
the DSE pass runs over all m twice (`:542`), the FT update scans all m spike entries twice
(`basis_factor.cpp:371`). Change: dense-backed vectors with index lists of their nonzeros
(Hall and McKinnon 2005; HiGHS HEkkDualRHS), produced by the solves and used by the updates; an
incrementally maintained list of primal-infeasible rows for CHUZR. Risk: a row that becomes
infeasible and is missed would end the search early: keep the full recomputation before
accepting optimality, and the verifier.

## 2. Adaptive refactorization (1–2 days; est. 1.0–2.5x)

`refactor_interval = 100` (`simplex.hpp:58`) means about 869 rebuilds over that run; at 0.1, 1
or 3 s per rebuild of a 467k basis that is 1, 10 or 30 ms per iteration, so step 0 decides
whether this matters. Change: a work-based policy (rebuild when the accumulated update cost
exceeds the rebuild cost), keeping the fill safeguard and the numerical checks.

## 3. Hyper-sparse FTRAN/BTRAN and FT update (5–10 days; est. 2–6x, 5–15x with 1 and 2)

The solves skip zero components but still traverse the whole factor order; BTRAN visits every
stored L entry (`basis_factor.cpp:344`, `:363`). Change: symbolic reach (DFS over the factor's
column/row graph, Gilbert–Peierls; HiGHS HFactor::solveHyper) when the right-hand side is
sparse, with the dense path kept when the reach is large; the same for the FT row elimination
(`basis_factor.cpp:393`). The generated LP is a chain of 11,128 periods, which favours small
reaches, but a basis with long inventory chains can still reach far: step 0's densities decide.
Risk: permutation and eta-order mistakes; differential FTRAN/BTRAN/update tests against the
current code are required.

## 4. Dual Devex for warm starts (2–4 days; est. 1.1–1.5x throughput)

DSE weights are reset to 1 at a warm start (`simplex.cpp:825`), exact only for a slack basis;
every iteration also pays the extra FTRAN for tau. A dual Devex framework (HiGHS HEkkDual uses
it for a near-optimal non-slack start) removes that FTRAN; it may raise the iteration count,
so judge by verified solve time.

## 5. Ratio test without the quadratic group scan (1–3 days)

`simplex.cpp:358` rescans the remaining candidates for each Harris group: O(C^2) with C
candidates. Harmless with short lists, pathological with long ones; step 0's candidate counts
decide (HiGHS HEkkDualRow::chooseFinal reduces candidates first).

## 6. A better crossover basis (2–4 days; full push phases longer)

`crossover_basis` ranks variables by distance to bound against reduced cost and takes the best
m without checking independence (`crossover.cpp:68`); rank repairs use logicals (`simplex.cpp:94`),
which on equality rows start infeasible; only statuses survive, not PDLP's small residual. A
rank-aware (triangular) basis construction and a quality check first; Megiddo push phases
later. How many of the 86,912 iterations were repair is unknown until step 0.

## Order

Step 0 first, on the large LP and Netlib, alone on the machine. Then 1 and 2 (they are cheap and
their value shows in step 0's table), then 3; 4–6 as the table says. Every step: the Netlib and
random-family differential tests, the verifier unchanged, the four presets, a planted bug.
Gate (from POST_SUBMISSION_PLAN W4): Netlib geomean from 0.37 s to 0.25 s or better with nothing
lost, and the large-LP crossover finishing within the hour.
