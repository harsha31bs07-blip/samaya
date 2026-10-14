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

### Step 0 results (14 October 2026, WSL on the laptop, alone on the machine)

Build d86b0fa; clock factor 1.016 (the WSL clock ran 1.6% fast against the Windows clock, so
seconds below are 1.6% high; ratios are unaffected). Raw data:
`bench/results/step0/20261014-0415-d86b0fa/`.

**Overhead gate: passed.** All 93 Netlib models, dual simplex, the instrumented release binary
against the uninstrumented 905c1a5 build, alternating, best of 2 rounds: geomean new/base
**1.0036** over the 27 models with base time ≥ 0.05 s (gate ≤ 1.03); the worst ratio is 1.052
(fit1p, 0.06 s). Iteration counts and objectives are identical on all 93; all optimal.

**Large LP, crossover** (`--lp-method pdlp --gpu --crossover --time-limit 600`, cuda build):
16,384 dual phase-2 iterations in 581 s, **35.4 ms per iteration**, no answer within 600 s.

| Component | Share | ms/iter |
|---|---:|---:|
| rebuild + recompute | **37.2%** | 13.18 |
| CHUZR | **14.9%** | 5.27 |
| BTRAN | 11.1% | 3.93 |
| entering-column FTRAN | 10.1% | 3.58 |
| DSE FTRAN (tau) | 9.7% | 3.44 |
| Forrest-Tomlin update | 6.0% | 2.11 |
| primal update | 4.3% | 1.54 |
| DSE weight update | 3.2% | 1.12 |
| pivot row | 1.8% | 0.65 |
| ratio test, bound flips, dual update, rest | 2.1% | 0.63 |

- 164 rebuilds, 163 of them for the 100-iteration interval, each about 1.26 s of factorization
  (L ≈ 0.33M, U ≈ 1.24M stored entries).
- Densities are low: rho 1.7%, alpha_col 1.5%, tau 2.5% of m on average.
- **110,681 rank repairs at the warm start**: almost a quarter of the 467,376 basic positions
  of the crossover basis were replaced by logicals.
- Ratio test: 4,867 candidates per iteration on average, 1.55 Harris groups.

**Netlib** (93 models, dual simplex, 41.3 profiled seconds; dfl001 alone is 20.7 s and pilot87
7.9 s):

| Component | Share |
|---|---:|
| pivot row | **19.5%** |
| ratio test | **14.5%** |
| rebuild + recompute | 14.3% |
| DSE FTRAN (tau) | 13.4% |
| entering-column FTRAN | 10.9% |
| BTRAN | 6.5% |
| dual/reduced-cost update | 6.1% |
| Forrest-Tomlin update | 4.7% |
| CHUZR | 3.3% |
| DSE weight update | 3.2% |
| bound-flip FTRAN, primal update, rest | 3.7% |

- Median densities: rho 1.2%, alpha_col 1.3%, tau 2.8% of m.
- Rebuild reasons: 1,092 interval, 118 fill limit, 152 phase, 9 confirmation; no numerical or
  mismatch rebuilds.
- Harris groups average about 1 (dfl001: 22,703 groups in 22,696 iterations): the quadratic
  rescan of step 5 does not occur. The ratio test costs time because the candidate lists are
  long (dfl001 2,285, pilot87 1,279 per iteration), not because of the group scan.

**What the table changes.**
- Large LP: step 2 comes first. The interval rebuilds are 37% of the time, and a work-based
  policy is the cheapest change in the plan. The factorization itself (1.26 s for this basis)
  is a separate question for step 3.
- Step 1 matters for the large LP (CHUZR 15%). On Netlib its share is small. There the pivot
  row (19.5%) leads, but it is already computed row by row over rho's nonzeros
  (`compute_pivot_row`, using the row-wise copy `At`). What remains there is the O(m) scan for
  those nonzeros, which needs rho's index list from a hyper-sparse BTRAN (step 3). Step 1
  therefore starts with CHUZR over a maintained list of infeasible rows (1b), with ties broken
  as the scan does, so iteration paths stay identical.
- Step 3: FTRAN and BTRAN together take 31% of the large LP and 31% of Netlib, at densities
  around 2%, so hyper-sparse solves are worth building.
- Step 6 is likely the largest lever for the hour: the crossover basis needs 110,681 repairs. At
  35 ms per iteration an hour is about 100,000 iterations. Even a 3x faster iteration does not
  finish if the repaired basis needs millions.
- Step 5 as written is dropped. A cheaper ratio test means shorter candidate lists (HiGHS's
  BFRT works on a pre-filtered list), which can be done later.

## 1. Indexed vectors through the iteration (2–4 days; est. 1.3–2.5x here, 1.0–1.3x Netlib)

Today: CHUZR scans all m basics (`simplex.cpp:312`), the BTRAN right-hand side and the entering
column are zeroed over m (`:479`, `:80`), the primal update writes every basic value (`:534`),
the DSE pass runs over all m twice (`:542`), the FT update scans all m spike entries twice
(`basis_factor.cpp:371`). Change: dense-backed vectors with index lists of their nonzeros
(Hall and McKinnon 2005; HiGHS HEkkDualRHS), produced by the solves and used by the updates; an
incrementally maintained list of primal-infeasible rows for CHUZR. Risk: a row that becomes
infeasible and is missed would end the search early: keep the full recomputation before
accepting optimality, and the verifier.

**1b done (2230a55): CHUZR over a list of primal-infeasible rows.** The large-LP crossover at
2230a55, still with the 100-update interval, 600 s
(`bench/results/step0/2230a55-crossover/`):
- CHUZR fell from 5.27 to 3.26 ms per iteration, and the iteration from 35.4 to 33.2 ms.
- The primal update rose from 1.54 to 1.85 ms, which is the cost of keeping the list.
- The gain is small because the list averages 215,885 entries, 46% of the 467,376 rows: about
  half the rows of this warm start are primal infeasible.
- On the MIPLIB A/B against 0a149c8 (60 s, four processes at once), simplex iterations per second
  had a geomean of 0.957. pk1 took the identical path 7% slower. The A/B also contains the
  presolve commits and its noise is large, so 1b's own overhead on small LPs is measured
  separately (`bench/results/step1b-overhead/`, sequential, 999caca against 2230a55).

## 2. Adaptive refactorization (done: 45487d0)

**Experiment (14 October 2026).** The large-LP crossover, 600 s, with the interval forced
through a throwaway environment override (raw data: `bench/results/step0/refactor-interval/`):

| Interval | Iterations | ms/iteration | Rebuild share | DSE FTRAN ms | R entries before a rebuild |
|---:|---:|---:|---:|---:|---:|
| 100 | 16,384 | 35.4 | 37% | 3.44 | 8k |
| 300 | 20,192 | 28.7 | 17% | 4.08 | 26k |
| 1000 | 21,856 | 26.4 | 5% | 4.43 | 122k |

No numerical, mismatch or unsafe-update rebuild at any interval. The solves get slowly dearer
as the Forrest-Tomlin etas grow, so the gain flattens. 1000 gives 1.34x iterations per hour.

**Policy (45487d0).** An amortized rule, deterministic (entries, not time):
- After each update, add 3 solves × (factor entries now − entries after the last
  factorization).
- Refactorize when that sum reaches `rebuild_cost` × the entries after factorization.
- `rebuild_cost` = 200 is calibrated from the step-0 profiles. The ratio of factorization
  cost per entry to solve cost per entry is about 70 on most Netlib models and 150–340 on
  pilot, pilot87, fit2p and the large LP. The optimal interval goes with its square root.
- `refactor_interval` is now a cap of 2000.

**Measured (`bench/results/step2/45487d0/`, against 2230a55, alone on the machine).**
- Large-LP crossover, 600 s: 25.3 ms per iteration, down from 33.2. That is 22,976 iterations
  instead of 17,408, about the same as the best fixed interval (26.4 ms at 1000).
- 21 rebuilds: 16 for update work, 4 at the 2000 cap, and the initial one. The rebuild share
  fell from 40% to 5.4%.
- Netlib, dual simplex, best of 2 rounds: all 93 models optimal with both builds. The time
  ratio geomean is 0.995, so the change is neutral on Netlib. 33 models take a different path,
  because rebuilding at a different point changes the rounding.
- Largest Netlib changes: pilot 2.21 → 1.74 s (0.79), nesm 0.85, greenbea 0.86; dfl001 1.05
  and greenbeb 1.04 are the slowest of the models above 1 s.

### Original plan

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

**Part 1 (1757a8d): repair inside `crossover_basis`.** The 110,681 dependent positions are now
replaced in 2 factorizations and 0.78 s, by logicals of the unpivoted rows. Each leaving variable
goes to the bound nearest its PDLP value. The warm start then needs no rank repair.

What it did not change: the dual simplex still starts with about half the rows primal
infeasible, and runs 17,408 iterations in 600 s without finishing. The final basis holds
110,681 equality-row logicals, the same number as the repairs. Those logicals are fixed, so a
basic one is infeasible whenever its row's activity is off the right-hand side. A basis that is independent from
the start (the triangular crash of part 2) is the next lever.

## Order

Step 0 is done (results above). Then 2 and 1b (the CHUZR list), then 3, which also gives the
index lists for the rest of step 1. Step 6 can run alongside, since it touches
`crossover.cpp`, not the iteration; its first part (repair at PDLP-nearest bounds) is in
1757a8d. Then 4; step 5 is dropped. Every step: the Netlib and
random-family differential tests, the verifier unchanged, the four presets, a planted bug.
Gate (from POST_SUBMISSION_PLAN W4): Netlib geomean from 0.37 s to 0.25 s or better with nothing
lost, and the large-LP crossover finishing within the hour.
