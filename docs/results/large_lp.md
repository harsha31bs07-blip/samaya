# Large LPs: PDLP on the GPU

The generated refinery LP `bench/instances/large/refinery_20c_10p_11128t.mps` (in
`bench/gpu_lp.test`): 467,376 rows, 1,001,520 columns, 4,117,330 nonzeros. Laptop: Dell G15,
Intel Core i5-13450HX, NVIDIA GeForce RTX 3050 6 GB Laptop GPU (driver 617.42), CUDA 12.0, WSL2,
`cmake --preset cuda`. Each run alone on the machine.

## Device-side evaluation (8 October 2026)

`samaya --lp-method pdlp --gpu --time-limit 600 --log-level 1`, before and after the commit
that computes PDLP's residuals, restarts and restart distances on the device (branch
`claude/w2-gpu-eval`, e1831a5) instead of copying the iterates back every 64 iterations.

| Build | Status at 600 s | Iterations | Iterations per second |
|---|---|---|---|
| fa76f94 (host evaluation) | time limit | 76,160 | 127 |
| e1831a5 (device evaluation) | time limit | 105,280 | 175 |

- **38% more iterations per second.** Neither run converges in 600 s: PDLP's relative
  tolerance here is 1e-8 (`kPdlpTolerance`), and in a 900 s run at log level 2 the KKT error at
  the last restart had fallen to 4.7e-8 (147,136 iterations), against 1.5e-7 (80,832) before.
- **Log level 2 is not usable on this model:** it prints a line per free column at every
  infeasibility check, 16 GB in 900 s before the change and 3.8 GB after. Those lines belong at
  level 3 (to do).
- **What limits the GPU now:** about 5.7 ms per iteration, far from what the memory bandwidth
  allows for two products with 4.1 million nonzeros. The SpMV kernel (`src/gpu/spmv.cu`) uses
  one warp per row, while rows of A average about 9 nonzeros and rows of A' about 4, so most
  lanes idle. A row-per-thread or merge-based kernel for short rows is the next step, measured
  the same way.

## GPU PDLP + crossover (13 October 2026)

`samaya --lp-method pdlp --gpu --crossover --time-limit 3600`, build f377d79 (tiled SpMV
3f63fb2, crossover 4ac6de6, time limits 22b767b), `cmake --preset cuda`, alone on the laptop,
clock factor 1.000 (`bench/clock_factor.py`).

| Phase | Where | Outcome |
|---|---|---|
| PDLP to its crossover tolerance 1e-4 | GPU | 5,312 iterations, within the first minute of the run (about 13 s at the ~400 iterations/s measured on 8 October) |
| Crossover: dual simplex from the crossover basis | CPU, one thread | 86,912 iterations in about 59 min, then the time limit |
| Total | | time limit at 3,602.6 s, no verified answer; peak memory 1.25 GB |

- **The GPU part is fast and the CPU simplex is the bottleneck.** At about 25 iterations per
  second on 467,376 rows the dual simplex cannot finish a basis of this size within an hour;
  this is what WP6/W4 (hyper-sparse FTRAN/BTRAN, profile first) targets, and a better crossover
  basis would cut the iterations.
- **No verified solve above 10⁶ columns yet.** The slides must not claim one.

## GPU PDLP to 1e-8, without crossover (13 October 2026)

`samaya --lp-method pdlp --gpu --time-limit 3600 --log-level 1`, build f377d79, alone on the
laptop, clock factor 1.000.

| Phase | Outcome |
|---|---|
| PDLP on the GPU, tolerance 1e-8 | stopped by its **automatic iteration cap of 1,000,000**, not by the tolerance or the time limit |
| Fallback: dual simplex from scratch (presolve removed nothing) | 86,400 iterations, then the time limit |
| Total | time limit at 3,601.0 s, no verified answer; peak memory 1.04 GB |

- **The run measured the cap, not PDLP.** `pdlp.cpp` capped automatic iterations at
  min(1000·(m + n), 10⁶) even when a time limit was given; at about 400 iterations per second
  that cap ends PDLP well inside the hour. The level-1 log did not say when, so the split
  between PDLP and the simplex is unknown. Fixed afterwards: with a time limit the 10⁶ ceiling
  no longer applies (the size-proportional cap, 1.5·10⁹ here, still does), and the level-1
  lines for the iteration and time limits give the elapsed time and the KKT error at the last
  restart.
- **Next:** the same run again with the fix, to learn whether 1e-8 is reachable within the hour.
  Then the barrier (AMD ordering, b212d4c) on the same LP.

## GPU PDLP to 1e-8, rerun after the cap fix (14 October 2026)

`samaya --lp-method pdlp --gpu --time-limit 3600 --log-level 1`, build d86b0fa, alone on the
laptop, clock factor 1.016 (WSL clock 1.6% fast). Raw data:
`bench/results/large_lp/20261014-0434-d86b0fa-pdlp1e8/`.

| Phase | Outcome |
|---|---|
| PDLP on the GPU, tolerance 1e-8 | 1,413,632 iterations (about 390 per second), then the **time limit**; not converged |
| Total | time limit at 3,600.9 s, no verified answer; peak memory 0.70 GB |

- **1e-8 is not reachable within the hour on this LP.** The weighted KKT error at the last
  restart was 5.5e-9. That number weights the primal residual by the primal weight ω and the
  dual one by 1/ω. The stopping test needs the relative primal residual, dual residual and gap
  each below 1e-8 without weights, and it was not met. This log did not show the unweighted
  values; the limit messages now print them from the last evaluation.
- Peak memory dropped from 1.04 GB to 0.70 GB because no simplex fallback ran this time.
- **Consequence:** PDLP alone does not give a verified solve above 10⁶ columns here. The path
  goes through crossover, and so through the simplex speed work in `docs/SIMPLEX_SPEED_PLAN.md`
  (its step 0 measured 35 ms per crossover iteration and 110,681 rank repairs of the crossover
  basis).

## GPU PDLP + crossover after simplex steps 1b and 2 (14 October 2026)

`samaya --lp-method pdlp --gpu --crossover --time-limit 3600 --log-level 2`, build f837aba,
alone on the laptop. Raw data: `bench/results/hour-f837aba/`.

| Build | Crossover iterations in the hour | ms per iteration | Outcome |
|---|---:|---:|---|
| f377d79 (13 October) | 86,912 | about 41 | time limit |
| f837aba | 128,576 | 27.8 | time limit |
| f837aba + structural repair (not committed) | 107,136 | 33.4 | time limit |

- The simplex is 1.48x faster per iteration than on 13 October, but crossover still does not
  finish.
- The structural repair takes the 110,681 dependent positions with structural columns instead
  of equality logicals. It made the basis denser and the iterations slower, and the start no
  better, so it stays out (the patch is kept outside the repository).

**Why the start is bad (59d4874, `bench/results/snap-diag/`).** PDLP's point is in the interior
of a large optimal face:
- 357,920 nonbasic variables are off their bounds, 323,335 of them by more than 1e-3 relative.
- `crossover_basis` returns statuses only, so all of them are snapped to a bound. Before the
  snap, 7 basic values violate their bounds; after it, 153,844 do. The basic values move by up
  to 6.2e9.
- The dual simplex then works with about 212k infeasible rows throughout. That count is the true
  number, not the list length.
- More simplex speed alone cannot fix this. A push phase keeps PDLP's values and moves each
  superbasic to a bound by ratio tests, so the basics stay feasible (Megiddo; the IPX
  crossover). It is in progress. At about 5 ms per FTRAN, 358k pushes would take about
  30 minutes.

## Still to run

- The fixed Mittelmann subset (`bench/mittelmann_lp.test`): dual simplex, barrier, PDLP on the
  CPU and on the GPU, against HiGHS.
- This LP to convergence (1,800 s), then crossover (branch `claude/crossover`) and the verifier.
