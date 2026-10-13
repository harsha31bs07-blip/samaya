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
- **Next, without crossover:** PDLP on the GPU to its own tolerance 1e-8. On 8 October (before
  the tiled SpMV, 2.26x fewer seconds per iteration since) the KKT error at the last restart was
  4.7e-8 after 900 s. An answer at 1e-8 can pass the verifier (1e-6 relative) on its own, which
  would give a verified interior solution above 10⁶ columns. Queued after the overnight
  attribution runs, at 3,600 s; then the barrier (AMD ordering, b212d4c) on the same LP.

## Still to run

- The fixed Mittelmann subset (`bench/mittelmann_lp.test`): dual simplex, barrier, PDLP on the
  CPU and on the GPU, against HiGHS.
- This LP to convergence (1,800 s), then crossover (branch `claude/crossover`) and the verifier.
