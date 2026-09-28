# Getting started: taking a work package with Claude Code

A step-by-step guide for a teammate joining the solver work, using the Claude Code CLI (any
model configured behind it). It covers **WP1, PDLP on the CPU** (Step 5) and then **WP2, PDLP on
the GPU** (Step 5b): the GPU story of the problem statement. The same steps apply to any
package in [ARCHITECTURE.md](ARCHITECTURE.md) §6.

The idea: **you** own the package and the decisions; the agent writes code in small,
test-checked steps; **you** verify every step by running the commands yourself.

---

## Step 0: prerequisites (once)

New to Git, WSL or GitHub? Do [ONBOARDING.md](ONBOARDING.md) first (repository access, WSL,
tools, GitHub login, the branch and pull-request workflow), then continue here at Step 2.

Linux or WSL2 (Ubuntu 22.04 or newer) with:

```sh
sudo apt update
sudo apt install -y git cmake ninja-build g++ clang python3-pip curl
pip install highspy          # comparison baseline for benchmarks only
```

You need CMake 3.22+, GCC 11+ or Clang 14+. For WP2 (GPU) later you also need an NVIDIA GPU,
a recent driver and the CUDA 12 toolkit (`nvcc --version`).

## Step 1: clone, build, run the tests

```sh
git clone https://github.com/harsha31bs07-blip/samaya.git
cd samaya
cmake --preset release && cmake --build --preset release && ctest --preset release
```

Expected: `100% tests passed`. If this fails, stop and fix the environment first; the agent
cannot help if the baseline is broken.

Quick smoke test:

```sh
build/release/apps/cli/samaya tests/instances/tiny_lp.mps    # "Optimal objective 11 ... (verified)"
```

## Step 2: get the benchmark data

```sh
bash bench/fetch_instances.sh netlib            # ~95 LPs into bench/instances/netlib/
python3 bench/harness.py bench/instances/netlib --baseline highspy --time-limit 60
```

Expected at the end: `samaya solved 93/93` and `disagreements with baselines: 0`. This is the
bar your package must not break.

## Step 3: read, then branch

Read (30 minutes, worth it):
1. `CLAUDE.md`: the rules the agent will follow (it loads automatically).
2. `docs/ARCHITECTURE.md` §3–5 and your package in §6.
3. `src/lp/simplex.hpp` (the `LpProblem` struct) and `src/lp/lp_solver.cpp` (scaling and result
   mapping you will reuse).

Then:

```sh
git checkout -b wp1-pdlp
```

One branch per package. Never work on `main`.

## Step 4: start the agent and check it understood

From the repository root:

```sh
claude
```

First prompt, to confirm it loaded the rules:

> Summarize the hard rules in CLAUDE.md and the WP1 section of docs/ARCHITECTURE.md. List the
> files WP1 may create and the files it may edit. Do not write any code yet.

If the summary misses the verifier rule, the "no external libraries" rule or the ownership
map, correct it before going further.

**Working pattern for every milestone:**
- One milestone per session. Run `/clear` between milestones so old context doesn't confuse it.
- Ask for a plan first (plan mode: press Shift+Tab until it shows plan mode, or say
  "plan only"). Read the plan. Only then let it implement.
- After it says "done", **run the checks yourself** (Step 6). Do not trust a summary.
- Commit after each milestone passes.

## Step 5: the milestone prompts (WP1)

Copy these one at a time. Each is small on purpose.

### M1: plain PDHG on tiny LPs

> Implement WP1 milestone M1 from docs/ARCHITECTURE.md.
>
> Create `src/lp/pdlp.hpp` and `src/lp/pdlp.cpp` with the `PdlpOptions`, `PdlpResult` and
> `solve_pdlp` interface given there. Solve the `LpProblem` form directly:
> `min costᵀv  s.t.  K v = 0,  lower ≤ v ≤ upper` with `K = [A  -I]` (`LpProblem::A` and its
> row-wise copy `At`). The primal step projects onto the box `[lower, upper]` (infinite bounds
> allowed); the dual variable is free because all constraints are equalities. Use plain PDHG
> with fixed steps `τ = σ = 0.9 / ‖K‖₂`, estimating `‖K‖₂` by 30 steps of power iteration. No
> restarts yet.
>
> Termination: relative primal residual `‖Kv‖ / (1 + ‖v‖)`, dual residual, and relative gap all
> ≤ `tol`; iteration limit otherwise.
>
> Add `tests/test_pdlp.cpp`: 30 random feasible LPs from `tests/lp_generators.hpp`
> (`LpFamily::kFeasible`, at most 8 rows and 8 columns, only models where
> `ReferenceLp` says optimal and every column has a finite bound). Build the `LpProblem` the
> way `tests/test_lp_random.cpp` does in `lp_dual_steepest_edge_weights_stay_exact`. Compare
> the objective with the reference to 1e-4 relative, with `tol = 1e-6` and up to 200,000
> iterations.
>
> Register the files in `CMakeLists.txt` and `tests/CMakeLists.txt`. Do not edit any other
> file. Run the debug and release presets and show me the test output.

### M2: restarts, adaptive steps, primal weight, preconditioning

> Implement WP1 milestone M2 in `src/lp/pdlp.cpp`, following Applegate et al. (NeurIPS 2021):
>
> - Ruiz equilibration (10 passes) plus Pock–Chambolle (α = 1) diagonal preconditioning of K,
>   undone on the result.
> - Adaptive step size.
> - Primal weight updates at restarts.
> - Adaptive restarts on the normalized duality gap (restart to the average or the current
>   iterate, whichever has the smaller gap).
>
> Keep the M1 test and add a test over `kFeasible`, `kDegenerate` and `kBoxed` families
> (up to 20×20, 200 models) comparing with `ReferenceLp` to 1e-6 relative at `tol = 1e-8`.
> Also report the total iteration count before and after this change in the test output: the
> restarts must cut iterations substantially. Show me the test output.

### M3: infeasibility and unboundedness detection

> Implement WP1 milestone M3: detect primal and dual infeasibility from the difference of
> consecutive restart iterates, as in the PDLP paper. Return `SimplexStatus::kInfeasible` /
> `kUnbounded` with a candidate Farkas ray (dual) or primal ray in `PdlpResult`. Add tests on
> `LpFamily::kRandom` (often infeasible) and `kLoose` (often unbounded): every claim must pass
> `verify_infeasibility` / `verify_unbounded_ray` from `samaya/verify.hpp` (build a `Model` for
> the check or map the ray back to the model the test generated), and the status must agree
> with `ReferenceLp`. Show me the test output with the counts per status.

### M4: wire it into the solver

> Implement WP1 milestone M4.
>
> In `src/core/solver.cpp`, `solve_lp_model`: when `params.lp_method == LpMethod::kPdlp`, build
> the scaled `LpProblem` the same way `src/lp/lp_solver.cpp` does (reuse or factor out its
> helper; do not duplicate the unscaling code), run `solve_pdlp`, map `x`, `y`, rays back with
> the existing unscaling, and pass the result through the existing verifier.
>
> PDLP is less accurate than the simplex. Do not loosen the default verifier. Instead: if the
> PDLP point fails the strict verifier, polish it by calling `Simplex::solve(statuses)` from a
> basis guessed from the PDLP point (columns strictly inside their bounds basic, up to m of
> them, the rest at their nearest bound), and verify the polished result. Record
> `result.barrier_iterations` as the PDLP iteration count.
>
> Keep the change in `solver.cpp` small. Add a test that solves the random families through
> `Solver` with `lp_method = kPdlp` and checks status and objective against `ReferenceLp`, with
> `verified == true`. Run all three presets and show me the output. Then run
> `build/release/apps/cli/samaya --lp-method pdlp --json` on five Netlib instances.

### M5: Netlib

> Run `build/release/apps/cli/samaya --lp-method pdlp --json --log-level 0` on every file in
> `bench/instances/netlib/` with a 300 s time limit (a shell loop is fine; put it in
> `bench/pdlp_netlib.sh`). Make a table of instance, status, objective, iterations, seconds and
> whether the simplex polish was needed. Compare objectives with `bench/results/netlib.csv`.
> Write the table to `docs/results/pdlp_netlib.md` with a short summary at the top.

### After M5

Open a pull request (Step 7), get it merged, then continue with WP2 below.

## Step 5b: WP2, PDLP on the GPU (RTX 3060 laptop)

### CUDA setup (once)

**On Windows with WSL2** (recommended):
1. Install the latest NVIDIA Game Ready or Studio driver **on Windows**.
2. Inside WSL, do **not** install a Linux NVIDIA driver; the Windows driver is shared.
3. Install only the toolkit, from NVIDIA's "WSL-Ubuntu" CUDA 12 repository (developer.nvidia.com
   → CUDA Toolkit → Linux → x86_64 → WSL-Ubuntu), e.g. `sudo apt install cuda-toolkit-12-6`.

**On native Linux:** the NVIDIA driver plus `cuda-toolkit-12-x` from NVIDIA's repository.

Check:

```sh
nvidia-smi                         # shows the RTX 3060 and the driver version
nvcc --version                     # CUDA 12.x (add /usr/local/cuda/bin to PATH if missing)
nvidia-smi --query-gpu=name,compute_cap,memory.total --format=csv   # compute_cap 8.6
```

The `cuda` preset builds for compute capability 8.6 (`CMAKE_CUDA_ARCHITECTURES=86`); change it
in `CMakePresets.json` if you use a different GPU. The GPU build and tests:

```sh
cmake --preset cuda && cmake --build --preset cuda && ctest --preset cuda
```

The CPU presets must keep working without CUDA installed; CI has no GPU.

### WP2 milestone prompts

#### G0: build wiring

> Wire CUDA into the build for WP2 (docs/ARCHITECTURE.md). In `CMakeLists.txt`, inside the
> existing `if(SAMAYA_CUDA)` block, add the sources in `src/gpu/` to the `samaya` library, link
> `CUDA::cudart`, and define `SAMAYA_HAVE_CUDA`. Create `src/gpu/device.hpp` / `device.cu`
> with a function `bool gpu_available()` (uses `cudaGetDeviceCount`) and a `SAMAYA_CUDA_CHECK`
> macro that turns any CUDA error into an exception with the file, line and error string.
> Without `SAMAYA_HAVE_CUDA`, `gpu_available()` must still exist and return false (a `.cpp`
> stub). Add `tests/test_gpu.cpp` with one test that adds two vectors on the device and
> compares with the CPU; when no GPU is present it prints `SKIPPED (no GPU)` and passes. Show
> me: the `cuda` preset build and test output, and that the `debug` preset still builds and
> passes without CUDA.

#### G1: kernels

> Implement the WP2 kernels in `src/gpu/`, double precision throughout: (1) CSR SpMV
> `y = A x` with one warp per row, used for both products: the CSC arrays of `LpProblem::At`
> are `A` in CSR form, and the CSC arrays of `LpProblem::A` are `Aᵀ` in CSR form, so upload both
> as-is (the `-I` block of `K = [A -I]` is handled inside the vector kernels); (2) fused vector
> kernels: `axpby`, projection onto `[lower, upper]` with infinite bounds, and the PDLP primal
> and dual update steps; (3) reductions for dot products and 2-norms (two-stage block
> reduction, no atomics on doubles). Do **not** use cuSPARSE or cuBLAS. Tests in
> `tests/test_gpu.cpp`: random sparse matrices (including empty rows and one very long row) and
> vectors, every kernel compared with a CPU implementation to 1e-12 relative. Show me the test
> output.

#### G2: the PDLP loop on the device

> Implement WP2 milestone M2: a GPU path in `solve_pdlp` used when
> `PdlpOptions::use_gpu && gpu_available()`, otherwise the CPU path. Upload the problem and the
> iterates once; keep every vector on the device during iterations; copy only the scalars
> needed for restart and termination checks back to the host (at most once per check, not per
> iteration). Reuse the CPU implementation's restart and step-size logic on the host. Test:
> on 50 random LPs, run 100 iterations on the CPU and on the GPU with identical settings and
> check the iterates agree to 1e-10 relative; then check full solves agree in status and
> objective with the CPU path. Wire `--gpu` (already parsed into `Params::use_gpu`) to
> `PdlpOptions::use_gpu`. Show me the test output and `nvidia-smi` during a long solve.

#### G3: benchmark

> Generate large LPs with `python3 bench/generate_lps.py --set lp --scale 10` and
> `--scale 20` (write to `bench/instances/generated-large/`). Run each with
> `--lp-method pdlp` on the CPU and with `--lp-method pdlp --gpu`, same tolerance (1e-4 and
> 1e-6), and record: status, objective, iterations, total seconds, seconds excluding the
> host-to-device upload, and whether the result verified. Also run the dual simplex on the same
> models for reference. Write the table to `docs/results/gpu_pdlp.md` with the GPU name, driver
> and CUDA version, and a one-paragraph honest summary (speedup, where the GPU does not help).

### GPU red flags (reject these)

- cuSPARSE/cuBLAS/cuSOLVER in the solver path (allowed only in a separate benchmark tool).
- `float` anywhere in the solver path, or "fast math" flags.
- Kernel launches without error checks, or `cudaDeviceSynchronize` after every kernel.
- A host-device copy of a full vector inside the iteration loop.
- GPU tests that silently pass without a GPU instead of printing `SKIPPED`.
- Speedups quoted without the same tolerance, or without verified results.

**What to expect on a laptop RTX 3060:** PDLP is limited by memory bandwidth (about 336 GB/s
on the GPU against roughly 60–80 GB/s from DDR5 for the CPU), and consumer GPUs run double
precision at a small fraction of their single-precision speed. A verified **2–4×** on large
LPs is a good, honest result; small LPs will be faster on the CPU (launch and transfer
overhead). Say so in the report.

### Who owns what

The MILP work (`src/mip/`, MILP parts of `bench/`) is being done in parallel by the lead. WP1/WP2
own `src/lp/pdlp.*`, `src/gpu/`, the CUDA block of `CMakeLists.txt`, `tests/test_pdlp.cpp` and
`tests/test_gpu.cpp`. Both sides touch `src/core/solver.cpp` only in small dispatch hunks:
pull `main` before starting each milestone and rebase before opening a PR.

## Step 6: verify every step yourself

After each milestone, run these yourself. This is the most important step.

```sh
git status && git diff --stat                     # only the files you expected?
cmake --preset debug   && cmake --build --preset debug   && ctest --preset debug
cmake --preset release && cmake --build --preset release && ctest --preset release
cmake --preset asan    && cmake --build --preset asan    && ctest --preset asan
python3 bench/harness.py bench/instances/netlib --baseline highspy --time-limit 60   # still 93/93
```

Red flags to reject:
- Tests deleted, skipped, commented out, or their tolerances loosened.
- `#pragma` or compiler flags that silence warnings.
- Edits outside the package's files (see the ownership map), especially `src/verify/`.
- An external library added to `CMakeLists.txt`.
- "Done" without the test output shown.

**Planted-bug check** (do it once per milestone, it takes two minutes): flip a sign in the new
code (for example the dual step), rebuild, and confirm your new tests **fail**. Revert. A test
that can't fail proves nothing.

## Step 7: commit, push, pull request

```sh
git add -A
git commit -m "WP1 M1: plain PDHG on the LpProblem form with reference tests"
git push -u origin wp1-pdlp
```

Open a PR into `main` on GitHub. The description should state what works, the test results
(paste the counts), the planted-bug check, and the known gaps. CI runs GCC and Clang in debug,
release and ASan; all six jobs must be green. Tag the lead for review before merging.

## Step 8: tips for working with the agent

- **Small steps win.** If a prompt produces a large, tangled change, `git checkout .` and ask
  for half of it.
- **Name files and functions** in your prompts; don't make it search.
- **When it loops** (same error three times), `/clear`, paste the exact error and the relevant
  20 lines, and ask for a diagnosis before a fix.
- **Numerical bugs**: ask it to add a tiny hand-checkable LP (2×2) as a test and print the
  iterates; compare with a hand calculation.
- **Never run benchmarks while rebuilding** the binary they use; results become meaningless.
- **Ask for evidence, not claims**: "show me the ctest output", "show me the diff of
  solver.cpp".

## Troubleshooting

| Symptom | Fix |
|---|---|
| `cmake --preset` unknown | CMake older than 3.22: install a newer one (`pip install cmake`) |
| `ninja: command not found` | `sudo apt install ninja-build` |
| `bench/fetch_instances.sh` fails on netlib.org | network or proxy; try again, or copy `bench/instances/` from a teammate |
| `highspy` import error | `pip install --upgrade highspy` in the same Python the harness uses |
| Debug build fails on warnings | that is intended (`-Werror`); fix the warning, never the flag |
