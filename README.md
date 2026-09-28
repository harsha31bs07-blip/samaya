# samaya

[![Download samaya for Windows](https://img.shields.io/badge/Download-samaya%20for%20Windows%20(v0.1.0)-1F4E80?style=for-the-badge&logo=windows&logoColor=white)](https://github.com/harsha31bs07-blip/samaya/releases/latest/download/samaya-windows-0.1.0.zip)
[![Release notes](https://img.shields.io/badge/Release-v0.1.0%20notes-3F8A5A?style=for-the-badge)](https://github.com/harsha31bs07-blip/samaya/releases/latest)

> **Download: [samaya-windows-0.1.0.zip](https://github.com/harsha31bs07-blip/samaya/releases/latest/download/samaya-windows-0.1.0.zip)**
> (Windows 10/11, 64-bit). Extract it and double-click `setup.bat`. It installs **samaya Studio**,
> the desktop app, and the `samaya` command-line solver, with no build tools or admin rights.
> An NVIDIA GPU is optional. Install steps and requirements:
> [Releases](https://github.com/harsha31bs07-blip/samaya/releases/latest).

A sovereign LP / MILP / QP optimization solver core, built from mathematical foundations for
SIH problem statement 26119 (MRPL). No external solver library is used. See [PLAN.md](PLAN.md)
for the architecture, algorithms, benchmarks and timeline, and
[docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) for how the code fits together and the open work
packages (PDLP, GPU, barrier, QP, case studies). New to the project? Start with
[docs/ONBOARDING.md](docs/ONBOARDING.md): machine setup, build, and the Git/pull-request workflow.

## Status

| Component | State |
|---|---|
| Model, sparse matrix (CSC), model statistics | done |
| MPS / QPS reader (free + fixed with spaces in names, ranges, bound types, integer markers, QUADOBJ/QMATRIX) | done |
| C++ API, C API, `samaya` CLI (`--stats`, `--json`, `--solution`) | done |
| Scaling: geometric mean + equilibration, powers of two | done |
| Sparse LU: Markowitz + threshold pivoting, Forrest–Tomlin updates | done |
| Dual simplex: dual steepest edge, bound-flipping Harris ratio test, perturbation, phase 1 | done |
| Primal simplex (Devex) for cleanup and unboundedness | done |
| Warm start from a given basis; cleanup of unscaled infeasibilities | done |
| Independent verifier: optimality, Farkas certificates, unbounded rays | done |
| LP presolve + postsolve (primal and dual), verified on the original model | done |
| Hyper-sparse solves, barrier, PDLP (CPU + GPU) | phase 2–3 |
| MILP branch-and-bound: warm-started dual simplex, propagation, reliability branching, plunging, rounding heuristics | done |
| Root cutting planes: Gomory mixed-integer, c-MIR, lifted knapsack covers | done |
| Cuts in the tree, flow covers, parallel tree search, QP | phase 3–4 |

**Netlib: all 93 feasible instances solved, every objective matching HiGHS; 28 of the 29
infeasible instances proven infeasible with a verified Farkas certificate** (the remaining one,
`cplex2`, is infeasible by less than the tolerance and is reported as unproven). Per-instance
results and timings: [docs/results/netlib.md](docs/results/netlib.md). MILP:
[docs/results/miplib.md](docs/results/miplib.md) (MIPLIB 2017 subset, no wrong answers) and
[docs/results/generated-mip.md](docs/results/generated-mip.md) (refinery scheduling, knapsack,
facility location; all match HiGHS).

LP models are solved by the dual simplex. Every optimal solution, infeasibility certificate and
unbounded ray is checked by the independent verifier before it is reported; an outcome that does
not verify becomes `numerical_error`. MILP models are solved by branch-and-bound; the returned solution is checked against the original model (bounds, rows, integrality). QP models return `not_implemented` for now.

## Build

Requires CMake 3.22+, Ninja and a C++20 compiler (GCC 11+ or Clang 14+).

```sh
cmake --preset release
cmake --build --preset release
ctest --preset release
```

Presets: `debug` (warnings as errors), `release`, `asan` (AddressSanitizer + UBSan), `cuda`
(reserved for the Phase 3 GPU kernels).

**Download for Windows:** the ready-to-run package is on the
[Releases page](https://github.com/harsha31bs07-blip/samaya/releases). Download
`samaya-windows-<version>.zip`, extract it, and double-click `setup.bat`. It installs
samaya Studio and the `samaya` command-line solver, with no build tools and no administrator
rights needed. Requirements: 64-bit Windows 10 (1809 or later) or Windows 11. An NVIDIA GPU is
optional: with driver 580 or newer, Studio can run PDLP on it.

**Windows from source (native, no WSL):** double-click `setup.bat`. It checks the PC (Visual Studio 2022
Build Tools, Python, the WebView2 Runtime), installs what is missing after asking, builds with
MSVC, runs the tests and installs **samaya Studio**, the desktop app: drag and drop models, see
verified results, download them as Excel. Presets `windows-debug`, `windows-release`,
`windows-asan` and `windows-cuda` (with the GPU kernels). Details and the release package: [docs/WINDOWS.md](docs/WINDOWS.md).

## Usage

```sh
build/release/apps/cli/samaya --stats model.mps              # sizes and coefficient ranges
build/release/apps/cli/samaya model.mps                      # solve and verify
build/release/apps/cli/samaya --json model.mps               # JSON summary on the last line
build/release/apps/cli/samaya --solution sol.txt model.mps   # primal/dual values, certificates
build/release/apps/cli/samaya --help
```

C++:

```cpp
#include "samaya.hpp"

samaya::Model model = samaya::read_mps("model.mps");
samaya::Result result = samaya::Solver().solve(model);
```

C: see [`include/samaya_c.h`](include/samaya_c.h).

## Benchmarks

```sh
bench/fetch_instances.sh netlib netlib-infeas miplib   # downloads into bench/instances/
bench/generate_lps.py --scale 1                     # transportation, refinery planning, sparse
bench/harness.py bench/instances/netlib --baseline highspy --time-limit 300
```

Baseline solvers (the `highs` executable or the `highspy` Python module) are run only for
comparison and are never linked into samaya. The harness reports status/objective agreement and
the shifted geometric mean of solve times.

## Layout

```
include/        public API (samaya.hpp, samaya/*.hpp, samaya_c.h)
src/core/       model, solver dispatch, status, logging, C API
src/io/         file readers
src/linalg/     sparse matrices, scaling, basis LU with Forrest–Tomlin updates
src/lp/         dual and primal simplex, LP driver (scaling, unscaling)
src/presolve/   LP/MILP presolve and postsolve
src/mip/        branch-and-bound
src/verify/     independent solution and certificate checks
apps/cli/       samaya command-line tool
tests/          unit tests (self-contained framework), dense reference solvers, random LP
                generators, small instances
bench/          benchmark harness, instance generator and download script
```

Later phases add `src/qp` and `src/gpu`, as described in PLAN.md §3.
