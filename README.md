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
| Barrier (interior point) for LP and convex QP; PDLP on the CPU and on NVIDIA GPUs (own CUDA kernels) | done |
| Hyper-sparse solves | planned |
| MILP branch-and-bound: warm-started dual simplex, propagation, reliability branching, plunging, rounding heuristics | done |
| Root cutting planes: Gomory mixed-integer, c-MIR, lifted knapsack covers | done |
| Cuts in the tree, flow covers, conflict analysis, Feasibility Jump, parallel tree search | done |
| Convex QP (interior point) and MIQP (branch and bound over QP relaxations) | done |
| Windows: native MSVC build, samaya Studio desktop app, `setup.bat` installer | done |

**Netlib: all 93 feasible instances solved, every objective matching HiGHS; 28 of the 29
infeasible instances proven infeasible with a verified Farkas certificate** (the remaining one,
`cplex2`, is infeasible by less than the tolerance and is reported as unproven). Per-instance
results and timings: [docs/results/netlib.md](docs/results/netlib.md). MILP:
[docs/results/miplib.md](docs/results/miplib.md) (MIPLIB 2017 subset, no wrong answers) and
[docs/results/generated-mip.md](docs/results/generated-mip.md) (refinery scheduling, knapsack,
facility location; all match HiGHS).

LP models are solved by the dual simplex. Every optimal solution, infeasibility certificate and
unbounded ray is checked by the independent verifier before it is reported; an outcome that does
not verify becomes `numerical_error`. MILP models are solved by branch-and-bound; the returned solution is checked against the original model (bounds, rows, integrality). Convex QP models are solved by the interior-point method (106 of the 138 Maros–Mészáros models solved and verified: [docs/results/qp.md](docs/results/qp.md)), and MIQP models by branch and bound over QP relaxations.

## Build and install

### Windows: one step, `setup.bat`

**Ready-to-run (no build):** download `samaya-windows-<version>.zip` from the
[Releases page](https://github.com/harsha31bs07-blip/samaya/releases/latest), extract it, and
double-click **`setup.bat`**.
- It checks the PC and installs the Microsoft Edge WebView2 Runtime if it is missing (asking
  first).
- It installs **samaya Studio** and the `samaya` command-line solver for the current user.
- It needs no build tools and no administrator rights.

Requirements:
- 64-bit Windows 10 (1809 or later) or Windows 11.
- An NVIDIA GPU is optional (driver 580 or newer); without one, everything runs on the CPU.

**From source:** clone the repository and double-click **`setup.bat`** in its folder. It
does the whole build:
1. **Checks** Windows, winget, the Visual Studio 2022 Build Tools (C++ tools, CMake, Ninja),
   Python 3 and the WebView2 Runtime.
2. **Installs what is missing** through winget, asking first (the Build Tools need one
   administrator prompt).
3. **Builds** samaya with MSVC (`windows-release`) and **runs the tests**.
4. **Generates** the sample refinery cases, **installs** Studio and the solver, and **opens**
   Studio.

| Command | What it does |
|---|---|
| `setup.bat` | check, install what is missing (asking first), build, install, open Studio |
| `setup.bat /check` | report what is present and what is missing; change nothing |
| `setup.bat /yes` | install missing components without asking |
| `setup.bat /nolaunch` | do not open Studio at the end |
| `setup.bat /uninstall` | remove samaya Studio for this user |

Studio and `samaya.exe` go to `%LOCALAPPDATA%\Programs\samaya`, with Start menu and desktop
shortcuts. The log is in `%LOCALAPPDATA%\samaya\setup.log`.

**Developers on Windows:** `tools\windows\build.ps1 -Preset <preset> [-Test]` builds one preset
from any shell. It finds Visual Studio and enters its developer shell.
- Presets: `windows-release`, `windows-debug`, `windows-asan`, and `windows-cuda`.
- `windows-cuda` has the GPU kernels and needs the NVIDIA CUDA Toolkit (12 or 13).
- `tools\windows\package.ps1` makes the release zip.
- Details: [docs/WINDOWS.md](docs/WINDOWS.md).

### Linux

Requires CMake 3.22+, Ninja and a C++20 compiler (GCC 11+ or Clang 14+).

```sh
cmake --preset release
cmake --build --preset release
ctest --preset release
```

Presets:
- `debug`: warnings as errors.
- `release`.
- `asan`: AddressSanitizer + UBSan.
- `cuda`: release with the GPU kernels for PDLP; needs the CUDA Toolkit.

## Usage

### samaya Studio (Windows)

Open **samaya Studio** from the Start menu or the desktop.
1. **Add a model:** drop `.mps` or `.qps` files on the window, use **Open models** (Ctrl+O), or
   press **Solve** next to a bundled MRPL-style sample.
2. **Choose the settings:** start date, time limit, threads, MIP gap, LP method, presolve, and
   **Compute**.
   - Compute is CPU or the NVIDIA GPU. The GPU runs PDLP and is meant for very large linear
     programs; small and medium models solve faster on the CPU.
   - Verification is always on.
3. **Read the results:**
   - the status and the verified checks (samaya's verifier and Studio's independent re-check);
   - the schedule or plan on calendar dates (a berth Gantt chart, tank stocks, weekly plan,
     hourly unit commitment);
   - every variable and constraint, and the solver log.
4. **Export** to Excel (real date cells), CSV, or samaya's solution file.
5. **Run again or compare:**
   - **Back** or **Run again** solves the same model with other settings.
   - **Compare cases** puts runs of the same model side by side.

### Command line

On Windows after `setup.bat`, the solver is
`%LOCALAPPDATA%\Programs\samaya\samaya.exe`. On Linux it is
`build/release/apps/cli/samaya`.

```sh
samaya --stats model.mps                      # sizes and coefficient ranges
samaya model.mps                              # solve and verify
samaya --json model.mps                       # JSON summary on the last line
samaya --solution sol.txt model.mps           # primal/dual values, certificates
samaya --time-limit 60 --threads 8 --mip-gap 1e-4 model.mps
samaya --lp-method pdlp --gpu model.mps       # PDLP on an NVIDIA GPU (CUDA builds)
samaya --gpu-info                             # does this build have CUDA, and which GPU
samaya --help
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
