# samaya on native Windows: build, tests and Linux equivalence (28 Sep 2026)

## Setup

- **Machine:** the Dell G15 laptop (Windows 11, build 26200; 16 logical cores; RTX 3050).
- **Toolchain:** MSVC 19.44 (Visual Studio 2022 Build Tools 17.14), Windows SDK 10.0.26100 and
  CMake/Ninja from Visual Studio, installed by winget. Everything builds with `/W4 /WX`
  `/permissive- /utf-8 /fp:precise` and the static runtime (`/MT`).
- **Linux reference:** the same commit built with GCC in WSL2 on the same laptop.

## Build and tests

- The whole tree (library, CLI, tests, samaya Studio) compiles with MSVC at `/W4 /WX` with
  **zero warnings**. The only change needed was `strncpy` in the C API, which MSVC deprecates.
- Three problems showed up only with MSVC; all are fixed:
  - **A test helper asked `std::bernoulli_distribution` for a probability above 1**
    (`3.0 / m` for m < 3). That is undefined behaviour; libstdc++ treats it as 1, while MSVC's
    debug library asserts. The helper now caps it at 1, which gives the same matrices on Linux.
  - **The barrier's starting point called `std::clamp` with lo > hi** (`src/lp/ipm.cpp`). When
    a box is narrower than twice the shift, rounding can put `u - margin` one ulp below
    `l + margin` (in about 3% of random boxes). That is undefined behaviour: GCC's library
    tolerated it silently, while MSVC's debug library asserted in every barrier and QP test. The
    upper end of the clamp is now never below the lower end.
  - **MSVC at `/O2` folds the round-trip check `(int64)(double)i == i` to true** even when i is
    not representable. It only mattered to a test's filter, which now takes its exact values
    from the doubles themselves.
- Debug-library assertions now go to stderr instead of a dialog (`tests/test_main.cpp`), so an
  unattended run fails rather than waits.

The test results for the commit are in its message.

## Linux and Windows compute the same answers

- **Method:** both builds solved all 93 Netlib LPs and the six MRPL cases (small and medium),
  with one thread and a 120 s limit (`--json --log-level 0 --threads 1`).
- **Result:** **99 of 99 have the same status, and all 99 are verified on both systems.** 97 give
  the same objective to 1e-9 relative. The other two:

| Model | Linux | Windows | Relative difference | Why |
|---|---:|---:|---:|---|
| mrpl_crude_medium (MILP) | 1,380,665.090 | 1,380,594.942 | 5.1e-05 | Both within the default MIP gap (1e-4) of the optimum; the searches take different paths (for example, `std::sort` orders ties differently in libstdc++ and MSVC's library) |
| pilot87 (LP) | 301.71036566 | 301.71035496 | 3.5e-08 | An ill-conditioned LP: two optimal points within the verifier's 1e-6 |

- **Precision:** the same on both systems.
  - Extended-precision sums (the verifier, and the reported objective) use `Xreal`
    (`src/core/xreal.hpp`), a double-double built from IEEE +, − and × only, with about 106
    significant bits.
  - It replaces `long double`, which was 80-bit with GCC but only 64-bit with MSVC.
  - Floating-point contraction is off on both compilers.
  - The verifier's tolerances are unchanged.
- **Speed:** on these small models (shifted geometric mean of solve time) Linux took 0.016 s and
  Windows 0.019 s. Other jobs were running at the same time, so this is only a rough comparison.

## samaya Studio and setup

- **The installed app:** `setup.bat` installed Studio from the release package, per user and
  without admin rights, to a path with spaces (`C:\Users\HARSHA B S\AppData\Local\Programs\samaya`).
  The installed app then solved the three bundled MRPL samples. Each run was optimal and
  verified, and Studio's independent re-check agrees with `cases/report.py` on the crude case:
  worst row error 2.25e-08, worst bound error 2.96e-08, integrality 0.
- **The Excel export:** the page's workbook writer (`apps/studio/ui/xlsx.js`) produced a workbook
  that Microsoft Excel opens with all sheets and values.
- **`setup.bat /check`:**
  - From the release package, it found Windows 11, the WebView2 Runtime 153 and the RTX 3050.
  - From the source tree, it also found winget, the Build Tools and Python 3.13.
