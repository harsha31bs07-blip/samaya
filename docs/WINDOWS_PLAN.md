# Plan: native Windows, one-click setup and samaya Studio (from Sun 27 Sep 2026, 23:50 IST)

The owner asked for three things:

1. samaya running on **Windows without WSL**, with **the same precision** as on Linux.
2. A **`.bat` or `.exe` that checks the PC** for everything samaya needs and installs what is
   missing.
3. A **clean, premium GUI**: attach or drag and drop models, see the results **in the app** (not
   in a separate HTML page), with an **Excel download** button.

This is a plan only; nothing here is started. It is meant for after the demo of 28 September.
The rules of `CLAUDE.md` hold throughout:
- differential tests with planted bugs;
- all presets green with zero warnings before every commit;
- the verifier on every result, its tolerances untouched.

## Status (28 Sep 2026, night)

The owner asked for all of it by the morning of 28 September. What was built that night:

**Done:**
- **W1, native Windows build with the same precision:**
  - MSVC presets (`windows-debug`, `windows-release`, `windows-asan`), zero warnings at `/W4 /WX`.
  - `Xreal` (double-double) replaces `long double`.
  - UTF-8 and Unicode paths in the reader and the CLI.
  - A Linux–Windows equivalence run on 99 models, in `docs/results/windows.md`.
- **W3, samaya Studio:** the WebView2 app with drag and drop, results in the app, Excel, CSV and
  solution export, the independent re-check, samples, light and dark themes.
- **CUDA on Windows (added on 28 Sep):**
  - the `windows-cuda` preset;
  - a Compute setting in Studio (CPU or the NVIDIA GPU, for PDLP);
  - `samaya --gpu-info`;
  - fallback to the CPU on unsupported GPUs and on GPU errors (docs/results/windows.md).
- **W4, setup:** `setup.bat` covers both the release package and building from source, with
  `/check`, `/yes` and `/uninstall`. `tools/windows/package.ps1` makes the release zip.

**Deferred:**
- **W2:** the log and progress callbacks, graceful cancel and the C API v2. Studio's Stop ends the
  worker process at once, which D4 makes safe.
- **CI:** the GitHub Actions workflow (`.github/workflows/ci.yml`).
- **An installer `.exe`:** a signed Inno Setup installer. For now the release is a zip with
  `setup.bat`.
- **Reading `.mps.gz` directly.**

**The acceptance criteria (§5), checked that night:**
- **1, developer setup:** `setup.bat /check` works from the source tree. The full source-mode
  install and build on a clean machine was not run (this laptop already had the tools).
- **2, user install:** met with the zip and `setup.bat` instead of `samaya-setup.exe`. Studio
  installed per user, solved the bundled samples, and Excel opened its workbook.
- **3, same precision:**
  - Met: all extended-precision sums use `Xreal` on both systems, and no tolerance changed.
  - The 99-model run shows the same status on both systems, verified on both. Verifier
    outputs were not compared bit for bit on shared solutions (the CLI has no verify-only mode).
- **4, Studio:**
  - Met: several files at once, results in the app, the Excel download.
  - Malformed files show the reader's message with the line number, for example
    `line 6: invalid number 'abc'`.
  - Not measured: the time Stop takes (it ends the worker process with `TerminateProcess`), the
    million-column LP, and a forced solver crash. The out-of-process worker is what keeps the app
    alive and responsive in those cases.
- **5, CI:** there is no CI yet. Every preset was run by hand (Linux and Windows: debug, release,
  asan), and the results are in the commit message.

## 1. What the code already has (checked on 27 Sep)

**Already portable:**
- Standard C++20 only: `std::thread`, no POSIX calls, no `<unistd.h>`.
- No recursion (tree searches use explicit stacks), so Windows' 1 MB default stack is not a risk.
- `long` is never used for sizes, so Windows' 32-bit `long` is not a trap (one test `printf`
  aside).
- The MPS reader already strips `\r` (Windows line endings).
- The one GCC attribute (`src/core/log.hpp`) is guarded by `#if defined(__GNUC__)`.
- `CMakeLists.txt` already has an MSVC branch (`/W4 /permissive-`).
- All time checks go through `Timer` (`src/core/log.hpp`); the only exception is the CLI's read
  timer.

**Gaps:**

| # | Gap | Where |
|---|---|---|
| G1 | Extended precision is `long double`: 80-bit (64-bit significand) with GCC on Linux, but plain `double` (53-bit) with MSVC. | `src/verify/verify.cpp` (`using Real = long double`, 22 uses); objective sums in `src/core/solver.cpp:137,323` and `src/lp/lp_solver.cpp:128` |
| G2 | `read_mps(const std::string&)` opens a narrow path, so file names outside the ANSI code page (e.g. Hindi folder names) fail on Windows. No `.mps.gz` support. | `include/samaya/io.hpp`, `src/io/mps_reader.cpp` |
| G3 | The logger writes straight to stdout. There is no log callback, no progress callback and no way to cancel a solve. | `src/core/log.*`, `include/samaya/params.hpp` |
| G4 | The C API can only read, solve with default parameters and return the objective. It has no parameters, no solution access, no callbacks, and no DLL export macro. | `include/samaya_c.h`, `src/core/c_api.cpp` |
| G5 | Never compiled with MSVC. Expect a `/W4` warning clean-up (mostly `size_t`→`int` and `double`→`float` conversions) before `/WX` holds. | whole tree |
| G6 | The demo and benchmark scripts are bash. Fine for developers (Git Bash or WSL); users get the GUI instead. | `cases/demo.sh`, `bench/*.sh` |
| G7 | On Windows, `nvcc` (CUDA) supports only MSVC as its host compiler, not MinGW. This decides the toolchain. | `src/gpu/*` |

## 2. Decisions (recommended defaults; the owner can change them)

**D1. Toolchain: MSVC (Visual Studio 2022 Build Tools) with CMake and Ninja.**
- It is the only Windows compiler CUDA supports (G7), and the standard one on Windows.
- It has AddressSanitizer (`/fsanitize=address`), and the WebView2 SDK targets it.
- *Rejected:* MinGW-w64 GCC. It would keep the 80-bit `long double` without code changes, but
  it cannot build the GPU code, and precision would still depend on the compiler.

**D2. Precision: a portable double-double accumulator replaces `long double` everywhere.**
- A new `samaya::Xreal` type (`src/core/xreal.hpp`) is a pair of doubles.
- Sums use TwoSum. Products use Dekker's TwoProduct with Veltkamp splitting (constant 2²⁷+1), so
  no FMA is needed.
- That gives a **106-bit significand, more than the 64 bits of today's Linux `long double`**,
  using only `+ − ×` on IEEE doubles.
- The verifier otherwise calls only `std::fabs`, `std::max`, `std::round` and `std::isfinite`,
  all exact. So **verification becomes bit-for-bit identical on Linux and Windows, with any
  compiler**. Linux switches too, so both systems run the same arithmetic.
- FP contraction must be off for these files: `-ffp-contract=off` for GCC/Clang, and `/fp:precise`
  without `/fp:contract` for MSVC.
- A unit test fails if a build ever enables fast-math.

**D3. GUI: "samaya Studio", a native Win32 window hosting Microsoft Edge WebView2.**
- The interface is HTML, CSS and JS embedded in the `.exe`. WebView2 is part of Windows 11 and
  of current Windows 10.
- This gives a premium, modern look with the same design language as the run report and the
  video, plus native file dialogs and drag and drop with real file paths.
- It is one `.exe` with no Python, .NET or Qt runtime, and the solver core stays untouched.
- *Alternatives:*
  - WPF on .NET 8 with a Fluent theme: very native, but a second language and a 70–150 MB
    runtime.
  - Qt 6/QML: a heavy SDK, and Qt's chart modules are GPL.
  - Dear ImGui: tiny, but looks like a developer tool.

**D4. Solves run out of process.**
- Studio starts `samaya.exe` as a worker for each run and talks to it by JSON lines over its
  stdout and stdin.
- A solver crash can never take the app down, a stuck solve can always be stopped, and memory
  goes back to Windows after a big model.
- Studio also uses exactly the CLI binary that the tests check.

**D5. Two ways to install.**
- **`samaya-setup.exe`** for users: prebuilt, no compilers needed, built with Inno Setup.
- **`setup.bat`** for developers: checks for and installs the toolchain, then builds, tests and
  launches.

**D6. GPU on Windows is optional.**
- It is built in the MSVC + CUDA 12 configuration and detected at run time, with a CPU fallback.
- Setup never installs GPU drivers itself.

**D7. Code signing.**
- Unsigned installers trigger Windows SmartScreen ("unknown publisher").
- A code-signing certificate (paid, yearly) is needed before a public release. For the
  hackathon we ship unsigned, with instructions.

## 3. Work packages

### WP-W1: native Windows build with the same precision (2–3 days)

**Creates:**
- `src/core/xreal.hpp`
- `tests/test_xreal.cpp`
- Windows presets in `CMakePresets.json`
- `.github/workflows/ci.yml`
- `docs/results/windows.md`

**Edits (small hunks):**
- `CMakeLists.txt` (MSVC options)
- `src/verify/verify.cpp`
- `src/core/solver.cpp` (two objective sums)
- `src/lp/lp_solver.cpp` (one)
- `include/samaya/io.hpp` and `src/io/mps_reader.cpp` (path type)

1. **Presets** `windows-debug`, `windows-release` and `windows-asan`: Ninja with MSVC, and
   `/W4 /WX /permissive- /utf-8 /EHsc /fp:precise`. Release links the static runtime (`/MT`), so
   users need no VC++ redistributable. Asan uses `/fsanitize=address`; MSVC has no UBSan, but it
   stays on in the Linux asan preset. There is also `windows-cuda`.
2. **Zero warnings under `/W4 /WX`.** Narrowing that is intended gets an explicit `static_cast`
   and a comment; no blanket pragmas.
3. **`Xreal`** replaces `long double` (G1).
   - Huge magnitudes (> 2⁹⁹⁵, where splitting would overflow) are pre-scaled by an exact power
     of two.
   - Tests:
     - exact small-integer dot products;
     - catastrophic cancellation (`1e16 + 1 − 1e16 = 1`);
     - random dot products against GCC's 80-bit `long double` (Linux-only reference) and against
       exact 64-bit integer arithmetic (all platforms);
     - a golden file with the verifier's outputs (as bits) on a fixed corpus (Netlib and the
       MRPL cases), which must match on Windows.
   - Planted bugs:
     - drop TwoSum's error term;
     - compile the tests with contraction on.
4. **Unicode paths (G2):** `read_mps(const std::filesystem::path&)`. The C API takes UTF-8 and
   converts through `std::u8string`. Tests use file names with spaces and Devanagari characters.
5. **Safety margins:** `/STACK:8388608` on the executables (worker threads inherit it), and
   `SetConsoleOutputCP(CP_UTF8)` in the CLI.
6. **CI:** GitHub Actions `windows-latest` (VS 2022 preinstalled) runs the three Windows presets.
   The Linux presets run on every push too.
7. **Equivalence run.** Linux (WSL) and native Windows on the same laptop, over three sets:
   Netlib (93), the MRPL cases (small, medium, large) and the MIPLIB 60 s subset. The report goes
   to `docs/results/windows.md`. Both systems must give:
   - the same statuses and `verified` flags;
   - objectives within 1e-9 relative;
   - bit-identical verifier outputs on the golden corpus.

   Iterates may still differ in the last bits: `exp` and `log` (in PDLP and the IPM) come from
   different maths libraries, and parallel MILP is nondeterministic on Linux already.

**Done when:**
- all six presets (three Linux, three Windows) are green with zero warnings;
- the equivalence report is committed;
- no verifier tolerance has changed.

### WP-W2: log, progress and cancel for applications (1–1.5 days)

**Edits:**
- `include/samaya/params.hpp`
- `src/core/log.{hpp,cpp}`
- `include/samaya/status.hpp` (a new code appended; existing values unchanged)
- `apps/cli/main.cpp`
- `include/samaya_c.h`, `src/core/c_api.cpp`

1. **Log sink:** `Params::log_callback`. The logger calls it when set and prints to stdout
   otherwise.
2. **Progress:** `Params::progress_callback` receives time, phase, primal bound, dual bound, gap,
   nodes and iterations. It is throttled to 10 per second and fed from the existing log points.
3. **Cancel:** a `StopToken` (a shared atomic flag) in `Params`, and `Timer` takes an optional
   pointer to it. Every existing time-limit check reads the clock through `Timer`, so every one
   of them then also stops on cancel, with no new check sites. The new status is `kInterrupted` (C code 12). The best solution so far is
   returned and verified.
4. **CLI `--events`:** JSON lines on stdout (`log`, `progress`, `result`), and `stop` on stdin
   cancels. This is Studio's worker protocol (D4). Large vectors go through the `--solution`
   file, never through JSON.
5. **C API v2** for integrators (C, C#, Python via ctypes):
   - parameter setters;
   - solve with callbacks and cancel;
   - access to x, duals, activities, certificates and the verification numbers;
   - a `SAMAYA_API` export macro, a `samaya.dll` option, and a test written in C.

**Tests:**
- cancel a long MILP and check it stops within 0.2 s;
- every method calls the callbacks;
- the `--events` stream parses and ends with `result`.

**Planted bug:** `Timer` ignores the token, so the cancel test fails.

### WP-W3: samaya Studio, the GUI (4–5 days)

**Creates:**
- `apps/studio/main.cpp`: the Win32 window and WebView2 host.
- `apps/studio/worker.{hpp,cpp}`: runs `samaya.exe --events` and parses the stream.
- `apps/studio/bridge.{hpp,cpp}`: JSON messages to and from the page.
- `apps/studio/recheck.{hpp,cpp}`: the independent re-check, with its own small MPS parser and
  `Xreal` sums. It deliberately shares no code with `src/verify` or `src/io`.
- `apps/studio/xlsx.{hpp,cpp}`: the workbook writer, a C++ port of `cases/report.py`.
- `apps/studio/ui/`: `index.html`, `app.css`, `app.js`, `charts.js` and icons, embedded as
  resources.

**Third-party:**
- the WebView2 SDK (NuGet, pinned version and SHA-256, static loader, so no extra DLL);
- Fluent UI System Icons (MIT), as inline SVG.

**Screens:**

```
┌─ ◆ samaya Studio ─────────────────────────────────────────── ⚙ Settings  ◐ Theme ─┐
│ RUNS              │                                                              │
│ ● crude_small     │     ┌──────────────────────────────────────────────────┐     │
│   ✓ verified 0.1s │     │              ⇩  Drop models here                 │     │
│ ● plan_small      │     │        .mps  ·  .qps      or   [ Browse… ]       │     │
│   ⟳ solving 3.2s  │     └──────────────────────────────────────────────────┘     │
│ ● utility_small   │     Try a sample:  [ Crude schedule ] [ Refinery plan ]     │
│   ⋯ queued        │     Time limit 60 s · Threads auto · MIP gap 0.01% · CPU    │
└───────────────────┴──────────────────────────────────────────────────────────────┘

┌─ MRPL_CRUDE_8k_14d   [OPTIMAL]  [✓ VERIFIED]        [⬇ Excel]  CSV  .sol  PDF ─┐
│ Objective 708,749.06 │ Gap 0.00% │ Time 0.10 s │ 164 × 251 │ 83 integer          │
│ Overview │ Plan │ Verification │ Solution │ Log                                   │
│ (charts and tables, as in the run report, inside the app)                        │
└───────────────────────────────────────────────────────────────────────────────────┘
```

1. **Input:**
   - a drop zone for one or many files, plus Browse (the native file dialog);
   - bundled sample cases;
   - a run queue with live status.

   Drops arrive with their real paths through WebView2's `postMessageWithAdditionalObjects`, so
   there is no copying of large files through JavaScript.
2. **Settings:** time limit, threads, MIP gap, LP method, GPU (greyed out without an NVIDIA
   driver) and presolve. Verification is shown as always on; the GUI cannot turn it off.
3. **During a solve:**
   - elapsed time and phase;
   - for MILP, a live chart of the best plan against the bound;
   - the log tail;
   - **Cancel** (graceful through `stop`; the process is killed after 2 s as a last resort).
4. **Results, in the app, in tabs:**
   - **Overview:** status pills and KPI cards.
   - **Plan:** MRPL-aware charts and tables (berth Gantt, CDU charge, tank stocks, unit feeds,
     marginal value of capacity). Other models get variable and constraint summaries and the
     largest duals.
   - **Verification:** samaya's verifier beside Studio's independent re-check, each check with
     its error against its tolerance. Infeasibility and unboundedness certificates are explained
     in words.
   - **Solution:** virtualized, searchable tables, paged over the bridge so a million rows never
     load at once.
   - **Log.**
5. **Export bar**, always visible:
   - **Download Excel (.xlsx)** as the primary button: the native Save dialog, then the same
     sheets as `cases/report.py`, then a toast with Open and Show in folder;
   - CSV (zip), the solution file, the JSON summary, and PDF (WebView2 `PrintToPdf`).
6. **Design:**
   - Windows 11 Mica title bar (plain on Windows 10);
   - Segoe UI Variable, an 8 px grid, cards with 12 px radius;
   - the palette of the video and report (navy `#1B3A6B`, orange `#F28C28`, green `#2E9E4F`);
   - light and dark themes following Windows;
   - 150–250 ms ease-out motion, per-monitor DPI;
   - keyboard shortcuts (Ctrl+O, Ctrl+Enter, Esc, Ctrl+S);
   - visible focus, AA contrast and screen-reader labels.
7. **History:** each run is kept in `%LOCALAPPDATA%\samaya\runs\<id>\` (parameters, summary,
   solution), so past results reopen instantly.

**Tests:**
- unit tests for the bridge and the worker-stream parser;
- the re-check differentially against `src/verify` and `cases/report.py` on the corpus;
- `.xlsx` structure checks, plus opening in Excel by COM on a developer machine;
- an end-to-end UI script (Playwright attached to WebView2 through its debugging port): drop a
  file, run, expect VERIFIED, export the workbook.

**Planted bug:** the re-check skips one row, and the differential test must fail.

### WP-W4: setup and installer (1.5–2 days)

**Creates:**
- `setup.bat`
- `tools/windows/setup.ps1`
- `tools/windows/components.json` (pinned versions, winget IDs, URLs, SHA-256)
- `installer/samaya.iss` (Inno Setup)
- `docs/WINDOWS.md` (user guide)

**A. `setup.bat` (developers, from source).**
- It starts `setup.ps1` (`-ExecutionPolicy Bypass`), which elevates itself only when an install
  needs it.
- It logs to `%LOCALAPPDATA%\samaya\setup.log`.
- Modes: `setup.bat` (check, install, build, test, launch), `setup.bat /check` (report only),
  `/nogpu`, `/offline <cache>`, `/full` (all three presets).

| Component | Needed for | Detect | Install | Size |
|---|---|---|---|---|
| Windows 10 1809+ or 11, x64 | everything | OS version, architecture | stop with a clear message | – |
| winget (App Installer) | the installs below | `winget --version` | signed `.msixbundle` from Microsoft's winget-cli releases | small |
| VS 2022 Build Tools: MSVC, Windows SDK, CMake, Ninja | building | `vswhere -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64` | `winget install Microsoft.VisualStudio.2022.BuildTools --override "--quiet --wait --norestart --add Microsoft.VisualStudio.Workload.VCTools --includeRecommended --add Microsoft.VisualStudio.Component.VC.CMake.Project"` | several GB |
| Git | getting the source | `git --version` | `winget install Git.Git` | ~300 MB |
| Python 3.12 | case generator, `report.py`, benchmarks | `py -3 --version` | `winget install Python.Python.3.12` | ~100 MB |
| WebView2 Runtime | Studio | registry key `EdgeUpdate\Clients\{F3017226-FE2A-4295-8BDF-00C3A9A7E4C5}` | `winget install Microsoft.EdgeWebView2Runtime` (preinstalled on Windows 11) | small |
| NVIDIA driver | GPU (optional) | `nvidia-smi`, `Win32_VideoController` | **not installed automatically**; setup links to NVIDIA | – |
| CUDA Toolkit 12.x | GPU build (optional) | `nvcc --version`, `CUDA_PATH` | `winget install Nvidia.CUDA`, only with an NVIDIA GPU and after asking | ~3 GB |

After installing, it:
- refreshes `PATH` in the session;
- checks versions (CMake ≥ 3.25);
- enters the VS developer environment (found by `vswhere`);
- configures and builds `windows-release`, and runs its tests (all three presets with `/full`);
- generates the sample cases;
- creates the Start-menu and desktop shortcuts, and opens Studio.

Safety:
- it asks before every large download and shows its size;
- it installs only from winget's signed manifests or official URLs with SHA-256 and Authenticode
  checks;
- it never installs drivers, and it is idempotent (a re-run skips what is present);
- it resumes after a reboot: the VS installer's exit code 3010 schedules the continuation
  through `RunOnce`.

**B. `samaya-setup.exe` (users, no compilers).** CI builds it with Inno Setup from the
`windows-release` artifacts.
- Contents:
  - `samaya.exe` (CLI and Studio's worker) and `samaya-studio.exe`;
  - `samaya.dll` with its headers, for integrators;
  - the sample cases, the docs and the licence.
- Checks:
  - Windows version and architecture, and disk space;
  - the **WebView2 Runtime**: if missing, it silently runs the 2 MB Evergreen bootstrapper (the
    offline variant bundles the standalone installer);
  - whether an **NVIDIA driver** is present, which enables the GPU option;
  - no VC++ runtime is needed (static `/MT`).
- Installs per user (no admin) by default, with an all-users option.
- Adds Start-menu and desktop shortcuts, an optional file association ("Open .mps with samaya
  Studio"), optional `PATH`, and an uninstaller.
- A portable `.zip` ships as well.

### WP-W5: testing and release (1 day)

**The QA matrix:**
- Windows 11 (23H2, 24H2) and Windows 10 22H2, on clean machines. Windows Sandbox is enough for
  installer tests on Pro editions; Home needs a Hyper-V or VirtualBox VM.
- A standard user without admin rights, and no internet (the offline installer).
- Paths with spaces and non-ASCII characters.
- 100–200% display scaling, light and dark themes.
- With and without an NVIDIA GPU.
- The million-column refinery LP (memory, and Studio staying responsive).
- Cancel during a solve.
- Malformed MPS files (the error must give the line number); infeasible and unbounded models
  (certificates shown).
- SmartScreen and antivirus behaviour of unsigned binaries.

**Documentation:** `docs/WINDOWS.md` with screenshots, the README, and the slide line changed to
"Windows (native) and Linux".

## 4. Schedule

| Days (one developer) | Work | With three people |
|---|---|---|
| 1–3 | WP-W1: MSVC build, `Xreal`, Unicode paths, CI, equivalence run | A: W1, then W2 |
| 3–4 | WP-W2: callbacks, cancel, `--events`, C API v2 | B: W3, starting on a mock worker stream |
| 4–8 | WP-W3: Studio | C: W4, then QA |
| 8–9 | WP-W4: `setup.bat`, `samaya-setup.exe` | |
| 10 | WP-W5: QA matrix, docs, release | all |

**Total:** about 10 working days alone, or 4–5 days with three people.

**Critical path:** the MSVC build (W1), then the worker protocol (W2), then Studio integration
(W3).

## 5. Acceptance criteria

1. **Developer setup:** `setup.bat` on a clean Windows 11 machine installs what is missing (with
   consent), builds, passes the tests and opens Studio. The only manual steps are UAC prompts.
2. **User install:** `samaya-setup.exe` installs as a standard user in under 2 minutes; Studio
   then solves a bundled sample, and its Excel file opens in Excel.
3. **Same precision:**
   - verifier outputs are bit-identical on Linux and Windows over the golden corpus;
   - all extended-precision sums are double-double (106-bit) on both systems;
   - no tolerance has changed.
4. **Studio:**
   - several files can be dropped at once, and results show in the app;
   - the Excel download works;
   - cancel takes at most 0.5 s;
   - the app stays responsive during the million-column LP;
   - a solver crash never closes the app.
5. **CI:** every preset (Linux and Windows: debug, release, asan) is green with zero warnings.
   Planted-bug checks are recorded in the commit messages.

## 6. Risks

- **The MSVC warning clean-up may be larger than expected.**
  - Fix it module by module. `/WX` stays on (the zero-warning rule), so budget the time rather
    than relax it.
- **MSVC's code can be slower than GCC's on numeric loops.**
  - Measure it in the equivalence run.
  - The fallback is clang-cl, which ships with VS: same ABI, and CUDA still compiles with `cl`.
- **Some machines lack parts of the tooling.**
  - Old Windows 10 builds without WebView2 are handled by the installer's bootstrapper.
  - LTSC editions without winget use the direct-download fallback.
- **Unsigned binaries trigger SmartScreen and antivirus heuristics.**
  - Get a code-signing certificate before any public release (D7).
- **GPU on Windows adds WDDM launch overhead and timeout detection.** Our kernels are short, but
  it needs measuring. The GPU stays optional (D6).

## 7. Not in this plan

- macOS and Windows on ARM64. `Xreal` makes both easy later.
- Reading `.mps.gz` directly. This needs our own inflate and is useful for MIPLIB files; it goes
  into W1 only if time allows.

## 8. Decisions for the owner

1. **GUI technology:** WebView2 with C++ (recommended), or WPF on .NET 8.
2. **GPU in the first Windows release:** now, or in the next one.
3. **What ships first:** the user installer (`samaya-setup.exe`, recommended for demos) or the
   developer `setup.bat`. Both are in the plan.
4. **A code-signing certificate** before a public release.
