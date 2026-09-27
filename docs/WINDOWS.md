# samaya on Windows

samaya runs natively on Windows 10 (1809 or later) and Windows 11, 64-bit; WSL is not needed.
The release has two programs:

- **samaya Studio** (`samaya-studio.exe`), the desktop app: drop models, see the results, download
  them as an Excel workbook.
- **samaya** (`samaya.exe`), the command-line solver, the same as on Linux.

## Install (release package)

1. Unzip `samaya-windows-<version>.zip`.
2. Double-click **`setup.bat`**.

`setup.bat` checks the PC, installs what is missing (asking first), and installs samaya Studio
for the current user, with no administrator rights needed:

- Windows version and architecture.
- The Microsoft Edge WebView2 Runtime, which Studio uses to draw its window. It is part of
  Windows 11 and current Windows 10; if it is missing, setup installs it through winget or
  Microsoft's signed installer.
- An NVIDIA GPU, reported only. This release runs on the CPU.

Studio is installed to `%LOCALAPPDATA%\Programs\samaya`, with shortcuts in the Start menu and
on the desktop. The installer has three more modes:

- `setup.bat /check` reports what is present and changes nothing.
- `setup.bat /yes` installs without asking.
- `setup.bat /uninstall` removes Studio.

The log is written to `%LOCALAPPDATA%\samaya\setup.log`.

## Using samaya Studio

Studio is laid out like other Windows planning tools: a toolbar with the commands, the list of
**cases** on the left, and a status bar. Each solve is a case: a model file, its settings and its
result.

1. **New case.** Drag MPS or QPS files onto the window, use **Open models** (Ctrl+O), or press
   **Solve** next to one of the bundled MRPL-style samples. Several files are solved one after
   another.
2. **Settings:**
   - **Start date:** day 1 of the model's horizon. The crude schedule's days, the plan's weeks
     and the utility case's hours are shown as calendar dates from there, in the charts, the
     tables and the Excel file.
   - Time limit, threads, MIP gap, LP method and presolve.
   - Verification is always on.
3. **Results** appear in the app:
   - **Overview:** the status, whether the answer is verified, the objective, the bound and gap,
     the horizon on the calendar, the time, and the case's settings.
   - **Schedule** (crude and utility cases) or **Plan** (planning case):
     - the berth schedule as a Gantt chart on dates, with weekends shaded;
     - tank stocks, CDU charge, crude processed and unit feeds by date or week;
     - boilers and the gas turbine on and off by hour;
     - the marginal value of each binding capacity.
     Other models get their largest values and duals.
   - **Verification:** samaya's own verifier next to an **independent re-check**. Studio reads the
     model file with its own code and recomputes every row, bound and the objective from the
     solution.
   - **Solution:** every variable and constraint, with a filter.
   - **Log:** the solver's output.
4. **Export:**
   - **Export to Excel** (Ctrl+S) saves a workbook. It has sheets for the summary, the
     verification, each schedule or plan table (dates as real Excel dates), all variables, all
     constraints and the log.
   - **CSV** saves the same tables as a zip of CSV files.
   - **Solution file** saves samaya's `.sol` file.
5. **Run again with other settings:**
   - **Back** (Alt+Left), above a case's name, or **Run again** in the toolbar, returns to the
     settings with that case's values filled in.
   - Change any of them and press **Run**. The new case is listed as, for example,
     `mrpl_plan_small (2)`.
6. **Compare cases:** with two or more finished cases of the same model, **Compare cases**
   shows them side by side, as planners compare a base case with alternatives. It shows the
   status, objective and change from the base case, the time, and the settings (those that
   differ are in bold), and then the variables that differ most. Export to Excel saves the
   comparison.

**Stop** (Esc) ends a running solve at once. Every solve runs in its own `samaya.exe` process,
so a problem in one solve cannot close the app. The **Dark theme** switch in the status bar
changes the colours; Studio keeps the choice.

## Command line

```bat
samaya.exe --json --solution plan.sol model.mps
samaya.exe --help
```

File names may contain spaces and any characters; on Windows they are opened through the
Unicode file APIs.

## Precision

Linux and Windows compute the same answers, verified the same way:

- Extended-precision sums (in the verifier and in the reported objective) use a double-double
  type, `src/core/xreal.hpp`, with about 106 significant bits.
- It is built from IEEE double +, − and × only, so the arithmetic does not depend on the
  compiler. `long double`, which it replaces, is 80-bit with GCC but plain double with Microsoft's
  compiler.
- Floating-point contraction is off on every compiler (`-ffp-contract=off`, `/fp:precise`).
- The verifier's tolerances are the same on both systems.

## Building from source

Double-click `setup.bat` in a source checkout. It checks for, and offers to install:
- Visual Studio 2022 Build Tools (the C++ compiler, Windows SDK, CMake and Ninja);
- Python 3 (for the sample cases);
- the WebView2 Runtime.

It then builds with MSVC, runs the tests, and installs Studio. By hand:

```bat
powershell -ExecutionPolicy Bypass -File tools\windows\build.ps1 -Preset windows-release -Test
powershell -ExecutionPolicy Bypass -File tools\windows\package.ps1
```

The presets are `windows-debug`, `windows-release` and `windows-asan` (AddressSanitizer). All
build with `/W4 /WX`.
