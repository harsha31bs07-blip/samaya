# samaya Python solver plugins

Install PuLP or Pyomo in your Python environment (`pip install pulp pyomo`),
and add this directory to `PYTHONPATH`; these modules require no other runtime
dependencies. The plugins target PuLP 4.0 and Pyomo 6.10. Each framework is
optional: import only its plugin. Install/build the samaya CLI separately.

Binary lookup: explicit `path=` (PuLP) / `executable=` (Pyomo), then `SAMAYA_BIN`,
then `samaya` on `PATH`, then this repository's `build/{windows-release,
windows-debug,release,debug}/apps/cli/` directories.

For an existing PuLP model:

```python
from samaya_pulp import SAMAYA_CMD

solver = SAMAYA_CMD(msg=False)
stats = model.solve(solver)
print(solver.verified, solver.samaya_result)
```

For an existing Pyomo model:

```python
import samaya_pyomo  # registers the shell plugin
from pyomo.opt import SolverFactory

results = SolverFactory("samaya").solve(model)
print(results.solver.message, results.solver.samaya)
```

PuLP options: `timeLimit` (seconds), `gapRel` (relative MIP gap), `threads`,
`msg`, `keepFiles`, `path`, and `options`, a list of CLI argument tokens
(e.g. `["--lp-method", "dual", "--no-presolve"]`). PuLP 4 returns `LpSolveStats`;
the plugin also sets legacy numeric `model.status` and `model.sol_status`.
A limit with a verified incumbent has legacy status not solved (0), solution
status integer feasible (2). Raw status, objective and bound are in
`solver.samaya_result`; verification is in `solver.verified`.

Pyomo options: `solver.options["time_limit"]`, `["mip_gap"]`, `["threads"]`,
and other CLI names with underscores (e.g. `lp_method`, `no_presolve=None`).
`solve(tee=True, keepfiles=True)` uses Pyomo's usual logging/file lifecycle;
`solve(timelimit=30)` also forwards the limit. Verification and the full JSON
summary are stored in `results.solver.samaya`; the message says "verified"
or "not verified". LP duals and reduced costs load into import-direction
`model.dual` and `model.rc` suffixes. PuLP LPs receive `constraint.pi` and
`variable.dj`. Both plugins preserve original objective sense and dual signs
via MPS `OBJSENSE`; MILPs do not receive LP duals. No solverModel is created.

Run tests from the repository root (pytest is only a test dependency):

```powershell
$env:PYTHONPATH = "$PWD\python"
$env:SAMAYA_BIN = 'C:\dev\samaya\mrpl-pro\build\windows-debug\apps\cli\samaya.exe'
C:\dev\samaya\.venv-py\Scripts\python.exe -m pytest python\tests -v
```

Tests skip a module when its optional framework is missing, or when the CLI
cannot be located. They cover min/max LP objectives, dual signs, primal name
mapping, infeasibility, and MILP optimality against exhaustive enumeration.
