"""PuLP command-line interface to samaya (PuLP is an optional dependency)."""

import pulp
from pulp.apis.core import PulpSolverError

from samaya_cli import find_binary, read_solution, run_solver


class SAMAYA_CMD(pulp.LpSolver_CMD):
    """Solve LP/MILP models with samaya; expose verified and samaya_result.

    Options are timeLimit (seconds), gapRel, threads, msg, keepFiles, path,
    and options (a list of individual extra CLI arguments).
    """

    name = "SAMAYA_CMD"

    def __init__(self, timeLimit=None, gapRel=None, threads=None, msg=True,
                 keepFiles=False, path=None, options=None):
        self._explicit_path = path
        super().__init__(path=path, keepFiles=keepFiles, msg=msg,
                         timeLimit=timeLimit, gapRel=gapRel, threads=threads,
                         options=options)
        self.verified = False
        self.samaya_result = {}

    def defaultPath(self):
        try:
            return find_binary()
        except FileNotFoundError:
            return "samaya"

    def available(self):
        try:
            return find_binary(self._explicit_path)
        except FileNotFoundError:
            return False

    def actualSolve(self, lp, **kwargs):
        self.verified = False
        self.samaya_result = {}
        modern = hasattr(self, "buildStats")
        if modern:
            from pulp.apis.core import clocks
            start = clocks()
        model, solution = self.create_tmp_files(lp.name, "mps", "sol")
        try:
            binary = find_binary(self._explicit_path)
            exported = lp.writeMPS(model, rename=True, with_objsense=True)
            if modern:
                var_names, row_names, _, original_vars = exported
                variables = dict(zip(var_names, original_vars))
                constraints = dict(zip(row_names, (row.name for row in lp.constraints())))
            else:
                _, variables, constraints, _ = exported
                variables = {v: k for k, v in variables.items()}
                constraints = {v: k for k, v in constraints.items()}
            arguments = []
            for flag, value in (("--time-limit", self.timeLimit),
                                ("--mip-gap", self.optionsDict.get("gapRel")),
                                ("--threads", self.optionsDict.get("threads"))):
                if value is not None:
                    arguments.extend([flag, str(value)])
            if not self.msg:
                arguments.extend(["--log-level", "0"])
            arguments.extend(self.options)
            result = run_solver(binary, model, solution, arguments, self.msg)
            self.samaya_result = result
            self.verified = bool(result.get("verified", False))
            values = read_solution(solution)
            status = result["status"]
            has_solution = bool(values.columns) and (
                status == "optimal" or (self.verified and status in
                    {"time_limit", "iteration_limit", "node_limit", "not_solved"})
            )
            if has_solution:
                lp.assignVarsVals({variables[k]: v for k, v in values.columns.items()
                                   if k in variables})
                lp.assignConsSlack({constraints[k]: v for k, v in values.rows.items()
                                    if k in constraints}, activity=True)
                if not lp.isMIP():
                    lp.assignVarsDj({variables[k]: v for k, v in values.reduced_costs.items()
                                     if k in variables})
                    lp.assignConsPi({constraints[k]: v for k, v in values.duals.items()
                                     if k in constraints})
            # PuLP 4 replaces the legacy status API with LpSolveStats. Preserve
            # conventional problem.status / sol_status for existing callers too.
            legacy = {"optimal": (1, 1), "infeasible": (-1, -1),
                      "unbounded": (-2, -2)}.get(status, (0, 2) if has_solution else (0, 0))
            if status in {"infeasible_or_unbounded", "numerical_error", "invalid_model",
                          "not_implemented", "not_convex"}:
                legacy = (-3, 0)
            if hasattr(lp, "assignStatus"):
                lp.assignStatus(*legacy)
            else:
                lp.status, lp.sol_status = legacy
            if modern:
                states = pulp.constants.LpSolveStatus
                mapped = {"optimal": states.Optimal, "infeasible": states.Infeasible,
                          "unbounded": states.Unbounded, "time_limit": states.TimeLimit,
                          "iteration_limit": states.IterationLimit, "node_limit": states.NodeLimit,
                          "numerical_error": states.NumericalError,
                          "infeasible_or_unbounded": states.Undefined,
                          "invalid_model": states.Undefined, "not_convex": states.Undefined,
                          "not_implemented": states.Undefined}.get(status, states.NotSolved)
                return self.buildStats(lp, mapped, has_solution, start=start,
                                       best_bound=result.get("dual_bound"))
            return legacy[0]
        except (OSError, RuntimeError, ValueError) as exc:
            raise PulpSolverError(f"samaya: {exc}") from exc
        finally:
            self.delete_tmp_files(model, solution)
