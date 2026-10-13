"""Import this module to register SolverFactory('samaya') with Pyomo."""

import re

from pyomo.common.collections import Bunch
from pyomo.common.tempfiles import TempfileManager
from pyomo.opt import (ProblemFormat, ResultsFormat, SolverFactory, SolverResults,
                       SolverStatus, TerminationCondition, SolutionStatus)
from pyomo.opt.solver import SystemCallSolver

from samaya_cli import command_line, find_binary, parse_summary, read_solution


@SolverFactory.register("samaya", doc="samaya verified LP/MILP shell solver")
class SAMAYA(SystemCallSolver):
    """MPS shell solver; results.solver.samaya holds samaya's JSON summary.

    Set options['time_limit'], options['mip_gap'], options['threads'], or
    other CLI option names using underscores. Standard tee, keepfiles and
    executable arguments are supported.
    """

    def __init__(self, **kwargs):
        kwargs["type"] = "samaya"
        super().__init__(**kwargs)
        self._valid_problem_formats = [ProblemFormat.mps]
        self._valid_result_formats = {ProblemFormat.mps: [ResultsFormat.soln]}
        self.set_problem_format(ProblemFormat.mps)
        self._capabilities = Bunch(linear=True, integer=True)

    def _default_executable(self):
        try:
            return find_binary()
        except FileNotFoundError:
            return None

    def available(self, exception_flag=False):
        available = super().available(exception_flag=False)
        if not available and exception_flag:
            raise RuntimeError("samaya executable not found; set SAMAYA_BIN or executable=")
        return available

    def create_command_line(self, executable, problem_files):
        self._soln_file = TempfileManager.create_tempfile(suffix=".samaya.sol")
        if self._log_file is None:
            self._log_file = TempfileManager.create_tempfile(suffix=".samaya.log")
        arguments = []
        for key, value in self.options.items():
            arguments.append("--" + key.replace("_", "-"))
            if value is not None and value != "":
                arguments.append(str(value))
        if self._timelimit is not None:
            arguments.extend(["--time-limit", str(self._timelimit)])
        if not self._tee and "log_level" not in self.options:
            arguments.extend(["--log-level", "0"])
        return Bunch(cmd=command_line(executable, problem_files[0], self._soln_file,
                                     arguments), log_file=self._log_file, env=None)

    def process_logfile(self):
        summary = parse_summary(self._log)
        results = SolverResults()
        results.solver.samaya = summary
        verified = "verified" if summary.get("verified") else "not verified"
        results.solver.message = f"{verified}: {summary.get('message', '')}"
        status = summary["status"]
        conditions = {
            "optimal": TerminationCondition.optimal,
            "infeasible": TerminationCondition.infeasible,
            "unbounded": TerminationCondition.unbounded,
            "infeasible_or_unbounded": TerminationCondition.infeasibleOrUnbounded,
            "time_limit": TerminationCondition.maxTimeLimit,
            "iteration_limit": TerminationCondition.maxIterations,
            "node_limit": TerminationCondition.maxEvaluations,
            "invalid_model": TerminationCondition.invalidProblem,
            "not_convex": TerminationCondition.invalidProblem,
            "not_implemented": TerminationCondition.error,
            "numerical_error": TerminationCondition.solverFailure,
            "not_solved": TerminationCondition.other,
        }
        results.solver.termination_condition = conditions.get(status, TerminationCondition.error)
        if status in {"optimal", "infeasible", "unbounded", "infeasible_or_unbounded"}:
            results.solver.status = SolverStatus.ok
        elif status in {"time_limit", "iteration_limit", "node_limit", "not_solved"}:
            results.solver.status = SolverStatus.warning
        else:
            results.solver.status = SolverStatus.error
        results.problem.number_of_variables = summary.get("cols", 0)
        results.problem.number_of_constraints = summary.get("rows", 0)
        return results

    def process_soln_file(self, results):
        summary = results.solver.samaya
        values = read_solution(self._soln_file)
        status = summary["status"]
        if not values.columns or not (status == "optimal" or
                (summary.get("verified") and status in
                 {"time_limit", "iteration_limit", "node_limit", "not_solved"})):
            return
        solution = results.solution.add()
        solution.status = SolutionStatus.optimal if status == "optimal" else SolutionStatus.feasible
        solution.objective["__default_objective__"] = {"Value": summary.get("objective")}
        if status == "optimal":
            results.problem.lower_bound = results.problem.upper_bound = summary.get("objective")
        dual = any(re.fullmatch(pattern, "dual") for pattern in self._suffixes)
        rc = any(re.fullmatch(pattern, "rc") for pattern in self._suffixes)
        for name, value in values.columns.items():
            solution.variable[name] = {"Value": value}
            if rc and name in values.reduced_costs and not summary.get("integers"):
                solution.variable[name]["Rc"] = values.reduced_costs[name]
        if dual and not summary.get("integers"):
            # Pyomo emits two aliases for ranged constraints. Both map to the
            # same component, so combine their multipliers before loading.
            ranged = {}
            for name, value in values.duals.items():
                if name.startswith(("r_l_", "r_u_")):
                    key = "r_l_" + name[4:]
                    ranged[key] = ranged.get(key, 0.0) + value
                else:
                    solution.constraint[name] = {"Dual": value}
            for name, value in ranged.items():
                solution.constraint[name] = {"Dual": value}
