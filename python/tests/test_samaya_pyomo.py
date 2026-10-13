"""Pyomo integration tests against the real samaya CLI."""

import itertools
from pathlib import Path
import sys

import pytest

pytest.importorskip("pyomo")
import pyomo.environ as pyo
from pyomo.opt import SolverStatus, TerminationCondition

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from samaya_cli import find_binary
import samaya_pyomo

try:
    BINARY = find_binary()
except FileNotFoundError as exc:
    pytest.skip(str(exc), allow_module_level=True)


def solve(model, **kwargs):
    solver = pyo.SolverFactory("samaya", executable=BINARY)
    solver.options.update({"time_limit": 30, "threads": 1, "mip_gap": 0})
    return solver.solve(model, **kwargs)


@pytest.mark.parametrize("maximize", [True, False])
def test_wyndor(maximize):
    model = pyo.ConcreteModel()
    model.x = pyo.Var(domain=pyo.NonNegativeReals)
    model.y = pyo.Var(domain=pyo.NonNegativeReals)
    model.obj = pyo.Objective(expr=(1 if maximize else -1)*(3*model.x + 5*model.y),
                              sense=pyo.maximize if maximize else pyo.minimize)
    model.first = pyo.Constraint(expr=model.x <= 4)
    model.second = pyo.Constraint(expr=2*model.y <= 12)
    model.third = pyo.Constraint(expr=3*model.x + 2*model.y <= 18)
    model.dual = pyo.Suffix(direction=pyo.Suffix.IMPORT)
    model.rc = pyo.Suffix(direction=pyo.Suffix.IMPORT)
    results = solve(model)
    assert results.solver.status == SolverStatus.ok
    assert results.solver.termination_condition == TerminationCondition.optimal
    assert results.solver.samaya["verified"]
    assert "verified" in results.solver.message
    assert pyo.value(model.obj) == pytest.approx(36 if maximize else -36, abs=1e-6)
    assert [pyo.value(model.x), pyo.value(model.y)] == pytest.approx([2, 6], abs=1e-6)
    sign = 1 if maximize else -1
    assert [model.dual[row] for row in (model.first, model.second, model.third)] == pytest.approx(
        [0, sign*1.5, sign], abs=1e-6)
    assert [model.rc[v] for v in (model.x, model.y)] == pytest.approx([0, 0], abs=1e-6)


def test_knapsack():
    profits, weights = [10, 13, 7, 8], [5, 7, 4, 3]
    optimum = max(sum(p*v for p, v in zip(profits, point))
                  for point in itertools.product((0, 1), repeat=4)
                  if sum(w*v for w, v in zip(weights, point)) <= 14)
    assert optimum == 28
    model = pyo.ConcreteModel()
    model.x = pyo.Var(range(4), domain=pyo.Binary)
    model.obj = pyo.Objective(expr=sum(profits[i]*model.x[i] for i in range(4)), sense=pyo.maximize)
    model.capacity = pyo.Constraint(expr=sum(weights[i]*model.x[i] for i in range(4)) <= 14)
    results = solve(model)
    assert results.solver.termination_condition == TerminationCondition.optimal
    assert results.solver.samaya["verified"]
    assert pyo.value(model.obj) == pytest.approx(optimum, abs=1e-6)
    assert [pyo.value(model.x[i]) for i in range(4)] == pytest.approx([0, 1, 1, 1], abs=1e-6)


def test_infeasible():
    model = pyo.ConcreteModel()
    model.x = pyo.Var(domain=pyo.NonNegativeReals)
    model.obj = pyo.Objective(expr=model.x)
    model.lower = pyo.Constraint(expr=model.x >= 2)
    model.upper = pyo.Constraint(expr=model.x <= 1)
    results = solve(model)
    assert results.solver.termination_condition == TerminationCondition.infeasible
    assert results.solver.status == SolverStatus.ok


@pytest.mark.parametrize("maximize", [True, False])
def test_nonzero_reduced_cost(maximize):
    model = pyo.ConcreteModel()
    model.x = pyo.Var(domain=pyo.NonNegativeReals)
    model.y = pyo.Var(domain=pyo.NonNegativeReals)
    sign = 1 if maximize else -1
    model.obj = pyo.Objective(expr=sign*(2*model.x + model.y),
                              sense=pyo.maximize if maximize else pyo.minimize)
    model.capacity = pyo.Constraint(expr=model.x + model.y <= 2)
    model.rc = pyo.Suffix(direction=pyo.Suffix.IMPORT)
    solve(model)
    assert [pyo.value(model.x), pyo.value(model.y)] == pytest.approx([2, 0], abs=1e-6)
    assert [model.rc[model.x], model.rc[model.y]] == pytest.approx([0, -sign], abs=1e-6)


def test_names_and_ranged_duals():
    model = pyo.ConcreteModel()
    model.x = pyo.Var([(1, 2)], domain=pyo.NonNegativeReals)
    model.add_component("other-var", pyo.Var(domain=pyo.NonNegativeReals))
    y = model.component("other-var")
    model.obj = pyo.Objective(expr=2*model.x[1, 2] + y, sense=pyo.maximize)
    model.add_component("row-with-dash", pyo.Constraint(expr=(0, model.x[1, 2], 2)))
    model.rows = pyo.Constraint([(1, 2)], rule=lambda m, i, j: y <= 3)
    model.dual = pyo.Suffix(direction=pyo.Suffix.IMPORT)
    results = solve(model, symbolic_solver_labels=True)
    assert results.solver.samaya["verified"]
    assert [pyo.value(model.x[1, 2]), pyo.value(y)] == pytest.approx([2, 3], abs=1e-6)
    assert model.dual[model.component("row-with-dash")] == pytest.approx(2, abs=1e-6)
    assert model.dual[model.rows[1, 2]] == pytest.approx(1, abs=1e-6)
