"""PuLP integration tests against the real samaya CLI."""

import itertools
from pathlib import Path
import sys

import pytest

pulp = pytest.importorskip("pulp")
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from samaya_cli import find_binary
from samaya_pulp import SAMAYA_CMD

try:
    BINARY = find_binary()
except FileNotFoundError as exc:
    pytest.skip(str(exc), allow_module_level=True)


def solve(model):
    solver = SAMAYA_CMD(path=BINARY, msg=False, timeLimit=30, threads=1, gapRel=0)
    model.solve(solver)
    return solver


@pytest.mark.parametrize("maximize", [True, False])
def test_wyndor(maximize):
    model = pulp.LpProblem("wyndor", pulp.LpMaximize if maximize else pulp.LpMinimize)
    x, y = model.add_variable("x", lowBound=0), model.add_variable("y", lowBound=0)
    model += (1 if maximize else -1) * (3*x + 5*y)
    model += x <= 4, "first"
    model += 2*y <= 12, "second"
    model += 3*x + 2*y <= 18, "third"
    solver = solve(model)
    assert model.status == 1
    assert model.sol_status == 1
    assert solver.samaya_result["status"] == "optimal"
    assert solver.verified
    assert pulp.value(model.objective) == pytest.approx(36 if maximize else -36, abs=1e-6)
    assert [x.value(), y.value()] == pytest.approx([2, 6], abs=1e-6)
    # Shadow prices differentiate the original objective with respect to RHS.
    sign = 1 if maximize else -1
    assert [row.pi for row in model.constraints()] == pytest.approx(
        [0, sign*1.5, sign], abs=1e-6)
    assert [x.dj, y.dj] == pytest.approx([0, 0], abs=1e-6)


def test_knapsack():
    profits, weights = [10, 13, 7, 8], [5, 7, 4, 3]
    optimum = max(sum(p*v for p, v in zip(profits, point))
                  for point in itertools.product((0, 1), repeat=4)
                  if sum(w*v for w, v in zip(weights, point)) <= 14)
    assert optimum == 28
    model = pulp.LpProblem("knapsack", pulp.LpMaximize)
    variables = [model.add_variable(name, cat="Binary") for name in "abcd"]
    model += pulp.lpSum(p*v for p, v in zip(profits, variables))
    model += pulp.lpSum(w*v for w, v in zip(weights, variables)) <= 14
    solver = solve(model)
    assert model.status == 1
    assert solver.verified
    assert pulp.value(model.objective) == pytest.approx(optimum, abs=1e-6)
    assert [v.value() for v in variables] == pytest.approx([0, 1, 1, 1], abs=1e-6)


def test_infeasible():
    model = pulp.LpProblem("infeasible")
    x = model.add_variable("x", lowBound=0)
    model += x
    model += x >= 2
    model += x <= 1
    solver = solve(model)
    assert model.status == -1
    assert model.sol_status == -1
    assert solver.samaya_result["status"] == "infeasible"


@pytest.mark.parametrize("maximize", [True, False])
def test_nonzero_reduced_cost(maximize):
    model = pulp.LpProblem("reduced_cost", pulp.LpMaximize if maximize else pulp.LpMinimize)
    x = model.add_variable("x", lowBound=0)
    y = model.add_variable("y", lowBound=0)
    sign = 1 if maximize else -1
    model += sign * (2*x + y)
    model += x + y <= 2
    solve(model)
    assert [x.value(), y.value()] == pytest.approx([2, 0], abs=1e-6)
    assert [x.dj, y.dj] == pytest.approx([0, -sign], abs=1e-6)


def test_names_and_keepfiles(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    model = pulp.LpProblem("names", pulp.LpMaximize)
    x = model.add_variable("x[1,2]", lowBound=0)
    y = model.add_variable("other-var", lowBound=0)
    model += 2*x + y
    model += x <= 2, "row[1,2]"
    model += y <= 3, "row-with-dash"
    solver = SAMAYA_CMD(path=BINARY, msg=False, keepFiles=True,
                        options=["--lp-method", "dual", "--no-presolve"])
    model.solve(solver)
    assert solver.verified
    assert [x.value(), y.value()] == pytest.approx([2, 3], abs=1e-6)
    assert [row.pi for row in model.constraints()] == pytest.approx([2, 1], abs=1e-6)
    assert (tmp_path / "names-pulp.mps").exists()
    assert (tmp_path / "names-pulp.sol").exists()
