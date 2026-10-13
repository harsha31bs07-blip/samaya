#!/usr/bin/env python3
"""Independent, standard-library checks of pooling MPS files and feasible starts.

Run from any directory with: python tests/pooling_check.py
This checks feasibility and objective values, not stationarity or global optimality.
"""

from __future__ import annotations

import json
import math
import shutil
import subprocess
import sys
import uuid
from contextlib import contextmanager
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TOLERANCE = 1e-9  # Absolute tolerance requested for the analytical Haverly fixtures.


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def read_mps(path: Path) -> dict:
    """Read only ROWS, COLUMNS, RHS, BOUNDS and full-symmetric QCMATRIX."""
    model = {"rows": {}, "columns": {}, "rhs": {}, "bounds": {}, "quadratic": {}}
    section = None
    quadratic_row = None
    ended = False
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if not fields or fields[0].startswith("*"):
            continue
        head = fields[0]
        if head == "NAME":
            continue
        if head == "ENDATA":
            ended = True
            break
        if head in ("ROWS", "COLUMNS", "RHS", "BOUNDS", "QCMATRIX"):
            section = head
            if head == "QCMATRIX":
                require(len(fields) == 2, "QCMATRIX needs a row name")
                quadratic_row = fields[1]
                require(quadratic_row in model["rows"], "unknown quadratic row")
                require(model["rows"][quadratic_row] != "N", "quadratic objective unsupported")
                require(quadratic_row not in model["quadratic"], "duplicate QCMATRIX block")
                model["quadratic"][quadratic_row] = {}
            continue
        if section == "ROWS":
            kind, row = fields
            require(kind in ("N", "E", "L", "G"), f"unsupported row sense {kind}")
            require(row not in model["rows"], f"duplicate row {row}")
            model["rows"][row] = kind
        elif section in ("COLUMNS", "RHS"):
            require(len(fields) >= 3 and len(fields) % 2 == 1, "invalid row/value pairs")
            name = fields[0]
            if section == "COLUMNS":
                entries = model["columns"].setdefault(name, {})
                model["bounds"].setdefault(name, [0.0, math.inf])
            else:
                entries = model["rhs"]
            for row, value in zip(fields[1::2], fields[2::2]):
                require(row in model["rows"], f"unknown row {row}")
                entries[row] = entries.get(row, 0.0) + float(value)
        elif section == "BOUNDS":
            kind, _, col, *values = fields
            require(col in model["columns"], f"unknown bound column {col}")
            require(kind in ("LO", "UP", "FX"), f"unsupported bound {kind}")
            require(len(values) == 1, "bound needs one value")
            value = float(values[0])
            if kind in ("LO", "FX"):
                model["bounds"][col][0] = value
            if kind in ("UP", "FX"):
                model["bounds"][col][1] = value
        elif section == "QCMATRIX":
            first, second, value = fields
            require(first in model["columns"] and second in model["columns"],
                    "unknown quadratic column")
            entries = model["quadratic"][quadratic_row]
            require((first, second) not in entries, "duplicate quadratic entry")
            entries[first, second] = float(value)
        else:
            raise AssertionError(f"unsupported MPS content: {line}")
    require(ended, "missing ENDATA")
    require(sum(kind == "N" for kind in model["rows"].values()) == 1,
            "expected one objective row")
    for row, entries in model["quadratic"].items():
        for (first, second), value in entries.items():
            require(math.isfinite(value), f"nonfinite matrix entry in {row}")
            require(entries.get((second, first)) == value, f"missing/unequal mirror in {row}")
    return model


def evaluate(model: dict, point: dict[str, float], expected: float | None = None) -> float:
    require(set(point) == set(model["columns"]), "point must specify every column exactly")
    for col, value in point.items():
        lo, up = model["bounds"][col]
        require(math.isfinite(value) and lo - TOLERANCE <= value <= up + TOLERANCE,
                f"bound violation: {col} = {value}, bounds [{lo}, {up}]")
    objective = None
    for row, kind in model["rows"].items():
        # Sum every Q entry separately: off-diagonal pairs already contain c/2 each.
        terms = [entries.get(row, 0.0) * point[col]
                 for col, entries in model["columns"].items()]
        terms.extend(value * point[first] * point[second]
                     for (first, second), value in model["quadratic"].get(row, {}).items())
        activity = math.fsum(terms)
        rhs = model["rhs"].get(row, 0.0)
        require(math.isfinite(activity), f"nonfinite activity in {row}")
        if kind == "N":
            objective = activity
        else:
            violation = {"E": abs(activity - rhs), "L": activity - rhs,
                         "G": rhs - activity}[kind]
            require(violation <= TOLERANCE,
                    f"row violation: {row}, activity {activity}, {kind} RHS {rhs}")
    if expected is not None:
        require(abs(objective - expected) <= TOLERANCE,
                f"objective {objective}, expected {expected}")
    return objective


def generate(*args: str) -> None:
    result = subprocess.run([sys.executable, str(ROOT / "cases" / "pooling.py"), *args],
                            capture_output=True, text=True)
    require(result.returncode == 0, f"generator failed: {result.stdout}{result.stderr}")


@contextmanager
def scratch_directory():
    # mkdir's default mode inherits Windows permissions; tempfile's 0700 mode can
    # deny subprocess access under restricted Windows tokens (Python 3.13).
    path = ROOT / "tests" / f"samaya_pooling_{uuid.uuid4().hex}"
    path.mkdir()
    try:
        yield path
    finally:
        require(path.resolve().parent == (ROOT / "tests").resolve(), "unsafe scratch path")
        shutil.rmtree(path)


def main() -> None:
    fixtures = ROOT / "tests" / "instances"
    points = json.loads((fixtures / "haverly_points.json").read_text(encoding="utf-8"))
    require(points["columns"] == ["a", "b", "u", "v", "z1", "z2", "p"],
            "incorrect fixture column order")
    require(set(points["instances"]) == {"haverly1", "haverly2", "haverly3"},
            "expected all three Haverly fixtures")
    with scratch_directory() as temporary:
        for number in (1, 2, 3):
            name = f"haverly{number}"
            model = read_mps(fixtures / f"{name}.mps")
            require(list(model["columns"]) == points["columns"], "MPS column order changed")
            require(list(model["rows"]) == ["obj", "pool_balance", "pool_quality",
                                           "prod1_quality", "prod2_quality",
                                           "prod1_demand", "prod2_demand"],
                    "Haverly row names/order changed")
            data = points["instances"][name]
            for point_key, value_key in (("global_solution", "global_min"),
                                         ("trap_start", "trap_start_value"),
                                         ("local_point", "local_value")):
                values = data[point_key]
                require(len(values) == len(points["columns"]), "wrong point dimension")
                point = dict(zip(points["columns"], values))
                objective = evaluate(model, point, data[value_key])
                print(f"{name} {point_key}: feasible, objective {objective:g}")
            path = temporary / f"{name}.mps"
            generate("--haverly", str(number), "--out", str(path))
            # Line endings aside: a Windows checkout with core.autocrlf turns the committed LF
            # fixtures into CRLF, while the generator writes LF.
            require(path.read_bytes().replace(b"\r\n", b"\n")
                    == (fixtures / path.name).read_bytes().replace(b"\r\n", b"\n"),
                    f"{name} regeneration differs")
            print(f"{name} regeneration: identical up to line endings")

        # A dropped bilinear term must fail the independent primal check.
        model = read_mps(fixtures / "haverly1.mps")
        model["quadratic"]["pool_quality"] = {}
        point = dict(zip(points["columns"], points["instances"]["haverly1"]["global_solution"]))
        try:
            evaluate(model, point, -400)
        except AssertionError as error:
            require("pool_quality" in str(error), f"unexpected mutation failure: {error}")
        else:
            raise AssertionError("dropped bilinear terms escaped detection")
        print("planted missing-bilinear-term check: rejected")

        for components, pools, products, seed in ((4, 2, 3, 26119), (1, 1, 1, 0),
                                                  (6, 3, 4, 42)):
            path = temporary / "blend.mps"
            args = ("--components", str(components), "--pools", str(pools),
                    "--products", str(products), "--seed", str(seed), "--out", str(path))
            generate(*args)
            model = read_mps(path)
            start = json.loads(path.with_suffix(".start.json").read_text(encoding="utf-8"))
            evaluate(model, start)
            for k in range(pools):
                require(sum(start[f"y_{k}_{j}"] for j in range(products)) > 0,
                        "generated start has an empty pool")
            print(f"generated {components}c/{pools} pools/{products} products seed {seed}: start feasible")
            original = (path.read_bytes(), path.with_suffix(".start.json").read_bytes())
            generate(*args)
            require(original == (path.read_bytes(), path.with_suffix(".start.json").read_bytes()),
                    "seeded generation is not deterministic")
            print(f"generated seed {seed}: deterministic MPS and start")
    print("pooling check passed")


if __name__ == "__main__":
    main()
