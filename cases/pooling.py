#!/usr/bin/env python3
"""Refinery-style pooling/blending models in the p-formulation (standard library only).

Component i has cost c_i and sulfur q_i. Flow x_i_k enters pool k, y_k_j
leaves it for product j, and z_i_j goes directly to that product. Pool quality
p_k obeys sum_i q_i*x_i_k = p_k*sum_j y_k_j, alongside the pool mass balance.
Product sulfur is limited by sum_k p_k*y_k_j + sum_i q_i*z_i_j <=
spec_j*(sum_k y_k_j + sum_i z_i_j). Minimize feed costs minus product revenue,
subject to component availability, pool capacities and product demand caps.
Flows are nonnegative; pool qualities lie between component quality extremes.

QCMATRIX contains the full symmetric matrix for x'Qx, without a factor 1/2.
Generated data are built around positive pool flows; the known feasible point
is saved as <name>.start.json beside the MPS. Figures are synthetic, not MRPL data.

Usage:
  python cases/pooling.py --components 4 --pools 2 --products 3 --seed 26119 --out blend.mps
  python cases/pooling.py --haverly 1 --out tests/instances/haverly1.mps
"""

from __future__ import annotations

import argparse
import json
import math
import random
from pathlib import Path


class PoolingWriter:
    """Small free-MPS writer for continuous models with quadratic constraint rows."""

    def __init__(self, name: str) -> None:
        self.name = name
        self.cols: dict[str, tuple[float, float, float | None]] = {}
        self.rows: dict[str, tuple[str, float, dict[str, float]]] = {}
        self.quadratic: dict[str, list[tuple[str, str, float]]] = {}

    def col(self, name: str, cost: float = 0.0, lo: float = 0.0,
            up: float | None = None) -> None:
        self.cols[name] = (cost, lo, up)

    def row(self, name: str, kind: str, rhs: float, terms: dict[str, float]) -> None:
        self.rows[name] = (kind, rhs, terms)

    def bilinear(self, row: str, first: str, second: str, value: float) -> None:
        # Store polynomial coefficients; write both matrix entries at half strength.
        self.quadratic.setdefault(row, []).append((first, second, value))

    def write(self, path: Path) -> None:
        lines = [f"NAME {self.name}", "ROWS", " N obj"]
        lines.extend(f" {kind} {name}" for name, (kind, _, _) in self.rows.items())
        lines.append("COLUMNS")
        for col, (cost, _, _) in self.cols.items():
            # Explicit zero costs declare columns used only in quadratic rows.
            lines.append(f" {col} obj {cost:.17g}")
            for row, (_, _, terms) in self.rows.items():
                if terms.get(col, 0.0):
                    lines.append(f" {col} {row} {terms[col]:.17g}")
        lines.append("RHS")
        for row, (_, rhs, _) in self.rows.items():
            if rhs:
                lines.append(f" rhs {row} {rhs:.17g}")
        lines.append("BOUNDS")
        for col, (_, lo, up) in self.cols.items():
            if lo != 0.0:
                lines.append(f" LO bnd {col} {lo:.17g}")
            if up is not None:
                lines.append(f" UP bnd {col} {up:.17g}")
        for row, terms in self.quadratic.items():
            lines.append(f"QCMATRIX {row}")
            for first, second, value in terms:
                if first == second:
                    lines.append(f" {first} {second} {value:.17g}")
                else:
                    lines.append(f" {first} {second} {value / 2:.17g}")
                    lines.append(f" {second} {first} {value / 2:.17g}")
        lines.append("ENDATA")
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(("\n".join(lines) + "\n").encode("utf-8"))


def haverly(number: int) -> PoolingWriter:
    demand, cost_b = {1: (100, 16), 2: (600, 16), 3: (100, 13)}[number]
    model = PoolingWriter(f"HAVERLY{number}")
    for name, cost in zip(("a", "b", "u", "v", "z1", "z2"),
                          (6, cost_b, -9, -15, 1, -5)):
        model.col(name, cost)
    model.col("p", lo=1, up=3)
    model.row("pool_balance", "E", 0, {"a": 1, "b": 1, "u": -1, "v": -1})
    model.row("pool_quality", "E", 0, {"a": 3, "b": 1})
    model.row("prod1_quality", "L", 0, {"u": -2.5, "z1": -0.5})
    model.row("prod2_quality", "L", 0, {"v": -1.5, "z2": 0.5})
    model.row("prod1_demand", "L", demand, {"u": 1, "z1": 1})
    model.row("prod2_demand", "L", 200, {"v": 1, "z2": 1})
    model.bilinear("pool_quality", "p", "u", -1)
    model.bilinear("pool_quality", "p", "v", -1)
    model.bilinear("prod1_quality", "p", "u", 1)
    model.bilinear("prod2_quality", "p", "v", 1)
    return model


def pooling(components: int, pools: int, products: int,
            seed: int) -> tuple[PoolingWriter, dict[str, float]]:
    if min(components, pools, products) < 1:
        raise ValueError("components, pools and products must all be positive")
    rng = random.Random(seed)
    quality = [rng.uniform(0.1, 3.5) for _ in range(components)]
    cost = [rng.uniform(40, 70) for _ in range(components)]
    price = [rng.uniform(75, 110) for _ in range(products)]
    model = PoolingWriter(f"POOLING_{components}c_{pools}k_{products}p_{seed}")
    start: dict[str, float] = {}
    for i in range(components):
        for k in range(pools):
            name = f"x_{i}_{k}"
            model.col(name, cost[i])
            start[name] = float(rng.randint(5, 30))
    for k in range(pools):
        total = math.fsum(start[f"x_{i}_{k}"] for i in range(components))
        weights = [rng.randint(1, 10) for _ in range(products)]
        for j in range(products):
            name = f"y_{k}_{j}"
            model.col(name, -price[j])
            start[name] = total * weights[j] / sum(weights)
    for i in range(components):
        for j in range(products):
            name = f"z_{i}_{j}"
            model.col(name, cost[i] - price[j])
            start[name] = float(rng.randint(1, 10))
    for k in range(pools):
        name = f"p_{k}"
        model.col(name, lo=min(quality), up=max(quality))
        inflows = {f"x_{i}_{k}": 1 for i in range(components)}
        outflows = {f"y_{k}_{j}": 1 for j in range(products)}
        total = math.fsum(start[col] for col in inflows)
        start[name] = math.fsum(quality[i] * start[f"x_{i}_{k}"]
                                for i in range(components)) / total
        model.row(f"pool_balance_{k}", "E", 0,
                  {**inflows, **{col: -1 for col in outflows}})
        model.row(f"pool_capacity_{k}", "L", total * 1.25, inflows)
        row = f"pool_quality_{k}"
        model.row(row, "E", 0, {f"x_{i}_{k}": quality[i] for i in range(components)})
        for col in outflows:
            model.bilinear(row, name, col, -1)
    for i in range(components):
        terms = {f"x_{i}_{k}": 1 for k in range(pools)}
        terms.update({f"z_{i}_{j}": 1 for j in range(products)})
        model.row(f"component_supply_{i}", "L",
                  1.25 * math.fsum(start[col] for col in terms), terms)
    for j in range(products):
        flows = {f"y_{k}_{j}": 1 for k in range(pools)}
        flows.update({f"z_{i}_{j}": 1 for i in range(components)})
        total = math.fsum(start[col] for col in flows)
        sulfur = math.fsum([start[f"p_{k}"] * start[f"y_{k}_{j}"] for k in range(pools)]
                           + [quality[i] * start[f"z_{i}_{j}"] for i in range(components)])
        # A small quality margin and 25% spare capacity keep the planted start feasible.
        spec = sulfur / total + 0.05
        model.row(f"product_demand_{j}", "L", total * 1.25, flows)
        row = f"product_quality_{j}"
        model.row(row, "L", 0, {**{f"y_{k}_{j}": -spec for k in range(pools)},
                                **{f"z_{i}_{j}": quality[i] - spec
                                   for i in range(components)}})
        for k in range(pools):
            model.bilinear(row, f"p_{k}", f"y_{k}_{j}", 1)
    return model, start


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--components", type=int, default=4)
    parser.add_argument("--pools", type=int, default=2)
    parser.add_argument("--products", type=int, default=3)
    parser.add_argument("--seed", type=int, default=26119)
    parser.add_argument("--haverly", type=int, choices=(1, 2, 3))
    parser.add_argument("--out", type=Path, help="output MPS path")
    args = parser.parse_args()
    if args.haverly is not None:
        path = args.out or Path(f"haverly{args.haverly}.mps")
        haverly(args.haverly).write(path)
    else:
        if min(args.components, args.pools, args.products) < 1:
            parser.error("components, pools and products must all be positive")
        path = args.out or Path("pooling.mps")
        model, start = pooling(args.components, args.pools, args.products, args.seed)
        model.write(path)
        path.with_suffix(".start.json").write_text(json.dumps(start, indent=2) + "\n",
                                                  encoding="utf-8")
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
