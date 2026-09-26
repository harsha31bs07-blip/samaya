#!/usr/bin/env python3
"""Runs samaya (and optionally baseline solvers) over a set of instances and writes a CSV.

Examples:
    bench/harness.py bench/instances/netlib --time-limit 60
    bench/harness.py bench/instances/netlib --baseline highs --out bench/results/netlib.csv

Baselines are external executables used only for comparison; they are never linked into samaya.
The summary reports solved counts and the shifted geometric mean of solve time (shift 10 s), the
standard metric in the Mittelmann benchmarks.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DEFAULT_SAMAYA = REPO / "build" / "release" / "apps" / "cli" / "samaya"
SOLVED = {"optimal", "infeasible", "unbounded"}


def run_samaya(exe: Path, instance: Path, time_limit: float, threads: int, gap: float) -> dict:
    cmd = [str(exe), "--json", "--log-level", "0", "--time-limit", str(time_limit),
           "--threads", str(threads), "--mip-gap", str(gap), str(instance)]
    start = time.perf_counter()
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=time_limit + 60)
    except subprocess.TimeoutExpired:
        return {"status": "killed", "wall_seconds": time.perf_counter() - start}
    wall = time.perf_counter() - start
    lines = [line for line in proc.stdout.splitlines() if line.startswith("{")]
    if not lines:
        return {"status": "crash", "wall_seconds": wall, "message": proc.stderr.strip()[-200:]}
    result = json.loads(lines[-1])
    result["wall_seconds"] = wall
    return result


def run_highs(exe: str, instance: Path, time_limit: float, threads: int, gap: float) -> dict:
    cmd = [exe, "--time_limit", str(time_limit), "--model_file", str(instance)]
    start = time.perf_counter()
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=time_limit + 60)
    except subprocess.TimeoutExpired:
        return {"status": "killed", "wall_seconds": time.perf_counter() - start}
    wall = time.perf_counter() - start
    out = proc.stdout
    status_match = re.search(r"Model status\s*:\s*(.+)", out)
    obj_match = re.search(r"Objective value\s*:\s*(\S+)", out)
    raw = status_match.group(1).strip().lower() if status_match else "unknown"
    status = {"optimal": "optimal", "infeasible": "infeasible", "unbounded": "unbounded",
              "time limit reached": "time_limit"}.get(raw, raw.replace(" ", "_"))
    return {"status": status, "objective": float(obj_match.group(1)) if obj_match else None,
            "solve_seconds": wall, "wall_seconds": wall}


def run_highspy(instance: Path, time_limit: float, threads: int, gap: float) -> dict:
    import highspy  # Optional dependency, used only as a comparison baseline.

    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.setOptionValue("time_limit", float(time_limit))
    h.setOptionValue("threads", int(threads))
    h.setOptionValue("mip_rel_gap", float(gap))
    h.readModel(str(instance))
    start = time.perf_counter()
    h.run()
    wall = time.perf_counter() - start
    raw = h.modelStatusToString(h.getModelStatus()).lower()
    status = {"optimal": "optimal", "infeasible": "infeasible", "unbounded": "unbounded",
              "time limit reached": "time_limit"}.get(raw, raw.replace(" ", "_"))
    info = h.getInfo()
    # The best solution's objective, also when a limit stopped the search.
    has_solution = info.primal_solution_status == 2  # kSolutionStatusFeasible
    objective = info.objective_function_value if status == "optimal" or has_solution else None
    is_mip = any(t != highspy.HighsVarType.kContinuous for t in h.getLp().integrality_)
    return {"status": status, "objective": objective, "solve_seconds": wall, "wall_seconds": wall,
            "nodes": info.mip_node_count if is_mip else None,
            "class": "MILP" if is_mip else "LP"}


def run_scip(instance: Path, time_limit: float, threads: int, gap: float) -> dict:
    import pyscipopt  # Optional dependency (pip install pyscipopt), a comparison baseline only.

    m = pyscipopt.Model()
    m.hideOutput()
    m.readProblem(str(instance))
    m.setParam("limits/time", float(time_limit))
    m.setParam("limits/gap", float(gap))
    start = time.perf_counter()
    m.optimize()
    wall = time.perf_counter() - start
    raw = m.getStatus()
    status = {"optimal": "optimal", "infeasible": "infeasible", "unbounded": "unbounded",
              "inforunbd": "infeasible_or_unbounded", "timelimit": "time_limit",
              "gaplimit": "optimal"}.get(raw, raw)  # Stopped at the requested gap.
    objective = m.getObjVal() if m.getNSols() > 0 else None
    is_mip = any(v.vtype() != "CONTINUOUS" for v in m.getVars(transformed=False))
    return {"status": status, "objective": objective, "solve_seconds": wall, "wall_seconds": wall,
            "nodes": m.getNNodes() if is_mip else None, "class": "MILP" if is_mip else "LP"}


def to_lp_file(instance: Path, tmp: Path) -> Path:
    """CBC reads free MPS as fixed format and GLPK rejects OBJSENSE, so both get the model in
    CPLEX LP format, written by SCIP (pip install pyscipopt). Conversion time is not measured."""
    import pyscipopt

    m = pyscipopt.Model()
    m.hideOutput()
    m.readProblem(str(instance))
    # An objective constant (e226) is written as a bare number that GLPK rejects and CBC drops:
    # carry it on a column fixed at 1 instead.
    offset = m.getObjoffset(original=True)
    if offset != 0.0:
        m.addVar(name="objective_constant", lb=1.0, ub=1.0, obj=offset)
        m.addObjoffset(-offset)
    out = tmp / "model.lp"
    m.writeProblem(str(out), trans=False, genericnames=True, verbose=False)
    # GLPK rejects rows without terms ("c1: = +0"); drop those that hold trivially (0 = 0,
    # 0 <= b with b >= 0, 0 >= b with b <= 0) and keep any that could make the model infeasible.
    kept = []
    for line in out.read_text().splitlines():
        empty = re.match(r"^\s*\S+:\s*(<=|>=|=)\s*([+-]?[0-9.eE+-]+)\s*$", line)
        if empty:
            op, rhs = empty.group(1), float(empty.group(2))
            if (op == "=" and rhs == 0) or (op == "<=" and rhs >= 0) or (op == ">=" and rhs <= 0):
                continue
        kept.append(line)
    out.write_text("\n".join(kept) + "\n")
    return out

def run_cbc(exe: str, instance: Path, time_limit: float, threads: int, gap: float) -> dict:
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        model = to_lp_file(instance, Path(tmp))
        cmd = [exe, str(model), "-seconds", str(time_limit), "-ratioGap", str(gap), "-threads", str(threads),
            "-solve", "-quit"]
        start = time.perf_counter()
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=time_limit + 60)
        except subprocess.TimeoutExpired:
            return {"status": "killed", "wall_seconds": time.perf_counter() - start}
        wall = time.perf_counter() - start
    out = proc.stdout
    status = "unknown"
    if re.search(r"Result - Optimal solution found|^Optimal - objective value", out, re.M):
        status = "optimal"
    elif re.search(r"infeasible", out, re.I) and re.search(r"Result - |Problem is infeasible",
                                                            out):
        status = "infeasible"
    elif re.search(r"unbounded", out, re.I) and re.search(r"Result - |Problem is unbounded", out):
        status = "unbounded"
    elif re.search(r"Stopped on time", out):
        status = "time_limit"
    obj = re.search(r"Objective value:\s*(\S+)", out) or re.search(
        r"Optimal - objective value\s*(\S+)", out)
    objective = None
    if obj and status in ("optimal", "time_limit"):
        try:
            objective = float(obj.group(1))
        except ValueError:
            objective = None
    if status == "time_limit" and re.search(r"No feasible solution found", out):
        objective = None
    nodes = re.search(r"Enumerated nodes:\s*(\d+)", out)
    return {"status": status, "objective": objective, "solve_seconds": wall, "wall_seconds": wall,
            "nodes": int(nodes.group(1)) if nodes else None}


def run_glpk(exe: str, instance: Path, time_limit: float, threads: int, gap: float) -> dict:
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        model = to_lp_file(instance, Path(tmp))
        report = Path(tmp) / "out.txt"
        cmd = [exe, "--lp", str(model), "--tmlim", str(max(1, int(time_limit))),
               "--mipgap", str(gap), "-o", str(report)]
        start = time.perf_counter()
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=time_limit + 60)
        except subprocess.TimeoutExpired:
            return {"status": "killed", "wall_seconds": time.perf_counter() - start}
        wall = time.perf_counter() - start
        text = report.read_text() if report.exists() else ""
    out = proc.stdout
    raw = re.search(r"^Status:\s*(.+)$", text, re.M)
    raw = raw.group(1).strip() if raw else ""
    status = {"OPTIMAL": "optimal", "INTEGER OPTIMAL": "optimal", "INFEASIBLE (FINAL)": "infeasible",
              "INTEGER EMPTY": "infeasible", "UNBOUNDED": "unbounded",
              "INTEGER NON-OPTIMAL": "time_limit", "UNDEFINED": "time_limit",
              "INTEGER UNDEFINED": "time_limit"}.get(raw, "unknown")
    if "TIME LIMIT EXCEEDED" in out:
        status = "time_limit"
    elif "RELATIVE MIP GAP TOLERANCE REACHED" in out and status == "time_limit":
        status = "optimal"  # Stopped at the requested gap, like every other solver.
    if re.search(r"PROBLEM HAS NO (PRIMAL )?FEASIBLE", out):
        status = "infeasible"
    obj = re.search(r"^Objective:\s*\S+\s*=\s*(\S+)", text, re.M)
    objective = float(obj.group(1)) if obj and status in ("optimal", "time_limit") and \
        raw not in ("UNDEFINED", "INTEGER UNDEFINED") else None
    return {"status": status, "objective": objective, "solve_seconds": wall, "wall_seconds": wall}


def read_solu(path: Path) -> dict[str, tuple[str, float | None]]:
    """Reads a MIPLIB .solu file: =opt= / =best= name value, =inf= name."""
    known: dict[str, tuple[str, float | None]] = {}
    for line in path.read_text().splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[0] in ("=opt=", "=best=", "=inf=", "=unbd="):
            value = float(parts[2]) if len(parts) > 2 and parts[0] in ("=opt=", "=best=") else None
            known[parts[1]] = (parts[0].strip("="), value)
    return known


def shifted_geomean(values: list[float], shift: float = 10.0) -> float:
    if not values:
        return float("nan")
    return math.exp(sum(math.log(v + shift) for v in values) / len(values)) - shift


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("instances", nargs="+", type=Path,
                        help="instance files or directories containing *.mps / *.qps")
    parser.add_argument("--samaya", type=Path, default=DEFAULT_SAMAYA)
    parser.add_argument("--baseline", action="append", default=[],
                        choices=["highs", "highspy", "scip", "cbc", "glpk"],
                        help="also run a baseline solver: the highs, cbc or glpsol executable on "
                             "PATH, or the highspy / pyscipopt Python module (scip)")
    parser.add_argument("--time-limit", type=float, default=300.0)
    parser.add_argument("--samaya-clock-factor", type=float, default=1.0,
                        help="multiply samaya's time limit by this: samaya times its limit on the "
                             "monotonic clock, the baselines on the real-time clock, and on a "
                             "machine whose monotonic clock runs fast (WSL2 on some laptops) "
                             "samaya would stop early; bench/clock_factor.py measures it")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--mip-gap", type=float, default=1e-4,
                        help="relative MIP gap for every solver (default 1e-4)")
    parser.add_argument("--solu", type=Path,
                        help="MIPLIB .solu file of known optimal values to check samaya against")
    parser.add_argument("--out", type=Path, default=REPO / "bench" / "results" / "results.csv")
    args = parser.parse_args()

    files: list[Path] = []
    for path in args.instances:
        if path.is_dir():
            files += sorted(p for p in path.iterdir() if p.suffix.lower() in {".mps", ".qps"})
        else:
            files.append(path)
    if not files:
        print("no instances found", file=sys.stderr)
        return 1
    if not args.samaya.exists():
        print(f"samaya binary not found at {args.samaya}; build with "
              "`cmake --preset release && cmake --build --preset release`", file=sys.stderr)
        return 1

    gap = args.mip_gap
    samaya_limit = args.time_limit * args.samaya_clock_factor
    if args.samaya_clock_factor != 1.0:
        print(f"samaya time limit {samaya_limit:.2f} s on its own clock (factor "
              f"{args.samaya_clock_factor}), {args.time_limit} s of real time", flush=True)
    solvers = {"samaya": lambda f: run_samaya(args.samaya, f, samaya_limit, args.threads, gap)}
    for name in args.baseline:
        if name == "highspy":
            solvers[name] = lambda f: run_highspy(f, args.time_limit, args.threads, gap)
            continue
        if name == "scip":
            solvers[name] = lambda f: run_scip(f, args.time_limit, args.threads, gap)
            continue
        exe = shutil.which({"glpk": "glpsol"}.get(name, name))
        if exe is None:
            print(f"baseline '{name}' not found on PATH", file=sys.stderr)
            return 1
        runner = {"highs": run_highs, "cbc": run_cbc, "glpk": run_glpk}[name]
        solvers[name] = lambda f, exe=exe, runner=runner: runner(
            exe, f, args.time_limit, args.threads, gap)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for instance in files:
        for solver, run in solvers.items():
            r = run(instance)
            name = instance.name
            for suffix in (".gz", ".mps", ".qps"):
                name = name.removesuffix(suffix)
            row = {"instance": name, "solver": solver, "status": r.get("status"),
                   "objective": r.get("objective"), "solve_seconds": r.get("solve_seconds"),
                   "wall_seconds": round(r["wall_seconds"], 4), "nodes": r.get("nodes"),
                   "class": r.get("class"), "verified": r.get("verified")}
            rows.append(row)
            print(f"{instance.stem:<20} {solver:<8} {row['status']:<16} obj={row['objective']}"
                  f" t={row['wall_seconds']:.3f}s", flush=True)

    with args.out.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nWrote {args.out}")
    if len(solvers) > 1:
        by_instance: dict[str, dict[str, dict]] = {}
        for r in rows:
            by_instance.setdefault(r["instance"], {})[r["solver"]] = r
        disagreements = 0
        for inst, results in by_instance.items():
            base = [r for s, r in results.items() if s != "samaya"]
            mine = results.get("samaya")
            for other in base:
                # Only final answers can contradict each other; a limit is not a disagreement.
                if mine["status"] not in SOLVED or other["status"] not in SOLVED:
                    continue
                same_status = mine["status"] == other["status"]
                same_obj = True
                if same_status and mine["status"] == "optimal":
                    a, b = mine["objective"], other["objective"]
                    # Each MIP solver may stop anywhere within its relative gap.
                    tol = 2 * gap if mine.get("class") == "MILP" else 1e-6
                    same_obj = abs(a - b) <= tol * max(1.0, abs(b))
                if not (same_status and same_obj):
                    disagreements += 1
                    print(f"  DISAGREE {inst}: samaya {mine['status']} {mine['objective']} vs "
                          f"{other['solver']} {other['status']} {other['objective']}")
        print(f"Status/objective disagreements with baselines: {disagreements}")
    if args.solu:
        known = read_solu(args.solu)
        correct = wrong = unknown = 0
        for r in rows:
            if r["solver"] != "samaya":
                continue
            kind, value = known.get(r["instance"], (None, None))
            if kind is None or r["status"] not in SOLVED:
                unknown += kind is None
                continue
            if kind == "inf":
                ok = r["status"] == "infeasible"
            elif kind == "opt" and r["status"] == "optimal":
                ok = abs(r["objective"] - value) <= 2 * gap * max(1.0, abs(value))
            else:
                continue
            correct += ok
            wrong += not ok
            if not ok:
                print(f"  WRONG {r['instance']}: samaya {r['status']} {r['objective']} vs "
                      f"known {kind} {value}")
        print(f"Known solutions: {correct} solved correctly, {wrong} wrong, "
              f"{unknown} not in {args.solu.name}")
    for solver in solvers:
        mine = [r for r in rows if r["solver"] == solver]
        solved = [r for r in mine if r["status"] in SOLVED]
        # Unsolved instances count at the time limit, as in the Mittelmann benchmarks.
        times = [r["wall_seconds"] if r["status"] in SOLVED else args.time_limit for r in mine]
        print(f"{solver:<8} solved {len(solved)}/{len(mine)}, "
              f"shifted geomean {shifted_geomean(times):.3f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
