#!/usr/bin/env python3
"""Step 0 of docs/IMPROVEMENT_PLAN.md: which component of HiGHS and SCIP solves each MIPLIB
instance that samaya loses?

Each baseline runs on the instances it solves (IMPROVEMENT_PLAN.md section 1), once with its
defaults and once with each component switched off. An instance that stops solving when a
component is off points at the technique samaya is missing there. The defaults are re-run on the
same machine, under the same load, so that "stops solving" compares like with like.

    bench/attribute.py [--time-limit 600] [--shards 4] [--out bench/results/attribution]
    bench/attribute.py --configs highs:presolve_rule_* --only-lost-by highs:presolve_off
    bench/attribute.py --configs scip:separator_* --instances n5-3,mc11 --list

Fine-grained configs are opt-in. --instance-dir selects the MPS directory; the legacy
--instances DIRECTORY form is also accepted. Quote globs if your shell expands them.

Every run is single-threaded with the 1e-4 gap, like bench/compare.sh. Results are appended to
OUT/runs.csv as they finish, so an interrupted run resumes where it stopped; OUT/summary.md is the
instance x component table. Run it alone on the machine (no build or other benchmark meanwhile).
"""

from __future__ import annotations

import argparse
import csv
import fnmatch
import multiprocessing
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
GAP = 1e-4

# The 18 instances HiGHS or SCIP solve at 600 s and samaya doesn't (IMPROVEMENT_PLAN.md section 1),
# with which baseline solved each there.
LOST = {
    "neos-3381206-awhea": ("highs", "scip"),
    "neos-911970": ("highs", "scip"),
    "n5-3": ("highs", "scip"),
    "mc11": ("highs", "scip"),
    "beasleyC3": ("highs", "scip"),
    "timtab1": ("highs", "scip"),
    "neos-3627168-kasai": ("highs", "scip"),
    "rococoC10-001000": ("highs",),
    "csched008": ("highs",),
    "gmu-35-40": ("highs",),
    "neos-4738912-atrato": ("highs",),
    "pg5_34": ("highs",),
    "tr12-30": ("highs",),
    "gen-ip054": ("highs",),
    "fastxgemm-n2r6s0t2": ("highs",),
    "ran14x18-disj-8": ("highs",),
    "graph20-20-1rand": ("scip",),
    "mcsched": ("scip",),
}

# Component switches, as named in IMPROVEMENT_PLAN.md section 2.
CONFIGS = {
    "highs": ["default", "presolve_off", "symmetry_off", "heuristics_off"],
    "scip": ["default", "presolve_off", "separating_off", "heuristics_off", "symmetry_off",
             "conflict_off"],
}

# HiGHS 1.15.1 highs/lp_data/HConst.h: PresolveRuleType. Bits 0--5 cannot
# be disabled; HPresolveAnalysis::setup uses 1 << rule_type, not a shifted index.
# Confirmed against the installed solver's "Presolve rules not allowed" log.
HIGHS_PRESOLVE_RULE_BITS = {
    "forcing_row": 6,
    "forcing_col": 7,
    "free_col_substitution": 8,
    "doubleton_equation": 9,
    "dependent_equations": 10,
    "dependent_free_cols": 11,
    "aggregator": 12,
    "parallel_rows_and_cols": 13,
    "sparsify": 14,
    "probing": 15,
    "enumeration": 16,
    "dual_fixing": 17,
    "col_stuffing": 18,
    "initial_sweep": 19,
}
HIGHS_RULE_CONFIGS = {f"presolve_rule_{name}_off": 1 << bit
                      for name, bit in HIGHS_PRESOLVE_RULE_BITS.items()}


def separator_configs() -> dict[str, str]:
    """Discover only separator plugins present in this SCIP installation."""
    import pyscipopt

    model = pyscipopt.Model()
    return {f"separator_{param.split('/')[1]}_off": param
            for param in sorted(model.getParams())
            if param.startswith("separating/") and param.endswith("/freq")
            and len(param.split('/')) == 3}


def select_configs(patterns: list[str] | None) -> dict[str, list[str]]:
    if patterns is None:
        return {solver: list(configs) for solver, configs in CONFIGS.items()}
    available = {"highs": CONFIGS["highs"] + list(HIGHS_RULE_CONFIGS),
                 "scip": list(CONFIGS["scip"])}
    # Default jobs should not need either solver installed merely to list them.
    if any(p.startswith("scip:") for p in patterns):
        available["scip"] += list(separator_configs())
    selected = {solver: [] for solver in CONFIGS}
    for pattern in patterns:
        solver, sep, glob = pattern.partition(":")
        if not sep or solver not in available or not glob:
            raise ValueError(f"expected solver:glob, got {pattern!r}")
        matches = [c for c in available[solver] if fnmatch.fnmatchcase(c, glob)]
        if not matches:
            raise ValueError(f"no configs match {pattern!r}")
        selected[solver].extend(c for c in matches if c not in selected[solver])
    return selected


def run_highs(instance: Path, config: str, time_limit: float) -> dict:
    import highspy

    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.setOptionValue("time_limit", float(time_limit))
    h.setOptionValue("threads", 1)
    h.setOptionValue("mip_rel_gap", GAP)
    if config == "presolve_off":
        h.setOptionValue("presolve", "off")
    elif config == "symmetry_off":
        h.setOptionValue("mip_detect_symmetry", False)
    elif config == "heuristics_off":
        h.setOptionValue("mip_heuristic_effort", 0.0)
    elif config in HIGHS_RULE_CONFIGS:
        if h.setOptionValue("presolve_rule_off", HIGHS_RULE_CONFIGS[config]) != highspy.HighsStatus.kOk:
            raise ValueError(f"HiGHS rejected presolve_rule_off for {config}")
    h.readModel(str(instance))
    start = time.perf_counter()
    h.run()
    seconds = time.perf_counter() - start
    raw = h.modelStatusToString(h.getModelStatus()).lower()
    status = {"optimal": "optimal", "infeasible": "infeasible",
              "time limit reached": "time_limit"}.get(raw, raw.replace(" ", "_"))
    return {"status": status, "seconds": seconds, "nodes": h.getInfo().mip_node_count}


def run_scip(instance: Path, config: str, time_limit: float) -> dict:
    import pyscipopt

    m = pyscipopt.Model()
    m.hideOutput()
    m.readProblem(str(instance))
    m.setParam("limits/time", float(time_limit))
    m.setParam("limits/gap", GAP)
    off = pyscipopt.SCIP_PARAMSETTING.OFF
    if config == "presolve_off":
        m.setPresolve(off)
    elif config == "separating_off":
        m.setSeparating(off)
    elif config == "heuristics_off":
        m.setHeuristics(off)
    elif config == "symmetry_off":
        m.setParam("misc/usesymmetry", 0)
    elif config == "conflict_off":
        m.setParam("conflict/enable", False)
    elif config.startswith("separator_") and config.endswith("_off"):
        # SCIP's freq=-1 means never, including the root; maxroundsroot is
        # only a cap on calls and does not override this scheduling switch.
        name = config[len("separator_"):-len("_off")]
        m.setParam(f"separating/{name}/freq", -1)
    start = time.perf_counter()
    m.optimize()
    seconds = time.perf_counter() - start
    raw = m.getStatus()
    status = {"optimal": "optimal", "infeasible": "infeasible", "timelimit": "time_limit",
              "gaplimit": "optimal"}.get(raw, raw)
    return {"status": status, "seconds": seconds, "nodes": m.getNNodes()}


def run_one(job: tuple[str, str, str, str, float]) -> dict:
    name, path, solver, config, time_limit = job
    runner = run_highs if solver == "highs" else run_scip
    try:
        result = runner(Path(path), config, time_limit)
    except Exception as exc:  # A baseline failing must not stop the other runs.
        result = {"status": f"error: {exc}", "seconds": float("nan"), "nodes": ""}
    return {"instance": name, "solver": solver, "config": config, **result}


def solved(row: dict) -> bool:
    return row["status"] in ("optimal", "infeasible")


def write_summary(rows: list[dict], out: Path, time_limit: float) -> None:
    by = {(r["instance"], r["solver"], r["config"]): r for r in rows}
    lines = ["# Step 0: attribution runs", "",
             f"Each cell: solve time in seconds, or TL (not solved in {time_limit:g} s). A component "
             "whose switch-off turns a solved default into TL is what samaya is likely missing.", ""]
    for solver, configs in CONFIGS.items():
        names = [n for n, s in LOST.items() if solver in s]
        lines += [f"## {solver}", "", "| Instance | " + " | ".join(configs) + " |",
                  "|---|" + "---|" * len(configs)]
        for name in names:
            cells = []
            for config in configs:
                r = by.get((name, solver, config))
                if r is None:
                    cells.append("–")
                elif solved(r):
                    cells.append(f"{float(r['seconds']):.0f}")
                else:
                    cells.append("TL" if r["status"] == "time_limit" else r["status"])
            lines.append(f"| {name} | " + " | ".join(cells) + " |")
        lines.append("")
    # Keep the original component tables intact. Add only families represented
    # in the CSV, so a legacy/default run produces the same summary as before.
    for solver, prefix, title in (("highs", "presolve_rule_", "presolve rules"),
                                   ("scip", "separator_", "separators")):
        configs = sorted({r["config"] for r in rows
                          if r["solver"] == solver and r["config"].startswith(prefix)})
        if not configs:
            continue
        names = [n for n, s in LOST.items() if solver in s]
        lines += [f"## {solver} {title}", "", "| Instance | " + " | ".join(configs) + " |",
                  "|---|" + "---|" * len(configs)]
        for name in names:
            cells = []
            for config in configs:
                r = by.get((name, solver, config))
                if r is None:
                    cells.append("–")
                elif solved(r):
                    cells.append(f"{float(r['seconds']):.0f}")
                else:
                    cells.append("TL" if r["status"] == "time_limit" else r["status"])
            lines.append(f"| {name} | " + " | ".join(cells) + " |")
        lines.append("")
    (out / "summary.md").write_text("\n".join(lines))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--time-limit", type=float, default=600.0)
    parser.add_argument("--shards", type=int, default=4)
    parser.add_argument("--instances", help="comma-separated LOST names (or a legacy MPS directory)")
    parser.add_argument("--instance-dir", type=Path,
                        default=REPO / "bench" / "instances" / "miplib_small")
    parser.add_argument("--configs", action="append", nargs="+", metavar="SOLVER:GLOB",
                        help="select configs; repeat or supply multiple solver:glob patterns")
    parser.add_argument("--only-lost-by", metavar="SOLVER:CONFIG",
                        help="keep instances solved by default but unsolved with this config in OUT/runs.csv")
    parser.add_argument("--list", action="store_true", help="list pending jobs without solving or requiring MPS files")
    parser.add_argument("--out", type=Path, default=REPO / "bench" / "results" / "attribution")
    args = parser.parse_args()
    try:
        configs = select_configs([p for group in args.configs for p in group]
                                 if args.configs is not None else None)
    except (ValueError, ImportError) as exc:
        parser.error(str(exc))
    names = set(LOST)
    if args.instances is not None:
        # Preserve the old --instances directory form, including missing paths
        # containing a directory separator; bare names are the new filter.
        if Path(args.instances).is_dir() or "/" in args.instances or "\\" in args.instances:
            args.instance_dir = Path(args.instances)
        else:
            names = set(args.instances.split(","))
            unknown = names - LOST.keys()
            if unknown:
                parser.error("unknown instances: " + ", ".join(sorted(unknown)))
    runs_csv = args.out / "runs.csv"
    fields = ["instance", "solver", "config", "status", "seconds", "nodes"]
    if not args.list:
        args.out.mkdir(parents=True, exist_ok=True)

    done: list[dict] = []
    if runs_csv.exists():
        with runs_csv.open() as f:
            done = list(csv.DictReader(f))
    finished = {(r["instance"], r["solver"], r["config"]) for r in done}

    if args.only_lost_by:
        solver, sep, config = args.only_lost_by.partition(":")
        if not sep or solver not in CONFIGS or not config:
            parser.error("--only-lost-by requires solver:config")
        if not runs_csv.exists():
            parser.error(f"--only-lost-by requires previous results in {runs_csv}")
        by = {(r["instance"], r["solver"], r["config"]): r for r in done}
        if not any(r["solver"] == solver and r["config"] == config for r in done):
            parser.error(f"no previous results for {args.only_lost_by}")
        names &= {name for name in LOST
                  if (name, solver, "default") in by
                  and solved(by[name, solver, "default"])
                  and (name, solver, config) in by
                  and not solved(by[name, solver, config])
                  and not by[name, solver, config]["status"].startswith("error:")}

    jobs = []
    for name, solvers in LOST.items():
        if name not in names:
            continue
        path = args.instance_dir / f"{name}.mps"
        if not args.list and any(configs[solver] for solver in solvers) and not path.exists():
            print(f"missing instance {path}", file=sys.stderr)
            return 1
        for solver in solvers:
            for config in configs[solver]:
                if (name, solver, config) not in finished:
                    jobs.append((name, str(path), solver, config, args.time_limit))
    print(f"{len(jobs)} runs to do ({len(finished)} already done), {args.shards} at a time",
          flush=True)
    if args.list:
        for name, _, solver, config, _ in jobs:
            print(f"{solver} {name} {config}")
        return 0

    new_file = not runs_csv.exists()
    with runs_csv.open("a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        if new_file:
            writer.writeheader()
        # One process per run (maxtasksperchild=1), so a solver's state never leaks into the next.
        with multiprocessing.Pool(args.shards, maxtasksperchild=1) as pool:
            for k, row in enumerate(pool.imap_unordered(run_one, jobs), 1):
                writer.writerow({key: row.get(key, "") for key in fields})
                f.flush()
                done.append(row)
                print(f"[{k}/{len(jobs)}] {row['instance']} {row['solver']} {row['config']}: "
                      f"{row['status']} {float(row['seconds']):.1f} s", flush=True)
    write_summary(done, args.out, args.time_limit)
    print(f"summary: {args.out / 'summary.md'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
