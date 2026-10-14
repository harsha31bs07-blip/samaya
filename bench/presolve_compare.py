#!/usr/bin/env python3
"""Compare samaya and HiGHS presolve sizes on bench/miplib_small.test.

    bench/presolve_compare.py --instances pk1,neos5 --jobs 2
    bench/presolve_compare.py --self-test

Fetch MPS files with bench/fetch_instances.sh miplib-list bench/miplib_small.test.
HiGHS uses highspy (imported only when running a baseline), one thread and presolve().
samaya has no presolve-only CLI: a one-second solve limit bounds subsequent search;
its logged presolve_seconds excludes that search and model reading. Presolve itself
is not interrupted by the solve limit. Each samaya subprocess has a 61-second
wall timeout; partial logs are retained. Results are runs.csv and summary.md in OUT.
"""

from __future__ import annotations

import argparse
from contextlib import nullcontext
import csv
import json
import math
import multiprocessing
import re
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DEFAULT_SAMAYA = REPO / "build" / "release" / "apps" / "cli" / "samaya"
COUNTERS = ["passes", "rows_removed", "cols_removed", "empty_rows", "singleton_rows",
            "forcing_rows", "redundant_rows", "duplicate_rows", "fixed_cols", "empty_cols",
            "dominated_cols", "free_col_singletons", "bounds_tightened",
            "implied_free_col_singletons", "implied_free_rejected_by_tolerance"]
DOUBLETON_COUNTERS = ["doubleton_continuous_cols", "doubleton_integer_cols",
                      "doubleton_rejected_pivot", "doubleton_rejected_integer_data",
                      "doubleton_rejected_divisibility", "doubleton_rejected_integer_bound",
                      "doubleton_rejected_crossed_bounds", "doubleton_rejected_residual",
                      "doubleton_rejected_nonfinite"]
COUNTERS += DOUBLETON_COUNTERS
SUBSTITUTION_COUNTERS = ["substitution_continuous_cols", "substitution_integer_cols",
                         "substitution_rejected_row_length", "substitution_rejected_col_length",
                         "substitution_rejected_not_implied_free", "substitution_rejected_by_tolerance",
                         "substitution_rejected_pivot", "substitution_rejected_integer",
                         "substitution_rejected_nonfinite", "substitution_rejected_residual",
                         "substitution_rejected_fill"]
COUNTERS += SUBSTITUTION_COUNTERS
SIZES = ["rows", "cols", "integers", "nonzeros"]
FIELDS = (["instance"] + [f"original_{s}" for s in SIZES]
          + ["samaya_status", "samaya_solve_status", "samaya_seconds", "samaya_wall_seconds"]
          + [f"samaya_{s}" for s in SIZES] + [f"samaya_{c}" for c in COUNTERS]
          + ["samaya_error", "highs_version", "highs_status", "highs_seconds"]
          + [f"highs_{s}" for s in SIZES] + ["highs_error", "column_ratio"])


def parse_samaya(text: str) -> dict:
    """Use the first presolve result; a later LP fallback must not replace it."""
    result = {"samaya_status": "failed"}
    original = re.search(r"samaya [^\n]+ with (\d+) rows, (\d+) cols, (\d+) nonzeros", text)
    if original:
        result.update(zip(["original_rows", "original_cols", "original_nonzeros"],
                          map(int, original.groups())))
    integers = re.search(r"cols \d+ \((\d+) integer", text)
    if integers:
        result["original_integers"] = int(integers[1])
    line = next((line for line in text.splitlines()
                 if re.search(r"Presolve: (?:\d+ rows|infeasib)", line)), "")
    seconds = re.search(r"presolve_seconds ([\d.eE+-]+)", line)
    if seconds:
        result["samaya_seconds"] = float(seconds[1])
    counters = re.search(r"Presolve: rows: (\d+) empty, (\d+) singleton, (\d+) forcing, "
                         r"(\d+) redundant, (\d+) duplicate; cols: (\d+) fixed, (\d+) empty, "
                         r"(\d+) dominated, (\d+) free singleton; (\d+) bounds tightened; "
                         r"implied_free_col_singletons (\d+); "
                         r"implied_free_rejected_by_tolerance (\d+)", text)
    if counters:
        result.update(zip([f"samaya_{c}" for c in COUNTERS[3:]], map(int, counters.groups())))
    # Binaries before the doubleton rule do not print this line; their fields stay blank.
    doubletons = re.search(r"Presolve doubletons: (\d+) continuous, (\d+) integer; rejected: "
                           r"(\d+) pivot, (\d+) integer_data, (\d+) divisibility, "
                           r"(\d+) integer_bound, (\d+) crossed_bounds, (\d+) residual, "
                           r"(\d+) nonfinite", text)
    if doubletons:
        result.update(zip([f"samaya_{c}" for c in DOUBLETON_COUNTERS],
                          map(int, doubletons.groups())))
    # Older binaries omit this line; their substitution fields stay blank.
    substitutions = re.search(r"Presolve substitutions: (\d+) continuous, (\d+) integer; rejected: "
                              r"(\d+) row_length, (\d+) col_length, (\d+) not_implied_free, "
                              r"(\d+) by_tolerance, (\d+) pivot, (\d+) integer, "
                              r"(\d+) nonfinite, (\d+) residual, (\d+) fill", text)
    if substitutions:
        result.update(zip([f"samaya_{c}" for c in SUBSTITUTION_COUNTERS],
                          map(int, substitutions.groups())))
    if "Presolve: infeasib" in line:
        result["samaya_status"] = "infeasible"
        passes = re.search(r"passes (\d+)", line)
        if passes:
            result["samaya_passes"] = int(passes[1])
        return result
    mip = re.search(r"Presolve: (\d+) rows, (\d+) cols \((\d+) integer\) -> "
                    r"(\d+) rows, (\d+) cols \((\d+) integer\); nonzeros (\d+) -> (\d+); "
                    r"passes (\d+)", line)
    lp = re.search(r"Presolve: (\d+) rows, (\d+) cols -> (\d+) rows, (\d+) cols, "
                   r"(\d+) nonzeros \((\d+) passes\); original_nonzeros (\d+)", line)
    if mip:
        keys = ["original_rows", "original_cols", "original_integers", "samaya_rows",
                "samaya_cols", "samaya_integers", "original_nonzeros", "samaya_nonzeros",
                "samaya_passes"]
        result.update(zip(keys, map(int, mip.groups())))
    elif lp:
        keys = ["original_rows", "original_cols", "samaya_rows", "samaya_cols",
                "samaya_nonzeros", "samaya_passes", "original_nonzeros"]
        result.update(zip(keys, map(int, lp.groups())))
        result.update(original_integers=0, samaya_integers=0)
    else:
        result["samaya_error"] = "missing size log (rebuild samaya with extended presolve logs)"
        return result
    if not counters or not seconds:
        result["samaya_error"] = "missing counter or timing log"
        return result
    result["samaya_rows_removed"] = result["original_rows"] - result["samaya_rows"]
    result["samaya_cols_removed"] = result["original_cols"] - result["samaya_cols"]
    result["samaya_status"] = "reduced"
    return result


def run_samaya(exe: Path, instance: Path) -> dict:
    start = time.perf_counter()
    try:
        proc = subprocess.run([str(exe), "--json", "--time-limit", "1", "--threads", "1",
                               "--log-level", "2", str(instance)],
                              capture_output=True, text=True, timeout=61)
    except subprocess.TimeoutExpired as exc:
        captured = exc.stdout or ""
        if isinstance(captured, bytes):
            captured = captured.decode(errors="replace")
        result = parse_samaya(captured)
        result.update(samaya_solve_status="killed", samaya_wall_seconds=time.perf_counter() - start)
        if result["samaya_status"] == "failed":
            result["samaya_error"] = "wall timeout before complete presolve logs"
        return result
    result = parse_samaya(proc.stdout)
    result["samaya_wall_seconds"] = time.perf_counter() - start
    for line in reversed(proc.stdout.splitlines()):
        if line.startswith("{"):
            try:
                result["samaya_solve_status"] = json.loads(line).get("status", "unknown")
            except json.JSONDecodeError:
                result["samaya_solve_status"] = "invalid_json"
            break
    if proc.returncode:
        result.update(samaya_status="failed",
                      samaya_error=f"exit {proc.returncode}: {proc.stderr.strip()[-300:]}")
    return result


def lp_size(lp, highspy) -> dict:
    return dict(zip(SIZES, [lp.num_row_, lp.num_col_,
                           sum(t != highspy.HighsVarType.kContinuous for t in lp.integrality_),
                           len(lp.a_matrix_.value_)]))


def run_highs(instance: Path) -> dict:
    import highspy  # Optional baseline: --help and --self-test work without it.

    h = highspy.Highs()
    result = {"highs_version": h.version(), "highs_status": "failed"}
    try:
        for option, value in (("output_flag", False), ("threads", 1), ("presolve", "on"),
                              ("mip_rel_gap", 1e-4)):
            if h.setOptionValue(option, value) != highspy.HighsStatus.kOk:
                raise RuntimeError(f"HiGHS rejected {option}")
        if h.readModel(str(instance)) == highspy.HighsStatus.kError:
            raise RuntimeError("HiGHS could not read model")
        original = lp_size(h.getLp(), highspy)
        result.update({f"original_{s}": v for s, v in original.items()})
        start = time.perf_counter()
        status = h.presolve()
        result["highs_seconds"] = time.perf_counter() - start
        raw = h.getModelPresolveStatus().name.removeprefix("k")
        result["highs_status"] = {"NotReduced": "not_reduced", "Reduced": "reduced",
                                  "ReducedToEmpty": "reduced_to_empty",
                                  "Infeasible": "infeasible",
                                  "UnboundedOrInfeasible": "unbounded_or_infeasible"}.get(raw, raw)
        if status == highspy.HighsStatus.kError:
            raise RuntimeError(f"HiGHS presolve failed: {raw}")
        if raw in ("NotReduced", "Reduced", "ReducedToEmpty"):
            # NotReduced has no populated presolved model in some HiGHS versions.
            size = original if raw == "NotReduced" else lp_size(h.getPresolvedLp(), highspy)
            result.update({f"highs_{s}": v for s, v in size.items()})
    except Exception as exc:
        result.update(highs_status="failed", highs_error=str(exc))
    return result


def run_one(job: tuple[str, str, str]) -> dict:
    name, path, exe = job
    row = {"instance": name}
    for solver, runner in (("samaya", lambda: run_samaya(Path(exe), Path(path))),
                           ("highs", lambda: run_highs(Path(path)))):
        try:
            row.update(runner())
        except Exception as exc:  # Keep failed instances in the CSV and summary.
            row.update({f"{solver}_status": "failed", f"{solver}_error": str(exc)})
    if row.get("samaya_status") == "reduced" and row.get("highs_cols", 0) > 0:
        row["column_ratio"] = row["samaya_cols"] / row["highs_cols"]
    return row


def write_summary(rows: list[dict], out: Path) -> None:
    paired = [r for r in rows if r.get("samaya_status") == "reduced" and "highs_cols" in r]
    ratios = [r["column_ratio"] for r in paired if r.get("column_ratio", 0) > 0]
    geomean = f"{math.exp(sum(map(math.log, ratios)) / len(ratios)):.6g}" if ratios else "n/a"
    versions = ", ".join(sorted({r["highs_version"] for r in rows if "highs_version" in r})) or "unavailable"
    lines = ["# Presolve comparison", "", f"HiGHS version: {versions}. Threads: 1 per solver.", "",
             "samaya uses --time-limit 1 --log-level 2. Presolve times exclude reading and solving; "
             "samaya wall time includes both. Blank sizes/counters mean unavailable, not zero.", "",
             f"Paired successful size results: {len(paired)}/{len(rows)}. "
             f"Geometric mean columns samaya/HiGHS: {geomean} ({len(ratios)} positive finite ratios).",
             "Zero-column cases are excluded from the geometric mean and shown separately below.", "",
             "## Size totals (paired successful results)", "",
             "| Model | Rows | Columns | Nonzeros |", "|---|---:|---:|---:|"]
    for prefix in ("original", "samaya", "highs"):
        totals = [sum(r[f"{prefix}_{s}"] for r in paired) for s in ("rows", "cols", "nonzeros")]
        lines.append(f"| {prefix} | " + " | ".join(map(str, totals)) + " |")
    lines += ["", "Total presolve seconds on paired results: "
              f"samaya {sum(r['samaya_seconds'] for r in paired):.6g}; "
              f"HiGHS {sum(r['highs_seconds'] for r in paired):.6g}."]
    lines += ["", "## samaya rule totals (paired successful results)", "",
              "| Counter | Total |", "|---|---:|"]
    lines += [f"| {c} | {sum(r.get('samaya_' + c, 0) for r in paired)} |" for c in COUNTERS]
    lines += ["", "## Most columns retained relative to HiGHS", "",
              "| Instance | samaya columns | HiGHS columns | Ratio |", "|---|---:|---:|---:|"]
    ranked = sorted(paired, key=lambda r: (math.inf if r["highs_cols"] == 0 and r["samaya_cols"] > 0
                                          else r.get("column_ratio", 0)), reverse=True)
    for r in ranked[:10]:
        ratio = f"{r['column_ratio']:.6g}" if "column_ratio" in r else ("infinity" if r["samaya_cols"] else "0/0")
        lines.append(f"| {r['instance']} | {r['samaya_cols']} | {r['highs_cols']} | {ratio} |")
    lines += ["", "## All instances", "",
              "Sizes are rows / columns / nonzeros. Times are presolve seconds.", "",
              "| Instance | Original | samaya status | samaya size | Seconds | HiGHS status | HiGHS size | Seconds |",
              "|---|---|---|---|---:|---|---|---:|"]
    errors = []
    for r in rows:
        size = lambda p: " / ".join(str(r.get(f"{p}_{s}", "")) for s in ("rows", "cols", "nonzeros"))
        lines.append(f"| {r['instance']} | {size('original')} | {r['samaya_status']} | {size('samaya')} | "
                     f"{r.get('samaya_seconds', '')} | {r['highs_status']} | {size('highs')} | {r.get('highs_seconds', '')} |")
        for solver in ("samaya", "highs"):
            if r.get(f"{solver}_error"):
                error = str(r[f"{solver}_error"]).replace("\n", " ").replace("|", "\\|")
                errors.append(f"- {r['instance']} {solver}: {error}")
    if errors:
        lines += ["", "## Run errors", ""] + errors
    (out / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def self_test() -> None:
    # Captured from tests/instances/{tiny_lp,mrpl_crude_small}.mps, Release/MSVC.
    captured_counters = ("Presolve: rows: 0 empty, 0 singleton, 0 forcing, 0 redundant, "
                         "1 duplicate; cols: 0 fixed, 0 empty, 0 dominated, 0 free singleton; "
                         "0 bounds tightened; implied_free_col_singletons 0; "
                         "implied_free_rejected_by_tolerance 0\n")
    captured_mip = parse_samaya(
        "Presolve: 164 rows, 251 cols (83 integer) -> 163 rows, 251 cols (83 integer); "
        "nonzeros 793 -> 787; passes 2; presolve_seconds 0.0003044\n" + captured_counters)
    assert captured_mip["samaya_status"] == "reduced" and captured_mip["samaya_nonzeros"] == 787
    captured_lp = parse_samaya(
        "Presolve: 2 rows, 2 cols -> 2 rows, 2 cols, 4 nonzeros (1 passes); "
        "original_nonzeros 4; presolve_seconds 8.25e-05\n"
        + captured_counters.replace("1 duplicate", "0 duplicate"))
    assert captured_lp["samaya_status"] == "reduced" and captured_lp["samaya_cols"] == 2
    # Synthetic nonzero counts ensure every counter maps to the correct field.
    counters = ("Presolve: rows: 1 empty, 2 singleton, 3 forcing, 4 redundant, 5 duplicate; "
                "cols: 6 fixed, 7 empty, 8 dominated, 9 free singleton; 10 bounds tightened; "
                "implied_free_col_singletons 11; implied_free_rejected_by_tolerance 12\n")
    mip = parse_samaya("Presolve: 20 rows, 30 cols (15 integer) -> 4 rows, 5 cols (3 integer); "
                        "nonzeros 80 -> 9; passes 2; presolve_seconds 0.001\n" + counters)
    assert mip["samaya_status"] == "reduced" and mip["samaya_nonzeros"] == 9
    assert mip["samaya_rows_removed"] == 16 and mip["samaya_cols_removed"] == 25
    assert mip["samaya_implied_free_col_singletons"] == 11
    assert mip["samaya_implied_free_rejected_by_tolerance"] == 12
    assert "samaya_doubleton_integer_cols" not in mip
    doubletons = ("Presolve doubletons: 13 continuous, 14 integer; rejected: 15 pivot, "
                  "16 integer_data, 17 divisibility, 18 integer_bound, 19 crossed_bounds, "
                  "20 residual, 21 nonfinite\n")
    mip = parse_samaya("Presolve: 20 rows, 30 cols (15 integer) -> 4 rows, 5 cols (3 integer); "
                        "nonzeros 80 -> 9; passes 2; presolve_seconds 0.001\n" + counters
                        + doubletons)
    assert [mip[f"samaya_{c}"] for c in DOUBLETON_COUNTERS] == list(range(13, 22))
    assert mip["samaya_implied_free_rejected_by_tolerance"] == 12
    assert all(f"samaya_{c}" not in mip for c in SUBSTITUTION_COUNTERS)
    substitutions = ("Presolve substitutions: 22 continuous, 23 integer; rejected: 24 row_length, "
                     "25 col_length, 26 not_implied_free, 27 by_tolerance, 28 pivot, 29 integer, "
                     "30 nonfinite, 31 residual, 32 fill\n")
    for mode in ("Presolve: 20 rows, 30 cols (15 integer) -> 4 rows, 5 cols (3 integer); "
                 "nonzeros 80 -> 9; passes 2; presolve_seconds 0.001\n",
                 "Presolve: 20 rows, 30 cols -> 4 rows, 5 cols, 9 nonzeros (2 passes); "
                 "original_nonzeros 80; presolve_seconds 0.001\n",
                 "Presolve: infeasible; presolve_seconds 0.001; passes 2\n"):
        parsed = parse_samaya(mode + counters + doubletons + substitutions)
        assert [parsed[f"samaya_{c}"] for c in SUBSTITUTION_COUNTERS] == list(range(22, 33))
        assert [parsed[f"samaya_{c}"] for c in DOUBLETON_COUNTERS] == list(range(13, 22))
        assert parsed["samaya_implied_free_rejected_by_tolerance"] == 12
    lp = parse_samaya("Presolve: 20 rows, 30 cols -> 4 rows, 5 cols, 9 nonzeros (2 passes); "
                       "original_nonzeros 80; presolve_seconds 2e-3\n" + counters)
    assert lp["samaya_status"] == "reduced" and lp["samaya_integers"] == 0
    assert lp["samaya_seconds"] == 0.002 and lp["original_nonzeros"] == 80
    for line in ("Presolve: infeasible; presolve_seconds 0.1",
                 "Presolve: infeasibility detected; solving the original model for a certificate; presolve_seconds 0.1"):
        assert parse_samaya(line)["samaya_status"] == "infeasible"
    infeasible = parse_samaya("samaya 0.1.0: MILP with 20 rows, 30 cols, 80 nonzeros\n"
                              "  rows 20 (0 equality, 0 ranged), cols 30 (15 integer, 0 binary, 0 free), nonzeros 80\n"
                              "Presolve: infeasible; presolve_seconds 0.1; passes 2\n" + counters)
    assert infeasible["original_nonzeros"] == 80 and infeasible["original_integers"] == 15
    assert infeasible["samaya_passes"] == 2 and infeasible["samaya_fixed_cols"] == 6
    assert parse_samaya("")["samaya_status"] == "failed"
    assert parse_samaya("Presolve: 2 rows, 3 cols (1 integer) -> 1 rows, 2 cols (1 integer)")["samaya_status"] == "failed"
    assert parse_samaya("Presolve: 20 rows, 30 cols -> 4 rows, 5 cols, 9 nonzeros (2 passes); "
                        "original_nonzeros 80; presolve_seconds 0.1")["samaya_status"] == "failed"
    print("self-test: PASS (LP, MIP, C1/C2/C3 counters, infeasibility, missing and legacy logs)")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--instances", help="comma-separated names from bench/miplib_small.test")
    parser.add_argument("--jobs", type=int, default=1, help="parallel processes (default 1)")
    parser.add_argument("--out", type=Path, default=REPO / "bench" / "results" / "presolve_compare")
    parser.add_argument("--samaya", type=Path, default=DEFAULT_SAMAYA)
    parser.add_argument("--instance-dir", type=Path, default=REPO / "bench" / "instances" / "miplib_small")
    parser.add_argument("--self-test", action="store_true", help="test log parsing without solvers or MPS files")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return 0
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    names = [line.strip() for line in (REPO / "bench" / "miplib_small.test").read_text().splitlines()
             if line.strip() and not line.lstrip().startswith("#")]
    if args.instances is not None:
        selected = {n.strip() for n in args.instances.split(",")}
        unknown = selected - set(names)
        if unknown:
            parser.error("unknown instances: " + ", ".join(sorted(unknown)))
        names = [n for n in names if n in selected]
    jobs = [(n, str(args.instance_dir / f"{n}.mps"), str(args.samaya)) for n in names]
    args.out.mkdir(parents=True, exist_ok=True)
    rows = []
    with (args.out / "runs.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        writer.writeheader()
        context = (multiprocessing.Pool(args.jobs, maxtasksperchild=1)
                   if args.jobs > 1 else nullcontext())
        with context as pool:
            results = pool.imap_unordered(run_one, jobs) if pool else map(run_one, jobs)
            for k, row in enumerate(results, 1):
                rows.append(row)
                writer.writerow({key: row.get(key, "") for key in FIELDS})
                f.flush()
                print(f"[{k}/{len(jobs)}] {row['instance']}: samaya {row['samaya_status']}, "
                      f"HiGHS {row['highs_status']}", flush=True)
    order = {name: k for k, name in enumerate(names)}
    write_summary(sorted(rows, key=lambda r: order[r["instance"]]), args.out)
    print(f"results: {args.out / 'runs.csv'}; summary: {args.out / 'summary.md'}")
    return int(any(r.get(f"{s}_status") == "failed" for r in rows for s in ("samaya", "highs")))


if __name__ == "__main__":
    sys.exit(main())
