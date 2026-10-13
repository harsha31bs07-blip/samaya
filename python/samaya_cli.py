"""Dependency-free helpers for samaya's command-line solver plugins."""

import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from typing import Any, Sequence


def find_binary(path: str | None = None) -> str:
    """Find samaya, respecting explicit paths and SAMAYA_BIN before defaults."""
    requested = path or os.environ.get("SAMAYA_BIN")
    if requested:
        resolved = shutil.which(str(requested))
        if resolved:
            return resolved
        raise FileNotFoundError(f"Cannot execute samaya at {requested!r}")
    for name in ("samaya", "samaya.exe"):
        resolved = shutil.which(name)
        if resolved:
            return resolved
    root = Path(__file__).resolve().parent.parent
    for folder in ("windows-release", "windows-debug", "release", "debug"):
        for name in ("samaya.exe", "samaya"):
            candidate = root / "build" / folder / "apps" / "cli" / name
            if candidate.is_file() and os.access(candidate, os.X_OK):
                return str(candidate)
    raise FileNotFoundError(
        "samaya executable not found; pass path=/executable=, set SAMAYA_BIN, "
        "put samaya on PATH, or build the CLI in a repository build folder"
    )


@dataclass
class Solution:
    """Named values from a samaya --solution file."""

    columns: dict[str, float] = field(default_factory=dict)
    reduced_costs: dict[str, float] = field(default_factory=dict)
    rows: dict[str, float] = field(default_factory=dict)
    duals: dict[str, float] = field(default_factory=dict)


def parse_summary(stdout: str) -> dict[str, Any]:
    """Read the last nonempty stdout line as the required JSON run summary."""
    try:
        result = json.loads(stdout.rstrip().splitlines()[-1])
        if not isinstance(result, dict) or "status" not in result:
            raise ValueError("missing status")
        return result
    except (IndexError, ValueError) as exc:
        raise RuntimeError("samaya did not produce a valid JSON summary") from exc


def read_solution(filename: str) -> Solution:
    """Parse primal/dual blocks, ignoring certificate and ray blocks."""
    solution = Solution()
    section = None
    remaining = 0
    with open(filename, encoding="utf-8") as stream:
        for line in stream:
            fields = line.split()
            if not fields or fields[0] == "#":
                continue
            if remaining:
                if section == "columns":
                    solution.columns[fields[0]] = float(fields[1])
                    if len(fields) > 2:
                        solution.reduced_costs[fields[0]] = float(fields[2])
                elif section == "rows":
                    solution.rows[fields[0]] = float(fields[1])
                    if len(fields) > 2:
                        solution.duals[fields[0]] = float(fields[2])
                remaining -= 1
            elif fields[0] in (
                "columns", "rows", "infeasibility_certificate", "unbounded_ray"
            ):
                section, remaining = fields[0], int(fields[1])
    if remaining:
        raise RuntimeError("Truncated samaya solution file")
    return solution


def command_line(binary: str, model: str, solution: str,
                 options: Sequence[str] = ()) -> list[str]:
    """Build an argv list without shell interpolation."""
    return [binary, *map(str, options), "--solution", solution, "--json", model]


def run_solver(binary: str, model: str, solution: str,
               options: Sequence[str] = (), msg: bool = False) -> dict[str, Any]:
    """Execute samaya and return its summary, raising on CLI failures."""
    completed = subprocess.run(
        command_line(binary, model, solution, options), capture_output=True,
        text=True, encoding="utf-8", errors="replace", check=False,
    )
    if msg:
        print(completed.stdout, end="")
        print(completed.stderr, end="", file=sys.stderr)
    if completed.returncode:
        raise RuntimeError(
            f"samaya exited with code {completed.returncode}: "
            f"{completed.stderr.strip() or completed.stdout.strip()}"
        )
    return parse_summary(completed.stdout)
