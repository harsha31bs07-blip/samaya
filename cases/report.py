#!/usr/bin/env python3
"""Run report for a samaya solve, in one of two formats:

    html   a self-contained web page: run summary, verification, the plan as charts and tables,
                  the solver log, and buttons to download the tables (CSV) and the Excel workbook;
    excel  a workbook (.xlsx) with the same content, one sheet per table.

    python3 cases/report.py MODEL.mps SOLUTION.sol [--run RUN.log] [--format html|excel] [--out FILE]

SOLUTION is the file written by `samaya --solution`; RUN is samaya's output with --json (the log,
whose last line is the JSON run summary). The report also re-checks the solution on its own: it
reads MODEL, recomputes every row activity, bound and the objective from SOLUTION, and shows the
largest violations next to samaya's verifier. Only the Python standard library is used.
"""
from __future__ import annotations

import argparse
import base64
import datetime as dt
import html
import io
import json
import math
import sys
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from mrpl import read_solution  # noqa: E402

# The verifier's default tolerances (include/samaya/verify.hpp, VerifyTolerances).
kPrimalTol = 1e-6
kIntegralityTol = 1e-6
kDualTol = 1e-6
# Values below this are shown as zero in tables and dropped from the "nonzero" lists.
kZero = 1e-9

PALETTE = ["#F28C28", "#1B3A6B", "#2E9E4F", "#7B4DBF", "#1565C0", "#D64545", "#0F9D9A", "#8D6E63",
                      "#C2185B", "#607D8B"]
COLORS: dict[str, str] = {}  # one colour per series name, the same in every chart of a report


def color(name: str) -> str:
    if name not in COLORS:
        COLORS[name] = PALETTE[len(COLORS) % len(PALETTE)]
    return COLORS[name]


# ------------------------------------------------------------------------------------------------
# Reading the model, the solution and the run.

class Mps:
    def __init__(self, path: Path):
        self.name = path.stem
        self.sense = "MIN"
        self.obj = None
        self.row_type: dict[str, str] = {}
        self.cols: dict[str, dict[str, float]] = {}
        self.integer: set[str] = set()
        self.rhs: dict[str, float] = {}
        self.ranges: dict[str, float] = {}
        self.lo: dict[str, float] = {}
        self.up: dict[str, float] = {}
        section, intmode = None, False
        for raw in path.read_text().splitlines():
            if not raw.strip() or raw.startswith("*"):
                continue
            tok = raw.split()
            if not raw[0].isspace():
                section = tok[0]
                if section == "NAME" and len(tok) > 1:
                    self.name = tok[1]
                if section == "OBJSENSE" and len(tok) > 1:
                    self.sense = tok[1]
                continue
            if section == "OBJSENSE":
                self.sense = tok[0]
            elif section == "ROWS":
                self.row_type[tok[1]] = tok[0]
                if tok[0] == "N" and self.obj is None:
                    self.obj = tok[1]
            elif section == "COLUMNS":
                if len(tok) >= 3 and tok[1] == "'MARKER'":
                    intmode = tok[2] == "'INTORG'"
                    continue
                col = self.cols.setdefault(tok[0], {})
                if intmode:
                    self.integer.add(tok[0])
                for k in range(1, len(tok) - 1, 2):
                    col[tok[k]] = float(tok[k + 1])
            elif section in ("RHS", "RANGES"):
                pairs = tok[1:] if len(tok) % 2 == 1 else tok
                target = self.rhs if section == "RHS" else self.ranges
                for k in range(0, len(pairs) - 1, 2):
                    target[pairs[k]] = float(pairs[k + 1])
            elif section == "BOUNDS":
                kind, col = tok[0], tok[2]
                val = float(tok[3]) if len(tok) > 3 else None
                if kind in ("UP", "UI"):
                    self.up[col] = val
                elif kind in ("LO", "LI"):
                    self.lo[col] = val
                elif kind == "FX":
                    self.lo[col] = self.up[col] = val
                elif kind == "FR":
                    self.lo[col], self.up[col] = -math.inf, math.inf
                elif kind == "MI":
                    self.lo[col] = -math.inf
                elif kind == "PL":
                    self.up[col] = math.inf
                elif kind == "BV":
                    self.lo[col], self.up[col] = 0.0, 1.0
                    self.integer.add(col)

    def row_bounds(self, r: str) -> tuple[float, float]:
        t, b = self.row_type[r], self.rhs.get(r, 0.0)
        rng = self.ranges.get(r)
        if t == "E":
            if rng is None:
                return b, b
            return (b, b + abs(rng)) if rng > 0 else (b - abs(rng), b)
        if t == "L":
            return (b - abs(rng) if rng is not None else -math.inf), b
        if t == "G":
            return b, (b + abs(rng) if rng is not None else math.inf)
        return -math.inf, math.inf


def recheck(model: Mps, x: dict[str, float]) -> dict:
    """Recomputes feasibility and the objective from the model file alone."""
    act: dict[str, float] = {r: 0.0 for r in model.row_type}
    mag: dict[str, float] = {r: 0.0 for r in model.row_type}
    missing = 0
    for c, entries in model.cols.items():
        if c not in x:
            missing += 1
        v = x.get(c, 0.0)
        for r, a in entries.items():
            act[r] += a * v
            mag[r] = max(mag[r], abs(a * v))
    worst_row, worst_row_name = 0.0, ""
    for r, t in model.row_type.items():
        if t == "N":
            continue
        lo, up = model.row_bounds(r)
        viol = max(lo - act[r], act[r] - up, 0.0)
        rel = viol / (1.0 + max(abs(lo) if math.isfinite(lo) else 0.0,
                                                        abs(up) if math.isfinite(up) else 0.0, mag[r]))
        if rel > worst_row:
            worst_row, worst_row_name = rel, r
    worst_bound, worst_int = 0.0, 0.0
    for c in model.cols:
        v = x.get(c, 0.0)
        lo = model.lo.get(c, 0.0)
        up = model.up.get(c, math.inf)
        worst_bound = max(worst_bound, (lo - v) / (1 + abs(lo)) if v < lo else 0.0,
                                            (v - up) / (1 + abs(up)) if v > up else 0.0)
        if c in model.integer:
            worst_int = max(worst_int, abs(v - round(v)))
    objective = act[model.obj] - model.rhs.get(model.obj, 0.0) if model.obj else float("nan")
    return {"rows": sum(1 for t in model.row_type.values() if t != "N"), "cols": len(model.cols),
                    "integers": len(model.integer), "missing": missing, "row": worst_row,
                    "row_name": worst_row_name, "bound": worst_bound, "integrality": worst_int,
                    "objective": objective, "activity": act}


def read_run(path: Path | None) -> tuple[dict, str]:
    if path is None or not path.exists():
        return {}, ""
    text = path.read_text()
    summary = {}
    lines = text.rstrip("\n").split("\n")
    if lines and lines[-1].startswith("{"):
        try:
            summary = json.loads(lines[-1])
            lines = lines[:-1]
        except json.JSONDecodeError:
            pass
    return summary, "\n".join(lines)


# ------------------------------------------------------------------------------------------------
# The plan, in the MRPL case's terms.

def series(x: dict[str, float], prefix: str) -> dict[str, list[float]]:
    out: dict[str, list[float]] = {}
    for k, v in x.items():
        if k.startswith(prefix):
            key, t = k[len(prefix):].rsplit("_", 1)
            if not t.isdigit():
                continue
            row = out.setdefault(key, [])
            idx = int(t)
            while len(row) <= idx:
                row.append(0.0)
            row[idx] = v
    return {k: v for k, v in out.items() if any(abs(a) > 1e-6 for a in v)}


class Table:
    def __init__(self, title: str, header: list[str], rows: list[list], note: str = ""):
        self.title, self.header, self.rows, self.note = title, header, rows, note


def period_table(title: str, data: dict[str, list[float]], unit: str, first: str) -> Table:
    n = max((len(v) for v in data.values()), default=0)
    rows = [[k] + [round(v[t], 3) + 0.0 if t < len(v) else 0.0 for t in range(n)] for k, v in sorted(data.items())]
    return Table(title, [first] + [f"{unit} {t + 1}" for t in range(n)], rows)


def plan_content(model: Mps, x: dict[str, float], y: dict[str, float]) -> tuple[list, list[Table]]:
    """Returns (charts, tables) for the recognized MRPL case families."""
    name = model.name.upper()
    charts, tables = [], []
    if "CRUDE" in name:
        arrivals = []
        for k, v in x.items():
            if (k.startswith("arrive_") or k.startswith("spot_")) and v > 0.5:
                parts = k.split("_")
                kind = "term" if parts[0] == "arrive" else "spot"
                arrivals.append((int(parts[-1]) + 1, f"{kind} cargo {parts[1]}", "_".join(parts[2:-1])))
        arrivals.sort()
        days = max((len(v) for v in series(x, "tank_").values()), default=14)
        skipped = sorted("cargo {} ({})".format(*k[len("skip_"):].split("_", 1))
                         for k, v in x.items() if k.startswith("skip_") and v > 0.5)
        charge = series(x, "charge_")
        tank = series(x, "tank_")
        charts.append(("Berth schedule: which cargo berths on which day", svg_gantt(arrivals, days)))
        charts.append(("CDU charge by crude (kt per day)", svg_stacked(charge, days, "day")))
        charts.append(("Crude tank stock (kt, end of day)", svg_lines(tank, days, "day")))
        tables.append(Table("Berth schedule", ["Day", "Cargo", "Crude"], [list(a) for a in arrivals],
                                                ("Not lifted: " + ", ".join(skipped)) if skipped else ""))
        tables.append(period_table("CDU charge (kt per day)", charge, "Day", "Crude"))
        tables.append(period_table("Tank stock (kt, end of day)", tank, "Day", "Crude"))
    elif "PLAN" in name:
        run = series(x, "run_")
        weeks = max((len(v) for v in run.values()), default=4)
        charts.append(("Crude processed (kt per week)", svg_stacked(run, weeks, "week")))
        caps = sorted(((k, abs(v)) for k, v in y.items() if "capacity" in k and abs(v) > 1e-6),
                                    key=lambda kv: -kv[1])
        if caps:
            charts.append(("Marginal value of capacity (lakh Rs per extra kt), binding limits",
                                          svg_hbars(caps[:12])))
        tables.append(period_table("Crude processed (kt per week)", run, "Week", "Crude"))
        tables.append(period_table("Unit feed (kt per week)", series(x, "feed_"), "Week", "Unit"))
        tables.append(period_table("Domestic sales (kt per week)", series(x, "sell_"), "Week", "Product"))
        tables.append(period_table("Exports (kt per week)", series(x, "export_"), "Week", "Product"))
        if caps:
            tables.append(Table("Marginal value of capacity", ["Constraint", "Lakh Rs per extra kt"],
                                                    [[k, round(v, 2)] for k, v in caps]))
    elif "UTILITY" in name:
        allv = series(x, "")
        on = {k[:-3]: v for k, v in allv.items() if k.endswith("_on")}
        hours = max((len(v) for v in on.values()), default=24)
        charts.append(("Unit commitment (on/off by hour)", svg_onoff(on, hours)))
        flows = {k: v for k, v in allv.items() if k.endswith(("_steam", "_mw")) and not k.startswith("stg")}
        charts.append(("Steam (t/h) and power (MW) by hour", svg_lines(flows, hours, "hour")))
        tables.append(period_table("Unit commitment (1 = on)", on, "Hour", "Unit"))
        tables.append(period_table("Steam and power", flows, "Hour", "Flow"))
    return charts, tables


# ------------------------------------------------------------------------------------------------
# SVG charts (inline, no scripts).

def _fmt(v: float) -> str:
    return f"{v:,.0f}" if abs(v) >= 100 else f"{v:,.1f}"


def _axes(w, h, pad, ymax, n, label):
    out = []
    for k in range(5):
        y = pad[1] + (h - pad[1] - pad[3]) * (1 - k / 4)
        out.append(f'<line x1="{pad[0]}" y1="{y:.1f}" x2="{w - pad[2]}" y2="{y:.1f}" class="grid"/>')
        out.append(f'<text x="{pad[0] - 8}" y="{y + 4:.1f}" class="ax" text-anchor="end">{_fmt(ymax * k / 4)}</text>')
    step = max(1, n // 14)
    for t in range(0, n, step):
        x = pad[0] + (w - pad[0] - pad[2]) * (t + 0.5) / n
        out.append(f'<text x="{x:.1f}" y="{h - pad[3] + 18}" class="ax" text-anchor="middle">{t + 1}</text>')
    out.append(f'<text x="{(w + pad[0] - pad[2]) / 2}" y="{h - 6}" class="ax" text-anchor="middle">{label}</text>')
    return out


def _legend(names, w):
    out, x = [], 10
    for i, nm in enumerate(names):
        out.append(f'<rect x="{x}" y="6" width="12" height="12" rx="3" fill="{color(nm)}"/>')
        out.append(f'<text x="{x + 18}" y="16" class="lg">{html.escape(nm)}</text>')
        x += 30 + 7.2 * len(nm)
    return out


def _svg(w, h, body):
    return f'<svg viewBox="0 0 {w} {h}" class="chart" role="img">' + "".join(body) + "</svg>"


def svg_stacked(data, n, label, w=900, h=330):
    pad = (64, 34, 16, 44)
    names = sorted(data)
    tot = [sum(data[k][t] if t < len(data[k]) else 0 for k in names) for t in range(n)]
    ymax = max(tot + [1e-9]) * 1.1
    body = _axes(w, h, pad, ymax, n, label) + _legend(names, w)
    bw = (w - pad[0] - pad[2]) / n * 0.7
    for t in range(n):
        x = pad[0] + (w - pad[0] - pad[2]) * (t + 0.5) / n - bw / 2
        base = 0.0
        for i, k in enumerate(names):
            v = max(0.0, data[k][t] if t < len(data[k]) else 0.0)
            if v <= 0:
                continue
            y0 = pad[1] + (h - pad[1] - pad[3]) * (1 - (base + v) / ymax)
            hh = (h - pad[1] - pad[3]) * v / ymax
            body.append(f'<rect x="{x:.1f}" y="{y0:.1f}" width="{bw:.1f}" height="{hh:.1f}" '
                                    f'fill="{color(k)}"><title>{html.escape(k)}, {label} {t + 1}: {v:,.1f}</title></rect>')
            base += v
    return _svg(w, h, body)


def svg_lines(data, n, label, w=900, h=330):
    pad = (64, 34, 16, 44)
    names = sorted(data)
    ymax = max([max(v) for v in data.values()] + [1e-9]) * 1.1
    body = _axes(w, h, pad, ymax, n, label) + _legend(names, w)
    for i, k in enumerate(names):
        pts = []
        for t, v in enumerate(data[k]):
            x = pad[0] + (w - pad[0] - pad[2]) * (t + 0.5) / n
            y = pad[1] + (h - pad[1] - pad[3]) * (1 - max(0.0, v) / ymax)
            pts.append(f"{x:.1f},{y:.1f}")
        body.append(f'<polyline points="{" ".join(pts)}" fill="none" stroke="{color(k)}" '
                                f'stroke-width="2.5" stroke-linejoin="round"/>')
    return _svg(w, h, body)


def svg_hbars(items, w=900):
    h = 40 + 30 * len(items)
    vmax = max(v for _, v in items) * 1.1
    body = []
    for i, (k, v) in enumerate(items):
        y = 20 + i * 30
        bw = (w - 330) * v / vmax
        body.append(f'<text x="250" y="{y + 15}" class="ax" text-anchor="end">{html.escape(k)}</text>')
        body.append(f'<rect x="260" y="{y + 2}" width="{bw:.1f}" height="18" rx="4" fill="{PALETTE[0]}"/>')
        body.append(f'<text x="{266 + bw:.1f}" y="{y + 16}" class="ax">{_fmt(v)}</text>')
    return _svg(w, h, body)


def svg_gantt(arrivals, days, w=900, h=210):
    pad_l, pad_r = 40, 20
    crudes = sorted({c for _, _, c in arrivals})
    body = _legend(crudes, w)
    y = 110
    body.append(f'<rect x="{pad_l}" y="{y - 5}" width="{w - pad_l - pad_r}" height="10" rx="5" fill="#E2E8F0"/>')
    for d in range(1, days + 1):
        x = pad_l + (w - pad_l - pad_r) * (d - 0.5) / days
        body.append(f'<text x="{x:.1f}" y="{y + 40}" class="ax" text-anchor="middle">{d}</text>')
    for d, cargo, crude in arrivals:
        x = pad_l + (w - pad_l - pad_r) * (d - 0.5) / days
        body.append(f'<circle cx="{x:.1f}" cy="{y}" r="15" fill="{color(crude)}" stroke="#fff" stroke-width="3">'
                                f'<title>day {d}: {html.escape(cargo)} ({html.escape(crude)})</title></circle>')
        body.append(f'<text x="{x:.1f}" y="{y - 26}" class="ax" text-anchor="middle">{html.escape(cargo.split()[-1])}</text>')
    body.append(f'<text x="{w / 2}" y="{h - 8}" class="ax" text-anchor="middle">day (one mooring, at most one cargo a day)</text>')
    return _svg(w, h, body)


def svg_onoff(on, n, w=900):
    names = sorted(on)
    h = 50 + 34 * len(names)
    cw = (w - 120) / n
    body = []
    for i, k in enumerate(names):
        y = 20 + i * 34
        body.append(f'<text x="110" y="{y + 18}" class="ax" text-anchor="end">{html.escape(k)}</text>')
        for t in range(n):
            v = on[k][t] if t < len(on[k]) else 0.0
            body.append(f'<rect x="{120 + t * cw:.1f}" y="{y}" width="{cw - 2:.1f}" height="26" rx="3" '
                                    f'fill="{"#2E9E4F" if v > 0.5 else "#E2E8F0"}"/>')
    return _svg(w, h, body)


# ------------------------------------------------------------------------------------------------
# Excel (.xlsx) and CSV, from the standard library.

def _col(i: int) -> str:
    s = ""
    i += 1
    while i:
        i, r = divmod(i - 1, 26)
        s = chr(65 + r) + s
    return s


def _cell(ref: str, v, bold: bool = False) -> str:
    style = ' s="1"' if bold else ""
    if isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v):
        return f'<c r="{ref}"{style}><v>{v!r}</v></c>'
    return f'<c r="{ref}" t="inlineStr"{style}><is><t xml:space="preserve">{html.escape(str(v))}</t></is></c>'


def xlsx_bytes(sheets: list[tuple[str, list[list]]]) -> bytes:
    """A minimal workbook: one sheet per (name, rows); the first row of each sheet is bold."""
    buf = io.BytesIO()
    used = set()
    names = []
    for name, _ in sheets:
        base = "".join(ch for ch in name if ch not in '[]:*?/\\')[:31] or "Sheet"
        nm, k = base, 2
        while nm.lower() in used:
            nm = f"{base[:28]} {k}"
            k += 1
        used.add(nm.lower())
        names.append(nm)
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("[Content_Types].xml",
                              '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                              '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
                              '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
                              '<Default Extension="xml" ContentType="application/xml"/>'
                              '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
                              '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>'
                              + "".join(f'<Override PartName="/xl/worksheets/sheet{i + 1}.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
                                                  for i in range(len(sheets))) + "</Types>")
        z.writestr("_rels/.rels",
                              '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                              '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
                              '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>'
                              "</Relationships>")
        z.writestr("xl/workbook.xml",
                              '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                              '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
                              'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets>'
                              + "".join(f'<sheet name="{html.escape(n)}" sheetId="{i + 1}" r:id="rId{i + 1}"/>' for i, n in enumerate(names))
                              + "</sheets></workbook>")
        z.writestr("xl/_rels/workbook.xml.rels",
                              '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                              '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
                              + "".join(f'<Relationship Id="rId{i + 1}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet{i + 1}.xml"/>'
                                                  for i in range(len(sheets)))
                              + f'<Relationship Id="rId{len(sheets) + 1}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>'
                              "</Relationships>")
        z.writestr("xl/styles.xml",
                              '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                              '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
                              '<fonts count="2"><font><sz val="11"/><name val="Calibri"/></font>'
                              '<font><b/><sz val="11"/><name val="Calibri"/></font></fonts>'
                              '<fills count="2"><fill><patternFill patternType="none"/></fill><fill><patternFill patternType="gray125"/></fill></fills>'
                              '<borders count="1"><border/></borders>'
                              '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>'
                              '<cellXfs count="2"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>'
                              '<xf numFmtId="0" fontId="1" fillId="0" borderId="0" xfId="0" applyFont="1"/></cellXfs>'
                              "</styleSheet>")
        for i, (_, rows) in enumerate(sheets):
            widths = {}
            for r in rows:
                for c, v in enumerate(r):
                    widths[c] = min(60, max(widths.get(c, 8), len(str(v)) + 2))
            cols = "".join(f'<col min="{c + 1}" max="{c + 1}" width="{w}" customWidth="1"/>' for c, w in sorted(widths.items()))
            data = "".join(f'<row r="{ri + 1}">' + "".join(_cell(f"{_col(ci)}{ri + 1}", v, ri == 0) for ci, v in enumerate(r)) + "</row>"
                                          for ri, r in enumerate(rows))
            z.writestr(f"xl/worksheets/sheet{i + 1}.xml",
                                  '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                                  '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
                                  + (f"<cols>{cols}</cols>" if cols else "") + f"<sheetData>{data}</sheetData></worksheet>")
    return buf.getvalue()


def csv_text(rows: list[list]) -> str:
    def q(v):
        s = str(v)
        return '"' + s.replace('"', '""') + '"' if any(ch in s for ch in ',"\n') else s
    return "\n".join(",".join(q(v) for v in r) for r in rows) + "\n"


# ------------------------------------------------------------------------------------------------
# Assembling the report.

def build(model_path: Path, sol_path: Path, run_path: Path | None) -> dict:
    model = Mps(model_path)
    header, x, y = read_solution(sol_path)
    summary, log = read_run(run_path)
    chk = recheck(model, x)
    charts, tables = plan_content(model, x, y)
    reported = float(header.get("objective", "nan"))
    obj_diff = abs(chk["objective"] - reported) / (1 + abs(reported)) if math.isfinite(reported) else float("nan")
    bound = summary.get("dual_bound")
    gap = (abs(reported - bound) / max(1e-9, abs(reported))) if bound is not None and math.isfinite(reported) else None
    status = header.get("status", summary.get("status", "unknown"))
    verified = header.get("verified", "yes" if summary.get("verified") else "no") in ("yes", "true", "True")
    checks_samaya = []
    if summary:
        checks_samaya = [
                ("Primal feasibility (rows and bounds)", summary.get("max_primal_violation", float("nan")), kPrimalTol),
                ("Dual feasibility (LP certificate)" if model.integer == set() else "Dual feasibility (final LP)",
                  summary.get("max_dual_violation", float("nan")), kDualTol)]
    checks_own = [("Rows: worst relative violation", chk["row"], kPrimalTol),
                                ("Bounds: worst relative violation", chk["bound"], kPrimalTol),
                                ("Integrality: worst distance to an integer", chk["integrality"], kIntegralityTol),
                                ("Objective recomputed from the model (relative difference)", obj_diff, kPrimalTol)]
    variables = [[k, v] for k, v in x.items() if abs(v) > kZero]
    constraints = [[r, chk["activity"][r], y.get(r, 0.0)] for r, t in model.row_type.items() if t != "N"]
    info = {
            "model": model.name, "file": model_path.name, "sense": "maximize" if model.sense.upper().startswith("MAX") else "minimize",
            "status": status, "verified": verified, "objective": reported, "bound": bound, "gap": gap,
            "solve_seconds": summary.get("solve_seconds"), "read_seconds": summary.get("read_seconds"),
            "iterations": summary.get("simplex_iterations"), "nodes": summary.get("nodes"),
            "rows": chk["rows"], "cols": chk["cols"], "integers": chk["integers"],
            "generated": dt.datetime.now().strftime("%d %b %Y, %H:%M"),
    }
    summary_rows = [["Item", "Value"], ["Model", info["model"]], ["Model file", info["file"]],
                                    ["Sense", info["sense"]], ["Status", status], ["Verified by samaya", "yes" if verified else "no"],
                                    ["Objective", reported], ["Best bound", bound if bound is not None else ""],
                                    ["Gap", gap if gap is not None else ""], ["Solve time (s)", info["solve_seconds"] or ""],
                                    ["Rows", info["rows"]], ["Columns", info["cols"]], ["Integer columns", info["integers"]],
                                    ["Simplex iterations", info["iterations"] or ""], ["Branch-and-bound nodes", info["nodes"] or ""],
                                    ["Report generated", info["generated"]]]
    verify_rows = [["Check", "Checked by", "Largest error", "Tolerance", "Passed"]]
    for name, v, tol in checks_samaya:
        verify_rows.append([name, "samaya verifier", v, tol, "yes" if v <= tol else "no"])
    for name, v, tol in checks_own:
        verify_rows.append([name, "this report (from the MPS file)", v, tol, "yes" if v <= tol else "no"])
    sheets = [("Summary", summary_rows), ("Verification", verify_rows)]
    for t in tables:
        rows = [t.header] + t.rows + ([[t.note]] if t.note else [])
        sheets.append((t.title, rows))
    sheets.append(("Variables", [["Column", "Value"]] + variables))
    sheets.append(("Constraints", [["Row", "Activity", "Dual"]] + constraints))
    if log:
        sheets.append(("Solver log", [["samaya output"]] + [[ln] for ln in log.split("\n")]))
    return {"info": info, "checks_samaya": checks_samaya, "checks_own": checks_own, "charts": charts,
                    "tables": tables, "variables": variables, "constraints": constraints, "log": log,
                    "sheets": sheets, "row_name": chk["row_name"], "missing": chk["missing"]}


CSS = """
:root{--bg:#F6F9FD;--card:#fff;--ink:#1B1B3A;--navy:#1B3A6B;--grey:#6B7A90;--line:#E2E8F0;
--green:#2E9E4F;--lgreen:#DDF2E3;--red:#D64545;--orange:#F28C28}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);
font:15px/1.5 "Segoe UI",system-ui,-apple-system,Roboto,Arial,sans-serif}
header{background:var(--navy);color:#fff;padding:18px 32px;display:flex;align-items:center;gap:18px}
header .brand{font-weight:700;font-size:26px;letter-spacing:.5px}header .sub{opacity:.8}
header .when{margin-left:auto;opacity:.75;font-size:13px}
main{max-width:1180px;margin:0 auto;padding:24px 20px 60px}
h2{color:var(--navy);font-size:20px;margin:30px 0 12px}
.hero{display:flex;flex-wrap:wrap;gap:12px;align-items:center;margin-bottom:6px}
.hero h1{font-size:24px;margin:0 12px 0 0}
.pill{border-radius:999px;padding:4px 14px;font-weight:600;font-size:13px}
.ok{background:var(--lgreen);color:var(--green)}.bad{background:#FBE3E3;color:var(--red)}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:14px}
.card{background:var(--card);border-radius:14px;box-shadow:0 2px 10px rgba(27,58,107,.07);padding:16px 18px}
.kpi .v{font-size:24px;font-weight:700;color:var(--navy)}.kpi .l{color:var(--grey);font-size:13px}
.two{display:grid;grid-template-columns:1fr 1fr;gap:14px}@media(max-width:860px){.two{grid-template-columns:1fr}}
table{border-collapse:collapse;width:100%;font-size:13.5px}
th,td{padding:6px 10px;border-bottom:1px solid var(--line);text-align:right;white-space:nowrap}
th:first-child,td:first-child,td.t{text-align:left}th{color:var(--grey);font-weight:600;background:#FAFCFE}
.check td:first-child{white-space:normal}.tick{color:var(--green);font-weight:700}.cross{color:var(--red);font-weight:700}
.chart{width:100%;height:auto}.chart .grid{stroke:#EDF1F6}.chart .ax{fill:#6B7A90;font-size:12px}
.chart .lg{fill:#1B1B3A;font-size:12px}
.scroll{overflow:auto;max-height:420px}
pre.log{background:#0F1720;color:#DCE3EC;border-radius:12px;padding:16px 18px;overflow:auto;font:13px/1.45 Consolas,"Cascadia Mono",monospace}
.btn{display:inline-block;background:var(--navy);color:#fff;text-decoration:none;border-radius:10px;
padding:9px 16px;margin:4px 8px 4px 0;font-weight:600;font-size:14px}.btn.alt{background:#fff;color:var(--navy);border:1.5px solid var(--navy)}
.note{color:var(--grey);font-size:13px;margin-top:8px}
input.search{width:100%;padding:8px 12px;border:1px solid var(--line);border-radius:8px;margin:6px 0 10px;font:inherit}
details summary{cursor:pointer;color:var(--navy);font-weight:600}
footer{color:var(--grey);font-size:12.5px;text-align:center;margin-top:40px}
"""

JS = """
function filt(inp,id){const q=inp.value.toLowerCase();for(const tr of document.getElementById(id).tBodies[0].rows){
tr.style.display=tr.cells[0].textContent.toLowerCase().includes(q)?'':'none'}}
"""


def _num(v, digits=6) -> str:
    if v is None or (isinstance(v, float) and not math.isfinite(v)):
        return "–"
    if isinstance(v, float):
        if v != 0 and (abs(v) < 1e-3 or abs(v) >= 1e9):
            return f"{v:.2e}"
        return f"{v:,.{digits}g}" if abs(v) < 1e5 else f"{v:,.2f}"
    return f"{v:,}" if isinstance(v, int) else html.escape(str(v))


def _table_html(t: Table, tid: str = "", limit: int | None = None) -> str:
    rows = t.rows if limit is None else t.rows[:limit]
    head = "".join(f"<th>{html.escape(str(h))}</th>" for h in t.header)
    body = "".join("<tr>" + "".join(f"<td>{_num(v, 4)}</td>" if isinstance(v, (int, float))
                                             else f'<td class="t">{html.escape(str(v))}</td>' for v in r) + "</tr>"
                                  for r in rows)
    return f'<table id="{tid}"><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>'


def _data_uri(mime: str, data: bytes) -> str:
    return f"data:{mime};base64," + base64.b64encode(data).decode()


def render_html(rep: dict, xlsx: bytes, stem: str) -> str:
    i = rep["info"]
    ok = i["verified"] and i["status"] == "optimal"
    pills = (f'<span class="pill {"ok" if i["status"] == "optimal" else "bad"}">{html.escape(i["status"]).upper()}</span>'
                      f'<span class="pill {"ok" if i["verified"] else "bad"}">{"✓ VERIFIED" if i["verified"] else "NOT VERIFIED"}</span>')
    kpis = [("Objective (" + i["sense"] + ")", _num(i["objective"], 10)),
                    ("Best bound · gap", (_num(i["bound"], 10) + " · " + (f"{i['gap']:.2%}" if i["gap"] is not None else "–"))
                      if i["bound"] is not None else "–"),
                    ("Solve time", f"{i['solve_seconds']:.3f} s" if i["solve_seconds"] is not None else "–"),
                    ("Model size", f"{i['rows']:,} rows · {i['cols']:,} cols"),
                    ("Integer decisions", f"{i['integers']:,}"),
                    ("Iterations · nodes", f"{_num(i['iterations'])} · {_num(i['nodes'])}")]
    kpi_html = "".join(f'<div class="card kpi"><div class="l">{html.escape(l)}</div><div class="v">{v}</div></div>'
                                          for l, v in kpis)

    def check_rows(items):
        return "".join(f"<tr><td>{html.escape(n)}</td><td>{_num(v)}</td><td>{tol:.0e}</td>"
                                      f'<td class="{"tick" if v <= tol else "cross"}">{"✓" if v <= tol else "✗"}</td></tr>'
                                      for n, v, tol in items)
    head = "<tr><th>Check</th><th>Largest error</th><th>Tolerance</th><th></th></tr>"
    verify = (f'<div class="two"><div class="card"><b>samaya’s built-in verifier</b>'
                        f'<table class="check">{head}{check_rows(rep["checks_samaya"])}</table>'
                        f'<div class="note">Recorded in the run summary (samaya --json).</div></div>'
                        f'<div class="card"><b>Independent re-check by this report</b>'
                        f'<table class="check">{head}{check_rows(rep["checks_own"])}</table>'
                        f'<div class="note">Recomputed here from the model file ({html.escape(i["file"])}) and the solution file alone.</div></div></div>')
    charts = "".join(f'<div class="card"><b>{html.escape(t)}</b>{svg}</div>' for t, svg in rep["charts"])
    tables = "".join(f'<div class="card"><b>{html.escape(t.title)}</b><div class="scroll">{_table_html(t)}</div>'
                                      + (f'<div class="note">{html.escape(t.note)}</div>' if t.note else "") + "</div>"
                                      for t in rep["tables"])
    downloads = [f'<a class="btn" download="{html.escape(stem)}.xlsx" href="{_data_uri("application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", xlsx)}">⬇ Excel workbook (.xlsx)</a>']
    for name, rows in rep["sheets"]:
        downloads.append(f'<a class="btn alt" download="{html.escape(stem)}_{html.escape(name.lower().replace(" ", "_"))}.csv" '
                                          f'href="{_data_uri("text/csv", csv_text(rows).encode())}">{html.escape(name)} (CSV)</a>')
    var_t = Table("Variables", ["Column", "Value"], rep["variables"])
    con_t = Table("Constraints", ["Row", "Activity", "Dual"], rep["constraints"])
    log = f'<pre class="log">{html.escape(rep["log"])}</pre>' if rep["log"] else '<div class="note">No run log given (--run).</div>'
    return f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>samaya report · {html.escape(i["model"])}</title><style>{CSS}</style></head><body>
<header><div class="brand">samaya</div><div class="sub">Run report</div><div class="when">{html.escape(i["generated"])}</div></header>
<main>
<div class="hero"><h1>{html.escape(i["model"])}</h1>{pills}</div>
<div class="note">{"The answer is optimal and passed every check below." if ok else "The answer did not pass every check; see below."}</div>
<h2>Summary</h2><div class="grid">{kpi_html}</div>
<h2>Verification</h2>{verify}
{"<h2>Plan</h2>" + charts + tables if charts or tables else ""}
<h2>Downloads</h2><div class="card">{"".join(downloads)}<div class="note">Every table in this report, as an Excel workbook or as CSV files.</div></div>
<h2>Solver log</h2>{log}
<h2>Full solution</h2>
<div class="card"><details><summary>Variables ({len(rep["variables"]):,} nonzero)</summary>
<input class="search" placeholder="Filter by name" oninput="filt(this,'vars')"><div class="scroll">{_table_html(var_t, "vars")}</div></details></div>
<div class="card" style="margin-top:14px"><details><summary>Constraints ({len(rep["constraints"]):,} rows)</summary>
<input class="search" placeholder="Filter by name" oninput="filt(this,'cons')"><div class="scroll">{_table_html(con_t, "cons")}</div></details></div>
<footer>Generated by cases/report.py from samaya’s output files: {html.escape(i["file"])}, the solution file and the run log.</footer>
</main><script>{JS}</script></body></html>"""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("model", type=Path)
    ap.add_argument("solution", type=Path)
    ap.add_argument("--run", type=Path, help="samaya's output with --json (log plus JSON summary)")
    ap.add_argument("--format", choices=["html", "excel"], default="html",
                                    help="html: a web page report (option 1); excel: a downloadable workbook (option 2)")
    ap.add_argument("--out", type=Path)
    a = ap.parse_args()
    rep = build(a.model, a.solution, a.run)
    stem = a.model.stem
    xlsx = xlsx_bytes(rep["sheets"])
    if a.format == "excel":
        out = a.out or a.solution.with_suffix(".xlsx")
        out.write_bytes(xlsx)
    else:
        out = a.out or a.solution.with_suffix(".html")
        out.write_text(render_html(rep, xlsx, stem), encoding="utf-8")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
