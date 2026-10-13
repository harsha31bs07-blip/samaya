"""Smoke test for cases/report.py: solve a generated MRPL case, write both report formats and
check that every verification row passed (samaya's verifier and the report's own re-check).
Portable (Linux and Windows).

  python3 tests/report_smoke.py SAMAYA SOURCE_DIR WORK_DIR
"""
import subprocess
import json
import math
import sys
import zipfile
from pathlib import Path


def check_ranging(samaya: Path, src: Path, work: Path) -> None:
    sys.path.insert(0, str(src / "cases"))
    from report import build, render_html, xlsx_bytes
    from mrpl import read_solution

    model = work / "tiny.mps"
    model.write_text("NAME TINY\nROWS\n N obj\n G demand\nCOLUMNS\n x obj 1 demand 1\n"
                     " y obj 2 demand 1\nRHS\n rhs demand 4\nBOUNDS\n UP bnd x 3\nENDATA\n")
    for ranged in (False, True):
        stem = "tiny_ranging" if ranged else "tiny_old"
        sol, log = work / f"{stem}.sol", work / f"{stem}.log"
        args = [str(samaya), "--lp-method", "dual", "--threads", "1", "--json", "--solution", str(sol)]
        if ranged:
            args.append("--ranging")
        output = subprocess.run(args + [str(model)], check=True, capture_output=True, text=True).stdout
        log.write_text(output, encoding="utf-8")
        summary = json.loads(output.splitlines()[-1])
        details = {}
        header, values, duals = read_solution(sol, details=details)
        assert values == {"x": 3.0, "y": 1.0} and duals == {"demand": 2.0}
        assert header["verified"] == "yes"
        if ranged:
            assert summary["ranging"]["available"] and details["ranging"]["available"]
            assert details["columns"]["x"]["cost_lower"] == -math.inf
            assert details["columns"]["y"]["cost_upper"] == math.inf
            assert math.isnan(details["columns"]["y"]["level_lower"])
            assert details["rows"]["demand"]["level_lower"] == 3.0
        else:
            assert summary["ranging"] is None and details["ranging"] is None
            assert sol.read_text() == ("status optimal\nverified yes\nobjective 5\n\ncolumns 2\n"
                                       "# name value reduced_cost\nx 3 -1\ny 1 0\n\nrows 1\n"
                                       "# name activity dual\ndemand 4 2\n")
        rep = build(model, sol, log)
        if ranged:
            assert len(rep["sensitivity"][1].rows) == 2  # x and demand; y is basic.
        for fmt in ("html", "excel"):
            target = work / f"{stem}.{'html' if fmt == 'html' else 'xlsx'}"
            subprocess.run([sys.executable, str(src / "cases" / "report.py"), str(model), str(sol),
                            "--run", str(log), "--format", fmt, "--out", str(target)], check=True)
            if fmt == "html":
                page = target.read_text(encoding="utf-8")
                assert ('<h2>Sensitivity</h2>' in page) == ranged
                assert 'class="cross"' not in page
                if ranged:
                    # Inspect the tables, not the embedded workbook's arbitrary base64 text.
                    sensitivity = page.split('<h2>Sensitivity</h2>', 1)[1].split('<h2>Downloads</h2>', 1)[0]
                    assert "−∞" in sensitivity and ">inf<" not in sensitivity and ">nan<" not in sensitivity
            else:
                with zipfile.ZipFile(target) as z:
                    assert z.testzip() is None
                    workbook = z.read("xl/workbook.xml").decode()
                    assert ('name="Sensitivity"' in workbook) == ranged
                    if ranged:
                        sheet_index = next(i for i, (name, _) in enumerate(rep["sheets"], 1) if name == "Sensitivity")
                        sheet = z.read(f"xl/worksheets/sheet{sheet_index}.xml").decode()
                        assert ">-inf<" in sheet and ">inf<" in sheet and "nan" not in sheet

    reordered = work / "reordered.sol"
    reordered.write_text("status optimal\nverified yes\nranging available alternative optimal vertex\n"
                         "objective 5\ncolumns 2\n# value name cost_upper reduced_cost cost_lower level_upper level_lower\n"
                         "3 x 2 -1 -inf 4 0\n1 y inf 0 1 nan nan\n"
                         "rows 1\n# dual level_upper name level_lower activity\n2 inf demand 3 4\n")
    rep = build(model, reordered, None)
    assert rep["ranging"]["message"] == "alternative optimal vertex"
    assert rep["sensitivity"][0].rows[0] == ["x", 3.0, -1.0, -math.inf, 2.0]
    assert rep["sensitivity"][1].rows[1] == ["demand", "Row", 4.0, 2.0, 3.0, math.inf]
    assert "alternative optimal vertex" in render_html(rep, xlsx_bytes(rep["sheets"]), "tiny")
    no_rows = work / "no_rows.mps"
    no_rows.write_text("NAME NO_ROWS\nROWS\n N obj\nCOLUMNS\n x obj 1\nENDATA\n")
    sol = work / "no_rows.sol"
    subprocess.run([str(samaya), "--ranging", "--log-level", "0", "--solution", str(sol), str(no_rows)], check=True)
    assert "# name value reduced_cost cost_lower cost_upper level_lower level_upper" in sol.read_text()
    assert "# name activity dual level_lower level_upper" in sol.read_text()
    print("ranging report formats and backward compatibility passed")


def main() -> int:
    samaya, src, work = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    work.mkdir(parents=True, exist_ok=True)
    py = sys.executable
    check_ranging(samaya, src, work)
    subprocess.run([py, str(src / "cases" / "mrpl.py"), "generate", "--out", str(work)],
                   check=True, stdout=subprocess.DEVNULL)
    model = work / "mrpl_crude_small.mps"
    sol, log = work / "crude.sol", work / "crude.log"
    with open(log, "w") as out:
        subprocess.run([str(samaya), "--ranging", "--json", "--solution", str(sol), str(model)], check=True,
                       stdout=out)
    summary = json.loads(log.read_text().splitlines()[-1])
    assert not summary["ranging"]["available"] and summary["ranging"]["message"]
    assert "ranging unavailable " in sol.read_text()
    for fmt, target in (("html", work / "crude.html"), ("excel", work / "crude.xlsx")):
        subprocess.run([py, str(src / "cases" / "report.py"), str(model), str(sol), "--run", str(log),
                        "--format", fmt, "--out", str(target)], check=True)
    page = (work / "crude.html").read_text(encoding="utf-8")
    if "VERIFIED" not in page:
        print("the report does not say VERIFIED", file=sys.stderr)
        return 1
    if 'class="cross"' in page:
        print("a verification check failed", file=sys.stderr)
        return 1
    assert "Sensitivity ranging is unavailable." in page
    with zipfile.ZipFile(work / "crude.xlsx") as z:
        if z.testzip() is not None or "xl/worksheets/sheet1.xml" not in z.namelist():
            print("the workbook is not a valid .xlsx", file=sys.stderr)
            return 1
        assert 'name="Sensitivity"' in z.read("xl/workbook.xml").decode()
    print("report smoke test passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
