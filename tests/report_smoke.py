"""Smoke test for cases/report.py: solve a generated MRPL case, write both report formats and
check that every verification row passed (samaya's verifier and the report's own re-check).
Portable (Linux and Windows).

  python3 tests/report_smoke.py SAMAYA SOURCE_DIR WORK_DIR
"""
import subprocess
import sys
import zipfile
from pathlib import Path


def main() -> int:
    samaya, src, work = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    work.mkdir(parents=True, exist_ok=True)
    py = sys.executable
    subprocess.run([py, str(src / "cases" / "mrpl.py"), "generate", "--out", str(work)],
                   check=True, stdout=subprocess.DEVNULL)
    model = work / "mrpl_crude_small.mps"
    sol, log = work / "crude.sol", work / "crude.log"
    with open(log, "w") as out:
        subprocess.run([str(samaya), "--json", "--solution", str(sol), str(model)], check=True,
                       stdout=out)
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
    with zipfile.ZipFile(work / "crude.xlsx") as z:
        if z.testzip() is not None or "xl/worksheets/sheet1.xml" not in z.namelist():
            print("the workbook is not a valid .xlsx", file=sys.stderr)
            return 1
    print("report smoke test passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
