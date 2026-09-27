#!/usr/bin/env bash
# Smoke test for cases/report.py: solve a generated MRPL case, write both report formats and
# check that every verification row passed (samaya's and the report's own re-check).
#   tests/report_smoke.sh SAMAYA SOURCE_DIR WORK_DIR
set -euo pipefail
samaya="$1"
src="$2"
dir="$3"
mkdir -p "${dir}"
python3 "${src}/cases/mrpl.py" generate --out "${dir}" > /dev/null
model="${dir}/mrpl_crude_small.mps"
"${samaya}" --json --solution "${dir}/crude.sol" "${model}" > "${dir}/crude.log"
python3 "${src}/cases/report.py" "${model}" "${dir}/crude.sol" --run "${dir}/crude.log" \
  --format html --out "${dir}/crude.html"
python3 "${src}/cases/report.py" "${model}" "${dir}/crude.sol" --run "${dir}/crude.log" \
  --format excel --out "${dir}/crude.xlsx"
grep -q "VERIFIED" "${dir}/crude.html"
if grep -q 'class="cross"' "${dir}/crude.html"; then
  echo "a verification check failed" >&2
  exit 1
fi
python3 - "${dir}/crude.xlsx" <<'PY'
import sys, zipfile
z = zipfile.ZipFile(sys.argv[1])
assert z.testzip() is None
assert "xl/worksheets/sheet1.xml" in z.namelist()
PY
echo "report smoke test passed"
