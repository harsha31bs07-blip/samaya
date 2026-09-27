#!/usr/bin/env bash
# End-to-end MRPL demo: generate the case-study models, solve each with samaya (every result is
# checked by the independent verifier), and print a readable plan.
#
#   cases/demo.sh [SIZE] [SAMAYA] [--output html|excel]
#
#   SIZE      small (default), medium or large
#   SAMAYA    solver binary (default build/release/apps/cli/samaya)
#   --output  also write a run report per model (cases/report.py): html is a web page with the
#             verification, the plan as charts and download buttons; excel is a workbook
set -euo pipefail

cd "$(dirname "$0")/.."
output=""
args=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --output) output="${2:-}"; shift 2 ;;
    --output=*) output="${1#--output=}"; shift ;;
    *) args+=("$1"); shift ;;
  esac
done
if [[ -n "${output}" && "${output}" != "html" && "${output}" != "excel" ]]; then
  echo "--output must be html or excel" >&2
  exit 1
fi
size="${args[0]:-small}"
samaya="${args[1]:-build/release/apps/cli/samaya}"
out="cases/instances"
if [[ ! -x "${samaya}" ]]; then
  echo "solver not found: ${samaya}" >&2
  echo "build it with: cmake --preset release && cmake --build --preset release" >&2
  exit 1
fi

python3 cases/mrpl.py generate --out "${out}" > /dev/null
for family in plan crude utility; do
  model="${out}/mrpl_${family}_${size}.mps"
  echo "=================================================================================="
  echo "${model}"
  echo "----------------------------------------------------------------------------------"
  start=$(date +%s.%N)
  if [[ -n "${output}" ]]; then
    # The log and its JSON summary feed the report.
    "${samaya}" --time-limit 60 --json --solution "${out}/${family}_${size}.sol" "${model}" \
      > "${out}/${family}_${size}.log"
  else
    "${samaya}" --time-limit 60 --log-level 0 --solution "${out}/${family}_${size}.sol" "${model}"
  fi
  end=$(date +%s.%N)
  python3 cases/mrpl.py report "${model}" "${out}/${family}_${size}.sol"
  if [[ -n "${output}" ]]; then
    python3 cases/report.py "${model}" "${out}/${family}_${size}.sol" \
      --run "${out}/${family}_${size}.log" --format "${output}"
  fi
  printf '\nsolved in %.2f s\n\n' "$(echo "${end} - ${start}" | bc)"
done
