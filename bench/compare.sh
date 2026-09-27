#!/usr/bin/env bash
# Side-by-side comparison of samaya with the open-source solvers it competes with (HiGHS, SCIP,
# CBC, GLPK), with the same 1e-4 MIP gap and thread count, on one machine.
#
#   bench/compare.sh SET [SHARDS]
#
#   SET     netlib         Netlib LP (93), 300 s, samaya HiGHS SCIP CBC GLPK
#           cases          MRPL cases (9) + generated MILPs (7), 300 s, all five solvers
#           miplib60       MIPLIB small (62), 60 s, samaya HiGHS SCIP CBC
#           miplib600      MIPLIB small (62), 600 s, samaya HiGHS SCIP
#           miplib600-t6   MIPLIB small (62), 600 s, samaya and HiGHS with 6 threads each
#                          (HiGHS 1.15 has a parallel MIP solver; use SHARDS 2)
#           miplib600-samaya  samaya alone on the 600 s set (same instance-to-shard split)
#           maros          Maros-Meszaros convex QP (138), 300 s, samaya and HiGHS
#           netlib-barrier Netlib LP (93), 300 s, samaya with --lp-method barrier alone
#           mittelmann     Mittelmann LP subset (13, bench/mittelmann_lp.test), 1000 s, samaya
#                          and HiGHS
#   SHARDS  instances run at the same time (default 4; keep SHARDS x threads <= cores)
#
# Needs: the release build, python3 with highspy and pyscipopt, cbc and glpsol on PATH, and the
# instances (bench/fetch_instances.sh netlib; bench/fetch_instances.sh miplib-list
# bench/miplib_small.test; bench/generate_lps.py --set mip --out bench/instances/generated-mip;
# bench/fetch_instances.sh maros; bench/fetch_instances.sh mittelmann).
# SAMAYA_CLOCK_FACTOR (default 1): see bench/clock_factor.py; it gives samaya the same real time as
# the baselines on a machine whose monotonic clock runs fast.
# Writes bench/results/compare/SET/shard*.csv and prints each shard's summary; the combined
# summary is at the end.
set -euo pipefail

cd "$(dirname "$0")/.."
set_name="${1:?usage: bench/compare.sh SET [SHARDS]}"
shards="${2:-4}"
samaya="build/release/apps/cli/samaya"
out="bench/results/compare/${set_name}"
threads=1
solu=()
case "${set_name}" in
  netlib)
    files=(bench/instances/netlib/*.mps*)
    limit=300
    baselines=(highspy scip cbc glpk)
    ;;
  cases)
    python3 cases/mrpl.py generate --out cases/instances > /dev/null
    files=(cases/instances/*.mps bench/instances/generated-mip/*.mps)
    limit=300
    baselines=(highspy scip cbc glpk)
    ;;
  miplib60 | miplib600 | miplib600-t6 | miplib600-samaya)
    dir="bench/instances/miplib_small"
    mapfile -t files < <(grep -v '^#' bench/miplib_small.test | grep -v '^$' |
                         sed "s|^|${dir}/|; s|$|.mps|")
    solu=(--solu "${dir}/miplib2017.solu")
    limit=600
    baselines=(highspy scip)
    if [[ "${set_name}" == miplib60 ]]; then
      limit=60
      baselines=(highspy scip cbc)
    elif [[ "${set_name}" == miplib600-samaya ]]; then
      baselines=()  # samaya alone, to pair with the baselines of an earlier miplib600 run.
    elif [[ "${set_name}" == miplib600-t6 ]]; then
      # Same thread count for both: HiGHS 1.15 searches the tree in parallel too.
      baselines=(highspy)
      threads=6
    fi
    ;;
  maros)
    files=(bench/instances/maros/*.mps)
    limit=300
    baselines=(highspy)
    ;;
  mittelmann)
    dir="bench/instances/mittelmann"
    mapfile -t files < <(grep -v '^#' bench/mittelmann_lp.test | grep -v '^$' |
                         sed "s|^|${dir}/|; s|$|.mps|")
    limit=1000
    baselines=(highspy)
    ;;
  netlib-barrier)
    files=(bench/instances/netlib/*.mps*)
    limit=300
    baselines=()
    samaya="bench/samaya_barrier.sh"
    ;;
  *)
    echo "unknown set: ${set_name}" >&2
    exit 1
    ;;
esac
[[ -x "${samaya}" ]] || { echo "build first: cmake --preset release && cmake --build --preset release" >&2; exit 1; }
[[ -e "${files[0]}" ]] || { echo "instances missing: ${files[0]} (see the header of this script)" >&2; exit 1; }

base_args=()
for b in "${baselines[@]}"; do base_args+=(--baseline "${b}"); done
mkdir -p "${out}"
for ((k = 0; k < shards; ++k)); do
  shard=()
  for ((i = k; i < ${#files[@]}; i += shards)); do shard+=("${files[i]}"); done
  [[ ${#shard[@]} -gt 0 ]] || continue
  python3 bench/harness.py "${shard[@]}" --samaya "${samaya}" --time-limit "${limit}" \
    --samaya-clock-factor "${SAMAYA_CLOCK_FACTOR:-1}" \
    --threads "${threads}" "${base_args[@]}" "${solu[@]}" --out "${out}/shard${k}.csv" \
    > "${out}/shard${k}.log" 2>&1 &
done
wait
grep -hE "DISAGREE|WRONG" "${out}"/shard*.log || true
python3 - "${out}" "${limit}" <<'EOF'
import csv, glob, math, sys
out, limit = sys.argv[1], float(sys.argv[2])
rows = [r for f in sorted(glob.glob(out + "/shard*.csv")) for r in csv.DictReader(open(f))]
solved_states = {"optimal", "infeasible", "unbounded"}
print(f"{'solver':<10} {'solved':>8} {'shifted geomean (s)':>20}")
for solver in dict.fromkeys(r["solver"] for r in rows):
    mine = [r for r in rows if r["solver"] == solver]
    solved = sum(r["status"] in solved_states for r in mine)
    times = [float(r["wall_seconds"]) if r["status"] in solved_states else limit for r in mine]
    geo = math.exp(sum(math.log(t + 10) for t in times) / len(times)) - 10
    print(f"{solver:<10} {solved:>4}/{len(mine):<3} {geo:>20.2f}")
EOF
