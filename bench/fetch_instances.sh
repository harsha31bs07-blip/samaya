#!/usr/bin/env bash
# Downloads public benchmark instances into bench/instances/<set>/ as plain .mps files.
#
#   bench/fetch_instances.sh netlib           # Netlib LP (decoded from netlib's compressed EMPS format)
#   bench/fetch_instances.sh netlib-infeas    # Netlib infeasible LPs
#   bench/fetch_instances.sh miplib           # MIPLIB 2017 benchmark set (330 MB)
#   bench/fetch_instances.sh miplib-list FILE # MIPLIB 2017 instances named in FILE, one per line
#   bench/fetch_instances.sh maros            # Maros-Meszaros convex QP set (138, QPS format)
#   bench/fetch_instances.sh mittelmann       # Mittelmann LP subset (bench/mittelmann_lp.test)
#
# Instances are never committed (see .gitignore). Requires curl, a C compiler, gunzip and unzip.
set -euo pipefail

here="$(cd "$(dirname "$0")" && pwd)"
dest_root="${here}/instances"

netlib_names=(
  25fv47 80bau3b adlittle afiro agg agg2 agg3 bandm beaconfd blend bnl1 bnl2 boeing1 boeing2
  bore3d brandy capri cycle czprob d2q06c d6cube degen2 degen3 dfl001 e226 etamacro fffff800
  finnis fit1d fit1p fit2d fit2p forplan ganges gfrd-pnc greenbea greenbeb grow15 grow22 grow7
  israel kb2 lotfi maros maros-r7 modszk1 nesm perold pilot pilot.ja pilot.we pilot4 pilot87
  pilotnov recipe sc105 sc205 sc50a sc50b scagr25 scagr7 scfxm1 scfxm2 scfxm3
  scorpion scrs8 scsd1 scsd6 scsd8 sctap1 sctap2 sctap3 seba share1b share2b shell ship04l
  ship04s ship08l ship08s ship12l ship12s sierra stair standata standgub standmps stocfor1
  stocfor2 tuff vtp.base wood1p woodw
)

netlib_infeas_names=(
  bgdbg1 bgetam bgindy bgprtr box1 ceria3d chemcom cplex1 cplex2 ex72a ex73a forest6 galenet
  gosh gran greenbea itest2 itest6 klein1 klein2 klein3 mondou2 pang pilot4i qual reactor
  refinery vol1 woodinfe
)

# The Maros-Meszaros set from the mirror of its SVN repository (files are QPS despite the .SIF
# name). Saved as .mps: HiGHS picks its reader by extension and reads QUADOBJ from MPS files.
maros_names=(
  AUG2D AUG2DC AUG2DCQP AUG2DQP AUG3D AUG3DC AUG3DCQP AUG3DQP BOYD1 BOYD2 CONT-050 CONT-100
  CONT-101 CONT-200 CONT-201 CONT-300 CVXQP1_L CVXQP1_M CVXQP1_S CVXQP2_L CVXQP2_M CVXQP2_S
  CVXQP3_L CVXQP3_M CVXQP3_S DPKLO1 DTOC3 DUAL1 DUAL2 DUAL3 DUAL4 DUALC1 DUALC2 DUALC5 DUALC8
  EXDATA GENHS28 GOULDQP2 GOULDQP3 HS118 HS21 HS268 HS35 HS35MOD HS51 HS52 HS53 HS76 HUES-MOD
  HUESTIS KSIP LASER LISWET1 LISWET10 LISWET11 LISWET12 LISWET2 LISWET3 LISWET4 LISWET5 LISWET6
  LISWET7 LISWET8 LISWET9 LOTSCHD MOSARQP1 MOSARQP2 POWELL20 PRIMAL1 PRIMAL2 PRIMAL3 PRIMAL4
  PRIMALC1 PRIMALC2 PRIMALC5 PRIMALC8 Q25FV47 QADLITTL QAFIRO QBANDM QBEACONF QBORE3D QBRANDY
  QCAPRI QE226 QETAMACR QFFFFF80 QFORPLAN QGFRDXPN QGROW15 QGROW22 QGROW7 QISRAEL QPCBLEND
  QPCBOEI1 QPCBOEI2 QPCSTAIR QPILOTNO QPTEST QRECIPE QSC205 QSCAGR25 QSCAGR7 QSCFXM1 QSCFXM2
  QSCFXM3 QSCORPIO QSCRS8 QSCSD1 QSCSD6 QSCSD8 QSCTAP1 QSCTAP2 QSCTAP3 QSEBA QSHARE1B QSHARE2B
  QSHELL QSHIP04L QSHIP04S QSHIP08L QSHIP08S QSHIP12L QSHIP12S QSIERRA QSTAIR QSTANDAT S268
  STADAT1 STADAT2 STADAT3 STCQP1 STCQP2 TAME UBH1 VALUES YAO ZECEVIC2
)

fetch_maros() {
  local dest="${dest_root}/maros"
  mkdir -p "${dest}"
  for name in "${maros_names[@]}"; do
    local out="${dest}/${name}.mps"
    [[ -s "${out}" ]] && continue
    echo "maros: ${name}"
    local url="https://raw.githubusercontent.com/optimizers/maros-meszaros-mirror/master"
    if curl -fsSL "${url}/${name}.SIF" -o "${out}.tmp"; then
      mv "${out}.tmp" "${out}"
    else
      rm -f "${out}.tmp"
      echo "maros: ${name} not available, skipped" >&2
    fi
  done
}

fetch_mittelmann() {
  local dest="${dest_root}/mittelmann"
  mkdir -p "${dest}"
  local name
  while read -r name; do
    [[ -z "${name}" || "${name}" == \#* || -s "${dest}/${name}.mps" ]] && continue
    echo "mittelmann: ${name}"
    local file="${name}.mps.bz2"
    [[ "${name}" == L1_sixm250obs ]] && file="${name}.bz2"  # Named without .mps upstream.
    if curl -fsSL "https://plato.asu.edu/ftp/lptestset/${file}" -o "${dest}/${name}.mps.bz2"; then
      bunzip2 -f "${dest}/${name}.mps.bz2"
    else
      rm -f "${dest}/${name}.mps.bz2"
      echo "mittelmann: ${name} not available, skipped" >&2
    fi
  done < "${here}/mittelmann_lp.test"
}

# fetch_emps <netlib directory> <destination> <names...>
fetch_emps() {
  local dir="$1" dest="$2"
  shift 2
  mkdir -p "${dest}"
  local emps="${dest}/.emps"
  if [[ ! -x "${emps}" ]]; then
    curl -fsSL https://www.netlib.org/lp/data/emps.c -o "${dest}/.emps.c"
    cc -O2 -o "${emps}" "${dest}/.emps.c"
  fi
  for name in "$@"; do
    local out="${dest}/${name}.mps"
    [[ -s "${out}" ]] && continue
    echo "netlib: ${name}"
    if curl -fsSL "https://www.netlib.org/lp/${dir}/${name}" | "${emps}" > "${out}.tmp" \
        && [[ -s "${out}.tmp" ]]; then
      mv "${out}.tmp" "${out}"
    else
      rm -f "${out}.tmp"
      echo "netlib: ${name} not available, skipped" >&2
    fi
  done
}

fetch_miplib() {
  local dest="${dest_root}/miplib2017"
  mkdir -p "${dest}"
  local zip="${dest}/.benchmark.zip"
  [[ -s "${zip}" ]] || curl -fSL https://miplib.zib.de/downloads/benchmark.zip -o "${zip}"
  unzip -oq "${zip}" -d "${dest}"
  find "${dest}" -name '*.mps.gz' -exec gunzip -f {} +
  curl -fsSL https://miplib.zib.de/downloads/miplib2017-v31.solu -o "${dest}/miplib2017.solu" || true
}

# fetch_miplib_list <file>: instances listed in <file> (names, '#' comments) into
# bench/instances/<file stem>/, with the known solutions.
fetch_miplib_list() {
  local list="$1"
  local dest="${dest_root}/$(basename "${list%.*}")"
  mkdir -p "${dest}"
  grep -v '^#' "${list}" | while read -r name; do
    [[ -z "${name}" || -s "${dest}/${name}.mps" ]] && continue
    echo "miplib: ${name}"
    if curl -fsSL "https://miplib.zib.de/WebData/instances/${name}.mps.gz" -o "${dest}/${name}.mps.gz"; then
      gunzip -f "${dest}/${name}.mps.gz"
    else
      echo "miplib: ${name} not available, skipped" >&2
    fi
  done
  curl -fsSL https://miplib.zib.de/downloads/miplib2017-v31.solu -o "${dest}/miplib2017.solu" || true
}

if [[ $# -eq 0 ]]; then
  echo "usage: $0 netlib|netlib-infeas|miplib|miplib-list FILE|maros|mittelmann ..." >&2
  exit 2
fi
while [[ $# -gt 0 ]]; do
  set="$1"
  shift
  case "${set}" in
    miplib-list)
      [[ $# -gt 0 ]] || { echo "miplib-list needs a file" >&2; exit 2; }
      fetch_miplib_list "$1"
      shift
      ;;
    netlib) fetch_emps data "${dest_root}/netlib" "${netlib_names[@]}" ;;
    netlib-infeas) fetch_emps infeas "${dest_root}/netlib-infeas" "${netlib_infeas_names[@]}" ;;
    miplib) fetch_miplib ;;
    maros) fetch_maros ;;
    mittelmann) fetch_mittelmann ;;
    *) echo "unknown instance set '${set}'" >&2; exit 2 ;;
  esac
done
