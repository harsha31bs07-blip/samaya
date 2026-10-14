#!/usr/bin/env bash
# Print each Simplex::solve profile (including an unscaled cleanup or a crossover).
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "usage: $0 MODEL [extra samaya args]" >&2
  exit 2
fi

here="$(cd "$(dirname "$0")" && pwd)"
samaya="${SAMAYA_BIN:-${here}/../build/release/apps/cli/samaya}"
if [[ ! -x "${samaya}" ]]; then
  echo "release CLI not found: ${samaya} (build release or set SAMAYA_BIN)" >&2
  exit 2
fi

model="$1"
shift
# pipefail propagates solver failures. Keep log level last so it cannot be overridden.
"${samaya}" "${model}" "$@" --log-level 2 | awk '
  /^simplex profile begin/ { printing = 1; found = 1 }
  printing { print }
  /^simplex profile end/ { printing = 0 }
  END { if (!found) { print "no simplex profile (no simplex solve ran)" > "/dev/stderr"; exit 1 } }
'
