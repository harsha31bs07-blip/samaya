#!/usr/bin/env bash
# samaya with the interior-point LP method and crossover, for bench/compare.sh netlib-crossover.
exec "$(dirname "$0")/../build/release/apps/cli/samaya" --lp-method barrier --crossover "$@"
