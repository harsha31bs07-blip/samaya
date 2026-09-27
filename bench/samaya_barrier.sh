#!/usr/bin/env bash
# samaya with the interior-point LP method, for bench/compare.sh netlib-barrier.
exec "$(dirname "$0")/../build/release/apps/cli/samaya" --lp-method barrier "$@"
