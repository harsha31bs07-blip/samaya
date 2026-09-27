#!/usr/bin/env bash
# samaya with PDLP on the CPU, from the CUDA build (the same binary as bench/samaya_pdlp_gpu.sh, so
# the two differ only in where the iterations run).
exec "$(dirname "$0")/../build/cuda/apps/cli/samaya" --lp-method pdlp "$@"
