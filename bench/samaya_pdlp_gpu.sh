#!/usr/bin/env bash
# samaya with PDLP, iterations on the GPU (see bench/samaya_pdlp.sh for the CPU twin).
exec "$(dirname "$0")/../build/cuda/apps/cli/samaya" --lp-method pdlp --gpu "$@"
