#pragma once

// SpMV kernels for GPU.
// CSR format with one warp per row.

#include "device.hpp"

namespace samaya {
namespace gpu {

// Launches y = A * x where A is m x n in CSR format.
// row_ptr: size m+1, col_idx: size nnz, values: size nnz
// x: size n, y: size m (output)
void csr_spmv(int m, const int* row_ptr, const int* col_idx, const double* values, const double* x,
              double* y, cudaStream_t stream = 0);

}  // namespace gpu
}  // namespace samaya