#pragma once

// SpMV kernels for GPU.
// CSR format with one warp per row.

#include "device.hpp"

namespace samaya {
namespace gpu {

// Launches y = A * x where A is m x n in CSR format.
// row_ptr: size m+1, col_idx: size nnz, values: size nnz
// x: size n, y: size m (output)
// threads_per_row: 2, 4, 8, 16 or 32 threads share a row (spmv_threads_per_row picks it from the
// average row length); 32 is one warp per row.
void csr_spmv(int m, const int* row_ptr, const int* col_idx, const double* values, const double* x,
              double* y, cudaStream_t stream = 0, int threads_per_row = 32);

// Threads per row for csr_spmv: the smallest power of two from 2 to 32 not below the average
// number of nonzeros per row.
int spmv_threads_per_row(long long nnz, int rows);

}  // namespace gpu
}  // namespace samaya