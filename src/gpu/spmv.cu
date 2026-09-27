// CSR SpMV kernel with one warp per row.
// Double precision throughout.

#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>
#include <cuda_runtime.h>

#include "device.hpp"
#include "spmv.hpp"

namespace cg = cooperative_groups;

namespace samaya {
namespace gpu {

// CSR SpMV: y = A * x
// A is in CSR format: row_ptr (size m+1), col_idx (size nnz), values (size nnz)
// One warp per row.
__global__ void csr_spmv_kernel(const int m, const int* __restrict__ row_ptr,
                                const int* __restrict__ col_idx, const double* __restrict__ values,
                                const double* __restrict__ x, double* __restrict__ y) {
  const int row = blockIdx.x * blockDim.y + threadIdx.y;
  if (row >= m) return;

  const int lane_id = threadIdx.x;

  const int row_start = row_ptr[row];
  const int row_end = row_ptr[row + 1];
  const int row_nnz = row_end - row_start;

  double sum = 0.0;

  // Each thread in the warp processes a subset of the row's nonzeros
  for (int i = lane_id; i < row_nnz; i += 32) {
    const int col = col_idx[row_start + i];
    const double val = values[row_start + i];
    sum += val * x[col];
  }

  // Warp-level reduction
  sum = warp_reduce_sum(sum);

  // First thread in warp writes the result
  if (lane_id == 0) {
    y[row] = sum;
  }
}

// Launch configuration for SpMV kernels.
// Uses 32 threads per warp (x dimension), multiple warps per block (y dimension).
inline dim3 spmv_grid_dim(int rows, int warps_per_block) {
  const int blocks = (rows + warps_per_block - 1) / warps_per_block;
  return dim3(blocks, 1, 1);
}

inline dim3 spmv_block_dim(int warps_per_block) {
  return dim3(32, warps_per_block, 1);
}

// Host launch functions
void csr_spmv(int m, const int* row_ptr, const int* col_idx, const double* values, const double* x,
              double* y, cudaStream_t stream) {
  constexpr int warps_per_block = 8;  // 8 warps = 256 threads per block
  dim3 grid = spmv_grid_dim(m, warps_per_block);
  dim3 block = spmv_block_dim(warps_per_block);
  csr_spmv_kernel<<<grid, block, 0, stream>>>(m, row_ptr, col_idx, values, x, y);
  SAMAYA_CUDA_CHECK(cudaGetLastError());
}

}  // namespace gpu
}  // namespace samaya