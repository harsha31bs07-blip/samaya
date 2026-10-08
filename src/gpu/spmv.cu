// CSR SpMV kernels: T threads per row (T = 2 ... 32), one warp per row when T = 32.
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

// T threads per row (a tile of a warp), for matrices with short rows: with one warp per row a
// row of 4 nonzeros keeps 28 of the 32 lanes idle. The tile sums its products in a fixed
// shuffle order, so the result does not depend on timing.
template <int T>
__global__ void csr_spmv_tile_kernel(const int m, const int* __restrict__ row_ptr,
                                     const int* __restrict__ col_idx,
                                     const double* __restrict__ values,
                                     const double* __restrict__ x, double* __restrict__ y) {
  const cg::thread_block_tile<T> tile = cg::tiled_partition<T>(cg::this_thread_block());
  const int row = static_cast<int>((static_cast<long long>(blockIdx.x) * blockDim.x +
                                    threadIdx.x) / T);
  const int lane = static_cast<int>(tile.thread_rank());
  double sum = 0.0;
  // No early return: every lane of the tile takes part in the shuffles below.
  if (row < m) {
    const int end = row_ptr[row + 1];
    for (int k = row_ptr[row] + lane; k < end; k += T) sum += values[k] * x[col_idx[k]];
  }
  for (int offset = T / 2; offset > 0; offset /= 2) sum += tile.shfl_down(sum, offset);
  if (row < m && lane == 0) y[row] = sum;
}

int spmv_threads_per_row(long long nnz, int rows) {
  // The smallest tile at least as long as the average row (CUSP's vector-CSR rule), from 2 to 32.
  const double mean = rows > 0 ? static_cast<double>(nnz) / rows : 0.0;
  int t = 2;
  while (t < 32 && t < mean) t *= 2;
  return t;
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
              double* y, cudaStream_t stream, int threads_per_row) {
  if (threads_per_row < 32 && m > 0) {
    constexpr int kThreads = 256;
    const long long threads = static_cast<long long>(m) * threads_per_row;
    const int blocks = static_cast<int>((threads + kThreads - 1) / kThreads);
    switch (threads_per_row) {
      case 2:
        csr_spmv_tile_kernel<2><<<blocks, kThreads, 0, stream>>>(m, row_ptr, col_idx, values, x, y);
        break;
      case 4:
        csr_spmv_tile_kernel<4><<<blocks, kThreads, 0, stream>>>(m, row_ptr, col_idx, values, x, y);
        break;
      case 8:
        csr_spmv_tile_kernel<8><<<blocks, kThreads, 0, stream>>>(m, row_ptr, col_idx, values, x, y);
        break;
      default:
        csr_spmv_tile_kernel<16><<<blocks, kThreads, 0, stream>>>(m, row_ptr, col_idx, values, x,
                                                                  y);
        break;
    }
    SAMAYA_CUDA_CHECK(cudaGetLastError());
    return;
  }
  constexpr int warps_per_block = 8;  // 8 warps = 256 threads per block
  dim3 grid = spmv_grid_dim(m, warps_per_block);
  dim3 block = spmv_block_dim(warps_per_block);
  csr_spmv_kernel<<<grid, block, 0, stream>>>(m, row_ptr, col_idx, values, x, y);
  SAMAYA_CUDA_CHECK(cudaGetLastError());
}

}  // namespace gpu
}  // namespace samaya