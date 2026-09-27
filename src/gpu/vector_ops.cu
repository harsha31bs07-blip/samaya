// Fused vector operations for PDLP on GPU.
// Double precision throughout.

#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>
#include <cuda_runtime.h>

#include "device.hpp"
#include "vector_ops.hpp"

namespace cg = cooperative_groups;

namespace samaya {
namespace gpu {

// Block size for vector kernels
constexpr int VECTOR_BLOCK_SIZE = 256;

// Projection onto bounds: x = clamp(x, lower, upper)
__global__ void project_bounds_kernel(int n, double* x, const double* lower, const double* upper) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= n) return;
  double val = x[idx];
  double lo = lower[idx];
  double up = upper[idx];
  if (val < lo) val = lo;
  if (val > up) val = up;
  x[idx] = val;
}

void project_bounds(int n, double* x, const double* lower, const double* upper,
                    cudaStream_t stream) {
  int blocks = (n + VECTOR_BLOCK_SIZE - 1) / VECTOR_BLOCK_SIZE;
  project_bounds_kernel<<<blocks, VECTOR_BLOCK_SIZE, 0, stream>>>(n, x, lower, upper);
  SAMAYA_CUDA_CHECK(cudaGetLastError());
}

// axpby: y = a*x + b*y
__global__ void axpby_kernel(int n, double a, const double* x, double b, double* y) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= n) return;
  y[idx] = a * x[idx] + b * y[idx];
}

void axpby(int n, double a, const double* x, double b, double* y, cudaStream_t stream) {
  int blocks = (n + VECTOR_BLOCK_SIZE - 1) / VECTOR_BLOCK_SIZE;
  axpby_kernel<<<blocks, VECTOR_BLOCK_SIZE, 0, stream>>>(n, a, x, b, y);
  SAMAYA_CUDA_CHECK(cudaGetLastError());
}

// axpy: y += a*x
__global__ void axpy_kernel(int n, double a, const double* x, double* y) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= n) return;
  y[idx] += a * x[idx];
}

void axpy(int n, double a, const double* x, double* y, cudaStream_t stream) {
  int blocks = (n + VECTOR_BLOCK_SIZE - 1) / VECTOR_BLOCK_SIZE;
  axpy_kernel<<<blocks, VECTOR_BLOCK_SIZE, 0, stream>>>(n, a, x, y);
  SAMAYA_CUDA_CHECK(cudaGetLastError());
}

// Scale: x *= a
__global__ void scale_kernel(int n, double a, double* x) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= n) return;
  x[idx] *= a;
}

void scale(int n, double a, double* x, cudaStream_t stream) {
  int blocks = (n + VECTOR_BLOCK_SIZE - 1) / VECTOR_BLOCK_SIZE;
  scale_kernel<<<blocks, VECTOR_BLOCK_SIZE, 0, stream>>>(n, a, x);
  SAMAYA_CUDA_CHECK(cudaGetLastError());
}

// Copy: y = x
__global__ void copy_kernel(int n, const double* x, double* y) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= n) return;
  y[idx] = x[idx];
}

void copy(int n, const double* x, double* y, cudaStream_t stream) {
  int blocks = (n + VECTOR_BLOCK_SIZE - 1) / VECTOR_BLOCK_SIZE;
  copy_kernel<<<blocks, VECTOR_BLOCK_SIZE, 0, stream>>>(n, x, y);
  SAMAYA_CUDA_CHECK(cudaGetLastError());
}

// Reduction kernels for dot, norm2, norm_inf, sum
// Using a two-stage reduction: block reduction then atomic add to global

// Block reduction for sum
__inline__ __device__ double block_reduce_sum(double val) {
  __shared__ double shared[32];
  cg::thread_block_tile<32> warp = cg::tiled_partition<32>(cg::this_thread_block());
  double sum = warp_reduce_sum(val);
  if (warp.thread_rank() == 0) {
    shared[warp.meta_group_rank()] = sum;
  }
  __syncthreads();
  // First warp reduces the partial sums
  if (warp.meta_group_rank() == 0) {
    double v = warp.thread_rank() < (blockDim.x / 32) ? shared[warp.thread_rank()] : 0.0;
    sum = warp_reduce_sum(v);
  }
  return warp.thread_rank() == 0 && warp.meta_group_rank() == 0 ? sum : 0.0;
}

// Block reduction for max
__inline__ __device__ double block_reduce_max(double val) {
  __shared__ double shared[32];
  cg::thread_block_tile<32> warp = cg::tiled_partition<32>(cg::this_thread_block());
  double max_val = val;
  for (int offset = 16; offset > 0; offset >>= 1) {
    double tmp = warp.shfl_down(max_val, offset);
    if (tmp > max_val) max_val = tmp;
  }
  if (warp.thread_rank() == 0) {
    shared[warp.meta_group_rank()] = max_val;
  }
  __syncthreads();
  if (warp.meta_group_rank() == 0) {
    double v = warp.thread_rank() < (blockDim.x / 32) ? shared[warp.thread_rank()] : -1.0;
    for (int offset = 16; offset > 0; offset >>= 1) {
      double tmp = warp.shfl_down(v, offset);
      if (tmp > v) v = tmp;
    }
    max_val = v;
  }
  return warp.thread_rank() == 0 && warp.meta_group_rank() == 0 ? max_val : 0.0;
}

// Dot product kernel
__global__ void dot_kernel(int n, const double* x, const double* y, double* result) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  double sum = 0.0;
  for (int i = idx; i < n; i += gridDim.x * blockDim.x) {
    sum += x[i] * y[i];
  }
  sum = block_reduce_sum(sum);
  if (threadIdx.x == 0) {
    atomicAdd(result, sum);
  }
}

double dot(int n, const double* x, const double* y, cudaStream_t stream) {
  double* d_result;
  SAMAYA_CUDA_CHECK(cudaMalloc(&d_result, sizeof(double)));
  SAMAYA_CUDA_CHECK(cudaMemsetAsync(d_result, 0, sizeof(double), stream));
  int blocks = min(65535, (n + VECTOR_BLOCK_SIZE - 1) / VECTOR_BLOCK_SIZE);
  dot_kernel<<<blocks, VECTOR_BLOCK_SIZE, 0, stream>>>(n, x, y, d_result);
  SAMAYA_CUDA_CHECK(cudaGetLastError());
  double result;
  SAMAYA_CUDA_CHECK(
      cudaMemcpyAsync(&result, d_result, sizeof(double), cudaMemcpyDeviceToHost, stream));
  SAMAYA_CUDA_CHECK(cudaStreamSynchronize(stream));
  SAMAYA_CUDA_CHECK(cudaFree(d_result));
  return result;
}

// Norm squared kernel
__global__ void norm2_squared_kernel(int n, const double* x, double* result) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  double sum = 0.0;
  for (int i = idx; i < n; i += gridDim.x * blockDim.x) {
    double val = x[i];
    sum += val * val;
  }
  sum = block_reduce_sum(sum);
  if (threadIdx.x == 0) {
    atomicAdd(result, sum);
  }
}

double norm2_squared(int n, const double* x, cudaStream_t stream) {
  double* d_result;
  SAMAYA_CUDA_CHECK(cudaMalloc(&d_result, sizeof(double)));
  SAMAYA_CUDA_CHECK(cudaMemsetAsync(d_result, 0, sizeof(double), stream));
  int blocks = min(65535, (n + VECTOR_BLOCK_SIZE - 1) / VECTOR_BLOCK_SIZE);
  norm2_squared_kernel<<<blocks, VECTOR_BLOCK_SIZE, 0, stream>>>(n, x, d_result);
  SAMAYA_CUDA_CHECK(cudaGetLastError());
  double result;
  SAMAYA_CUDA_CHECK(
      cudaMemcpyAsync(&result, d_result, sizeof(double), cudaMemcpyDeviceToHost, stream));
  SAMAYA_CUDA_CHECK(cudaStreamSynchronize(stream));
  SAMAYA_CUDA_CHECK(cudaFree(d_result));
  return result;
}

double norm2(int n, const double* x, cudaStream_t stream) {
  return std::sqrt(norm2_squared(n, x, stream));
}

// Inf norm kernel
__global__ void norm_inf_kernel(int n, const double* x, double* result) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  double max_val = 0.0;
  for (int i = idx; i < n; i += gridDim.x * blockDim.x) {
    double val = fabs(x[i]);
    if (val > max_val) max_val = val;
  }
  max_val = block_reduce_max(max_val);
  if (threadIdx.x == 0) {
    atomicMax(reinterpret_cast<unsigned long long*>(result), __double_as_longlong(max_val));
  }
}

double norm_inf(int n, const double* x, cudaStream_t stream) {
  double* d_result;
  SAMAYA_CUDA_CHECK(cudaMalloc(&d_result, sizeof(double)));
  SAMAYA_CUDA_CHECK(cudaMemsetAsync(d_result, 0, sizeof(double), stream));
  int blocks = min(65535, (n + VECTOR_BLOCK_SIZE - 1) / VECTOR_BLOCK_SIZE);
  norm_inf_kernel<<<blocks, VECTOR_BLOCK_SIZE, 0, stream>>>(n, x, d_result);
  SAMAYA_CUDA_CHECK(cudaGetLastError());
  double result;
  SAMAYA_CUDA_CHECK(
      cudaMemcpyAsync(&result, d_result, sizeof(double), cudaMemcpyDeviceToHost, stream));
  SAMAYA_CUDA_CHECK(cudaStreamSynchronize(stream));
  SAMAYA_CUDA_CHECK(cudaFree(d_result));
  return result;
}

// Sum kernel
__global__ void sum_kernel(int n, const double* x, double* result) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  double sum = 0.0;
  for (int i = idx; i < n; i += gridDim.x * blockDim.x) {
    sum += x[i];
  }
  sum = block_reduce_sum(sum);
  if (threadIdx.x == 0) {
    atomicAdd(result, sum);
  }
}

double sum(int n, const double* x, cudaStream_t stream) {
  double* d_result;
  SAMAYA_CUDA_CHECK(cudaMalloc(&d_result, sizeof(double)));
  SAMAYA_CUDA_CHECK(cudaMemsetAsync(d_result, 0, sizeof(double), stream));
  int blocks = min(65535, (n + VECTOR_BLOCK_SIZE - 1) / VECTOR_BLOCK_SIZE);
  sum_kernel<<<blocks, VECTOR_BLOCK_SIZE, 0, stream>>>(n, x, d_result);
  SAMAYA_CUDA_CHECK(cudaGetLastError());
  double result;
  SAMAYA_CUDA_CHECK(
      cudaMemcpyAsync(&result, d_result, sizeof(double), cudaMemcpyDeviceToHost, stream));
  SAMAYA_CUDA_CHECK(cudaStreamSynchronize(stream));
  SAMAYA_CUDA_CHECK(cudaFree(d_result));
  return result;
}

// Primal residual norm: ||Ax - r||
__global__ void primal_residual_kernel(int m, const double* Ax, const double* r, double* result) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  double sum = 0.0;
  for (int i = idx; i < m; i += gridDim.x * blockDim.x) {
    double diff = Ax[i] - r[i];
    sum += diff * diff;
  }
  sum = block_reduce_sum(sum);
  if (threadIdx.x == 0) {
    atomicAdd(result, sum);
  }
}

double primal_residual_norm(int m, const double* Ax, const double* r, cudaStream_t stream) {
  double* d_result;
  SAMAYA_CUDA_CHECK(cudaMalloc(&d_result, sizeof(double)));
  SAMAYA_CUDA_CHECK(cudaMemsetAsync(d_result, 0, sizeof(double), stream));
  int blocks = min(65535, (m + VECTOR_BLOCK_SIZE - 1) / VECTOR_BLOCK_SIZE);
  primal_residual_kernel<<<blocks, VECTOR_BLOCK_SIZE, 0, stream>>>(m, Ax, r, d_result);
  SAMAYA_CUDA_CHECK(cudaGetLastError());
  double result;
  SAMAYA_CUDA_CHECK(
      cudaMemcpyAsync(&result, d_result, sizeof(double), cudaMemcpyDeviceToHost, stream));
  SAMAYA_CUDA_CHECK(cudaStreamSynchronize(stream));
  SAMAYA_CUDA_CHECK(cudaFree(d_result));
  return std::sqrt(result);
}

// Dual residual norm: ||A^T y + d - c||
__global__ void dual_residual_kernel(int n, const double* ATy, const double* d, const double* c,
                                     double* result) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  double sum = 0.0;
  for (int i = idx; i < n; i += gridDim.x * blockDim.x) {
    double diff = ATy[i] + d[i] - c[i];
    sum += diff * diff;
  }
  sum = block_reduce_sum(sum);
  if (threadIdx.x == 0) {
    atomicAdd(result, sum);
  }
}

double dual_residual_norm(int n, const double* ATy, const double* d, const double* c,
                          cudaStream_t stream) {
  double* d_result;
  SAMAYA_CUDA_CHECK(cudaMalloc(&d_result, sizeof(double)));
  SAMAYA_CUDA_CHECK(cudaMemsetAsync(d_result, 0, sizeof(double), stream));
  int blocks = min(65535, (n + VECTOR_BLOCK_SIZE - 1) / VECTOR_BLOCK_SIZE);
  dual_residual_kernel<<<blocks, VECTOR_BLOCK_SIZE, 0, stream>>>(n, ATy, d, c, d_result);
  SAMAYA_CUDA_CHECK(cudaGetLastError());
  double result;
  SAMAYA_CUDA_CHECK(
      cudaMemcpyAsync(&result, d_result, sizeof(double), cudaMemcpyDeviceToHost, stream));
  SAMAYA_CUDA_CHECK(cudaStreamSynchronize(stream));
  SAMAYA_CUDA_CHECK(cudaFree(d_result));
  return std::sqrt(result);
}

// Gap numerator: |c^T x - b^T y|
// We compute c^T x and b^T y separately
double gap_numerator(int n, int m, const double* c, const double* x, const double* b,
                     const double* y, cudaStream_t stream) {
  double cx = dot(n, c, x, stream);
  double by = dot(m, b, y, stream);
  return fabs(cx - by);
}

}  // namespace gpu
}  // namespace samaya