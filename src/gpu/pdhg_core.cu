#include <cuda_runtime.h>

#include <algorithm>
#include <climits>
#include <utility>

#include "gpu/device.hpp"
#include "gpu/pdhg_core.hpp"
#include "gpu/spmv.hpp"

namespace samaya {
namespace gpu {

namespace {

constexpr int kBlock = 256;

int blocks(int n) {
  return (n + kBlock - 1) / kBlock;
}

// v_bar = v + theta (v - v_prev)
__global__ void extrapolate_kernel(int nt, const double* v, const double* v_prev, double theta,
                                   double* v_bar) {
  const int k = blockIdx.x * blockDim.x + threadIdx.x;
  if (k < nt) v_bar[k] = v[k] + theta * (v[k] - v_prev[k]);
}

// y_new = y + sigma (A v_bar_x - v_bar_r), with ax = A v_bar_x.
__global__ void dual_kernel(int m, int n, const double* y, const double* ax, const double* v_bar,
                            double sigma, double* y_new) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < m) y_new[i] = y[i] + sigma * (ax[i] - v_bar[n + i]);
}

// v_new = clamp(v - tau (cost + [A -I]' y_new)), with aty = A' y_new for the first n entries.
__global__ void primal_kernel(int n, int m, const double* v, const double* cost, const double* aty,
                              const double* y_new, const double* lower, const double* upper,
                              double tau, double* v_new) {
  const int k = blockIdx.x * blockDim.x + threadIdx.x;
  if (k >= n + m) return;
  const double grad = cost[k] + (k < n ? aty[k] : -y_new[k - n]);
  v_new[k] = fmin(fmax(v[k] - tau * grad, lower[k]), upper[k]);
}

// avg = (1 - alpha) avg + alpha x
__global__ void average_kernel(int len, const double* x, double alpha, double* avg) {
  const int k = blockIdx.x * blockDim.x + threadIdx.x;
  if (k < len) avg[k] = (1.0 - alpha) * avg[k] + alpha * x[k];
}

// Fixed reduction grid: every sum is split over the same blocks in the same order, so the sums
// (and with them the algorithm's decisions) are repeatable from run to run, unlike atomicAdd.
constexpr int kReduceBlocks = 256;
constexpr int kSums = 7;  // PdhgDevice::Sums.

// Sums the K per-thread values of a block in a fixed tree order; thread 0 writes them to
// partial[q * gridDim.x + blockIdx.x].
template <int K>
__device__ void block_sums(const double (&local)[K], double* partial) {
  __shared__ double shared[K][kBlock];
  const int t = threadIdx.x;
  for (int q = 0; q < K; ++q) shared[q][t] = local[q];
  for (int stride = kBlock / 2; stride > 0; stride >>= 1) {
    __syncthreads();
    if (t < stride) {
      for (int q = 0; q < K; ++q) shared[q][t] += shared[q][t + stride];
    }
  }
  if (t == 0) {
    for (int q = 0; q < K; ++q) partial[q * gridDim.x + blockIdx.x] = shared[q][0];
  }
}

// out[q] = sum over the blocks of partial[q * blocks + b], in a fixed order (one block).
template <int K>
__global__ void finish_sums_kernel(int blocks, const double* partial, double* out) {
  double local[K];
  for (int q = 0; q < K; ++q) {
    local[q] = 0.0;
    for (int b = threadIdx.x; b < blocks; b += blockDim.x) local[q] += partial[q * blocks + b];
  }
  __shared__ double shared[K][kBlock];
  const int t = threadIdx.x;
  for (int q = 0; q < K; ++q) shared[q][t] = local[q];
  for (int stride = kBlock / 2; stride > 0; stride >>= 1) {
    __syncthreads();
    if (t < stride) {
      for (int q = 0; q < K; ++q) shared[q][t] += shared[q][t + stride];
    }
  }
  if (t == 0) {
    for (int q = 0; q < K; ++q) out[q] = shared[q][0];
  }
}

// The terms of PdhgDevice::Sums for the iterate (v, y), in the original space: x = cs v_x,
// r = v_r / rs, y_o = rs y. ax = A_eq v_x and aty = A_eq' y (equilibrated), so
// A x = ax / rs and A' y_o = aty / cs. The dual residual and objective follow solve_pdlp's host
// evaluation: the reduced cost g = cost + [A -I]' y_o of each variable is split into the part
// lambda its bounds can absorb and the rest (Applegate et al. 2021).
__global__ void measure_kernel(int n, int m, const double* v, const double* y, const double* ax,
                               const double* aty, const double* cs, const double* rs,
                               const double* cost, const double* lower, const double* upper,
                               double* partial) {
  double local[kSums] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < n + m; k += gridDim.x * blockDim.x) {
    double g = 0.0;
    double value = 0.0;
    if (k < n) {
      value = v[k] * cs[k];
      g = cost[k] + (m > 0 ? aty[k] / cs[k] : 0.0);
      local[4] += value * value;
    } else {
      const int i = k - n;
      value = v[k] / rs[i];
      const double y_o = y[i] * rs[i];
      g = cost[k] - y_o;
      const double e = ax[i] / rs[i] - value;
      local[0] += e * e;
      local[5] += value * value;
      local[6] += y_o * y_o;
    }
    local[2] += cost[k] * value;
    const bool lo = lower[k] > -INFINITY;
    const bool up = upper[k] < INFINITY;
    const double lambda = lo && up ? g : lo ? fmax(g, 0.0) : up ? fmin(g, 0.0) : 0.0;
    local[1] += (g - lambda) * (g - lambda);
    if (lambda > 0.0) local[3] += lambda * lower[k];
    if (lambda < 0.0) local[3] += lambda * upper[k];
  }
  block_sums<kSums>(local, partial);
}

// ||v - v0||^2 (n + m entries) and ||y - y0||^2 (m entries).
__global__ void distance_kernel(int nt, int m, const double* v, const double* v0, const double* y,
                                const double* y0, double* partial) {
  double local[2] = {0.0, 0.0};
  for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < nt; k += gridDim.x * blockDim.x) {
    local[0] += (v[k] - v0[k]) * (v[k] - v0[k]);
    if (k < m) local[1] += (y[k] - y0[k]) * (y[k] - y0[k]);
  }
  block_sums<2>(local, partial);
}

// Uploads a host array; returns the device pointer.
template <typename T>
T* upload(const T* data, std::size_t count) {
  T* d = nullptr;
  SAMAYA_CUDA_CHECK(cudaMalloc(&d, std::max<std::size_t>(count, 1) * sizeof(T)));
  if (count > 0) SAMAYA_CUDA_CHECK(cudaMemcpy(d, data, count * sizeof(T), cudaMemcpyHostToDevice));
  return d;
}

double* zeros(std::size_t count) {
  double* d = nullptr;
  SAMAYA_CUDA_CHECK(cudaMalloc(&d, std::max<std::size_t>(count, 1) * sizeof(double)));
  SAMAYA_CUDA_CHECK(cudaMemset(d, 0, std::max<std::size_t>(count, 1) * sizeof(double)));
  return d;
}

// CSR arrays with 32-bit offsets (the SpMV kernel's format) from a CSC matrix read as the CSR
// of its transpose.
struct DeviceCsr {
  int rows = 0;
  int threads_per_row = 32;  // csr_spmv's tile width for this matrix's average row length.
  int* ptr = nullptr;
  int* idx = nullptr;
  double* val = nullptr;
};

DeviceCsr upload_csr_of_transpose(const SparseMatrix& csc) {
  DeviceCsr d;
  d.rows = csc.cols();
  const auto start = csc.col_start();
  // The device kernels index with 32-bit ints; a matrix too large for them stays on the CPU.
  if (!start.empty() && start.back() > INT_MAX) {
    throw std::runtime_error("matrix too large for the GPU kernels");
  }
  std::vector<int> ptr(start.size());
  for (std::size_t k = 0; k < start.size(); ++k) ptr[k] = static_cast<int>(start[k]);
  d.threads_per_row =
      spmv_threads_per_row(start.empty() ? 0 : static_cast<long long>(start.back()), d.rows);
  d.ptr = upload(ptr.data(), ptr.size());
  d.idx = upload(csc.row_index().data(), csc.row_index().size());
  d.val = upload(csc.values().data(), csc.values().size());
  return d;
}

}  // namespace

struct PdhgDevice::Impl {
  int n = 0;
  int m = 0;
  DeviceCsr a_rows;   // A by rows (the CSC of At read as CSR): A v.
  DeviceCsr at_rows;  // A' by rows (the CSC of A read as CSR): A' y.
  double* cost = nullptr;
  double* lower = nullptr;
  double* upper = nullptr;
  double* v = nullptr;
  double* v_prev = nullptr;
  double* v_new = nullptr;
  double* v_bar = nullptr;
  double* y = nullptr;
  double* y_new = nullptr;
  double* ax = nullptr;
  double* aty = nullptr;
  double* v_erg = nullptr;
  double* y_erg = nullptr;
  // For evaluate(): the equilibration's scales and the cost and bounds before it.
  double* col_scale = nullptr;
  double* row_scale = nullptr;
  double* cost_original = nullptr;
  double* lower_original = nullptr;
  double* upper_original = nullptr;
  double* v_restart = nullptr;
  double* y_restart = nullptr;
  double* partial = nullptr;  // kSums x kReduceBlocks.
  double* sums = nullptr;     // kSums.

  ~Impl() {
    for (void* p : {static_cast<void*>(a_rows.ptr), static_cast<void*>(a_rows.idx),
                    static_cast<void*>(a_rows.val), static_cast<void*>(at_rows.ptr),
                    static_cast<void*>(at_rows.idx), static_cast<void*>(at_rows.val),
                    static_cast<void*>(cost), static_cast<void*>(lower), static_cast<void*>(upper),
                    static_cast<void*>(v), static_cast<void*>(v_prev), static_cast<void*>(v_new),
                    static_cast<void*>(v_bar), static_cast<void*>(y), static_cast<void*>(y_new),
                    static_cast<void*>(ax), static_cast<void*>(aty), static_cast<void*>(v_erg),
                    static_cast<void*>(y_erg), static_cast<void*>(col_scale),
                    static_cast<void*>(row_scale), static_cast<void*>(cost_original),
                    static_cast<void*>(lower_original), static_cast<void*>(upper_original),
                    static_cast<void*>(v_restart), static_cast<void*>(y_restart),
                    static_cast<void*>(partial), static_cast<void*>(sums)}) {
      if (p != nullptr) cudaFree(p);
    }
  }
};

PdhgDevice::PdhgDevice(std::unique_ptr<Impl> impl) : impl_(std::move(impl)) {}
PdhgDevice::~PdhgDevice() = default;

std::unique_ptr<PdhgDevice> PdhgDevice::create(const SparseMatrix& A, const SparseMatrix& At,
                                               const std::vector<double>& cost,
                                               const std::vector<double>& lower,
                                               const std::vector<double>& upper) {
  if (!gpu_available()) return nullptr;
  if (A.nnz() >= static_cast<NnzIndex>(INT_MAX)) return nullptr;  // 32-bit CSR offsets.
  auto impl = std::make_unique<Impl>();
  impl->m = A.rows();
  impl->n = A.cols();
  const auto nt = static_cast<std::size_t>(impl->n + impl->m);
  const auto m = static_cast<std::size_t>(impl->m);
  impl->a_rows = upload_csr_of_transpose(At);
  impl->at_rows = upload_csr_of_transpose(A);
  impl->cost = upload(cost.data(), nt);
  impl->lower = upload(lower.data(), nt);
  impl->upper = upload(upper.data(), nt);
  impl->v = zeros(nt);
  impl->v_prev = zeros(nt);
  impl->v_new = zeros(nt);
  impl->v_bar = zeros(nt);
  impl->v_erg = zeros(nt);
  impl->y = zeros(m);
  impl->y_new = zeros(m);
  impl->y_erg = zeros(m);
  impl->ax = zeros(m);
  impl->aty = zeros(static_cast<std::size_t>(impl->n));
  impl->v_restart = zeros(nt);
  impl->y_restart = zeros(m);
  impl->partial = zeros(static_cast<std::size_t>(kSums * kReduceBlocks));
  impl->sums = zeros(kSums);
  return std::unique_ptr<PdhgDevice>(new PdhgDevice(std::move(impl)));
}

void PdhgDevice::step(double tau, double sigma, double theta, long long iter) {
  Impl& d = *impl_;
  const int nt = d.n + d.m;
  extrapolate_kernel<<<blocks(nt), kBlock>>>(nt, d.v, d.v_prev, theta, d.v_bar);
  if (d.m > 0) {
    csr_spmv(d.m, d.a_rows.ptr, d.a_rows.idx, d.a_rows.val, d.v_bar, d.ax, 0,
             d.a_rows.threads_per_row);
    dual_kernel<<<blocks(d.m), kBlock>>>(d.m, d.n, d.y, d.ax, d.v_bar, sigma, d.y_new);
    csr_spmv(d.n, d.at_rows.ptr, d.at_rows.idx, d.at_rows.val, d.y_new, d.aty, 0,
             d.at_rows.threads_per_row);
  }
  primal_kernel<<<blocks(nt), kBlock>>>(d.n, d.m, d.v, d.cost, d.aty, d.y_new, d.lower, d.upper,
                                        tau, d.v_new);
  const double alpha = 1.0 / static_cast<double>(iter + 1);
  average_kernel<<<blocks(nt), kBlock>>>(nt, d.v_new, alpha, d.v_erg);
  if (d.m > 0) average_kernel<<<blocks(d.m), kBlock>>>(d.m, d.y_new, alpha, d.y_erg);
  SAMAYA_CUDA_CHECK(cudaGetLastError());
}

void PdhgDevice::advance() {
  Impl& d = *impl_;
  std::swap(d.v_prev, d.v);  // v_prev <- v
  std::swap(d.v, d.v_new);   // v <- v_new (v_new now holds scratch)
  std::swap(d.y, d.y_new);
}

void PdhgDevice::restart() {
  Impl& d = *impl_;
  const auto bytes = static_cast<std::size_t>(d.n + d.m) * sizeof(double);
  SAMAYA_CUDA_CHECK(cudaMemcpy(d.v_prev, d.v, bytes, cudaMemcpyDeviceToDevice));
}

void PdhgDevice::set_original(const std::vector<double>& col_scale,
                              const std::vector<double>& row_scale,
                              const std::vector<double>& cost, const std::vector<double>& lower,
                              const std::vector<double>& upper) {
  Impl& d = *impl_;
  const auto nt = static_cast<std::size_t>(d.n + d.m);
  d.col_scale = upload(col_scale.data(), static_cast<std::size_t>(d.n));
  d.row_scale = upload(row_scale.data(), static_cast<std::size_t>(d.m));
  d.cost_original = upload(cost.data(), nt);
  d.lower_original = upload(lower.data(), nt);
  d.upper_original = upload(upper.data(), nt);
}

PdhgDevice::Sums PdhgDevice::evaluate(Iterate which) {
  Impl& d = *impl_;
  const double* v = which == Iterate::kCurrent ? d.v_new : d.v_erg;
  const double* y = which == Iterate::kCurrent ? d.y_new : d.y_erg;
  // ax and aty are scratch between steps: step() recomputes both from v_bar and y_new.
  if (d.m > 0) {
    csr_spmv(d.m, d.a_rows.ptr, d.a_rows.idx, d.a_rows.val, v, d.ax, 0, d.a_rows.threads_per_row);
    csr_spmv(d.n, d.at_rows.ptr, d.at_rows.idx, d.at_rows.val, y, d.aty, 0,
             d.at_rows.threads_per_row);
  }
  measure_kernel<<<kReduceBlocks, kBlock>>>(d.n, d.m, v, y, d.ax, d.aty, d.col_scale, d.row_scale,
                                            d.cost_original, d.lower_original, d.upper_original,
                                            d.partial);
  finish_sums_kernel<kSums><<<1, kBlock>>>(kReduceBlocks, d.partial, d.sums);
  SAMAYA_CUDA_CHECK(cudaGetLastError());
  double out[kSums];
  SAMAYA_CUDA_CHECK(cudaMemcpy(out, d.sums, sizeof(out), cudaMemcpyDeviceToHost));
  Sums s;
  s.primal_sq = out[0];
  s.dual_sq = out[1];
  s.primal_obj = out[2];
  s.dual_obj = out[3];
  s.x_sq = out[4];
  s.r_sq = out[5];
  s.y_sq = out[6];
  return s;
}

void PdhgDevice::restart_distance(Iterate which, double& dv_sq, double& dy_sq) {
  Impl& d = *impl_;
  const double* v = which == Iterate::kCurrent ? d.v_new : d.v_erg;
  const double* y = which == Iterate::kCurrent ? d.y_new : d.y_erg;
  distance_kernel<<<kReduceBlocks, kBlock>>>(d.n + d.m, d.m, v, d.v_restart, y, d.y_restart,
                                             d.partial);
  finish_sums_kernel<2><<<1, kBlock>>>(kReduceBlocks, d.partial, d.sums);
  SAMAYA_CUDA_CHECK(cudaGetLastError());
  double out[2];
  SAMAYA_CUDA_CHECK(cudaMemcpy(out, d.sums, sizeof(out), cudaMemcpyDeviceToHost));
  dv_sq = out[0];
  dy_sq = out[1];
}

void PdhgDevice::restart_to(Iterate which) {
  Impl& d = *impl_;
  const auto vbytes = static_cast<std::size_t>(d.n + d.m) * sizeof(double);
  const auto ybytes = static_cast<std::size_t>(d.m) * sizeof(double);
  const double* v = which == Iterate::kCurrent ? d.v_new : d.v_erg;
  const double* y = which == Iterate::kCurrent ? d.y_new : d.y_erg;
  const auto copy = [](double* to, const double* from, std::size_t bytes) {
    if (bytes > 0 && to != from) {
      SAMAYA_CUDA_CHECK(cudaMemcpy(to, from, bytes, cudaMemcpyDeviceToDevice));
    }
  };
  // The restart point first: v_new or y_new may be the source.
  copy(d.v_restart, v, vbytes);
  copy(d.y_restart, y, ybytes);
  copy(d.v, d.v_restart, vbytes);
  copy(d.v_new, d.v_restart, vbytes);
  copy(d.y, d.y_restart, ybytes);
  copy(d.y_new, d.y_restart, ybytes);
}

void PdhgDevice::download(std::vector<double>* v, std::vector<double>* v_new,
                          std::vector<double>* y, std::vector<double>* y_new,
                          std::vector<double>* v_erg, std::vector<double>* y_erg) const {
  const Impl& d = *impl_;
  const auto nt = static_cast<std::size_t>(d.n + d.m);
  const auto m = static_cast<std::size_t>(d.m);
  const auto get = [](std::vector<double>* out, const double* src, std::size_t count) {
    if (out == nullptr) return;
    out->resize(count);
    if (count > 0) {
      SAMAYA_CUDA_CHECK(
          cudaMemcpy(out->data(), src, count * sizeof(double), cudaMemcpyDeviceToHost));
    }
  };
  get(v, d.v, nt);
  get(v_new, d.v_new, nt);
  get(y, d.y, m);
  get(y_new, d.y_new, m);
  get(v_erg, d.v_erg, nt);
  get(y_erg, d.y_erg, m);
}

}  // namespace gpu
}  // namespace samaya
