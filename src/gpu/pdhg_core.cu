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
  int* ptr = nullptr;
  int* idx = nullptr;
  double* val = nullptr;
};

DeviceCsr upload_csr_of_transpose(const SparseMatrix& csc) {
  DeviceCsr d;
  d.rows = csc.cols();
  const auto start = csc.col_start();
  std::vector<int> ptr(start.begin(), start.end());
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

  ~Impl() {
    for (void* p : {static_cast<void*>(a_rows.ptr), static_cast<void*>(a_rows.idx),
                    static_cast<void*>(a_rows.val), static_cast<void*>(at_rows.ptr),
                    static_cast<void*>(at_rows.idx), static_cast<void*>(at_rows.val),
                    static_cast<void*>(cost), static_cast<void*>(lower), static_cast<void*>(upper),
                    static_cast<void*>(v), static_cast<void*>(v_prev), static_cast<void*>(v_new),
                    static_cast<void*>(v_bar), static_cast<void*>(y), static_cast<void*>(y_new),
                    static_cast<void*>(ax), static_cast<void*>(aty), static_cast<void*>(v_erg),
                    static_cast<void*>(y_erg)}) {
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
  return std::unique_ptr<PdhgDevice>(new PdhgDevice(std::move(impl)));
}

void PdhgDevice::step(double tau, double sigma, double theta, long long iter) {
  Impl& d = *impl_;
  const int nt = d.n + d.m;
  extrapolate_kernel<<<blocks(nt), kBlock>>>(nt, d.v, d.v_prev, theta, d.v_bar);
  if (d.m > 0) {
    csr_spmv(d.m, d.a_rows.ptr, d.a_rows.idx, d.a_rows.val, d.v_bar, d.ax);
    dual_kernel<<<blocks(d.m), kBlock>>>(d.m, d.n, d.y, d.ax, d.v_bar, sigma, d.y_new);
    csr_spmv(d.n, d.at_rows.ptr, d.at_rows.idx, d.at_rows.val, d.y_new, d.aty);
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
