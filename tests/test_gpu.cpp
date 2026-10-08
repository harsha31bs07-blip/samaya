// GPU kernel tests: validates all GPU kernels against CPU reference implementations.
// When no GPU is available, prints "SKIPPED (no GPU)" and passes.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <memory>
#include <random>
#include <span>
#include <vector>

#include "core/log.hpp"
#include "gpu/device.hpp"
#include "gpu/pdhg_core.hpp"
#include "gpu/spmv.hpp"
#include "gpu/vector_ops.hpp"
#include "lp/lp_solver.hpp"
#include "lp/pdlp.hpp"
#include "lp_generators.hpp"
#include "reference_lp.hpp"
#include "samaya/verify.hpp"
#include "test_framework.hpp"

using samaya::LpProblem;
using samaya::PdlpOptions;
using samaya::PdlpResult;
using samaya::Scaling;
using samaya::SimplexStatus;
using samaya::test::LpFamily;
using samaya::test::ReferenceLp;
using samaya::test::ReferenceResult;

namespace {

// CPU reference implementations for comparison

// Vector operations
void cpu_project_bounds(int n, double* x, const double* lower, const double* upper) {
  for (int i = 0; i < n; ++i) {
    double val = x[i];
    if (val < lower[i]) val = lower[i];
    if (val > upper[i]) val = upper[i];
    x[i] = val;
  }
}

void cpu_axpby(int n, double a, const double* x, double b, double* y) {
  for (int i = 0; i < n; ++i) y[i] = a * x[i] + b * y[i];
}

void cpu_axpy(int n, double a, const double* x, double* y) {
  for (int i = 0; i < n; ++i) y[i] += a * x[i];
}

void cpu_scale(int n, double a, double* x) {
  for (int i = 0; i < n; ++i) x[i] *= a;
}

void cpu_copy(int n, const double* x, double* y) {
  for (int i = 0; i < n; ++i) y[i] = x[i];
}

double cpu_dot(int n, const double* x, const double* y) {
  double sum = 0.0;
  for (int i = 0; i < n; ++i) sum += x[i] * y[i];
  return sum;
}

double cpu_norm2(int n, const double* x) {
  double sum = 0.0;
  for (int i = 0; i < n; ++i) sum += x[i] * x[i];
  return std::sqrt(sum);
}

double cpu_norm2_squared(int n, const double* x) {
  double sum = 0.0;
  for (int i = 0; i < n; ++i) sum += x[i] * x[i];
  return sum;
}

double cpu_norm_inf(int n, const double* x) {
  double max_val = 0.0;
  for (int i = 0; i < n; ++i) {
    double abs_val = std::fabs(x[i]);
    if (abs_val > max_val) max_val = abs_val;
  }
  return max_val;
}

double cpu_sum(int n, const double* x) {
  double sum = 0.0;
  for (int i = 0; i < n; ++i) sum += x[i];
  return sum;
}

// CSR SpMV on CPU
void cpu_csr_spmv(int m, const int* row_ptr, const int* col_idx, const double* values,
                  const double* x, double* y) {
  for (int i = 0; i < m; ++i) {
    double sum = 0.0;
    for (int p = row_ptr[i]; p < row_ptr[i + 1]; ++p) {
      sum += values[p] * x[col_idx[p]];
    }
    y[i] = sum;
  }
}

// Relative error check
bool check_relative(const char* name, const double* gpu, const double* cpu, int n,
                    double tol = 1e-12) {
  double max_rel = 0.0;
  int max_idx = -1;
  for (int i = 0; i < n; ++i) {
    double cpu_val = cpu[i];
    double gpu_val = gpu[i];
    double denom = std::max(std::fabs(cpu_val), 1.0);
    const double rel = std::fabs(gpu_val - cpu_val) / denom;
    if (rel > max_rel) {
      max_rel = rel;
      max_idx = i;
    }
  }
  if (max_rel > tol) {
    std::fprintf(stderr, "  %s: largest relative error %.2e at %d (tolerance %.2e)\n", name,
                 max_rel, max_idx, tol);
    return false;
  }
  return true;
}

// Relative error check of one reduction result (reductions sum in a different order on the
// GPU, so the tolerance allows round-off).
bool check_scalar(const char* name, double gpu, double cpu, double tol = 1e-10) {
  const double rel = std::fabs(gpu - cpu) / std::max(std::fabs(cpu), 1.0);
  if (rel > tol) {
    std::fprintf(stderr, "  %s: gpu %.17g cpu %.17g relative error %.2e\n", name, gpu, cpu, rel);
    return false;
  }
  return true;
}

// Test vector operations
void test_vector_ops() {
  std::printf("Testing vector operations...\n");
  std::mt19937 rng(42);
  std::uniform_real_distribution<double> dist(-10.0, 10.0);

  const int sizes[] = {1, 10, 100, 1000, 10000, 100000};

  for (int n : sizes) {
    std::printf("  Size n=%d\n", n);

    std::vector<double> x(n), y(n), lower(n), upper(n);
    for (int i = 0; i < n; ++i) {
      x[i] = dist(rng);
      y[i] = dist(rng);
      lower[i] = dist(rng);
      upper[i] = lower[i] + std::fabs(dist(rng)) + 1.0;
    }

    // GPU buffers
    double *d_x, *d_y, *d_lower, *d_upper;
    SAMAYA_CUDA_CHECK(cudaMalloc(&d_x, n * sizeof(double)));
    SAMAYA_CUDA_CHECK(cudaMalloc(&d_y, n * sizeof(double)));
    SAMAYA_CUDA_CHECK(cudaMalloc(&d_lower, n * sizeof(double)));
    SAMAYA_CUDA_CHECK(cudaMalloc(&d_upper, n * sizeof(double)));

    SAMAYA_CUDA_CHECK(cudaMemcpy(d_x, x.data(), n * sizeof(double), cudaMemcpyHostToDevice));
    SAMAYA_CUDA_CHECK(cudaMemcpy(d_y, y.data(), n * sizeof(double), cudaMemcpyHostToDevice));
    SAMAYA_CUDA_CHECK(
        cudaMemcpy(d_lower, lower.data(), n * sizeof(double), cudaMemcpyHostToDevice));
    SAMAYA_CUDA_CHECK(
        cudaMemcpy(d_upper, upper.data(), n * sizeof(double), cudaMemcpyHostToDevice));

    // project_bounds
    std::vector<double> x_proj_cpu = x;
    cpu_project_bounds(n, x_proj_cpu.data(), lower.data(), upper.data());
    samaya::gpu::project_bounds(n, d_x, d_lower, d_upper);
    SAMAYA_CUDA_CHECK(cudaMemcpy(x.data(), d_x, n * sizeof(double), cudaMemcpyDeviceToHost));
    CHECK(check_relative("project_bounds", x.data(), x_proj_cpu.data(), n));

    // Restore x
    SAMAYA_CUDA_CHECK(
        cudaMemcpy(d_x, x_proj_cpu.data(), n * sizeof(double), cudaMemcpyHostToDevice));

    // axpby
    double a = dist(rng), b = dist(rng);
    std::vector<double> y_axpby_cpu = y;
    cpu_axpby(n, a, x_proj_cpu.data(), b, y_axpby_cpu.data());
    samaya::gpu::axpby(n, a, d_x, b, d_y);
    SAMAYA_CUDA_CHECK(cudaMemcpy(y.data(), d_y, n * sizeof(double), cudaMemcpyDeviceToHost));
    CHECK(check_relative("axpby", y.data(), y_axpby_cpu.data(), n));

    // Restore y
    SAMAYA_CUDA_CHECK(cudaMemcpy(d_y, y.data(), n * sizeof(double), cudaMemcpyHostToDevice));

    // axpy
    double alpha = dist(rng);
    std::vector<double> y_axpy_cpu = y;
    cpu_axpy(n, alpha, x_proj_cpu.data(), y_axpy_cpu.data());
    samaya::gpu::axpy(n, alpha, d_x, d_y);
    SAMAYA_CUDA_CHECK(cudaMemcpy(y.data(), d_y, n * sizeof(double), cudaMemcpyDeviceToHost));
    CHECK(check_relative("axpy", y.data(), y_axpy_cpu.data(), n));

    // scale
    double scale_factor = dist(rng);
    std::vector<double> x_scale_cpu = x_proj_cpu;
    cpu_scale(n, scale_factor, x_scale_cpu.data());
    samaya::gpu::scale(n, scale_factor, d_x);
    SAMAYA_CUDA_CHECK(cudaMemcpy(x.data(), d_x, n * sizeof(double), cudaMemcpyDeviceToHost));
    CHECK(check_relative("scale", x.data(), x_scale_cpu.data(), n));

    // copy
    std::vector<double> x_copy_cpu(n);
    cpu_copy(n, x_scale_cpu.data(), x_copy_cpu.data());
    samaya::gpu::copy(n, d_x, d_y);  // copy x -> y
    SAMAYA_CUDA_CHECK(cudaMemcpy(y.data(), d_y, n * sizeof(double), cudaMemcpyDeviceToHost));
    CHECK(check_relative("copy", y.data(), x_copy_cpu.data(), n));

    // dot
    double dot_gpu = samaya::gpu::dot(n, d_x, d_y);
    double dot_cpu = cpu_dot(n, x_scale_cpu.data(), x_copy_cpu.data());
    CHECK(check_scalar("dot", dot_gpu, dot_cpu));

    // norm2
    double norm2_gpu = samaya::gpu::norm2(n, d_x);
    double norm2_cpu = cpu_norm2(n, x_scale_cpu.data());
    CHECK(check_scalar("norm2", norm2_gpu, norm2_cpu));

    // norm2_squared
    double norm2_sq_gpu = samaya::gpu::norm2_squared(n, d_x);
    double norm2_sq_cpu = cpu_norm2_squared(n, x_scale_cpu.data());
    CHECK(check_scalar("norm2_squared", norm2_sq_gpu, norm2_sq_cpu));

    // norm_inf
    double norm_inf_gpu = samaya::gpu::norm_inf(n, d_x);
    double norm_inf_cpu = cpu_norm_inf(n, x_scale_cpu.data());
    CHECK(check_scalar("norm_inf", norm_inf_gpu, norm_inf_cpu));

    // sum
    double sum_gpu = samaya::gpu::sum(n, d_x);
    double sum_cpu = cpu_sum(n, x_scale_cpu.data());
    CHECK(check_scalar("sum", sum_gpu, sum_cpu));

    cudaFree(d_x);
    cudaFree(d_y);
    cudaFree(d_lower);
    cudaFree(d_upper);
  }
}

// Generate random sparse matrix in CSR format
struct SparseMatrixCSR {
  int m, n, nnz;
  std::vector<int> row_ptr;
  std::vector<int> col_idx;
  std::vector<double> values;
};

SparseMatrixCSR generate_random_csr(int m, int n, double density, std::mt19937& rng) {
  SparseMatrixCSR A;
  A.m = m;
  A.n = n;
  A.row_ptr.assign(static_cast<std::size_t>(m) + 1, 0);
  std::uniform_real_distribution<double> dist_val(-1.0, 1.0);
  std::uniform_real_distribution<double> coin(0.0, 1.0);
  // One pass: each entry is drawn once, row by row.
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      if (coin(rng) < density) {
        A.col_idx.push_back(j);
        A.values.push_back(dist_val(rng));
      }
    }
    A.row_ptr[static_cast<std::size_t>(i) + 1] = static_cast<int>(A.col_idx.size());
  }
  A.nnz = static_cast<int>(A.col_idx.size());
  return A;
}

// Test SpMV kernels
void test_spmv() {
  std::printf("Testing SpMV kernels...\n");
  std::mt19937 rng(123);

  struct TestCase {
    int m, n;
    double density;
    const char* name;
  };

  TestCase cases[] = {
      {10, 10, 0.5, "small_dense"},
      {100, 50, 0.1, "medium"},
      {1000, 500, 0.01, "large_sparse"},
      {100, 200, 0.05, "wide"},
      {200, 100, 0.05, "tall"},
      // Edge cases
      {1, 10, 1.0, "single_row"},
      {10, 1, 1.0, "single_col"},
      {100, 100, 0.0, "empty"},
  };

  for (const auto& tc : cases) {
    std::printf("  %s: m=%d n=%d density=%.2f\n", tc.name, tc.m, tc.n, tc.density);

    SparseMatrixCSR A = generate_random_csr(tc.m, tc.n, tc.density, rng);

    // Allocate GPU memory
    int *d_row_ptr, *d_col_idx;
    double *d_values, *d_x, *d_y;
    SAMAYA_CUDA_CHECK(cudaMalloc(&d_row_ptr, (tc.m + 1) * sizeof(int)));
    SAMAYA_CUDA_CHECK(cudaMalloc(&d_col_idx, A.nnz * sizeof(int)));
    SAMAYA_CUDA_CHECK(cudaMalloc(&d_values, A.nnz * sizeof(double)));
    SAMAYA_CUDA_CHECK(cudaMalloc(&d_x, tc.n * sizeof(double)));
    SAMAYA_CUDA_CHECK(cudaMalloc(&d_y, tc.m * sizeof(double)));

    SAMAYA_CUDA_CHECK(
        cudaMemcpy(d_row_ptr, A.row_ptr.data(), (tc.m + 1) * sizeof(int), cudaMemcpyHostToDevice));
    SAMAYA_CUDA_CHECK(
        cudaMemcpy(d_col_idx, A.col_idx.data(), A.nnz * sizeof(int), cudaMemcpyHostToDevice));
    SAMAYA_CUDA_CHECK(
        cudaMemcpy(d_values, A.values.data(), A.nnz * sizeof(double), cudaMemcpyHostToDevice));

    // Test y = A * x
    std::vector<double> x(tc.n);
    std::vector<double> y_cpu(tc.m), y_gpu(tc.m);
    std::uniform_real_distribution<double> dist(-10.0, 10.0);
    for (int i = 0; i < tc.n; ++i) x[i] = dist(rng);

    SAMAYA_CUDA_CHECK(cudaMemcpy(d_x, x.data(), tc.n * sizeof(double), cudaMemcpyHostToDevice));

    cpu_csr_spmv(tc.m, A.row_ptr.data(), A.col_idx.data(), A.values.data(), x.data(), y_cpu.data());
    samaya::gpu::csr_spmv(tc.m, d_row_ptr, d_col_idx, d_values, d_x, d_y);
    SAMAYA_CUDA_CHECK(cudaMemcpy(y_gpu.data(), d_y, tc.m * sizeof(double), cudaMemcpyDeviceToHost));
    CHECK(check_relative("csr_spmv", y_gpu.data(), y_cpu.data(), tc.m));

    cudaFree(d_row_ptr);
    cudaFree(d_col_idx);
    cudaFree(d_values);
    cudaFree(d_x);
    cudaFree(d_y);
  }

  // Test with empty rows and one very long row
  std::printf("  Edge case: empty rows and long row\n");
  {
    int m = 100, n = 50;
    SparseMatrixCSR A;
    A.m = m;
    A.n = n;
    A.row_ptr.resize(m + 1, 0);

    // Row 0: empty
    // Row 1: very long (all columns)
    // Row 2..99: empty
    for (int i = 0; i < n; ++i) {
      A.row_ptr[2]++;  // row 1 has n entries
    }
    for (int i = 2; i < m; ++i) {
      A.row_ptr[i + 1] = A.row_ptr[i];
    }
    A.nnz = A.row_ptr[m];

    A.col_idx.resize(A.nnz);
    A.values.resize(A.nnz);
    for (int i = 0; i < n; ++i) {
      A.col_idx[i] = i;
      A.values[i] = 1.0;
    }

    int *d_row_ptr, *d_col_idx;
    double *d_values, *d_x, *d_y;
    SAMAYA_CUDA_CHECK(cudaMalloc(&d_row_ptr, (m + 1) * sizeof(int)));
    SAMAYA_CUDA_CHECK(cudaMalloc(&d_col_idx, A.nnz * sizeof(int)));
    SAMAYA_CUDA_CHECK(cudaMalloc(&d_values, A.nnz * sizeof(double)));
    SAMAYA_CUDA_CHECK(cudaMalloc(&d_x, n * sizeof(double)));
    SAMAYA_CUDA_CHECK(cudaMalloc(&d_y, m * sizeof(double)));

    SAMAYA_CUDA_CHECK(
        cudaMemcpy(d_row_ptr, A.row_ptr.data(), (m + 1) * sizeof(int), cudaMemcpyHostToDevice));
    SAMAYA_CUDA_CHECK(
        cudaMemcpy(d_col_idx, A.col_idx.data(), A.nnz * sizeof(int), cudaMemcpyHostToDevice));
    SAMAYA_CUDA_CHECK(
        cudaMemcpy(d_values, A.values.data(), A.nnz * sizeof(double), cudaMemcpyHostToDevice));

    std::vector<double> x(n, 1.0);
    std::vector<double> y_cpu(m), y_gpu(m);

    SAMAYA_CUDA_CHECK(cudaMemcpy(d_x, x.data(), n * sizeof(double), cudaMemcpyHostToDevice));

    cpu_csr_spmv(m, A.row_ptr.data(), A.col_idx.data(), A.values.data(), x.data(), y_cpu.data());
    samaya::gpu::csr_spmv(m, d_row_ptr, d_col_idx, d_values, d_x, d_y);
    SAMAYA_CUDA_CHECK(cudaMemcpy(y_gpu.data(), d_y, m * sizeof(double), cudaMemcpyDeviceToHost));
    CHECK(check_relative("csr_spmv_empty_long", y_gpu.data(), y_cpu.data(), m));

    cudaFree(d_row_ptr);
    cudaFree(d_col_idx);
    cudaFree(d_values);
    cudaFree(d_x);
    cudaFree(d_y);
  }
}

// Test PDLP GPU solver against CPU reference on small problems
void test_pdlp_gpu_vs_cpu() {
  std::printf("Testing PDLP GPU vs CPU...\n");

  if (!samaya::gpu::gpu_available()) {
    std::printf("  SKIPPED (no GPU)\n");
    return;
  }

  std::mt19937 rng(42);
  PdlpOptions opts;
  opts.tol = 1e-6;
  opts.max_iterations = 10000;
  opts.use_gpu = true;

  samaya::Logger quiet(0);

  int passed = 0, total = 0;

  for (int k = 0; k < 20; ++k) {
    samaya::Model model = samaya::test::random_lp(LpFamily::kFeasible, 6, 6, rng);
    if (!model.Q.empty()) continue;  // Skip QP

    ReferenceResult ref = ReferenceLp(model).solve();
    if (ref.status != ReferenceResult::Status::kOptimal) continue;

    Scaling scaled = samaya::compute_scaling(model.A);
    LpProblem scaled_lp = samaya::make_problem(model, scaled);

    total++;

    // Solve with CPU PDLP
    PdlpOptions cpu_opts = opts;
    cpu_opts.use_gpu = false;
    PdlpResult cpu_result = samaya::solve_pdlp(scaled_lp, cpu_opts, quiet);

    // Solve with GPU PDLP
    PdlpResult gpu_result = samaya::solve_pdlp(scaled_lp, opts, quiet);

    if (cpu_result.status == SimplexStatus::kOptimal &&
        gpu_result.status == SimplexStatus::kOptimal) {
      // Compare solutions (in original space)
      bool ok = true;
      for (int j = 0; j < model.num_cols(); ++j) {
        double rel =
            std::fabs(cpu_result.x[j] - gpu_result.x[j]) /
            std::max(1.0, std::max(std::fabs(cpu_result.x[j]), std::fabs(gpu_result.x[j])));
        if (rel > 1e-4) {
          std::fprintf(stderr, "  Mismatch col %d: cpu=%.6e gpu=%.6e rel=%.2e\n", j,
                       cpu_result.x[j], gpu_result.x[j], rel);
          ok = false;
        }
      }
      if (ok) {
        std::printf("  PASS test %d: both optimal, solutions match\n", k);
        passed++;
      } else {
        std::fprintf(stderr, "  FAIL test %d: solution mismatch\n", k);
      }
    } else if (cpu_result.status != gpu_result.status) {
      std::fprintf(stderr, "  FAIL test %d: status mismatch cpu=%d gpu=%d\n", k,
                   static_cast<int>(cpu_result.status), static_cast<int>(gpu_result.status));
    } else {
      std::printf("  PASS test %d: both status=%d\n", k, static_cast<int>(cpu_result.status));
      passed++;
    }
  }

  std::printf("  PDLP GPU vs CPU: %d/%d passed\n", passed, total);
  CHECK(passed > 0);                     // At least some should pass
  CHECK(passed == total || total == 0);  // All should pass if GPU works
}

}  // namespace

// Main test entry points
TEST(gpu_device_available) {
  bool avail = samaya::gpu::gpu_available();
  std::printf("GPU available: %s\n", avail ? "true" : "false");
  if (!avail) {
    std::printf("  SKIPPED (no GPU)\n");
  }
  // Always passes - just reports availability
}

TEST(gpu_vector_ops) {
  if (!samaya::gpu::gpu_available()) {
    std::printf("  SKIPPED (no GPU)\n");
    return;
  }
  test_vector_ops();
}

TEST(gpu_spmv) {
  if (!samaya::gpu::gpu_available()) {
    std::printf("  SKIPPED (no GPU)\n");
    return;
  }
  test_spmv();
}

TEST(gpu_pdlp_vs_cpu) {
  if (!samaya::gpu::gpu_available()) {
    std::printf("  SKIPPED (no GPU)\n");
    return;
  }
  test_pdlp_gpu_vs_cpu();
}
namespace {

// A small LP in computational form with a random equilibration (row scales rs, column scales cs)
// applied as solve_pdlp applies it, and a device built from it.
struct EquilibratedLp {
  LpProblem lp;
  std::vector<double> rs, cs;
  std::unique_ptr<samaya::gpu::PdhgDevice> device;
};

EquilibratedLp make_device_lp(std::mt19937& rng) {
  EquilibratedLp e;
  const samaya::Model model = samaya::test::random_lp(LpFamily::kFeasible, 20, 20, rng);
  e.lp = samaya::make_problem(model, Scaling::identity(model.num_rows(), model.num_cols()));
  const int m = e.lp.m;
  const int n = e.lp.n;
  std::uniform_real_distribution<double> scale(0.5, 2.0);
  e.rs.resize(m);
  e.cs.resize(n);
  for (double& r : e.rs) r = scale(rng);
  for (double& c : e.cs) c = scale(rng);
  // A_eq = diag(rs) A diag(cs); costs and bounds of x' = x / cs and r' = rs r.
  std::vector<samaya::Triplet> t;
  for (int j = 0; j < n; ++j) {
    for (samaya::NnzIndex p = e.lp.A.col_start()[j]; p < e.lp.A.col_start()[j + 1]; ++p) {
      const int i = e.lp.A.row_index()[p];
      t.push_back({i, j, e.rs[i] * e.lp.A.values()[p] * e.cs[j]});
    }
  }
  const samaya::SparseMatrix a_eq = samaya::SparseMatrix::from_triplets(m, n, std::move(t));
  std::vector<double> cost = e.lp.cost, lower = e.lp.lower, upper = e.lp.upper;
  for (int j = 0; j < n; ++j) {
    cost[j] *= e.cs[j];
    lower[j] /= e.cs[j];
    upper[j] /= e.cs[j];
  }
  for (int i = 0; i < m; ++i) {
    cost[n + i] *= e.rs[i];
    lower[n + i] *= e.rs[i];
    upper[n + i] *= e.rs[i];
  }
  e.device = samaya::gpu::PdhgDevice::create(a_eq, a_eq.transpose(), cost, lower, upper);
  if (e.device) e.device->set_original(e.cs, e.rs, e.lp.cost, e.lp.lower, e.lp.upper);
  return e;
}

// The sums of PdhgDevice::Sums for the scaled iterate (v, y), computed independently on the host
// from the unequilibrated A.
samaya::gpu::PdhgDevice::Sums host_sums(const EquilibratedLp& e, const std::vector<double>& v,
                                        const std::vector<double>& y) {
  const int m = e.lp.m;
  const int n = e.lp.n;
  std::vector<double> value(n + m), y_o(m), ax(m, 0.0), aty(n, 0.0);
  for (int j = 0; j < n; ++j) value[j] = v[j] * e.cs[j];
  for (int i = 0; i < m; ++i) value[n + i] = v[n + i] / e.rs[i];
  for (int i = 0; i < m; ++i) y_o[i] = y[i] * e.rs[i];
  e.lp.A.multiply(std::span<const double>(value.data(), n), ax);
  e.lp.A.multiply_transpose(y_o, aty);
  samaya::gpu::PdhgDevice::Sums s;
  for (int k = 0; k < n + m; ++k) {
    const double g = e.lp.cost[k] + (k < n ? aty[k] : -y_o[k - n]);
    const bool lo = e.lp.lower[k] > -samaya::kInf;
    const bool up = e.lp.upper[k] < samaya::kInf;
    const double lambda = lo && up ? g : lo ? std::max(g, 0.0) : up ? std::min(g, 0.0) : 0.0;
    s.dual_sq += (g - lambda) * (g - lambda);
    if (lambda > 0.0) s.dual_obj += lambda * e.lp.lower[k];
    if (lambda < 0.0) s.dual_obj += lambda * e.lp.upper[k];
    s.primal_obj += e.lp.cost[k] * value[k];
    if (k < n) {
      s.x_sq += value[k] * value[k];
    } else {
      const int i = k - n;
      s.primal_sq += (ax[i] - value[k]) * (ax[i] - value[k]);
      s.r_sq += value[k] * value[k];
      s.y_sq += y_o[i] * y_o[i];
    }
  }
  return s;
}

bool same_sums(const samaya::gpu::PdhgDevice::Sums& a, const samaya::gpu::PdhgDevice::Sums& b,
               double tol) {
  const auto near = [tol](double x, double y) {
    return std::fabs(x - y) <= tol * (1.0 + std::fabs(y));
  };
  return near(a.primal_sq, b.primal_sq) && near(a.dual_sq, b.dual_sq) &&
         near(a.primal_obj, b.primal_obj) && near(a.dual_obj, b.dual_obj) &&
         near(a.x_sq, b.x_sq) && near(a.r_sq, b.r_sq) && near(a.y_sq, b.y_sq);
}

}  // namespace

TEST(gpu_evaluate_matches_host_sums) {
  // The residual terms summed on the device equal a host computation from the unequilibrated
  // matrix, for the current iterate and the average, and are identical when evaluated twice.
  if (!samaya::gpu::gpu_available()) {
    std::printf("  SKIPPED (no GPU)\n");
    return;
  }
  using Iterate = samaya::gpu::PdhgDevice::Iterate;
  std::mt19937 rng(611);
  int compared = 0;
  for (int k = 0; k < 20; ++k) {
    EquilibratedLp e = make_device_lp(rng);
    REQUIRE(e.device != nullptr);
    for (long long it = 0; it < 7; ++it) {
      e.device->step(0.3, 0.3, 1.0, it);
      e.device->advance();
    }
    e.device->step(0.3, 0.3, 1.0, 7);
    std::vector<double> v_new, y_new, v_erg, y_erg;
    e.device->download(nullptr, &v_new, nullptr, &y_new, &v_erg, &y_erg);
    for (const Iterate which : {Iterate::kCurrent, Iterate::kAverage}) {
      const auto device = e.device->evaluate(which);
      const auto again = e.device->evaluate(which);
      const auto host = which == Iterate::kCurrent ? host_sums(e, v_new, y_new)
                                                   : host_sums(e, v_erg, y_erg);
      CHECK(same_sums(device, host, 1e-11));
      CHECK(device.primal_sq == again.primal_sq && device.dual_sq == again.dual_sq &&
            device.dual_obj == again.dual_obj);
      ++compared;
    }
  }
  CHECK_EQ(compared, 40);
}

TEST(gpu_restart_distance_and_restart_to) {
  // restart_distance measures from the last restart point (zero at first); restart_to makes the
  // iterate the current point and the new restart point, so the distance to it becomes zero.
  if (!samaya::gpu::gpu_available()) {
    std::printf("  SKIPPED (no GPU)\n");
    return;
  }
  using Iterate = samaya::gpu::PdhgDevice::Iterate;
  std::mt19937 rng(612);
  EquilibratedLp e = make_device_lp(rng);
  REQUIRE(e.device != nullptr);
  for (long long it = 0; it < 5; ++it) {
    e.device->step(0.3, 0.3, 1.0, it);
    e.device->advance();
  }
  e.device->step(0.3, 0.3, 1.0, 5);
  std::vector<double> v_erg, y_erg;
  e.device->download(nullptr, nullptr, nullptr, nullptr, &v_erg, &y_erg);
  double dv = 0.0, dy = 0.0, host_dv = 0.0, host_dy = 0.0;
  e.device->restart_distance(Iterate::kAverage, dv, dy);
  for (const double x : v_erg) host_dv += x * x;
  for (const double x : y_erg) host_dy += x * x;
  CHECK(std::fabs(dv - host_dv) <= 1e-12 * (1.0 + host_dv));
  CHECK(std::fabs(dy - host_dy) <= 1e-12 * (1.0 + host_dy));
  CHECK(host_dv > 0.0);

  e.device->restart_to(Iterate::kAverage);
  std::vector<double> v, y, v_new, y_new;
  e.device->download(&v, &v_new, &y, &y_new, nullptr, nullptr);
  CHECK(v == v_erg && v_new == v_erg && y == y_erg && y_new == y_erg);
  e.device->restart_distance(Iterate::kCurrent, dv, dy);
  CHECK_EQ(dv, 0.0);
  CHECK_EQ(dy, 0.0);
}

TEST(gpu_pdlp_runs_are_repeatable) {
  // The device sums in a fixed order, so two GPU runs make the same decisions: the same
  // iterations and bit-identical answers.
  if (!samaya::gpu::gpu_available()) {
    std::printf("  SKIPPED (no GPU)\n");
    return;
  }
  std::mt19937 rng(613);
  const samaya::Logger quiet(0);
  PdlpOptions opts;
  opts.tol = 1e-6;
  opts.max_iterations = 20000;
  opts.use_gpu = true;
  for (int k = 0; k < 5; ++k) {
    const samaya::Model model = samaya::test::random_lp(LpFamily::kFeasible, 30, 30, rng);
    const LpProblem lp = samaya::make_problem(model, samaya::compute_scaling(model.A));
    const PdlpResult a = samaya::solve_pdlp(lp, opts, quiet);
    const PdlpResult b = samaya::solve_pdlp(lp, opts, quiet);
    CHECK(a.status == b.status);
    CHECK_EQ(a.iterations, b.iterations);
    CHECK(a.x == b.x && a.y == b.y);
  }
}

TEST(gpu_spmv_tile_widths_match_the_cpu) {
  // Every tile width gives the CPU product on rows that are empty, shorter than the tile, and
  // much longer than it, and the same bits when run twice.
  if (!samaya::gpu::gpu_available()) {
    std::printf("  SKIPPED (no GPU)\n");
    return;
  }
  std::mt19937 rng(614);
  const int m = 3000;
  const int n = 2000;
  std::vector<int> ptr{0}, idx;
  std::vector<double> val;
  std::uniform_int_distribution<int> col(0, n - 1);
  std::uniform_real_distribution<double> value(-1.0, 1.0);
  for (int i = 0; i < m; ++i) {
    const int len = i % 97 == 0 ? 0 : i % 41 == 0 ? 300 : 1 + static_cast<int>(rng() % 12);
    for (int k = 0; k < len; ++k) {
      idx.push_back(col(rng));
      val.push_back(value(rng));
    }
    ptr.push_back(static_cast<int>(idx.size()));
  }
  std::vector<double> x(n);
  for (double& v : x) v = value(rng);
  std::vector<double> y_cpu(m, 0.0);
  for (int i = 0; i < m; ++i) {
    for (int k = ptr[i]; k < ptr[i + 1]; ++k) y_cpu[i] += val[k] * x[idx[k]];
  }
  int *d_ptr = nullptr, *d_idx = nullptr;
  double *d_val = nullptr, *d_x = nullptr, *d_y = nullptr;
  SAMAYA_CUDA_CHECK(cudaMalloc(&d_ptr, ptr.size() * sizeof(int)));
  SAMAYA_CUDA_CHECK(cudaMalloc(&d_idx, idx.size() * sizeof(int)));
  SAMAYA_CUDA_CHECK(cudaMalloc(&d_val, val.size() * sizeof(double)));
  SAMAYA_CUDA_CHECK(cudaMalloc(&d_x, x.size() * sizeof(double)));
  SAMAYA_CUDA_CHECK(cudaMalloc(&d_y, m * sizeof(double)));
  SAMAYA_CUDA_CHECK(
      cudaMemcpy(d_ptr, ptr.data(), ptr.size() * sizeof(int), cudaMemcpyHostToDevice));
  SAMAYA_CUDA_CHECK(
      cudaMemcpy(d_idx, idx.data(), idx.size() * sizeof(int), cudaMemcpyHostToDevice));
  SAMAYA_CUDA_CHECK(
      cudaMemcpy(d_val, val.data(), val.size() * sizeof(double), cudaMemcpyHostToDevice));
  SAMAYA_CUDA_CHECK(cudaMemcpy(d_x, x.data(), x.size() * sizeof(double), cudaMemcpyHostToDevice));
  for (const int t : {2, 4, 8, 16, 32}) {
    std::vector<double> a(m), b(m);
    samaya::gpu::csr_spmv(m, d_ptr, d_idx, d_val, d_x, d_y, 0, t);
    SAMAYA_CUDA_CHECK(cudaMemcpy(a.data(), d_y, m * sizeof(double), cudaMemcpyDeviceToHost));
    samaya::gpu::csr_spmv(m, d_ptr, d_idx, d_val, d_x, d_y, 0, t);
    SAMAYA_CUDA_CHECK(cudaMemcpy(b.data(), d_y, m * sizeof(double), cudaMemcpyDeviceToHost));
    CHECK(check_relative("csr_spmv tile", a.data(), y_cpu.data(), m));
    CHECK(a == b);
  }
  CHECK_EQ(samaya::gpu::spmv_threads_per_row(4117330, 467376), 16);
  CHECK_EQ(samaya::gpu::spmv_threads_per_row(4117330, 1001520), 8);
  for (void* p : {static_cast<void*>(d_ptr), static_cast<void*>(d_idx), static_cast<void*>(d_val),
                  static_cast<void*>(d_x), static_cast<void*>(d_y)}) {
    cudaFree(p);
  }
}
