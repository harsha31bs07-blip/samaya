#include "lp/pdlp.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <numeric>
#include <random>

#include "core/log.hpp"
#include "linalg/scaling.hpp"
#include "samaya/sparse_matrix.hpp"

#ifdef SAMAYA_HAVE_CUDA
#include "gpu/device.hpp"
#include "gpu/pdhg_core.hpp"
#endif

namespace samaya {

namespace {

// Iterates evaluated (residuals, restarts, detection) once every this many iterations, PDLP's
// termination_evaluation_frequency (Applegate et al. 2021).
constexpr long long kCheckFrequency = 64;

// Power iteration to estimate ||A|| (largest singular value).
// Returns an upper bound on the operator norm of [A -I].
// A is m x n, At is its transpose (n x m).
double estimate_operator_norm(const SparseMatrix& A, const SparseMatrix& At, int max_iters = 20) {
  const int m = A.rows();
  const int n = A.cols();

  // Handle case with no rows (no constraints)
  if (m == 0) return 1.0;

  // We need to estimate the norm of [A -I], which is sqrt(||A||^2 + 1)
  // First estimate ||A||
  std::vector<double> u(static_cast<std::size_t>(m));
  std::vector<double> v(static_cast<std::size_t>(n));
  std::mt19937_64 rng(12345);
  std::uniform_real_distribution<double> dist(-1.0, 1.0);

  // Initialize random vector
  for (int i = 0; i < m; ++i) u[i] = dist(rng);
  double u_norm = 0.0;
  for (double val : u) u_norm += val * val;
  u_norm = std::sqrt(u_norm);
  for (double& val : u) val /= u_norm;

  double sigma = 0.0;
  for (int iter = 0; iter < max_iters; ++iter) {
    // v = A^T u
    std::fill(v.begin(), v.end(), 0.0);
    const auto at_start = At.col_start();
    const auto at_index = At.row_index();
    const auto at_value = At.values();
    for (int j = 0; j < m; ++j) {
      for (NnzIndex p = at_start[j]; p < at_start[j + 1]; ++p) {
        int row = at_index[p];
        v[row] += at_value[p] * u[j];
      }
    }

    double v_norm = 0.0;
    for (double val : v) v_norm += val * val;
    v_norm = std::sqrt(v_norm);
    if (v_norm == 0.0) return 1.0;  // Norm of [A -I] is at least 1

    // u = A v
    std::fill(u.begin(), u.end(), 0.0);
    const auto a_start = A.col_start();
    const auto a_index = A.row_index();
    const auto a_value = A.values();
    for (int j = 0; j < n; ++j) {
      if (v[j] == 0.0) continue;
      for (NnzIndex p = a_start[j]; p < a_start[j + 1]; ++p) {
        int row = a_index[p];
        u[row] += a_value[p] * v[j];
      }
    }

    u_norm = 0.0;
    for (double val : u) u_norm += val * val;
    u_norm = std::sqrt(u_norm);
    if (u_norm == 0.0) return 1.0;

    sigma = u_norm / v_norm;

    // Normalize u for next iteration
    for (double& val : u) val /= u_norm;
  }

  // Norm of [A -I] is sqrt(||A||^2 + 1)
  double L_A = sigma * 1.01;  // Slight overestimate for safety
  return std::sqrt(L_A * L_A + 1.0);
}

// Ruiz equilibration on A only (not the -I part).
// Returns (row_scale, col_scale) where A' = D_r * A * D_c
// row_scale multiplies rows: A' = diag(row_scale) * A
// col_scale multiplies cols: A' = A * diag(col_scale)
void ruiz_equilibration(const SparseMatrix& A, const SparseMatrix& At,
                        std::vector<double>& row_scale, std::vector<double>& col_scale,
                        int max_iters, const Logger& log) {
  const int m = A.rows();
  const int n = A.cols();

  row_scale.assign(static_cast<std::size_t>(m), 1.0);
  col_scale.assign(static_cast<std::size_t>(n), 1.0);

  if (m == 0 || n == 0) return;

  const auto a_start = A.col_start();
  const auto a_index = A.row_index();
  const auto a_value = A.values();

  const auto at_start = At.col_start();
  const auto at_index = At.row_index();
  const auto at_value = At.values();

  std::vector<double> row_norm(m);
  std::vector<double> col_norm(n);

  for (int iter = 0; iter < max_iters; ++iter) {
    // Compute column 2-norms of D_r * A * D_c (current scaling)
    for (int j = 0; j < n; ++j) {
      double sum = 0.0;
      for (NnzIndex p = a_start[j]; p < a_start[j + 1]; ++p) {
        int i = a_index[p];
        double val = a_value[p] * row_scale[i] * col_scale[j];
        sum += val * val;
      }
      col_norm[j] = std::sqrt(sum);
    }

    // Compute row 2-norms of D_r * A * D_c (current scaling)
    for (int i = 0; i < m; ++i) {
      double sum = 0.0;
      for (NnzIndex p = at_start[i]; p < at_start[i + 1]; ++p) {
        int j = at_index[p];
        double val = at_value[p] * row_scale[i] * col_scale[j];
        sum += val * val;
      }
      row_norm[i] = std::sqrt(sum);
    }

    // Update scales: target norm = 1
    // new_scale = old_scale / norm (since norm = |old_scale| * ||unscaled||)
    bool converged = true;
    constexpr double kMaxScale = 1e8;
    constexpr double kMinScale = 1e-8;
    for (int j = 0; j < n; ++j) {
      if (col_norm[j] > 0) {
        double new_scale = col_scale[j] / col_norm[j];
        new_scale = std::max(kMinScale, std::min(kMaxScale, new_scale));
        if (std::fabs(new_scale - col_scale[j]) / std::max(std::fabs(col_scale[j]), 1.0) > 0.1)
          converged = false;
        col_scale[j] = new_scale;
      }
    }
    for (int i = 0; i < m; ++i) {
      if (row_norm[i] > 0) {
        double new_scale = row_scale[i] / row_norm[i];
        new_scale = std::max(kMinScale, std::min(kMaxScale, new_scale));
        if (std::fabs(new_scale - row_scale[i]) / std::max(std::fabs(row_scale[i]), 1.0) > 0.1)
          converged = false;
        row_scale[i] = new_scale;
      }
    }

    if (converged) break;
  }

  log.log(3,
          "pdlp: Ruiz equilibration done, row_scale min/max=%.2e/%.2e col_scale min/max=%.2e/%.2e",
          *std::min_element(row_scale.begin(), row_scale.end()),
          *std::max_element(row_scale.begin(), row_scale.end()),
          *std::min_element(col_scale.begin(), col_scale.end()),
          *std::max_element(col_scale.begin(), col_scale.end()));
}

// Apply equilibration to A: A' = D_r * A * D_c
SparseMatrix apply_equilibration(const SparseMatrix& A, const std::vector<double>& row_scale,
                                 const std::vector<double>& col_scale) {
  std::vector<Triplet> triplets;
  triplets.reserve(A.nnz());

  const auto& a_start = A.col_start();
  const auto& a_index = A.row_index();
  const auto& a_value = A.values();
  const int m = A.rows();
  const int n = A.cols();

  for (int j = 0; j < n; ++j) {
    double col_factor = col_scale[j];
    for (NnzIndex p = a_start[j]; p < a_start[j + 1]; ++p) {
      double scaled_val = a_value[p] * row_scale[a_index[p]] * col_factor;
      triplets.push_back({a_index[p], j, scaled_val});
    }
  }
  return SparseMatrix::from_triplets(m, n, std::move(triplets));
}

// Projection onto box constraints [lower, upper].
// For free variables (infinite bounds), this is identity.
void project_box(const std::vector<double>& lower, const std::vector<double>& upper,
                 std::vector<double>& x) {
  const int nt = static_cast<int>(x.size());
  for (int j = 0; j < nt; ++j) {
    double val = x[j];
    if (lower[j] > -kInf && val < lower[j]) val = lower[j];
    if (upper[j] < kInf && val > upper[j]) val = upper[j];
    x[j] = val;
  }
}

// Check if we should restart based on normalized duality gap.
bool should_restart(double gap, double prev_gap, int iterations_since_restart) {
  if (iterations_since_restart < 200) return false;
  // Restart if gap is significantly increasing in magnitude with same sign (divergence)
  if (gap * prev_gap > 0) {  // Same sign
    if (std::fabs(gap) > 5.0 * std::fabs(prev_gap)) return true;
  }
  // Restart if gap is negative and consistently becoming more negative (strong unbounded signal)
  // Only trigger if gap is significantly negative (< -1e-4) and decreasing
  if (gap < -1e-4 && gap < prev_gap) {
    return true;
  }
  if (std::fabs(gap) < 1e-12 && std::fabs(prev_gap) < 1e-12) return false;
  return false;
}
}  // namespace

PdlpResult solve_pdlp(const LpProblem& lp, const PdlpOptions& options, const Logger& log) {
  const int m = lp.m;
  const int n = lp.n;
  const int nt = n + m;
  PdlpResult result;
  log.log(2, "pdlp: m=%d n=%d nnz=%lld tol=%.2e", m, n, static_cast<long long>(lp.A.nnz()),
          options.tol);

  // --- Preconditioning: Ruiz equilibration on A only ---
  std::vector<double> row_scale, col_scale;
  ruiz_equilibration(lp.A, lp.At, row_scale, col_scale, 10, log);

  // Apply equilibration to A
  SparseMatrix A_eq = apply_equilibration(lp.A, row_scale, col_scale);
  SparseMatrix At_eq = apply_equilibration(lp.At, col_scale, row_scale);

  // Original problem matrix accessors (for certificate verification in original space)
  // Use lp.At (row-compressed, size m+1) for row-wise access, NOT lp.A (column-compressed, size
  // n+1)
  const auto at_start = lp.At.col_start();
  const auto at_index = lp.At.row_index();
  const auto at_value = lp.At.values();
  // Also keep lp.A (column-compressed) for column-wise access in dual infeasibility check
  const auto a_start = lp.A.col_start();
  const auto a_index = lp.A.row_index();
  const auto a_value = lp.A.values();

  // Scaled problem data
  // Primal variables: v' = (x', r') where x' = D_c^{-1} x, r' = D_r r
  // Dual variables: y' = D_r^{-1} y
  // Cost: c'_x = D_c c_x, c'_r = D_r c_r
  // Bounds: l'_x = D_c^{-1} l_x, u'_x = D_c^{-1} u_x, l'_r = D_r l_r, u'_r = D_r u_r
  std::vector<double> cost_eq = lp.cost;
  std::vector<double> lower_eq = lp.lower;
  std::vector<double> upper_eq = lp.upper;

  // Apply scaling to structural variables (0 to n-1): x' = x / col_scale
  for (int j = 0; j < n; ++j) {
    cost_eq[j] *= col_scale[j];
    lower_eq[j] /= col_scale[j];
    upper_eq[j] /= col_scale[j];
  }
  // Row activity variables (n to n+m-1): r' = row_scale * r
  for (int i = 0; i < m; ++i) {
    cost_eq[n + i] *= row_scale[i];
    lower_eq[n + i] *= row_scale[i];
    upper_eq[n + i] *= row_scale[i];
  }

  log.log(3, "pdlp: after scaling, nt=%d, cost_eq size=%d", nt, static_cast<int>(cost_eq.size()));

  // 2. Estimate operator norm for step size: L = ||[A_eq -I]||
  double L = estimate_operator_norm(A_eq, At_eq);
  log.log(3, "pdlp: estimated ||[A -I]|| = %.3e", L);

  // Step sizes: tau * sigma * L^2 <= 1. Use 0.99/L for tau and sigma with theta=1.0 (standard
  // PDHG).
  double tau = 0.99 / L;
  double sigma = 0.99 / L;
  double theta = 1.0;  // Standard PDHG/Chambolle-Pock extrapolation

  log.log(3, "pdlp: tau=%.3e, sigma=%.3e, theta=%.3e, m=%d, n=%d, nt=%d", tau, sigma, theta, m, n,
          nt);

  long long iter_since_restart = 0;

  // --- Initialization ---
  std::vector<double> v(nt, 0.0);       // Primal: (x, r) in SCALED space
  std::vector<double> y(m, 0.0);        // Dual in SCALED space
  std::vector<double> v_prev(nt, 0.0);  // Previous primal

  // Ergodic averages for gap computation (running averages)
  std::vector<double> v_erg(nt, 0.0);  // Ergodic primal average
  std::vector<double> y_erg(m, 0.0);   // Ergodic dual average

  // For unbounded detection: store iterate at restart
  std::vector<double> v_at_restart(nt, 0.0);
  bool has_restart_point = false;
  // Workspace vectors
  std::vector<double> residual(std::max(m, 1));
  std::vector<double> ATy(nt);
  std::vector<double> v_new(nt);
  std::vector<double> v_bar(nt);
  std::vector<double> y_new(m);

  // The iteration arithmetic on the GPU when asked for and available; the algorithm's decisions
  // below stay on the host either way.
#ifdef SAMAYA_HAVE_CUDA
  std::unique_ptr<gpu::PdhgDevice> device;
  if (options.use_gpu) {
    device = gpu::PdhgDevice::create(A_eq, At_eq, cost_eq, lower_eq, upper_eq);
    log.log(1, "pdlp: %s",
            device ? ("iterations on the GPU (" + gpu::device_name() + ")").c_str()
                   : "no usable GPU, iterations on the CPU");
  }
#else
  if (options.use_gpu) log.log(1, "pdlp: built without CUDA, iterations on the CPU");
#endif

  long long max_iter = options.max_iterations;
  if (max_iter <= 0) {
    // Default: proportional to problem size, minimum 10000
    max_iter = std::max<long long>(10000LL, std::min<long long>(10000LL * (m + n) / 10, 1000000LL));
  }

  Timer timer;
  // Largest column norm of A and largest row norm (column of At), for the relative residuals;
  // they do not change during the solve.
  double a_norm_max_col = 0.0;
  double at_norm_max_row = 0.0;
  if (m > 0) {
    for (int j = 0; j < lp.n; ++j) {
      double col_norm = 0.0;
      for (NnzIndex p = lp.A.col_start()[j]; p < lp.A.col_start()[j + 1]; ++p) {
        col_norm += lp.A.values()[p] * lp.A.values()[p];
      }
      a_norm_max_col = std::max(a_norm_max_col, std::sqrt(col_norm));
    }
    for (int i = 0; i < lp.m; ++i) {
      double row_norm = 0.0;
      for (NnzIndex p = lp.At.col_start()[i]; p < lp.At.col_start()[i + 1]; ++p) {
        row_norm += lp.At.values()[p] * lp.At.values()[p];
      }
      at_norm_max_row = std::max(at_norm_max_row, std::sqrt(row_norm));
    }
  }
  double prev_gap = std::numeric_limits<double>::infinity();
  double prev_prev_gap = std::numeric_limits<double>::infinity();
  bool should_restart_next = false;

  log.log(3, "pdlp: after initialization, before main loop");

  // --- Main PDHG loop ---
  // Standard PDHG with extrapolation:
  // v_bar = v + theta * (v - v_prev)
  // y_new = y + sigma * (A_eq * v_bar_x - v_bar_r)
  // v_new = proj(v - tau * (cost + A_eq^T y_new))
  log.log(3, "pdlp: entering main loop, max_iter=%lld, m=%d, n=%d, nt=%d", max_iter, m, n, nt);
  bool done = false;
  for (long long iter = 0; iter < max_iter && !done; ++iter) {
    if (options.time_limit < kInf && timer.seconds() > options.time_limit) {
      result.status = SimplexStatus::kTimeLimit;
      result.iterations = iter;
      break;
    }

    // Adaptive restart at the BEGINNING of the iteration (affects current iteration's
    // extrapolation) Decision is based on gaps from previous two iterations
    if (should_restart_next) {
      log.log(2, "pdlp: RESTART at iter %lld (triggered by prev gap=%.2e, prev_prev_gap=%.2e)",
              iter, prev_gap, prev_prev_gap);
      // Save current iterate as restart point for unbounded detection (only first restart)
      if (!has_restart_point) {
        v_at_restart = v;
        has_restart_point = true;
      }
      // Restart: reset momentum by setting v_prev = v (so v_bar = v, no extrapolation)
      // Also restart from the current iterate (which is the best so far)
      v_prev = v;
#ifdef SAMAYA_HAVE_CUDA
      if (device) device->restart();
#endif
      iter_since_restart = 0;
      should_restart_next = false;
    }

    // Only every kCheckFrequency-th iterate is evaluated (termination, restarts, detection), as
    // in PDLP: the evaluation costs about as much as an iteration.
    const bool evaluate = (iter + 1) % kCheckFrequency == 0 || iter + 1 == max_iter;
#ifdef SAMAYA_HAVE_CUDA
    if (device) {
      device->step(tau, sigma, theta, iter);
      if (!evaluate) {
        device->advance();
        continue;
      }
      device->download(&v, &v_new, nullptr, &y_new, &v_erg, &y_erg);
    } else {
#endif
      // Standard PDHG extrapolation: v_bar = v + theta * (v - v_prev)
      for (int j = 0; j < nt; ++j) {
        v_bar[j] = v[j] + theta * (v[j] - v_prev[j]);
      }

      // Dual update: y_new = y + sigma * (A_eq * v_bar_x - v_bar_r)
      // Residual is A_eq * v_bar_x - v_bar_r (this is [A_eq -I] * v_bar)
      const auto a_eq_start = A_eq.col_start();
      const auto a_eq_index = A_eq.row_index();
      const auto a_eq_value = A_eq.values();

      residual.assign(static_cast<std::size_t>(m), 0.0);
      if (m > 0) {
        for (int j = 0; j < n; ++j) {
          double vj = v_bar[j];
          if (vj == 0.0) continue;
          for (NnzIndex p = a_eq_start[j]; p < a_eq_start[j + 1]; ++p) {
            residual[a_eq_index[p]] += a_eq_value[p] * vj;
          }
        }
        for (int i = 0; i < m; ++i) {
          residual[i] -= v_bar[n + i];  // -I * v_bar_r
        }
      }

      // y_new = y + sigma * residual (dual variables are FREE)
      if (m > 0) {
        for (int i = 0; i < m; ++i) {
          y_new[i] = y[i] + sigma * residual[i];
        }
      }

      // Compute A_eq^T y_new
      const auto at_eq_start = At_eq.col_start();
      const auto at_eq_index = At_eq.row_index();
      const auto at_eq_value = At_eq.values();

      ATy.assign(static_cast<std::size_t>(nt), 0.0);
      if (m > 0) {
        for (int i = 0; i < m; ++i) {
          double yi = y_new[i];
          if (yi == 0.0) continue;
          for (NnzIndex p = at_eq_start[i]; p < at_eq_start[i + 1]; ++p) {
            ATy[at_eq_index[p]] += at_eq_value[p] * yi;
          }
          ATy[n + i] = -yi;  // -I^T * y
        }
      }

      // Primal update: v_new = proj(v - tau * (cost + A_eq^T y_new))
      for (int j = 0; j < nt; ++j) {
        double grad = cost_eq[j] + ATy[j];
        v_new[j] = v[j] - tau * grad;
      }
      project_box(lower_eq, upper_eq, v_new);

      // Update ergodic averages (running averages)
      double alpha = 1.0 / (iter + 1);
      for (int j = 0; j < nt; ++j) {
        v_erg[j] = (1.0 - alpha) * v_erg[j] + alpha * v_new[j];
      }
      if (m > 0) {
        for (int i = 0; i < m; ++i) {
          y_erg[i] = (1.0 - alpha) * y_erg[i] + alpha * y_new[i];
        }
      }

#ifdef SAMAYA_HAVE_CUDA
    }
#endif
    if (!evaluate) {
      v_prev.swap(v);
      v.swap(v_new);
      y.swap(y_new);
      continue;
    }

    // Map current iterate to original space for residual computation
    // x_orig[j] = v_new[j] * col_scale[j] for j in 0..n-1
    // r_orig[i] = v_new[n+i] / row_scale[i] for i in 0..m-1
    // y_orig[i] = y_new[i] * row_scale[i] for i in 0..m-1
    std::vector<double> x_orig(lp.n);
    std::vector<double> r_orig(lp.m);
    std::vector<double> y_orig(lp.m);
    for (int j = 0; j < lp.n; ++j) x_orig[j] = v_new[j] * col_scale[j];
    for (int i = 0; i < lp.m; ++i) r_orig[i] = v_new[lp.n + i] / row_scale[i];
    for (int i = 0; i < lp.m; ++i) y_orig[i] = y_new[i] * row_scale[i];

    // Primal residual in original space: ||A * x_orig - r_orig||
    double primal_res = 0.0;
    if (m > 0) {
      std::vector<double> Ax(lp.m, 0.0);
      for (int j = 0; j < lp.n; ++j) {
        double xj = x_orig[j];
        if (xj == 0.0) continue;
        for (NnzIndex p = a_start[j]; p < a_start[j + 1]; ++p) {
          Ax[a_index[p]] += a_value[p] * xj;
        }
      }
      for (int i = 0; i < lp.m; ++i) {
        double diff = Ax[i] - r_orig[i];
        primal_res += diff * diff;
      }
      primal_res = std::sqrt(primal_res);
    }

    // Dual residual in original space: based on A^T y_orig and original cost/bounds
    double dual_res = 0.0;
    std::vector<double> ATy_orig(lp.n + lp.m, 0.0);
    if (m > 0) {
      // Compute A^T y_orig using original lp.At
      for (int i = 0; i < lp.m; ++i) {
        double yi = y_orig[i];
        if (yi == 0.0) continue;
        for (NnzIndex p = lp.At.col_start()[i]; p < lp.At.col_start()[i + 1]; ++p) {
          ATy_orig[lp.At.row_index()[p]] += lp.At.values()[p] * yi;
        }
      }
      // Row part: -I^T y_orig
      for (int i = 0; i < lp.m; ++i) {
        ATy_orig[lp.n + i] = -y_orig[i];
      }
      // Dual residual (Applegate et al. 2021): the reduced cost g = cost + [A -I]'y of each
      // variable is split into the part its bounds can absorb, lambda (any sign with both bounds
      // finite, >= 0 with only a lower bound, <= 0 with only an upper bound, 0 if free), and the
      // rest, whose norm is the dual residual.
      for (int k = 0; k < lp.n + lp.m; ++k) {
        const double g = lp.cost[k] + ATy_orig[k];
        const bool lo = lp.lower[k] > -kInf;
        const bool up = lp.upper[k] < kInf;
        const double lambda = lo && up ? g : lo ? std::max(g, 0.0) : up ? std::min(g, 0.0) : 0.0;
        dual_res += (g - lambda) * (g - lambda);
      }
      dual_res = std::sqrt(dual_res);
    }

    // Primal/dual objective in original space (scale-invariant)
    double primal_obj = 0.0;
    for (int j = 0; j < lp.n; ++j) primal_obj += lp.cost[j] * x_orig[j];
    for (int i = 0; i < lp.m; ++i) primal_obj += lp.cost[lp.n + i] * r_orig[i];

    // Dual objective: min_{lower <= v <= upper} (cost + A^T y)^T v
    // Note: ATy_orig is zero-initialized above and only populated when m > 0, so this
    // computation is correct (and necessary) even for row-less (pure box-constrained)
    // problems -- it must NOT be gated on m > 0, or dual_obj silently stays 0 forever.
    double dual_obj = 0.0;
    for (int k = 0; k < lp.n + lp.m; ++k) {
      const double g = lp.cost[k] + ATy_orig[k];
      const bool lo = lp.lower[k] > -kInf;
      const bool up = lp.upper[k] < kInf;
      const double lambda = lo && up ? g : lo ? std::max(g, 0.0) : up ? std::min(g, 0.0) : 0.0;
      if (lambda > 0.0) dual_obj += lambda * lp.lower[k];
      if (lambda < 0.0) dual_obj += lambda * lp.upper[k];
    }

    double gap = primal_obj - dual_obj;

    // Compute scaling for relative residuals using ORIGINAL problem data
    // Primal scale: 1 + ||A|| * ||x_orig|| + ||r_orig||
    const double A_norm = a_norm_max_col;  // Constant: computed once before the loop.
    double x_norm = 0.0, r_norm = 0.0;
    for (int j = 0; j < lp.n; ++j) x_norm += x_orig[j] * x_orig[j];
    for (int i = 0; i < lp.m; ++i) r_norm += r_orig[i] * r_orig[i];
    x_norm = std::sqrt(x_norm);
    r_norm = std::sqrt(r_norm);
    double primal_scale = 1.0 + A_norm * x_norm + r_norm;

    // Dual scale: 1 + ||A^T|| * ||y_orig|| + ||c||
    const double At_norm = at_norm_max_row;  // Constant: computed once before the loop.
    double y_norm = 0.0, c_norm = 0.0;
    for (int i = 0; i < lp.m; ++i) y_norm += y_orig[i] * y_orig[i];
    for (int j = 0; j < nt; ++j) c_norm += lp.cost[j] * lp.cost[j];
    y_norm = std::sqrt(y_norm);
    c_norm = std::sqrt(c_norm);
    double dual_scale = 1.0 + At_norm * y_norm + c_norm;

    // Gap scale: 1 + |primal_obj| + |dual_obj|
    double gap_scale = 1.0 + std::fabs(primal_obj) + std::fabs(dual_obj);

    double rel_primal = (m > 0) ? (primal_res / primal_scale) : 0.0;
    double rel_dual = (m > 0) ? (dual_res / dual_scale) : 0.0;
    double rel_gap = gap / gap_scale;

    // Check convergence - early termination when all residuals are below tolerance
    // The gap must be small in magnitude: a negative gap (dual objective above the primal one)
    // means the iterates are not yet consistent, not that they converged.
    if (rel_primal < options.tol && rel_dual < options.tol && std::fabs(rel_gap) < options.tol) {
      log.log(2, "pdlp: converged at iter %lld: primal=%.2e dual=%.2e gap=%.2e", iter, rel_primal,
              rel_dual, rel_gap);
      result.status = SimplexStatus::kOptimal;
      result.iterations = iter + 1;
      result.primal_residual = rel_primal;
      result.dual_residual = rel_dual;
      result.gap = rel_gap;
      v = v_new;
      y = y_new;
      break;
    }

    // --- Infeasibility / Unboundedness detection (M3) ---
    // Check every 100 iterations after the first 1000 iterations
    // (Give ergodic averages time to converge; avoid checking after restarts)
    if (iter > 1000 && m > 0) {
      // --- Check for primal infeasibility (Farkas certificate) ---
      // Unscale the dual certificate to original space: y_orig = y_erg * row_scale
      std::vector<double> y_farkas(m, 0.0);
      double y_max = 0.0;
      for (int i = 0; i < m; ++i) {
        y_farkas[i] = y_erg[i] * row_scale[i];
        y_max = std::max(y_max, std::fabs(y_farkas[i]));
      }
      if (y_max > 1e-6) {  // Non-trivial certificate
        // Compute A^T y_farkas in original space using lp.At (row-compressed)
        std::vector<double> aty_farkas(lp.n + m, 0.0);
        for (int i = 0; i < m; ++i) {
          double yi = y_farkas[i];
          if (yi == 0.0) continue;
          for (NnzIndex p = at_start[i]; p < at_start[i + 1]; ++p) {
            aty_farkas[at_index[p]] += at_value[p] * yi;
          }
        }
        // Row part: -I^T y_farkas
        for (int i = 0; i < m; ++i) {
          aty_farkas[lp.n + i] = -y_farkas[i];
        }
        // Check Farkas certificate against original bounds
        // f(x, r) = (A^T y)^T x - y^T r = aty_farkas[0:n]^T x + aty_farkas[n:n+m]^T r
        // Handle infinite bounds like verify_infeasibility
        auto add_term = [](double coef, double lower, double upper, double& lo, double& hi) {
          if (coef > 0) {
            lo += (lower > -kInf) ? coef * lower : -INFINITY;
            hi += (upper < kInf) ? coef * upper : INFINITY;
          } else if (coef < 0) {
            lo += (upper < kInf) ? coef * upper : -INFINITY;
            hi += (lower > -kInf) ? coef * lower : INFINITY;
          }
        };
        double min_f = 0.0, max_f = 0.0;

        // Compute w_scale for relative tolerance (like verify_infeasibility)
        double w_scale = 0.0;
        for (int j = 0; j < lp.n; ++j) {
          w_scale += std::fabs(aty_farkas[j]);
        }
        for (int i = 0; i < m; ++i) {
          w_scale += std::fabs(aty_farkas[lp.n + i]);
        }
        double tol = 1e-9 * w_scale;  // Relative tolerance

        // Debug: log non-zero coefficients for free variables
        for (int j = 0; j < lp.n; ++j) {
          if (std::fabs(aty_farkas[j]) > tol &&
              (lp.lower[j] <= -kInf / 2 || lp.upper[j] >= kInf / 2)) {
            log.log(2, "pdlp: infeas cert free col j=%d coef=%.2e tol=%.2e lower=%.2e upper=%.2e",
                    j, aty_farkas[j], tol, lp.lower[j], lp.upper[j]);
          }
        }
        for (int i = 0; i < m; ++i) {
          if (std::fabs(aty_farkas[lp.n + i]) > tol &&
              (lp.lower[lp.n + i] <= -kInf / 2 || lp.upper[lp.n + i] >= kInf / 2)) {
            log.log(2, "pdlp: infeas cert free row i=%d coef=%.2e tol=%.2e lower=%.2e upper=%.2e",
                    i, aty_farkas[lp.n + i], tol, lp.lower[lp.n + i], lp.upper[lp.n + i]);
          }
        }

        for (int j = 0; j < lp.n; ++j) {
          double coef = std::fabs(aty_farkas[j]) <= tol ? 0.0 : aty_farkas[j];
          add_term(coef, lp.lower[j], lp.upper[j], min_f, max_f);
        }
        for (int i = 0; i < m; ++i) {
          double coef = aty_farkas[lp.n + i];
          coef = std::fabs(coef) <= tol ? 0.0 : coef;
          add_term(coef, lp.lower[lp.n + i], lp.upper[lp.n + i], min_f, max_f);
        }
        // Use margin matching verify_infeasibility (tol=1e-8)
        // Compute scale similar to verify_infeasibility
        double scale = 0.0;
        for (int j = 0; j < lp.n; ++j) {
          double w = aty_farkas[j] / y_max;
          if (std::fabs(w) > 1e-9) {  // Non-negligible (absolute threshold)
            if (lp.lower[j] > -kInf) scale = std::max(scale, std::fabs(w * lp.lower[j]));
            if (lp.upper[j] < kInf) scale = std::max(scale, std::fabs(w * lp.upper[j]));
          }
        }
        for (int i = 0; i < m; ++i) {
          double coef = -aty_farkas[lp.n + i] / y_max;  // -y_farkas[i]
          if (std::fabs(coef) > 1e-9) {                 // Non-negligible (absolute threshold)
            if (lp.lower[lp.n + i] > -kInf)
              scale = std::max(scale, std::fabs(coef * lp.lower[lp.n + i]));
            if (lp.upper[lp.n + i] < kInf)
              scale = std::max(scale, std::fabs(coef * lp.upper[lp.n + i]));
          }
        }
        if (!std::isfinite(scale)) scale = 0.0;
        double margin = 1e-8 * (1.0 + scale);

        log.log(
            2,
            "pdlp: infeas check iter=%lld y_max=%.2e min_f=%.2e max_f=%.2e margin=%.2e scale=%.2e",
            iter, y_max, min_f, max_f, margin, scale);
        if ((std::isfinite(min_f) && min_f > margin) || (std::isfinite(max_f) && max_f < -margin)) {
          log.log(1, "pdlp: detected primal infeasibility at iter %lld (Farkas certificate)", iter);
          result.status = SimplexStatus::kInfeasible;
          result.iterations = iter + 1;
          result.primal_residual = rel_primal;
          result.dual_residual = rel_dual;
          result.gap = rel_gap;
          // Store Farkas certificate in scaled space (solver.cpp will unscale)
          result.dual_ray.assign(m, 0.0);
          for (int i = 0; i < m; ++i) {
            result.dual_ray[i] = y_erg[i];
          }
          v = v_new;
          y = y_new;
          done = true;
          break;
        }
      }

      // --- Check for dual infeasibility (improving ray) ---
      // Use current iterate as ray candidate (ergodic average v_erg includes
      // early transients that decay as O(1/iter), causing high primal_feas).
      // The current iterate v converges to the ray direction directly.
      // The difference of consecutive iterates: when the LP is unbounded the PDHG iterates
      // diverge along the ray, v_k ~ k d + c, so v_k - v_(k-1) tends to d much faster than
      // v_k / |v_k| does (Applegate, Diaz, Lu, Lubin, "Infeasibility detection with primal-dual
      // hybrid gradient for large-scale linear programming", 2021).
      std::vector<double> ray_candidate(v_new.size());
      for (std::size_t k = 0; k < v_new.size(); ++k) ray_candidate[k] = v_new[k] - v[k];

      double v_max = 0.0;
      for (double vj : ray_candidate) v_max = std::max(v_max, std::fabs(vj));
      if (v_max > 1e-6) {
        // An improving ray is the *direction* v/||v|| converges to, not the raw
        // iterate v itself. Normalize by v_max up front so every downstream
        // check (residual, cost, bounds) operates in scale-invariant ray space.
        // Without this, a bounded/fixed component of v sits at a fixed value
        // (e.g. its bound), which only vanishes relative to the unbounded
        // components as v_max grows -- checking the raw iterate against a
        // tolerance that also grows with v_max can never catch that.
        std::vector<double> ray_dir(ray_candidate.size());
        for (std::size_t k = 0; k < ray_candidate.size(); ++k) {
          ray_dir[k] = ray_candidate[k] / v_max;
        }

        // Check primal feasibility of ray in original space: A * x_orig - r_orig ≈ 0
        std::vector<double> res_check(m, 0.0);
        for (int j = 0; j < lp.n; ++j) {
          double xj = ray_dir[j] * col_scale[j];
          if (xj == 0.0) continue;
          for (NnzIndex p = a_start[j]; p < a_start[j + 1]; ++p) {
            res_check[a_index[p]] += a_value[p] * xj;
          }
        }
        for (int i = 0; i < m; ++i) {
          res_check[i] -=
              ray_dir[lp.n + i] / row_scale[i];  // r' = row_scale * r, so r = r' / row_scale
        }
        double primal_feas = 0.0;
        for (double r : res_check) primal_feas += r * r;
        primal_feas = std::sqrt(primal_feas);

        // Check improving direction: cost^T * (x, r) < 0
        double cost_dir = 0.0;
        for (int j = 0; j < lp.n; ++j) {
          cost_dir += lp.cost[j] * (ray_dir[j] * col_scale[j]);
        }
        for (int i = 0; i < m; ++i) {
          cost_dir += lp.cost[lp.n + i] * (ray_dir[lp.n + i] / row_scale[i]);
        }

        // Check ray bound conditions in original space, using the normalized
        // ray direction. The tolerance is a FIXED constant: ray_dir is
        // already scale-invariant, so it must NOT be scaled up by v_max.
        //
        // Box-constrained variables (lower < upper, both finite): ray component
        // must be 0 in the limit, but PDHG's finite-iterate approximation only gets there with
        // O(1/iter) transient error, so kRayBoundTol has to exceed that error
        // at the iteration budget in use (~2e-5 at 50k iters) rather than sit
        // near machine precision.
        // Fixed variables/rows (lower == upper, i.e. equality constraints): ray
        // component must be exactly 0 in theory, but PDHG's current iterate v
        // converges poorly for equality-constrained components (stagnates at
        // ~1e-4 even at 50k iters). Use a practical tolerance that still
        // rejects non-ray directions but allows the algorithm's limitation.
        // The final verification step (verify_unbounded_ray) will catch any
        // false positives by checking the ray in original space with exact math.
        constexpr double kRayBoundTol = 1e-4;
        constexpr double kRayEqualityTol = 1e-3;
        double ray_bound_tol = kRayBoundTol;
        bool ray_bounds_ok = true;

        // Polish columns first: snap to exact bounds (0 for box/fixed, correct
        // sign for one-sided). This ensures the ray we validate is exactly
        // the ray we will store and return.
        std::vector<double> ray_x(lp.n);
        for (int j = 0; j < lp.n && ray_bounds_ok; ++j) {
          double xj = ray_dir[j] * col_scale[j];
          bool lower_bounded = lp.lower[j] > -kInf / 2;
          bool upper_bounded = lp.upper[j] < kInf / 2;
          bool is_fixed = lower_bounded && upper_bounded &&
                          (lp.upper[j] - lp.lower[j] <= 1e-10 * (1.0 + std::fabs(lp.lower[j])));

          if (lower_bounded && upper_bounded) {
            // Box-constrained (or fixed/equality) variable: ray component must be ~0
            double tol = is_fixed ? kRayEqualityTol : ray_bound_tol;
            if (std::fabs(xj) > tol) {
              ray_bounds_ok = false;
              if (iter % 100 == 0) {
                log.log(3, "pdlp: BOUND FAIL %s col %d: xj=%.6e tol=%.6e bounds=[%.6e, %.6e]",
                        is_fixed ? "fixed" : "box", j, xj, tol, lp.lower[j], lp.upper[j]);
              }
            }
            // Polish to exact 0
            xj = 0.0;
          } else if (lower_bounded) {
            // Only lower bounded: ray component must be >= 0 (up to tol)
            if (xj < -ray_bound_tol) {
              ray_bounds_ok = false;
              if (iter % 100 == 0) {
                log.log(3, "pdlp: BOUND FAIL lb col %d: xj=%.6e tol=%.6e lower=%.6e", j, xj,
                        ray_bound_tol, lp.lower[j]);
              }
            }
            // Polish to >= 0
            xj = std::max(xj, 0.0);
          } else if (upper_bounded) {
            // Only upper bounded: ray component must be <= 0 (up to tol)
            if (xj > ray_bound_tol) {
              ray_bounds_ok = false;
              if (iter % 100 == 0) {
                log.log(3, "pdlp: BOUND FAIL ub col %d: xj=%.6e tol=%.6e upper=%.6e", j, xj,
                        ray_bound_tol, lp.upper[j]);
              }
            }
            // Polish to <= 0
            xj = std::min(xj, 0.0);
          }
          ray_x[j] = xj;
        }

        // Row bound check: verify A * ray_x (row activity) respects row bounds,
        // using the POLISHED ray_x. The verifier checks A * ray_x.
        std::vector<double> av(m, 0.0);
        for (int j = 0; j < lp.n; ++j) {
          double xj = ray_x[j];
          if (xj == 0.0) continue;
          for (NnzIndex p = a_start[j]; p < a_start[j + 1]; ++p) {
            av[a_index[p]] += a_value[p] * xj;
          }
        }
        for (int i = 0; i < m && ray_bounds_ok; ++i) {
          double avi = av[i];
          bool lower_bounded = lp.lower[lp.n + i] > -kInf / 2;
          bool upper_bounded = lp.upper[lp.n + i] < kInf / 2;
          bool is_fixed = lower_bounded && upper_bounded &&
                          (lp.upper[lp.n + i] - lp.lower[lp.n + i] <=
                           1e-10 * (1.0 + std::fabs(lp.lower[lp.n + i])));

          if (lower_bounded && upper_bounded) {
            // Equality constraint: A * ray_x must be ~0
            double tol = is_fixed ? kRayEqualityTol : ray_bound_tol;
            if (std::fabs(avi) > tol) {
              ray_bounds_ok = false;
              if (iter % 100 == 0) {
                log.log(3, "pdlp: BOUND FAIL %s row %d: av=%.6e tol=%.6e bounds=[%.6e, %.6e]",
                        is_fixed ? "fixed" : "box", i, avi, tol, lp.lower[lp.n + i],
                        lp.upper[lp.n + i]);
              }
            }
          } else if (lower_bounded) {
            // A * ray_x >= 0
            if (avi < -ray_bound_tol) {
              ray_bounds_ok = false;
              if (iter % 100 == 0) {
                log.log(3, "pdlp: BOUND FAIL lb row %d: av=%.6e tol=%.6e lower=%.6e", i, avi,
                        ray_bound_tol, lp.lower[lp.n + i]);
              }
            }
          } else if (upper_bounded) {
            // A * ray_x <= 0
            if (avi > ray_bound_tol) {
              ray_bounds_ok = false;
              if (iter % 100 == 0) {
                log.log(3, "pdlp: BOUND FAIL ub row %d: av=%.6e tol=%.6e upper=%.6e", i, avi,
                        ray_bound_tol, lp.upper[lp.n + i]);
              }
            }
          }
        }

        // Debug: log unbounded check
        if (iter % 100 == 0) {
          log.log(2,
                  "pdlp: unbounded check iter=%lld v_max=%.2e primal_feas=%.2e cost_dir=%.2e "
                  "ray_bounds=%d gap=%.2e",
                  iter, v_max, primal_feas, cost_dir, ray_bounds_ok, gap);
        }

        if (primal_feas < 1e-6 && cost_dir < -1e-8 && ray_bounds_ok) {
          log.log(1, "pdlp: detected dual infeasibility (unbounded) at iter %lld (improving ray)",
                  iter);
          result.status = SimplexStatus::kUnbounded;
          result.iterations = iter + 1;
          result.primal_residual = rel_primal;
          result.dual_residual = rel_dual;
          result.gap = rel_gap;
          // ray_x is already polished and validated above; just store it.
          result.primal_ray.assign(lp.n, 0.0);
          for (int j = 0; j < lp.n; ++j) {
            result.primal_ray[j] = ray_x[j];
          }
          v = v_new;
          y = y_new;
          done = true;
          break;
        }
      }
    }

    // Log progress (only at level 3 to avoid excessive output)
    if (log.level() >= 3 && (iter % 100 == 0 || iter < 10)) {
      log.log(3, "pdlp iter %lld: primal=%.2e dual=%.2e gap=%.2e tau=%.2e sigma=%.2e", iter,
              rel_primal, rel_dual, rel_gap, tau, sigma);
    }

    // Adaptive step sizes (like GPU kernel): reduce tau/sigma if residuals increase
    // DISABLED for now - fixed step sizes work better
    // if (iter > 0 && iter % 50 == 0) {
    //   bool reduce = false;
    //   if (rel_primal > prev_primal_res) reduce = true;
    //   if (rel_dual > prev_dual_res) reduce = true;
    //   if (reduce) {
    //     tau *= 0.9;
    //     sigma *= 0.9;
    //     // Don't let step sizes go below a minimum
    //     double min_step = 1e-6 / L;
    //     if (tau < min_step) tau = min_step;
    //     if (sigma < min_step) sigma = min_step;
    //   }
    //   prev_primal_res = rel_primal;
    //   prev_dual_res = rel_dual;
    // }

    // Compute restart decision for NEXT iteration (based on current and previous gaps)
    if (should_restart(gap, prev_gap, iter_since_restart)) {
      should_restart_next = true;
      log.log(2, "pdlp: Will RESTART at next iter (gap=%.2e prev_gap=%.2e iter_since_restart=%lld)",
              gap, prev_gap, iter_since_restart);
    }
    iter_since_restart++;

    // Debug: log y_erg norm for infeasibility detection
    if (iter > 1000 && m > 0) {
      double y_erg_norm = 0.0;
      for (int i = 0; i < m; ++i) y_erg_norm += y_erg[i] * y_erg[i];
      y_erg_norm = std::sqrt(y_erg_norm);
      log.log(2, "pdlp: iter=%lld y_erg_norm=%.2e primal=%.2e dual=%.2e gap=%.2e", iter, y_erg_norm,
              rel_primal, rel_dual, rel_gap);
    }

    // theta is fixed at 1.0 for standard PDHG/Chambolle-Pock

    // Update gap history for next iteration's restart decision
    prev_prev_gap = prev_gap;
    prev_gap = gap;

    v_prev.swap(v);
    v.swap(v_new);
    y.swap(y_new);
#ifdef SAMAYA_HAVE_CUDA
    if (device) device->advance();
#endif
  }
#ifdef SAMAYA_HAVE_CUDA
  // A limit stopped the loop between evaluations: the host copies are behind the device.
  if (device && result.status != SimplexStatus::kOptimal &&
      result.status != SimplexStatus::kInfeasible && result.status != SimplexStatus::kUnbounded) {
    device->download(&v, nullptr, &y, nullptr, nullptr, nullptr);
  }
#endif

  if (result.status == SimplexStatus::kNumericalError) {
    result.status = SimplexStatus::kIterationLimit;
    result.iterations = max_iter;
    log.log(1, "pdlp: iteration limit reached (%lld)", max_iter);
  }

  // --- Map solution back to original scaling ---
  log.log(3, "pdlp: mapping solution back, nt=%d", nt);
  result.x.assign(nt, 0.0);
  log.log(3, "pdlp: after assign, result.x size=%d", static_cast<int>(result.x.size()));

  // Primal variables (0 to n-1): x = x' * col_scale (since x' = x / col_scale)
  for (int j = 0; j < n; ++j) {
    result.x[j] = v[j] * col_scale[j];
  }
  // Row activity variables (n to n+m-1): r = r' / row_scale (since r' = row_scale * r)
  for (int i = 0; i < m; ++i) {
    result.x[n + i] = v[n + i] / row_scale[i];
  }

  // Dual variables: y = y' * row_scale (since y' = y / row_scale)
  if (m > 0) {
    result.y.assign(m, 0.0);
    for (int i = 0; i < m; ++i) {
      result.y[i] = y[i] * row_scale[i];
    }
  } else {
    result.y.assign(0, 0.0);
  }

  result.stats.dual_iterations = result.iterations;
  return result;
}

}  // namespace samaya