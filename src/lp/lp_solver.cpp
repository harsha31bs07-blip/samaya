#include "lp/lp_solver.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <memory>

#include "core/xreal.hpp"
#include "linalg/scaling.hpp"

namespace samaya {

LpSolveOptions lp_options_from_params(const Params& params) {
  LpSolveOptions options;
  options.simplex.primal_tol = params.primal_feasibility_tol;
  options.simplex.dual_tol = params.dual_feasibility_tol;
  options.simplex.time_limit = params.time_limit;
  return options;
}

namespace {

constexpr double kRayDropTolerance = 1e-12;

}  // namespace

// The model in computational form (minimization), scaled by `scaling`.
LpProblem make_problem(const Model& model, const Scaling& scaling) {
  const Index m = model.num_rows();
  const Index n = model.num_cols();
  const double sense = model.sense == ObjSense::kMaximize ? -1.0 : 1.0;
  LpProblem lp;
  lp.m = m;
  lp.n = n;
  lp.A = scaling.apply(model.A);
  lp.At = lp.A.transpose();
  const auto nt = static_cast<std::size_t>(n + m);
  lp.cost.assign(nt, 0.0);
  lp.lower.resize(nt);
  lp.upper.resize(nt);
  for (Index j = 0; j < n; ++j) {
    lp.cost[j] = sense * model.obj[j] * scaling.col[j];
    lp.lower[j] = model.col_lower[j] / scaling.col[j];
    lp.upper[j] = model.col_upper[j] / scaling.col[j];
  }
  for (Index i = 0; i < m; ++i) {
    lp.lower[n + i] = model.row_lower[i] * scaling.row[i];
    lp.upper[n + i] = model.row_upper[i] * scaling.row[i];
  }
  return lp;
}

namespace {

void add_stats(SimplexStats& total, const SimplexStats& s) {
  total.dual_iterations += s.dual_iterations;
  total.primal_iterations += s.primal_iterations;
  total.bound_flips += s.bound_flips;
  total.cost_shifts += s.cost_shifts;
  total.refactorizations += s.refactorizations;
  total.rebuilds_on_mismatch += s.rebuilds_on_mismatch;
}

}  // namespace

LpResult solve_lp(const Model& model, const LpSolveOptions& options, const Logger& log) {
  assert(model.Q.empty());
  const Index m = model.num_rows();
  const Index n = model.num_cols();
  const double sense = model.sense == ObjSense::kMaximize ? -1.0 : 1.0;

  const Scaling identity = Scaling::identity(m, n);
  const Scaling scaled = options.scale ? compute_scaling(model.A) : identity;
  const LpProblem scaled_lp = make_problem(model, scaled);
  LpProblem unscaled_lp;
  // The problem, scaling and simplex whose solution is reported.
  const LpProblem* lp_ptr = &scaled_lp;
  const Scaling* scaling_ptr = &scaled;

  auto simplex = std::make_unique<Simplex>(scaled_lp, options.simplex, log);
  LpResult result;
  result.status = simplex->solve();
  result.iterations = simplex->iterations();
  result.stats = simplex->stats();

  if (options.scale && result.status == SimplexStatus::kOptimal) {
    // Tolerances met in the scaled space can be violated once unscaled (reduced costs are divided
    // by the column scale factors). Continue on the unscaled problem from the optimal basis; this
    // usually needs no iterations. If it fails, the scaled solution is kept.
    unscaled_lp = make_problem(model, identity);
    auto polish = std::make_unique<Simplex>(unscaled_lp, options.simplex, log);
    const SimplexStatus status = polish->solve(simplex->status());
    result.iterations += polish->iterations();
    add_stats(result.stats, polish->stats());
    log.log(2, "lp: unscaled cleanup %s after %lld iterations", to_string(status),
            polish->iterations());
    if (status == SimplexStatus::kOptimal) {
      lp_ptr = &unscaled_lp;
      scaling_ptr = &identity;
      simplex = std::move(polish);
    }
  }
  const LpProblem& lp = *lp_ptr;
  const Scaling& scaling = *scaling_ptr;

  // Map the scaled solution back to the original model.
  const std::vector<double>& v = simplex->values();
  result.col_value.assign(v.begin(), v.begin() + n);
  scaling.unscale_cols(result.col_value);
  for (Index j = 0; j < n; ++j) {
    // Snap values that sit on a bound in scaled space exactly onto the original bound.
    if (v[j] == lp.lower[j]) result.col_value[j] = model.col_lower[j];
    if (v[j] == lp.upper[j]) result.col_value[j] = model.col_upper[j];
  }
  result.row_activity.assign(static_cast<std::size_t>(m), 0.0);
  model.A.multiply(result.col_value, result.row_activity);

  result.row_dual = simplex->duals();
  scaling.unscale_row_duals(result.row_dual);
  for (double& y : result.row_dual) y *= sense;
  result.col_dual.assign(static_cast<std::size_t>(n), 0.0);
  model.A.multiply_transpose(result.row_dual, result.col_dual);
  for (Index j = 0; j < n; ++j) result.col_dual[j] = model.obj[j] - result.col_dual[j];
  // Normalize signed zeros for cleaner output.
  for (double& value : result.row_dual) value += 0.0;
  for (double& value : result.col_dual) value += 0.0;
  for (double& value : result.col_value) value += 0.0;

  Xreal objective = model.obj_offset;  // double-double: the same bits on every platform
  for (Index j = 0; j < n; ++j) {
    objective += Xreal(model.obj[j]) * result.col_value[j];
  }
  result.objective = static_cast<double>(objective);

  if (result.status == SimplexStatus::kInfeasible) {
    result.dual_ray = simplex->dual_ray();
    scaling.unscale_row_duals(result.dual_ray);
    // Round-off entries (~1e-16 of the largest) on rows with an infinite bound on the wrong side
    // would make the certificate's range infinite; they carry no information, so drop them.
    double y_max = 0.0;
    for (const double y : result.dual_ray) y_max = std::max(y_max, std::fabs(y));
    for (double& y : result.dual_ray) {
      if (std::fabs(y) <= kRayDropTolerance * y_max) y = 0.0;
    }
  } else if (result.status == SimplexStatus::kUnbounded) {
    const std::vector<double>& ray = simplex->primal_ray();
    result.primal_ray.assign(ray.begin(), ray.begin() + n);
    scaling.unscale_cols(result.primal_ray);
  }
  log.log(2, "lp: %s after %lld simplex iterations", to_string(result.status), result.iterations);
  return result;
}

}  // namespace samaya
