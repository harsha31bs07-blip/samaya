#include "lp/ranging.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

#include "lp/crossover.hpp"
#include "lp/lp_solver.hpp"
#include "samaya/verify.hpp"

namespace samaya {
namespace {

constexpr Index kRangingMaxRows = 100000;
// Tiny tableau entries would give unreliable, very large sensitivity limits.
constexpr double kRangingPivotTol = 1e-9;
// Crossover may select a different vertex, but must preserve the reported optimum.
constexpr double kRangingAgreementTol = 1e-6;

Ranging unavailable(const std::string& message) {
  Ranging ranging;
  ranging.message = message;
  return ranging;
}

// Slack slightly below zero within the feasibility tolerance is roundoff, not a negative step.
double feasible_slack(double slack, double tolerance) {
  return slack >= -tolerance ? std::max(0.0, slack) : slack;
}

}  // namespace

Ranging compute_ranging(const Model& model, Result& result, const Params& params,
                        const Logger& log) {
  if (model.problem_class() != ProblemClass::kLP) {
    return unavailable("ranging is available for LPs only");
  }
  if (result.status != Status::kOptimal) return unavailable("ranging requires an optimal LP");
  const Index n = model.num_cols();
  const Index m = model.num_rows();
  if (result.col_value.size() != static_cast<std::size_t>(n) ||
      result.row_dual.size() != static_cast<std::size_t>(m)) {
    return unavailable("ranging requires column values and row duals of the original model");
  }
  if (m > kRangingMaxRows) return unavailable("ranging row limit exceeded (100000)");
  if (params.time_limit <= 0.0) return unavailable("ranging time limit used up");
  const Timer timer;
  if (!verify_lp_optimality(model, result.col_value, result.row_dual).ok) {
    return unavailable("ranging requires a verified optimal LP answer");
  }

  const LpProblem lp = make_problem(model, Scaling::identity(m, n));
  const auto start = crossover_basis(model, result.col_value, result.row_dual);
  SimplexOptions options = lp_options_from_params(params).simplex;
  options.perturb = false;
  options.time_limit -= timer.seconds();
  if (options.time_limit <= 0.0) return unavailable("ranging time limit used up");
  const Logger quiet(0);
  Simplex simplex(lp, options, quiet);
  const SimplexStatus status = simplex.solve(start);
  if (status != SimplexStatus::kOptimal) {
    return unavailable(std::string("ranging basis simplex: ") + to_string(status));
  }
  const double sense = model.sense == ObjSense::kMaximize ? -1.0 : 1.0;
  std::vector<double> x(simplex.values().begin(), simplex.values().begin() + n);
  std::vector<double> y = simplex.duals();
  for (double& value : y) value *= sense;
  const VerifyReport report = verify_lp_optimality(model, x, y);
  if (!report.ok) return unavailable("ranging basis failed optimality verification: " + report.message);
  if (!std::isfinite(result.objective) ||
      std::fabs(report.objective - result.objective) >
          kRangingAgreementTol * (1.0 + std::fabs(result.objective))) {
    return unavailable("ranging basis objective differs from the reported optimum");
  }
  Ranging ranging;
  for (Index j = 0; j < n; ++j) {
    if (std::fabs(x[j] - result.col_value[j]) >
        kRangingAgreementTol * (1.0 + std::fabs(result.col_value[j]))) {
      ranging.message = "the reported solution is the ranging basis's vertex, an alternative "
                        "optimum to the one first found";
      break;
    }
  }
  const double nan = std::numeric_limits<double>::quiet_NaN();
  ranging.cost_lower.resize(static_cast<std::size_t>(n));
  ranging.cost_upper.resize(static_cast<std::size_t>(n));
  ranging.col_level_lower.assign(static_cast<std::size_t>(n), nan);
  ranging.col_level_upper.assign(static_cast<std::size_t>(n), nan);
  ranging.row_level_lower.assign(static_cast<std::size_t>(m), nan);
  ranging.row_level_upper.assign(static_cast<std::size_t>(m), nan);

  const auto& statuses = simplex.status();
  const auto& basic = simplex.basic();
  const auto& values = simplex.values();
  std::vector<double> reduced = simplex.reduced_costs();
  std::vector<Index> position(static_cast<std::size_t>(n + m), -1);
  for (Index r = 0; r < m; ++r) position[basic[r]] = r;
  for (Index k = 0; k < n + m; ++k) {
    if (statuses[k] == VarStatus::kAtLower && reduced[k] < 0.0 &&
        reduced[k] >= -params.dual_feasibility_tol) reduced[k] = 0.0;
    if (statuses[k] == VarStatus::kAtUpper && reduced[k] > 0.0 &&
        reduced[k] <= params.dual_feasibility_tol) reduced[k] = 0.0;
  }
  std::vector<double> alpha;
  for (Index j = 0; j < n; ++j) {
    if (timer.seconds() >= params.time_limit) return unavailable("ranging time limit used up");
    double lo = -kInf;
    double hi = kInf;
    if (statuses[j] == VarStatus::kBasic) {
      simplex.tableau_row(position[j], alpha);
      double down = kInf;
      double up = kInf;
      for (Index k = 0; k < n + m; ++k) {
        if (statuses[k] == VarStatus::kBasic || lp.lower[k] == lp.upper[k] ||
            std::fabs(alpha[k]) <= kRangingPivotTol) continue;
        if (statuses[k] == VarStatus::kAtZero) {
          down = up = 0.0;
          continue;
        }
        // d_k changes by -delta * alpha_k. Convert upper-bound dual signs to lower signs.
        const double side = statuses[k] == VarStatus::kAtLower ? 1.0 : -1.0;
        const double slope = side * alpha[k];
        const double slack = std::max(0.0, side * reduced[k]);
        if (slope > 0.0) up = std::min(up, slack / slope);
        else down = std::min(down, slack / -slope);
      }
      lo = lp.cost[j] - down;
      hi = lp.cost[j] + up;
    } else if (lp.lower[j] != lp.upper[j]) {
      if (statuses[j] == VarStatus::kAtLower) lo = lp.cost[j] - reduced[j];
      else if (statuses[j] == VarStatus::kAtUpper) hi = lp.cost[j] - reduced[j];
      else lo = hi = lp.cost[j];
    }
    ranging.cost_lower[j] = sense > 0.0 ? lo : -hi;
    ranging.cost_upper[j] = sense > 0.0 ? hi : -lo;
  }

  BasisFactor factor(lp.A);
  if (factor.factorize(basic) != 0) return unavailable("ranging basis factorization failed");
  std::vector<double> column(static_cast<std::size_t>(m));
  const auto col_start = lp.A.col_start();
  const auto row_index = lp.A.row_index();
  const auto coefficients = lp.A.values();
  for (Index k = 0; k < n + m; ++k) {
    if (timer.seconds() >= params.time_limit) return unavailable("ranging time limit used up");
    // Basic variables have no binding bound, and neither does a free one nonbasic at zero.
    if (statuses[k] == VarStatus::kBasic || statuses[k] == VarStatus::kAtZero) continue;
    std::fill(column.begin(), column.end(), 0.0);
    if (k < n) {
      for (NnzIndex p = col_start[k]; p < col_start[k + 1]; ++p) {
        column[row_index[p]] = coefficients[p];
      }
    } else {
      column[k - n] = -1.0;
    }
    factor.ftran(column);
    double down = kInf;
    double up = kInf;
    for (Index r = 0; r < m; ++r) {
      if (std::fabs(column[r]) <= kRangingPivotTol) continue;
      const Index b = basic[r];
      const double lower_slack = feasible_slack(values[b] - lp.lower[b], params.primal_feasibility_tol);
      const double upper_slack = feasible_slack(lp.upper[b] - values[b], params.primal_feasibility_tol);
      if (lower_slack < 0.0 || upper_slack < 0.0) {
        return unavailable("ranging basis violates primal bounds");
      }
      // A positive transformed coefficient makes the basic value decrease as the bound rises.
      if (column[r] > 0.0) {
        up = std::min(up, lower_slack / column[r]);
        down = std::min(down, upper_slack / column[r]);
      } else {
        up = std::min(up, upper_slack / -column[r]);
        down = std::min(down, lower_slack / -column[r]);
      }
    }
    // A fixed variable moves both bounds together; otherwise keep the opposite bound intact.
    if (lp.lower[k] != lp.upper[k]) {
      if (statuses[k] == VarStatus::kAtLower) up = std::min(up, lp.upper[k] - values[k]);
      if (statuses[k] == VarStatus::kAtUpper) down = std::min(down, values[k] - lp.lower[k]);
    }
    auto& lower = k < n ? ranging.col_level_lower : ranging.row_level_lower;
    auto& upper = k < n ? ranging.col_level_upper : ranging.row_level_upper;
    const Index index = k < n ? k : k - n;
    lower[index] = values[k] - down;
    upper[index] = values[k] + up;
  }
  ranging.available = true;

  // Report the basis's own solution (verified above, same objective), so the duals and binding
  // bounds the caller reads are the ones the ranges are for. Reduced costs back in the original
  // sense: d = c - A'y, the min form's d times the sense.
  const std::vector<double>& min_reduced = simplex.reduced_costs();
  result.col_value = std::move(x);
  result.row_activity.assign(values.begin() + n, values.end());
  result.row_dual = std::move(y);
  result.col_dual.resize(static_cast<std::size_t>(n));
  for (Index j = 0; j < n; ++j) result.col_dual[j] = sense * min_reduced[j];
  result.objective = result.dual_bound = report.objective;
  result.max_primal_violation = std::max(report.max_bound_violation, report.max_row_violation);
  result.max_dual_violation = report.max_dual_violation;

  log.log(1, "Ranging: %d columns, %d rows, %lld simplex iterations", n, m, simplex.iterations());
  return ranging;
}

}  // namespace samaya
