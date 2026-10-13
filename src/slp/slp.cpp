#include "slp/slp.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "core/lp_verified.hpp"
#include "lp/lp_solver.hpp"
#include "samaya/verify.hpp"

namespace samaya {
namespace {

// Trust region radius delta, as a fraction of each column's scale (docs/SLP_PLAN.md 4.2).
constexpr double kInitialRadius = 0.1;
constexpr double kMinRadius = 1e-10;  // Below this the loop has stalled.
constexpr double kMaxRadius = 1e3;
// Ratio test (4.5): accept above kAcceptRatio, expand above kExpandRatio when the step reached
// the trust region, shrink below kShrinkRatio and on rejection.
constexpr double kAcceptRatio = 0.1;
constexpr double kExpandRatio = 0.75;
constexpr double kShrinkRatio = 0.25;
constexpr double kExpandFactor = 2.0;
constexpr double kShrinkFactor = 0.25;
// Penalty (4.5, 4.6): start at max(1, 10 max|c|), raise by 10 up to an operational cap.
constexpr double kPenaltyGrowth = 10.0;
constexpr double kMaxPenalty = 1e8;
constexpr int kMaxPenaltyRaisesPerIteration = 8;
// Steering: raise the penalty while the step achieves less than this share of the violation
// reduction the feasibility LP shows is possible, if that is at least kSteeringMinReduction of
// the current violation.
constexpr double kSteeringShare = 0.1;
constexpr double kSteeringMinReduction = 0.1;
constexpr int kMaxIterations = 200;
// Radius of the second-order correction (4.7) as a share of the step's radius: the residual it
// removes is second order in the step, so a much smaller region suffices.
constexpr double kCorrectionRadiusShare = 0.1;
// Componentwise step, relative to the column scales, below which the iterate has stopped moving.
constexpr double kStepTol = 1e-7;
// A finite column range larger than this many times max(1, |x0|) is not a useful scale.
constexpr double kRangeScaleCap = 1e6;
// Weighted violation below which elastics count as zero (an absolute floor; the final verdict
// is the verifier's, with its own relative tolerances).
constexpr double kViolationFloor = 1e-9;
// Predicted reductions below this (relative to the merit) are treated as zero: no division.
constexpr double kPredictedFloor = 1e-12;

class Slp {
 public:
  Slp(const Model& model, const Params& params, const Logger& log)
      : model_(model),
        params_(params),
        log_(log),
        lp_log_(params.log_level >= 3 ? params.log_level : 0),
        m_(model.num_rows()),
        n_(model.num_cols()),
        sigma_(model.sense == ObjSense::kMaximize ? -1.0 : 1.0) {
    nl_index_.assign(static_cast<std::size_t>(m_), -1);
    nonlinear_col_.assign(static_cast<std::size_t>(n_), 0);
    for (const RowQuadraticTerm& t : model.row_quadratic) {
      if (nl_index_[t.row] < 0) {
        nl_index_[t.row] = static_cast<Index>(nl_rows_.size());
        nl_rows_.push_back(t.row);
      }
      nonlinear_col_[t.col1] = 1;
      nonlinear_col_[t.col2] = 1;
    }
    cost_.resize(static_cast<std::size_t>(n_));
    for (Index j = 0; j < n_; ++j) cost_[j] = sigma_ * model.obj[j];
  }

  void run(Result& result);
  Model linearization_at(const std::vector<double>& center, double radius) {
    compute_scales_and_weights(center);
    return linearization(center, radius, 1.0, false);
  }

 private:
  Index num_nl() const { return static_cast<Index>(nl_rows_.size()); }
  double remaining() const { return params_.time_limit - timer_.seconds(); }

  std::vector<double> activities(const std::vector<double>& x) const;
  double violation(const std::vector<double>& x) const;  // V(x), weighted.
  double merit(const std::vector<double>& x, double mu) const;
  double objective(const std::vector<double>& x) const;  // sigma c'x

  Model linearization(const std::vector<double>& center, double radius, double mu,
                      bool feasibility_only) const;
  bool solve(const Model& lp, std::vector<VarStatus>* basis, LpResult& out);
  double elastic_sum(const LpResult& lp) const;  // Σ w_i (p_i + n_i)

  bool project_start(std::vector<double>& x, Result& result);
  void prove_linear_infeasibility(Result& result);
  void compute_scales_and_weights(const std::vector<double>& x0);
  bool stationary(const std::vector<double>& x, const LpResult& lp, Result& result);
  void report_point(const std::vector<double>& x, Status status, const std::string& message,
                    Result& result) const;

  const Model& model_;
  const Params& params_;
  const Logger& log_;
  const Logger lp_log_;
  const Timer timer_;
  Index m_;
  Index n_;
  double sigma_;
  std::vector<Index> nl_rows_;      // Rows with quadratic terms.
  std::vector<Index> nl_index_;     // Position in nl_rows_, or -1.
  std::vector<char> nonlinear_col_;
  std::vector<double> cost_;        // sigma * c: the minimization form.
  std::vector<double> scale_;       // s_j (4.2).
  // Columns without a usable finite range scale with their value: s_j = max(base, |x_j|), so a
  // flow that grows from 1 to 100 keeps steps proportional to itself (Haverly 3 from the
  // default start crawled 200 iterations at 0.2 units per step with a fixed scale).
  std::vector<double> base_scale_;
  std::vector<char> relative_scale_;
  std::vector<double> weight_;      // w_i per nonlinear row (4.4).
  long long simplex_iterations_ = 0;
};

std::vector<double> Slp::activities(const std::vector<double>& x) const {
  std::vector<double> h(static_cast<std::size_t>(m_), 0.0);
  model_.A.multiply(x, h);
  for (const RowQuadraticTerm& t : model_.row_quadratic) h[t.row] += t.value * x[t.col1] * x[t.col2];
  return h;
}

double Slp::violation(const std::vector<double>& x) const {
  const std::vector<double> h = activities(x);
  double v = 0.0;
  for (Index k = 0; k < num_nl(); ++k) {
    const Index i = nl_rows_[k];
    const double below = model_.row_lower[i] - h[i];
    const double above = h[i] - model_.row_upper[i];
    v += weight_[k] * std::max({0.0, below, above});
  }
  return v;
}

double Slp::objective(const std::vector<double>& x) const {
  double f = 0.0;
  for (Index j = 0; j < n_; ++j) f += cost_[j] * x[j];
  return f;
}

double Slp::merit(const std::vector<double>& x, double mu) const {
  return objective(x) + mu * violation(x);
}

// LP_k (4.3, 4.4): columns x then (p_k, n_k) per nonlinear row; each nonlinear row linearized
// at `center` (v x_a x_b -> v (c_b x_a + c_a x_b - c_a c_b)) with the elastics added; trust
// bounds on the columns in quadratic terms. The objective is sigma c'x plus mu w_k (p_k + n_k),
// or only Σ w_k (p_k + n_k) for the feasibility LP of the penalty steering.
Model Slp::linearization(const std::vector<double>& center, double radius, double mu,
                         bool feasibility_only) const {
  const Index nl = num_nl();
  const Index cols = n_ + 2 * nl;
  Model lp;
  lp.sense = ObjSense::kMinimize;
  lp.obj.assign(static_cast<std::size_t>(cols), 0.0);
  lp.col_lower.assign(static_cast<std::size_t>(cols), 0.0);
  lp.col_upper.assign(static_cast<std::size_t>(cols), kInf);
  lp.col_type.assign(static_cast<std::size_t>(cols), VarType::kContinuous);
  for (Index j = 0; j < n_; ++j) {
    if (!feasibility_only) lp.obj[j] = cost_[j];
    double lo = model_.col_lower[j];
    double up = model_.col_upper[j];
    if (nonlinear_col_[j] && lo < up) {
      lo = std::max(lo, center[j] - radius * scale_[j]);
      up = std::min(up, center[j] + radius * scale_[j]);
    }
    lp.col_lower[j] = lo;
    lp.col_upper[j] = up;
  }
  for (Index k = 0; k < nl; ++k) {
    const double price = (feasibility_only ? 1.0 : mu) * weight_[k];
    lp.obj[n_ + 2 * k] = price;
    lp.obj[n_ + 2 * k + 1] = price;
  }
  lp.row_lower = model_.row_lower;
  lp.row_upper = model_.row_upper;

  std::vector<Triplet> triplets;
  triplets.reserve(static_cast<std::size_t>(model_.A.nnz()) + 2 * model_.row_quadratic.size() +
                   2 * static_cast<std::size_t>(nl));
  const auto start = model_.A.col_start();
  const auto index = model_.A.row_index();
  const auto value = model_.A.values();
  for (Index j = 0; j < n_; ++j) {
    for (NnzIndex p = start[j]; p < start[j + 1]; ++p) triplets.push_back({index[p], j, value[p]});
  }
  for (const RowQuadraticTerm& t : model_.row_quadratic) {
    // On the diagonal (col1 == col2) the two entries sum to 2 v c_a, the derivative of v x_a^2.
    triplets.push_back({t.row, t.col1, t.value * center[t.col2]});
    triplets.push_back({t.row, t.col2, t.value * center[t.col1]});
    const double shift = t.value * center[t.col1] * center[t.col2];
    if (lp.row_lower[t.row] > -kInf) lp.row_lower[t.row] += shift;
    if (lp.row_upper[t.row] < kInf) lp.row_upper[t.row] += shift;
  }
  for (Index k = 0; k < nl; ++k) {
    triplets.push_back({nl_rows_[k], n_ + 2 * k, 1.0});
    triplets.push_back({nl_rows_[k], n_ + 2 * k + 1, -1.0});
  }
  lp.A = SparseMatrix::from_triplets(m_, cols, std::move(triplets));
  return lp;
}

// Solves an LP through solve_lp_verified with verification on; true only for a verified optimum
// (the helper also passes limit statuses, which claim nothing).
bool Slp::solve(const Model& lp, std::vector<VarStatus>* basis, LpResult& out) {
  Params p = params_;
  p.verify = true;
  p.time_limit = remaining();
  Result scratch;
  bool verified = false;
  const std::vector<VarStatus> none;
  const bool warm = basis != nullptr &&
                    basis->size() == static_cast<std::size_t>(lp.num_cols() + lp.num_rows());
  out = solve_lp_verified(lp, p, lp_log_, scratch, verified, warm ? *basis : none);
  simplex_iterations_ += out.iterations;
  const bool ok = verified && out.status == SimplexStatus::kOptimal;
  if (ok && basis != nullptr) *basis = out.basis;
  return ok;
}

double Slp::elastic_sum(const LpResult& lp) const {
  double v = 0.0;
  for (Index k = 0; k < num_nl(); ++k) {
    v += weight_[k] * (lp.col_value[n_ + 2 * k] + lp.col_value[n_ + 2 * k + 1]);
  }
  return v;
}

void Slp::compute_scales_and_weights(const std::vector<double>& x0) {
  scale_.resize(static_cast<std::size_t>(n_));
  base_scale_.resize(static_cast<std::size_t>(n_));
  relative_scale_.assign(static_cast<std::size_t>(n_), 0);
  for (Index j = 0; j < n_; ++j) {
    const double base = std::max(1.0, std::fabs(x0[j]));
    const double range = model_.col_upper[j] - model_.col_lower[j];
    const bool use_range = std::isfinite(range) && range > 0.0 && range <= kRangeScaleCap * base;
    scale_[j] = base_scale_[j] = use_range ? range : base;
    relative_scale_[j] = use_range ? 0 : 1;
  }
  // w_k = 1 / max(1, largest |derivative| of row k at x0): rows in different units count alike.
  std::vector<double> largest(static_cast<std::size_t>(num_nl()), 0.0);
  const auto start = model_.A.col_start();
  const auto index = model_.A.row_index();
  const auto value = model_.A.values();
  std::vector<std::vector<std::pair<Index, double>>> grad(static_cast<std::size_t>(num_nl()));
  for (Index j = 0; j < n_; ++j) {
    for (NnzIndex p = start[j]; p < start[j + 1]; ++p) {
      if (nl_index_[index[p]] >= 0) grad[nl_index_[index[p]]].push_back({j, value[p]});
    }
  }
  for (const RowQuadraticTerm& t : model_.row_quadratic) {
    grad[nl_index_[t.row]].push_back({t.col1, t.value * x0[t.col2]});
    grad[nl_index_[t.row]].push_back({t.col2, t.value * x0[t.col1]});
  }
  weight_.resize(static_cast<std::size_t>(num_nl()));
  for (Index k = 0; k < num_nl(); ++k) {
    auto& g = grad[k];
    std::sort(g.begin(), g.end());
    double best = 0.0;
    for (std::size_t s = 0; s < g.size();) {
      double sum = 0.0;
      const Index j = g[s].first;
      for (; s < g.size() && g[s].first == j; ++s) sum += g[s].second;
      best = std::max(best, std::fabs(sum));
    }
    weight_[k] = 1.0 / std::max(1.0, best);
  }
}

// 4.1: the projection LP, min Σ_j |x_j - x0_j| / s_j over the nonlinear columns subject to the
// linear rows and all bounds (the nonlinear rows are left out). Its answer is the first iterate.
bool Slp::project_start(std::vector<double>& x, Result& result) {
  std::vector<Index> linear_rows;
  for (Index i = 0; i < m_; ++i) {
    if (nl_index_[i] < 0) linear_rows.push_back(i);
  }
  // Pulled toward the start: the columns in quadratic terms, and every column the caller gave
  // a value for (a column left free here lands on an arbitrary vertex and starts the loop with
  // a violation the start did not have).
  const bool has_start = params_.mip_start.size() == static_cast<std::size_t>(n_);
  std::vector<Index> targets;
  for (Index j = 0; j < n_; ++j) {
    const bool given = has_start && !std::isnan(params_.mip_start[j]);
    if ((nonlinear_col_[j] || given) && model_.col_lower[j] < model_.col_upper[j]) {
      targets.push_back(j);
    }
  }
  const auto r = static_cast<Index>(linear_rows.size());
  const auto t = static_cast<Index>(targets.size());
  Model lp;
  lp.sense = ObjSense::kMinimize;
  const Index cols = n_ + 2 * t;
  lp.obj.assign(static_cast<std::size_t>(cols), 0.0);
  lp.col_lower = model_.col_lower;
  lp.col_upper = model_.col_upper;
  lp.col_lower.resize(static_cast<std::size_t>(cols), 0.0);
  lp.col_upper.resize(static_cast<std::size_t>(cols), kInf);
  lp.col_type.assign(static_cast<std::size_t>(cols), VarType::kContinuous);
  std::vector<Index> new_row(static_cast<std::size_t>(m_), -1);
  for (Index k = 0; k < r; ++k) {
    new_row[linear_rows[k]] = k;
    lp.row_lower.push_back(model_.row_lower[linear_rows[k]]);
    lp.row_upper.push_back(model_.row_upper[linear_rows[k]]);
  }
  std::vector<Triplet> triplets;
  const auto start = model_.A.col_start();
  const auto index = model_.A.row_index();
  const auto value = model_.A.values();
  for (Index j = 0; j < n_; ++j) {
    for (NnzIndex p = start[j]; p < start[j + 1]; ++p) {
      if (new_row[index[p]] >= 0) triplets.push_back({new_row[index[p]], j, value[p]});
    }
  }
  for (Index k = 0; k < t; ++k) {
    // x_j - dp_k + dn_k = x0_j
    const Index j = targets[k];
    const Index row = r + k;
    triplets.push_back({row, j, 1.0});
    triplets.push_back({row, n_ + 2 * k, -1.0});
    triplets.push_back({row, n_ + 2 * k + 1, 1.0});
    lp.row_lower.push_back(x[j]);
    lp.row_upper.push_back(x[j]);
    lp.obj[n_ + 2 * k] = 1.0 / scale_[j];
    lp.obj[n_ + 2 * k + 1] = 1.0 / scale_[j];
  }
  lp.A = SparseMatrix::from_triplets(r + t, cols, std::move(triplets));
  LpResult out;
  if (solve(lp, nullptr, out)) {
    x.assign(out.col_value.begin(), out.col_value.begin() + n_);
    return true;
  }
  if (out.status == SimplexStatus::kInfeasible) {
    prove_linear_infeasibility(result);
  } else {
    report_point(x, out.status == SimplexStatus::kTimeLimit ? Status::kTimeLimit
                                                             : Status::kNumericalError,
                 "the starting-point projection LP failed", result);
    result.col_value.clear();
  }
  return false;
}

// The linear rows and bounds alone are infeasible: solve that LP for a Farkas certificate and
// report it on the original model (zero on the nonlinear rows, which it does not use).
void Slp::prove_linear_infeasibility(Result& result) {
  Model lp = model_;
  lp.row_quadratic.clear();
  for (const Index i : nl_rows_) {
    lp.row_lower[i] = -kInf;
    lp.row_upper[i] = kInf;
  }
  std::fill(lp.obj.begin(), lp.obj.end(), 0.0);
  LpResult out;
  solve(lp, nullptr, out);
  result.status = Status::kNumericalError;
  result.message = "the linear rows look infeasible, but no certificate was found";
  if (out.status == SimplexStatus::kInfeasible) {
    const VerifyReport report = verify_infeasibility(model_, out.dual_ray);
    if (report.ok) {
      result.status = Status::kInfeasible;
      result.verified = true;
      result.message = "the linear rows and bounds are infeasible";
      result.infeasibility_certificate = out.dual_ray;
    }
  }
}

// The final check (section 3): x feasible for the original model, and the duals of the LP at x
// (which reached no better point, so x is optimal for it with zero elastics and its trust bounds
// inactive at the center) satisfying the original model's KKT conditions.
bool Slp::stationary(const std::vector<double>& x, const LpResult& lp, Result& result) {
  std::vector<double> y(lp.row_dual.begin(), lp.row_dual.begin() + m_);
  for (double& v : y) v *= sigma_;
  const VerifyReport primal = verify_primal(model_, x);
  if (!primal.ok) return false;
  const VerifyReport kkt = verify_qcp_kkt(model_, x, y);
  if (!kkt.ok) {
    log_.log(2, "SLP: point not certified: %s", kkt.message.c_str());
    return false;
  }
  report_point(x, Status::kFirstOrderStationary,
               "feasible, first-order stationary (KKT) point; local optimality not certified",
               result);
  result.verified = true;
  result.max_dual_violation = kkt.max_dual_violation;
  result.row_dual = y;
  // d = c - A'y - J(x)'y in the model's own sense, as verify_qcp_kkt checked it.
  result.col_dual.assign(static_cast<std::size_t>(n_), 0.0);
  model_.A.multiply_transpose(y, result.col_dual);
  for (Index j = 0; j < n_; ++j) result.col_dual[j] = model_.obj[j] - result.col_dual[j];
  for (const RowQuadraticTerm& t : model_.row_quadratic) {
    result.col_dual[t.col1] -= t.value * x[t.col2] * y[t.row];
    result.col_dual[t.col2] -= t.value * x[t.col1] * y[t.row];
  }
  return true;
}

void Slp::report_point(const std::vector<double>& x, Status status, const std::string& message,
                       Result& result) const {
  const VerifyReport primal = verify_primal(model_, x);
  result.status = status;
  result.message = message;
  result.col_value = x;
  result.row_activity = activities(x);
  result.objective = primal.objective;
  result.dual_bound = std::numeric_limits<double>::quiet_NaN();
  result.max_primal_violation = std::max(primal.max_bound_violation, primal.max_row_violation);
  // The worst quadratic row, measured as the verifier measures rows.
  double worst = 0.0;
  const std::vector<double>& h = result.row_activity;
  for (const Index i : nl_rows_) {
    const double lo = model_.row_lower[i];
    const double up = model_.row_upper[i];
    double v = 0.0;
    if (h[i] < lo) v = (lo - h[i]) / (1.0 + std::fabs(lo));
    if (h[i] > up) v = (h[i] - up) / (1.0 + std::fabs(up));
    worst = std::max(worst, v);
  }
  result.max_nonlinear_violation = worst;
  result.verified = status != Status::kFirstOrderStationary && primal.ok;
  result.simplex_iterations = simplex_iterations_;
}

void Slp::run(Result& result) {
  // 4.1: the starting point.
  std::vector<double> x(static_cast<std::size_t>(n_), 0.0);
  const bool has_start = params_.mip_start.size() == static_cast<std::size_t>(n_);
  for (Index j = 0; j < n_; ++j) {
    const double lo = model_.col_lower[j];
    const double up = model_.col_upper[j];
    double v = has_start ? params_.mip_start[j] : std::numeric_limits<double>::quiet_NaN();
    if (std::isnan(v)) {
      // A column in a quadratic term starts inside its bounds: at a bound, a pool's flows are
      // zero and its quality has no first-order effect, so the LP sees no way to move it (the
      // empty-pool trap, docs/SLP_PLAN.md 4.9). One unit above a lone lower bound (below a lone
      // upper one) gives the projection LP flows to balance.
      if (std::isfinite(lo) && std::isfinite(up)) v = nonlinear_col_[j] ? 0.5 * (lo + up) : lo;
      else if (std::isfinite(lo)) v = nonlinear_col_[j] ? lo + 1.0 : lo;
      else if (std::isfinite(up)) v = nonlinear_col_[j] ? up - 1.0 : up;
      else v = 0.0;
    }
    x[j] = std::clamp(v, lo, up);
  }
  compute_scales_and_weights(x);
  if (!project_start(x, result)) return;
  compute_scales_and_weights(x);

  double mu = 1.0;
  for (const double c : cost_) mu = std::max(mu, 10.0 * std::fabs(c));
  double radius = kInitialRadius;
  std::vector<VarStatus> basis;
  long long iteration = 0;
  for (; iteration < kMaxIterations; ++iteration) {
    if (remaining() <= 0.0) {
      report_point(x, Status::kTimeLimit, "time limit before stationarity", result);
      result.slp_iterations = iteration;
      return;
    }
    const double v_now = violation(x);
    LpResult step;
    bool solved = solve(linearization(x, radius, mu, false), &basis, step);
    // 4.6: steer the penalty while the step leaves violation the feasibility LP could remove.
    for (int raise = 0; solved && raise < kMaxPenaltyRaisesPerIteration && mu < kMaxPenalty &&
                        v_now > kViolationFloor && elastic_sum(step) > kViolationFloor;
         ++raise) {
      LpResult feasibility;
      std::vector<VarStatus> feasibility_basis = basis;
      if (!solve(linearization(x, radius, mu, true), &feasibility_basis, feasibility)) break;
      const double possible = v_now - elastic_sum(feasibility);
      const double achieved = v_now - elastic_sum(step);
      if (possible < kSteeringMinReduction * v_now || achieved >= kSteeringShare * possible) break;
      mu = std::min(kMaxPenalty, mu * kPenaltyGrowth);
      solved = solve(linearization(x, radius, mu, false), &basis, step);
    }
    if (!solved) {
      if (step.status == SimplexStatus::kUnbounded) {
        // Unbounded along columns outside the quadratic terms (those have trust bounds): the
        // quadratic rows do not change along the ray, so the original model is unbounded too.
        std::vector<double> ray(step.primal_ray.begin(), step.primal_ray.begin() + n_);
        bool linear_only = true;
        for (Index j = 0; j < n_; ++j) linear_only = linear_only && (!nonlinear_col_[j] || ray[j] == 0.0);
        if (linear_only && verify_unbounded_ray(model_, ray).ok) {
          report_point(x, Status::kUnbounded, "unbounded along columns outside quadratic terms",
                       result);
          result.unbounded_ray = std::move(ray);
          result.verified = true;
          result.slp_iterations = iteration;
          return;
        }
      }
      report_point(x, step.status == SimplexStatus::kTimeLimit ? Status::kTimeLimit
                                                                : Status::kNumericalError,
                   std::string("an LP of the SLP loop ended ") + to_string(step.status), result);
      result.slp_iterations = iteration;
      return;
    }

    std::vector<double> trial(step.col_value.begin(), step.col_value.begin() + n_);
    const double phi = merit(x, mu);
    const double predicted = objective(x) - objective(trial) + mu * (v_now - elastic_sum(step));
    double moved = 0.0;
    double moved_nonlinear = 0.0;
    for (Index j = 0; j < n_; ++j) {
      const double s = std::fabs(trial[j] - x[j]) / scale_[j];
      moved = std::max(moved, s);
      if (nonlinear_col_[j]) moved_nonlinear = std::max(moved_nonlinear, s);
    }
    if (predicted <= kPredictedFloor * (1.0 + std::fabs(phi)) || moved <= kStepTol) {
      // No predicted progress: x is optimal for its own LP. Certify it, or shrink and retry.
      if (stationary(x, step, result)) {
        result.slp_iterations = iteration + 1;
        log_.log(1, "SLP: first-order stationary point after %lld iterations, objective %.12g",
                 iteration + 1, result.objective);
        return;
      }
      radius *= kShrinkFactor;
    } else {
      double rho = (phi - merit(trial, mu)) / predicted;
      // 4.7: second-order correction. A step that satisfies the linearized rows leaves the
      // curvature's residual (v s_a s_b) in them, which the next LP again predicts to remove:
      // the ratio stays low and the radius cannot grow (Haverly 3 from the default start
      // crawled at a ratio near 0.63). One LP linearized at the trial point, over a small
      // radius, removes that residual; its point replaces the trial when the ratio improves.
      if (rho < kExpandRatio && violation(trial) > kViolationFloor) {
        LpResult correction;
        std::vector<VarStatus> correction_basis = basis;
        if (solve(linearization(trial, kCorrectionRadiusShare * radius, mu, false),
                  &correction_basis, correction)) {
          std::vector<double> corrected(correction.col_value.begin(),
                                        correction.col_value.begin() + n_);
          const double rho_corrected = (phi - merit(corrected, mu)) / predicted;
          if (rho_corrected > rho) {
            trial = std::move(corrected);
            rho = rho_corrected;
          }
        }
      }
      if (rho > kAcceptRatio) {
        const bool at_boundary = moved_nonlinear >= 0.99 * radius;
        x = trial;
        for (Index j = 0; j < n_; ++j) {
          if (relative_scale_[j]) scale_[j] = std::max(base_scale_[j], std::fabs(x[j]));
        }
        if (rho > kExpandRatio && at_boundary) {
          radius = std::min(kMaxRadius, radius * kExpandFactor);
        } else if (rho < kShrinkRatio) {
          radius *= kShrinkFactor;
        }
      } else {
        radius *= kShrinkFactor;
      }
      log_.log(2, "SLP %lld: objective %.10g violation %.3e radius %.2e penalty %.1e ratio %.3f",
               iteration + 1, sigma_ * objective(x), violation(x), radius, mu, rho);
    }
    if (radius < kMinRadius) break;
  }
  const bool feasible = verify_primal(model_, x).ok;
  report_point(x, Status::kIterationLimit,
               feasible ? "stalled before stationarity (feasible point)"
                        : "no feasible point found from this start",
               result);
  result.slp_iterations = iteration;
  log_.log(1, "SLP: %s after %lld iterations", result.message.c_str(), iteration);
}

}  // namespace

void solve_qcp(const Model& model, const Params& params, const Logger& log, Result& result) {
  Slp(model, params, log).run(result);
}

Model slp_linearization_for_test(const Model& model, const std::vector<double>& center,
                                 double radius) {
  const Params params;
  const Logger quiet(0);
  return Slp(model, params, quiet).linearization_at(center, radius);
}

}  // namespace samaya
