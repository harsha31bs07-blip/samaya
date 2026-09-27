#include "qp/miqp.hpp"

#include <algorithm>
#include <cmath>
#include <queue>
#include <utility>

#include "lp/ipm.hpp"
#include "samaya/solver.hpp"
#include "samaya/verify.hpp"

namespace samaya {

namespace {

// The rounding heuristic runs at the root and every this many nodes.
constexpr long long kRoundingFrequency = 20;
// Relaxation values come from interior-point solutions converged to a relative gap of about
// 1e-9; a node's bound is its value less this relative margin, so rounding in the relaxation can
// never prune a node that holds a better solution.
constexpr double kBoundMargin = 1e-7;

struct Node {
  std::vector<double> lower;
  std::vector<double> upper;
  double bound = -kInf;  // Minimization sense.
  long long id = 0;
};

struct Worse {
  bool operator()(const Node& a, const Node& b) const {
    return a.bound > b.bound || (a.bound == b.bound && a.id > b.id);
  }
};

}  // namespace

MiqpOutcome solve_miqp(const Model& model, const MiqpOptions& options, const Logger& log) {
  const Timer timer;
  MiqpOutcome out;
  const Index n = model.num_cols();
  const double sense = model.sense == ObjSense::kMaximize ? -1.0 : 1.0;
  double incumbent = kInf;  // Minimization sense, offset included.
  bool incomplete = false;  // A node could be neither solved nor proven infeasible.

  // Relaxation of the model with the given column bounds, integers relaxed.
  Model relaxed = model;
  relaxed.col_type.assign(static_cast<std::size_t>(n), VarType::kContinuous);
  const auto remaining = [&] { return options.time_limit - timer.seconds(); };

  // Accepts x as the incumbent if it passes verify_primal and improves.
  const auto try_incumbent = [&](const std::vector<double>& x) {
    VerifyTolerances tol;
    tol.integrality = options.integrality_tol;
    const VerifyReport report = verify_primal(model, x, tol);
    if (!report.ok) return false;
    const double value = sense * report.objective;
    if (!(value < incumbent)) return false;
    incumbent = value;
    out.x = x;
    log.log(2, "MIQP: incumbent %.10g after %lld nodes", report.objective, out.nodes);
    return true;
  };

  // Solves the relaxation under the bounds; false if it has no solution (then `proven` says
  // whether infeasibility is certified).
  const auto solve_relaxation = [&](const std::vector<double>& lower,
                                    const std::vector<double>& upper, std::vector<double>& x,
                                    double& value, bool& proven, IpmStatus& status) {
    relaxed.col_lower = lower;
    relaxed.col_upper = upper;
    proven = false;
    for (Index j = 0; j < n; ++j) {
      if (lower[j] > upper[j]) {
        proven = true;  // Crossed bounds.
        return false;
      }
    }
    IpmOptions ipm_options;
    ipm_options.time_limit = remaining();
    const IpmResult ipm = solve_ipm(relaxed, ipm_options, Logger(0));
    out.ipm_iterations += ipm.iterations;
    status = ipm.status;
    if (ipm.status == IpmStatus::kOptimal) {
      x = ipm.x;
      const VerifyReport report = verify_primal(relaxed, x);
      value = sense * report.objective;
      value -= kBoundMargin * (1.0 + std::fabs(value));
      return true;
    }
    if (ipm.status == IpmStatus::kNotConvex || ipm.status == IpmStatus::kTimeLimit) return false;
    // Not converged: the constraints alone, by the dual simplex, decide feasibility.
    Model feasibility = relaxed;
    feasibility.Q = SparseMatrix();
    std::fill(feasibility.obj.begin(), feasibility.obj.end(), 0.0);
    Params params;
    params.log_level = 0;
    params.lp_method = LpMethod::kDualSimplex;
    params.time_limit = std::max(0.0, remaining());
    const Result lp = Solver(params).solve(feasibility);
    proven = lp.status == Status::kInfeasible && lp.verified;
    return false;
  };

  // Rounds the integers of x, fixes them and solves the continuous rest.
  const auto round_and_solve = [&](const std::vector<double>& x, const Node& node) {
    std::vector<double> lower = node.lower;
    std::vector<double> upper = node.upper;
    for (Index j = 0; j < n; ++j) {
      if (model.col_type[j] != VarType::kInteger) continue;
      const double v = std::clamp(std::round(x[j]), lower[j], upper[j]);
      lower[j] = upper[j] = v;
    }
    std::vector<double> y;
    double value = 0.0;
    bool proven = false;
    IpmStatus status = IpmStatus::kNumericalError;
    if (solve_relaxation(lower, upper, y, value, proven, status)) {
      for (Index j = 0; j < n; ++j) {
        if (model.col_type[j] == VarType::kInteger) y[j] = lower[j];
      }
      try_incumbent(y);
    }
  };

  const auto gap_closed = [&](double bound) {
    return bound >= incumbent - std::max(options.abs_gap,
                                         options.rel_gap * std::fabs(incumbent));
  };

  std::priority_queue<Node, std::vector<Node>, Worse> open;
  long long next_id = 0;
  open.push({model.col_lower, model.col_upper, -kInf, next_id++});
  double best_bound = -kInf;
  while (!open.empty()) {
    if (remaining() <= 0.0) {
      out.status = Status::kTimeLimit;
      break;
    }
    if (options.node_limit >= 0 && out.nodes >= options.node_limit) {
      out.status = Status::kNodeLimit;
      break;
    }
    Node node = open.top();
    open.pop();
    if (incumbent < kInf && gap_closed(node.bound)) continue;
    ++out.nodes;

    std::vector<double> x;
    double value = 0.0;
    bool proven = false;
    IpmStatus status = IpmStatus::kNumericalError;
    if (!solve_relaxation(node.lower, node.upper, x, value, proven, status)) {
      if (status == IpmStatus::kNotConvex) {
        out.status = Status::kNotConvex;
        return out;
      }
      if (status == IpmStatus::kTimeLimit) {
        out.status = Status::kTimeLimit;
        open.push(node);  // Its bound still counts.
        break;
      }
      if (!proven) incomplete = true;
      continue;
    }
    if (incumbent < kInf && gap_closed(value)) continue;

    // Branch on the most fractional integer column.
    Index branch = -1;
    double best_frac = options.integrality_tol;
    for (Index j = 0; j < n; ++j) {
      if (model.col_type[j] != VarType::kInteger) continue;
      const double frac = std::fabs(x[j] - std::round(x[j]));
      if (frac > best_frac) {
        best_frac = frac;
        branch = j;
      }
    }
    if (branch < 0) {
      // Integral relaxation: fix the rounded integers and re-solve the continuous rest, so the
      // point is exactly integral.
      round_and_solve(x, node);
      continue;
    }
    if (out.nodes == 1 || out.nodes % kRoundingFrequency == 0) round_and_solve(x, node);
    Node down{node.lower, node.upper, value, next_id++};
    down.upper[branch] = std::floor(x[branch]);
    Node up{node.lower, node.upper, value, next_id++};
    up.lower[branch] = std::ceil(x[branch]);
    open.push(std::move(down));
    open.push(std::move(up));
  }
  if (open.empty() && out.status == Status::kNotSolved) {
    if (incomplete) {
      out.status = Status::kNumericalError;
    } else {
      out.status = incumbent < kInf ? Status::kOptimal : Status::kInfeasible;
    }
    best_bound = incumbent;
  } else {
    best_bound = incumbent;
    while (!open.empty()) {
      best_bound = std::min(best_bound, open.top().bound);
      open.pop();
    }
  }
  if (incomplete && out.status != Status::kNumericalError) {
    best_bound = -kInf;  // A dropped node has no bound.
  }
  out.objective = incumbent < kInf ? sense * incumbent : sense * kInf;
  out.bound = sense * best_bound;
  log.log(1, "MIQP: %s after %lld nodes, %lld interior-point iterations, %.2f s",
          to_string(out.status), out.nodes, out.ipm_iterations, timer.seconds());
  return out;
}

}  // namespace samaya
