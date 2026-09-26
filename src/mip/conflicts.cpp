// Conflict analysis from node LPs (dual proofs; Witzig, Berthold and Heinz, "Experiments with
// conflict analysis in mixed integer programming", CPAIOR 2017).
//
// With the row activities s = A x inside the row bounds, any multipliers y give the identity
// (A'y)'x = y's, so  (A'y)'x  lies within [min y's, max y's] over the row bounds. For the Farkas
// ray of an infeasible node LP that range and the range of (A'y)'x over the node's column bounds
// are disjoint. For a node LP whose value exceeds the cutoff, the LP duals give
// c'x = y's + d'x with d = c - A'y, so every improving solution satisfies
// d'x <= cutoff - min y's. Both use only the rows and their bounds, never the node's column
// bounds, so they hold at every node; stored, they prune other nodes whose bounds violate them.

#include <algorithm>
#include <cmath>

#include "mip/branch_and_bound.hpp"
#include "samaya/verify.hpp"

namespace samaya {

namespace {

// A proof may involve at most MipOptions::conflict_density x columns plus this many.
constexpr std::size_t kConflictMinSize = 10;
constexpr std::size_t kMaxConflicts = 1000;
// Coefficients below this, relative to the largest, are dropped by relaxing the side through the
// column's global bounds; a proof needing an infinite bound for that is dropped.
constexpr double kConflictDropTol = 1e-9;
// Relative safety margin on each side, for the round-off of the aggregation.
constexpr double kConflictSafety = 1e-9;
// Farkas ray or dual entries up to this share of the largest are round-off (a row whose logical
// is basic has 0 in theory) and are zeroed, as lp_solver.cpp does for LP certificates.
constexpr double kMultiplierDropTol = 1e-12;
// With an integral objective an improving solution's objective is an integer below the cutoff;
// this slack matches effective_bound's rounding.
constexpr double kIntegralObjectiveSlack = 1e-6;

}  // namespace

double BranchAndBound::conflict_limit() const {
  const double limit = cutoff();
  if (!integral_objective_ || !(limit < kInf)) return limit;
  return offset_ + std::ceil(limit - offset_) - 1.0 + kIntegralObjectiveSlack;
}

bool BranchAndBound::add_conflict(const std::vector<double>& g, double lower, double upper,
                                  bool relative) {
  if (!(lower > -kInf) && !(upper < kInf)) return false;
  double largest = 0.0;
  for (Index j = 0; j < n_; ++j) largest = std::max(largest, std::fabs(g[j]));
  if (!(largest > 0.0) || !std::isfinite(largest)) return false;
  Conflict c;
  c.lower = lower;
  c.upper = upper;
  c.cutoff_relative = relative;
  double magnitude = (std::isfinite(lower) ? std::fabs(lower) : 0.0) +
                     (std::isfinite(upper) ? std::fabs(upper) : 0.0);
  for (Index j = 0; j < n_; ++j) {
    const double a = g[j];
    if (a == 0.0) continue;
    if (std::fabs(a) < kConflictDropTol * largest) {
      // Relax each finite side by the column's contribution over its global bounds. A side needs
      // only the bound that limits it: lower <= g'x loses the largest a x_j, g'x <= upper the
      // smallest.
      const double lo = root_lower_[j];
      const double up = root_upper_[j];
      if (c.lower > -kInf) {
        const double max_term = a > 0.0 ? a * up : a * lo;
        if (!std::isfinite(max_term)) return false;
        c.lower -= max_term;
      }
      if (c.upper < kInf) {
        const double min_term = a > 0.0 ? a * lo : a * up;
        if (!std::isfinite(min_term)) return false;
        c.upper -= min_term;
      }
      continue;
    }
    c.index.push_back(j);
    c.value.push_back(a);
    const double bound = std::max(std::fabs(root_lower_[j]), std::fabs(root_upper_[j]));
    if (std::isfinite(bound)) magnitude += std::fabs(a) * bound;
  }
  const std::size_t limit =
      kConflictMinSize +
      static_cast<std::size_t>(options_.conflict_density * static_cast<double>(n_));
  if (c.index.empty() || c.index.size() > limit) return false;
  const double safety = kConflictSafety * (1.0 + magnitude);
  if (c.lower > -kInf) c.lower -= safety;
  if (c.upper < kInf) c.upper += safety;

  if (!options_.debug_solution.empty()) {
    const std::vector<double>& d = options_.debug_solution;
    double debug_value = offset_;
    for (Index j = 0; j < n_; ++j) debug_value += cost_[j] * d[j];
    // A cutoff proof only speaks about improving solutions.
    if (!relative || debug_value < cutoff()) {
      double activity = 0.0;
      for (std::size_t k = 0; k < c.index.size(); ++k) activity += c.value[k] * d[c.index[k]];
      const double upper_now = relative ? c.upper + conflict_limit() : c.upper;
      const double tol = 1e-6 * (1.0 + std::fabs(activity));
      if (activity < c.lower - tol || activity > upper_now + tol) {
        ++outcome_.debug_conflict_violations;
      }
    }
  }
  c.last_used = outcome_.nodes;
  if (conflicts_.size() >= kMaxConflicts) {
    // Replace the proof that has been idle longest.
    auto oldest = std::min_element(conflicts_.begin(), conflicts_.end(),
                                   [](const Conflict& a, const Conflict& b) {
                                     return a.last_used < b.last_used;
                                   });
    *oldest = std::move(c);
  } else {
    conflicts_.push_back(std::move(c));
  }
  ++outcome_.conflicts_found;
  return true;
}

void BranchAndBound::analyze_infeasible_lp() {
  ++outcome_.conflict_infeasible_lps;
  const bool debug = !options_.debug_solution.empty();
  const std::vector<double>& ray = simplex_->dual_ray();
  if (ray.size() != static_cast<std::size_t>(m_)) {
    if (debug) ++outcome_.debug_farkas_failures;  // An infeasible LP must come with its ray.
    return;
  }
  // A round-off entry on a row whose bound on that side is infinite would make [min y's, max y's]
  // infinite; any multipliers give a valid proof, so zeroing it is safe. The threshold is taken
  // in the scaled space, where the round-off of btran is relative to the largest entry: after
  // unscaling, very different row factors could lift it above the threshold.
  double ray_max = 0.0;
  for (const double r : ray) ray_max = std::max(ray_max, std::fabs(r));
  // Original-space multipliers: the scaled rows are row_i times the original ones.
  std::vector<double> y(static_cast<std::size_t>(m_), 0.0);
  double min_ys = 0.0;
  double max_ys = 0.0;
  for (Index i = 0; i < m_; ++i) {
    if (std::fabs(ray[i]) <= kMultiplierDropTol * ray_max) continue;
    y[i] = ray[i] * scaling_.row[i];
    if (y[i] == 0.0) continue;
    const double lo = y[i] > 0.0 ? model_.row_lower[i] : model_.row_upper[i];
    const double up = y[i] > 0.0 ? model_.row_upper[i] : model_.row_lower[i];
    min_ys += y[i] * lo;  // -inf when the bound is infinite.
    max_ys += y[i] * up;
  }
  if (debug) {
    // Tests: the very multipliers the proof uses must certify, by the independent verifier, that
    // no point within this node's bounds satisfies the rows.
    Model node = model_;
    node.col_lower = lower_;
    node.col_upper = upper_;
    if (!verify_infeasibility(node, y).ok) ++outcome_.debug_farkas_failures;
  }
  std::vector<double> g(static_cast<std::size_t>(n_), 0.0);
  model_.A.multiply_transpose(y, g);
  // Which side is violated at this node decides the proof's direction.
  double min_act = 0.0;
  double max_act = 0.0;
  for (Index j = 0; j < n_; ++j) {
    if (g[j] == 0.0) continue;
    min_act += g[j] > 0.0 ? g[j] * lower_[j] : g[j] * upper_[j];
    max_act += g[j] > 0.0 ? g[j] * upper_[j] : g[j] * lower_[j];
  }
  bool kept = false;
  if (std::isfinite(min_ys) && max_act < min_ys) {
    kept = add_conflict(g, min_ys, kInf, false);
  } else if (std::isfinite(max_ys) && min_act > max_ys) {
    kept = add_conflict(g, -kInf, max_ys, false);
  }
  if (kept) ++outcome_.farkas_proofs;
}

void BranchAndBound::analyze_cutoff_lp() {
  if (!(cutoff() < kInf)) return;
  const std::vector<double>& dual = simplex_->duals();
  // Round-off duals are zeroed as the ray's entries are: any multipliers give a valid proof.
  double dual_max = 0.0;
  for (Index i = 0; i < m_; ++i) dual_max = std::max(dual_max, std::fabs(dual[i]));
  std::vector<double> y(static_cast<std::size_t>(m_), 0.0);
  double min_ys = 0.0;
  for (Index i = 0; i < m_; ++i) {
    if (std::fabs(dual[i]) <= kMultiplierDropTol * dual_max) continue;
    y[i] = dual[i] * scaling_.row[i];
    if (y[i] == 0.0) continue;
    min_ys += y[i] * (y[i] > 0.0 ? model_.row_lower[i] : model_.row_upper[i]);
  }
  if (!std::isfinite(min_ys)) return;
  // d = c - A'y; improving solutions satisfy  d'x <= limit - offset - min y's.
  std::vector<double> d(static_cast<std::size_t>(n_), 0.0);
  model_.A.multiply_transpose(y, d);
  for (Index j = 0; j < n_; ++j) d[j] = cost_[j] - d[j];
  // Keep the proof only if it excludes the LP it came from. One that does not rarely prunes
  // elsewhere, yet every node would pass over it.
  double min_act = 0.0;
  for (Index j = 0; j < n_; ++j) {
    if (d[j] == 0.0) continue;
    min_act += d[j] > 0.0 ? d[j] * lower_[j] : d[j] * upper_[j];
  }
  const double rhs = -offset_ - min_ys;
  if (!(min_act > rhs + conflict_limit())) return;
  add_conflict(d, -kInf, rhs, true);
}

bool BranchAndBound::propagate_conflicts(std::vector<BoundChange>* record) {
  std::vector<Index> changed;
  const double limit = conflict_limit();
  for (Conflict& c : conflicts_) {
    const double rl = c.lower;
    double ru = c.upper;
    if (c.cutoff_relative) {
      if (!(limit < kInf)) continue;
      ru += limit;
    }
    double min_act = 0.0;
    double max_act = 0.0;
    int min_inf = 0;
    int max_inf = 0;
    for (std::size_t k = 0; k < c.index.size(); ++k) {
      const Index j = c.index[k];
      const double a = c.value[k];
      const double lo = a > 0.0 ? lower_[j] : upper_[j];
      const double up = a > 0.0 ? upper_[j] : lower_[j];
      if (std::isfinite(lo)) {
        min_act += a * lo;
      } else {
        ++min_inf;
      }
      if (std::isfinite(up)) {
        max_act += a * up;
      } else {
        ++max_inf;
      }
    }
    if ((min_inf == 0 && min_act > ru) || (max_inf == 0 && max_act < rl)) {
      c.last_used = outcome_.nodes;
      return false;
    }
    // Tighten the integer columns (as propagate() does for the rows).
    bool used = false;
    for (std::size_t k = 0; k < c.index.size(); ++k) {
      const Index j = c.index[k];
      if (model_.col_type[j] != VarType::kInteger) continue;
      const double a = c.value[k];
      double new_lo = lower_[j];
      double new_up = upper_[j];
      if (ru < kInf && min_inf == 0) {
        const double residual = min_act - (a > 0.0 ? a * lower_[j] : a * upper_[j]);
        const double b = (ru - residual) / a;
        const double t = options_.integrality_tol * (1.0 + std::fabs(b));
        if (a > 0.0) {
          new_up = std::min(new_up, std::floor(b + t));
        } else {
          new_lo = std::max(new_lo, std::ceil(b - t));
        }
      }
      if (rl > -kInf && max_inf == 0) {
        const double residual = max_act - (a > 0.0 ? a * upper_[j] : a * lower_[j]);
        const double b = (rl - residual) / a;
        const double t = options_.integrality_tol * (1.0 + std::fabs(b));
        if (a > 0.0) {
          new_lo = std::max(new_lo, std::ceil(b - t));
        } else {
          new_up = std::min(new_up, std::floor(b + t));
        }
      }
      if (new_lo > new_up) {
        c.last_used = outcome_.nodes;
        return false;
      }
      if (new_lo > lower_[j] || new_up < upper_[j]) {
        set_bound(j, new_lo, new_up);
        if (record != nullptr) record->push_back({j, new_lo, new_up});
        changed.push_back(j);
        used = true;
      }
    }
    if (used) c.last_used = outcome_.nodes;
  }
  // The tightened columns may let the rows tighten more.
  return changed.empty() || propagate(std::move(changed), record);
}

}  // namespace samaya
