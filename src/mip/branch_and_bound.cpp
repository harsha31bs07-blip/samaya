#include "mip/branch_and_bound.hpp"

#include "mip/cuts.hpp"
#include "mip/search_constants.hpp"
#include "presolve/presolve.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <optional>
#include <utility>

namespace samaya {

namespace {

constexpr double kScoreFloor = 1e-6;
// The search restarts after the root when at least this share of the integer columns was fixed
// there (by propagation, strong branching and reduced costs): the fixings are presolved away and
// the smaller model gets its own root cuts and heuristics.
constexpr double kRestartFraction = 0.2;
// Reduced-cost fixing ignores reduced costs below this and allows this relative slack on the
// cutoff for the LP's tolerances.
constexpr double kReducedCostTol = 1e-7;
constexpr double kReducedCostSlack = 1e-6;
constexpr double kStrongIterationShare = 0.5;
constexpr double kStrongIterationOffset = 20000.0;
// A solution replaces the incumbent only if it is better by this much (relative).
constexpr double kMinImprovement = 1e-9;
// Releasing the stored nodes when the search ends takes about 1 us per node (0.85 us measured on
// gen-ip002 with 175k open nodes). The time limit reserves this, with a margin, so a search with
// millions of open nodes still finishes within the limit.
constexpr double kReleaseSecondsPerNode = 2e-6;
// Rounding followed by an LP over the continuous columns runs at the root and every this many
// nodes (with search::kRoundAndSolveIterations).
constexpr long long kRoundAndSolveFrequency = 50;
// Primal and dual tolerance of the re-solve that repairs a rounded solution.
constexpr double kTightPrimalTol = 1e-9;
// Stored warm-start bases may use at most this many bytes; beyond it nodes start from the root
// basis.
constexpr double kBasisMemoryBudget = 256.0 * 1024 * 1024;
// A node re-solves after strong branching tightens it at most this many times.
constexpr int kMaxTightenRounds = 5;
// Propagation stops after this many row visits per call and ignores huge derived bounds.
constexpr double kMaxDerivedBound = 1e12;
// A derived bound on a continuous column is applied only if it moves the bound by this fraction
// of min(domain width, |bound|), at least kMinContinuousStep, or moves it across zero (SCIP's
// SCIPsetIsLbBetter with numerics/boundstreps = 0.05): without it, rows that chain continuous
// columns shrink their domains by ever smaller steps until the budget runs out.
constexpr double kContinuousBoundStep = 0.05;
constexpr double kMinContinuousStep = 1e-3;
// Relative round-off of a row's activity sum, relaxing continuous bounds derived from it.
constexpr double kActivityRoundoff = 1e-12;
// Tests: debug_solution is inside a box if within this relative distance of its bounds (the
// reference solutions satisfy the rows to about 1e-9).
constexpr double kDebugBoxTolerance = 1e-7;
// Cut selection and the root cut loop.
constexpr double kMinCutEfficacy = 1e-5;
constexpr double kMaxCutParallelism = 0.995;
constexpr double kCutStallImprovement = 1e-4;  // Relative bound gain that counts as progress.
constexpr int kCutStallRounds = 3;
// The cut loop may take this multiple of the root LP time (at least kCutLoopMinSeconds, at most
// kCutLoopTimeFraction of the time limit); cuts are dropped if they improve the root bound by
// less than kCutMinTotalGain relative.
constexpr double kCutLoopRootFactor = 5.0;
constexpr double kCutLoopMinSeconds = 1.0;
constexpr double kCutLoopTimeFraction = 0.1;
constexpr double kCutMinTotalGain = 1e-6;
// Cuts in the tree (fresh separation every MipOptions::tree_separation_frequency depths, the pool
// at every node). A tree cut must be violated this much (efficacy, as SCIP's
// separating/minefficacy), a node adds at most kMaxTreeCutsPerNode, and tree cut rounds may take
// at most kTreeCutTimeFraction of the solve time.
constexpr double kMinTreeCutEfficacy = 1e-4;
constexpr int kMaxTreeCutsPerNode = 20;
constexpr double kTreeCutTimeFraction = 0.1;
constexpr std::size_t kMaxPoolCuts = 5000;
// A node's tree cuts stay in the LP only if they close this share of the node's gap.
constexpr double kTreeCutMinGapShare = 0.01;
constexpr long long kMaxTreeCutInterval = 10000;

}  // namespace

BranchAndBound::BranchAndBound(const Model& model, const MipOptions& options, const Logger& log)
    : model_(model),
      options_(options),
      log_(log),
      m_(model.num_rows()),
      n_(model.num_cols()),
      original_rows_(model.num_rows()),
      sense_(model.sense == ObjSense::kMaximize ? -1.0 : 1.0),
      offset_(sense_ * model.obj_offset) {
  cost_.resize(static_cast<std::size_t>(n_));
  integral_objective_ = true;
  bool any_cost = false;
  for (Index j = 0; j < n_; ++j) {
    cost_[j] = sense_ * model.obj[j];
    if (model.col_type[j] == VarType::kInteger) {
      integers_.push_back(j);
      if (cost_[j] != std::round(cost_[j])) integral_objective_ = false;
      any_cost = any_cost || cost_[j] != 0.0;
    } else if (cost_[j] != 0.0) {
      integral_objective_ = false;
    }
  }
  integral_objective_ = integral_objective_ && any_cost;

  // Locks: rows that may become violated when a column moves down (up).
  down_locks_.assign(static_cast<std::size_t>(n_), 0);
  up_locks_.assign(static_cast<std::size_t>(n_), 0);
  const auto start = model.A.col_start();
  const auto index = model.A.row_index();
  const auto val = model.A.values();
  for (Index j = 0; j < n_; ++j) {
    for (NnzIndex p = start[j]; p < start[j + 1]; ++p) {
      const bool has_lower = model.row_lower[index[p]] > -kInf;
      const bool has_upper = model.row_upper[index[p]] < kInf;
      if (val[p] > 0.0) {
        down_locks_[j] += has_lower;
        up_locks_[j] += has_upper;
      } else {
        down_locks_[j] += has_upper;
        up_locks_[j] += has_lower;
      }
    }
  }

  root_lower_ = model.col_lower;
  root_upper_ = model.col_upper;
  lower_ = root_lower_;
  upper_ = root_upper_;
  is_touched_.assign(static_cast<std::size_t>(n_), 0);
  x_.assign(static_cast<std::size_t>(n_), 0.0);
  for (int d = 0; d < 2; ++d) {
    pc_sum_[d].assign(static_cast<std::size_t>(n_), 0.0);
    pc_count_[d].assign(static_cast<std::size_t>(n_), 0);
  }
  build_relaxation();
}

// (Re)builds the scaled LP, its simplex and the row-wise matrix from model_ and the current
// column bounds; called at construction and whenever cuts change the rows.
void BranchAndBound::build_relaxation() {
  m_ = model_.num_rows();
  At_ = model_.A.transpose();
  row_mark_.assign(static_cast<std::size_t>(m_), 0);
  scaling_ = compute_scaling(model_.A);
  lp_ = LpProblem{};
  lp_.m = m_;
  lp_.n = n_;
  lp_.A = scaling_.apply(model_.A);
  lp_.At = lp_.A.transpose();
  const auto nt = static_cast<std::size_t>(n_ + m_);
  lp_.cost.assign(nt, 0.0);
  lp_.lower.resize(nt);
  lp_.upper.resize(nt);
  for (Index j = 0; j < n_; ++j) {
    lp_.cost[j] = cost_[j] * scaling_.col[j];
    lp_.lower[j] = lower_[j] / scaling_.col[j];
    lp_.upper[j] = upper_[j] / scaling_.col[j];
  }
  for (Index i = 0; i < m_; ++i) {
    lp_.lower[n_ + i] = model_.row_lower[i] * scaling_.row[i];
    lp_.upper[n_ + i] = model_.row_upper[i] * scaling_.row[i];
  }
  simplex_ = std::make_unique<Simplex>(lp_, SimplexOptions{}, lp_log_);
}

BranchAndBound::~BranchAndBound() = default;

// ---------------------------------------------------------------------------------------------
// Relaxation

void BranchAndBound::set_bound(Index j, double lower, double upper) {
  if (logging_undo_) undo_log_.push_back({j, lower_[j], upper_[j]});
  lower_[j] = lower;
  upper_[j] = upper;
  lp_.lower[j] = lower / scaling_.col[j];
  lp_.upper[j] = upper / scaling_.col[j];
  if (!is_touched_[j] && (lower != root_lower_[j] || upper != root_upper_[j])) {
    is_touched_[j] = 1;
    touched_.push_back(j);
  }
}

void BranchAndBound::apply_node_bounds(const Node& node) {
  for (const Index j : touched_) {
    lower_[j] = root_lower_[j];
    upper_[j] = root_upper_[j];
    lp_.lower[j] = lower_[j] / scaling_.col[j];
    lp_.upper[j] = upper_[j] / scaling_.col[j];
    is_touched_[j] = 0;
  }
  touched_.clear();
  std::vector<const PathSegment*> segments;
  for (const PathSegment* seg = node.path.get(); seg != nullptr; seg = seg->parent.get()) {
    segments.push_back(seg);
  }
  for (auto it = segments.rbegin(); it != segments.rend(); ++it) {
    for (const BoundChange& c : (*it)->changes) set_bound(c.col, c.lower, c.upper);
  }
  for (const BoundChange& c : node.local) set_bound(c.col, c.lower, c.upper);
}

SimplexStatus BranchAndBound::solve_relaxation(const std::vector<VarStatus>* start,
                                               long long iteration_limit) {
  simplex_->set_iteration_limit(iteration_limit);
  simplex_->set_time_limit(std::max(0.0, remaining_time()));
  // A basis stored before tree cuts were appended lacks their rows: their logicals are basic.
  std::vector<VarStatus> padded;
  const auto total = static_cast<std::size_t>(n_ + m_);
  if (start != nullptr && !start->empty() && start->size() < total) {
    padded = *start;
    padded.resize(total, VarStatus::kBasic);
    start = &padded;
  }
  lp_status_ = start != nullptr && !start->empty() ? simplex_->solve(*start) : simplex_->solve();
  outcome_.lp_iterations += simplex_->iterations();
  const std::vector<double>& v = simplex_->values();
  for (Index j = 0; j < n_; ++j) x_[j] = v[j] * scaling_.col[j];
  return lp_status_;
}

double BranchAndBound::relaxation_objective() const {
  double obj = offset_;
  for (Index j = 0; j < n_; ++j) obj += cost_[j] * x_[j];
  return obj;
}

std::vector<VarStatus> BranchAndBound::current_basis() const { return simplex_->status(); }

// ---------------------------------------------------------------------------------------------
// Bounds, pseudocosts and limits

double BranchAndBound::cutoff() const {
  if (incumbent_value_ == kInf) return kInf;
  return incumbent_value_ -
         std::max(options_.abs_gap, options_.rel_gap * std::fabs(incumbent_value_));
}

double BranchAndBound::effective_bound(double objective) const {
  if (!integral_objective_ || !std::isfinite(objective)) return objective;
  return std::ceil(objective - offset_ - 1e-6) + offset_;
}

double BranchAndBound::best_bound(double extra) const {
  double bound = extra;
  if (!open_.empty()) bound = std::min(bound, open_.begin()->first);
  for (const Node& node : dive_stack_) bound = std::min(bound, node.bound);
  return effective_bound(bound);
}

double BranchAndBound::remaining_time() const {
  const double stored =
      static_cast<double>(open_.size() + dive_stack_.size() + shared_stored_);
  return options_.time_limit - timer_.seconds() - kReleaseSecondsPerNode * stored;
}

bool BranchAndBound::time_up() const { return remaining_time() <= 0.0; }

// Strong branching may use at most kStrongIterationShare of the node LP iterations plus
// kStrongIterationOffset (as SCIP's reliability branching does); beyond that unreliable columns
// are scored by their (average) pseudocosts. Without the cap strong branching took 88% of the
// LP iterations on beasleyC3 and 77% on mc11.
bool BranchAndBound::strong_budget_left() const {
  const double node_iterations = static_cast<double>(
      outcome_.lp_iterations - outcome_.strong_branching_iterations -
      outcome_.heuristic_lp_iterations);
  return static_cast<double>(outcome_.strong_branching_iterations) <
         kStrongIterationShare * node_iterations + kStrongIterationOffset;
}

// The child to plunge into: the one with the smaller pseudocost estimate of the bound
// degradation (children[0] is the down branch, children[1] the up branch); ties go to the side
// the LP value is closer to.
std::size_t BranchAndBound::plunge_child(const std::vector<Node>& children) const {
  const Index j = children[0].branch_col;
  const double down = pseudocost(j, false) * children[0].branch_distance;
  const double up = pseudocost(j, true) * children[1].branch_distance;
  if (down < up) return 0;
  if (up < down) return 1;
  return children[0].branch_distance >= 0.5 ? 1 : 0;
}

double BranchAndBound::pseudocost(Index j, bool up) const {
  const int d = up ? 1 : 0;
  if (pc_count_[d][j] > 0) return pc_sum_[d][j] / pc_count_[d][j];
  if (pc_total_count_[d] > 0) return pc_total_sum_[d] / static_cast<double>(pc_total_count_[d]);
  return 1.0;
}

void BranchAndBound::record_pseudocost(Index j, bool up, double gain, double distance) {
  if (distance < 1e-9 || !std::isfinite(gain)) return;
  const int d = up ? 1 : 0;
  const double unit = std::max(gain, 0.0) / distance;
  pc_sum_[d][j] += unit;
  ++pc_count_[d][j];
  pc_total_sum_[d] += unit;
  ++pc_total_count_[d];
}

bool BranchAndBound::is_fractional(double v) const {
  return std::fabs(v - std::round(v)) > options_.integrality_tol;
}

// ---------------------------------------------------------------------------------------------
// Propagation

namespace {

// Whether a derived bound on a continuous column is worth applying (kContinuousBoundStep).
bool improves_lower(double new_lower, double lower, double upper) {
  if (!(new_lower > lower)) return false;
  if (!std::isfinite(lower) || (lower < 0.0 && new_lower >= 0.0)) return true;
  const double scale = std::max(std::min(upper - lower, std::fabs(lower)), kMinContinuousStep);
  return new_lower > lower + kContinuousBoundStep * scale;
}

bool improves_upper(double new_upper, double lower, double upper) {
  return improves_lower(-new_upper, -upper, -lower);
}

}  // namespace

bool BranchAndBound::propagate(std::vector<Index> changed, std::vector<BoundChange>* record) {
  // Tests: whether debug_solution lies in the current box, which propagation must keep.
  const bool debug = !options_.debug_solution.empty() && contains_debug_solution();
  const bool feasible = propagate_rows(std::move(changed), record);
  if (debug && (!feasible || !contains_debug_solution())) ++outcome_.debug_propagation_violations;
  return feasible;
}

bool BranchAndBound::contains_debug_solution() const {
  const std::vector<double>& d = options_.debug_solution;
  for (Index j = 0; j < n_; ++j) {
    const double tol = kDebugBoxTolerance * (1.0 + std::fabs(d[j]));
    if (d[j] < lower_[j] - tol || d[j] > upper_[j] + tol) return false;
  }
  return true;
}

bool BranchAndBound::propagate_rows(std::vector<Index> changed,
                                    std::vector<BoundChange>* record) {
  const auto cstart = model_.A.col_start();
  const auto cindex = model_.A.row_index();
  const auto rstart = At_.col_start();
  const auto rindex = At_.row_index();
  const auto rval = At_.values();
  std::vector<Index> queue;
  const auto enqueue_col = [&](Index j) {
    for (NnzIndex p = cstart[j]; p < cstart[j + 1]; ++p) {
      const Index i = cindex[p];
      if (!row_mark_[i]) {
        row_mark_[i] = 1;
        queue.push_back(i);
      }
    }
  };
  for (const Index j : changed) enqueue_col(j);

  const std::size_t budget = 20 * static_cast<std::size_t>(m_) + 1000;
  bool feasible = true;
  std::size_t head = 0;
  for (; head < queue.size() && head < budget && feasible; ++head) {
    const Index i = queue[head];
    row_mark_[i] = 0;
    const double rl = model_.row_lower[i];
    const double ru = model_.row_upper[i];
    double min_fin = 0.0;
    double max_fin = 0.0;
    int min_inf = 0;
    int max_inf = 0;
    double max_range = 0.0;  // Largest |a_j| (u_j - l_j) over the columns.
    double magnitude = 0.0;  // Sum of the finite |a_j l_j| and |a_j u_j|, for round-off.
    for (NnzIndex p = rstart[i]; p < rstart[i + 1]; ++p) {
      const Index j = rindex[p];
      const double a = rval[p];
      const double lo = a > 0.0 ? lower_[j] : upper_[j];
      const double up = a > 0.0 ? upper_[j] : lower_[j];
      max_range = std::max(max_range, std::fabs(a) * (upper_[j] - lower_[j]));
      if (std::isfinite(lo)) {
        min_fin += a * lo;
        magnitude += std::fabs(a * lo);
      } else {
        ++min_inf;
      }
      if (std::isfinite(up)) {
        max_fin += a * up;
        magnitude += std::fabs(a * up);
      } else {
        ++max_inf;
      }
    }
    const double tol_u = options_.feasibility_tol * (1.0 + std::fabs(ru));
    const double tol_l = options_.feasibility_tol * (1.0 + std::fabs(rl));
    if ((min_inf == 0 && min_fin > ru + tol_u) || (max_inf == 0 && max_fin < rl - tol_l)) {
      feasible = false;
      break;
    }
    // A side tightens column j only if its slack is below |a_j| (u_j - l_j): skip the row when
    // both sides have more slack than the largest such range.
    const bool upper_slack = ru == kInf || (min_inf == 0 && ru - min_fin >= max_range);
    const bool lower_slack = rl == -kInf || (max_inf == 0 && max_fin - rl >= max_range);
    if (upper_slack && lower_slack) continue;
    for (NnzIndex p = rstart[i]; p < rstart[i + 1]; ++p) {
      const Index j = rindex[p];
      const bool integer = model_.col_type[j] == VarType::kInteger;
      if (!integer && !options_.propagate_continuous) continue;
      const double a = rval[p];
      double new_lo = lower_[j];
      double new_up = upper_[j];
      // A continuous bound is relaxed by the bound tolerance, by the row tolerance an accepted
      // solution may use (divided by |a|, as it moves x_j that much) and by the round-off of the
      // activity sum, so that it never cuts off a solution the tolerances accept.
      const auto slack = [&](double b, double row_tol) {
        return options_.feasibility_tol * (1.0 + std::fabs(b)) +
               (row_tol + kActivityRoundoff * magnitude) / std::fabs(a);
      };
      // a x_j <= ru - (minimum activity of the other columns), and the same for rl.
      if (ru < kInf) {
        const double contrib = a > 0.0 ? a * lower_[j] : a * upper_[j];
        double residual = kInf;
        if (!std::isfinite(contrib)) {
          if (min_inf == 1) residual = min_fin;
        } else if (min_inf == 0) {
          residual = min_fin - contrib;
        }
        const double b = (ru - residual) / a;
        if (std::isfinite(residual) && std::fabs(b) < kMaxDerivedBound) {
          if (integer) {
            const double t = options_.integrality_tol * (1.0 + std::fabs(b));
            if (a > 0.0) {
              new_up = std::min(new_up, std::floor(b + t));
            } else {
              new_lo = std::max(new_lo, std::ceil(b - t));
            }
          } else if (a > 0.0) {
            new_up = std::min(new_up, b + slack(b, tol_u));
          } else {
            new_lo = std::max(new_lo, b - slack(b, tol_u));
          }
        }
      }
      if (rl > -kInf) {
        const double contrib = a > 0.0 ? a * upper_[j] : a * lower_[j];
        double residual = -kInf;
        if (!std::isfinite(contrib)) {
          if (max_inf == 1) residual = max_fin;
        } else if (max_inf == 0) {
          residual = max_fin - contrib;
        }
        const double b = (rl - residual) / a;
        if (std::isfinite(residual) && std::fabs(b) < kMaxDerivedBound) {
          if (integer) {
            const double t = options_.integrality_tol * (1.0 + std::fabs(b));
            if (a > 0.0) {
              new_lo = std::max(new_lo, std::ceil(b - t));
            } else {
              new_up = std::min(new_up, std::floor(b + t));
            }
          } else if (a > 0.0) {
            new_lo = std::max(new_lo, b - slack(b, tol_l));
          } else {
            new_up = std::min(new_up, b + slack(b, tol_l));
          }
        }
      }
      if (!integer) {
        if (!improves_lower(new_lo, lower_[j], upper_[j])) new_lo = lower_[j];
        if (!improves_upper(new_up, lower_[j], upper_[j])) new_up = upper_[j];
        if (new_lo <= new_up && (new_lo > lower_[j] || new_up < upper_[j])) {
          ++outcome_.continuous_tightened;
        }
      }
      if (new_lo > new_up) {
        feasible = false;
        break;
      }
      if (new_lo > lower_[j] || new_up < upper_[j]) {
        set_bound(j, new_lo, new_up);
        if (record != nullptr) record->push_back({j, new_lo, new_up});
        enqueue_col(j);
      }
    }
  }
  for (std::size_t k = head; k < queue.size(); ++k) row_mark_[queue[k]] = 0;
  return feasible;
}

// ---------------------------------------------------------------------------------------------
// Solutions

bool BranchAndBound::improves(double value) const {
  if (incumbent_value_ == kInf) return value < kInf;
  return value < incumbent_value_ - kMinImprovement * (1.0 + std::fabs(incumbent_value_));
}

bool BranchAndBound::try_solution(std::vector<double> x) {
  for (const Index j : integers_) {
    if (is_fractional(x[j])) return false;
    x[j] = std::round(x[j]);
  }
  const double tol = options_.feasibility_tol;
  // The root bounds (after propagation and probing) hold for every feasible point, and the
  // tightened rows are equivalent to the model's only inside them.
  for (Index j = 0; j < n_; ++j) {
    if (x[j] < root_lower_[j] - tol * (1.0 + std::fabs(root_lower_[j])) ||
        x[j] > root_upper_[j] + tol * (1.0 + std::fabs(root_upper_[j]))) {
      return false;
    }
  }
  std::vector<double> activity(static_cast<std::size_t>(m_), 0.0);
  model_.A.multiply(x, activity);
  // Only the model's own rows: a solution is never rejected because of a cut's round-off.
  for (Index i = 0; i < original_rows_; ++i) {
    const double rl = model_.row_lower[i];
    const double ru = model_.row_upper[i];
    if (activity[i] < rl - tol * (1.0 + std::fabs(rl)) ||
        activity[i] > ru + tol * (1.0 + std::fabs(ru))) {
      return false;
    }
  }
  double value = offset_;
  for (Index j = 0; j < n_; ++j) value += cost_[j] * x[j];
  if (!improves(value)) return false;
  incumbent_value_ = value;
  incumbent_ = std::move(x);
  if (shared_ != nullptr) publish_incumbent();
  log_.log(2, "mip: new incumbent %.12g after %lld nodes, %.2f s", sense_ * value,
           outcome_.nodes, timer_.seconds());
  // Prune open nodes the incumbent cuts off.
  const double limit = cutoff();
  for (auto it = open_.begin(); it != open_.end();) {
    const double b = effective_bound(it->first);
    if (b >= limit) {
      pruned_bound_ = std::min(pruned_bound_, b);
      it = open_.erase(it);
    } else {
      ++it;
    }
  }
  return true;
}

void BranchAndBound::simple_rounding(const std::vector<double>& x) {
  std::vector<double> y = x;
  for (const Index j : integers_) {
    const double v = y[j];
    if (!is_fractional(v)) {
      y[j] = std::round(v);
      continue;
    }
    if (down_locks_[j] == 0) {
      y[j] = std::floor(v);
    } else if (up_locks_[j] == 0) {
      y[j] = std::ceil(v);
    } else {
      return;
    }
    y[j] = std::clamp(y[j], lower_[j], upper_[j]);
  }
  if (try_solution(std::move(y))) ++outcome_.heuristic_solutions;
}

bool BranchAndBound::round_and_solve(const std::vector<double>& x,
                                     const std::vector<VarStatus>& basis,
                                     long long iteration_limit) {
  std::vector<std::pair<double, double>> saved;
  saved.reserve(integers_.size());
  for (const Index j : integers_) {
    saved.emplace_back(lower_[j], upper_[j]);
    const double v = std::clamp(std::round(x[j]), lower_[j], upper_[j]);
    set_bound(j, v, v);
  }
  const SimplexStatus status = solve_relaxation(&basis, iteration_limit);
  bool found = status == SimplexStatus::kOptimal && try_solution(x_);
  if (status == SimplexStatus::kOptimal && !found && relaxation_objective() < incumbent_value_) {
    // The scaled tolerance (1e-7) can leave a row violated by more than the acceptance
    // tolerance once unscaled; re-solve this LP from its optimal basis with tight tolerances.
    SimplexOptions tight;
    tight.primal_tol = kTightPrimalTol;
    tight.dual_tol = kTightPrimalTol;
    tight.perturb = false;
    Simplex exact(lp_, tight, lp_log_);
    exact.set_iteration_limit(iteration_limit);
    exact.set_time_limit(std::max(0.0, remaining_time()));
    if (exact.solve(simplex_->status()) == SimplexStatus::kOptimal) {
      outcome_.lp_iterations += exact.iterations();
      std::vector<double> y(static_cast<std::size_t>(n_));
      for (Index j = 0; j < n_; ++j) y[j] = exact.values()[j] * scaling_.col[j];
      found = try_solution(std::move(y));
    }
  }
  if (found) ++outcome_.heuristic_solutions;
  for (std::size_t k = 0; k < integers_.size(); ++k) {
    set_bound(integers_[k], saved[k].first, saved[k].second);
  }
  return found;
}

// ---------------------------------------------------------------------------------------------
// Cuts

void BranchAndBound::remove_rows(const std::vector<Index>& rows) {
  if (rows.empty()) return;
  std::vector<Index> new_index(static_cast<std::size_t>(m_), 0);
  for (const Index i : rows) new_index[i] = -1;
  Index next = 0;
  for (Index i = 0; i < m_; ++i) {
    if (new_index[i] < 0) continue;
    new_index[i] = next;
    model_.row_lower[next] = model_.row_lower[i];
    model_.row_upper[next] = model_.row_upper[i];
    ++next;
  }
  model_.row_lower.resize(static_cast<std::size_t>(next));
  model_.row_upper.resize(static_cast<std::size_t>(next));
  model_.row_names.clear();
  std::vector<Triplet> t;
  const auto start = model_.A.col_start();
  const auto index = model_.A.row_index();
  const auto val = model_.A.values();
  for (Index j = 0; j < n_; ++j) {
    for (NnzIndex p = start[j]; p < start[j + 1]; ++p) {
      if (new_index[index[p]] >= 0) t.push_back({new_index[index[p]], j, val[p]});
    }
  }
  model_.A = SparseMatrix::from_triplets(next, n_, std::move(t));
}

void BranchAndBound::pool_rows(const std::vector<Index>& rows) {
  const auto start = At_.col_start();
  const auto index = At_.row_index();
  const auto val = At_.values();
  for (const Index i : rows) {
    if (i < original_rows_ || cut_pool_.size() >= kMaxPoolCuts) continue;
    Cut cut;
    for (NnzIndex p = start[i]; p < start[i + 1]; ++p) {
      cut.index.push_back(index[p]);
      cut.value.push_back(val[p]);
    }
    cut.lower = model_.row_lower[i];
    cut_pool_.push_back(std::move(cut));
  }
}

void BranchAndBound::separate_gomory(const CutContext& ctx, std::vector<Cut>& candidates) {
  const std::vector<double>& x = ctx.x;
  // Gomory mixed-integer cuts from the rows of fractional basic integer columns, most fractional
  // first. The scaled tableau row is mapped to original units: structural j gets
  // a_j col_k / col_j, the logical of row i a_{n+i} col_k row_i.
  const std::vector<Index>& basic = simplex_->basic();
  const std::vector<VarStatus> status = simplex_->status();
  std::vector<std::pair<double, Index>> rows;
  for (Index r = 0; r < m_; ++r) {
    const Index k = basic[r];
    if (k >= n_ || model_.col_type[k] != VarType::kInteger || !is_fractional(x[k])) continue;
    rows.emplace_back(std::fabs(x[k] - std::floor(x[k]) - 0.5), r);
  }
  std::sort(rows.begin(), rows.end());
  if (rows.size() > 2 * static_cast<std::size_t>(options_.max_cuts_per_round)) {
    rows.resize(2 * static_cast<std::size_t>(options_.max_cuts_per_round));
  }
  std::vector<double> scaled_row;
  std::vector<double> row(static_cast<std::size_t>(n_ + m_));
  for (const auto& [score, r] : rows) {
    const Index k = basic[r];
    simplex_->tableau_row(r, scaled_row);
    const double ck = scaling_.col[k];
    for (Index j = 0; j < n_; ++j) row[j] = scaled_row[j] * ck / scaling_.col[j];
    for (Index i = 0; i < m_; ++i) row[n_ + i] = scaled_row[n_ + i] * ck * scaling_.row[i];
    Cut cut;
    if (gomory_mixed_integer_cut(ctx, k, row, status, cut) && finalize_cut(ctx, cut)) {
      candidates.push_back(std::move(cut));
    }
  }
}

int BranchAndBound::separate_and_add_cuts(bool tree, bool fresh, int tree_max_cuts) {
  const std::vector<double> x = x_;
  const CutContext ctx{model_, At_, tree ? root_lower_ : lower_, tree ? root_upper_ : upper_, x,
                       original_rows_};
  std::vector<Cut> candidates;
  const double min_efficacy = tree ? kMinTreeCutEfficacy : kMinCutEfficacy;
  const int max_cuts = tree ? tree_max_cuts : options_.max_cuts_per_round;
  if (tree) {
    for (std::size_t k = 0; k < cut_pool_.size(); ++k) {
      const Cut& c = cut_pool_[k];
      double activity = 0.0;
      double norm = 0.0;
      for (std::size_t q = 0; q < c.index.size(); ++q) {
        activity += c.value[q] * x[c.index[q]];
        norm += c.value[q] * c.value[q];
      }
      const double efficacy = norm > 0.0 ? (c.lower - activity) / std::sqrt(norm) : 0.0;
      if (efficacy < min_efficacy) continue;
      Cut candidate = c;
      candidate.efficacy = efficacy;
      candidate.pool_index = static_cast<int>(k);
      candidates.push_back(std::move(candidate));
    }
  }
  if (!tree) separate_gomory(ctx, candidates);
  if (fresh) {
    separate_mir(ctx, candidates);
    aggregation_start_ = separate_aggregated_mir(ctx, candidates, aggregation_start_);
    separate_knapsack_covers(ctx, candidates);
    separate_cliques(ctx, candidates);
  }

  // Select the most efficacious cuts, skipping near-parallel ones.
  for (Cut& c : candidates) {
    std::vector<std::size_t> order(c.index.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
              [&](std::size_t a, std::size_t b) { return c.index[a] < c.index[b]; });
    Cut sorted;
    for (const std::size_t k : order) {
      sorted.index.push_back(c.index[k]);
      sorted.value.push_back(c.value[k]);
    }
    c.index = std::move(sorted.index);
    c.value = std::move(sorted.value);
  }
  std::sort(candidates.begin(), candidates.end(),
            [](const Cut& a, const Cut& b) { return a.efficacy > b.efficacy; });
  const auto norm = [](const Cut& c) {
    double s = 0.0;
    for (const double v : c.value) s += v * v;
    return std::sqrt(s);
  };
  const auto cosine = [&](const Cut& a, const Cut& b) {
    double dot = 0.0;
    std::size_t p = 0;
    std::size_t q = 0;
    while (p < a.index.size() && q < b.index.size()) {
      if (a.index[p] == b.index[q]) {
        dot += a.value[p++] * b.value[q++];
      } else if (a.index[p] < b.index[q]) {
        ++p;
      } else {
        ++q;
      }
    }
    return dot / (norm(a) * norm(b));
  };
  std::vector<const Cut*> chosen;
  for (const Cut& c : candidates) {
    if (static_cast<int>(chosen.size()) >= max_cuts) break;
    if (c.efficacy < min_efficacy) break;
    if (std::any_of(chosen.begin(), chosen.end(),
                    [&](const Cut* o) { return cosine(c, *o) > kMaxCutParallelism; })) {
      continue;
    }
    chosen.push_back(&c);
  }
  if (chosen.empty()) return 0;

  if (!options_.debug_solution.empty()) {
    for (const Cut* c : chosen) {
      double activity = 0.0;
      for (std::size_t k = 0; k < c->index.size(); ++k) {
        activity += c->value[k] * options_.debug_solution[c->index[k]];
      }
      if (activity < c->lower - 1e-6 * (1.0 + std::fabs(c->lower))) {
        ++outcome_.debug_cut_violations;
        log_.log(1, "mip: cut separates the debug solution (%.9g < %.9g)", activity, c->lower);
      }
    }
  }

  // Append the cuts as rows  lower <= pi'x <= +inf.
  std::vector<Triplet> t;
  const auto start = model_.A.col_start();
  const auto index = model_.A.row_index();
  const auto val = model_.A.values();
  for (Index j = 0; j < n_; ++j) {
    for (NnzIndex p = start[j]; p < start[j + 1]; ++p) t.push_back({index[p], j, val[p]});
  }
  Index row_id = m_;
  for (const Cut* c : chosen) {
    for (std::size_t k = 0; k < c->index.size(); ++k) t.push_back({row_id, c->index[k], c->value[k]});
    model_.row_lower.push_back(c->lower);
    model_.row_upper.push_back(kInf);
    ++row_id;
  }
  model_.row_names.clear();
  model_.A = SparseMatrix::from_triplets(row_id, n_, std::move(t));
  // Pool cuts now in the LP leave the pool.
  std::vector<char> taken(cut_pool_.size(), 0);
  bool any_taken = false;
  for (const Cut* c : chosen) {
    if (c->pool_index >= 0) {
      taken[static_cast<std::size_t>(c->pool_index)] = 1;
      any_taken = true;
    }
  }
  if (any_taken) {
    std::size_t next = 0;
    for (std::size_t k = 0; k < cut_pool_.size(); ++k) {
      if (!taken[k]) cut_pool_[next++] = std::move(cut_pool_[k]);
    }
    cut_pool_.resize(next);
  }
  return static_cast<int>(chosen.size());
}

int BranchAndBound::tree_cut_round(int depth) {
  if (!options_.cuts || !options_.tree_cuts || shared_ != nullptr || depth == 0) return 0;
  const int room =
      static_cast<int>(options_.tree_cut_row_fraction * static_cast<double>(original_rows_)) -
      outcome_.tree_cuts;
  if (room <= 0) return 0;
  if (tree_cut_seconds_ > kTreeCutTimeFraction * timer_.seconds()) return 0;
  if (outcome_.nodes < next_tree_cut_node_) return 0;
  // Fresh separation also waits after finding nothing (on mas76 it finds nothing at hundreds of
  // thousands of nodes); the pool check is a few dot products and runs at every node.
  const bool fresh = options_.tree_separation_frequency > 0 &&
                     depth % options_.tree_separation_frequency == 0 &&
                     outcome_.nodes >= next_fresh_node_;
  if (!fresh && cut_pool_.empty()) return 0;
  const Timer timer;
  const int added = separate_and_add_cuts(true, fresh, std::min(kMaxTreeCutsPerNode, room));
  if (added > 0) {
    outcome_.tree_cuts += added;
    outcome_.tree_cuts_separated += added;
    build_relaxation();
  } else if (fresh) {
    fresh_interval_ = std::min(kMaxTreeCutInterval, 2 * fresh_interval_);
    next_fresh_node_ = outcome_.nodes + fresh_interval_;
  }
  tree_cut_seconds_ += timer.seconds();
  return added;
}

void BranchAndBound::root_cut_loop() {
  if (!options_.cuts || integers_.empty()) return;
  const Timer loop_timer;
  if (solve_relaxation(nullptr, -1) != SimplexStatus::kOptimal) return;
  const double root_seconds = loop_timer.seconds();
  double bound = relaxation_objective();
  const double initial_bound = bound;
  outcome_.root_bound = sense_ * bound;
  // Rounds re-solve an ever larger LP; on big models they can cost more than they save.
  double budget = std::max(kCutLoopMinSeconds, kCutLoopRootFactor * root_seconds);
  if (std::isfinite(options_.time_limit)) {
    budget = std::min(budget, kCutLoopTimeFraction * options_.time_limit);
  }
  int stalls = 0;
  for (int round = 0; round < options_.max_cut_rounds && !time_up(); ++round) {
    if (loop_timer.seconds() > budget) break;
    if (std::none_of(integers_.begin(), integers_.end(),
                     [&](Index j) { return is_fractional(x_[j]); })) {
      break;
    }
    std::vector<VarStatus> basis = simplex_->status();
    const Index old_m = m_;
    const int added = separate_and_add_cuts();
    if (added == 0) break;
    basis.insert(basis.end(), static_cast<std::size_t>(added), VarStatus::kBasic);
    build_relaxation();
    if (solve_relaxation(&basis, -1) != SimplexStatus::kOptimal) {
      // Numerical trouble with this round's cuts: drop them.
      std::vector<Index> added_rows;
      for (Index i = old_m; i < m_; ++i) added_rows.push_back(i);
      remove_rows(added_rows);
      build_relaxation();
      basis.resize(basis.size() - static_cast<std::size_t>(added));
      solve_relaxation(&basis, -1);
      break;
    }
    ++outcome_.cut_rounds;
    const double next = relaxation_objective();
    if (next - bound <= kCutStallImprovement * std::max(1.0, std::fabs(bound))) {
      if (++stalls >= kCutStallRounds) {
        bound = next;
        break;
      }
    } else {
      stalls = 0;
    }
    bound = next;
  }

  // Cuts that did not move the bound only slow down every node LP: drop them all.
  if (lp_status_ == SimplexStatus::kOptimal && m_ > original_rows_ &&
      relaxation_objective() - initial_bound <=
          kCutMinTotalGain * std::max(1.0, std::fabs(initial_bound))) {
    std::vector<VarStatus> basis = simplex_->status();
    basis.resize(static_cast<std::size_t>(n_ + original_rows_));
    std::vector<Index> rows;
    for (Index i = original_rows_; i < m_; ++i) rows.push_back(i);
    pool_rows(rows);
    remove_rows(rows);
    build_relaxation();
    solve_relaxation(&basis, -1);
  }
  // Keep only the cuts that are binding at the final root LP.
  if (lp_status_ == SimplexStatus::kOptimal && m_ > original_rows_) {
    std::vector<double> activity(static_cast<std::size_t>(m_), 0.0);
    model_.A.multiply(x_, activity);
    std::vector<Index> slack_rows;
    for (Index i = original_rows_; i < m_; ++i) {
      const double lo = model_.row_lower[i];
      if (activity[i] > lo + 1e-6 * (1.0 + std::fabs(lo))) slack_rows.push_back(i);
    }
    if (!slack_rows.empty()) {
      std::vector<VarStatus> basis = simplex_->status();
      std::vector<VarStatus> kept;
      std::vector<char> drop(static_cast<std::size_t>(m_), 0);
      for (const Index i : slack_rows) drop[i] = 1;
      for (Index v = 0; v < n_ + m_; ++v) {
        if (v >= n_ && drop[v - n_]) continue;
        kept.push_back(basis[v]);
      }
      pool_rows(slack_rows);
      remove_rows(slack_rows);
      build_relaxation();
      solve_relaxation(&kept, -1);
    }
  }
  outcome_.cuts_added = m_ - original_rows_;
  if (lp_status_ == SimplexStatus::kOptimal) {
    root_basis_ = simplex_->status();
    outcome_.root_bound_cuts = sense_ * relaxation_objective();
  }
  log_.log(1, "MIP root bound %.10g -> %.10g after %d cut rounds, %d cuts kept, %.2f s",
           outcome_.root_bound, outcome_.root_bound_cuts, outcome_.cut_rounds,
           outcome_.cuts_added, loop_timer.seconds());
}

// ---------------------------------------------------------------------------------------------
// Branching

Index BranchAndBound::select_branching(const std::vector<Index>& fractional,
                                       const std::vector<double>& x, double objective,
                                       const std::vector<VarStatus>& basis,
                                       bool& node_infeasible, BoundChange& tighten,
                                       bool& has_tighten) {
  struct Candidate {
    Index col;
    double score;
  };
  std::vector<Candidate> candidates;
  candidates.reserve(fractional.size());
  for (const Index j : fractional) {
    const double down = x[j] - std::floor(x[j]);
    const double up = std::ceil(x[j]) - x[j];
    const double score = std::max(pseudocost(j, false) * down, kScoreFloor) *
                         std::max(pseudocost(j, true) * up, kScoreFloor);
    candidates.push_back({j, score});
  }
  std::sort(candidates.begin(), candidates.end(),
            [](const Candidate& a, const Candidate& b) { return a.score > b.score; });

  Index best = -1;
  double best_score = -1.0;
  int strong = 0;
  // A strong-branching child found infeasible or cut off yields a proof just as a node LP does:
  // its bounds are still set, and the proof itself holds everywhere.
  const auto learn_from_child = [&](bool infeasible) {
    if (!options_.conflicts) return;
    const long long before = outcome_.conflicts_found;
    if (infeasible) {
      analyze_infeasible_lp();
    } else {
      analyze_cutoff_lp();
    }
    outcome_.strong_branching_conflicts += outcome_.conflicts_found - before;
  };
  for (const Candidate& c : candidates) {
    const Index j = c.col;
    double score = c.score;
    const bool reliable = std::min(pc_count_[0][j], pc_count_[1][j]) >= options_.reliability;
    if (!reliable && strong < options_.max_strong_branching && strong_budget_left() &&
        !time_up()) {
      ++strong;
      const double lo = lower_[j];
      const double up = upper_[j];
      const double down_value = std::floor(x[j]);
      const double up_value = std::ceil(x[j]);
      double gain[2] = {0.0, 0.0};
      bool cut[2] = {false, false};
      for (int d = 0; d < 2; ++d) {
        if (d == 0) {
          set_bound(j, lo, down_value);
        } else {
          set_bound(j, up_value, up);
        }
        const SimplexStatus st = solve_relaxation(&basis, options_.strong_iterations);
        outcome_.strong_branching_iterations += simplex_->iterations();
        const double distance = d == 0 ? x[j] - down_value : up_value - x[j];
        if (st == SimplexStatus::kInfeasible) {
          cut[d] = true;
          learn_from_child(true);
        } else if (st == SimplexStatus::kOptimal || st == SimplexStatus::kIterationLimit) {
          const double child = relaxation_objective();
          gain[d] = std::max(0.0, child - objective);
          // Record the gain even when the child is cut off: otherwise, once a good incumbent
          // exists, such columns never become reliable and are strong-branched at every node
          // (mas76: 25k strong-branching LPs instead of 1.5k).
          record_pseudocost(j, d == 1, gain[d], distance);
          if (st == SimplexStatus::kOptimal && effective_bound(child) >= cutoff()) {
            cut[d] = true;
            learn_from_child(false);
          }
        } else {
          gain[d] = pseudocost(j, d == 1) * distance;
        }
      }
      set_bound(j, lo, up);
      if (cut[0] && cut[1]) {
        node_infeasible = true;
        return -1;
      }
      if (cut[0] || cut[1]) {
        tighten = cut[0] ? BoundChange{j, up_value, up} : BoundChange{j, lo, down_value};
        has_tighten = true;
        return j;
      }
      score = std::max(gain[0], kScoreFloor) * std::max(gain[1], kScoreFloor);
    }
    if (score > best_score) {
      best_score = score;
      best = j;
    }
  }
  return best;
}

// Reduced-cost fixing: an integer column nonbasic at a bound with reduced cost d_j can move k
// units into its domain only if objective + |d_j| k stays below the cutoff (LP duality), so its
// other bound tightens to the largest such k. The changes hold in the node's subtree and are
// recorded with the node. The current LP point stays feasible and optimal.
void BranchAndBound::reduced_cost_fixing(Node& node, double objective,
                                         const std::vector<double>& reduced,
                                         const std::vector<VarStatus>& basis) {
  const double limit = cutoff();
  if (!(limit < kInf)) return;
  // Slack for the LP's own tolerances, so round-off never fixes a column wrongly.
  const double room = limit - objective + kReducedCostSlack * (1.0 + std::fabs(limit));
  if (room < 0.0) return;
  for (const Index j : integers_) {
    if (lower_[j] == upper_[j]) continue;
    const double d = reduced[j];
    if (basis[j] == VarStatus::kAtLower && d > kReducedCostTol) {
      const double up = lower_[j] + std::floor(room / d);
      if (up < upper_[j]) {
        set_bound(j, lower_[j], up);
        node.local.push_back({j, lower_[j], up});
        ++outcome_.reduced_cost_fixings;
      }
    } else if (basis[j] == VarStatus::kAtUpper && d < -kReducedCostTol) {
      const double lo = upper_[j] - std::floor(room / -d);
      if (lo > lower_[j]) {
        set_bound(j, lo, upper_[j]);
        node.local.push_back({j, lo, upper_[j]});
        ++outcome_.reduced_cost_fixings;
      }
    }
  }
}

// ---------------------------------------------------------------------------------------------
// Nodes

BranchAndBound::NodeResult BranchAndBound::process_node(Node& node, std::vector<Node>& children) {
  apply_node_bounds(node);
  if (node.branch_col >= 0 && !propagate({node.branch_col}, &node.local)) {
    return NodeResult::kPruned;
  }
  if (!conflicts_.empty() && !propagate_conflicts(&node.local)) {
    ++outcome_.conflict_prunes;
    return NodeResult::kPruned;
  }
  std::shared_ptr<const std::vector<VarStatus>> start = node.basis;
  bool cuts_tried = false;
  double before_cuts = kInf;  // The objective before this node's cut round, until re-solved.
  int cut_rows_added = 0;
  std::shared_ptr<const std::vector<VarStatus>> pre_cut_basis;
  for (int round = 0;; ++round) {
    const std::vector<VarStatus>* start_basis =
        start ? start.get() : (root_basis_.empty() ? nullptr : &root_basis_);
    SimplexStatus status = solve_relaxation(start_basis, -1);
    if (status == SimplexStatus::kNumericalError) status = solve_relaxation(nullptr, -1);
    switch (status) {
      case SimplexStatus::kOptimal: break;
      case SimplexStatus::kInfeasible:
        if (options_.conflicts) analyze_infeasible_lp();
        return NodeResult::kPruned;
      case SimplexStatus::kUnbounded:
        return node.depth == 0 ? NodeResult::kUnbounded : NodeResult::kFailed;
      case SimplexStatus::kTimeLimit: return NodeResult::kStopped;
      case SimplexStatus::kIterationLimit:
      case SimplexStatus::kNumericalError: return NodeResult::kFailed;
    }
    const double objective = relaxation_objective();
    if (std::isfinite(before_cuts)) {
      // The node's cut round: keep its rows only if they closed enough of the node's gap (or, with
      // no incumbent, moved the bound as a root round must). They are the LP's last rows and no
      // stored basis covers them yet, so removing them restores the previous LP exactly.
      const double gain = objective - before_cuts;
      const bool has_gap = cutoff() < kInf;
      const double needed = has_gap ? kTreeCutMinGapShare * (cutoff() - before_cuts)
                                    : kCutStallImprovement * std::max(1.0, std::fabs(before_cuts));
      ++tree_cut_rounds_;
      before_cuts = kInf;
      if (gain < needed && effective_bound(objective) < cutoff()) {
        const Timer undo_timer;
        std::vector<Index> rows;
        for (Index i = m_ - cut_rows_added; i < m_; ++i) rows.push_back(i);
        remove_rows(rows);
        build_relaxation();
        outcome_.tree_cuts -= cut_rows_added;
        tree_cut_seconds_ += undo_timer.seconds();
        tree_cut_interval_ = std::min(kMaxTreeCutInterval, 2 * tree_cut_interval_);
        next_tree_cut_node_ = outcome_.nodes + tree_cut_interval_;
        start = std::move(pre_cut_basis);
        continue;
      }
      ++tree_cut_kept_rounds_;
      tree_cut_interval_ = 1;
      fresh_interval_ = 1;
    }
    if (round == 0 && node.branch_col >= 0) {
      record_pseudocost(node.branch_col, node.branch_up, objective - node.bound,
                        node.branch_distance);
    }
    const double bound = effective_bound(objective);
    if (bound >= cutoff()) {
      pruned_bound_ = std::min(pruned_bound_, bound);
      if (options_.conflicts) analyze_cutoff_lp();
      return NodeResult::kPruned;
    }
    node.bound = std::max(node.bound, objective);
    const std::vector<double> x = x_;
    // The node LP's reduced costs (original units), kept for reduced-cost fixing after the
    // heuristics, which re-solve other LPs on the same simplex.
    std::vector<double> reduced(static_cast<std::size_t>(n_));
    for (Index j = 0; j < n_; ++j) reduced[j] = simplex_->reduced_costs()[j] / scaling_.col[j];
    std::vector<Index> fractional;
    for (const Index j : integers_) {
      if (is_fractional(x[j])) fractional.push_back(j);
    }
    auto basis = std::make_shared<const std::vector<VarStatus>>(current_basis());
    if (fractional.empty()) {
      if (try_solution(x)) return NodeResult::kPruned;
      // No better than the incumbent: the node's LP bound rules out anything better.
      if (!improves(objective)) return NodeResult::kPruned;
      // Integral within the tolerance, yet rounding the integers broke a row (a large
      // coefficient times a tiny rounding): re-solve the continuous columns with the integers
      // fixed. If that fails, branch on the columns that are not exactly integral (x <= k or
      // x >= k + 1 is a valid split at any value); only an exactly integral point that still
      // fails leaves the node unresolved.
      if (round_and_solve(x, *basis, search::kRoundAndSolveIterations) || bound >= cutoff()) {
        return NodeResult::kPruned;
      }
      for (const Index j : integers_) {
        // Both children must be strictly smaller: a value just outside a bound (within the LP
        // tolerance) gives no split.
        if (x[j] != std::round(x[j]) && std::floor(x[j]) >= lower_[j] &&
            std::ceil(x[j]) <= upper_[j]) {
          fractional.push_back(j);
        }
      }
      if (fractional.empty()) {
        // The LP point is its integer rounding up to the LP tolerances. If that rounding does
        // not beat the incumbent, the difference is round-off and nothing better is left here.
        double rounded = offset_;
        for (Index k = 0; k < n_; ++k) {
          const bool integer = model_.col_type[k] == VarType::kInteger;
          rounded += cost_[k] * (integer ? std::round(x[k]) : x[k]);
        }
        return improves(rounded) ? NodeResult::kFailed : NodeResult::kPruned;
      }
    }
    if (node.depth == 0 && round == 0) {
      root_basis_ = *basis;
      log_.log(1, "MIP root relaxation %.12g, %zu fractional of %zu integers, %lld iterations",
               sense_ * objective, fractional.size(), integers_.size(), outcome_.lp_iterations);
    }

    simple_rounding(x);
    if (node.depth == 0 || outcome_.nodes % kRoundAndSolveFrequency == 0) {
      round_and_solve(x, *basis, search::kRoundAndSolveIterations);
    }
    if (options_.heuristics && round == 0) run_heuristics(node, x, *basis);
    if (bound >= cutoff()) {
      pruned_bound_ = std::min(pruned_bound_, bound);
      return NodeResult::kPruned;
    }
    reduced_cost_fixing(node, objective, reduced, *basis);
    if (node.depth == 0 && round == 0) {
      std::size_t fixed = 0;
      for (const Index j : integers_) fixed += lower_[j] == upper_[j];
      log_.log(1, "MIP root: %zu of %zu integers fixed, %lld reduced-cost fixings", fixed,
               integers_.size(), outcome_.reduced_cost_fixings);
    }

    if (!cuts_tried) {
      cuts_tried = true;
      const int added = tree_cut_round(node.depth);
      if (added > 0) {
        before_cuts = objective;
        cut_rows_added = added;
        pre_cut_basis = basis;
        auto padded = std::make_shared<std::vector<VarStatus>>(*basis);
        padded->resize(static_cast<std::size_t>(n_ + m_), VarStatus::kBasic);
        start = std::move(padded);
        continue;
      }
    }

    bool infeasible = false;
    bool has_tighten = false;
    BoundChange tighten{};
    Index j = select_branching(fractional, x, objective, *basis, infeasible, tighten,
                               has_tighten);
    if (infeasible) return NodeResult::kPruned;
    if (has_tighten && round < kMaxTightenRounds) {
      set_bound(tighten.col, tighten.lower, tighten.upper);
      node.local.push_back(tighten);
      if (!propagate({tighten.col}, &node.local)) return NodeResult::kPruned;
      start = basis;
      continue;
    }
    if (j < 0) j = fractional.front();
    if (has_tighten) {
      // Out of tightening rounds: branch on the tightened column itself.
      j = tighten.col;
    }

    const double open_bytes = static_cast<double>(open_.size() + 2) *
                              static_cast<double>(n_ + m_) * sizeof(VarStatus);
    std::shared_ptr<const std::vector<VarStatus>> child_basis =
        open_bytes <= kBasisMemoryBudget ? basis : nullptr;
    const double v = x[j];
    Node down;
    down.bound = node.bound;
    down.depth = node.depth + 1;
    // The children share one segment holding this node's changes.
    auto segment = std::make_shared<PathSegment>();
    segment->parent = node.path;
    segment->changes = std::move(node.local);
    down.path = std::move(segment);
    down.local = {{j, lower_[j], std::floor(v)}};
    down.basis = child_basis;
    down.branch_col = j;
    down.branch_up = false;
    down.branch_distance = v - std::floor(v);
    Node up = down;
    up.local = {{j, std::ceil(v), upper_[j]}};
    up.branch_up = true;
    up.branch_distance = std::ceil(v) - v;
    children.push_back(std::move(down));
    children.push_back(std::move(up));
    return NodeResult::kBranched;
  }
}

// Solves the model under the root node's bounds (valid everywhere: they come from the root) as a
// new, presolved search with the incumbent as its cutoff, and maps the result back.
MipOutcome BranchAndBound::restart() {
  const Timer timer;
  Model sub = model_;
  sub.col_lower = lower_;
  sub.col_upper = upper_;
  PresolveOptions presolve_options;
  presolve_options.mip = true;
  presolve_options.integrality_tol = options_.integrality_tol;
  Presolve presolve(sub, presolve_options);
  const bool infeasible = presolve.run() == PresolveStatus::kInfeasible;
  MipOutcome out;
  if (infeasible) {
    out.status = Status::kInfeasible;
  } else if (presolve.reduced().num_cols() == 0) {
    out.status = Status::kOptimal;
    out.x.clear();
    std::vector<double> x;
    std::vector<double> unused_y;
    presolve.postsolve({}, {}, x, unused_y);
    if (try_solution(std::move(x))) ++outcome_.heuristic_solutions;
  } else {
    MipOptions options = options_;
    options.restart = false;
    options.time_limit = std::max(0.0, remaining_time());
    if (options.node_limit >= 0) {
      options.node_limit = std::max(0LL, options.node_limit - outcome_.nodes);
    }
    options.debug_solution.clear();
    options.start.clear();  // The incumbent below already carries whatever the start gave.
    options.objective_cutoff.reset();
    if (incumbent_value_ < kInf) options.objective_cutoff = sense_ * incumbent_value_;
    BranchAndBound search(presolve.reduced(), options, log_);
    out = search.solve();
    if (!out.x.empty()) {
      std::vector<double> x;
      std::vector<double> unused_y;
      presolve.postsolve(out.x, {}, x, unused_y);
      try_solution(std::move(x));
    }
    outcome_.nodes += out.nodes;
    outcome_.lp_iterations += out.lp_iterations;
    outcome_.strong_branching_iterations += out.strong_branching_iterations;
    outcome_.heuristic_solutions += out.heuristic_solutions;
    outcome_.heuristic_lp_iterations += out.heuristic_lp_iterations;
    outcome_.reduced_cost_fixings += out.reduced_cost_fixings;
    outcome_.conflicts_found += out.conflicts_found;
    outcome_.conflict_prunes += out.conflict_prunes;
    outcome_.conflict_infeasible_lps += out.conflict_infeasible_lps;
    outcome_.farkas_proofs += out.farkas_proofs;
    outcome_.strong_branching_conflicts += out.strong_branching_conflicts;
    outcome_.threads_used = out.threads_used;
  }
  outcome_.restarted = true;
  // The restarted search either completed (optimal / infeasible: nothing better than the
  // incumbent exists) or stopped with a bound valid for the whole model.
  double bound = kInf;
  switch (out.status) {
    case Status::kOptimal:
    case Status::kInfeasible:
      outcome_.status = incumbent_.empty() ? Status::kInfeasible : Status::kOptimal;
      bound = incumbent_value_;
      break;
    case Status::kTimeLimit:
    case Status::kNodeLimit:
      outcome_.status = out.status;
      bound = std::min(sense_ * out.bound, incumbent_value_);
      break;
    default:
      outcome_.status = out.status;
      bound = -kInf;
      break;
  }
  outcome_.bound = sense_ * bound;
  if (!incumbent_.empty()) {
    outcome_.x = incumbent_;
    outcome_.objective = sense_ * incumbent_value_;
  } else {
    outcome_.objective = sense_ * kInf;
  }
  log_.log(1, "MIP: restarted search %s after %lld nodes, %.2f s", to_string(out.status),
           out.nodes, timer.seconds());
  return outcome_;
}

MipOutcome BranchAndBound::solve() {
  timer_ = Timer();
  outcome_ = MipOutcome{};
  if (options_.objective_cutoff) incumbent_value_ = sense_ * *options_.objective_cutoff;

  // Root propagation becomes part of the root bounds.
  std::vector<Index> all_columns(static_cast<std::size_t>(n_));
  std::iota(all_columns.begin(), all_columns.end(), Index{0});
  if (!propagate(std::move(all_columns), nullptr)) {
    outcome_.status = Status::kInfeasible;
    return outcome_;
  }
  if (options_.probing && !integers_.empty()) {
    outcome_.coefficients_tightened = tighten_coefficients();
    if (outcome_.coefficients_tightened > 0) {
      log_.log(1, "MIP: %d coefficients tightened", outcome_.coefficients_tightened);
    }
    if (!probe()) {
      outcome_.status = Status::kInfeasible;
      return outcome_;
    }
  }
  root_lower_ = lower_;
  root_upper_ = upper_;
  if (!options_.debug_solution.empty()) {
    const std::vector<double>& d = options_.debug_solution;
    const double tol = options_.feasibility_tol;
    for (Index j = 0; j < n_; ++j) {
      if (d[j] < root_lower_[j] - tol || d[j] > root_upper_[j] + tol) {
        ++outcome_.debug_reduction_violations;
      }
    }
    std::vector<double> activity(static_cast<std::size_t>(m_), 0.0);
    model_.A.multiply(d, activity);
    for (Index i = 0; i < m_; ++i) {
      if (activity[i] < model_.row_lower[i] - tol * (1.0 + std::fabs(model_.row_lower[i])) ||
          activity[i] > model_.row_upper[i] + tol * (1.0 + std::fabs(model_.row_upper[i]))) {
        ++outcome_.debug_reduction_violations;
      }
    }
  }
  for (const Index j : touched_) is_touched_[j] = 0;
  touched_.clear();
  const bool started = !options_.start.empty() && use_start();
  if (options_.heuristics && !integers_.empty() && !started) run_feasibility_jump();
  root_cut_loop();

  std::optional<Node> current = Node{};
  bool unbounded = false;
  bool stopped = false;
  bool gap_closed = false;
  double last_log = 0.0;
  for (;;) {
    if (!current && !dive_stack_.empty()) {
      current = std::move(dive_stack_.back());
      dive_stack_.pop_back();
    }
    if (!current) {
      if (open_.empty()) break;
      auto it = open_.begin();
      current = std::move(it->second);
      open_.erase(it);
    }
    const double node_bound = effective_bound(current->bound);
    if (node_bound >= cutoff()) {
      pruned_bound_ = std::min(pruned_bound_, node_bound);
      current.reset();
      continue;
    }
    if (incumbent_value_ < kInf &&
        incumbent_value_ - best_bound(current->bound) <=
            std::max(options_.abs_gap, options_.rel_gap * std::fabs(incumbent_value_))) {
      open_.emplace(current->bound, std::move(*current));
      current.reset();
      gap_closed = true;
      break;
    }
    if (open_.size() >= options_.max_open_nodes && !open_limit_logged_) {
      open_limit_logged_ = true;
      log_.log(1, "MIP: %zu open nodes, stopping to bound memory use", open_.size());
    }
    if (time_up() || open_.size() >= options_.max_open_nodes ||
        (options_.node_limit >= 0 && outcome_.nodes >= options_.node_limit)) {
      open_.emplace(current->bound, std::move(*current));
      current.reset();
      stopped = true;
      break;
    }

    if (options_.threads > 1 && outcome_.nodes >= options_.parallel_start_nodes &&
        !open_.empty()) {
      open_.emplace(current->bound, std::move(*current));
      current.reset();
      run_parallel(unbounded, stopped, gap_closed);
      break;
    }

    ++outcome_.nodes;
    std::vector<Node> children;
    const NodeResult result = process_node(*current, children);
    if (result == NodeResult::kBranched && current->depth == 0 && options_.restart &&
        shared_ == nullptr) {
      std::size_t fixed = 0;
      for (const Index j : integers_) {
        fixed += lower_[j] == upper_[j] && root_lower_[j] < root_upper_[j];
      }
      if (static_cast<double>(fixed) >=
          kRestartFraction * static_cast<double>(integers_.size())) {
        log_.log(1, "MIP: restart after the root, %zu of %zu integers fixed", fixed,
                 integers_.size());
        return restart();
      }
    }
    if (log_.level() >= 1 && timer_.seconds() - last_log >= 5.0) {
      last_log = timer_.seconds();
      const double bound = best_bound(current->bound);
      log_.log(1, "MIP nodes %lld, open %zu, bound %.10g, incumbent %.10g, %.1f s",
               outcome_.nodes, open_.size(), sense_ * bound, sense_ * incumbent_value_,
               timer_.seconds());
    }
    if (result == NodeResult::kUnbounded) {
      unbounded = true;
      break;
    }
    if (result == NodeResult::kStopped) {
      open_.emplace(current->bound, std::move(*current));
      current.reset();
      stopped = true;
      break;
    }
    if (result == NodeResult::kFailed) incomplete_ = true;
    if (result != NodeResult::kBranched) {
      current.reset();
      continue;
    }
    // Plunge into the child on the side the value is closer to; queue the other.
    const std::size_t dive = plunge_child(children);
    Node& other = children[1 - dive];
    Node& next = children[dive];
    if (!dive_stack_.empty() || open_.size() >= options_.max_open_nodes_soft) {
      // Memory mode: depth first with a LIFO stack, so the stored nodes are bounded by the
      // depth instead of growing with the node count. Best-bound resumes when it empties.
      dive_stack_.push_back(std::move(other));
      current = std::move(next);
      continue;
    }
    open_.emplace(other.bound, std::move(other));
    bool keep = false;
    if (incumbent_value_ == kInf) {
      keep = next.depth < search::kMaxPlungeDepth;
    } else {
      const double best = open_.empty() ? next.bound : open_.begin()->first;
      keep = effective_bound(next.bound) <=
             best + search::kPlungeGapFraction * (incumbent_value_ - best);
    }
    if (keep) {
      current = std::move(next);
    } else {
      open_.emplace(next.bound, std::move(next));
      current.reset();
    }
  }

  if (outcome_.conflicts_found > 0) {
    log_.log(1,
             "MIP conflicts: %lld proofs stored (%lld from %lld infeasible LPs, %lld from "
             "strong branching), %lld nodes pruned by them",
             outcome_.conflicts_found, outcome_.farkas_proofs, outcome_.conflict_infeasible_lps,
             outcome_.strong_branching_conflicts, outcome_.conflict_prunes);
  }
  if (tree_cut_rounds_ > 0) {
    log_.log(1, "MIP tree cuts: %d kept from %lld of %lld rounds, %zu left in the pool, %.2f s",
             outcome_.tree_cuts, tree_cut_kept_rounds_, tree_cut_rounds_, cut_pool_.size(),
             tree_cut_seconds_);
  }
  // Final status and bound (minimization internally).
  double bound = kInf;
  if (unbounded) {
    outcome_.status = incumbent_.empty() ? Status::kInfeasibleOrUnbounded : Status::kUnbounded;
    bound = -kInf;
  } else if (stopped) {
    outcome_.status = time_up() ? Status::kTimeLimit : Status::kNodeLimit;
    bound = std::min(best_bound(kInf), pruned_bound_);
    if (outcome_.nodes == 0 || !std::isfinite(bound)) bound = -kInf;
  } else if (gap_closed) {
    outcome_.status = Status::kOptimal;
    bound = std::min({best_bound(kInf), pruned_bound_, incumbent_value_});
  } else if (incomplete_) {
    outcome_.status = Status::kNumericalError;
    bound = -kInf;
  } else if (incumbent_.empty()) {
    outcome_.status = Status::kInfeasible;
  } else {
    outcome_.status = Status::kOptimal;
    bound = std::min(incumbent_value_, pruned_bound_);
  }
  outcome_.bound = sense_ * bound;
  if (!incumbent_.empty()) {
    outcome_.x = incumbent_;
    outcome_.objective = sense_ * incumbent_value_;
  } else {
    outcome_.objective = sense_ * kInf;
  }
  return outcome_;
}

}  // namespace samaya
