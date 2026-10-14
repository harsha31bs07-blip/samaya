#include "lp/simplex.hpp"

#include <algorithm>
#include <cmath>

namespace samaya {

const char* to_string(SimplexStatus status) {
  switch (status) {
    case SimplexStatus::kOptimal: return "optimal";
    case SimplexStatus::kInfeasible: return "infeasible";
    case SimplexStatus::kUnbounded: return "unbounded";
    case SimplexStatus::kIterationLimit: return "iteration_limit";
    case SimplexStatus::kTimeLimit: return "time_limit";
    case SimplexStatus::kNumericalError: return "numerical_error";
  }
  return "unknown";
}

namespace {

using Component = SimplexStats::Component;
using RebuildReason = BasisFactor::RebuildReason;

// One clock read per boundary; unclassified work is obtained from the solve's wall time.
// Segments marked rest do not accumulate, so nested rebuild/flip timers cannot double count.
class ProfileTimer {
 public:
  ProfileTimer(bool enabled, SimplexStats& stats, Component component,
               double* maximum = nullptr, bool iteration = false)
      : enabled_(enabled), stats_(stats), component_(component), maximum_(maximum),
        iteration_(iteration) {
    if (enabled_) {
      if (iteration_) initial_measured_ = measured_seconds();
      start_ = first_start_ = std::chrono::steady_clock::now();
    }
  }
  ~ProfileTimer() {
    switch_to(Component::kRest);
    if (enabled_ && iteration_) {
      // The enclosing iteration owns unclassified time. Subtract its measured segments and
      // nested flip/update/rebuild timers, so refreshes cannot be counted twice.
      const double wall = std::chrono::duration<double>(start_ - first_start_).count();
      stats_.seconds[static_cast<std::size_t>(Component::kRest)] +=
          std::max(0.0, wall - (measured_seconds() - initial_measured_));
    }
    if (enabled_ && maximum_ != nullptr) *maximum_ = std::max(*maximum_, elapsed_);
  }
  void switch_to(Component component) {
    if (!enabled_) return;
    const auto now = std::chrono::steady_clock::now();
    if (component_ != Component::kRest) {
      const double seconds = std::chrono::duration<double>(now - start_).count();
      stats_.seconds[static_cast<std::size_t>(component_)] += seconds;
      elapsed_ += seconds;
    }
    start_ = now;
    component_ = component;
  }

 private:
  double measured_seconds() const {
    double sum = 0.0;
    for (std::size_t k = 0; k < static_cast<std::size_t>(Component::kRest); ++k) {
      sum += stats_.seconds[k];
    }
    return sum;
  }
  bool enabled_;
  SimplexStats& stats_;
  Component component_;
  double* maximum_;
  bool iteration_;
  double elapsed_ = 0.0;
  double initial_measured_ = 0.0;
  std::chrono::steady_clock::time_point start_{};
  std::chrono::steady_clock::time_point first_start_{};
};

constexpr double kMinDseWeight = 1e-4;
constexpr double kPhase1FreeBound = 1000.0;
constexpr double kPerturbationBase = 5e-7;
constexpr double kFeasibilityCostBase = 1e-3;
// Relative disagreement between the pivot element computed from the column (ftran) and the row
// (btran) that triggers a refactorization.
constexpr double kPivotConsistencyTol = 1e-7;

BasisFactor::Options factor_options(const SimplexOptions& options, const Logger& log) {
  BasisFactor::Options f;
  f.max_updates = options.refactor_interval;
  f.collect_stats = options.collect_stats || log.enabled(2);
  return f;
}

}  // namespace

void Simplex::set_iteration_limit(long long limit) {
  options_.max_iterations = limit;
  max_iterations_ = limit >= 0 ? limit : std::max<long long>(100000, 100LL * static_cast<long long>(nt_));
}

Simplex::Simplex(const LpProblem& lp, const SimplexOptions& options, const Logger& log)
    : lp_(lp),
      options_(options),
      log_(log),
      m_(lp.m),
      n_(lp.n),
      nt_(lp.n + lp.m),
      factor_(lp.A, factor_options(options, log)),
      rng_(options.seed) {
  set_iteration_limit(options.max_iterations);
  const auto m = static_cast<std::size_t>(m_);
  const auto nt = static_cast<std::size_t>(nt_);
  basic_.resize(m);
  position_.assign(nt, -1);
  status_.assign(nt, VarStatus::kAtLower);
  x_.assign(nt, 0.0);
  y_.assign(m, 0.0);
  d_.assign(nt, 0.0);
  dse_weight_.assign(m, 1.0);
  devex_weight_.assign(nt, 1.0);
  alpha_row_.assign(nt, 0.0);
  in_row_nz_.assign(nt, 0);
  alpha_col_.assign(m, 0.0);
  in_infeasible_.assign(m, 0);
}

template <typename F>
void Simplex::for_column(Index j, F&& f) const {
  if (j < n_) {
    const auto start = lp_.A.col_start();
    const auto index = lp_.A.row_index();
    const auto value = lp_.A.values();
    for (NnzIndex p = start[j]; p < start[j + 1]; ++p) f(index[p], value[p]);
  } else {
    f(j - n_, -1.0);
  }
}

void Simplex::load_column(Index j, std::vector<double>& out) const {
  out.assign(static_cast<std::size_t>(m_), 0.0);
  for_column(j, [&](Index i, double a) { out[i] = a; });
}

void Simplex::set_value_from_status(Index j) {
  switch (status_[j]) {
    case VarStatus::kAtLower: x_[j] = lower_[j]; break;
    case VarStatus::kAtUpper: x_[j] = upper_[j]; break;
    case VarStatus::kAtZero: x_[j] = 0.0; break;
    case VarStatus::kBasic: break;
  }
}

bool Simplex::rebuild(RebuildReason reason) {
  if (profile_) ++stats_.rebuild_calls;
  const ProfileTimer measurement(profile_, stats_, Component::kRebuild,
                                 &stats_.max_rebuild_seconds);
  for (int attempt = 0; attempt < 4; ++attempt) {
    ++stats_.refactorizations;
    if (factor_.factorize(basic_, attempt == 0 ? reason : RebuildReason::kRankRepair) == 0) {
      compute_primal();
      compute_dual();
      factor_valid_ = true;
      return true;
    }
    // Replace dependent basis columns with the logicals of the rows left unpivoted.
    const auto& positions = factor_.singular_positions();
    const auto& rows = factor_.unpivoted_rows();
    log_.log(3, "simplex: basis rank deficiency %zu, repairing", positions.size());
    for (std::size_t k = 0; k < positions.size(); ++k) {
      const Index pos = positions[k];
      const Index entering = n_ + rows[k];
      const Index leaving = basic_[pos];
      if (position_[entering] != -1) return false;
      if (profile_) {
        ++stats_.rank_repairs;
        if (warm_start_setup_) ++stats_.warm_start_rank_repairs;
      }
      const bool has_lower = lower_[leaving] > -kInf;
      const bool has_upper = upper_[leaving] < kInf;
      if (has_lower && (!has_upper || std::fabs(x_[leaving] - lower_[leaving]) <=
                                          std::fabs(x_[leaving] - upper_[leaving]))) {
        status_[leaving] = VarStatus::kAtLower;
      } else if (has_upper) {
        status_[leaving] = VarStatus::kAtUpper;
      } else {
        status_[leaving] = VarStatus::kAtZero;
      }
      set_value_from_status(leaving);
      position_[leaving] = -1;
      basic_[pos] = entering;
      position_[entering] = pos;
      status_[entering] = VarStatus::kBasic;
      dse_weight_[pos] = 1.0;
    }
  }
  factor_valid_ = false;
  return false;
}

void Simplex::compute_primal() {
  std::vector<double> rhs(static_cast<std::size_t>(m_), 0.0);
  for (Index j = 0; j < nt_; ++j) {
    if (status_[j] == VarStatus::kBasic) continue;
    const double v = x_[j];
    if (v == 0.0) continue;
    for_column(j, [&](Index i, double a) { rhs[i] -= a * v; });
  }
  factor_.ftran(rhs);
  for (Index k = 0; k < m_; ++k) x_[basic_[k]] = rhs[k];
  rebuild_infeasible_list();
}

bool Simplex::primal_infeasible(Index j) const {
  const double tol = options_.primal_tol;
  return x_[j] < lower_[j] - tol || x_[j] > upper_[j] + tol;
}

void Simplex::rebuild_infeasible_list() {
  for (const Index k : infeasible_) in_infeasible_[k] = 0;
  infeasible_.clear();
  for (Index k = 0; k < m_; ++k) {
    if (primal_infeasible(basic_[k])) {
      in_infeasible_[k] = 1;
      infeasible_.push_back(k);
    }
  }
}

void Simplex::note_basic_change(Index position) {
  if (!in_infeasible_[position] && primal_infeasible(basic_[position])) {
    in_infeasible_[position] = 1;
    infeasible_.push_back(position);
  }
}

void Simplex::compute_dual() {
  std::vector<double> cb(static_cast<std::size_t>(m_));
  for (Index k = 0; k < m_; ++k) cb[k] = cost_[basic_[k]];
  factor_.btran(cb);
  y_ = std::move(cb);
  for (Index j = 0; j < nt_; ++j) {
    if (status_[j] == VarStatus::kBasic) {
      d_[j] = 0.0;
      continue;
    }
    double s = cost_[j];
    for_column(j, [&](Index i, double a) { s -= a * y_[i]; });
    d_[j] = s;
  }
}

void Simplex::compute_pivot_row() {
  // alpha_row = rho' [A -I], accumulated row-wise so only rows with rho_i != 0 are touched.
  for (const Index j : row_nz_) {
    alpha_row_[j] = 0.0;
    in_row_nz_[j] = 0;
  }
  row_nz_.clear();
  const auto start = lp_.At.col_start();
  const auto index = lp_.At.row_index();
  const auto value = lp_.At.values();
  long long rho_nnz = 0;
  for (Index i = 0; i < m_; ++i) {
    const double r = rho_[i];
    if (r == 0.0) continue;
    ++rho_nnz;
    for (NnzIndex p = start[i]; p < start[i + 1]; ++p) {
      const Index j = index[p];
      if (!in_row_nz_[j]) {
        in_row_nz_[j] = 1;
        row_nz_.push_back(j);
      }
      alpha_row_[j] += r * value[p];
    }
    in_row_nz_[n_ + i] = 1;
    row_nz_.push_back(n_ + i);
    alpha_row_[n_ + i] = -r;
  }
  if (profile_) stats_.rho_nnz.add(rho_nnz);
}

void Simplex::place_nonbasic_dual_feasible() {
  for (Index j = 0; j < nt_; ++j) {
    if (status_[j] == VarStatus::kBasic) continue;
    const bool has_lower = lower_[j] > -kInf;
    const bool has_upper = upper_[j] < kInf;
    if (has_lower && has_upper) {
      status_[j] = (lower_[j] == upper_[j] || d_[j] >= 0.0) ? VarStatus::kAtLower
                                                             : VarStatus::kAtUpper;
    } else if (has_lower) {
      status_[j] = VarStatus::kAtLower;
    } else if (has_upper) {
      status_[j] = VarStatus::kAtUpper;
    } else {
      status_[j] = VarStatus::kAtZero;
    }
    set_value_from_status(j);
  }
}

bool Simplex::correct_dual_infeasibilities() {
  bool flipped = false;
  bool changed = false;
  const double tol = options_.dual_tol;
  for (Index j = 0; j < nt_; ++j) {
    const VarStatus st = status_[j];
    if (st == VarStatus::kBasic || lower_[j] == upper_[j]) continue;
    const double dj = d_[j];
    const bool infeasible = (st == VarStatus::kAtLower && dj < -tol) ||
                            (st == VarStatus::kAtUpper && dj > tol) ||
                            (st == VarStatus::kAtZero && std::fabs(dj) > tol);
    if (!infeasible) continue;
    changed = true;
    if (lower_[j] > -kInf && upper_[j] < kInf) {
      status_[j] = dj >= 0.0 ? VarStatus::kAtLower : VarStatus::kAtUpper;
      set_value_from_status(j);
      flipped = true;
      if (profile_ && warm_start_setup_) ++stats_.warm_start_bound_flips;
    } else {
      // Shift the cost so that d_j = 0; shifts are removed before optimality is declared.
      ++stats_.cost_shifts;
      if (profile_ && warm_start_setup_) ++stats_.warm_start_cost_shifts;
      cost_[j] -= dj;
      d_[j] = 0.0;
    }
  }
  if (flipped) compute_primal();
  return changed;
}

void Simplex::perturb_costs(double base) {
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  for (Index j = 0; j < nt_; ++j) {
    const bool has_lower = lower_[j] > -kInf;
    const bool has_upper = upper_[j] < kInf;
    if (lower_[j] == upper_[j] || (!has_lower && !has_upper)) continue;
    const double xi = base * (1.0 + std::fabs(cost_[j])) * (1.0 + unit(rng_));
    switch (status_[j]) {
      case VarStatus::kAtLower: cost_[j] += xi; break;
      case VarStatus::kAtUpper: cost_[j] -= xi; break;
      case VarStatus::kAtZero: break;
      case VarStatus::kBasic:
        if (has_lower && !has_upper) {
          cost_[j] += xi;
        } else if (has_upper && !has_lower) {
          cost_[j] -= xi;
        } else {
          cost_[j] += cost_[j] >= 0.0 ? xi : -xi;
        }
        break;
    }
  }
}

Index Simplex::count_dual_infeasibilities_unboxed() const {
  Index count = 0;
  const double tol = options_.dual_tol;
  for (Index j = 0; j < nt_; ++j) {
    if (status_[j] == VarStatus::kBasic) continue;
    const bool has_lower = lower_[j] > -kInf;
    const bool has_upper = upper_[j] < kInf;
    if (has_lower && has_upper) continue;
    const double dj = d_[j];
    if ((has_lower && dj < -tol) || (has_upper && dj > tol) ||
        (!has_lower && !has_upper && std::fabs(dj) > tol)) {
      ++count;
    }
  }
  return count;
}

double Simplex::max_primal_infeasibility() const {
  double worst = 0.0;
  for (Index k = 0; k < m_; ++k) {
    const Index j = basic_[k];
    worst = std::max({worst, lower_[j] - x_[j], x_[j] - upper_[j]});
  }
  return worst;
}

double Simplex::max_dual_infeasibility() const {
  double worst = 0.0;
  for (Index j = 0; j < nt_; ++j) {
    if (lower_[j] == upper_[j]) continue;
    switch (status_[j]) {
      case VarStatus::kAtLower: worst = std::max(worst, -d_[j]); break;
      case VarStatus::kAtUpper: worst = std::max(worst, d_[j]); break;
      case VarStatus::kAtZero: worst = std::max(worst, std::fabs(d_[j])); break;
      case VarStatus::kBasic: break;
    }
  }
  return worst;
}

bool Simplex::limit_reached(SimplexStatus& status) const {
  if (iterations_ >= max_iterations_) {
    status = SimplexStatus::kIterationLimit;
    return true;
  }
  if ((iterations_ & 31) == 0 && timer_.seconds() > options_.time_limit) {
    status = SimplexStatus::kTimeLimit;
    return true;
  }
  return false;
}

Index Simplex::choose_leaving_row() {
  Index best = -1;
  double best_score = 0.0;
  const double tol = options_.primal_tol;
  // Dual steepest-edge score of position k, or -1 when it is feasible.
  const auto score_of = [&](Index k) {
    const Index j = basic_[k];
    double infeasibility = 0.0;
    if (x_[j] < lower_[j] - tol) {
      infeasibility = lower_[j] - x_[j];
    } else if (x_[j] > upper_[j] + tol) {
      infeasibility = x_[j] - upper_[j];
    } else {
      return -1.0;
    }
    // A non-finite weight would hide the row from pricing forever.
    const double w = std::isfinite(dse_weight_[k]) ? dse_weight_[k] : 1.0;
    return infeasibility * infeasibility / w;
  };
  if (options_.chuzr_full_scan) {
    for (Index k = 0; k < m_; ++k) {
      const double score = score_of(k);
      if (score > best_score) {
        best_score = score;
        best = k;
      }
    }
    return best;
  }
  if (profile_) stats_.chuzr_list.add(static_cast<long long>(infeasible_.size()));
  // The list is unordered: an equal score goes to the lower position, as in the scan, so both
  // choose the same row.
  for (std::size_t t = 0; t < infeasible_.size();) {
    const Index k = infeasible_[t];
    const double score = score_of(k);
    if (score < 0.0) {
      in_infeasible_[k] = 0;
      infeasible_[t] = infeasible_.back();
      infeasible_.pop_back();
      continue;
    }
    if (score > best_score || (score == best_score && best >= 0 && k < best)) {
      best_score = score;
      best = k;
    }
    ++t;
  }
  return best;
}

Index Simplex::dual_ratio_test(double direction, double slope) {
  // Along the dual step t >= 0, d_j(t) = d_j - t * a_j with a_j = direction * alpha_rj.
  //
  // Bound-flipping ("long step") ratio test: passing the breakpoint of a boxed variable only
  // requires flipping it to its other bound, which lowers the slope of the dual objective (the
  // primal infeasibility of the leaving row) by |a_j| * (u_j - l_j). Breakpoints are passed in
  // Harris groups while the slope stays positive; the entering variable is the largest pivot of
  // the group where it would turn non-positive or a variable cannot be flipped.
  const double tol = options_.dual_tol;
  const double pivot_tol = options_.pivot_tol;
  flips_.clear();
  candidates_.clear();
  for (const Index j : row_nz_) {
    const VarStatus st = status_[j];
    const double a = direction * alpha_row_[j];
    if (st == VarStatus::kBasic || lower_[j] == upper_[j] || std::fabs(a) < pivot_tol) continue;
    const bool candidate = a > 0.0 ? (st == VarStatus::kAtLower || st == VarStatus::kAtZero)
                                   : (st == VarStatus::kAtUpper || st == VarStatus::kAtZero);
    if (candidate) candidates_.push_back(j);
  }

  if (profile_) stats_.ratio_candidates.add(static_cast<long long>(candidates_.size()));
  long long groups = 0;
  const auto finish = [&](Index entering) {
    if (profile_) stats_.harris_groups.add(groups);
    return entering;
  };
  std::size_t remaining = candidates_.size();
  while (remaining > 0) {
    ++groups;
    // Harris pass 1: largest step keeping every remaining d_j feasible within the tolerance.
    double max_step = kInf;
    for (std::size_t k = 0; k < remaining; ++k) {
      const Index j = candidates_[k];
      const double a = direction * alpha_row_[j];
      max_step = std::min(max_step, a > 0.0 ? (d_[j] + tol) / a : (d_[j] - tol) / a);
    }
    // The group of breakpoints within that step; pass 2 picks its largest pivot.
    Index best = -1;
    double best_abs = 0.0;
    double reduction = 0.0;
    for (std::size_t k = 0; k < remaining; ++k) {
      const Index j = candidates_[k];
      const double a = direction * alpha_row_[j];
      if (d_[j] / a > max_step) continue;
      if (std::fabs(a) > best_abs) {
        best_abs = std::fabs(a);
        best = j;
      }
      reduction += std::fabs(a) * (upper_[j] - lower_[j]);  // Infinite if not boxed.
    }
    // Stop at this group once the remaining infeasibility would be within the tolerance: the row
    // becomes feasible at (or numerically at) this breakpoint.
    if (!(slope - reduction > options_.primal_tol)) return finish(best);

    // Pass the whole group: its variables flip to their other bounds.
    for (std::size_t k = 0; k < remaining;) {
      const Index j = candidates_[k];
      if (d_[j] / (direction * alpha_row_[j]) <= max_step) {
        flips_.push_back(j);
        candidates_[k] = candidates_[--remaining];
      } else {
        ++k;
      }
    }
    slope -= reduction;
  }
  // Even with every candidate flipped the leaving row stays infeasible: no entering variable.
  flips_.clear();
  return finish(-1);
}

void Simplex::apply_flips() {
  // Flipping boxed nonbasic variables changes x_B by -B^-1 (sum_j a_j * change_j).
  std::vector<double> column(static_cast<std::size_t>(m_), 0.0);
  long long rhs_nnz = 0;
  for (const Index j : flips_) {
    const double old_value = x_[j];
    status_[j] = status_[j] == VarStatus::kAtUpper ? VarStatus::kAtLower : VarStatus::kAtUpper;
    set_value_from_status(j);
    const double change = x_[j] - old_value;
    for_column(j, [&](Index i, double a) {
      const double old = column[i];
      column[i] += a * change;
      if (profile_) rhs_nnz += (column[i] != 0.0) - (old != 0.0);
    });
  }
  if (profile_) stats_.flip_rhs_nnz.add(rhs_nnz);
  ProfileTimer measurement(time_iteration_, stats_, Component::kFlipFtran);
  factor_.ftran(column);
  measurement.switch_to(Component::kPrimalUpdate);
  for (Index k = 0; k < m_; ++k) {
    x_[basic_[k]] -= column[k];
    if (column[k] != 0.0) note_basic_change(k);
  }
  stats_.bound_flips += static_cast<long long>(flips_.size());
}

Index Simplex::choose_entering_column() const {
  Index best = -1;
  double best_score = 0.0;
  const double tol = options_.dual_tol;
  for (Index j = 0; j < nt_; ++j) {
    if (lower_[j] == upper_[j]) continue;
    double score = 0.0;
    switch (status_[j]) {
      case VarStatus::kAtLower: score = -d_[j]; break;
      case VarStatus::kAtUpper: score = d_[j]; break;
      case VarStatus::kAtZero: score = std::fabs(d_[j]); break;
      case VarStatus::kBasic: continue;
    }
    if (score <= tol) continue;
    score = score * score / devex_weight_[j];
    if (score > best_score) {
      best_score = score;
      best = j;
    }
  }
  return best;
}

SimplexStatus Simplex::dual_loop() {
  // A warm start that just rebuilt (or reused) the factorization has current values: skip another
  // factorization here. The optimal point is still confirmed on a fresh one below.
  if (skip_rebuild_) {
    skip_rebuild_ = false;
  } else if (!rebuild()) {
    return SimplexStatus::kNumericalError;
  }
  // Bounds may have changed since the last primal values (phase 1 boxes, a warm start).
  rebuild_infeasible_list();
  correct_dual_infeasibilities();
  bool fresh = true;
  const auto refresh = [&](RebuildReason reason) {
    if (!rebuild(reason)) return false;
    correct_dual_infeasibilities();
    fresh = true;
    return true;
  };

  for (;;) {
    SimplexStatus limit;
    if (limit_reached(limit)) return limit;
    if (factor_.should_refactor() && !refresh(factor_.refactor_reason())) {
      return SimplexStatus::kNumericalError;
    }

    time_iteration_ = sample_iteration();
    ProfileTimer measurement(time_iteration_, stats_, Component::kChuzr, nullptr, true);
    const Index r = choose_leaving_row();
    measurement.switch_to(Component::kRest);
    if (r < 0) {
      if (fresh) return SimplexStatus::kOptimal;
      // Confirm optimality on values recomputed from scratch with the current factorization
      // (its updates are bounded by the refactorization interval) instead of refactorizing:
      // that was one factorization per node LP.
      compute_primal();
      compute_dual();
      correct_dual_infeasibilities();
      fresh = true;
      continue;
    }
    const Index p = basic_[r];
    const bool to_lower = x_[p] < lower_[p];
    const double bound = to_lower ? lower_[p] : upper_[p];
    const double direction = to_lower ? -1.0 : 1.0;

    measurement.switch_to(Component::kBtran);
    rho_.assign(static_cast<std::size_t>(m_), 0.0);
    rho_[r] = 1.0;
    factor_.btran(rho_);
    measurement.switch_to(Component::kPivotRow);
    compute_pivot_row();
    measurement.switch_to(Component::kRest);

    measurement.switch_to(Component::kRatioTest);
    const Index q = dual_ratio_test(direction, std::fabs(x_[p] - bound));
    measurement.switch_to(Component::kRest);
    if (q < 0) {
      if (!fresh) {
        if (!refresh(RebuildReason::kNumerical)) return SimplexStatus::kNumericalError;
        continue;
      }
      dual_ray_ = rho_;
      return SimplexStatus::kInfeasible;
    }

    measurement.switch_to(Component::kEnteringFtran);
    load_column(q, alpha_col_);
    factor_.ftran(alpha_col_, &spike_, profile_ ? &stats_.alpha_col_nnz : nullptr);
    measurement.switch_to(Component::kRest);
    const double a_col = alpha_col_[r];
    const double a_row = alpha_row_[q];
    if (std::fabs(a_col - a_row) > kPivotConsistencyTol * (1.0 + std::fabs(a_col))) {
      if (!fresh) {
        ++stats_.rebuilds_on_mismatch;
        if (!refresh(RebuildReason::kMismatch)) return SimplexStatus::kNumericalError;
        continue;
      }
      log_.log(3, "simplex: pivot mismatch %.3e vs %.3e after refactor", a_col, a_row);
    }
    if (std::fabs(a_col) < 1e-11) return SimplexStatus::kNumericalError;

    measurement.switch_to(Component::kDseFtran);
    tau_ = rho_;
    factor_.ftran(tau_, nullptr, profile_ ? &stats_.tau_nnz : nullptr);
    measurement.switch_to(Component::kRest);

    // Bound flips from the long-step ratio test move x_B; the leaving variable stays infeasible
    // on the same side, so its distance to the bound is taken afterwards.
    if (!flips_.empty()) apply_flips();
    const double delta = x_[p] - bound;

    measurement.switch_to(Component::kDualUpdate);
    // Dual step. A slightly infeasible d_q would make the step go the wrong way; shift its cost
    // instead so the step is zero.
    double theta_d = d_[q] / a_row;
    if (direction * theta_d < 0.0) {
      ++stats_.cost_shifts;
      cost_[q] -= d_[q];
      theta_d = 0.0;
    }
    if (theta_d != 0.0) {
      for (const Index j : row_nz_) {
        if (status_[j] != VarStatus::kBasic) d_[j] -= theta_d * alpha_row_[j];
      }
    }
    d_[q] = 0.0;
    d_[p] = -theta_d;

    measurement.switch_to(Component::kPrimalUpdate);
    // Primal step: p moves exactly to the violated bound.
    const double theta_p = delta / a_col;
    for (Index k = 0; k < m_; ++k) {
      x_[basic_[k]] -= theta_p * alpha_col_[k];
      if (alpha_col_[k] != 0.0) note_basic_change(k);
    }
    x_[q] += theta_p;
    x_[p] = bound;

    measurement.switch_to(Component::kDseUpdate);
    // Dual steepest-edge weights (Forrest–Goldfarb update). The pivot row's weight is known
    // exactly as ||rho||^2; using it instead of the updated value keeps rounding errors from
    // accumulating in the recurrence.
    double w_r = 0.0;
    for (const double v : rho_) w_r += v * v;
    for (Index k = 0; k < m_; ++k) {
      if (k == r || alpha_col_[k] == 0.0) continue;
      const double ratio = alpha_col_[k] / a_col;
      dse_weight_[k] =
          std::max(dse_weight_[k] + ratio * (ratio * w_r - 2.0 * tau_[k]), kMinDseWeight);
    }
    dse_weight_[r] = std::max(w_r / (a_col * a_col), kMinDseWeight);
    measurement.switch_to(Component::kRest);

    basic_[r] = q;
    position_[q] = r;
    status_[q] = VarStatus::kBasic;
    position_[p] = -1;
    status_[p] = (to_lower || lower_[p] == upper_[p]) ? VarStatus::kAtLower : VarStatus::kAtUpper;
    note_basic_change(r);
    ++iterations_;
    ++stats_.dual_iterations;
    if (profile_) {
      if (dual_phase1_) ++stats_.dual_phase1_iterations;
      else ++stats_.dual_phase2_iterations;
    }
    if (update_factor(r)) {
      fresh = false;
    } else if (!refresh(RebuildReason::kUnsafeUpdate)) {
      return SimplexStatus::kNumericalError;
    }
  }
}

SimplexStatus Simplex::primal_loop(LoopResult& result) {
  result = LoopResult::kDone;
  if (!rebuild()) return SimplexStatus::kNumericalError;
  // Devex reference framework: the current nonbasic variables.
  std::fill(devex_weight_.begin(), devex_weight_.end(), 1.0);
  bool fresh = true;
  const auto refresh = [&](RebuildReason reason) {
    if (!rebuild(reason)) return false;
    fresh = true;
    return true;
  };
  const double ptol = options_.primal_tol;
  // Rounding drift of a few tolerances is absorbed by the Harris ratio test during the loop; only
  // the final point must be feasible within the tolerance (the dual simplex repairs it if not).
  const double drift_tol = 100.0 * ptol;
  if (max_primal_infeasibility() > drift_tol) {
    result = LoopResult::kLostFeasibility;
    return SimplexStatus::kOptimal;
  }

  for (;;) {
    SimplexStatus limit;
    if (limit_reached(limit)) return limit;
    if (factor_.should_refactor()) {
      if (!refresh(factor_.refactor_reason())) return SimplexStatus::kNumericalError;
      if (max_primal_infeasibility() > drift_tol) {
        result = LoopResult::kLostFeasibility;
        return SimplexStatus::kOptimal;
      }
    }

    time_iteration_ = sample_iteration();
    ProfileTimer measurement(time_iteration_, stats_, Component::kRest, nullptr, true);
    const Index q = choose_entering_column();
    if (q < 0) {
      if (!fresh) {
        if (!refresh(RebuildReason::kConfirmation)) return SimplexStatus::kNumericalError;
        continue;
      }
      if (max_primal_infeasibility() > ptol) result = LoopResult::kLostFeasibility;
      return SimplexStatus::kOptimal;
    }
    const double dir = d_[q] < 0.0 ? 1.0 : -1.0;
    measurement.switch_to(Component::kEnteringFtran);
    load_column(q, alpha_col_);
    factor_.ftran(alpha_col_, &spike_, profile_ ? &stats_.alpha_col_nnz : nullptr);
    measurement.switch_to(Component::kRest);

    measurement.switch_to(Component::kRatioTest);
    // Harris two-pass ratio test on x_B(t) = x_B - t * dir * alpha.
    double max_step = kInf;
    for (Index k = 0; k < m_; ++k) {
      const double a = dir * alpha_col_[k];
      if (std::fabs(a) < options_.pivot_tol) continue;
      const Index j = basic_[k];
      if (a > 0.0 && lower_[j] > -kInf) {
        max_step = std::min(max_step, (x_[j] - lower_[j] + ptol) / a);
      } else if (a < 0.0 && upper_[j] < kInf) {
        max_step = std::min(max_step, (x_[j] - upper_[j] - ptol) / a);
      }
    }
    Index r = -1;
    double step = kInf;
    double best_abs = 0.0;
    if (max_step < kInf) {
      for (Index k = 0; k < m_; ++k) {
        const double a = dir * alpha_col_[k];
        if (std::fabs(a) < options_.pivot_tol) continue;
        const Index j = basic_[k];
        double ratio;
        if (a > 0.0 && lower_[j] > -kInf) {
          ratio = (x_[j] - lower_[j]) / a;
        } else if (a < 0.0 && upper_[j] < kInf) {
          ratio = (x_[j] - upper_[j]) / a;
        } else {
          continue;
        }
        if (ratio <= max_step && std::fabs(a) > best_abs) {
          best_abs = std::fabs(a);
          r = k;
          step = std::max(ratio, 0.0);
        }
      }
    }
    const double range = upper_[q] - lower_[q];
    measurement.switch_to(Component::kRest);

    if (r < 0 && !(range < kInf)) {
      if (!fresh) {
        if (!refresh(RebuildReason::kNumerical)) return SimplexStatus::kNumericalError;
        continue;
      }
      primal_ray_.assign(static_cast<std::size_t>(nt_), 0.0);
      primal_ray_[q] = dir;
      for (Index k = 0; k < m_; ++k) primal_ray_[basic_[k]] = -dir * alpha_col_[k];
      return SimplexStatus::kUnbounded;
    }

    if (range <= step) {
      measurement.switch_to(Component::kPrimalUpdate);
      // Bound flip of the entering variable; the basis is unchanged.
      for (Index k = 0; k < m_; ++k) x_[basic_[k]] -= dir * range * alpha_col_[k];
      status_[q] = status_[q] == VarStatus::kAtLower ? VarStatus::kAtUpper : VarStatus::kAtLower;
      set_value_from_status(q);
      ++iterations_;
      ++stats_.primal_iterations;
      if (profile_) ++stats_.primal_cleanup_iterations;
      measurement.switch_to(Component::kRest);
      ++stats_.bound_flips;
      continue;
    }

    const Index p = basic_[r];
    const double a_col = alpha_col_[r];
    const bool p_to_lower = dir * a_col > 0.0;

    measurement.switch_to(Component::kBtran);
    rho_.assign(static_cast<std::size_t>(m_), 0.0);
    rho_[r] = 1.0;
    factor_.btran(rho_);
    measurement.switch_to(Component::kPivotRow);
    compute_pivot_row();
    measurement.switch_to(Component::kRest);
    const double a_row = alpha_row_[q];
    if (std::fabs(a_col - a_row) > kPivotConsistencyTol * (1.0 + std::fabs(a_col)) && !fresh) {
      ++stats_.rebuilds_on_mismatch;
      if (!refresh(RebuildReason::kMismatch)) return SimplexStatus::kNumericalError;
      continue;
    }

    measurement.switch_to(Component::kDualUpdate);
    const double theta_d = d_[q] / a_row;
    const double w_q = devex_weight_[q];
    for (const Index j : row_nz_) {
      if (status_[j] != VarStatus::kBasic) {
        d_[j] -= theta_d * alpha_row_[j];
        const double ratio = alpha_row_[j] / a_row;
        devex_weight_[j] = std::max(devex_weight_[j], ratio * ratio * w_q);
      }
    }
    d_[q] = 0.0;
    d_[p] = -theta_d;
    devex_weight_[p] = std::max(w_q / (a_row * a_row), 1.0);

    measurement.switch_to(Component::kPrimalUpdate);
    for (Index k = 0; k < m_; ++k) x_[basic_[k]] -= step * dir * alpha_col_[k];
    x_[q] += step * dir;
    x_[p] = p_to_lower ? lower_[p] : upper_[p];
    measurement.switch_to(Component::kRest);

    basic_[r] = q;
    position_[q] = r;
    status_[q] = VarStatus::kBasic;
    position_[p] = -1;
    status_[p] = (p_to_lower || lower_[p] == upper_[p]) ? VarStatus::kAtLower : VarStatus::kAtUpper;
    dse_weight_[r] = 1.0;
    ++iterations_;
    ++stats_.primal_iterations;
    if (profile_) ++stats_.primal_cleanup_iterations;
    if (update_factor(r)) {
      fresh = false;
    } else if (!refresh(RebuildReason::kUnsafeUpdate)) {
      return SimplexStatus::kNumericalError;
    }
  }
}

void Simplex::tableau_row(Index position, std::vector<double>& row) const {
  std::vector<double> rho(static_cast<std::size_t>(m_), 0.0);
  rho[position] = 1.0;
  factor_.btran(rho);
  row.assign(static_cast<std::size_t>(nt_), 0.0);
  const auto start = lp_.At.col_start();
  const auto index = lp_.At.row_index();
  const auto value = lp_.At.values();
  for (Index i = 0; i < m_; ++i) {
    const double r = rho[i];
    if (r == 0.0) continue;
    for (NnzIndex p = start[i]; p < start[i + 1]; ++p) row[index[p]] += r * value[p];
    row[n_ + i] = -r;
  }
}

std::vector<double> Simplex::exact_dse_weights() const {
  std::vector<double> weights(static_cast<std::size_t>(m_));
  std::vector<double> row;
  for (Index r = 0; r < m_; ++r) {
    row.assign(static_cast<std::size_t>(m_), 0.0);
    row[r] = 1.0;
    factor_.btran(row);
    double sum = 0.0;
    for (const double v : row) sum += v * v;
    weights[r] = sum;
  }
  return weights;
}

SimplexStatus Simplex::dual_phase1() {
  for (Index j = 0; j < nt_; ++j) {
    const bool has_lower = lp_.lower[j] > -kInf;
    const bool has_upper = lp_.upper[j] < kInf;
    if (has_lower && has_upper) {
      lower_[j] = upper_[j] = 0.0;
    } else if (has_lower) {
      lower_[j] = 0.0;
      upper_[j] = 1.0;
    } else if (has_upper) {
      lower_[j] = -1.0;
      upper_[j] = 0.0;
    } else {
      lower_[j] = -kPhase1FreeBound;
      upper_[j] = kPhase1FreeBound;
    }
  }
  place_nonbasic_dual_feasible();
  dual_phase1_ = true;
  SimplexStatus status = dual_loop();
  dual_phase1_ = false;
  lower_ = lp_.lower;
  upper_ = lp_.upper;
  cost_ = lp_.cost;
  log_.log(2, "simplex: dual phase 1 %s after %lld iterations", to_string(status), iterations_);
  // The auxiliary problem is always feasible (v = 0), so "infeasible" means numerical trouble.
  if (status == SimplexStatus::kInfeasible) return SimplexStatus::kNumericalError;
  if (status != SimplexStatus::kOptimal) return status;
  compute_dual();
  return SimplexStatus::kOptimal;
}

SimplexStatus Simplex::phase2() {
  for (int round = 0; round < 8; ++round) {
    SimplexStatus status = dual_loop();
    log_.log(2, "simplex: dual phase 2 %s after %lld iterations", to_string(status), iterations_);
    if (status != SimplexStatus::kOptimal) return status;

    // Remove perturbations and shifts; the basis stays primal feasible.
    if (cost_ != lp_.cost) {
      cost_ = lp_.cost;
      compute_dual();
    }
    if (max_dual_infeasibility() <= options_.dual_tol) return SimplexStatus::kOptimal;

    LoopResult result;
    status = primal_loop(result);
    log_.log(2, "simplex: primal cleanup %s after %lld iterations", to_string(status), iterations_);
    if (result == LoopResult::kLostFeasibility) {
      std::fill(dse_weight_.begin(), dse_weight_.end(), 1.0);
      continue;
    }
    return status;
  }
  return SimplexStatus::kNumericalError;
}

SimplexStatus Simplex::solve_dual_infeasible() {
  // Decide primal feasibility with a zero objective (trivially dual feasible) ... An all-zero
  // objective makes every ratio test a tie and the dual simplex stalls, so use small random costs
  // oriented to keep the starting basis dual feasible; the feasible region is unchanged.
  std::fill(cost_.begin(), cost_.end(), 0.0);
  compute_dual();
  place_nonbasic_dual_feasible();
  compute_primal();
  perturb_costs(kFeasibilityCostBase);
  dual_phase1_ = true;
  SimplexStatus status = dual_loop();
  dual_phase1_ = false;
  log_.log(2, "simplex: feasibility search %s after %lld iterations", to_string(status),
           iterations_);
  if (status != SimplexStatus::kOptimal) return status;

  // ... then the primal simplex from the feasible basis finds the unbounded ray.
  cost_ = lp_.cost;
  compute_dual();
  LoopResult result;
  status = primal_loop(result);
  if (result == LoopResult::kLostFeasibility) return phase2();
  return status;
}

void Simplex::reset() {
  timer_ = Timer();
  iterations_ = 0;
  cost_ = lp_.cost;
  lower_ = lp_.lower;
  upper_ = lp_.upper;
  std::fill(dse_weight_.begin(), dse_weight_.end(), 1.0);
  std::fill(devex_weight_.begin(), devex_weight_.end(), 1.0);
}

bool Simplex::update_factor(Index position) {
  const ProfileTimer measurement(time_iteration_, stats_, Component::kFactorUpdate);
  return factor_.update(position, spike_);
}

bool Simplex::sample_iteration() {
  if (!profile_) return false;
  const bool sample = stats_.timing_attempts++ % stats_.timing_sample_interval == 0;
  if (sample) ++stats_.timing_samples;
  return sample;
}

void Simplex::begin_profile(bool warm_start) {
  profile_ = options_.collect_stats || log_.enabled(2);
  warm_start_ = warm_start;
  warm_start_setup_ = warm_start;
  dual_phase1_ = false;
  // Keep the original lifetime work counters, while starting a fresh per-solve profile.
  SimplexStats fresh;
  fresh.dual_iterations = stats_.dual_iterations;
  fresh.primal_iterations = stats_.primal_iterations;
  fresh.bound_flips = stats_.bound_flips;
  fresh.cost_shifts = stats_.cost_shifts;
  fresh.refactorizations = stats_.refactorizations;
  fresh.rebuilds_on_mismatch = stats_.rebuilds_on_mismatch;
  stats_ = fresh;
  stats_.collected = profile_;
  stats_.timing_sample_interval = std::max(1, options_.timing_sample_interval);
  time_iteration_ = false;
  factor_.reset_stats();
  if (profile_) profile_start_ = std::chrono::steady_clock::now();
}

void Simplex::finish_profile() {
  if (!profile_) return;
  stats_.total_seconds = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - profile_start_).count();
  stats_.factor = factor_.stats();
  if (warm_start_) stats_.warm_start_iterations = iterations_;
  stats_.sampled_seconds = stats_.seconds;
  const auto rebuild_slot = static_cast<std::size_t>(Component::kRebuild);
  const auto rest_slot = static_cast<std::size_t>(Component::kRest);
  double sampled = stats_.seconds[rest_slot];
  for (std::size_t k = 0; k < rebuild_slot; ++k) sampled += stats_.seconds[k];
  double scale = stats_.timing_samples > 0
      ? static_cast<double>(stats_.timing_attempts) / static_cast<double>(stats_.timing_samples)
      : 1.0;
  // Deterministic stride sampling can overestimate (e.g. a very short/irregular solve). Keep
  // estimates within the unsampled wall-time budget, preserving component proportions.
  if (sampled > 0.0) {
    scale = std::min(scale, std::max(0.0, stats_.total_seconds - stats_.seconds[rebuild_slot]) /
                                sampled);
  }
  for (std::size_t k = 0; k < rebuild_slot; ++k) stats_.seconds[k] *= scale;
  stats_.seconds[rest_slot] = 0.0;
  double measured = 0.0;
  for (double seconds : stats_.seconds) measured += seconds;
  stats_.seconds[static_cast<std::size_t>(Component::kRest)] =
      std::max(0.0, stats_.total_seconds - measured);
  if (!log_.enabled(2)) return;

  log_.log(2, "simplex profile begin (m=%d n=%d iterations=%lld)", m_, n_, iterations_);
  log_.log(2, "  timing: every %d iteration attempts; sampled=%lld/%lld scale=%.6f; "
              "component estimates, rebuild unsampled, counters/densities exact",
           stats_.timing_sample_interval, stats_.timing_samples, stats_.timing_attempts, scale);
  log_.log(2, "  component                    seconds     share     ms/iter");
  const auto time_row = [&](const char* name, double seconds) {
    log_.log(2, "  %-27s %10.6f %8.2f%% %11.6f", name, seconds,
             stats_.total_seconds > 0.0 ? 100.0 * seconds / stats_.total_seconds : 0.0,
             iterations_ > 0 ? 1000.0 * seconds / static_cast<double>(iterations_) : 0.0);
  };
  constexpr const char* names[] = {
      "CHUZR", "BTRAN", "pivot row", "ratio test", "entering-column FTRAN", "DSE FTRAN (tau)",
      "bound-flip FTRAN", "primal update", "dual/reduced-cost update", "DSE weight update",
      "Forrest-Tomlin update", "rebuild + recompute", "rest"};
  for (std::size_t k = 0; k < stats_.seconds.size(); ++k) time_row(names[k], stats_.seconds[k]);
  time_row("total", stats_.total_seconds);
  log_.log(2, "  phases: dual2=%lld dual1/feasibility=%lld primal_cleanup=%lld",
           stats_.dual_phase2_iterations, stats_.dual_phase1_iterations,
           stats_.primal_cleanup_iterations);
  log_.log(2, "  warm start: iterations=%lld (overlaps phases; repair not separable) "
              "rank_repairs=%lld flips=%lld shifts=%lld",
           stats_.warm_start_iterations, stats_.warm_start_rank_repairs,
           stats_.warm_start_bound_flips, stats_.warm_start_cost_shifts);
  log_.log(2, "  density/vector                 solves      nnz sum      nnz max    avg nnz/m");
  const auto density_row = [&](const char* name, const MeasurementCount& count) {
    const double density = count.count > 0 && m_ > 0
        ? static_cast<double>(count.sum) / static_cast<double>(count.count) / m_ : 0.0;
    log_.log(2, "  %-27s %10lld %12lld %12lld %12.6f", name, count.count, count.sum,
             count.maximum, density);
  };
  density_row("rho (after BTRAN)", stats_.rho_nnz);
  density_row("alpha_col (after FTRAN)", stats_.alpha_col_nnz);
  density_row("tau (after FTRAN)", stats_.tau_nnz);
  density_row("flip RHS (before FTRAN)", stats_.flip_rhs_nnz);
  const auto count_row = [&](const char* name, const MeasurementCount& count) {
    log_.log(2, "  %-27s count=%lld sum=%lld max=%lld avg=%.3f", name, count.count, count.sum,
             count.maximum, count.count > 0
                 ? static_cast<double>(count.sum) / static_cast<double>(count.count) : 0.0);
  };
  count_row("dual ratio candidates", stats_.ratio_candidates);
  count_row("dual Harris groups", stats_.harris_groups);
  count_row("CHUZR infeasibility list", stats_.chuzr_list);
  constexpr const char* reasons[] = {
      "initial", "phase", "interval", "fill_limit", "unsafe_update", "mismatch", "numerical",
      "confirmation", "rank_repair"};
  log_.log(2, "  factor rebuilds=%lld factorize_s=%.6f max_factorize_s=%.6f "
              "max_rebuild_recompute_s=%.6f rank_repairs=%lld",
           stats_.factor.rebuilds, stats_.factor.rebuild_seconds,
           stats_.factor.max_rebuild_seconds, stats_.max_rebuild_seconds, stats_.rank_repairs);
  log_.log(2, "  rebuild calls=%lld avg_rebuild_recompute_s=%.6f avg_factorize_s=%.6f",
           stats_.rebuild_calls, stats_.rebuild_calls > 0
               ? stats_.seconds[rebuild_slot] / static_cast<double>(stats_.rebuild_calls) : 0.0,
           stats_.factor.rebuilds > 0
               ? stats_.factor.rebuild_seconds / static_cast<double>(stats_.factor.rebuilds) : 0.0);
  for (std::size_t k = 0; k < stats_.factor.reasons.size(); ++k) {
    log_.log(2, "  rebuild reason %-12s %lld", reasons[k], stats_.factor.reasons[k]);
  }
  log_.log(2, "  factor nnz: stored L/R eta entries; U includes diagonal; pre-state includes retained updates");
  count_row("L before rebuild", stats_.factor.old_l_nnz);
  count_row("U before rebuild", stats_.factor.old_u_nnz);
  count_row("R before rebuild", stats_.factor.old_r_nnz);
  count_row("L after rebuild", stats_.factor.l_nnz);
  count_row("U after rebuild", stats_.factor.u_nnz);
  count_row("R after rebuild", stats_.factor.r_nnz);
  log_.log(2, "simplex profile end");
}

SimplexStatus Simplex::solve() {
  begin_profile(false);
  const SimplexStatus result = solve_cold();
  finish_profile();
  return result;
}

SimplexStatus Simplex::solve(const std::vector<VarStatus>& start) {
  begin_profile(true);
  const SimplexStatus result = solve_warm(start);
  finish_profile();
  return result;
}

// Keeps a nonbasic status consistent with the bounds.
VarStatus Simplex::normalized_status(Index j, VarStatus st) const {
  const bool has_lower = lower_[j] > -kInf;
  const bool has_upper = upper_[j] < kInf;
  if (st == VarStatus::kAtLower && !has_lower) {
    return has_upper ? VarStatus::kAtUpper : VarStatus::kAtZero;
  }
  if (st == VarStatus::kAtUpper && !has_upper) {
    return has_lower ? VarStatus::kAtLower : VarStatus::kAtZero;
  }
  if (st == VarStatus::kAtZero && (has_lower || has_upper)) {
    return has_lower ? VarStatus::kAtLower : VarStatus::kAtUpper;
  }
  return st;
}

SimplexStatus Simplex::solve_warm(const std::vector<VarStatus>& start) {
  if (start.size() != static_cast<std::size_t>(nt_) ||
      std::count(start.begin(), start.end(), VarStatus::kBasic) != m_) {
    return solve_cold();
  }
  reset();
  // Same statuses as the factorized basis (a child solved right after its parent): keep the
  // factorization and its updates, only the bounds and costs changed.
  if (factor_valid_) {
    bool same = true;
    for (Index j = 0; j < nt_ && same; ++j) same = normalized_status(j, start[j]) == status_[j];
    if (same) {
      for (Index j = 0; j < nt_; ++j) {
        if (status_[j] != VarStatus::kBasic) set_value_from_status(j);
      }
      compute_primal();
      compute_dual();
      correct_dual_infeasibilities();
      skip_rebuild_ = true;
      warm_start_setup_ = false;
      return phase2();
    }
  }
  factor_valid_ = false;
  Index k = 0;
  for (Index j = 0; j < nt_; ++j) {
    const VarStatus st = normalized_status(j, start[j]);
    status_[j] = st;
    if (st == VarStatus::kBasic) {
      basic_[k] = j;
      position_[j] = k++;
    } else {
      position_[j] = -1;
      set_value_from_status(j);
    }
  }
  if (!rebuild(RebuildReason::kInitial)) return solve_cold();
  correct_dual_infeasibilities();
  skip_rebuild_ = true;
  warm_start_setup_ = false;
  return phase2();
}

SimplexStatus Simplex::solve_cold() {
  warm_start_ = false;
  warm_start_setup_ = false;
  reset();
  factor_valid_ = false;

  // Slack basis.
  for (Index i = 0; i < m_; ++i) {
    basic_[i] = n_ + i;
    position_[n_ + i] = i;
    status_[n_ + i] = VarStatus::kBasic;
  }
  for (Index j = 0; j < n_; ++j) {
    position_[j] = -1;
    if (lower_[j] > -kInf) {
      status_[j] = VarStatus::kAtLower;
    } else if (upper_[j] < kInf) {
      status_[j] = VarStatus::kAtUpper;
    } else {
      status_[j] = VarStatus::kAtZero;
    }
    set_value_from_status(j);
  }
  if (!rebuild(RebuildReason::kInitial)) return SimplexStatus::kNumericalError;

  if (count_dual_infeasibilities_unboxed() > 0) {
    const SimplexStatus status = dual_phase1();
    if (status != SimplexStatus::kOptimal) return status;
    if (count_dual_infeasibilities_unboxed() > 0) return solve_dual_infeasible();
  }
  place_nonbasic_dual_feasible();
  compute_primal();
  if (options_.perturb) perturb_costs(kPerturbationBase);
  return phase2();
}

}  // namespace samaya
