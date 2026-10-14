#pragma once

#include <array>
#include <cstdint>
#include <random>
#include <vector>

#include "core/log.hpp"
#include "linalg/basis_factor.hpp"
#include "samaya/sparse_matrix.hpp"
#include "samaya/types.hpp"

namespace samaya {

// LP in computational form:  min cost'v  s.t.  [A  -I] v = 0,  lower <= v <= upper,
// where v = (x, r) holds the n structural variables followed by the m row activities
// ("logicals"). Logical n+i carries the bounds of row i.
struct LpProblem {
  Index m = 0;
  Index n = 0;
  SparseMatrix A;   // m x n.
  SparseMatrix At;  // Transpose of A, i.e. A in row-wise form.
  std::vector<double> cost;   // n + m.
  std::vector<double> lower;  // n + m.
  std::vector<double> upper;  // n + m.
};

enum class SimplexStatus : std::uint8_t {
  kOptimal,
  kInfeasible,
  kUnbounded,
  kIterationLimit,
  kTimeLimit,
  kNumericalError,
};

const char* to_string(SimplexStatus status);

enum class VarStatus : std::uint8_t { kBasic, kAtLower, kAtUpper, kAtZero };

// Work counters, mainly for tests and logs: a healthy dual simplex does nearly all of its work in
// dual iterations, with few cost shifts and little primal cleanup.
struct SimplexStats {
  // Push work is counted separately from simplex iterations, but shares their limit.
  long long push_superbasics = 0;
  long long push_targets = 0;
  long long push_pivots = 0;
  long long push_degenerate_pivots = 0;
  long long push_refactorizations = 0;
  long long push_primal_recomputations = 0;
  long long push_remaining = 0;
  long long push_infeasible = 0;
  long long push_absorbed = 0;           // Snaps absorbed into the frozen right-hand side.
  long long push_infeasible_blocks = 0;  // Step-0 blocks by an already infeasible basic.
  double push_max_rhs = 0.0;             // Largest residual left at the end of the push.
  long long push_cost_shifts = 0;        // Dual infeasibilities shifted away at the handoff.
  double push_seconds = 0.0;
  double push_max_infeasibility = 0.0;
  double push_primal_error = 0.0;  // Incremental versus final recomputed basics, relative.
  bool push_used = false;
  bool push_primal = false;
  bool push_dual_fallback = false;
  long long dual_iterations = 0;
  long long primal_iterations = 0;
  long long bound_flips = 0;
  long long cost_shifts = 0;
  long long refactorizations = 0;
  long long rebuilds_on_mismatch = 0;
  // Detailed measurements are per solve; the legacy work counters above retain their lifetime
  // semantics. A warm start has no separate repair phase: its iterations are an overlapping subset.
  enum class Component : std::uint8_t {
    kChuzr, kBtran, kPivotRow, kRatioTest, kEnteringFtran, kDseFtran, kFlipFtran,
    kPrimalUpdate, kDualUpdate, kDseUpdate, kFactorUpdate, kRebuild, kRest, kCount
  };
  bool collected = false;
  long long dual_phase2_iterations = 0;
  long long dual_phase1_iterations = 0;  // Includes the auxiliary feasibility search.
  long long primal_cleanup_iterations = 0;
  long long warm_start_iterations = 0;
  long long warm_start_rank_repairs = 0;
  long long warm_start_bound_flips = 0;
  long long warm_start_cost_shifts = 0;
  long long rank_repairs = 0;
  long long rebuild_calls = 0;  // A call can contain several factorization attempts/rank repairs.
  MeasurementCount rho_nnz, alpha_col_nnz, tau_nnz, flip_rhs_nnz;
  MeasurementCount ratio_candidates, harris_groups;  // Dual ratio-test calls, including retries.
  MeasurementCount chuzr_list;  // Infeasibility list length at each CHUZR, stale entries included.
  MeasurementCount primal_infeasible_rows;  // Infeasible basic rows at each CHUZR, either mode.
  long long chuzr_scans = 0;    // CHUZRs that scanned all m rows instead (dense updates).
  BasisFactor::Stats factor;
  std::array<double, static_cast<std::size_t>(Component::kCount)> seconds{};
  // Raw sampled component seconds; seconds above contains the estimates at solve completion.
  std::array<double, static_cast<std::size_t>(Component::kCount)> sampled_seconds{};
  long long timing_attempts = 0;
  long long timing_samples = 0;
  int timing_sample_interval = 1;
  double total_seconds = 0.0;
  double max_rebuild_seconds = 0.0;
};

struct SimplexOptions {
  double primal_tol = 1e-7;
  double dual_tol = 1e-7;
  double pivot_tol = 1e-7;
  long long max_iterations = -1;  // -1: automatic, proportional to the problem size.
  double time_limit = kInf;
  bool perturb = true;
  // Upper bound on updates between factorizations; the amortized rule (BasisFactor::Options::
  // rebuild_cost) usually refactorizes earlier. On the 467k-row LP intervals of 100, 300 and
  // 1000 gave 35.4, 28.7 and 26.4 ms per iteration with no numerical rebuilds (step 2 in
  // SIMPLEX_SPEED_PLAN.md).
  int refactor_interval = 2000;
  std::uint64_t seed = 1;
  bool collect_stats = false;  // Detailed counters/timing; also enabled automatically at log level 2.
  int timing_sample_interval = 32;  // 1: time every iteration; counters/densities are always exact.
  // CHUZR scans all m basic rows instead of the infeasibility list. The two choose the same row
  // (ties go to the lowest position in both); the scan is the reference of the equivalence tests.
  bool chuzr_full_scan = false;
  // The list is kept only while the primal updates since the last recomputation averaged at most
  // this fraction of m nonzeros; denser updates scan instead (both choose the same row). 1 keeps
  // the list always (tests); see Simplex::rebuild_infeasible_list for the measurements.
  double chuzr_list_max_density = 0.1;
};

// Bounded dual simplex with dual steepest-edge pricing, a Harris two-pass ratio test, cost
// perturbation and shifting against degeneracy, and a primal simplex for cleanup.
//
// Phase 1 uses the auxiliary-problem approach: bounds are replaced by boxes ([0,1] for
// lower-bounded, [-1,0] for upper-bounded, [-1000,1000] for free variables, [0,0] otherwise) and
// the dual simplex minimizes the sum of dual infeasibilities. If the original problem is dual
// infeasible, a zero-cost dual simplex decides primal feasibility and the primal simplex then
// proves unboundedness.
class Simplex {
 public:
  Simplex(const LpProblem& lp, const SimplexOptions& options, const Logger& log);

  SimplexStatus solve();
  // Starts from the given variable statuses (n + m, exactly m basic) instead of the slack basis:
  // bounds are restored, dual infeasibilities are removed by bound flips and cost shifts, and the
  // dual simplex continues. Falls back to solve() if the statuses do not form a usable basis.
  SimplexStatus solve(const std::vector<VarStatus>& start);
  // Preserve the n structural values, push nonbasics to bounds, then finish by simplex.
  SimplexStatus solve(const std::vector<VarStatus>& start, const std::vector<double>& values);

  // All vectors are in the (scaled) space of the LpProblem.
  const std::vector<double>& values() const { return x_; }         // n + m.
  const std::vector<double>& duals() const { return y_; }          // m.
  const std::vector<double>& reduced_costs() const { return d_; }  // n + m.
  const std::vector<VarStatus>& status() const { return status_; }
  // Variable basic at each basis position (m) after a solve.
  const std::vector<Index>& basic() const { return basic_; }
  // Row `position` of the simplex tableau, e_r' B^-1 [A -I], over all n + m variables (dense,
  // scaled space). Used to derive Gomory cuts.
  void tableau_row(Index position, std::vector<double>& row) const;
  // Farkas ray y (length m) after kInfeasible: 0 is outside the range of
  // y'A x - y'r over the variable bounds.
  const std::vector<double>& dual_ray() const { return dual_ray_; }
  // Improving direction (length n + m) after kUnbounded.
  const std::vector<double>& primal_ray() const { return primal_ray_; }
  long long iterations() const { return iterations_; }
  // Limits for the next solve (strong branching and tree search re-solve with small budgets).
  void set_iteration_limit(long long limit);
  void set_time_limit(double seconds) { options_.time_limit = seconds; }
  const SimplexStats& stats() const { return stats_; }

  // Dual steepest-edge weights as maintained by the updates, and recomputed from scratch as
  // ||e_r' B^-1||^2 (m btran solves; for tests).
  const std::vector<double>& dse_weights() const { return dse_weight_; }
  std::vector<double> exact_dse_weights() const;

 private:
  enum class LoopResult : std::uint8_t { kDone, kLostFeasibility };

  void reset();
  void begin_profile(bool warm_start);
  void finish_profile();
  SimplexStatus solve_cold();
  SimplexStatus solve_warm(const std::vector<VarStatus>& start);
  SimplexStatus solve_push(const std::vector<VarStatus>& start, const std::vector<double>& values);
  SimplexStatus push_to_vertex();
  VarStatus nearest_status(Index j) const;
  double push_target(Index j) const;
  bool rebuild(BasisFactor::RebuildReason reason = BasisFactor::RebuildReason::kPhase,
               bool recompute = true);
  bool update_factor(Index position);
  bool sample_iteration();
  VarStatus normalized_status(Index j, VarStatus st) const;
  void compute_primal();
  void compute_dual();
  void compute_pivot_row();
  void set_value_from_status(Index j);
  void place_nonbasic_dual_feasible();
  // Removes dual infeasibilities by flipping boxed nonbasics to their other bound and shifting
  // the costs of the others; flip_boxed = false shifts every cost, keeping the primal values.
  bool correct_dual_infeasibilities(bool flip_boxed = true);
  void perturb_costs(double base);
  Index count_dual_infeasibilities_unboxed() const;
  double max_primal_infeasibility() const;
  double max_dual_infeasibility() const;
  bool limit_reached(SimplexStatus& status) const;

  SimplexStatus dual_phase1();
  SimplexStatus phase2();
  SimplexStatus solve_dual_infeasible();
  SimplexStatus dual_loop();
  SimplexStatus primal_loop(LoopResult& result);

  Index choose_leaving_row();
  bool primal_infeasible(Index j) const;
  void rebuild_infeasible_list();
  void note_basic_change(Index position);
  // x_[basic_[k]] -= scale * delta[k] for every position k; in list mode, lists the changed ones
  // that become primal infeasible. Counts the update's nonzeros for the mode choice.
  void update_basic_values(const std::vector<double>& delta, double scale);
  Index dual_ratio_test(double direction, double slope);
  void apply_flips();
  Index choose_entering_column() const;

  template <typename F>
  void for_column(Index j, F&& f) const;
  void load_column(Index j, std::vector<double>& out) const;

  const LpProblem& lp_;
  SimplexOptions options_;
  const Logger& log_;
  Index m_;
  Index n_;
  Index nt_;
  long long max_iterations_;
  Timer timer_;

  BasisFactor factor_;
  // factor_ holds the current basic_ (with its updates). skip_rebuild_: set by a warm start whose
  // values are current, so the first dual loop of phase 2 does not factorize again.
  bool factor_valid_ = false;
  bool skip_rebuild_ = false;
  bool pushing_ = false;  // Rank repair preserves leaving values during the push.
  // Right-hand side of [A -I] v = push_rhs_ during the push (the start point's residual plus the
  // absorbed snaps); empty otherwise.
  std::vector<double> push_rhs_;
  std::vector<Index> basic_;
  std::vector<Index> position_;
  std::vector<VarStatus> status_;
  std::vector<double> x_;
  std::vector<double> y_;
  std::vector<double> d_;
  std::vector<double> cost_;
  std::vector<double> lower_;
  std::vector<double> upper_;
  std::vector<double> dse_weight_;
  std::vector<double> devex_weight_;  // Primal Devex reference weights (n + m).

  std::vector<double> rho_;
  std::vector<double> tau_;
  std::vector<double> alpha_row_;
  std::vector<Index> row_nz_;       // Indices of the nonzeros of alpha_row_.
  std::vector<char> in_row_nz_;
  std::vector<Index> candidates_;   // Ratio test workspace.
  std::vector<Index> flips_;        // Boxed variables flipped by the long-step ratio test.
  // Basis positions that may be primal infeasible: every infeasible position is listed (entries
  // can be stale-feasible; CHUZR drops them). Rebuilt with the primal values, then kept up to
  // date for each position whose value changes in the dual loop.
  std::vector<Index> infeasible_;
  std::vector<char> in_infeasible_;
  bool track_infeasible_ = true;  // List mode; otherwise CHUZR scans all m rows.
  long long update_entries_ = 0;  // Positions and nonzeros of the primal updates since the list
  long long update_nonzeros_ = 0;  // was last rebuilt, to choose the mode.
  std::vector<double> alpha_col_;
  std::vector<double> spike_;
  std::vector<double> dual_ray_;
  std::vector<double> primal_ray_;

  long long iterations_ = 0;
  SimplexStats stats_;
  bool profile_ = false;
  bool time_iteration_ = false;
  bool warm_start_ = false;
  bool warm_start_setup_ = false;
  bool dual_phase1_ = false;
  std::chrono::steady_clock::time_point profile_start_{};
  std::mt19937_64 rng_;
};

}  // namespace samaya
