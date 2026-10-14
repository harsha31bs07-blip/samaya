#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>

#include "lp/crossover.hpp"
#include "lp/ipm.hpp"
#include "lp/lp_solver.hpp"
#include "lp_generators.hpp"
#include "reference_lp.hpp"
#include "samaya/verify.hpp"
#include "test_framework.hpp"

using namespace samaya;

namespace {

// Compare to the independent reference at the existing crossover objective tolerance.
constexpr double kObjectiveTolerance = 1e-6;
// Incremental push arithmetic must agree with a fresh solve to this relative tolerance.
constexpr double kPushTolerance = 1e-9;

Model single_column(double lower, double upper, double row_lower, double row_upper) {
  Model model;
  model.A = SparseMatrix::from_triplets(1, 1, {{0, 0, 1.0}});
  model.obj = {0.0};
  model.col_lower = {lower};
  model.col_upper = {upper};
  model.col_type = {VarType::kContinuous};
  model.row_lower = {row_lower};
  model.row_upper = {row_upper};
  return model;
}

void check_optimum(const Model& model, const LpResult& answer, double objective) {
  REQUIRE(answer.status == SimplexStatus::kOptimal);
  CHECK(verify_lp_optimality(model, answer.col_value, answer.row_dual).ok);
  CHECK_NEAR(answer.objective, objective, kObjectiveTolerance * (1.0 + std::fabs(objective)));
}

LpResult push_model(const Model& model, const std::vector<VarStatus>& start,
                    const std::vector<double>& values) {
  const auto ref = test::ReferenceLp(model).solve();
  REQUIRE(ref.status == test::ReferenceResult::Status::kOptimal);
  LpSolveOptions options;
  options.scale = false;
  options.start = start;
  options.start_values = values;
  const LpResult answer = solve_lp(model, options, Logger(0));
  check_optimum(model, answer, ref.objective);
  CHECK(answer.stats.push_used);
  CHECK_EQ(answer.stats.push_remaining, 0);
  CHECK(answer.stats.push_primal_error <= kPushTolerance);
  // The push starts from the point's own values (no recomputation) and recomputes with the frozen
  // right-hand side at least once, to check the incremental values at its end.
  CHECK(answer.stats.push_primal_recomputations >= 1);
  CHECK(answer.stats.push_primal);
  CHECK(!answer.stats.push_dual_fallback);
  return answer;
}

}  // namespace

TEST(crossover_push_reaches_bound_without_pivot) {
  const auto answer = push_model(single_column(0.0, 10.0, 0.0, 20.0),
      {VarStatus::kAtLower, VarStatus::kBasic}, {3.0});
  CHECK_EQ(answer.stats.push_superbasics, 1);
  CHECK_EQ(answer.stats.push_targets, 1);
  CHECK_EQ(answer.stats.push_pivots, 0);
}

TEST(crossover_push_blocked_by_basic) {
  Model model;
  model.A = SparseMatrix::from_triplets(1, 2, {{0, 0, 1.0}, {0, 1, 1.0}});
  model.obj = {0.0, 0.0};
  model.col_lower = {1.0, 0.0};
  model.col_upper = {10.0, 10.0};
  model.col_type.assign(2, VarType::kContinuous);
  model.row_lower = model.row_upper = {10.0};
  const auto answer = push_model(model,
      {VarStatus::kBasic, VarStatus::kAtUpper, VarStatus::kAtLower}, {3.0, 7.0});
  CHECK_EQ(answer.stats.push_pivots, 1);
  CHECK_EQ(answer.stats.push_targets, 0);
  CHECK_EQ(answer.stats.push_infeasible, 0);
  CHECK_NEAR(answer.col_value[0], 1.0, kPushTolerance);
  CHECK_NEAR(answer.col_value[1], 9.0, kPushTolerance);
}

TEST(crossover_push_zero_step_fixed_logical) {
  const auto answer = push_model(single_column(0.0, 10.0, 4.0, 4.0),
      {VarStatus::kAtLower, VarStatus::kBasic}, {4.0});
  CHECK_EQ(answer.stats.push_pivots, 1);
  CHECK_EQ(answer.stats.push_degenerate_pivots, 1);
}

TEST(crossover_push_order_prefers_distance_then_index) {
  Model model;
  model.A = SparseMatrix::from_triplets(1, 2, {{0, 0, 1.0}, {0, 1, 1.0}});
  model.obj = {0.0, 0.0};
  model.col_lower = {0.0, 0.0};
  model.col_upper = {10.0, 10.0};
  model.col_type.assign(2, VarType::kContinuous);
  model.row_lower = model.row_upper = {8.0};
  const std::vector<VarStatus> start = {VarStatus::kAtLower, VarStatus::kAtLower,
                                        VarStatus::kBasic};
  const auto tied = push_model(model, start, {4.0, 4.0});
  CHECK_NEAR(tied.col_value[0], 8.0, kPushTolerance);
  CHECK_NEAR(tied.col_value[1], 0.0, kPushTolerance);
  model.row_lower = model.row_upper = {7.0};
  const auto unequal = push_model(model, start, {3.0, 4.0});
  CHECK_NEAR(unequal.col_value[0], 0.0, kPushTolerance);
  CHECK_NEAR(unequal.col_value[1], 7.0, kPushTolerance);
}

TEST(crossover_push_free_superbasic) {
  const auto answer = push_model(single_column(-kInf, kInf, -kInf, kInf),
      {VarStatus::kAtZero, VarStatus::kBasic}, {3.0});
  CHECK_EQ(answer.stats.push_targets, 1);
  CHECK_NEAR(answer.col_value[0], 0.0, kPushTolerance);
}

TEST(crossover_push_upper_only) {
  Model model = single_column(-kInf, 10.0, -100.0, 100.0);
  model.obj = {-1.0};
  const auto answer = push_model(model, {VarStatus::kAtUpper, VarStatus::kBasic}, {3.0});
  CHECK_EQ(answer.stats.push_targets, 1);
  CHECK_NEAR(answer.col_value[0], 10.0, kPushTolerance);
}

TEST(crossover_push_superbasic_logical) {
  const auto answer = push_model(single_column(0.0, 10.0, 0.0, 20.0),
      {VarStatus::kBasic, VarStatus::kAtLower}, {3.0});
  CHECK_EQ(answer.stats.push_superbasics, 1);
  CHECK_EQ(answer.stats.push_targets, 1);
}

TEST(crossover_push_clamps_outside_point) {
  // Below the simplex tolerance, as in an approximate PDLP point.
  constexpr double kOutside = 1e-8;
  const auto answer = push_model(single_column(0.0, 10.0, 0.0, 20.0),
      {VarStatus::kAtLower, VarStatus::kBasic}, {-kOutside});
  CHECK_EQ(answer.stats.push_superbasics, 0);
  CHECK_NEAR(answer.col_value[0], 0.0, kPushTolerance);
}

TEST(crossover_push_primal_handoff_preserves_original_costs_and_point) {
  Model model = single_column(0.0, 10.0, 0.0, 20.0);
  model.obj = {-1.0};
  const auto answer = push_model(model, {VarStatus::kAtLower, VarStatus::kBasic}, {3.0});
  CHECK_EQ(answer.stats.push_targets, 1);
  // Dual correction would flip x from 0 to 10, hiding this required primal iteration.
  CHECK_EQ(answer.stats.primal_iterations, 1);
  CHECK_EQ(answer.stats.dual_iterations, 0);
  CHECK_EQ(answer.stats.cost_shifts, 0);
}

TEST(crossover_push_rank_repair_keeps_point) {
  Model model;
  model.A = SparseMatrix::from_triplets(2, 2,
      {{0, 0, 1.0}, {1, 0, 1.0}, {0, 1, 1.0}, {1, 1, 1.0}});
  model.obj = {1.0, 1.0};
  model.col_lower = {0.0, 0.0};
  model.col_upper = {10.0, 10.0};
  model.col_type.assign(2, VarType::kContinuous);
  model.row_lower = {18.0, 0.0};
  model.row_upper = {18.0, 18.0};
  const auto answer = push_model(model,
      {VarStatus::kBasic, VarStatus::kBasic, VarStatus::kAtLower, VarStatus::kAtUpper}, {9.0, 9.0});
  CHECK_EQ(answer.stats.push_superbasics, 1);
  CHECK_EQ(answer.stats.push_targets + answer.stats.push_pivots, 1);
}

TEST(crossover_push_limits_and_empty_values) {
  const Model model = single_column(0.0, 10.0, 4.0, 4.0);
  const LpProblem lp = make_problem(model, Scaling::identity(1, 1));
  const Logger quiet(0);
  const std::vector<VarStatus> start = {VarStatus::kAtLower, VarStatus::kBasic};
  SimplexOptions options;
  options.max_iterations = 0;
  Simplex limited(lp, options, quiet);
  CHECK(limited.solve(start, {4.0}) == SimplexStatus::kIterationLimit);
  options.max_iterations = -1;
  options.time_limit = 0.0;
  Simplex timed(lp, options, quiet);
  CHECK(timed.solve(start, {4.0}) == SimplexStatus::kTimeLimit);
  options.time_limit = kInf;
  Simplex old(lp, options, quiet), empty(lp, options, quiet);
  CHECK(old.solve(start) == empty.solve(start, {}));
  CHECK(old.values() == empty.values());
  CHECK(old.status() == empty.status());
  CHECK_EQ(old.iterations(), empty.iterations());
  CHECK(!empty.stats().push_used);
}

TEST(crossover_push_dual_fallback_for_infeasible_handoff) {
  const Model model = single_column(0.0, 10.0, 4.0, 4.0);
  LpSolveOptions options;
  options.scale = false;
  options.start = {VarStatus::kAtLower, VarStatus::kBasic};
  options.start_values = {0.0};
  const auto answer = solve_lp(model, options, Logger(0));
  const auto ref = test::ReferenceLp(model).solve();
  REQUIRE(ref.status == test::ReferenceResult::Status::kOptimal);
  check_optimum(model, answer, ref.objective);
  CHECK_EQ(answer.stats.push_superbasics, 0);
  CHECK_EQ(answer.stats.push_infeasible, 1);
  CHECK_NEAR(answer.stats.push_max_infeasibility, 4.0, kPushTolerance);
  CHECK(!answer.stats.push_primal);
  CHECK(answer.stats.push_dual_fallback);
}

TEST(crossover_push_dual_fallback_shifts_instead_of_flipping) {
  // The handoff is infeasible (row x0 = 4 with x0 at 0), so the dual simplex takes over. x1 (no
  // row entries, cost -1) sits at its lower bound with a wrong-sign reduced cost: the fallback
  // must shift its cost and keep the pushed point, not flip x1 to its upper bound.
  Model model;
  model.A = SparseMatrix::from_triplets(1, 2, {{0, 0, 1.0}});
  model.obj = {0.0, -1.0};
  model.col_lower = {0.0, 0.0};
  model.col_upper = {10.0, 10.0};
  model.col_type.assign(2, VarType::kContinuous);
  model.row_lower = model.row_upper = {4.0};
  LpSolveOptions options;
  options.scale = false;
  options.start = {VarStatus::kAtLower, VarStatus::kAtLower, VarStatus::kBasic};
  options.start_values = {0.0, 0.0};
  const auto answer = solve_lp(model, options, Logger(0));
  const auto ref = test::ReferenceLp(model).solve();
  REQUIRE(ref.status == test::ReferenceResult::Status::kOptimal);
  check_optimum(model, answer, ref.objective);
  CHECK(answer.stats.push_dual_fallback);
  CHECK(answer.stats.cost_shifts >= 1);
}

TEST(crossover_push_refactor_interval) {
  const Model model = single_column(0.0, 10.0, 4.0, 4.0);
  LpSolveOptions options;
  options.scale = false;
  options.simplex.refactor_interval = 1;
  options.start = {VarStatus::kAtLower, VarStatus::kBasic};
  options.start_values = {4.0};
  const auto answer = solve_lp(model, options, Logger(0));
  const auto ref = test::ReferenceLp(model).solve();
  REQUIRE(ref.status == test::ReferenceResult::Status::kOptimal);
  check_optimum(model, answer, ref.objective);
  CHECK_EQ(answer.stats.push_pivots, 1);
  CHECK_EQ(answer.stats.push_refactorizations, 3);
  CHECK_EQ(answer.stats.push_remaining, 0);
  CHECK(answer.stats.push_primal_error <= kPushTolerance);
}

TEST(crossover_push_starts_from_the_point_not_the_solved_basics) {
  Model model;
  model.A = SparseMatrix::from_triplets(1, 2, {{0, 0, 1.0}, {0, 1, 1.0}});
  model.obj = {0.0, 0.0};
  model.col_lower = {1.0, 0.0};
  model.col_upper = {10.0, 10.0};
  model.col_type.assign(2, VarType::kContinuous);
  model.row_lower = model.row_upper = {10.0};
  // The point's activity is 10.5, clamped to 10: the residual 0.5 is frozen. Solving the basis
  // would put x0 at 0.5, below its bound; the push keeps the point's x0 = 1 instead, so x1's push
  // toward 10 is a feasible zero step: x0 leaves at its bound and x1 enters.
  const auto answer = push_model(model,
      {VarStatus::kBasic, VarStatus::kAtUpper, VarStatus::kAtLower}, {1.0, 9.5});
  CHECK_EQ(answer.stats.push_degenerate_pivots, 1);
  CHECK_EQ(answer.stats.push_infeasible_blocks, 0);
  CHECK_EQ(answer.stats.push_targets, 0);
  CHECK_NEAR(answer.stats.push_max_rhs, 0.5, kPushTolerance);
  CHECK_NEAR(answer.col_value[0], 1.0, kPushTolerance);
  CHECK_NEAR(answer.col_value[1], 9.0, kPushTolerance);
}

TEST(crossover_push_absorbs_small_snaps_on_an_ill_conditioned_basis) {
  // Basis {x0, x1} with columns (1, 1) and (1, 1 + 1e-6): B^-1 has entries near 1e6. x2 is 5e-8
  // from its bound, within the tolerance; snapping it by moving the basics would shift them by
  // about 0.05, out of [0.99, 1.01], and every later push would be blocked by infeasible basics.
  // Absorbed into the frozen right-hand side, it leaves them in place, x3's push pivots, and the
  // residual only meets the final, well-conditioned basis.
  constexpr double kNearDependence = 1e-6;
  Model model;
  model.A = SparseMatrix::from_triplets(2, 4,
      {{0, 0, 1.0}, {0, 1, 1.0}, {1, 0, 1.0}, {1, 1, 1.0 + kNearDependence}, {1, 2, 1.0},
       {1, 3, 1.0}});
  model.obj.assign(4, 0.0);
  model.col_lower = {0.99, 0.99, 0.0, 0.0};
  model.col_upper = {1.01, 1.01, 10.0, 10.0};
  model.col_type.assign(4, VarType::kContinuous);
  const std::vector<double> point = {1.0, 1.0, 5e-8, 5.0};
  model.row_lower.assign(2, 0.0);
  model.A.multiply(point, model.row_lower);
  model.row_upper = model.row_lower;
  const auto answer = push_model(model,
      {VarStatus::kBasic, VarStatus::kBasic, VarStatus::kAtLower, VarStatus::kAtLower,
       VarStatus::kAtLower, VarStatus::kAtLower}, point);
  CHECK_EQ(answer.stats.push_superbasics, 1);
  CHECK(answer.stats.push_absorbed >= 1);
  CHECK(answer.stats.push_pivots >= 1);
  CHECK_EQ(answer.stats.push_infeasible_blocks, 0);
}

TEST(crossover_push_freezes_a_clamped_activity) {
  // Activity 5.00001 above the row's upper bound 5: the logical starts clamped at 5 and the
  // excess 1e-5 is the frozen residual. The push of x toward 10 is blocked at once by the
  // logical; on the true equations x = 5, which is optimal for the objective -x.
  Model model = single_column(0.0, 10.0, 0.0, 5.0);
  model.obj = {-1.0};
  constexpr double kExcess = 1e-5;
  const auto answer = push_model(model, {VarStatus::kAtLower, VarStatus::kBasic},
                                 {5.0 + kExcess});
  CHECK_EQ(answer.stats.push_superbasics, 1);
  CHECK_NEAR(answer.stats.push_max_rhs, kExcess, 1e-12);
  CHECK_NEAR(answer.col_value[0], 5.0, kPushTolerance);
}

TEST(crossover_push_interior_face_incremental_invariants) {
  std::mt19937 rng(616);
  std::uniform_int_distribution<int> entry(-3, 3);
  long long pivots = 0;
  for (int trial = 0; trial < 100; ++trial) {
    Model model;
    constexpr Index kRows = 4;
    constexpr Index kCols = 10;
    std::vector<Triplet> entries;
    for (Index i = 0; i < kRows; ++i) {
      for (Index j = 0; j < kCols; ++j) {
        const double value = i == j ? 12.0 : static_cast<double>(entry(rng));
        entries.push_back({i, j, value});
      }
    }
    model.A = SparseMatrix::from_triplets(kRows, kCols, std::move(entries));
    model.obj.assign(kCols, 0.0);
    model.col_lower.assign(kCols, 0.0);
    model.col_upper.assign(kCols, 10.0);
    // Narrow basic bounds force several blocking pivots in each dense basis.
    for (Index j = 0; j < kRows; ++j) {
      model.col_lower[j] = 3.0;
      model.col_upper[j] = 5.0;
    }
    model.col_type.assign(kCols, VarType::kContinuous);
    std::vector<double> point(kCols, 4.0);
    model.row_lower.assign(kRows, 0.0);
    model.A.multiply(point, model.row_lower);
    model.row_upper = model.row_lower;
    std::vector<VarStatus> start(kCols + kRows, VarStatus::kAtLower);
    for (Index j = 0; j < kRows; ++j) start[j] = VarStatus::kBasic;
    const auto answer = push_model(model, start, point);
    pivots += answer.stats.push_pivots;
    CHECK_EQ(answer.stats.push_superbasics, kCols - kRows);
    CHECK_EQ(answer.stats.push_targets + answer.stats.push_pivots, kCols - kRows);
  }
  CHECK(pivots > 50);
}

TEST(crossover_push_differential_reference_perturbed_and_interior) {
  std::mt19937 rng(615);
  const Logger quiet(0);
  // Approximate primal points exercise the fallback as well as near-feasible handoffs.
  constexpr double kPerturbation = 1e-5;
  std::uniform_real_distribution<double> perturb(-kPerturbation, kPerturbation);
  long long pushed_iterations = 0, snapped_iterations = 0, push_work = 0;
  int checked = 0, interiors = 0;
  for (const auto family : {test::LpFamily::kFeasible, test::LpFamily::kDegenerate,
                            test::LpFamily::kBoxed, test::LpFamily::kLoose}) {
    for (int k = 0; k < 60; ++k) {
      Model model = test::random_lp(family, 12, 12, rng);
      for (const auto sense : {ObjSense::kMinimize, ObjSense::kMaximize}) {
        model.sense = sense;
        const auto ref = test::ReferenceLp(model).solve();
        if (ref.status != test::ReferenceResult::Status::kOptimal) continue;
        const auto optimum = solve_lp(model, {}, quiet);
        REQUIRE(optimum.status == SimplexStatus::kOptimal);
        const auto ipm = solve_ipm(model, {}, quiet);
        for (int point = 0; point < 3; ++point) {
          if (point == 2 && ipm.status != IpmStatus::kOptimal) continue;
          std::vector<double> x = point == 2 ? ipm.x : ref.x;
          const auto& y = point == 2 ? ipm.y : optimum.row_dual;
          if (point == 1) {
            for (double& value : x) value += perturb(rng) * (1.0 + std::fabs(value));
          }
          for (Index j = 0; j < model.num_cols(); ++j) {
            x[j] = std::clamp(x[j], model.col_lower[j], model.col_upper[j]);
          }
          for (const bool scale : {false, true}) {
            LpSolveOptions options;
            options.scale = scale;
            options.start = crossover_basis(model, x, y);
            const auto snapped = solve_lp(model, options, quiet);
            options.start_values = x;
            const auto pushed = solve_lp(model, options, quiet);
            check_optimum(model, pushed, ref.objective);
            check_optimum(model, snapped, ref.objective);
            CHECK_EQ(pushed.stats.push_remaining, 0);
            CHECK(pushed.stats.push_primal_error <= kPushTolerance);
            pushed_iterations += pushed.iterations;
            snapped_iterations += snapped.iterations;
            push_work += pushed.stats.push_targets + pushed.stats.push_pivots;
            ++checked;
            if (point == 2) ++interiors;
          }
        }
      }
    }
  }
  std::printf("  push differential: %d solves (%d interior), simplex after push %lld (avg %.6f), "
              "statuses-only %lld (avg %.6f), push work %lld\n", checked, interiors,
              pushed_iterations, static_cast<double>(pushed_iterations) / checked,
              snapped_iterations, static_cast<double>(snapped_iterations) / checked, push_work);
  CHECK(checked > 1000);
  CHECK(interiors > 300);
  CHECK(pushed_iterations < snapped_iterations);
}
