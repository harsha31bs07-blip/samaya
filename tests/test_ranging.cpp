#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>

#include "lp/ranging.hpp"
#include "lp_generators.hpp"
#include "samaya/solver.hpp"
#include "test_framework.hpp"

using samaya::Index;
using samaya::Model;
using samaya::Result;
using samaya::Status;
using samaya::test::LpFamily;

namespace {

constexpr double kKnownTol = 1e-9;
// Match the requested relative finite-difference accuracy, without changing verifier settings.
constexpr double kObjectiveTol = 1e-6;
constexpr double kNondegenerateTol = 1e-7;
constexpr double kOutsideImprovementTol = 1e-7;
constexpr double kInsideFraction = 0.9;
constexpr double kOutsideFraction = 0.1;
constexpr double kMinimumOutsideStep = 1e-3;
constexpr double kInfiniteCostStep = 10.0;
constexpr int kMinimumOutsideColumns = 30;

samaya::Params ranging_params() {
  samaya::Params params;
  params.ranging = true;
  params.lp_method = samaya::LpMethod::kDualSimplex;
  params.log_level = 0;
  return params;
}

Model wyndor() {
  Model model;
  model.sense = samaya::ObjSense::kMaximize;
  model.obj = {3.0, 5.0};
  model.col_type.assign(2, samaya::VarType::kContinuous);
  model.col_lower = {0.0, 0.0};
  model.col_upper = {samaya::kInf, samaya::kInf};
  model.row_lower = {-samaya::kInf, -samaya::kInf, -samaya::kInf};
  model.row_upper = {4.0, 12.0, 18.0};
  model.row_names = {"x <= 4", "2y <= 12", "3x + 2y <= 18"};
  model.A = samaya::SparseMatrix::from_triplets(3, 2,
      {{0, 0, 1.0}, {1, 1, 2.0}, {2, 0, 3.0}, {2, 1, 2.0}});
  return model;
}

void check_wyndor(const Result& result, bool maximize) {
  REQUIRE(result.status == Status::kOptimal);
  REQUIRE(result.verified);
  REQUIRE(result.ranging.available);
  REQUIRE(result.ranging.cost_lower.size() == 2);
  REQUIRE(result.ranging.row_level_lower.size() == 3);
  CHECK_NEAR(result.col_value[0], 2.0, kKnownTol);
  CHECK_NEAR(result.col_value[1], 6.0, kKnownTol);
  CHECK_NEAR(result.objective, maximize ? 36.0 : -36.0, kKnownTol);
  const auto& range = result.ranging;
  if (maximize) {
    CHECK_NEAR(range.cost_lower[0], 0.0, kKnownTol);
    CHECK_NEAR(range.cost_upper[0], 7.5, kKnownTol);
    CHECK_NEAR(range.cost_lower[1], 2.0, kKnownTol);
    CHECK(range.cost_upper[1] == samaya::kInf);
  } else {
    CHECK_NEAR(range.cost_lower[0], -7.5, kKnownTol);
    CHECK_NEAR(range.cost_upper[0], 0.0, kKnownTol);
    CHECK(range.cost_lower[1] == -samaya::kInf);
    CHECK_NEAR(range.cost_upper[1], -2.0, kKnownTol);
  }
  CHECK(std::isnan(range.col_level_lower[0]));
  CHECK(std::isnan(range.col_level_upper[0]));
  CHECK(std::isnan(range.col_level_lower[1]));
  CHECK(std::isnan(range.col_level_upper[1]));
  CHECK(std::isnan(range.row_level_lower[0]));
  CHECK(std::isnan(range.row_level_upper[0]));
  CHECK_NEAR(range.row_level_lower[1], 6.0, kKnownTol);
  CHECK_NEAR(range.row_level_upper[1], 18.0, kKnownTol);
  CHECK_NEAR(range.row_level_lower[2], 12.0, kKnownTol);
  CHECK_NEAR(range.row_level_upper[2], 24.0, kKnownTol);
}

double old_objective(const Model& model, const Result& old) {
  double objective = model.obj_offset;
  for (Index j = 0; j < model.num_cols(); ++j) objective += model.obj[j] * old.col_value[j];
  return objective;
}

bool nondegenerate(const Model& model, const Result& result) {
  for (Index k = 0; k < model.num_cols() + model.num_rows(); ++k) {
    const bool col = k < model.num_cols();
    const Index index = col ? k : k - model.num_cols();
    const double lo = col ? model.col_lower[index] : model.row_lower[index];
    const double hi = col ? model.col_upper[index] : model.row_upper[index];
    const double value = col ? result.col_value[index] : result.row_activity[index];
    const double dual = col ? result.col_dual[index] : result.row_dual[index];
    const double range = col ? result.ranging.col_level_lower[index]
                             : result.ranging.row_level_lower[index];
    if (std::isnan(range)) {
      if (value - lo <= kNondegenerateTol || hi - value <= kNondegenerateTol) return false;
    } else if (std::fabs(dual) <= kNondegenerateTol) {
      return false;
    }
  }
  return true;
}

double outside_step(double lo, double hi, double current) {
  // An infinite interval has no width: use its finite distance from the current coefficient.
  double width = 0.0;
  if (std::isfinite(lo) && std::isfinite(hi)) width = hi - lo;
  else if (std::isfinite(lo)) width = current - lo;
  else if (std::isfinite(hi)) width = hi - current;
  return std::max(kMinimumOutsideStep, kOutsideFraction * width);
}

void move_bound(Model& model, bool column, Index index, bool at_lower, bool fixed, double level) {
  auto& lower = column ? model.col_lower : model.row_lower;
  auto& upper = column ? model.col_upper : model.row_upper;
  if (fixed || at_lower) lower[index] = level;
  if (fixed || !at_lower) upper[index] = level;
}

}  // namespace

TEST(ranging_wyndor_maximization) {
  check_wyndor(samaya::Solver(ranging_params()).solve(wyndor()), true);
}

TEST(ranging_wyndor_minimization) {
  Model model = wyndor();
  model.sense = samaya::ObjSense::kMinimize;
  model.obj = {-3.0, -5.0};
  check_wyndor(samaya::Solver(ranging_params()).solve(model), false);
}

TEST(ranging_presolve_and_barrier) {
  for (const bool presolve : {false, true}) {
    for (const auto method : {samaya::LpMethod::kDualSimplex, samaya::LpMethod::kBarrier}) {
      auto params = ranging_params();
      params.presolve = presolve;
      params.lp_method = method;
      check_wyndor(samaya::Solver(params).solve(wyndor()), true);
    }
  }
}

TEST(ranging_cost_finite_differences) {
  std::mt19937 rng(401);
  const samaya::Solver solver(ranging_params());
  int inside_checks = 0;
  int outside_columns = 0;
  int solved_models = 0;
  for (const auto family : {LpFamily::kFeasible, LpFamily::kBoxed}) {
    for (int sample = 0; sample < 120; ++sample) {
      const Model model = samaya::test::random_lp(family, 15, 15, rng);
      const Result old = solver.solve(model);
      if (old.status == Status::kUnbounded) {
        CHECK(old.verified);
        continue;
      }
      REQUIRE(old.status == Status::kOptimal);
      REQUIRE(old.verified);
      REQUIRE(old.ranging.available);
      ++solved_models;
      const bool strict = nondegenerate(model, old);
      for (Index j = 0; j < model.num_cols(); ++j) {
        const double current = model.obj[j];
        const double lo = old.ranging.cost_lower[j];
        const double hi = old.ranging.cost_upper[j];
        CHECK(lo <= current && current <= hi);
        bool checked_outside = false;
        for (const bool lower : {true, false}) {
          const double end = lower ? lo : hi;
          const double direction = lower ? -1.0 : 1.0;
          const double inside = std::isfinite(end) ? current + kInsideFraction * (end - current)
              : current + direction * kInfiniteCostStep * (1.0 + std::fabs(current));
          Model changed = model;
          changed.obj[j] = inside;
          const Result result = solver.solve(changed);
          REQUIRE(result.status == Status::kOptimal);
          CHECK(result.verified);
          const double prediction = old_objective(changed, old);
          CHECK_NEAR(result.objective, prediction, kObjectiveTol * (1.0 + std::fabs(prediction)));
          ++inside_checks;
          if (!strict || !std::isfinite(end)) continue;
          changed.obj[j] = end + direction * outside_step(lo, hi, current);
          const Result outside = solver.solve(changed);
          const double old_cost = old_objective(changed, old);
          if (outside.status == Status::kUnbounded) {
            CHECK(outside.verified);
            continue;
          }
          REQUIRE(outside.status == Status::kOptimal);
          CHECK(outside.verified);
          const double sense = model.sense == samaya::ObjSense::kMinimize ? 1.0 : -1.0;
          CHECK(sense * (old_cost - outside.objective) >
                kOutsideImprovementTol * (1.0 + std::fabs(outside.objective)));
          checked_outside = true;
        }
        if (checked_outside) ++outside_columns;
      }
    }
  }
  std::printf("  ranging costs: %d models, %d inside checks, %d outside columns\n",
              solved_models, inside_checks, outside_columns);
  CHECK(solved_models >= 150);
  CHECK(outside_columns >= kMinimumOutsideColumns);
}

TEST(ranging_bound_finite_differences) {
  std::mt19937 rng(402);
  const samaya::Solver solver(ranging_params());
  int inside_checks = 0;
  int outside_variables = 0;
  int solved_models = 0;
  for (const auto family : {LpFamily::kFeasible, LpFamily::kBoxed}) {
    for (int sample = 0; sample < 120; ++sample) {
      const Model model = samaya::test::random_lp(family, 15, 15, rng);
      const Result old = solver.solve(model);
      if (old.status == Status::kUnbounded) {
        CHECK(old.verified);
        continue;
      }
      REQUIRE(old.status == Status::kOptimal);
      REQUIRE(old.verified);
      REQUIRE(old.ranging.available);
      ++solved_models;
      const bool strict = nondegenerate(model, old);
      for (Index k = 0; k < model.num_cols() + model.num_rows(); ++k) {
        const bool col = k < model.num_cols();
        const Index index = col ? k : k - model.num_cols();
        const double lower = col ? model.col_lower[index] : model.row_lower[index];
        const double upper = col ? model.col_upper[index] : model.row_upper[index];
        const double current = col ? old.col_value[index] : old.row_activity[index];
        const double dual = col ? old.col_dual[index] : old.row_dual[index];
        const double lo = col ? old.ranging.col_level_lower[index] : old.ranging.row_level_lower[index];
        const double hi = col ? old.ranging.col_level_upper[index] : old.ranging.row_level_upper[index];
        if (std::isnan(lo)) {
          CHECK(std::isnan(hi));
          continue;
        }
        CHECK(lo <= current + kKnownTol && current <= hi + kKnownTol);
        const bool fixed = lower == upper;
        const bool at_lower = std::fabs(current - lower) <= kNondegenerateTol;
        // A free nonbasic has no binding bound to edit in the original model.
        if (!std::isfinite(lower) && !std::isfinite(upper)) continue;
        if (!fixed) {
          if (at_lower) CHECK(hi <= upper + kKnownTol);
          else CHECK(lo >= lower - kKnownTol);
        }
        bool checked_outside = false;
        for (const bool lower_end : {true, false}) {
          const double end = lower_end ? lo : hi;
          if (!std::isfinite(end)) continue;
          const double inside = current + kInsideFraction * (end - current);
          Model changed = model;
          move_bound(changed, col, index, at_lower, fixed, inside);
          const Result result = solver.solve(changed);
          REQUIRE(result.status == Status::kOptimal);
          CHECK(result.verified);
          const double prediction = old.objective + dual * (inside - current);
          if (std::fabs(result.objective - prediction) > kObjectiveTol * (1.0 + std::fabs(prediction))) {
            std::printf("  bound mismatch: family %s sample %d variable %d, level %.12g -> %.12g, "
                        "dual %.12g, range [%.12g, %.12g], basis note '%s'\n",
                        samaya::test::to_string(family), sample, k, current, inside, dual, lo, hi,
                        old.ranging.message.c_str());
          }
          CHECK_NEAR(result.objective, prediction, kObjectiveTol * (1.0 + std::fabs(prediction)));
          ++inside_checks;
          if (!strict) continue;
          const double direction = lower_end ? -1.0 : 1.0;
          const double outside_level = end + direction * outside_step(lo, hi, current);
          move_bound(changed, col, index, at_lower, fixed, outside_level);
          const Result outside = solver.solve(changed);
          // Leaving the feasible RHS interval can make the LP infeasible instead of changing price.
          if (outside.status == Status::kInfeasible) {
            CHECK(outside.verified);
            continue;
          }
          REQUIRE(outside.status == Status::kOptimal);
          CHECK(outside.verified);
          const double outside_prediction = old.objective + dual * (outside_level - current);
          CHECK(std::fabs(outside.objective - outside_prediction) >
                kOutsideImprovementTol * (1.0 + std::fabs(outside.objective)));
          checked_outside = true;
        }
        if (checked_outside) ++outside_variables;
      }
    }
  }
  std::printf("  ranging bounds: %d models, %d inside checks, %d outside variables\n",
              solved_models, inside_checks, outside_variables);
  CHECK(solved_models >= 150);
  CHECK(outside_variables >= kMinimumOutsideColumns);
}

TEST(ranging_opposite_bounds_fixed_and_free_columns) {
  Model model;
  model.obj = {2.0, -3.0, 7.0, 0.0};
  model.col_type.assign(4, samaya::VarType::kContinuous);
  model.col_lower = {1.0, -2.0, 4.0, -samaya::kInf};
  model.col_upper = {5.0, 6.0, 4.0, samaya::kInf};
  model.A = samaya::SparseMatrix::from_triplets(0, 4, {});
  const auto result = samaya::Solver(ranging_params()).solve(model);
  REQUIRE(result.status == Status::kOptimal);
  REQUIRE(result.ranging.available);
  const auto& range = result.ranging;
  CHECK_NEAR(range.cost_lower[0], 0.0, kKnownTol);
  CHECK(range.cost_upper[0] == samaya::kInf);
  CHECK(range.cost_lower[1] == -samaya::kInf);
  CHECK_NEAR(range.cost_upper[1], 0.0, kKnownTol);
  CHECK(range.cost_lower[2] == -samaya::kInf);
  CHECK(range.cost_upper[2] == samaya::kInf);
  CHECK_NEAR(range.cost_lower[3], 0.0, kKnownTol);
  CHECK_NEAR(range.cost_upper[3], 0.0, kKnownTol);
  CHECK(range.col_level_lower[0] == -samaya::kInf);
  CHECK_NEAR(range.col_level_upper[0], 5.0, kKnownTol);
  CHECK_NEAR(range.col_level_lower[1], -2.0, kKnownTol);
  CHECK(range.col_level_upper[1] == samaya::kInf);
  CHECK(range.col_level_lower[2] == -samaya::kInf);
  CHECK(range.col_level_upper[2] == samaya::kInf);

  // Both sides of a ranged row must also stop at its unchanged opposite bound.
  model.obj = {1.0};
  model.col_type.assign(1, samaya::VarType::kContinuous);
  model.col_lower = {-samaya::kInf};
  model.col_upper = {samaya::kInf};
  model.row_lower = {2.0};
  model.row_upper = {8.0};
  model.A = samaya::SparseMatrix::from_triplets(1, 1, {{0, 0, 1.0}});
  const auto row = samaya::Solver(ranging_params()).solve(model);
  REQUIRE(row.ranging.available);
  CHECK(row.ranging.row_level_lower[0] == -samaya::kInf);
  CHECK_NEAR(row.ranging.row_level_upper[0], 8.0, kKnownTol);
  model.obj[0] = -1.0;
  const auto upper_row = samaya::Solver(ranging_params()).solve(model);
  REQUIRE(upper_row.ranging.available);
  CHECK_NEAR(upper_row.ranging.row_level_lower[0], 2.0, kKnownTol);
  CHECK(upper_row.ranging.row_level_upper[0] == samaya::kInf);
}

TEST(ranging_not_available) {
  auto params = ranging_params();
  Model model = wyndor();
  model.col_type = {samaya::VarType::kInteger, samaya::VarType::kContinuous};
  const auto mip = samaya::Solver(params).solve(model);
  CHECK(mip.status == Status::kOptimal);
  CHECK(!mip.ranging.available);
  CHECK(mip.ranging.message == "ranging is available for LPs only");
  model = wyndor();
  model.row_lower[0] = 5.0;
  const auto infeasible = samaya::Solver(params).solve(model);
  CHECK(infeasible.status == Status::kInfeasible);
  CHECK(!infeasible.ranging.available);
  CHECK(infeasible.ranging.cost_lower.empty());
  params.ranging = false;
  const auto disabled = samaya::Solver(params).solve(wyndor());
  CHECK(disabled.status == Status::kOptimal);
  CHECK(!disabled.ranging.available);
  CHECK(disabled.ranging.message.empty());
  CHECK(disabled.ranging.cost_lower.empty());
  CHECK(disabled.ranging.col_level_lower.empty());
  CHECK(disabled.ranging.row_level_lower.empty());
}

TEST(ranging_input_and_time_guards) {
  const auto params = ranging_params();
  const Model model = wyndor();
  Result result = samaya::Solver(params).solve(model);
  REQUIRE(result.ranging.available);
  const samaya::Logger quiet(0);
  auto limited = params;
  limited.time_limit = 0.0;
  const auto timed = samaya::compute_ranging(model, result, limited, quiet);
  CHECK(!timed.available);
  CHECK(timed.message.find("time limit") != std::string::npos);
  result.col_value.clear();
  CHECK(!samaya::compute_ranging(model, result, params, quiet).available);
  result = samaya::Solver(params).solve(model);
  result.row_dual.clear();
  CHECK(!samaya::compute_ranging(model, result, params, quiet).available);
  result.status = Status::kInfeasible;
  CHECK(!samaya::compute_ranging(model, result, params, quiet).available);
}

TEST(ranging_reconstructed_basis_preserves_reported_shadow_prices) {
  // Duplicate tight rows admit different dual optima at the same primal vertex. A bound range
  // advertised with Result::row_dual must still predict the changed objective on its interior.
  Model model;
  model.obj = {-1.0};
  model.col_type = {samaya::VarType::kContinuous};
  model.col_lower = {0.0};
  model.col_upper = {samaya::kInf};
  model.row_lower = {-samaya::kInf, -samaya::kInf};
  model.row_upper = {1.0, 1.0};
  model.A = samaya::SparseMatrix::from_triplets(2, 1, {{0, 0, 1.0}, {1, 0, 1.0}});
  const samaya::Solver solver(ranging_params());
  const Result old = solver.solve(model);
  REQUIRE(old.status == Status::kOptimal);
  REQUIRE(old.ranging.available);
  int checked = 0;
  for (Index i = 0; i < 2; ++i) {
    if (std::isnan(old.ranging.row_level_lower[i])) continue;
    const double level = 0.5;
    REQUIRE(old.ranging.row_level_lower[i] < level && level < old.ranging.row_level_upper[i]);
    Model changed = model;
    changed.row_upper[i] = level;
    const Result result = solver.solve(changed);
    REQUIRE(result.status == Status::kOptimal);
    const double prediction = old.objective + old.row_dual[i] * (level - old.row_activity[i]);
    std::printf("  duplicate row %d: reported dual %.12g, range [%.12g, %.12g], "
                "changed objective %.12g, prediction %.12g\n", i, old.row_dual[i],
                old.ranging.row_level_lower[i], old.ranging.row_level_upper[i],
                result.objective, prediction);
    CHECK_NEAR(result.objective, prediction, kKnownTol);
    ++checked;
  }
  CHECK_EQ(checked, 1);
}
