// Crossover tests: a basis guessed from a primal-dual point must lead the warm-started simplex to
// the reference optimum, and from an interior-point answer in far fewer iterations than a cold
// start; through the solver, --crossover must report verified basic optima.

#include <cmath>
#include <cstdio>
#include <random>

#include "core/log.hpp"
#include "lp/crossover.hpp"
#include "lp/ipm.hpp"
#include "lp/lp_solver.hpp"
#include "lp_generators.hpp"
#include "reference_lp.hpp"
#include "samaya/solver.hpp"
#include "samaya/verify.hpp"
#include "test_framework.hpp"

using samaya::Index;
using samaya::LpResult;
using samaya::LpSolveOptions;
using samaya::Model;
using samaya::SimplexStatus;
using samaya::VarStatus;
using samaya::test::LpFamily;
using samaya::test::ReferenceLp;
using samaya::test::ReferenceResult;

namespace {

bool close(double a, double b) { return std::fabs(a - b) <= 1e-6 * (1.0 + std::fabs(b)); }

}  // namespace

TEST(crossover_basis_has_m_basic_variables_on_valid_sides) {
  std::mt19937 rng(301);
  const samaya::Logger quiet(0);
  int checked = 0;
  for (const LpFamily family : {LpFamily::kFeasible, LpFamily::kDegenerate, LpFamily::kBoxed,
                                LpFamily::kLoose}) {
    for (int k = 0; k < 60; ++k) {
      const Model model = samaya::test::random_lp(family, 12, 12, rng);
      const samaya::IpmResult ipm = samaya::solve_ipm(model, {}, quiet);
      if (ipm.x.size() != static_cast<std::size_t>(model.num_cols())) continue;
      const std::vector<VarStatus> status = samaya::crossover_basis(model, ipm.x, ipm.y);
      const Index n = model.num_cols();
      const Index m = model.num_rows();
      REQUIRE(status.size() == static_cast<std::size_t>(n + m));
      Index basic = 0;
      for (Index j = 0; j < n + m; ++j) {
        const double lower = j < n ? model.col_lower[j] : model.row_lower[j - n];
        const double upper = j < n ? model.col_upper[j] : model.row_upper[j - n];
        switch (status[j]) {
          case VarStatus::kBasic: ++basic; break;
          case VarStatus::kAtLower: CHECK(lower > -samaya::kInf); break;
          case VarStatus::kAtUpper: CHECK(upper < samaya::kInf); break;
          case VarStatus::kAtZero: CHECK(lower == -samaya::kInf && upper == samaya::kInf); break;
        }
      }
      CHECK_EQ(basic, m);
      ++checked;
    }
  }
  CHECK(checked > 150);
}

TEST(crossover_from_interior_point_reaches_reference_optimum_quickly) {
  // From the interior-point answer, the warm-started simplex must reach the reference optimum
  // with a verified basic solution, and need far fewer iterations than a cold start: that is
  // what makes the crossover a cleanup rather than a second solve.
  std::mt19937 rng(302);
  const samaya::Logger quiet(0);
  int solved = 0;
  int mismatches = 0;
  long long warm_iterations = 0;
  long long cold_iterations = 0;
  for (const LpFamily family : {LpFamily::kFeasible, LpFamily::kDegenerate, LpFamily::kBoxed}) {
    for (int k = 0; k < 80; ++k) {
      const Model model = samaya::test::random_lp(family, 30, 30, rng);
      const ReferenceResult ref = ReferenceLp(model).solve();
      if (ref.status != ReferenceResult::Status::kOptimal) continue;
      const samaya::IpmResult ipm = samaya::solve_ipm(model, {}, quiet);
      if (ipm.status != samaya::IpmStatus::kOptimal) continue;

      LpSolveOptions warm_options;
      warm_options.start = samaya::crossover_basis(model, ipm.x, ipm.y);
      const LpResult warm = samaya::solve_lp(model, warm_options, quiet);
      const LpResult cold = samaya::solve_lp(model, LpSolveOptions{}, quiet);
      REQUIRE(cold.status == SimplexStatus::kOptimal);
      const bool ok = warm.status == SimplexStatus::kOptimal &&
                      close(warm.objective, ref.objective) &&
                      samaya::verify_lp_optimality(model, warm.col_value, warm.row_dual).ok;
      if (!ok) {
        ++mismatches;
        std::printf("  mismatch (%s #%d): crossover %s %.10g vs reference %.10g\n",
                    samaya::test::to_string(family), k, samaya::to_string(warm.status),
                    warm.objective, ref.objective);
      }
      warm_iterations += warm.iterations;
      cold_iterations += cold.iterations;
      ++solved;
    }
  }
  std::printf("  %d models: %lld simplex iterations after crossover, %lld from the slack basis\n",
              solved, warm_iterations, cold_iterations);
  CHECK_EQ(mismatches, 0);
  CHECK(solved > 150);
  CHECK(warm_iterations * 4 < cold_iterations);
}

TEST(crossover_through_the_solver_reports_verified_basic_optima) {
  // --lp-method barrier --crossover and --lp-method pdlp --crossover: the reported answer is the
  // simplex's, from the crossover basis, verified, and equal to the reference optimum; and the
  // simplex needs far fewer iterations than the dual simplex alone (so the basis is really used).
  std::mt19937 rng(303);
  int mismatches = 0;
  int crossed = 0;
  long long crossover_iterations = 0;
  long long dual_iterations = 0;
  for (const samaya::LpMethod method : {samaya::LpMethod::kBarrier, samaya::LpMethod::kPdlp}) {
    for (int k = 0; k < 60; ++k) {
      const Model model = samaya::test::random_lp(LpFamily::kFeasible, 10, 10, rng);
      const ReferenceResult ref = ReferenceLp(model).solve();
      samaya::Params params;
      params.log_level = 0;
      params.lp_method = method;
      params.crossover = true;
      const samaya::Result r = samaya::Solver(params).solve(model);
      using RS = ReferenceResult::Status;
      const bool same =
          (ref.status == RS::kOptimal && r.status == samaya::Status::kOptimal && r.verified &&
           close(r.objective, ref.objective)) ||
          (ref.status == RS::kInfeasible && r.status == samaya::Status::kInfeasible) ||
          (ref.status == RS::kUnbounded && r.status == samaya::Status::kUnbounded);
      if (!same) {
        ++mismatches;
        std::printf("  mismatch (%s #%d): %s %.10g vs reference %.10g\n",
                    method == samaya::LpMethod::kBarrier ? "barrier" : "pdlp", k,
                    samaya::to_string(r.status), r.objective, ref.objective);
      }
      if (r.status == samaya::Status::kOptimal && r.crossover) {
        ++crossed;
        crossover_iterations += r.simplex_iterations;
        samaya::Params dual = params;
        dual.lp_method = samaya::LpMethod::kDualSimplex;
        dual_iterations += samaya::Solver(dual).solve(model).simplex_iterations;
      }
    }
  }
  std::printf("  %d crossovers: %lld simplex iterations, %lld by the dual simplex alone\n", crossed,
              crossover_iterations, dual_iterations);
  CHECK_EQ(mismatches, 0);
  CHECK(crossed > 80);
  CHECK(crossover_iterations * 4 < dual_iterations);
}

TEST(crossover_that_cannot_finish_falls_back_to_the_dual_simplex) {
  // With no crossover iterations allowed, the crossover stops without an answer; the dual
  // simplex must then solve the model from the slack basis, and the reported answer is its,
  // verified (a crossover stopped at its iteration limit used to be reported as the answer).
  std::mt19937 rng(304);
  int mismatches = 0;
  int fell_back = 0;
  for (const samaya::LpMethod method : {samaya::LpMethod::kBarrier, samaya::LpMethod::kPdlp}) {
    for (int k = 0; k < 40; ++k) {
      const Model model = samaya::test::random_lp(LpFamily::kFeasible, 10, 10, rng);
      const ReferenceResult ref = ReferenceLp(model).solve();
      if (ref.status != ReferenceResult::Status::kOptimal) continue;
      samaya::Params params;
      params.log_level = 0;
      params.lp_method = method;
      params.crossover = true;
      params.crossover_max_iterations = 0;
      const samaya::Result r = samaya::Solver(params).solve(model);
      if (!(r.status == samaya::Status::kOptimal && r.verified &&
            close(r.objective, ref.objective))) {
        ++mismatches;
        std::printf("  mismatch (#%d): %s %.10g vs reference %.10g\n", k,
                    samaya::to_string(r.status), r.objective, ref.objective);
      }
      if (!r.crossover) ++fell_back;
    }
  }
  CHECK_EQ(mismatches, 0);
  CHECK(fell_back > 20);
}
