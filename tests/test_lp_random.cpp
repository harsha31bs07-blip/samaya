// Differential tests: the sparse dual simplex against the dense reference simplex on randomly
// generated LPs.

#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <random>
#include <string>

#include "core/log.hpp"
#include "lp/lp_solver.hpp"
#include "lp_generators.hpp"
#include "reference_lp.hpp"
#include "samaya/verify.hpp"
#include "test_framework.hpp"

using samaya::LpResult;
using samaya::LpSolveOptions;
using samaya::Model;
using samaya::SimplexStatus;
using samaya::test::LpFamily;
using samaya::test::ReferenceLp;
using samaya::test::ReferenceResult;

namespace {

bool same_status(SimplexStatus a, ReferenceResult::Status b) {
  switch (b) {
    case ReferenceResult::Status::kOptimal: return a == SimplexStatus::kOptimal;
    case ReferenceResult::Status::kInfeasible: return a == SimplexStatus::kInfeasible;
    case ReferenceResult::Status::kUnbounded: return a == SimplexStatus::kUnbounded;
  }
  return false;
}

const char* to_string(ReferenceResult::Status s) {
  switch (s) {
    case ReferenceResult::Status::kOptimal: return "optimal";
    case ReferenceResult::Status::kInfeasible: return "infeasible";
    case ReferenceResult::Status::kUnbounded: return "unbounded";
  }
  return "?";
}

double max_primal_violation(const Model& m, const LpResult& r) {
  double worst = 0.0;
  for (samaya::Index j = 0; j < m.num_cols(); ++j) {
    worst = std::max({worst, m.col_lower[j] - r.col_value[j], r.col_value[j] - m.col_upper[j]});
  }
  for (samaya::Index i = 0; i < m.num_rows(); ++i) {
    worst = std::max(
        {worst, m.row_lower[i] - r.row_activity[i], r.row_activity[i] - m.row_upper[i]});
  }
  return worst;
}

struct Tally {
  int optimal = 0;
  int infeasible = 0;
  int unbounded = 0;
  int mismatches = 0;
  samaya::SimplexStats work;  // Summed over the optimal solves.
};

// Solves `count` random LPs of the family and compares with the reference.
Tally cross_check(LpFamily family, int count, int max_rows, int max_cols, unsigned seed,
                  const LpSolveOptions& options) {
  std::mt19937 rng(seed);
  const samaya::Logger quiet(0);
  Tally tally;
  for (int k = 0; k < count; ++k) {
    const Model model = samaya::test::random_lp(family, max_rows, max_cols, rng);
    const ReferenceResult ref = ReferenceLp(model).solve();
    const LpResult got = samaya::solve_lp(model, options, quiet);
    bool ok = same_status(got.status, ref.status);
    if (ok && ref.status == ReferenceResult::Status::kOptimal) {
      ok = std::fabs(got.objective - ref.objective) <= 1e-6 * (1.0 + std::fabs(ref.objective)) &&
           max_primal_violation(model, got) <= 1e-6;
    }
    // Every outcome must also pass the independent verifier.
    std::string verify_message;
    if (ok) {
      samaya::VerifyReport report;
      if (got.status == SimplexStatus::kOptimal) {
        report = samaya::verify_lp_optimality(model, got.col_value, got.row_dual);
      } else if (got.status == SimplexStatus::kInfeasible) {
        report = samaya::verify_infeasibility(model, got.dual_ray);
      } else {
        report = samaya::verify_unbounded_ray(model, got.primal_ray);
      }
      ok = report.ok;
      verify_message = report.message;
    }
    if (got.status == SimplexStatus::kOptimal) {
      tally.work.dual_iterations += got.stats.dual_iterations;
      tally.work.primal_iterations += got.stats.primal_iterations;
      tally.work.cost_shifts += got.stats.cost_shifts;
      tally.work.bound_flips += got.stats.bound_flips;
      tally.work.rebuilds_on_mismatch += got.stats.rebuilds_on_mismatch;
    }
    switch (ref.status) {
      case ReferenceResult::Status::kOptimal: ++tally.optimal; break;
      case ReferenceResult::Status::kInfeasible: ++tally.infeasible; break;
      case ReferenceResult::Status::kUnbounded: ++tally.unbounded; break;
    }
    if (!ok) {
      ++tally.mismatches;
      std::fprintf(stderr,
                   "  %s #%d (%dx%d): reference %s %.9g, solver %s %.9g (viol %.2e) %s\n",
                   samaya::test::to_string(family), k, model.num_rows(), model.num_cols(),
                   to_string(ref.status), ref.objective, samaya::to_string(got.status),
                   got.objective,
                   got.status == SimplexStatus::kOptimal ? max_primal_violation(model, got) : 0.0,
                   verify_message.c_str());
    }
  }
  std::printf("  %-10s %4d LPs: %d optimal, %d infeasible, %d unbounded, %d mismatches | "
              "optimal solves: %lld dual it, %lld primal it, %lld shifts, %lld mismatch rebuilds\n",
              samaya::test::to_string(family), count, tally.optimal, tally.infeasible,
              tally.unbounded, tally.mismatches, tally.work.dual_iterations,
              tally.work.primal_iterations, tally.work.cost_shifts,
              tally.work.rebuilds_on_mismatch);
  // The fallbacks (cost shifting, primal cleanup) can mask a broken dual simplex, so also check
  // that the dual simplex does the work.
  CHECK(tally.work.primal_iterations <= tally.work.dual_iterations / 10 + 5);
  CHECK(tally.work.cost_shifts <= tally.work.dual_iterations / 5 + 5);
  CHECK_EQ(tally.work.rebuilds_on_mismatch, 0);
  return tally;
}

}  // namespace

TEST(lp_random_default_options) {
  const LpSolveOptions options;
  int optimal = 0, infeasible = 0, unbounded = 0;
  for (LpFamily f : {LpFamily::kFeasible, LpFamily::kRandom, LpFamily::kDegenerate,
                     LpFamily::kLoose}) {
    const Tally t = cross_check(f, 250, 8, 8, 1000 + static_cast<unsigned>(f), options);
    CHECK_EQ(t.mismatches, 0);
    optimal += t.optimal;
    infeasible += t.infeasible;
    unbounded += t.unbounded;
  }
  // The families must actually exercise every outcome.
  CHECK(optimal > 200);
  CHECK(infeasible > 50);
  CHECK(unbounded > 50);
}

TEST(lp_random_boxed_models_use_bound_flips) {
  const LpSolveOptions options;
  const Tally small = cross_check(LpFamily::kBoxed, 300, 8, 8, 5000, options);
  CHECK_EQ(small.mismatches, 0);
  const Tally large = cross_check(LpFamily::kBoxed, 60, 30, 30, 5001, options);
  CHECK_EQ(large.mismatches, 0);
  CHECK(small.work.bound_flips + large.work.bound_flips > 100);
  CHECK(small.optimal + large.optimal > 100);
}

TEST(lp_random_larger_models) {
  const LpSolveOptions options;
  for (LpFamily f : {LpFamily::kFeasible, LpFamily::kDegenerate}) {
    CHECK_EQ(cross_check(f, 40, 30, 30, 2000 + static_cast<unsigned>(f), options).mismatches, 0);
  }
}

TEST(lp_random_without_scaling_or_perturbation) {
  LpSolveOptions options;
  options.scale = false;
  options.simplex.perturb = false;
  for (LpFamily f : {LpFamily::kFeasible, LpFamily::kDegenerate, LpFamily::kLoose}) {
    CHECK_EQ(cross_check(f, 150, 8, 8, 3000 + static_cast<unsigned>(f), options).mismatches, 0);
  }
}

TEST(lp_random_frequent_refactorization) {
  // A tiny refactorization interval exercises factorize/update switching on every model.
  LpSolveOptions options;
  options.simplex.refactor_interval = 2;
  for (LpFamily f : {LpFamily::kFeasible, LpFamily::kRandom, LpFamily::kDegenerate}) {
    CHECK_EQ(cross_check(f, 100, 20, 20, 4000 + static_cast<unsigned>(f), options).mismatches, 0);
  }
}

TEST(lp_dual_steepest_edge_weights_stay_accurate_on_long_runs) {
  // Solve each model once to learn its iteration count, then stop a second run three quarters of
  // the way and compare the updated weights with exact ones.
  std::mt19937 rng(91);
  const samaya::Logger quiet(0);
  int compared = 0;
  double worst = 0.0;
  for (int k = 0; k < 16; ++k) {
    const Model model = samaya::test::random_lp(LpFamily::kBoxed, 150, 200, rng);
    samaya::LpProblem lp;
    lp.m = model.num_rows();
    lp.n = model.num_cols();
    lp.A = model.A;
    lp.At = model.A.transpose();
    lp.cost = model.obj;
    lp.cost.resize(static_cast<std::size_t>(lp.n + lp.m), 0.0);
    lp.lower = model.col_lower;
    lp.upper = model.col_upper;
    lp.lower.insert(lp.lower.end(), model.row_lower.begin(), model.row_lower.end());
    lp.upper.insert(lp.upper.end(), model.row_upper.begin(), model.row_upper.end());
    samaya::SimplexOptions options;
    samaya::Simplex full(lp, options, quiet);
    full.solve();
    if (full.iterations() < 40 || full.stats().primal_iterations > 0) continue;
    options.max_iterations = full.iterations() * 3 / 4;
    samaya::Simplex simplex(lp, options, quiet);
    if (simplex.solve() != SimplexStatus::kIterationLimit) continue;
    const std::vector<double> exact = simplex.exact_dse_weights();
    for (std::size_t r = 0; r < exact.size(); ++r) {
      if (exact[r] < 1e-3) continue;
      ++compared;
      worst = std::max(worst, std::fabs(simplex.dse_weights()[r] - exact[r]) / exact[r]);
    }
  }
  std::printf("  compared %d weights, worst relative error %.2e\n", compared, worst);
  CHECK(compared > 300);
  CHECK(worst < 1e-6);
}

TEST(lp_dual_steepest_edge_weights_stay_exact) {
  // Stop the dual simplex after a few iterations and compare the updated weights with weights
  // recomputed from the factorization.
  std::mt19937 rng(77);
  const samaya::Logger quiet(0);
  int compared = 0;
  int wrong = 0;
  for (int k = 0; k < 60; ++k) {
    const Model model = samaya::test::random_lp(LpFamily::kFeasible, 25, 25, rng);
    for (long long limit : {3LL, 8LL, 20LL}) {
      samaya::LpProblem lp;
      lp.m = model.num_rows();
      lp.n = model.num_cols();
      lp.A = model.A;
      lp.At = model.A.transpose();
      lp.cost = model.obj;
      lp.cost.resize(static_cast<std::size_t>(lp.n + lp.m), 0.0);
      lp.lower = model.col_lower;
      lp.upper = model.col_upper;
      lp.lower.insert(lp.lower.end(), model.row_lower.begin(), model.row_lower.end());
      lp.upper.insert(lp.upper.end(), model.row_upper.begin(), model.row_upper.end());
      samaya::SimplexOptions options;
      options.max_iterations = limit;
      samaya::Simplex simplex(lp, options, quiet);
      if (simplex.solve() != SimplexStatus::kIterationLimit) continue;
      if (simplex.stats().primal_iterations > 0) continue;  // Primal steps reset the weights.
      const std::vector<double> exact = simplex.exact_dse_weights();
      const std::vector<double>& kept = simplex.dse_weights();
      for (std::size_t r = 0; r < exact.size(); ++r) {
        if (exact[r] < 1e-3) continue;  // Updates clamp weights at 1e-4.
        ++compared;
        if (std::fabs(kept[r] - exact[r]) > 1e-6 * exact[r]) ++wrong;
      }
    }
  }
  CHECK(compared > 500);
  CHECK_EQ(wrong, 0);
}

TEST(lp_warm_start_from_optimal_and_nearby_bases) {
  // From its own optimal basis the simplex needs no iterations; from the optimal basis of a
  // perturbed objective it must reach the same optimum as a cold start.
  std::mt19937 rng(123);
  const samaya::Logger quiet(0);
  int optimal_starts = 0;
  int nearby_starts = 0;
  for (int k = 0; k < 150; ++k) {
    const Model model = samaya::test::random_lp(LpFamily::kFeasible, 30, 40, rng);
    samaya::LpProblem lp;
    lp.m = model.num_rows();
    lp.n = model.num_cols();
    lp.A = model.A;
    lp.At = model.A.transpose();
    lp.cost = model.obj;
    lp.cost.resize(static_cast<std::size_t>(lp.n + lp.m), 0.0);
    lp.lower = model.col_lower;
    lp.upper = model.col_upper;
    lp.lower.insert(lp.lower.end(), model.row_lower.begin(), model.row_lower.end());
    lp.upper.insert(lp.upper.end(), model.row_upper.begin(), model.row_upper.end());
    const auto objective = [&](const samaya::Simplex& s) {
      double sum = 0.0;
      for (samaya::Index j = 0; j < lp.n; ++j) sum += lp.cost[j] * s.values()[j];
      return sum;
    };
    const samaya::SimplexOptions options;
    samaya::Simplex cold(lp, options, quiet);
    if (cold.solve() != SimplexStatus::kOptimal) continue;

    samaya::Simplex again(lp, options, quiet);
    REQUIRE(again.solve(cold.status()) == SimplexStatus::kOptimal);
    CHECK_EQ(again.iterations(), 0);
    CHECK_NEAR(objective(again), objective(cold), 1e-6 * (1.0 + std::fabs(objective(cold))));
    ++optimal_starts;

    samaya::LpProblem other = lp;
    for (samaya::Index j = 0; j < lp.n; ++j) {
      other.cost[j] += std::uniform_real_distribution<double>(-2.0, 2.0)(rng);
    }
    samaya::Simplex nearby(other, options, quiet);
    if (nearby.solve() != SimplexStatus::kOptimal) continue;
    samaya::Simplex warm(lp, options, quiet);
    REQUIRE(warm.solve(nearby.status()) == SimplexStatus::kOptimal);
    CHECK_NEAR(objective(warm), objective(cold), 1e-6 * (1.0 + std::fabs(objective(cold))));
    ++nearby_starts;
  }
  CHECK(optimal_starts > 60);
  CHECK(nearby_starts > 40);
}

TEST(lp_random_large_models_verify_first_time) {
  // Larger sparse models where steepest-edge weight drift once produced false optima. The dense
  // reference is too slow at this size, so the verifier is the oracle: every outcome must be
  // proven on the first attempt.
  std::mt19937 rng(606);
  const samaya::Logger quiet(0);
  int proven = 0;
  for (int k = 0; k < 12; ++k) {
    const LpFamily family = k % 2 == 0 ? LpFamily::kBoxed : LpFamily::kFeasible;
    const Model model = samaya::test::random_lp(family, 150, 220, rng);
    const LpResult got = samaya::solve_lp(model, LpSolveOptions{}, quiet);
    samaya::VerifyReport report;
    if (got.status == SimplexStatus::kOptimal) {
      report = samaya::verify_lp_optimality(model, got.col_value, got.row_dual);
    } else if (got.status == SimplexStatus::kInfeasible) {
      report = samaya::verify_infeasibility(model, got.dual_ray);
    } else if (got.status == SimplexStatus::kUnbounded) {
      report = samaya::verify_unbounded_ray(model, got.primal_ray);
    } else {
      report.ok = false;
      report.message = samaya::to_string(got.status);
    }
    if (!report.ok) {
      std::fprintf(stderr, "  large #%d (%dx%d) %s: %s\n", k, model.num_rows(), model.num_cols(),
                   samaya::to_string(got.status), report.message.c_str());
    }
    proven += report.ok ? 1 : 0;
  }
  CHECK_EQ(proven, 12);
}

TEST(lp_random_profile_preserves_objective_bits_and_iterations) {
  std::mt19937 rng(20261014);
  const samaya::Logger quiet(0);
  long long profiled_iterations = 0;
  const auto compare = [&](const LpResult& plain, const LpResult& profiled) {
    CHECK(plain.status == profiled.status);
    CHECK_EQ(plain.iterations, profiled.iterations);
    CHECK_EQ(std::bit_cast<std::uint64_t>(plain.objective),
             std::bit_cast<std::uint64_t>(profiled.objective));
    CHECK(plain.basis == profiled.basis);
    const auto same_bits = [](const std::vector<double>& a, const std::vector<double>& b) {
      REQUIRE(a.size() == b.size());
      for (std::size_t j = 0; j < a.size(); ++j) {
        CHECK_EQ(std::bit_cast<std::uint64_t>(a[j]), std::bit_cast<std::uint64_t>(b[j]));
      }
    };
    same_bits(plain.col_value, profiled.col_value);
    same_bits(plain.row_dual, profiled.row_dual);
    same_bits(plain.col_dual, profiled.col_dual);
    same_bits(plain.dual_ray, profiled.dual_ray);
    same_bits(plain.primal_ray, profiled.primal_ray);
    profiled_iterations += profiled.iterations;
    CHECK(!plain.stats.collected);
    CHECK(profiled.stats.collected);
    CHECK(profiled.stats.total_seconds > 0.0);
    CHECK_EQ(profiled.stats.factor.rebuilds, profiled.stats.refactorizations);
    CHECK_EQ(profiled.stats.dual_phase1_iterations + profiled.stats.dual_phase2_iterations +
                 profiled.stats.primal_cleanup_iterations, profiled.iterations);
    CHECK_EQ(profiled.stats.ratio_candidates.count, profiled.stats.harris_groups.count);
    double seconds = 0.0;
    for (double component : profiled.stats.seconds) {
      CHECK(component >= 0.0);
      seconds += component;
    }
    CHECK_NEAR(seconds, profiled.stats.total_seconds, 1e-10);
  };
  for (LpFamily family : {LpFamily::kFeasible, LpFamily::kRandom, LpFamily::kDegenerate,
                          LpFamily::kLoose, LpFamily::kBoxed}) {
    for (int k = 0; k < 16; ++k) {
      const Model model = samaya::test::random_lp(family, 30, 40, rng);
      LpSolveOptions options;
      options.scale = false;  // Inspect a single Simplex::solve profile, without unscaled polishing.
      options.simplex.refactor_interval = 2;
      options.simplex.perturb = k % 2 == 0;
      options.simplex.timing_sample_interval = k % 2 == 0 ? 1 : 32;
      const LpResult plain = samaya::solve_lp(model, options, quiet);
      options.simplex.collect_stats = true;
      const LpResult profiled = samaya::solve_lp(model, options, quiet);
      compare(plain, profiled);
      CHECK_EQ(profiled.stats.timing_samples,
               (profiled.stats.timing_attempts + profiled.stats.timing_sample_interval - 1) /
                   profiled.stats.timing_sample_interval);
      if (plain.status != SimplexStatus::kOptimal) continue;
      options.start = plain.basis;
      options.simplex.collect_stats = false;
      const LpResult warm_plain = samaya::solve_lp(model, options, quiet);
      options.simplex.collect_stats = true;
      const LpResult warm_profiled = samaya::solve_lp(model, options, quiet);
      compare(warm_plain, warm_profiled);
      CHECK_EQ(warm_profiled.stats.warm_start_iterations, warm_profiled.iterations);
    }
  }
  CHECK(profiled_iterations > 100);
}

TEST(lp_profile_counts_rank_repairs_and_resets_on_reused_warm_start) {
  samaya::LpProblem lp;
  lp.m = lp.n = 2;
  lp.A = samaya::SparseMatrix::from_triplets(2, 2,
      {{0, 0, 1.0}, {1, 0, 1.0}, {0, 1, 1.0}, {1, 1, 1.0}});
  lp.At = lp.A.transpose();
  lp.cost = {1.0, 2.0, 0.0, 0.0};
  lp.lower = {0.0, 0.0, 1.0, 1.0};
  lp.upper = {10.0, 10.0, samaya::kInf, samaya::kInf};
  samaya::SimplexOptions options;
  options.collect_stats = true;
  const samaya::Logger quiet(0);
  samaya::Simplex simplex(lp, options, quiet);
  const std::vector<samaya::VarStatus> singular = {
      samaya::VarStatus::kBasic, samaya::VarStatus::kBasic,
      samaya::VarStatus::kAtLower, samaya::VarStatus::kAtLower};
  REQUIRE(simplex.solve(singular) == SimplexStatus::kOptimal);
  CHECK_EQ(simplex.stats().warm_start_rank_repairs, 1);
  CHECK_EQ(simplex.stats().rank_repairs, 1);
  CHECK_EQ(simplex.stats().factor.reasons[static_cast<std::size_t>(
      samaya::BasisFactor::RebuildReason::kRankRepair)], 1);
  const auto start = simplex.status();
  REQUIRE(simplex.solve(start) == SimplexStatus::kOptimal);
  CHECK_EQ(simplex.iterations(), 0);
  CHECK_EQ(simplex.stats().warm_start_rank_repairs, 0);
  CHECK_EQ(simplex.stats().factor.rebuilds, 0);
  CHECK_EQ(simplex.stats().rho_nnz.count, 0);
}

TEST(lp_profile_counts_warm_start_flips_and_cost_shifts) {
  samaya::LpProblem lp;
  lp.m = 1;
  lp.n = 2;
  lp.A = samaya::SparseMatrix::from_triplets(1, 2, {{0, 0, 1.0}, {0, 1, 1.0}});
  lp.At = lp.A.transpose();
  lp.cost = {-1.0, -2.0, 0.0};
  lp.lower = {0.0, 0.0, -samaya::kInf};
  lp.upper = {2.0, samaya::kInf, 1.0};
  samaya::SimplexOptions options;
  options.collect_stats = true;
  const samaya::Logger quiet(0);
  samaya::Simplex simplex(lp, options, quiet);
  REQUIRE(simplex.solve({samaya::VarStatus::kAtLower, samaya::VarStatus::kAtLower,
                         samaya::VarStatus::kBasic}) == SimplexStatus::kOptimal);
  CHECK_EQ(simplex.stats().warm_start_bound_flips, 1);
  CHECK_EQ(simplex.stats().warm_start_cost_shifts, 1);
  CHECK_EQ(simplex.stats().warm_start_iterations, simplex.iterations());
}

TEST(lp_profile_ftran_density_counts_outputs_after_updates) {
  const auto matrix = samaya::SparseMatrix::from_triplets(2, 1,
      {{0, 0, 1.0}, {1, 0, 1.0}});
  samaya::BasisFactor factor(matrix);
  const std::vector<samaya::Index> logicals = {1, 2};
  REQUIRE(factor.factorize(logicals) == 0);
  samaya::MeasurementCount count;
  std::vector<double> column = {1.0, 1.0};
  std::vector<double> spike;
  factor.ftran(column, &spike, &count);
  CHECK_EQ(count.sum, 2);
  REQUIRE(factor.update(0, spike));
  // In the updated basis [a, -e_1], this RHS has solution (2, 0). The moved pivot
  // leaves a tombstone in the factor order; it must not be counted a second time.
  column = {2.0, 2.0};
  factor.ftran(column, nullptr, &count);
  CHECK_EQ(column[0], 2.0);
  CHECK_EQ(column[1], 0.0);
  column = {0.0, 0.0};
  factor.ftran(column, nullptr, &count);
  CHECK_EQ(count.count, 3);
  CHECK_EQ(count.sum, 3);
  CHECK_EQ(count.maximum, 2);
}
