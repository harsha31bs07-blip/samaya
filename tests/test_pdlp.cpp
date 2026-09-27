// PDLP tests: basic functionality and cross-check against reference simplex.

#include <cmath>
#include <cstdio>
#include <random>

#include "core/log.hpp"
#include "lp/lp_solver.hpp"
#include "lp/pdlp.hpp"
#include "lp_generators.hpp"
#include "reference_lp.hpp"
#include "samaya/solver.hpp"
#include "samaya/verify.hpp"
#include "test_framework.hpp"

using samaya::LpProblem;
using samaya::LpResult;
using samaya::LpSolveOptions;
using samaya::Model;
using samaya::PdlpOptions;
using samaya::PdlpResult;
using samaya::Scaling;
using samaya::SimplexStatus;
using samaya::test::LpFamily;
using samaya::test::ReferenceLp;
using samaya::test::ReferenceResult;

namespace {

bool same_status(SimplexStatus a, ReferenceResult::Status b) {
  switch (b) {
    case ReferenceResult::Status::kOptimal:
      return a == SimplexStatus::kOptimal;
    case ReferenceResult::Status::kInfeasible:
      return a == SimplexStatus::kInfeasible;
    case ReferenceResult::Status::kUnbounded:
      return a == SimplexStatus::kUnbounded;
  }
  return false;
}

const char* to_string(ReferenceResult::Status s) {
  switch (s) {
    case ReferenceResult::Status::kOptimal:
      return "optimal";
    case ReferenceResult::Status::kInfeasible:
      return "infeasible";
    case ReferenceResult::Status::kUnbounded:
      return "unbounded";
  }
  return "?";
}

void test_pdlp_basic() {
  // Tiny LP: min x s.t. x >= 1, x <= 2
  Model model;
  model.sense = samaya::ObjSense::kMinimize;
  model.obj = {1.0};
  model.col_lower = {1.0};
  model.col_upper = {2.0};
  model.col_type = {samaya::VarType::kContinuous};
  model.row_lower = {};
  model.row_upper = {};
  model.A = samaya::SparseMatrix(0, 1);

  const Scaling scaled = compute_scaling(model.A);
  const LpProblem scaled_lp = make_problem(model, scaled);

  PdlpOptions opts;
  opts.tol = 1e-4;  // PDLP design tolerance (per ARCHITECTURE.md M5)
  opts.max_iterations = 10000;
  samaya::Logger quiet(0);

  PdlpResult result = solve_pdlp(scaled_lp, opts, quiet);
  CHECK(same_status(result.status, ReferenceResult::Status::kOptimal));
  CHECK(result.iterations > 0);
  CHECK(result.primal_residual < opts.tol);
  CHECK(result.dual_residual < opts.tol);
  CHECK(result.gap < opts.tol);

  // Check solution: should be x = 1.0
  CHECK_NEAR(result.x[0], 1.0, 1e-4);
}

void test_pdlp_random_feasible() {
  std::mt19937 rng(42);
  const samaya::Logger log(0);
  PdlpOptions opts;
  opts.tol = 1e-8;               // The strict end of ARCHITECTURE.md M5 (1e-4 and 1e-8)
  opts.max_iterations = 100000;  // More iterations for difficult unbounded cases

  int optimal = 0, infeasible = 0, unbounded = 0, mismatches = 0;
  int verified_optimal = 0, verified_infeasible = 0, verified_unbounded = 0;

  for (int k = 0; k < 50; ++k) {
    Model model = samaya::test::random_lp(LpFamily::kFeasible, 8, 8, rng);
    if (model.Q.empty() == false) continue;  // Skip if Q is not empty

    const ReferenceResult ref = ReferenceLp(model).solve();
    const Scaling scaled = compute_scaling(model.A);
    const LpProblem scaled_lp = make_problem(model, scaled);
    PdlpResult pdlp = solve_pdlp(scaled_lp, opts, log);

    // Map PDLP result back to original space for verification
    LpResult lp;
    lp.status = pdlp.status;
    const int n = model.num_cols();
    lp.col_value.assign(pdlp.x.begin(), pdlp.x.begin() + n);
    scaled.unscale_cols(lp.col_value);

    // Compute row_activity = A * x in original space
    lp.row_activity.assign(model.num_rows(), 0.0);
    const auto a_start = model.A.col_start();
    const auto a_index = model.A.row_index();
    const auto a_value = model.A.values();
    for (int j = 0; j < n; ++j) {
      const double xj = lp.col_value[j];
      if (xj == 0.0) continue;
      for (samaya::NnzIndex p = a_start[j]; p < a_start[j + 1]; ++p) {
        lp.row_activity[a_index[p]] += a_value[p] * xj;
      }
    }

    // Compute objective = c^T x + obj_offset in original space
    lp.objective = model.obj_offset;
    for (int j = 0; j < n; ++j) {
      lp.objective += model.obj[j] * lp.col_value[j];
    }

    // Copy certificates if present
    if (!pdlp.dual_ray.empty()) {
      lp.dual_ray = pdlp.dual_ray;
      scaled.unscale_row_duals(lp.dual_ray);
    }
    if (!pdlp.primal_ray.empty()) {
      lp.primal_ray = pdlp.primal_ray;
      // primal_ray is in scaled space matching x, unscale cols
      scaled.unscale_cols(lp.primal_ray);
    }

    // The independent verifier at its default (strict) tolerances.
    const samaya::VerifyTolerances vtol;

    bool verified = false;
    switch (pdlp.status) {
      case SimplexStatus::kOptimal: {
        // Map dual variables (y) back to original space
        // PDLP uses Lagrangian L = c^T x + y^T (A x - s) so reduced cost = c + A^T y.
        // Verifier expects reduced cost = c - A^T y, so negate y.
        std::vector<double> y_orig = pdlp.y;
        scaled.unscale_row_duals(y_orig);
        for (double& y : y_orig) y = -y;  // Sign flip for verifier convention
        // Adjust for objective sense (same as src/core/solver.cpp:180-181)
        const double sense = model.sense == samaya::ObjSense::kMaximize ? -1.0 : 1.0;
        for (double& y : y_orig) y *= sense;
        auto report = samaya::verify_lp_optimality(model, lp.col_value, y_orig, vtol);
        if (report.ok) {
          verified = true;
          ++verified_optimal;
        } else {
          std::fprintf(stderr, "  Verification failed (optimal): %s\n", report.message.c_str());
        }
        break;
      }
      case SimplexStatus::kInfeasible: {
        auto report = samaya::verify_infeasibility(model, lp.dual_ray, vtol);
        if (report.ok) {
          verified = true;
          ++verified_infeasible;
        } else {
          std::fprintf(stderr, "  Verification failed (infeasible): %s\n", report.message.c_str());
        }
        break;
      }
      case SimplexStatus::kUnbounded: {
        auto report = samaya::verify_unbounded_ray(model, lp.primal_ray, vtol);
        if (report.ok) {
          verified = true;
          ++verified_unbounded;
        } else {
          std::fprintf(stderr, "  Verification failed (unbounded): %s\n", report.message.c_str());
        }
        break;
      }
      default:
        break;
    }

    // Status should match reference (allowing for PDLP detecting infeasible/unbounded that
    // reference might miss)
    bool ok = same_status(pdlp.status, ref.status) || verified;
    // Both optimal: the objectives must agree to PDLP's tolerance (a status match alone let a
    // wrong optimum through).
    if (ok && pdlp.status == SimplexStatus::kOptimal &&
        ref.status == ReferenceResult::Status::kOptimal &&
        std::fabs(lp.objective - ref.objective) > 1e-6 * (1.0 + std::fabs(ref.objective))) {
      ok = false;
    }

    if (!ok) {
      ++mismatches;
      std::fprintf(stderr, "  %s #%d (%dx%d): ref %s %.9g, pdlp %s %.9g (verified=%d)\n",
                   samaya::test::to_string(LpFamily::kFeasible), k, model.num_rows(),
                   model.num_cols(), to_string(ref.status), ref.objective,
                   samaya::to_string(pdlp.status), lp.objective, verified);
    }
    switch (ref.status) {
      case ReferenceResult::Status::kOptimal:
        ++optimal;
        break;
      case ReferenceResult::Status::kInfeasible:
        ++infeasible;
        break;
      case ReferenceResult::Status::kUnbounded:
        ++unbounded;
        break;
    }
  }
  std::printf("  PDLP feasible: %d optimal, %d infeasible, %d unbounded, %d mismatches\n", optimal,
              infeasible, unbounded, mismatches);
  std::printf("  Verified: %d optimal, %d infeasible, %d unbounded\n", verified_optimal,
              verified_infeasible, verified_unbounded);
  CHECK(optimal > 10);  // Should solve most feasible problems
  CHECK(mismatches == 0);
}

void test_pdlp_against_simplex() {
  // Compare PDLP with dual simplex on random problems - verify solutions independently
  std::mt19937 rng(123);
  const samaya::Logger quiet(0);
  PdlpOptions pdlp_opts;
  pdlp_opts.tol = 1e-4;               // PDLP design tolerance (per ARCHITECTURE.md M5)
  pdlp_opts.max_iterations = 100000;  // More iterations for difficult cases
  LpSolveOptions simplex_opts;
  simplex_opts.simplex.primal_tol = 1e-7;
  simplex_opts.simplex.dual_tol = 1e-7;

  int mismatches = 0;
  int compared = 0;
  int pdlp_verified = 0, simplex_verified = 0;

  for (int k = 0; k < 30; ++k) {
    Model model = samaya::test::random_lp(LpFamily::kFeasible, 10, 10, rng);
    if (model.Q.empty() == false) continue;

    const Scaling scaled = compute_scaling(model.A);
    const LpProblem scaled_lp = make_problem(model, scaled);

    PdlpResult pdlp = solve_pdlp(scaled_lp, pdlp_opts, quiet);
    LpResult simplex = samaya::solve_lp(model, simplex_opts, quiet);

    // Map PDLP result back to original space for verification
    LpResult lp;
    lp.status = pdlp.status;
    const int n = model.num_cols();
    lp.col_value.assign(pdlp.x.begin(), pdlp.x.begin() + n);
    scaled.unscale_cols(lp.col_value);

    // Compute row_activity = A * x in original space
    lp.row_activity.assign(model.num_rows(), 0.0);
    const auto a_start = model.A.col_start();
    const auto a_index = model.A.row_index();
    const auto a_value = model.A.values();
    for (int j = 0; j < n; ++j) {
      const double xj = lp.col_value[j];
      if (xj == 0.0) continue;
      for (samaya::NnzIndex p = a_start[j]; p < a_start[j + 1]; ++p) {
        lp.row_activity[a_index[p]] += a_value[p] * xj;
      }
    }

    // Compute objective
    lp.objective = model.obj_offset;
    for (int j = 0; j < n; ++j) {
      lp.objective += model.obj[j] * lp.col_value[j];
    }

    // Copy certificates
    if (!pdlp.dual_ray.empty()) {
      lp.dual_ray = pdlp.dual_ray;
      scaled.unscale_row_duals(lp.dual_ray);
    }
    if (!pdlp.primal_ray.empty()) {
      lp.primal_ray = pdlp.primal_ray;
      scaled.unscale_cols(lp.primal_ray);
    }

    // Verify PDLP solution with PDLP's design tolerance
    // Using 2e-4 to account for numerical differences between PDLP residuals and verifier checks
    samaya::VerifyTolerances vtol;
    vtol.primal = 2e-4;
    vtol.dual = 2e-4;
    vtol.integrality = 2e-4;

    bool pdlp_ok = false;
    switch (pdlp.status) {
      case SimplexStatus::kOptimal: {
        // Map dual variables (y) back to original space
        // PDLP uses Lagrangian L = c^T x + y^T (A x - s) so reduced cost = c + A^T y.
        // Verifier expects reduced cost = c - A^T y, so negate y.
        std::vector<double> y_orig = pdlp.y;
        scaled.unscale_row_duals(y_orig);
        for (double& y : y_orig) y = -y;  // Sign flip for verifier convention
        // Adjust for objective sense (same as src/core/solver.cpp:180-181)
        const double sense = model.sense == samaya::ObjSense::kMaximize ? -1.0 : 1.0;
        for (double& y : y_orig) y *= sense;
        auto report = samaya::verify_lp_optimality(model, lp.col_value, y_orig, vtol);
        if (report.ok) {
          pdlp_ok = true;
          ++pdlp_verified;
        }
        break;
      }
      case SimplexStatus::kInfeasible: {
        auto report = samaya::verify_infeasibility(model, lp.dual_ray, vtol);
        if (report.ok) {
          pdlp_ok = true;
        }
        break;
      }
      case SimplexStatus::kUnbounded: {
        auto report = samaya::verify_unbounded_ray(model, lp.primal_ray, vtol);
        if (report.ok) {
          pdlp_ok = true;
        }
        break;
      }
      case SimplexStatus::kIterationLimit:
      case SimplexStatus::kTimeLimit:
      case SimplexStatus::kNumericalError:
        break;  // No claim to check.
    }

    // Verify simplex solution with tighter tolerance (simplex is exact)
    samaya::VerifyTolerances simplex_vtol;
    simplex_vtol.primal = 1e-7;
    simplex_vtol.dual = 1e-7;
    simplex_vtol.integrality = 1e-7;

    bool simplex_ok = false;
    switch (simplex.status) {
      case SimplexStatus::kOptimal: {
        // Use full verifier for simplex too
        auto report =
            samaya::verify_lp_optimality(model, simplex.col_value, simplex.row_dual, simplex_vtol);
        if (report.ok) {
          simplex_ok = true;
          ++simplex_verified;
        }
        break;
      }
      case SimplexStatus::kInfeasible:
        // Simplex infeasibility certificates not exposed, skip
        simplex_ok = true;
        break;
      case SimplexStatus::kUnbounded:
        // Simplex unbounded rays not exposed, skip
        simplex_ok = true;
        break;
      case SimplexStatus::kIterationLimit:
      case SimplexStatus::kTimeLimit:
      case SimplexStatus::kNumericalError:
        break;  // No claim to check.
    }

    // Check: when both verified, status should agree
    // PDLP may not verify on all problems (first-order method limitation on small problems)
    bool ok = true;
    if (pdlp_ok && simplex_ok && pdlp.status != simplex.status) {
      ok = false;
    }
    // Also check that PDLP doesn't claim optimal with a non-verified solution when simplex says
    // unbounded/infeasible
    if (!pdlp_ok && simplex_ok &&
        (simplex.status == SimplexStatus::kInfeasible ||
         simplex.status == SimplexStatus::kUnbounded)) {
      // PDLP failed to verify but simplex found a certificate - this is a mismatch
      ok = false;
    }

    if (!ok) {
      ++mismatches;
      std::fprintf(stderr, "  mismatch #%d: pdlp=%s (verified=%d) simplex=%s (verified=%d)\n", k,
                   samaya::to_string(pdlp.status), pdlp_ok, samaya::to_string(simplex.status),
                   simplex_ok);
    }
    ++compared;
  }
  std::printf(
      "  PDLP vs Simplex: %d compared, %d mismatches (PDLP verified: %d, Simplex verified: %d)\n",
      compared, mismatches, pdlp_verified, simplex_verified);
  CHECK(mismatches == 0);
}

void test_pdlp_infeasible_unbounded() {
  std::mt19937 rng(4242);
  const samaya::Logger log(0);
  PdlpOptions opts;
  opts.tol = 1e-4;  // PDLP design tolerance (per ARCHITECTURE.md M5)
  opts.max_iterations = 50000;

  int infeasible_detected = 0, unbounded_detected = 0;
  int infeasible_verified = 0, unbounded_verified = 0;

  // Test infeasible problems (LpFamily::kRandom often generates infeasible)
  for (int k = 0; k < 30; ++k) {
    Model model = samaya::test::random_lp(LpFamily::kRandom, 8, 8, rng);
    if (model.Q.empty() == false) continue;

    const ReferenceResult ref = ReferenceLp(model).solve();
    if (ref.status != ReferenceResult::Status::kInfeasible) continue;

    const Scaling scaled = compute_scaling(model.A);
    const LpProblem scaled_lp = make_problem(model, scaled);
    PdlpResult pdlp = solve_pdlp(scaled_lp, opts, log);

    if (pdlp.status == SimplexStatus::kInfeasible && !pdlp.dual_ray.empty()) {
      ++infeasible_detected;
      // Unscale the Farkas certificate to original space for verification
      std::vector<double> dual_ray_orig = pdlp.dual_ray;
      scaled.unscale_row_duals(dual_ray_orig);

      LpResult lp;
      lp.status = pdlp.status;
      lp.dual_ray = dual_ray_orig;

      auto report = samaya::verify_infeasibility(model, lp.dual_ray);
      if (report.ok) {
        ++infeasible_verified;
      } else {
        std::fprintf(stderr, "  Infeasible verification failed: %s\n", report.message.c_str());
      }
    }
  }

  // Test unbounded problems (LpFamily::kLoose often generates unbounded)
  for (int k = 0; k < 30; ++k) {
    Model model = samaya::test::random_lp(LpFamily::kLoose, 8, 8, rng);
    if (model.Q.empty() == false) continue;

    const ReferenceResult ref = ReferenceLp(model).solve();
    if (ref.status != ReferenceResult::Status::kUnbounded) continue;

    const Scaling scaled = compute_scaling(model.A);
    const LpProblem scaled_lp = make_problem(model, scaled);
    PdlpResult pdlp = solve_pdlp(scaled_lp, opts, log);

    if (pdlp.status == SimplexStatus::kUnbounded && !pdlp.primal_ray.empty()) {
      ++unbounded_detected;
      // Ray is already in original space (matching result.x), no unscaling needed.
      LpResult lp;
      lp.status = pdlp.status;
      lp.primal_ray = pdlp.primal_ray;

      auto report = samaya::verify_unbounded_ray(model, lp.primal_ray);
      if (report.ok) {
        ++unbounded_verified;
      } else {
        std::fprintf(stderr, "  Unbounded verification failed: %s\n", report.message.c_str());
      }
    }
  }

  std::printf("  PDLP infeasible: %d detected, %d verified\n", infeasible_detected,
              infeasible_verified);
  std::printf("  PDLP unbounded: %d detected, %d verified\n", unbounded_detected,
              unbounded_verified);
  // At least some should be detected and verified
  CHECK(infeasible_detected > 0 || unbounded_detected > 0);
  CHECK(infeasible_verified > 0 || unbounded_verified > 0);
}

}  // namespace

TEST(pdlp_basic) {
  test_pdlp_basic();
}

TEST(pdlp_random_feasible) {
  test_pdlp_random_feasible();
}

TEST(pdlp_vs_simplex) {
  test_pdlp_against_simplex();
}

TEST(pdlp_infeasible_unbounded) {
  test_pdlp_infeasible_unbounded();
}
TEST(pdlp_through_the_solver_is_verified) {
  // --lp-method pdlp through the public API: every answer must match the reference simplex, and
  // most optimal answers must be PDLP's own (verified without the dual simplex fallback), which
  // checks the mapping of PDLP's solution, duals and rays back to the model.
  std::mt19937 rng(515);
  int models = 0;
  int by_pdlp = 0;
  int optimal = 0;
  int mismatches = 0;
  for (int k = 0; k < 60; ++k) {
    const Model model = samaya::test::random_lp(LpFamily::kFeasible, 8, 8, rng);
    const ReferenceResult ref = samaya::test::ReferenceLp(model).solve();
    samaya::Params params;
    params.log_level = 0;
    params.presolve = false;
    params.lp_method = samaya::LpMethod::kPdlp;
    const samaya::Result r = samaya::Solver(params).solve(model);
    ++models;
    const bool same =
        (ref.status == ReferenceResult::Status::kOptimal && r.status == samaya::Status::kOptimal &&
         r.verified &&
         std::fabs(r.objective - ref.objective) <= 1e-6 * (1.0 + std::fabs(ref.objective))) ||
        (ref.status == ReferenceResult::Status::kInfeasible &&
         r.status == samaya::Status::kInfeasible) ||
        (ref.status == ReferenceResult::Status::kUnbounded &&
         r.status == samaya::Status::kUnbounded);
    if (!same) {
      ++mismatches;
      std::printf(
          "  mismatch #%d: %s %.10g (verified %d, %lld simplex iterations) vs reference "
          "%d %.10g\n",
          k, samaya::to_string(r.status), r.objective, r.verified, r.simplex_iterations,
          static_cast<int>(ref.status), ref.objective);
    }
    if (ref.status == ReferenceResult::Status::kOptimal) {
      ++optimal;
      by_pdlp += r.status == samaya::Status::kOptimal && r.simplex_iterations == 0;
    }
  }
  std::printf("  %d models, %d optimal, %d answered by PDLP alone, %d mismatches\n", models,
              optimal, by_pdlp, mismatches);
  CHECK_EQ(mismatches, 0);
  CHECK(by_pdlp * 10 >= optimal * 8);  // At least 80% of the optima come from PDLP itself.
}
