// Convex QP and the barrier (interior-point) LP method against dense references.

#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

#include "lp_generators.hpp"
#include "reference_lp.hpp"
#include "reference_qp.hpp"
#include "samaya/solver.hpp"
#include "samaya/verify.hpp"
#include "test_framework.hpp"

using samaya::Model;
using samaya::ObjSense;
using samaya::Params;
using samaya::Result;
using samaya::Solver;
using samaya::SparseMatrix;
using samaya::Status;
using samaya::Triplet;
using samaya::test::ReferenceQpResult;

namespace {

enum class QpFamily {
  kStrictlyConvex,  // Q = R'R + 0.1 I.
  kLowRank,         // Q = R'R with fewer rows than columns (positive semidefinite, singular).
  kNearlyLinear,    // One small diagonal entry: an LP with a little curvature.
  kInfeasible,      // A row that no point of the box satisfies.
  kNonConvex,       // A negative diagonal entry.
};

Model random_qp(QpFamily family, std::mt19937& rng) {
  const auto uniform_int = [&](int lo, int hi) {
    return std::uniform_int_distribution<int>(lo, hi)(rng);
  };
  const auto real = [&](double lo, double hi) {
    return std::uniform_real_distribution<double>(lo, hi)(rng);
  };
  const int n = uniform_int(1, 4);
  const int m = uniform_int(1, 3);
  Model model;
  model.sense = uniform_int(0, 9) < 3 ? ObjSense::kMaximize : ObjSense::kMinimize;
  const double sense = model.sense == ObjSense::kMaximize ? -1.0 : 1.0;
  model.obj_offset = uniform_int(-3, 3);
  std::vector<double> x0(static_cast<std::size_t>(n));
  for (int j = 0; j < n; ++j) {
    model.obj.push_back(uniform_int(-5, 5));
    const double lo = uniform_int(-5, 0);
    model.col_lower.push_back(lo);
    model.col_upper.push_back(lo + uniform_int(1, 6));
    model.col_type.push_back(samaya::VarType::kContinuous);
    x0[j] = real(model.col_lower[j], model.col_upper[j]);
  }
  std::vector<Triplet> t;
  int equalities = 0;
  for (int i = 0; i < m; ++i) {
    double activity = 0.0;
    double max_act = 0.0;
    for (int j = 0; j < n; ++j) {
      const int a = uniform_int(-3, 3);
      if (a == 0) continue;
      t.push_back({i, j, static_cast<double>(a)});
      activity += a * x0[j];
      max_act += a > 0 ? a * model.col_upper[j] : a * model.col_lower[j];
    }
    // At most n - 1 equality rows, so the equalities never over-determine the columns.
    int kind = uniform_int(0, 3);
    if (kind == 0 && equalities + 1 >= n) kind = 3;
    equalities += kind == 0;
    double rl = -samaya::kInf;
    double ru = samaya::kInf;
    if (kind == 0) {
      rl = ru = activity;
    } else if (kind == 1) {
      ru = activity + real(0.0, 2.0);
    } else if (kind == 2) {
      rl = activity - real(0.0, 2.0);
    } else {
      rl = activity - real(0.0, 2.0);
      ru = activity + real(0.0, 2.0);
    }
    if (family == QpFamily::kInfeasible && i == 0) {
      rl = max_act + 1.0;  // Beyond the largest activity over the box.
      ru = samaya::kInf;
    }
    model.row_lower.push_back(rl);
    model.row_upper.push_back(ru);
  }
  model.A = SparseMatrix::from_triplets(m, n, std::move(t));

  // Q in the minimization sense, then stored with the model's sense (lower triangle).
  std::vector<std::vector<double>> q(static_cast<std::size_t>(n),
                                     std::vector<double>(static_cast<std::size_t>(n), 0.0));
  if (family == QpFamily::kNearlyLinear) {
    const int j = uniform_int(0, n - 1);
    q[j][j] = real(0.01, 0.1);
  } else {
    const int k = family == QpFamily::kLowRank ? std::max(1, n - 1) : n;
    for (int r = 0; r < k; ++r) {
      std::vector<double> row(static_cast<std::size_t>(n));
      for (double& v : row) v = uniform_int(-2, 2);
      for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) q[i][j] += row[i] * row[j];
      }
    }
    if (family != QpFamily::kLowRank) {
      for (int i = 0; i < n; ++i) q[i][i] += 0.1;
    }
    if (family == QpFamily::kNonConvex) q[0][0] = -1.0 - q[0][0];
  }
  std::vector<Triplet> qt;
  for (int j = 0; j < n; ++j) {
    for (int i = j; i < n; ++i) {
      if (q[i][j] != 0.0) qt.push_back({i, j, sense * q[i][j]});
    }
  }
  model.Q = SparseMatrix::from_triplets(n, n, std::move(qt));
  if (model.Q.empty()) {  // Keep it a QP: a tiny curvature on the first column.
    model.Q = SparseMatrix::from_triplets(n, n, {{0, 0, sense * 0.05}});
  }
  return model;
}

Result solve(const Model& model, samaya::LpMethod method = samaya::LpMethod::kAuto) {
  Params params;
  params.log_level = 0;
  params.lp_method = method;
  params.presolve = false;
  return Solver(params).solve(model);
}

}  // namespace

TEST(qp_random_models_match_reference) {
  // Convex QPs of 1-4 columns and 1-3 rows (equalities, one-sided and ranged rows, boxed
  // columns, both senses): the interior-point answer must pass verify_qp_optimality and match
  // the dense active-set reference's objective.
  std::mt19937 rng(71);
  int models = 0;
  int mismatches = 0;
  int unverified = 0;
  int by_ipm = 0;
  for (const QpFamily family :
       {QpFamily::kStrictlyConvex, QpFamily::kLowRank, QpFamily::kNearlyLinear}) {
    for (int k = 0; k < 150; ++k) {
      const Model model = random_qp(family, rng);
      const ReferenceQpResult ref = samaya::test::reference_qp(model);
      CHECK(ref.status == ReferenceQpResult::Status::kOptimal);  // The rows hold at x0.
      const Result r = solve(model);
      ++models;
      if (r.status != Status::kOptimal || !r.verified) {
        ++unverified;
        continue;
      }
      if (std::fabs(r.objective - ref.objective) > 1e-6 * (1.0 + std::fabs(ref.objective))) {
        ++mismatches;
        std::printf("  mismatch: %.10g vs reference %.10g\n", r.objective, ref.objective);
      }
      by_ipm += r.barrier_iterations > 0;
    }
  }
  std::printf("  %d models, %d unverified, %d objective mismatches\n", models, unverified,
              mismatches);
  CHECK(models == 450);
  CHECK_EQ(unverified, 0);
  CHECK_EQ(mismatches, 0);
  CHECK(by_ipm == models);
}

TEST(qp_infeasible_models_get_verified_certificates) {
  // A row beyond the box: the interior point cannot converge, and the constraints alone,
  // solved by the simplex, give a Farkas certificate that the verifier accepts.
  std::mt19937 rng(72);
  int infeasible = 0;
  for (int k = 0; k < 60; ++k) {
    const Model model = random_qp(QpFamily::kInfeasible, rng);
    const Result r = solve(model);
    CHECK(r.status == Status::kInfeasible);
    CHECK(r.verified);
    CHECK(samaya::verify_infeasibility(model, r.infeasibility_certificate).ok);
    infeasible += r.status == Status::kInfeasible;
  }
  CHECK_EQ(infeasible, 60);
}

TEST(qp_nonconvex_models_are_refused) {
  std::mt19937 rng(73);
  for (int k = 0; k < 60; ++k) {
    const Model model = random_qp(QpFamily::kNonConvex, rng);
    CHECK(solve(model).status == Status::kNotConvex);
  }
}

TEST(barrier_lp_matches_reference) {
  // --lp-method barrier on the random LP families: the answer must match the reference simplex
  // whether the interior point's own point verified or the dual simplex finished the solve; and
  // the interior point must verify on its own for most feasible models.
  std::mt19937 rng(74);
  int models = 0;
  int by_barrier = 0;
  int mismatches = 0;
  for (const samaya::test::LpFamily family :
       {samaya::test::LpFamily::kFeasible, samaya::test::LpFamily::kDegenerate,
        samaya::test::LpFamily::kBoxed, samaya::test::LpFamily::kRandom}) {
    for (int k = 0; k < 150; ++k) {
      const Model model = samaya::test::random_lp(family, 8, 8, rng);
      const samaya::test::ReferenceResult ref = samaya::test::ReferenceLp(model).solve();
      const Result r = solve(model, samaya::LpMethod::kBarrier);
      ++models;
      using RS = samaya::test::ReferenceResult::Status;
      const bool same =
          (ref.status == RS::kOptimal && r.status == Status::kOptimal && r.verified &&
           std::fabs(r.objective - ref.objective) <= 1e-6 * (1.0 + std::fabs(ref.objective))) ||
          (ref.status == RS::kInfeasible && r.status == Status::kInfeasible) ||
          (ref.status == RS::kUnbounded && r.status == Status::kUnbounded);
      if (!same) {
        ++mismatches;
        std::printf("  mismatch (%s): status %s objective %.10g vs reference %d %.10g\n",
                    samaya::test::to_string(family), samaya::to_string(r.status), r.objective,
                    static_cast<int>(ref.status), ref.objective);
      }
      if (ref.status == RS::kOptimal && r.status == Status::kOptimal &&
          r.simplex_iterations == 0) {
        ++by_barrier;
      }
    }
  }
  std::printf("  %d models, %d optimal by the interior point alone, %d mismatches\n", models,
              by_barrier, mismatches);
  CHECK_EQ(mismatches, 0);
  CHECK(by_barrier > 200);
}

namespace {

// Brute force for small MIQPs: every assignment of the integer columns within their (finite)
// bounds, the continuous rest solved by the dense QP reference; the best feasible one.
ReferenceQpResult brute_force_miqp(const Model& model) {
  std::vector<int> ints;
  for (int j = 0; j < model.num_cols(); ++j) {
    if (model.col_type[j] == samaya::VarType::kInteger) ints.push_back(j);
  }
  ReferenceQpResult best;
  const bool maximize = model.sense == ObjSense::kMaximize;
  std::vector<double> value(ints.size());
  for (std::size_t k = 0; k < ints.size(); ++k) value[k] = model.col_lower[ints[k]];
  for (;;) {
    Model fixed = model;
    fixed.col_type.assign(fixed.col_type.size(), samaya::VarType::kContinuous);
    for (std::size_t k = 0; k < ints.size(); ++k) {
      fixed.col_lower[ints[k]] = fixed.col_upper[ints[k]] = value[k];
    }
    const ReferenceQpResult r = samaya::test::reference_qp(fixed);
    if (r.status == ReferenceQpResult::Status::kOptimal &&
        (best.status != ReferenceQpResult::Status::kOptimal ||
         (maximize ? r.objective > best.objective : r.objective < best.objective))) {
      best = r;
    }
    std::size_t k = 0;
    for (; k < ints.size(); ++k) {
      if (value[k] < model.col_upper[ints[k]]) {
        value[k] += 1.0;
        break;
      }
      value[k] = model.col_lower[ints[k]];
    }
    if (k == ints.size()) break;
  }
  return best;
}

}  // namespace

TEST(miqp_random_models_match_brute_force) {
  // Convex QPs with 1-2 integer columns (boxed, integral bounds): branch and bound over the QP
  // relaxations must find the brute-force optimum, or prove infeasibility when no integer
  // assignment is feasible; every solution passes the verifier with integrality.
  std::mt19937 rng(75);
  int optimal = 0;
  int infeasible = 0;
  int mismatches = 0;
  for (const QpFamily family : {QpFamily::kStrictlyConvex, QpFamily::kLowRank}) {
    for (int k = 0; k < 100; ++k) {
      Model model = random_qp(family, rng);
      const int ints = std::min(model.num_cols(), 1 + static_cast<int>(rng() % 2));
      for (int j = 0; j < ints; ++j) model.col_type[j] = samaya::VarType::kInteger;
      const ReferenceQpResult ref = brute_force_miqp(model);
      Params params;
      params.log_level = 0;
      params.presolve = false;
      params.mip_rel_gap = 0.0;
      params.mip_abs_gap = 1e-9;
      const Result r = Solver(params).solve(model);
      bool same = false;
      if (ref.status == ReferenceQpResult::Status::kOptimal) {
        ++optimal;
        same = r.status == Status::kOptimal && r.verified &&
               std::fabs(r.objective - ref.objective) <= 1e-6 * (1.0 + std::fabs(ref.objective));
      } else {
        ++infeasible;
        same = r.status == Status::kInfeasible;
      }
      if (!same) {
        ++mismatches;
        std::printf("  mismatch: status %s objective %.10g vs brute force %s %.10g\n",
                    samaya::to_string(r.status), r.objective,
                    ref.status == ReferenceQpResult::Status::kOptimal ? "optimal" : "infeasible",
                    ref.objective);
      }
    }
  }
  std::printf("  %d optimal, %d infeasible, %d mismatches\n", optimal, infeasible, mismatches);
  CHECK_EQ(mismatches, 0);
  CHECK(optimal > 100);
  CHECK(infeasible > 5);
}

TEST(lp_methods_stop_at_the_time_limit) {
  // When the interior point or PDLP stops at the time limit, the solve ends there: no dual
  // simplex fallback runs past the limit (it used to get the whole limit again).
  std::mt19937 rng(76);
  const Model model = samaya::test::random_lp(samaya::test::LpFamily::kFeasible, 300, 300, rng);
  for (const samaya::LpMethod method : {samaya::LpMethod::kBarrier, samaya::LpMethod::kPdlp}) {
    Params params;
    params.log_level = 0;
    params.presolve = false;
    params.lp_method = method;
    params.time_limit = 0.02;
    const Result r = Solver(params).solve(model);
    // Too little time for either method: it stops at the limit (or, on a fast machine, solves the
    // LP itself); in neither case may the dual simplex run.
    CHECK(r.status == Status::kTimeLimit || r.status == Status::kOptimal);
    CHECK_EQ(r.simplex_iterations, 0);
  }
}
