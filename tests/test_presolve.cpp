// Presolve and postsolve against the dense reference simplex, on random LPs with planted
// structure (singleton rows, forcing rows, duplicate rows, free column singletons, fixed and
// empty columns). Postsolve is exercised directly, so a solver fallback cannot hide a bug: the
// restored primal and dual solution must pass the verifier on the original model.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

#include "core/log.hpp"
#include "lp/lp_solver.hpp"
#include "presolve/presolve.hpp"
#include "reference_lp.hpp"
#include "samaya/solver.hpp"
#include "samaya/verify.hpp"
#include "test_framework.hpp"

using samaya::Index;
using samaya::kInf;
using samaya::Model;
using samaya::Presolve;
using samaya::PresolveStatus;
using samaya::test::ReferenceLp;
using samaya::test::ReferenceResult;

namespace {

// Random LP built around a point, with each presolve reduction planted a few times.
Model structured_lp(std::mt19937& rng, bool integer_columns = false) {
  const auto uniform_int = [&](int lo, int hi) {
    return std::uniform_int_distribution<int>(lo, hi)(rng);
  };
  const auto chance = [&](double p) { return std::bernoulli_distribution(p)(rng); };
  const auto real = [&](double lo, double hi) {
    return std::uniform_real_distribution<double>(lo, hi)(rng);
  };

  Model model;
  model.sense = chance(0.3) ? samaya::ObjSense::kMaximize : samaya::ObjSense::kMinimize;
  const int n = uniform_int(3, 24);
  std::vector<double> point(static_cast<std::size_t>(n));
  for (int j = 0; j < n; ++j) {
    const double p = integer_columns ? uniform_int(-3, 6) : std::round(real(-3, 6) * 4) / 4;
    point[j] = p;
    const double r = real(0, 1);
    double lo = -kInf;
    double up = kInf;
    if (r < 0.1) {
      lo = up = p;  // Fixed.
    } else if (r < 0.5) {
      lo = p - uniform_int(0, 3);
      up = p + uniform_int(0, 4);
    } else if (r < 0.75) {
      lo = p - uniform_int(0, 3);
    } else if (r < 0.9) {
      up = p + uniform_int(0, 3);
    }
    model.col_lower.push_back(lo);
    model.col_upper.push_back(up);
    model.obj.push_back(uniform_int(-5, 5));
    model.col_type.push_back(integer_columns && chance(0.5) ? samaya::VarType::kInteger
                                                            : samaya::VarType::kContinuous);
  }

  std::vector<samaya::Triplet> t;
  int m = 0;
  const auto add_row = [&](const std::vector<std::pair<int, double>>& entries, int kind,
                           double slack) {
    double act = 0.0;
    for (const auto& [j, a] : entries) {
      t.push_back({m, j, a});
      act += a * point[j];
    }
    double lo = -kInf;
    double up = kInf;
    switch (kind) {
      case 0: up = act + slack; break;
      case 1: lo = act - slack; break;
      case 2: lo = up = act; break;
      case 3: lo = act - slack; up = act + slack + 1; break;
      default: break;  // Free row.
    }
    model.row_lower.push_back(lo);
    model.row_upper.push_back(up);
    return m++;
  };
  const auto random_entries = [&](int count) {
    std::vector<std::pair<int, double>> e;
    std::vector<int> cols(static_cast<std::size_t>(n));
    for (int j = 0; j < n; ++j) cols[j] = j;
    std::shuffle(cols.begin(), cols.end(), rng);
    for (int k = 0; k < std::min(count, n); ++k) {
      int v = uniform_int(-4, 4);
      if (v == 0) v = 2;
      e.push_back({cols[k], static_cast<double>(v)});
    }
    std::sort(e.begin(), e.end());
    return e;
  };

  // Forcing rows: set the point on the bounds that minimize the activity, then make the row's
  // upper bound that minimum.
  for (int f = uniform_int(0, 2); f > 0; --f) {
    auto e = random_entries(uniform_int(2, 4));
    bool finite = true;
    for (const auto& [j, a] : e) {
      const double b = a > 0 ? model.col_lower[j] : model.col_upper[j];
      finite = finite && std::isfinite(b);
    }
    if (!finite) continue;
    for (const auto& [j, a] : e) point[j] = a > 0 ? model.col_lower[j] : model.col_upper[j];
    add_row(e, 0, 0.0);
  }
  // General rows.
  const int general = uniform_int(1, 14);
  std::vector<std::vector<std::pair<int, double>>> rows;
  for (int k = 0; k < general; ++k) {
    rows.push_back(random_entries(uniform_int(2, 6)));
    add_row(rows.back(), uniform_int(0, 4), chance(0.4) ? 0.0 : real(0, 3));
  }
  // Duplicate rows (scaled copies).
  for (int d = uniform_int(0, 3); d > 0; --d) {
    auto e = rows[uniform_int(0, static_cast<int>(rows.size()) - 1)];
    const double factor = chance(0.5) ? -2.0 : 0.5;
    for (auto& entry : e) entry.second *= factor;
    add_row(e, uniform_int(0, 3), chance(0.5) ? 0.0 : real(0, 2));
  }
  // Singleton rows.
  for (int s = uniform_int(0, 4); s > 0; --s) {
    add_row(random_entries(1), uniform_int(0, 3), chance(0.5) ? 0.0 : real(0, 2));
  }
  // Empty rows.
  if (chance(0.3)) {
    model.row_lower.push_back(chance(0.5) ? -1.0 : -kInf);
    model.row_upper.push_back(1.0);
    ++m;
  }
  // Free column singletons: a new free column in an equality row (any cost) or in an inequality
  // row (no cost).
  for (int s = uniform_int(0, 3); s > 0; --s) {
    const int j = static_cast<int>(model.obj.size());
    const bool equality = chance(0.6);
    const double v = real(-2, 2);
    point.push_back(v);
    model.col_lower.push_back(-kInf);
    model.col_upper.push_back(kInf);
    model.obj.push_back(equality ? uniform_int(-3, 3) : 0.0);
    model.col_type.push_back(samaya::VarType::kContinuous);
    auto e = random_entries(uniform_int(1, 3));
    e.push_back({j, chance(0.5) ? 1.0 : -3.0});
    add_row(e, equality ? 2 : uniform_int(0, 1), real(0, 2));
  }
  // Empty columns.
  if (chance(0.3)) {
    model.col_lower.push_back(-1.0);
    model.col_upper.push_back(chance(0.5) ? 2.0 : kInf);
    model.obj.push_back(uniform_int(0, 3));
    model.col_type.push_back(samaya::VarType::kContinuous);
  }
  model.A = samaya::SparseMatrix::from_triplets(m, static_cast<Index>(model.obj.size()),
                                                std::move(t));
  return model;
}

}  // namespace

TEST(presolve_postsolve_restores_verified_optimal_solutions) {
  std::mt19937 rng(2026);
  const samaya::Logger quiet(0);
  samaya::PresolveStats total;
  int optimal = 0;
  int infeasible = 0;
  int unbounded = 0;
  int failures = 0;
  Index rows_before = 0;
  Index rows_after = 0;
  for (int k = 0; k < 1500; ++k) {
    const Model model = structured_lp(rng);
    const ReferenceResult ref = ReferenceLp(model).solve();
    Presolve presolve(model);
    const PresolveStatus status = presolve.run();
    const samaya::PresolveStats& st = presolve.stats();
    total.singleton_rows += st.singleton_rows;
    total.forcing_rows += st.forcing_rows;
    total.redundant_rows += st.redundant_rows;
    total.duplicate_rows += st.duplicate_rows;
    total.empty_rows += st.empty_rows;
    total.fixed_cols += st.fixed_cols;
    total.empty_cols += st.empty_cols;
    total.dominated_cols += st.dominated_cols;
    total.free_col_singletons += st.free_col_singletons;

    if (ref.status == ReferenceResult::Status::kInfeasible) {
      ++infeasible;
      continue;  // Presolve may or may not notice; the solver falls back either way.
    }
    if (status == PresolveStatus::kInfeasible) {
      ++failures;
      std::fprintf(stderr, "  #%d: presolve claims infeasible, reference %s\n", k,
                   ref.status == ReferenceResult::Status::kOptimal ? "optimal" : "unbounded");
      continue;
    }
    if (ref.status == ReferenceResult::Status::kUnbounded) {
      ++unbounded;
      continue;
    }
    ++optimal;
    rows_before += model.num_rows();
    rows_after += presolve.reduced().num_rows();

    const Model& reduced = presolve.reduced();
    samaya::LpResult lp;
    lp.status = samaya::SimplexStatus::kOptimal;
    if (reduced.num_cols() > 0) {
      lp = samaya::solve_lp(reduced, samaya::LpSolveOptions{}, quiet);
    }
    if (lp.status != samaya::SimplexStatus::kOptimal) {
      ++failures;
      std::fprintf(stderr, "  #%d: reduced model %s, reference optimal\n", k,
                   samaya::to_string(lp.status));
      continue;
    }
    std::vector<double> x;
    std::vector<double> y;
    presolve.postsolve(lp.col_value, lp.row_dual, x, y);
    double objective = model.obj_offset;
    for (Index j = 0; j < model.num_cols(); ++j) objective += model.obj[j] * x[j];
    const samaya::VerifyReport report = samaya::verify_lp_optimality(model, x, y);
    const bool same = std::fabs(objective - ref.objective) <= 1e-6 * (1 + std::fabs(ref.objective));
    if (!report.ok || !same) {
      ++failures;
      std::fprintf(stderr, "  #%d (%dx%d -> %dx%d): objective %.9g vs reference %.9g, %s\n", k,
                   model.num_rows(), model.num_cols(), reduced.num_rows(), reduced.num_cols(),
                   objective, ref.objective, report.message.c_str());
    }
  }
  std::printf("  %d optimal, %d infeasible, %d unbounded; rows %d -> %d; reductions: "
              "%d singleton, %d forcing, %d redundant, %d duplicate, %d empty rows; "
              "%d fixed, %d empty, %d dominated, %d free singleton cols\n",
              optimal, infeasible, unbounded, rows_before, rows_after, total.singleton_rows,
              total.forcing_rows, total.redundant_rows, total.duplicate_rows, total.empty_rows,
              total.fixed_cols, total.empty_cols, total.dominated_cols,
              total.free_col_singletons);
  CHECK_EQ(failures, 0);
  CHECK(optimal > 600);
  CHECK(total.singleton_rows > 100);
  CHECK(total.forcing_rows > 50);
  CHECK(total.redundant_rows > 50);
  CHECK(total.duplicate_rows > 50);
  CHECK(total.empty_rows > 50);
  CHECK(total.fixed_cols > 100);
  CHECK(total.empty_cols > 50);
  CHECK(total.dominated_cols > 50);
  CHECK(total.free_col_singletons > 50);
}

TEST(presolve_solver_results_match_reference_with_and_without_presolve) {
  // Through the public solver: every status (including infeasible and unbounded, which fall back
  // to the original model) must agree with the reference and verify.
  std::mt19937 rng(7);
  int failures = 0;
  for (int k = 0; k < 400; ++k) {
    const Model model = structured_lp(rng);
    const ReferenceResult ref = ReferenceLp(model).solve();
    samaya::Params params;
    params.log_level = 0;
    const samaya::Result result = samaya::Solver(params).solve(model);
    samaya::Status expected = samaya::Status::kOptimal;
    if (ref.status == ReferenceResult::Status::kInfeasible) expected = samaya::Status::kInfeasible;
    if (ref.status == ReferenceResult::Status::kUnbounded) expected = samaya::Status::kUnbounded;
    bool ok = result.status == expected && result.verified;
    if (ok && expected == samaya::Status::kOptimal) {
      ok = std::fabs(result.objective - ref.objective) <= 1e-6 * (1 + std::fabs(ref.objective));
    }
    if (!ok) {
      ++failures;
      std::fprintf(stderr, "  #%d: solver %s %.9g (verified %d), reference objective %.9g\n", k,
                   samaya::to_string(result.status), result.objective, result.verified,
                   ref.objective);
    }
  }
  CHECK_EQ(failures, 0);
}

namespace {

Model forcing_snapshot_lp() {
  Model model;
  model.obj = {-2.0, -1.0, 1.0};
  model.col_lower = {0.0, 0.0, 0.0};
  model.col_upper = {2.0, 2.0, 2.0};
  model.col_type.assign(3, samaya::VarType::kContinuous);
  model.row_lower = {-kInf, 1.0};
  model.row_upper = {0.0, kInf};
  model.A = samaya::SparseMatrix::from_triplets(2, 3,
      {{0, 0, 1.0}, {0, 1, 1.0}, {1, 0, 1.0}, {1, 2, 1.0}});
  return model;
}

}  // namespace

TEST(presolve_forcing_snapshot_precedes_fixing) {
  Model model = forcing_snapshot_lp();
  Presolve p(model);
  REQUIRE(p.run() == PresolveStatus::kReduced);
  std::vector<double> x, y;
  p.postsolve({}, {}, x, y);
  CHECK_EQ(x[2], 1.0);
  CHECK_EQ(y[1], 1.0);
  CHECK_EQ(y[0], -3.0);
  CHECK(samaya::verify_lp_optimality(model, x, y).ok);
  // Postsolve must use the saved coefficients even if the source matrix changes afterwards.
  model.A = samaya::SparseMatrix(2, 3);
  std::vector<double> again_x, again_y;
  p.postsolve({}, {}, again_x, again_y);
  CHECK(again_x == x);
  CHECK(again_y == y);
}

TEST(presolve_singleton_snapshot_independent_of_source_matrix) {
  Model model;
  model.obj = {1.0};
  model.col_lower = {0.0};
  model.col_upper = {4.0};
  model.col_type = {samaya::VarType::kContinuous};
  model.row_lower = {1.0, 2.0};
  model.row_upper = {kInf, kInf};
  model.A = samaya::SparseMatrix::from_triplets(2, 1, {{0, 0, 1.0}, {1, 0, 1.0}});
  Presolve p(model);
  REQUIRE(p.run() == PresolveStatus::kReduced);
  model.A = samaya::SparseMatrix(2, 1);
  std::vector<double> x, y;
  p.postsolve({}, {}, x, y);
  CHECK_EQ(x[0], 2.0);
  CHECK_EQ(y[0], 0.0);
  CHECK_EQ(y[1], 1.0);
}
