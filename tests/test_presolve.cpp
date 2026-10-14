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
#include "lp_generators.hpp"
#include "presolve/presolve.hpp"
#include "reference_lp.hpp"
#include "reference_milp.hpp"
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

namespace samaya {

// Exercise C3 alone for atomic rejection and boundary cases, without earlier reductions.
struct PresolveSubstitutionTestAccess {
  static bool substitutions(Presolve& p) {
    bool changed = false;
    p.process_substitutions(changed);
    p.build_reduced();
    return changed;
  }
};

}  // namespace samaya

TEST(presolve_near_parallel_equations_remain_feasible) {
  const double eps = std::ldexp(1.0, -41);
  const double v = 1.0 / eps;
  Model model;
  model.obj = {0.0, 0.0};
  model.col_lower = {-2.0 * v, 0.0};
  model.col_upper = {0.0, 2.0 * v};
  model.col_type.assign(2, samaya::VarType::kContinuous);
  model.row_lower = model.row_upper = {0.0, 1.0};
  model.A = samaya::SparseMatrix::from_triplets(2, 2,
      {{0, 0, 1.0}, {0, 1, 1.0}, {1, 0, 1.0}, {1, 1, 1.0 + eps}});
  Presolve p(model);
  REQUIRE(p.run() == PresolveStatus::kReduced);
  CHECK_EQ(p.stats().duplicate_rows, 0);
  CHECK_EQ(p.reduced().num_rows(), 2);
  const std::vector<double> witness = {-v, v};
  CHECK(samaya::verify_primal(p.reduced(), witness).ok);
  // The reference's absolute 1e-9 pivot threshold needs equivalent row scaling to see eps.
  Model scaled = model;
  const double scale = std::ldexp(1.0, 20);
  scaled.row_lower = scaled.row_upper = {0.0, scale};
  scaled.A = samaya::SparseMatrix::from_triplets(2, 2,
      {{0, 0, scale}, {0, 1, scale}, {1, 0, scale}, {1, 1, scale * (1.0 + eps)}});
  const auto ref = ReferenceLp(scaled).solve();
  REQUIRE(ref.status == ReferenceResult::Status::kOptimal);
  CHECK_NEAR(ref.objective, p.reduced().obj_offset, 0.0);
  CHECK(samaya::verify_primal(model, ref.x).ok);
}

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

namespace {

Model implied_singleton_lp(double a = 1.0) {
  Model model;
  model.obj = {1.0, 0.0};
  model.col_lower = {0.0, 0.0};
  model.col_upper = {100.0, 4.0};
  model.col_type.assign(2, samaya::VarType::kContinuous);
  model.row_lower = model.row_upper = {a > 0.0 ? 10.0 : -10.0};
  model.A = samaya::SparseMatrix::from_triplets(1, 2, {{0, 0, a}, {0, 1, 1.0}});
  return model;
}

void check_implied_postsolve(const Model& model, const Presolve& p) {
  samaya::LpResult lp;
  lp.status = samaya::SimplexStatus::kOptimal;
  if (p.reduced().num_cols() > 0) {
    lp = samaya::solve_lp(p.reduced(), samaya::LpSolveOptions{}, samaya::Logger(0));
  }
  REQUIRE(lp.status == samaya::SimplexStatus::kOptimal);
  std::vector<double> x, y;
  p.postsolve(lp.col_value, lp.row_dual, x, y);
  CHECK(samaya::verify_lp_optimality(model, x, y).ok);
  const auto ref = ReferenceLp(model).solve();
  REQUIRE(ref.status == ReferenceResult::Status::kOptimal);
  double objective = model.obj_offset;
  for (Index j = 0; j < model.num_cols(); ++j) objective += model.obj[j] * x[j];
  CHECK_NEAR(objective, ref.objective, 1e-8);
}

// A small bounded family around a feasible point. The first three continuous columns occur
// only in their planted rows and have bounds strictly outside the hand-computed implied range.
Model redundant_singleton_family(std::mt19937& rng, bool mip, int variant) {
  const auto integer = [&](int lo, int hi) {
    return std::uniform_int_distribution<int>(lo, hi)(rng);
  };
  const auto real = [&](double lo, double hi) {
    return std::uniform_real_distribution<double>(lo, hi)(rng);
  };
  Model model;
  model.sense = integer(0, 1) ? samaya::ObjSense::kMinimize : samaya::ObjSense::kMaximize;
  model.obj_offset = real(-2.0, 2.0);
  model.obj.resize(6);
  model.col_lower.resize(6);
  model.col_upper.resize(6);
  model.col_type.assign(6, samaya::VarType::kContinuous);
  std::vector<double> point(6);
  for (Index j = 3; j < 6; ++j) {
    point[j] = integer(0, 3);
    model.col_lower[j] = point[j] - 2.0;
    model.col_upper[j] = point[j] + 2.0;
    model.obj[j] = real(-3.0, 3.0);
    if (mip && j != 5) model.col_type[j] = samaya::VarType::kInteger;
  }
  std::vector<samaya::Triplet> entries;
  for (Index i = 0; i < 3; ++i) {
    const double a = integer(0, 1) ? 2.0 : -3.0;
    point[i] = real(-3.0, 3.0);
    const bool equality = i != 2;
    model.obj[i] = equality ? real(-3.0, 3.0) : 0.0;
    double activity = a * point[i];
    double lower = 0.0, upper = 0.0;
    entries.push_back({i, i, a});
    for (Index j = 3; j < 6; ++j) {
      double b = integer(-3, 3);
      if (b == 0.0) b = 1.0;
      entries.push_back({i, j, b});
      activity += b * point[j];
      lower += b * (b > 0.0 ? model.col_lower[j] : model.col_upper[j]);
      upper += b * (b > 0.0 ? model.col_upper[j] : model.col_lower[j]);
    }
    const double rl = activity - (equality ? 0.0 : 1.0);
    const double ru = activity + (equality ? 0.0 : 2.0);
    model.row_lower.push_back(rl);
    model.row_upper.push_back(ru);
    double lo = (rl - upper) / a;
    double up = (ru - lower) / a;
    if (a < 0.0) std::swap(lo, up);
    model.col_lower[i] = lo - real(1.0, 3.0);
    model.col_upper[i] = up + real(1.0, 3.0);
  }
  // Retain a coupled core after singleton substitution; include infeasible cases as well as
  // feasible equalities and inequalities, with fractional LP optima in the integer core.
  const double activity = 2.0 * point[3] + point[4] - point[5];
  model.row_lower.push_back(variant == 1 ? activity : -kInf);
  model.row_upper.push_back(variant == 2 ? activity - 20.0 : activity + (variant == 1 ? 0.0 : 0.5));
  entries.insert(entries.end(), {{3, 3, 2.0}, {3, 4, 1.0}, {3, 5, -1.0}});
  model.A = samaya::SparseMatrix::from_triplets(4, 6, std::move(entries));
  return model;
}

}  // namespace

TEST(presolve_implied_free_equality_both_pivot_signs) {
  for (const double a : {1.0, -2.0}) {
    const Model model = implied_singleton_lp(a);  // [6, 10] or [5, 7].
    Presolve p(model);
    REQUIRE(p.run() == PresolveStatus::kReduced);
    CHECK_EQ(p.stats().implied_free_col_singletons, 1);
    CHECK_EQ(p.stats().free_col_singletons, 0);
    check_implied_postsolve(model, p);
  }
}

TEST(presolve_implied_free_rejects_nonredundant_bound) {
  Model model = implied_singleton_lp();
  model.col_lower[0] = 7.0;  // The row allows x = 6.
  Presolve p(model);
  REQUIRE(p.run() == PresolveStatus::kReduced);
  CHECK_EQ(p.stats().implied_free_col_singletons, 0);
  CHECK_EQ(p.stats().implied_free_rejected_by_tolerance, 0);
  CHECK(std::find(p.col_map().begin(), p.col_map().end(), 0) != p.col_map().end());
  check_implied_postsolve(model, p);
}

TEST(presolve_implied_free_inequality_slack) {
  for (const double a : {1.0, -2.0}) {
    Model model = implied_singleton_lp(a);
    model.obj = {0.0, -1.0};
    model.row_upper[0] += 2.0;  // Finite implied range, zero-cost slack.
    Presolve p(model);
    REQUIRE(p.run() == PresolveStatus::kReduced);
    CHECK_EQ(p.stats().implied_free_col_singletons, 1);
    check_implied_postsolve(model, p);
    std::vector<double> x, y;
    p.postsolve({}, {}, x, y);
    CHECK(x[0] >= model.col_lower[0]);
    CHECK(x[0] <= model.col_upper[0]);
  }
}

TEST(presolve_implied_free_infinite_contributions) {
  for (const int count : {1, 2}) {
    for (const double a : {1.0, -2.0}) {
      Model model = implied_singleton_lp(a);
      model.col_upper[1] = kInf;
      if (count == 2) {
        model.obj.push_back(0.0);
        model.col_lower.push_back(0.0);
        model.col_upper.push_back(kInf);
        model.col_type.push_back(samaya::VarType::kContinuous);
        model.A = samaya::SparseMatrix::from_triplets(1, 3,
            {{0, 0, a}, {0, 1, 1.0}, {0, 2, 1.0}});
      }
      Presolve p(model);
      REQUIRE(p.run() == PresolveStatus::kReduced);
      CHECK_EQ(p.stats().implied_free_col_singletons, 0);
      CHECK_EQ(p.stats().implied_free_rejected_by_tolerance, 0);
      check_implied_postsolve(model, p);
    }
  }
  // No proof is needed on an already infinite column bound, even when that activity side
  // has an infinite contribution. The finite upper side still has to be proved.
  Model model = implied_singleton_lp();
  model.col_lower[0] = -kInf;
  model.col_upper[1] = kInf;
  samaya::PresolveOptions options;
  options.max_passes = 1;
  Presolve p(model, options);
  REQUIRE(p.run() == PresolveStatus::kReduced);
  CHECK_EQ(p.stats().implied_free_col_singletons, 1);
  std::vector<double> x, y;
  p.postsolve({20.0}, {}, x, y);
  CHECK_EQ(x[0], -10.0);
  CHECK(x[0] <= model.col_upper[0]);
}

TEST(presolve_implied_free_later_bound_tightening) {
  Model model = implied_singleton_lp();
  model.col_upper[1] = 6.0;  // Avoid forcing before x is substituted.
  model.obj.push_back(0.0);
  model.col_lower.push_back(0.0);
  model.col_upper.push_back(0.0);
  model.col_type.push_back(samaya::VarType::kContinuous);
  model.obj.push_back(0.0);
  model.col_lower.push_back(0.0);
  model.col_upper.push_back(3.0);
  model.col_type.push_back(samaya::VarType::kContinuous);
  model.row_lower.push_back(4.0);
  model.row_upper.push_back(4.0);
  // Keep y in two rows after x's removal, so the next row pass tightens it to 4.
  model.row_lower.push_back(5.0);
  model.row_upper.push_back(5.0);
  model.A = samaya::SparseMatrix::from_triplets(3, 4,
      {{0, 0, 1.0}, {0, 1, 1.0}, {1, 1, 1.0}, {1, 2, 1.0},
       {2, 1, 1.0}, {2, 3, 1.0}});
  // C3 would substitute column 1 out of its equations first; this test isolates C1.
  samaya::PresolveOptions first_pass;
  first_pass.max_passes = 1;
  first_pass.substitute = false;
  Presolve early(model, first_pass);
  REQUIRE(early.run() == PresolveStatus::kReduced);
  CHECK_EQ(early.stats().implied_free_col_singletons, 1);
  CHECK_EQ(early.stats().bounds_tightened, 0);
  CHECK(std::find(early.col_map().begin(), early.col_map().end(), 0) == early.col_map().end());
  CHECK(std::find(early.col_map().begin(), early.col_map().end(), 1) != early.col_map().end());
  samaya::PresolveOptions c1_only;
  c1_only.substitute = false;
  Presolve p(model, c1_only);
  REQUIRE(p.run() == PresolveStatus::kReduced);
  CHECK_EQ(p.stats().implied_free_col_singletons, 2);  // x first, the last column after y is fixed.
  CHECK(p.stats().passes >= 2);
  CHECK(p.stats().bounds_tightened > 0);
  check_implied_postsolve(model, p);
}

TEST(presolve_implied_free_rechecks_tightened_candidate) {
  Model model = implied_singleton_lp();
  model.row_lower.push_back(7.0);
  model.row_upper.push_back(kInf);
  model.A = samaya::SparseMatrix::from_triplets(2, 2,
      {{0, 0, 1.0}, {0, 1, 1.0}, {1, 0, 1.0}});
  Presolve p(model);
  REQUIRE(p.run() == PresolveStatus::kReduced);
  CHECK(p.stats().bounds_tightened > 0);
  CHECK_EQ(p.stats().implied_free_col_singletons, 0);
  CHECK(std::find(p.col_map().begin(), p.col_map().end(), 0) != p.col_map().end());
  check_implied_postsolve(model, p);
}

TEST(presolve_implied_free_one_sided_row_and_column_bounds) {
  Model model = implied_singleton_lp();
  model.obj = {0.0, -1.0};
  model.col_lower[0] = -kInf;
  model.row_lower[0] = -kInf;
  Presolve p(model);
  REQUIRE(p.run() == PresolveStatus::kReduced);
  CHECK_EQ(p.stats().implied_free_col_singletons, 1);
  check_implied_postsolve(model, p);
  // With a negative pivot, the infinite row side instead leaves the upper bound unproved.
  model.col_lower[0] = -100.0;
  model.A = samaya::SparseMatrix::from_triplets(1, 2, {{0, 0, -1.0}, {0, 1, 1.0}});
  model.col_upper[0] = kInf;
  Presolve negative(model);
  REQUIRE(negative.run() == PresolveStatus::kReduced);
  // Here the row implies x >= -10, so the finite lower bound is redundant.
  CHECK_EQ(negative.stats().implied_free_col_singletons, 1);
  check_implied_postsolve(model, negative);
}

TEST(presolve_implied_free_error_allowance) {
  for (const bool lower : {true, false}) {
    Model model = implied_singleton_lp();
    if (lower) {
      model.col_lower[0] = 6.0 - 6e-13;
    } else {
      model.col_upper[0] = 10.0 + 1e-12;
    }
    Presolve p(model);
    REQUIRE(p.run() == PresolveStatus::kReduced);
    CHECK_EQ(p.stats().implied_free_col_singletons, 0);
    CHECK(p.stats().implied_free_rejected_by_tolerance > 0);
    check_implied_postsolve(model, p);
  }
}

namespace {
void check_doubleton_reference(const Model& model, const Presolve& p, PresolveStatus status);
}  // namespace

TEST(presolve_implied_free_integer_columns_excluded_in_mip) {
  Model model = implied_singleton_lp();
  model.col_type[0] = samaya::VarType::kInteger;
  samaya::PresolveOptions options;
  options.mip = true;
  Presolve p(model, options);
  REQUIRE(p.run() == PresolveStatus::kReduced);
  CHECK_EQ(p.stats().implied_free_col_singletons, 0);
  // C1 keeps the integer column; only C2 (exact integer data) may still eliminate it.
  const bool kept =
      std::find(p.col_map().begin(), p.col_map().end(), 0) != p.col_map().end();
  CHECK(kept || p.stats().doubleton_integer_cols + p.stats().doubleton_continuous_cols > 0);
  check_doubleton_reference(model, p, PresolveStatus::kReduced);
}

TEST(presolve_implied_free_differential_with_and_without_presolve) {
  std::mt19937 rng(130102026);
  Index lp_substitutions = 0, mip_substitutions = 0;
  for (const bool mip : {false, true}) {
    for (int k = 0; k < 300; ++k) {
      const Model model = redundant_singleton_family(rng, mip, k % 3);
      samaya::Status expected;
      double objective;
      if (mip) {
        const auto ref = samaya::test::reference_milp(model);
        REQUIRE(ref.status == samaya::test::ReferenceMilpResult::Status::kOptimal ||
                ref.status == samaya::test::ReferenceMilpResult::Status::kInfeasible);
        expected = ref.status == samaya::test::ReferenceMilpResult::Status::kOptimal
            ? samaya::Status::kOptimal : samaya::Status::kInfeasible;
        objective = ref.objective;
      } else {
        const auto ref = ReferenceLp(model).solve();
        REQUIRE(ref.status == ReferenceResult::Status::kOptimal ||
                ref.status == ReferenceResult::Status::kInfeasible);
        expected = ref.status == ReferenceResult::Status::kOptimal
            ? samaya::Status::kOptimal : samaya::Status::kInfeasible;
        objective = ref.objective;
      }
      samaya::PresolveOptions options;
      options.mip = mip;
      Presolve p(model, options);
      const auto status = p.run();
      (mip ? mip_substitutions : lp_substitutions) += p.stats().implied_free_col_singletons;
      if (expected == samaya::Status::kOptimal) {
        REQUIRE(status == PresolveStatus::kReduced);
        if (!mip) check_implied_postsolve(model, p);  // Bypass public-solver fallbacks.
      }
      for (const bool enabled : {false, true}) {
        samaya::Params params;
        params.log_level = 0;
        params.presolve = enabled;
        const auto result = samaya::Solver(params).solve(model);
        CHECK(result.status == expected);
        // MIP infeasibility is a search/presolve conclusion, without an LP certificate.
        if (!mip || expected == samaya::Status::kOptimal) CHECK(result.verified);
        if (expected == samaya::Status::kOptimal) {
          CHECK_NEAR(result.objective, objective, 1e-6 * (1.0 + std::fabs(objective)));
          if (!mip) {
            CHECK(samaya::verify_lp_optimality(model, result.col_value, result.row_dual).ok);
          }
        }
      }
    }
  }
  std::printf("  implied-free differential: %d LP, %d MILP substitutions (%d total)\n",
              lp_substitutions, mip_substitutions, lp_substitutions + mip_substitutions);
  CHECK(lp_substitutions >= 100);
  CHECK(mip_substitutions >= 100);
}

namespace {

Model doubleton_model(bool integer_x = false, double a = 1.0, double b = 2.0,
                      double c = 0.0) {
  Model model;
  model.obj = {1.0, -2.0, -0.5};
  model.obj_offset = 1.25;
  model.col_lower = {integer_x ? -7.0 : -2.5, -3.0, -2.0};
  model.col_upper = {integer_x ? 7.0 : 2.5, 3.0, 2.0};
  model.col_type = {integer_x ? samaya::VarType::kInteger : samaya::VarType::kContinuous,
                    samaya::VarType::kInteger, samaya::VarType::kContinuous};
  model.row_lower = {c, -kInf, -kInf, -3.0};
  model.row_upper = {c, 3.0, 3.0, kInf};
  // x has fewer other entries than y. No column singleton can mask C2.
  model.A = samaya::SparseMatrix::from_triplets(4, 3,
      {{0, 0, a}, {0, 1, b}, {1, 0, 1.0}, {1, 2, 1.0},
       {2, 1, 1.0}, {2, 2, 1.0}, {3, 1, 1.0}, {3, 2, -1.0}});
  return model;
}

samaya::PresolveOptions doubleton_options(int passes = 50) {
  samaya::PresolveOptions options;
  options.mip = true;
  options.max_passes = passes;
  return options;
}

void check_doubleton_reference(const Model& model, const Presolve& p, PresolveStatus status) {
  using RefStatus = samaya::test::ReferenceMilpResult::Status;
  const auto ref = samaya::test::reference_milp(model);
  REQUIRE(ref.status == RefStatus::kOptimal || ref.status == RefStatus::kInfeasible);
  if (status == PresolveStatus::kInfeasible) {
    CHECK(ref.status == RefStatus::kInfeasible);
    return;
  }
  const auto reduced_ref = samaya::test::reference_milp(p.reduced());
  REQUIRE(reduced_ref.status == ref.status);
  if (ref.status != RefStatus::kOptimal) return;
  CHECK_NEAR(reduced_ref.objective, ref.objective, 1e-7 * (1.0 + std::fabs(ref.objective)));
  std::vector<double> x, y;
  p.postsolve(reduced_ref.x, {}, x, y);
  CHECK(samaya::verify_primal(model, x).ok);
  double objective = model.obj_offset;
  for (Index j = 0; j < model.num_cols(); ++j) objective += model.obj[j] * x[j];
  CHECK_NEAR(objective, ref.objective, 1e-7 * (1.0 + std::fabs(ref.objective)));
}

}  // namespace

TEST(presolve_doubleton_both_signs_and_objective_senses) {
  for (const bool integer_x : {false, true}) {
    for (const double a : {-1.0, 1.0}) {
      for (const double b : {-2.0, 2.0}) {
        for (const auto sense : {samaya::ObjSense::kMinimize, samaya::ObjSense::kMaximize}) {
          const Model model = [&] {
            Model result = doubleton_model(integer_x, a, b, 2.0);
            result.sense = sense;
            return result;
          }();
          Presolve p(model, doubleton_options());
          const auto status = p.run();
          REQUIRE(status == PresolveStatus::kReduced);
          CHECK((integer_x ? p.stats().doubleton_integer_cols :
                            p.stats().doubleton_continuous_cols) > 0);
          check_doubleton_reference(model, p, status);
        }
      }
    }
  }
}

TEST(presolve_doubleton_exact_integer_divisibility) {
  for (const double b : {4.0, 3.0, 4.0 + 1e-10}) {
    Model model = doubleton_model(true, 2.0, b, 2.0);
    model.col_lower[0] = -8.0;
    model.col_upper[0] = 8.0;
    Presolve p(model, doubleton_options(1));
    const auto status = p.run();
    REQUIRE(status == PresolveStatus::kReduced);
    if (b == 4.0) CHECK_EQ(p.stats().doubleton_integer_cols, 1);
    else {
      CHECK_EQ(p.stats().doubleton_integer_cols, 0);
      CHECK((b == 3.0 ? p.stats().doubleton_rejected_divisibility :
                       p.stats().doubleton_rejected_integer_data) > 0);
    }
    check_doubleton_reference(model, p, status);
  }
  for (const double c : {1.0, 2.0 + 1e-10, 0x1p32}) {
    const Model model = doubleton_model(true, 2.0, 4.0, c);
    Presolve p(model, doubleton_options(1));
    const auto status = p.run();
    CHECK_EQ(p.stats().doubleton_integer_cols, 0);
    check_doubleton_reference(model, p, status);
  }
}

TEST(presolve_doubleton_integer_transfer_and_ambiguity) {
  // The review's ratio 1000 is refused by the pivot guard. Ratio 100 passes it and must
  // prove the same integer-transfer contradiction without set_bounds' 1e-6 allowance.
  for (const double b : {1000.0, 100.0}) {
    Model model = doubleton_model(false, 1.0, b, b - 0.0005);
    model.col_lower[0] = 0.0;
    model.col_upper[0] = 0.1;
    model.col_lower[1] = 0.0;
    model.col_upper[1] = 2.0;
    Presolve p(model, doubleton_options(1));
    const auto status = p.run();
    if (b == 1000.0) CHECK_EQ(p.stats().doubleton_rejected_pivot, 1);
    else CHECK(status == PresolveStatus::kInfeasible);
    check_doubleton_reference(model, p, status);
  }
  // x = -2y lands exactly on x's bound -6 at y = 3: computed exactly, so it is decided without
  // an allowance and y keeps 3.
  Model model = doubleton_model(true);
  model.col_lower[0] = -6.0;
  model.col_upper[0] = 6.0;
  Presolve exact(model, doubleton_options(1));
  REQUIRE(exact.run() == PresolveStatus::kReduced);
  CHECK_EQ(exact.stats().doubleton_rejected_integer_bound, 0);
  CHECK_EQ(exact.stats().doubleton_integer_cols, 1);
  check_doubleton_reference(model, exact, PresolveStatus::kReduced);
  // One ulp below -6, x = -6 is inside and y keeps 3. One ulp above, x = -6 is outside only
  // within the feasibility tolerance: cutting y = 3 off is not proven, so C2 skips.
  for (const double direction : {-kInf, kInf}) {
    model.col_lower[0] = std::nextafter(-6.0, direction);
    model.col_type[0] = samaya::VarType::kContinuous;
    Presolve p(model, doubleton_options(1));
    REQUIRE(p.run() == PresolveStatus::kReduced);
    if (direction < 0.0) {
      CHECK_EQ(p.stats().doubleton_rejected_integer_bound, 0);
      CHECK_EQ(p.stats().doubleton_continuous_cols, 1);
      REQUIRE(p.reduced().num_cols() == 2);
      CHECK_EQ(p.reduced().col_upper[0], 3.0);
    } else {
      CHECK_EQ(p.stats().doubleton_rejected_integer_bound, 1);
      CHECK_EQ(p.stats().doubleton_continuous_cols, 0);
    }
    check_doubleton_reference(model, p, PresolveStatus::kReduced);
  }
}

TEST(presolve_doubleton_integer_endpoint_on_exact_bound) {
  // x + y = 0, x continuous in [-2, 2], y integer: y's transferred range [-2, 2] has integral
  // endpoints where x sits exactly on its bounds. They must be kept (a rounding that started
  // at the computed endpoint could step past them), and the optimum (y = 2, x = -2) uses one.
  Model model = doubleton_model(false, 1.0, 1.0, 0.0);
  model.col_lower[0] = -2.0;
  model.col_upper[0] = 2.0;
  Presolve p(model, doubleton_options(1));
  REQUIRE(p.run() == PresolveStatus::kReduced);
  CHECK_EQ(p.stats().doubleton_continuous_cols, 1);
  CHECK_EQ(p.stats().doubleton_rejected_integer_bound, 0);
  REQUIRE(p.reduced().num_cols() == 2);
  CHECK_EQ(p.reduced().col_lower[0], -2.0);
  CHECK_EQ(p.reduced().col_upper[0], 2.0);
  check_doubleton_reference(model, p, PresolveStatus::kReduced);
}

TEST(presolve_doubleton_crossing_within_tolerance_is_not_infeasible) {
  // x + 2y = 1, x in [0, 1 - 1e-12], y integer in [0, 1]: no integer y fits in exact
  // arithmetic (y would be in [5e-13, 0.5]), but y = 0, x = 1 is feasible within the
  // tolerance. Presolve must skip it, not prove infeasibility. (A continuous crossing makes the
  // row forcing and never reaches C2; integer rounding does.)
  Model model = doubleton_model(false, 1.0, 2.0, 1.0);
  model.col_lower[0] = 0.0;
  model.col_upper[0] = 1.0 - 1e-12;
  model.col_lower[1] = 0.0;
  model.col_upper[1] = 1.0;
  Presolve p(model, doubleton_options(1));
  REQUIRE(p.run() == PresolveStatus::kReduced);
  CHECK_EQ(p.stats().doubleton_continuous_cols, 0);
  CHECK(p.stats().doubleton_rejected_integer_bound + p.stats().doubleton_rejected_crossed_bounds >
        0);
  check_doubleton_reference(model, p, PresolveStatus::kReduced);
}

TEST(presolve_doubleton_integer_data_magnitude_limit) {
  for (const double a : {-0x1p31, 0x1p31, 0x1p32}) {
    const Model model = doubleton_model(true, a, a);
    Presolve p(model, doubleton_options(1));
    const auto status = p.run();
    REQUIRE(status == PresolveStatus::kReduced);
    if (std::fabs(a) <= 0x1p31) CHECK_EQ(p.stats().doubleton_integer_cols, 1);
    else {
      CHECK_EQ(p.stats().doubleton_integer_cols, 0);
      CHECK_EQ(p.stats().doubleton_rejected_integer_data, 1);
    }
    check_doubleton_reference(model, p, status);
  }
}

TEST(presolve_doubleton_nonfinite_objective_is_rejected_atomically) {
  Model model = doubleton_model();
  model.obj[0] = 1e308;  // The staged cost transfer overflows; the matrix changes are finite.
  Presolve p(model, doubleton_options(1));
  REQUIRE(p.run() == PresolveStatus::kReduced);
  CHECK_EQ(p.stats().doubleton_rejected_nonfinite, 1);
  CHECK_EQ(p.stats().doubleton_continuous_cols, 0);
  CHECK(p.reduced().obj == model.obj);
  CHECK_EQ(p.reduced().obj_offset, model.obj_offset);
  CHECK_EQ(p.reduced().A.nnz(), model.A.nnz());
  CHECK(p.reduced().col_lower == model.col_lower);
  CHECK(p.reduced().col_upper == model.col_upper);
  CHECK(p.reduced().row_lower == model.row_lower);
  CHECK(p.reduced().row_upper == model.row_upper);
}

TEST(presolve_doubleton_pivot_choice) {
  // Continuous beats fewer entries and a larger integer coefficient.
  Model continuous = doubleton_model(false, 1.0, 2.0);
  continuous.A = samaya::SparseMatrix::from_triplets(4, 3,
      {{0, 0, 1.0}, {0, 1, 2.0}, {1, 0, 1.0}, {1, 2, 1.0},
       {2, 1, 1.0}, {2, 2, 1.0}, {3, 0, 1.0}, {3, 2, -1.0}});
  Presolve p(continuous, doubleton_options(1));
  REQUIRE(p.run() == PresolveStatus::kReduced);
  CHECK_EQ(p.stats().doubleton_continuous_cols, 1);
  CHECK(std::find(p.col_map().begin(), p.col_map().end(), 0) == p.col_map().end());
  check_doubleton_reference(continuous, p, PresolveStatus::kReduced);
  // Equal counts choose the larger coefficient. Eliminating y (coefficient 2) is indivisible;
  // eliminating x (coefficient 1) would pass, but is not the specified choice.
  Model integer = doubleton_model(true);
  integer.row_upper[3] = 3.0;
  integer.row_lower[3] = -kInf;
  integer.A = samaya::SparseMatrix::from_triplets(4, 3,
      {{0, 0, 1.0}, {0, 1, 2.0}, {1, 0, 1.0}, {1, 2, 1.0},
       {2, 1, 1.0}, {2, 2, 1.0}, {3, 0, 1.0}, {3, 1, -1.0}, {3, 2, -1.0}});
  Presolve q(integer, doubleton_options(1));
  REQUIRE(q.run() == PresolveStatus::kReduced);
  CHECK_EQ(q.stats().doubleton_integer_cols, 0);
  CHECK_EQ(q.stats().doubleton_rejected_divisibility, 1);
  check_doubleton_reference(integer, q, PresolveStatus::kReduced);
}

TEST(presolve_doubleton_eps_cancellation_is_rejected_atomically) {
  const double eps = std::ldexp(1.0, -41);
  Model model = doubleton_model(false, 1.0, -1.0);
  model.col_lower = {-1.0, 0.0, -2.0};
  model.col_upper = {0x1p43, 0x1p42, 2.0};
  model.row_lower = {0.0, 1.0, -kInf, -kInf};
  model.row_upper = {0.0, 1.0, 0x1p41, 0x1p41};
  model.A = samaya::SparseMatrix::from_triplets(4, 3,
      {{0, 0, 1.0}, {0, 1, -1.0}, {1, 0, 1.0}, {1, 1, -1.0 + eps},
       {1, 2, 1.0}, {2, 1, 1.0}, {2, 2, 1.0}, {3, 1, 1.0}, {3, 2, -1.0}});
  Presolve p(model, doubleton_options(1));
  REQUIRE(p.run() == PresolveStatus::kReduced);
  CHECK_EQ(p.stats().doubleton_rejected_residual, 1);
  CHECK_EQ(p.stats().doubleton_continuous_cols, 0);
  CHECK_EQ(p.reduced().A.nnz(), model.A.nnz());
  const std::vector<double> witness = {0x1p41, 0x1p41, 0.0};
  CHECK(samaya::verify_primal(model, witness).ok);
  CHECK(samaya::verify_primal(p.reduced(), witness).ok);
  CHECK_EQ(p.reduced().obj_offset, model.obj_offset);
  CHECK(p.reduced().obj == model.obj);
}

TEST(presolve_doubleton_lp_mode_does_not_fire) {
  Model model = doubleton_model();
  model.col_type.assign(3, samaya::VarType::kContinuous);
  samaya::PresolveOptions options;
  options.max_passes = 1;
  // C3 may substitute an implied-free column of this LP; this test isolates C2.
  options.substitute = false;
  Presolve p(model, options);
  REQUIRE(p.run() == PresolveStatus::kReduced);
  CHECK_EQ(p.stats().doubleton_continuous_cols, 0);
  CHECK_EQ(p.stats().doubleton_integer_cols, 0);
  CHECK_EQ(p.reduced().num_cols(), 3);
  check_implied_postsolve(model, p);
  // With C3 on, C2 still stays off in LP mode and the reduced model still postsolves.
  Presolve with_c3(model, samaya::PresolveOptions{.max_passes = 1});
  REQUIRE(with_c3.run() == PresolveStatus::kReduced);
  CHECK_EQ(with_c3.stats().doubleton_continuous_cols + with_c3.stats().doubleton_integer_cols, 0);
  check_implied_postsolve(model, with_c3);
}

TEST(presolve_doubleton_survivor_tightened_fixed_or_eliminated_later) {
  for (const int variant : {0, 1, 2}) {
    Model model = doubleton_model(false, 1.0, -2.0, 1.0);
    model.row_lower[1] = variant == 2 ? -1.0 : -kInf;
    model.row_upper[1] = variant == 2 ? -1.0 : 3.0;
    model.row_lower[2] = variant == 0 ? 0.5 : -kInf;
    model.row_upper[2] = variant == 0 ? 3.0 : -1.0;
    model.row_lower[3] = variant == 0 ? -3.0 : -1.0;
    model.row_upper[3] = variant == 0 ? kInf : -1.0;
    Presolve early(model, doubleton_options(1));
    REQUIRE(early.run() == PresolveStatus::kReduced);
    CHECK_EQ(early.stats().doubleton_continuous_cols, variant == 0 ? 1 : 2);
    Presolve p(model, doubleton_options());
    const auto status = p.run();
    REQUIRE(status == PresolveStatus::kReduced);
    CHECK(p.stats().passes >= 2);
    CHECK(p.stats().bounds_tightened > 0);
    if (variant != 0) CHECK(p.stats().cols_removed > early.stats().cols_removed);
    check_doubleton_reference(model, p, status);
  }
}

TEST(presolve_doubleton_forcing_and_duplicates_before_and_after) {
  for (const bool before : {false, true}) {
    for (const bool forcing : {false, true}) {
      Model model = doubleton_model(false, 1.0, -2.0, 1.0);
      std::vector<samaya::Triplet> entries = {
        {0, 0, 1.0}, {0, 1, -2.0}, {1, 0, 1.0}, {1, 2, 1.0},
        {2, 1, 1.0}, {2, 2, 1.0}, {3, 1, 1.0}, {3, 2, -1.0}};
      if (before && forcing) {
        for (int j = 0; j < 2; ++j) {
          model.obj.push_back(0.0);
          model.col_lower.push_back(0.0);
          model.col_upper.push_back(2.0);
          model.col_type.push_back(samaya::VarType::kContinuous);
        }
        model.row_lower.push_back(-kInf);
        model.row_upper.push_back(0.0);
        entries.insert(entries.end(), {{4, 3, 1.0}, {4, 4, 1.0}});
      } else if (before) {
        model.row_lower.push_back(2.0);
        model.row_upper.push_back(2.0);
        entries.insert(entries.end(), {{4, 0, 2.0}, {4, 1, -4.0}});
      } else if (forcing) {
        // Transfer fixes y = 1; the updated first inequality then forces z = 0 next round.
        model.col_lower[0] = 2.5;
        model.col_upper[0] = 4.5;
        model.col_lower[2] = 0.0;
      } else {
        // The first inequality becomes 2y + z <= 2, parallel to row 2 only after C2.
        entries[4].value = 2.0;
        model.col_upper[2] = 3.0;
        model.row_lower[2] = -1.0;
        model.row_upper[2] = kInf;
      }
      model.A = samaya::SparseMatrix::from_triplets(model.num_rows(), model.num_cols(), entries);
      Presolve p(model, doubleton_options());
      const auto status = p.run();
      check_doubleton_reference(model, p, status);
      CHECK(p.stats().doubleton_continuous_cols > 0);
      if (forcing) CHECK(p.stats().forcing_rows > 0);
      else CHECK(p.stats().duplicate_rows > 0);
    }
  }
}

TEST(presolve_doubleton_differential_with_and_without_presolve) {
  std::mt19937 rng(14102026);
  Index continuous = 0, integral = 0;
  int optimal = 0, infeasible = 0;
  for (const bool integer_x : {false, true}) {
    for (int k = 0; k < 300; ++k) {
      const auto draw = [&](int lo, int hi) {
        return std::uniform_int_distribution<int>(lo, hi)(rng);
      };
      const double a = draw(0, 1) ? 1.0 : -1.0;
      const double b = draw(0, 1) ? 2.0 : -2.0;
      const double px = draw(-2, 2), py = draw(-2, 2), pz = draw(-1, 1);
      Model model = doubleton_model(integer_x, a, b, a * px + b * py);
      model.sense = draw(0, 1) ? samaya::ObjSense::kMinimize : samaya::ObjSense::kMaximize;
      for (Index j = 0; j < 3; ++j) model.obj[j] = draw(-5, 5) * 0.5;
      model.col_lower[0] += px;
      model.col_upper[0] += px;
      model.col_lower[1] += py;
      model.col_upper[1] += py;
      model.col_lower[2] += pz;
      model.col_upper[2] += pz;
      model.row_upper[1] = px + pz + draw(0, 2) + 0.5;
      model.row_upper[2] = py + pz + draw(0, 2) + 0.5;
      model.row_lower[3] = py - pz - draw(0, 2) - 0.5;
      if (k % 5 == 0) model.row_lower[3] += 20.0;
      const auto ref = samaya::test::reference_milp(model);
      const bool feasible = ref.status == samaya::test::ReferenceMilpResult::Status::kOptimal;
      REQUIRE(feasible || ref.status == samaya::test::ReferenceMilpResult::Status::kInfeasible);
      if (feasible) ++optimal;
      else ++infeasible;
      Presolve p(model, doubleton_options());
      const auto status = p.run();
      continuous += p.stats().doubleton_continuous_cols;
      integral += p.stats().doubleton_integer_cols;
      check_doubleton_reference(model, p, status);  // Bypasses solver fallback.
      for (const bool enabled : {false, true}) {
        samaya::Params params;
        params.log_level = 0;
        params.presolve = enabled;
        const auto result = samaya::Solver(params).solve(model);
        CHECK(result.status == (feasible ? samaya::Status::kOptimal : samaya::Status::kInfeasible));
        if (feasible) {
          CHECK(result.verified);
          CHECK_NEAR(result.objective, ref.objective, 1e-6 * (1.0 + std::fabs(ref.objective)));
        }
      }
    }
  }
  std::printf("  doubleton differential: %d optimal, %d infeasible; %d continuous, %d integer "
              "eliminations (600 MILPs, presolve on/off)\n", optimal, infeasible, continuous, integral);
  CHECK(continuous >= 150);
  CHECK(integral >= 150);
}


namespace {

Model substitution_model(double a = 1.0, double b = 2.0, double c = 2.0,
                         double rhs = 2.0, bool integer_x = false) {
  Model model;
  model.obj = {1.0, -2.0, -0.5};
  model.obj_offset = 1.25;
  model.col_lower = {-20.0, -2.0, -2.0};
  model.col_upper = {20.0, 2.0, 2.0};
  model.col_type.assign(3, samaya::VarType::kContinuous);
  if (integer_x) model.col_type.assign(3, samaya::VarType::kInteger);
  model.row_lower = {rhs, -kInf, -kInf, -kInf};
  model.row_upper = {rhs, 1.0, 1.0, 1.0};
  model.A = samaya::SparseMatrix::from_triplets(4, 3,
      {{0, 0, a}, {0, 1, b}, {0, 2, c},
       {1, 0, 1.0}, {1, 1, 1.0}, {1, 2, 1.0},
       {2, 0, -1.0}, {2, 1, 2.0}, {2, 2, -1.0},
       {3, 1, 1.0}, {3, 2, -1.0}});
  return model;
}

samaya::PresolveOptions substitution_options(bool mip = false, int passes = 50) {
  samaya::PresolveOptions options;
  options.mip = mip;
  options.max_passes = passes;
  return options;
}

void check_substitution_unchanged(const Model& model, const Presolve& p) {
  REQUIRE(p.reduced().num_rows() == model.num_rows());
  REQUIRE(p.reduced().num_cols() == model.num_cols());
  REQUIRE(p.reduced().A.nnz() == model.A.nnz());
  CHECK(std::equal(model.A.values().begin(), model.A.values().end(), p.reduced().A.values().begin()));
  CHECK(std::equal(model.A.row_index().begin(), model.A.row_index().end(),
                   p.reduced().A.row_index().begin()));
  CHECK(p.reduced().row_lower == model.row_lower);
  CHECK(p.reduced().row_upper == model.row_upper);
  CHECK(p.reduced().obj == model.obj);
  CHECK_EQ(p.reduced().obj_offset, model.obj_offset);
}

// Add an equation defining x from 2--5 bounded neighbours of a generated random core.
// x's finite bounds strictly contain the whole implied interval, for either pivot sign.
Model substitution_family(std::mt19937& rng, bool mip, int variant) {
  const auto family = static_cast<samaya::test::LpFamily>(variant % 5);
  Model core = samaya::test::random_lp(family, 5, 5, rng);
  const Index neighbours = std::uniform_int_distribution<Index>(2, 5)(rng);
  const Index n = std::max(core.num_cols(), neighbours);
  core.obj.resize(n, 1.0);
  core.col_lower.resize(n, -2.0);
  core.col_upper.resize(n, 2.0);
  core.col_type.resize(n, samaya::VarType::kContinuous);
  // Bound the tiny MILP reference search; retain the generated row families, including
  // infeasible and degenerate cores. LP loose families still include unbounded instances.
  for (Index j = 0; j < n; ++j) {
    if (mip || j < neighbours) {
      if (!std::isfinite(core.col_lower[j])) core.col_lower[j] = -3.0;
      if (!std::isfinite(core.col_upper[j])) core.col_upper[j] = core.col_lower[j] + 4.0;
    }
    if (mip && (variant % 3 == 0 || j % 2 == 0)) core.col_type[j] = samaya::VarType::kInteger;
  }
  Model model = core;
  const Index x = n;
  const double a = variant % 2 == 0 ? 2.0 : -2.0;
  const bool integer_x = mip && variant % 3 == 0;
  const double rhs = 2.0 * std::uniform_int_distribution<int>(-2, 2)(rng);
  std::vector<samaya::Triplet> entries;
  // Put the planted equation first, before other generated equations compete for x.
  entries.push_back({0, x, a});
  double min_act = 0.0, max_act = 0.0;
  for (Index j = 0; j < neighbours; ++j) {
    double b = 2.0 * std::uniform_int_distribution<int>(-3, 3)(rng);
    if (b == 0.0) b = 2.0;
    entries.push_back({0, j, b});
    min_act += b * (b > 0.0 ? model.col_lower[j] : model.col_upper[j]);
    max_act += b * (b > 0.0 ? model.col_upper[j] : model.col_lower[j]);
  }
  double lo = (rhs - max_act) / a, up = (rhs - min_act) / a;
  if (a < 0.0) std::swap(lo, up);
  model.obj.push_back(std::uniform_int_distribution<int>(-3, 3)(rng));
  model.col_lower.push_back(lo - 2.0);
  model.col_upper.push_back(up + 2.0);
  model.col_type.push_back(integer_x ? samaya::VarType::kInteger : samaya::VarType::kContinuous);
  model.row_lower.insert(model.row_lower.begin(), rhs);
  model.row_upper.insert(model.row_upper.begin(), rhs);
  const auto starts = core.A.col_start();
  const auto indices = core.A.row_index();
  const auto values = core.A.values();
  for (Index j = 0; j < core.A.cols(); ++j) {
    for (auto k = starts[j]; k < starts[j + 1]; ++k) entries.push_back({indices[k] + 1, j, values[k]});
  }
  // Two extra rows couple x to the core and make some optima depend on their duals.
  for (Index r = 0; r < 2; ++r) {
    const Index i = static_cast<Index>(model.row_lower.size());
    entries.insert(entries.end(), {{i, x, r == 0 ? 1.0 : -1.0}, {i, 0, 1.0},
                                  {i, 1, r == 0 ? -1.0 : 2.0}});
    model.row_lower.push_back(-kInf);
    model.row_upper.push_back(std::uniform_int_distribution<int>(1, 8)(rng));
  }
  model.A = samaya::SparseMatrix::from_triplets(static_cast<Index>(model.row_lower.size()), n + 1,
                                               std::move(entries));
  return model;
}

}  // namespace

TEST(presolve_substitution_primal_dual_both_signs_and_senses) {
  for (const double a : {1.0, -2.0}) {
    for (const auto sense : {samaya::ObjSense::kMinimize, samaya::ObjSense::kMaximize}) {
      Model model = substitution_model(a);
      model.sense = sense;
      Presolve p(model);
      REQUIRE(p.run() == PresolveStatus::kReduced);
      CHECK(p.stats().substitution_continuous_cols > 0);
      check_implied_postsolve(model, p);
      // Also exercise the snapshot while both other rows still exist.
      Presolve early(model, substitution_options(false, 1));
      REQUIRE(early.run() == PresolveStatus::kReduced);
      CHECK(early.stats().substitution_continuous_cols > 0);
      check_implied_postsolve(model, early);
    }
  }
}

TEST(presolve_substitution_continuous_mip) {
  Model model = substitution_model();
  model.col_type[1] = model.col_type[2] = samaya::VarType::kInteger;
  Presolve p(model, substitution_options(true));
  const auto status = p.run();
  REQUIRE(status == PresolveStatus::kReduced);
  CHECK(p.stats().substitution_continuous_cols > 0);
  check_doubleton_reference(model, p, status);
}

TEST(presolve_substitution_exact_zero_removes_entry) {
  Model model = substitution_model();  // x + 2y + 2z = 2.
  model.A = samaya::SparseMatrix::from_triplets(4, 3,
      {{0, 0, 1.0}, {0, 1, 2.0}, {0, 2, 2.0},
       {1, 0, 1.0}, {1, 1, 2.0}, {1, 2, 3.0},
       {2, 0, -1.0}, {2, 1, 2.0}, {2, 2, -1.0}, {3, 1, 1.0}, {3, 2, -1.0}});
  Presolve p(model, substitution_options(false, 1));
  REQUIRE(p.run() == PresolveStatus::kReduced);
  CHECK_EQ(p.stats().substitution_continuous_cols, 1);
  CHECK_EQ(p.reduced().A.nnz(), 5);
  const auto starts = p.reduced().A.col_start();
  const auto rows = p.reduced().A.row_index();
  // y's entry in original row 1 is exactly zero and must be absent.
  for (auto k = starts[0]; k < starts[1]; ++k) CHECK(p.row_map()[rows[k]] != 1);
  check_implied_postsolve(model, p);
}

TEST(presolve_substitution_exact_integer_divisibility) {
  for (const double a : {2.0, -2.0}) {
    for (const double b : {4.0, 3.0, 4.5}) {
      const Model model = substitution_model(a, b, 2.0, 2.0, true);
      Presolve p(model, substitution_options(true, 1));
      REQUIRE(p.run() == PresolveStatus::kReduced);
      if (b == 4.0) CHECK_EQ(p.stats().substitution_integer_cols, 1);
      else {
        CHECK_EQ(p.stats().substitution_integer_cols, 0);
        CHECK(p.stats().substitution_rejected_integer > 0);
      }
      check_doubleton_reference(model, p, PresolveStatus::kReduced);
    }
    const Model model = substitution_model(a, 4.0, 2.0, 1.0, true);
    Presolve p(model, substitution_options(true, 1));
    REQUIRE(p.run() == PresolveStatus::kReduced);
    CHECK_EQ(p.stats().substitution_integer_cols, 0);
    CHECK(p.stats().substitution_rejected_integer > 0);
    check_doubleton_reference(model, p, PresolveStatus::kReduced);
  }
}

TEST(presolve_substitution_nonredundant_and_tolerance_bounds) {
  for (const double bound : {1.0, 6.0, 6.0 + 1e-13}) {
    Model model = substitution_model();  // x's implied range is [-6, 10].
    model.col_lower[0] = -bound;
    Presolve p(model, substitution_options(false, 1));
    REQUIRE(p.run() == PresolveStatus::kReduced);
    CHECK_EQ(p.stats().substitution_continuous_cols, 0);
    if (bound == 1.0) CHECK(p.stats().substitution_rejected_not_implied_free > 0);
    else CHECK(p.stats().substitution_rejected_by_tolerance > 0);
    check_implied_postsolve(model, p);
  }
}

TEST(presolve_substitution_residual_rejected_atomically) {
  Model model = substitution_model();
  const double eps = std::ldexp(1.0, -41);
  model.A = samaya::SparseMatrix::from_triplets(4, 3,
      {{0, 0, 1.0}, {0, 1, 2.0}, {0, 2, 2.0},
       {1, 0, 1.0}, {1, 1, 2.0 + eps}, {1, 2, 3.0},
       {2, 0, -1.0}, {2, 1, 2.0}, {2, 2, -1.0}, {3, 1, 1.0}, {3, 2, -1.0}});
  Presolve p(model, substitution_options(false, 1));
  REQUIRE(p.run() == PresolveStatus::kReduced);
  CHECK_EQ(p.stats().substitution_continuous_cols, 0);
  CHECK_EQ(p.stats().substitution_rejected_residual, 1);
  check_substitution_unchanged(model, p);
}

TEST(presolve_substitution_fill_limit_and_length_caps) {
  // Five new neighbours in each of three other rows gives fill 15, excluding x itself.
  Model model;
  model.obj.assign(8, 1.0);
  model.obj[6] = -1.0;  // Prevent a dominated-column reduction from consuming this fixture.
  model.col_lower.assign(8, -1.0);
  model.col_upper.assign(8, 1.0);
  model.col_lower[0] = -100.0;
  model.col_upper[0] = 100.0;
  model.col_type.assign(8, samaya::VarType::kContinuous);
  model.row_lower = {0.0, -kInf, -kInf, -kInf};
  model.row_upper = {0.0, 1.0, 1.0, 1.0};
  std::vector<samaya::Triplet> entries = {{0, 0, 1.0}, {1, 0, 1.0}, {2, 0, -1.0},
      {3, 0, 2.0}, {1, 6, 1.0}, {1, 7, 1.0}, {2, 6, 1.0}, {2, 7, -1.0},
      {3, 6, 2.0}, {3, 7, -1.0}};
  for (Index j = 1; j <= 5; ++j) {
    entries.insert(entries.end(), {{0, j, 1.0}, {j + 3, j, 1.0}, {j + 3, 6, 1.0}});
    model.row_lower.push_back(-kInf);
    model.row_upper.push_back(0.5);
  }
  model.A = samaya::SparseMatrix::from_triplets(9, 8, entries);
  Presolve p(model, substitution_options(false, 1));
  REQUIRE(p.run() == PresolveStatus::kReduced);
  CHECK_EQ(p.stats().substitution_rejected_fill, 1);
  CHECK_EQ(p.stats().substitution_continuous_cols, 0);
  check_substitution_unchanged(model, p);

  // Reducing to two other rows is exactly the accepted fill limit of 10.
  entries.erase(std::remove_if(entries.begin(), entries.end(),
                              [](const samaya::Triplet& t) { return t.row == 3; }), entries.end());
  model.A = samaya::SparseMatrix::from_triplets(9, 8, entries);
  Presolve boundary(model);
  CHECK(samaya::PresolveSubstitutionTestAccess::substitutions(boundary));
  CHECK_EQ(boundary.stats().substitution_continuous_cols, 1);
  CHECK_EQ(boundary.stats().substitution_rejected_fill, 0);
  check_implied_postsolve(model, boundary);

  Model long_col = substitution_model();
  std::vector<samaya::Triplet> long_entries = {{0, 0, 1.0}, {0, 1, 2.0}, {0, 2, 2.0}};
  long_col.row_lower.assign(12, -kInf);
  long_col.row_upper.assign(12, 1.0);
  long_col.row_lower[0] = long_col.row_upper[0] = 2.0;
  for (Index i = 1; i < 12; ++i) long_entries.insert(long_entries.end(),
      {{i, 0, i % 2 == 0 ? 1.0 : -1.0}, {i, 1, 1.0}, {i, 2, 1.0}});
  long_col.A = samaya::SparseMatrix::from_triplets(12, 3, long_entries);
  Presolve q(long_col);
  CHECK(!samaya::PresolveSubstitutionTestAccess::substitutions(q));
  CHECK(q.stats().substitution_rejected_col_length > 0);
  check_substitution_unchanged(long_col, q);

  Model long_row;
  long_row.obj.assign(11, 1.0);
  long_row.col_lower.assign(11, -1.0);
  long_row.col_upper.assign(11, 1.0);
  long_row.col_type.assign(11, samaya::VarType::kContinuous);
  long_row.row_lower = long_row.row_upper = {0.0};
  std::vector<samaya::Triplet> row_entries;
  for (Index j = 0; j < 11; ++j) row_entries.push_back({0, j, 1.0});
  long_row.A = samaya::SparseMatrix::from_triplets(1, 11, std::move(row_entries));
  Presolve r(long_row);
  CHECK(!samaya::PresolveSubstitutionTestAccess::substitutions(r));
  CHECK_EQ(r.stats().substitution_rejected_row_length, 1);
  check_substitution_unchanged(long_row, r);
}

TEST(presolve_substitution_nonfinite_objective_rejected_atomically) {
  Model model = substitution_model();
  model.obj[0] = 1e308;
  Presolve p(model, substitution_options(false, 1));
  REQUIRE(p.run() == PresolveStatus::kReduced);
  CHECK_EQ(p.stats().substitution_continuous_cols, 0);
  CHECK_EQ(p.stats().substitution_rejected_nonfinite, 1);
  check_substitution_unchanged(model, p);
}

TEST(presolve_substitution_differential_with_and_without_presolve) {
  std::mt19937 rng(14102026);
  Index lp_substitutions = 0, mip_substitutions = 0, integer_substitutions = 0;
  int lp_instances = 0, mip_instances = 0;
  int optimal = 0, infeasible = 0, unbounded = 0;
  for (const bool mip : {false, true}) {
    for (int k = 0; k < 600; ++k) {
      const Model model = substitution_family(rng, mip, k);
      samaya::Status expected;
      double objective = 0.0;
      if (mip) {
        const auto ref = samaya::test::reference_milp(model);
        REQUIRE(ref.status == samaya::test::ReferenceMilpResult::Status::kOptimal ||
                ref.status == samaya::test::ReferenceMilpResult::Status::kInfeasible);
        expected = ref.status == samaya::test::ReferenceMilpResult::Status::kOptimal
            ? samaya::Status::kOptimal : samaya::Status::kInfeasible;
        objective = ref.objective;
      } else {
        const auto ref = ReferenceLp(model).solve();
        expected = ref.status == ReferenceResult::Status::kOptimal ? samaya::Status::kOptimal
            : ref.status == ReferenceResult::Status::kInfeasible ? samaya::Status::kInfeasible
                                                                : samaya::Status::kUnbounded;
        objective = ref.objective;
      }
      if (expected == samaya::Status::kOptimal) ++optimal;
      else if (expected == samaya::Status::kInfeasible) ++infeasible;
      else ++unbounded;
      Presolve p(model, substitution_options(mip));
      const auto status = p.run();
      const Index fired = p.stats().substitution_continuous_cols + p.stats().substitution_integer_cols;
      (mip ? mip_substitutions : lp_substitutions) += fired;
      integer_substitutions += p.stats().substitution_integer_cols;
      if (fired > 0) ++(mip ? mip_instances : lp_instances);
      if (mip) check_doubleton_reference(model, p, status);  // No solver fallback.
      else if (expected == samaya::Status::kOptimal) {
        REQUIRE(status == PresolveStatus::kReduced);
        check_implied_postsolve(model, p);  // Original-model dual verification.
      } else if (status == PresolveStatus::kInfeasible) {
        CHECK(expected == samaya::Status::kInfeasible);
      } else {
        const auto reduced_ref = ReferenceLp(p.reduced()).solve();
        CHECK(reduced_ref.status == (expected == samaya::Status::kInfeasible
            ? ReferenceResult::Status::kInfeasible : ReferenceResult::Status::kUnbounded));
      }
      for (const bool enabled : {false, true}) {
        samaya::Params params;
        params.log_level = 0;
        params.presolve = enabled;
        const auto result = samaya::Solver(params).solve(model);
        CHECK(result.status == expected);
        if (!mip || expected == samaya::Status::kOptimal) CHECK(result.verified);
        if (expected == samaya::Status::kOptimal) {
          CHECK_NEAR(result.objective, objective, 1e-6 * (1.0 + std::fabs(objective)));
          if (!mip) CHECK(samaya::verify_lp_optimality(model, result.col_value, result.row_dual).ok);
        }
      }
    }
  }
  std::printf("  substitution differential: %d optimal, %d infeasible, %d unbounded; "
              "%d LP, %d MILP (%d integer) eliminations (600 LPs + 600 MILPs, presolve on/off)\n",
              optimal, infeasible, unbounded, lp_substitutions, mip_substitutions,
              integer_substitutions);
  CHECK(lp_substitutions >= 150);
  CHECK(mip_substitutions >= 150);
  CHECK(integer_substitutions >= 30);
  CHECK(lp_instances >= 150);
  CHECK(mip_instances >= 150);
  std::printf("  substitution fired on %d/600 LPs and %d/600 MILPs\n", lp_instances, mip_instances);
}


TEST(presolve_substitution_doubletons_lp_and_rejected_mip) {
  for (const bool mip : {false, true}) {
    Model model;
    model.obj = {1.0, -2.0, -0.5};
    model.col_lower = {-2.0, -20.0, -2.0};
    model.col_upper = {2.0, 20.0, 2.0};
    model.col_type.assign(3, mip ? samaya::VarType::kInteger : samaya::VarType::kContinuous);
    model.row_lower = {0.0, -kInf, -kInf, -kInf};
    model.row_upper = {0.0, 1.0, 1.0, 1.0};
    model.A = samaya::SparseMatrix::from_triplets(4, 3,
        {{0, 0, 2.0}, {0, 1, 1.0}, {1, 0, 1.0}, {1, 2, 1.0},
         {2, 1, 1.0}, {2, 2, 1.0}, {3, 1, -1.0}, {3, 2, 1.0}});
    Presolve p(model, substitution_options(mip, 1));
    REQUIRE(p.run() == PresolveStatus::kReduced);
    if (mip) {
      CHECK_EQ(p.stats().doubleton_integer_cols, 0);
      CHECK_EQ(p.stats().doubleton_rejected_divisibility, 1);
      CHECK_EQ(p.stats().substitution_integer_cols, 1);
      check_doubleton_reference(model, p, PresolveStatus::kReduced);
    } else {
      CHECK_EQ(p.stats().substitution_continuous_cols, 1);
      check_implied_postsolve(model, p);
    }
  }
}

TEST(presolve_substitution_pivot_and_markowitz_choice) {
  Model tiny = substitution_model(0.001);
  tiny.col_lower[0] = -20000.0;
  tiny.col_upper[0] = 20000.0;
  Presolve rejected(tiny);
  CHECK(!samaya::PresolveSubstitutionTestAccess::substitutions(rejected));
  CHECK(rejected.stats().substitution_rejected_pivot > 0);
  check_substitution_unchanged(tiny, rejected);
  for (const bool shorter : {false, true}) {
    Model model = substitution_model();
    model.obj.assign(3, 0.0);
    model.col_lower[0] = model.col_lower[1] = -kInf;
    model.col_upper[0] = model.col_upper[1] = kInf;
    model.A = samaya::SparseMatrix::from_triplets(4, 3,
        {{0, 0, 1.0}, {0, 1, 2.0}, {0, 2, 1.0},
         {1, 0, 1.0}, {1, 1, 1.0}, {1, 2, 1.0}, {2, shorter ? 0 : 1, 1.0}, {2, 2, -1.0},
         {3, shorter ? 0 : 1, 1.0}, {3, 2, -1.0}});
    Presolve p(model);
    CHECK(samaya::PresolveSubstitutionTestAccess::substitutions(p));
    CHECK_EQ(p.stats().substitution_continuous_cols, 1);
    const Index removed = shorter ? 1 : 0;
    CHECK(std::find(p.col_map().begin(), p.col_map().end(), removed) == p.col_map().end());
    check_implied_postsolve(model, p);
  }
  // Equal column lengths choose the larger coefficient.
  Model tied = substitution_model();
  tied.obj.assign(3, 0.0);
  tied.col_lower[0] = tied.col_lower[1] = -kInf;
  tied.col_upper[0] = tied.col_upper[1] = kInf;
  tied.A = samaya::SparseMatrix::from_triplets(3, 3,
      {{0, 0, 1.0}, {0, 1, 2.0}, {0, 2, 1.0},
       {1, 0, 1.0}, {1, 2, 1.0}, {2, 1, 1.0}, {2, 2, -1.0}});
  tied.row_lower.resize(3);
  tied.row_upper.resize(3);
  Presolve p(tied);
  CHECK(samaya::PresolveSubstitutionTestAccess::substitutions(p));
  CHECK(std::find(p.col_map().begin(), p.col_map().end(), 1) == p.col_map().end());
  check_implied_postsolve(tied, p);
}

// 3x + y = 0 and 0.9x + 0.30000000000000004y in [-1, 1]: substituting x leaves y the coefficient
// 0.30000000000000004 - 0.9/3, a few ulps of rounding noise. A row left with only such noise
// turns into a wrong bound (random LP #866 of presolve_postsolve_restores_verified_optimal_
// solutions was declared infeasible that way), so both C2 and C3 must skip the substitution.
Model cancellation_model() {
  Model model;
  model.obj = {0.0, 1.0};
  model.col_lower = {-kInf, -1.0};
  model.col_upper = {kInf, 1.0};
  model.col_type.assign(2, samaya::VarType::kContinuous);
  model.row_lower = {0.0, -1.0};
  model.row_upper = {0.0, 1.0};
  model.A = samaya::SparseMatrix::from_triplets(2, 2,
      {{0, 0, 3.0}, {0, 1, 1.0}, {1, 0, 0.9}, {1, 1, 0.1 + 0.2}});
  return model;
}

TEST(presolve_substitution_rejects_cancelled_coefficient) {
  const Model model = cancellation_model();
  REQUIRE(0.1 + 0.2 - 0.9 * 1.0 / 3.0 != 0.0);  // The noise this test is about.
  Presolve p(model);
  REQUIRE(p.run() == PresolveStatus::kReduced);
  CHECK_EQ(p.stats().substitution_continuous_cols, 0);
  CHECK(p.stats().substitution_rejected_residual > 0);
  check_implied_postsolve(model, p);
}

TEST(presolve_doubleton_rejects_cancelled_coefficient) {
  const Model model = cancellation_model();
  Presolve p(model, doubleton_options());
  const PresolveStatus status = p.run();
  REQUIRE(status == PresolveStatus::kReduced);
  CHECK_EQ(p.stats().doubleton_continuous_cols, 0);
  CHECK(p.stats().doubleton_rejected_residual > 0);
  check_doubleton_reference(model, p, status);
}