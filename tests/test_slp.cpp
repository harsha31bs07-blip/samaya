// SLP tests (docs/SLP_PLAN.md section 5): Haverly's pooling problems from near-global and trap
// starts, a stationary point that is not a local optimum, the starting-point projection, linear
// infeasibility, and a maximization.

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "samaya/io.hpp"
#include "samaya/solver.hpp"
#include "samaya/verify.hpp"
#include "slp/slp.hpp"
#include "test_framework.hpp"

using samaya::kInf;
using samaya::Model;
using samaya::Result;
using samaya::Status;

namespace {

constexpr double kTol = 1e-6;

bool near(double a, double b) { return std::fabs(a - b) <= kTol * (1.0 + std::fabs(b)); }

// Haverly's pooling problem in the p-formulation; columns (a, b, u, v, z1, z2, p): a, b into
// the pool, u, v from the pool to products 1 and 2, z1, z2 component C directly, p the pool's
// sulfur quality.
//   min 6a + cB b + 10(z1 + z2) - 9(u + z1) - 15(v + z2)
//   a + b - u - v = 0;  3a + b - p u - p v = 0;
//   p u + 2 z1 - 2.5 u - 2.5 z1 <= 0;  p v + 2 z2 - 1.5 v - 1.5 z2 <= 0;
//   u + z1 <= D1;  v + z2 <= 200;  flows >= 0;  1 <= p <= 3.
Model haverly(double d1, double cb) {
  enum { kA, kB, kU, kV, kZ1, kZ2, kP };
  Model m;
  m.obj = {6, cb, -9, -15, 10 - 9, 10 - 15, 0};
  m.col_lower = {0, 0, 0, 0, 0, 0, 1};
  m.col_upper = {kInf, kInf, kInf, kInf, kInf, kInf, 3};
  m.col_type.assign(7, samaya::VarType::kContinuous);
  m.row_lower = {0, 0, -kInf, -kInf, -kInf, -kInf};
  m.row_upper = {0, 0, 0, 0, d1, 200};
  m.A = samaya::SparseMatrix::from_triplets(
      6, 7,
      {{0, kA, 1}, {0, kB, 1}, {0, kU, -1}, {0, kV, -1},
       {1, kA, 3}, {1, kB, 1},
       {2, kU, -2.5}, {2, kZ1, 2 - 2.5},
       {3, kV, -1.5}, {3, kZ2, 2 - 1.5},
       {4, kU, 1}, {4, kZ1, 1},
       {5, kV, 1}, {5, kZ2, 1}});
  // col1 <= col2: p (6) after u (2) and v (3).
  m.row_quadratic = {{1, kU, kP, -1.0}, {1, kV, kP, -1.0}, {2, kU, kP, 1.0}, {3, kV, kP, 1.0}};
  return m;
}

Result solve_from(const Model& m, const std::vector<double>& start) {
  samaya::Params params;
  params.log_level = 0;
  params.mip_start = start;
  return samaya::Solver(params).solve(m);
}

struct Case {
  const char* name;
  double d1, cb;
  double global;
  std::vector<double> near_global;  // A feasible start next to a global solution.
  std::vector<double> trap;          // The reviewer's trap start ...
  std::vector<double> local;         // ... and the stationary point it leads to.
  double local_value;
};

std::vector<Case> haverly_cases() {
  return {
      {"Haverly 1", 100, 16, -400, {0, 99, 0, 99, 0, 99, 1}, {49, 0, 49, 0, 49, 0, 3},
       {50, 0, 50, 0, 50, 0, 3}, -100},
      {"Haverly 2", 600, 16, -600, {297, 0, 297, 0, 297, 0, 3}, {0, 99, 0, 99, 0, 99, 1},
       {0, 100, 0, 100, 0, 100, 1}, -400},
      {"Haverly 3", 100, 13, -750, {49.5, 148.5, 0, 198, 0, 0, 1.5},
       {74.25, 24.75, 99, 0, 0, 0, 2.5}, {75, 25, 100, 0, 0, 0, 2.5}, -125},
  };
}

}  // namespace

TEST(slp_haverly_global_from_near_global_starts) {
  for (const Case& c : haverly_cases()) {
    const Model m = haverly(c.d1, c.cb);
    const Result r = solve_from(m, c.near_global);
    std::printf("  %s near global: %s, objective %.9g (global %g), %lld SLP iterations\n",
                c.name, samaya::to_string(r.status), r.objective, c.global, r.slp_iterations);
    CHECK(r.status == Status::kFirstOrderStationary);
    CHECK(r.verified);
    CHECK(near(r.objective, c.global));
    CHECK(samaya::verify_primal(m, r.col_value).ok);
    CHECK(samaya::verify_qcp_kkt(m, r.col_value, r.row_dual).ok);
  }
}

TEST(slp_haverly_trap_starts_end_at_stationary_local_points) {
  for (const Case& c : haverly_cases()) {
    const Model m = haverly(c.d1, c.cb);
    const Result r = solve_from(m, c.trap);
    std::printf("  %s trap: %s, objective %.9g (local %g)\n", c.name,
                samaya::to_string(r.status), r.objective, c.local_value);
    CHECK(r.status == Status::kFirstOrderStationary);  // Stationary, never "optimal".
    CHECK(r.verified);
    CHECK(near(r.objective, c.local_value));
    for (std::size_t j = 0; j < c.local.size(); ++j) CHECK(near(r.col_value[j], c.local[j]));
  }
}

TEST(slp_haverly_from_the_default_start) {
  // No start: the columns in quadratic terms begin at the middle of their finite bounds (p = 2)
  // or at a bound, then the projection and the SLP. Whatever stationary point it reaches must be
  // verified, and no better than the global optimum.
  for (const Case& c : haverly_cases()) {
    const Model m = haverly(c.d1, c.cb);
    const Result r = solve_from(m, {});
    std::printf("  %s default start: %s, objective %.9g, %lld SLP iterations\n", c.name,
                samaya::to_string(r.status), r.objective, r.slp_iterations);
    CHECK(r.status == Status::kFirstOrderStationary);
    CHECK(r.verified);
    CHECK(r.objective >= c.global - kTol * (1.0 + std::fabs(c.global)));
  }
}

TEST(slp_haverly_from_mps_files) {
  // The same three problems read from tests/instances/haverly*.mps (QCMATRIX sections, written
  // by cases/pooling.py --haverly): the reader's terms must equal the hand-built model's, and
  // the solves must agree.
  const char* files[] = {"haverly1.mps", "haverly2.mps", "haverly3.mps"};
  const std::vector<Case> cases = haverly_cases();
  for (std::size_t k = 0; k < cases.size(); ++k) {
    const Model read = samaya::read_mps(std::string(SAMAYA_TEST_DATA_DIR) + "/" + files[k]);
    const Model built = haverly(cases[k].d1, cases[k].cb);
    REQUIRE(read.problem_class() == samaya::ProblemClass::kQCP);
    REQUIRE(read.num_cols() == built.num_cols());
    // Both list each (row, col1, col2) once; compare the activities they define at a point.
    const std::vector<double> x{3, 5, 7, 2, 4, 6, 1.7};
    double act_read = 0.0, act_built = 0.0;
    for (const auto& t : read.row_quadratic) act_read += (t.row + 1) * t.value * x[t.col1] * x[t.col2];
    for (const auto& t : built.row_quadratic) act_built += (t.row + 1) * t.value * x[t.col1] * x[t.col2];
    CHECK(near(act_read, act_built));
    const Result from_file = solve_from(read, {});
    const Result from_code = solve_from(built, {});
    std::printf("  %s: %s, objective %.9g (hand-built %.9g)\n", files[k],
                samaya::to_string(from_file.status), from_file.objective, from_code.objective);
    CHECK(from_file.status == Status::kFirstOrderStationary);
    CHECK(from_file.verified);
    CHECK(near(from_file.objective, from_code.objective));
  }
}

TEST(slp_linearization_has_trust_bounds_and_the_right_rows) {
  // Haverly 1 linearized at the trap start (a, b, u, v, z1, z2, p) = (49, 0, 49, 0, 49, 0, 3),
  // radius 0.1. Columns in quadratic terms are u, v (no finite range: scale max(1, |x|)) and p
  // (range 2): u in [44.1, 53.9], v in [0, 0.1], p in [2.8, 3]; the others keep their bounds.
  const Model m = haverly(100, 16);
  const std::vector<double> c{49, 0, 49, 0, 49, 0, 3};
  const Model lp = samaya::slp_linearization_for_test(m, c, 0.1);
  REQUIRE(lp.num_cols() == 7 + 2 * 3);  // Three quadratic rows, two elastics each.
  CHECK(near(lp.col_lower[2], 44.1));
  CHECK(near(lp.col_upper[2], 53.9));
  CHECK(near(lp.col_lower[3], 0.0));
  CHECK(near(lp.col_upper[3], 0.1));
  CHECK(near(lp.col_lower[6], 2.8));
  CHECK(near(lp.col_upper[6], 3.0));
  for (const int j : {0, 1, 4, 5}) {
    CHECK(lp.col_lower[j] == m.col_lower[j]);
    CHECK(lp.col_upper[j] == m.col_upper[j]);
  }
  // pool_quality, 3a + b - p u - p v = 0, linearized: 3a + b - 3u - 3v - 49p = -147.
  std::vector<double> row(7, 0.0);
  const auto start = lp.A.col_start();
  const auto index = lp.A.row_index();
  const auto value = lp.A.values();
  for (int j = 0; j < 7; ++j) {
    for (auto p = start[j]; p < start[j + 1]; ++p) {
      if (index[p] == 1) row[j] = value[p];
    }
  }
  CHECK(near(row[0], 3) && near(row[1], 1) && near(row[2], -3) && near(row[3], -3));
  CHECK(near(row[6], -49));
  CHECK(near(lp.row_lower[1], -147) && near(lp.row_upper[1], -147));
}

TEST(slp_stationary_is_not_called_optimal) {
  // min z s.t. z - u v = 0, -1 <= u, v <= 1, from the origin: a KKT point where the LP sees no
  // progress, though (t, -t, -t^2) is better nearby. The status says stationary.
  Model m;
  m.obj = {0, 0, 1};
  m.col_lower = {-1, -1, -10};
  m.col_upper = {1, 1, 10};
  m.col_type.assign(3, samaya::VarType::kContinuous);
  m.row_lower = {0};
  m.row_upper = {0};
  m.A = samaya::SparseMatrix::from_triplets(1, 3, {{0, 2, 1}});
  m.row_quadratic = {{0, 0, 1, -1.0}};
  const Result r = solve_from(m, {0, 0, 0});
  CHECK(r.status == Status::kFirstOrderStationary);
  CHECK(r.status != Status::kOptimal);
  CHECK(r.message.find("not certified") != std::string::npos);
}

TEST(slp_projection_repairs_a_start_that_breaks_linear_rows) {
  // u + v = 1 with u, v in [0, 10] and u v >= 0.16; the default start (5, 5) breaks the linear
  // row, so the first LP would be infeasible without the projection. max u: the KKT point is
  // u = 0.8, v = 0.2 (u v = 0.16 with u + v = 1).
  Model m;
  m.sense = samaya::ObjSense::kMaximize;
  m.obj = {1, 0};
  m.col_lower = {0, 0};
  m.col_upper = {10, 10};
  m.col_type.assign(2, samaya::VarType::kContinuous);
  m.row_lower = {1, 0.16};
  m.row_upper = {1, kInf};
  m.A = samaya::SparseMatrix::from_triplets(2, 2, {{0, 0, 1}, {0, 1, 1}});
  m.row_quadratic = {{1, 0, 1, 1.0}};
  samaya::Params params;
  params.log_level = 0;
  const Result r = samaya::Solver(params).solve(m);
  std::printf("  projection: %s, u = %.9g, v = %.9g\n", samaya::to_string(r.status),
              r.col_value.empty() ? NAN : r.col_value[0], r.col_value.empty() ? NAN : r.col_value[1]);
  CHECK(r.status == Status::kFirstOrderStationary);
  CHECK(r.verified);
  CHECK(near(r.objective, 0.8));
}

TEST(slp_linear_infeasibility_is_proven) {
  // x + y >= 5 with x, y <= 1 cannot hold whatever the quadratic row says.
  Model m;
  m.obj = {1, 1};
  m.col_lower = {0, 0};
  m.col_upper = {1, 1};
  m.col_type.assign(2, samaya::VarType::kContinuous);
  m.row_lower = {5, 0.1};
  m.row_upper = {kInf, kInf};
  m.A = samaya::SparseMatrix::from_triplets(2, 2, {{0, 0, 1}, {0, 1, 1}});
  m.row_quadratic = {{1, 0, 1, 1.0}};
  const Result r = solve_from(m, {});
  CHECK(r.status == Status::kInfeasible);
  CHECK(r.verified);
  CHECK(samaya::verify_infeasibility(m, r.infeasibility_certificate).ok);
}

TEST(slp_haverly_as_a_maximization) {
  // Profit maximization: the objective negated, max profit 400 from the near-global start.
  Model m = haverly(100, 16);
  m.sense = samaya::ObjSense::kMaximize;
  for (double& c : m.obj) c = -c;
  const Result r = solve_from(m, {0, 99, 0, 99, 0, 99, 1});
  CHECK(r.status == Status::kFirstOrderStationary);
  CHECK(r.verified);
  CHECK(near(r.objective, 400.0));
  CHECK(samaya::verify_qcp_kkt(m, r.col_value, r.row_dual).ok);
}
