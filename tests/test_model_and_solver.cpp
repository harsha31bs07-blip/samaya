#include <cmath>
#include <sstream>
#include <string>

#include "samaya.hpp"
#include "samaya_c.h"
#include "test_framework.hpp"

using samaya::kInf;
using samaya::Model;

namespace {

Model two_by_two() {
  Model m;
  m.obj = {1.0, 1.0};
  m.col_lower = {0.0, 0.0};
  m.col_upper = {kInf, 1.0};
  m.col_type = {samaya::VarType::kContinuous, samaya::VarType::kInteger};
  m.row_lower = {1.0, -kInf};
  m.row_upper = {1.0, 1e6};
  m.A = samaya::SparseMatrix::from_triplets(2, 2, {{0, 0, 1e-3}, {0, 1, 1.0}, {1, 1, 1e4}});
  return m;
}

}  // namespace

TEST(model_validate_accepts_consistent_model) {
  CHECK(two_by_two().validate().empty());
}

TEST(model_validate_rejects_bad_models) {
  Model m = two_by_two();
  m.obj.pop_back();
  CHECK(!m.validate().empty());

  m = two_by_two();
  m.col_lower[0] = 2.0;
  m.col_upper[0] = 1.0;
  CHECK(m.validate().empty());  // Crossed bounds are infeasible, not invalid.
  CHECK(m.crossed_bounds().find("column 0") != std::string::npos);

  m = two_by_two();
  m.row_lower[1] = kInf;
  CHECK(!m.validate().empty());

  m = two_by_two();
  m.obj[0] = std::nan("");
  CHECK(!m.validate().empty());

  m = two_by_two();
  m.Q = samaya::SparseMatrix::from_triplets(2, 2, {{0, 1, 1.0}});  // Upper triangle.
  CHECK(!m.validate().empty());
}

TEST(model_stats_ranges) {
  const Model m = two_by_two();
  const samaya::ModelStats s = samaya::compute_stats(m);
  CHECK_EQ(s.rows, 2);
  CHECK_EQ(s.cols, 2);
  CHECK_EQ(s.nnz, 3);
  CHECK_EQ(s.integers, 1);
  CHECK_EQ(s.binaries, 1);
  CHECK_EQ(s.equality_rows, 1);
  CHECK_EQ(s.matrix.min_abs, 1e-3);
  CHECK_EQ(s.matrix.max_abs, 1e4);
  CHECK_EQ(s.rhs.max_abs, 1e6);

  std::ostringstream out;
  samaya::print_stats(out, m, s);
  CHECK(out.str().find("MILP") != std::string::npos);
}

TEST(solver_reports_invalid_and_unimplemented) {
  samaya::Params params;
  params.log_level = 0;
  const samaya::Solver solver(params);

  Model bad = two_by_two();
  bad.col_upper.pop_back();
  CHECK(solver.solve(bad).status == samaya::Status::kInvalidModel);

  // A convex QP is solved by the interior point and verified; MIQP is not implemented yet.
  Model qp = two_by_two();
  qp.col_type.assign(qp.col_type.size(), samaya::VarType::kContinuous);
  qp.Q = samaya::SparseMatrix::from_triplets(2, 2, {{0, 0, 1.0}});
  const samaya::Result r = solver.solve(qp);
  CHECK(r.status == samaya::Status::kOptimal);
  CHECK(r.verified);
  Model miqp = two_by_two();
  miqp.Q = samaya::SparseMatrix::from_triplets(2, 2, {{0, 0, 1.0}});
  const samaya::Result r2 = solver.solve(miqp);
  CHECK(r2.status == samaya::Status::kNotImplemented);
  CHECK(std::isnan(r2.objective));

  // The mixed-integer model itself is solved: x1 = 1, x0 = 0 beats x1 = 0, x0 = 1000.
  const samaya::Result mip = solver.solve(two_by_two());
  CHECK(mip.status == samaya::Status::kOptimal);
  CHECK(mip.verified);
  CHECK_NEAR(mip.objective, 1.0, 1e-9);
}

TEST(c_api_round_trip) {
  samaya_model* model = nullptr;
  const std::string path = std::string(SAMAYA_TEST_DATA_DIR) + "/tiny_lp.mps";
  REQUIRE(samaya_read_mps(path.c_str(), &model, nullptr, 0) == 0);
  CHECK_EQ(samaya_model_num_rows(model), 2);
  CHECK_EQ(samaya_model_num_cols(model), 2);
  CHECK_EQ(samaya_model_num_nonzeros(model), 4LL);
  double objective = 0.0;
  CHECK_EQ(samaya_solve(model, &objective), static_cast<int>(SAMAYA_OPTIMAL));
  CHECK_NEAR(objective, 11.0, 1e-9);
  samaya_model_free(model);

  char err[128];
  samaya_model* missing = nullptr;
  CHECK(samaya_read_mps("/nonexistent.mps", &missing, err, sizeof err) != 0);
  CHECK(missing == nullptr);
  CHECK(std::string(err).find("cannot open") != std::string::npos);
  CHECK_EQ(std::string(samaya_status_string(SAMAYA_OPTIMAL)), "optimal");
  CHECK_EQ(std::string(samaya_version()), std::string(samaya::version()));
}

TEST(solver_solves_and_verifies_lp_outcomes) {
  samaya::Params params;
  params.log_level = 0;
  const samaya::Solver solver(params);

  const Model lp = samaya::read_mps(std::string(SAMAYA_TEST_DATA_DIR) + "/tiny_lp.mps");
  samaya::Result r = solver.solve(lp);
  REQUIRE(r.status == samaya::Status::kOptimal);
  CHECK(r.verified);
  CHECK_NEAR(r.objective, 11.0, 1e-9);
  CHECK_EQ(r.dual_bound, r.objective);
  CHECK_EQ(r.col_value.size(), 2u);
  CHECK_EQ(r.row_dual.size(), 2u);
  CHECK(r.simplex_iterations > 0);
  CHECK(r.max_primal_violation <= 1e-9);

  Model crossed = lp;
  crossed.row_lower[0] = 5.0;  // 5 <= x + y <= 4.
  r = solver.solve(crossed);
  CHECK(r.status == samaya::Status::kInfeasible);
  CHECK(r.verified);
  CHECK(r.message.find("row 0") != std::string::npos);

  Model infeasible = lp;
  infeasible.row_upper[1] = -1.0;  // x + 3y <= -1 with x, y >= 0.
  r = solver.solve(infeasible);
  CHECK(r.status == samaya::Status::kInfeasible);
  CHECK(r.verified);
  CHECK(samaya::verify_infeasibility(infeasible, r.infeasibility_certificate).ok);

  Model unbounded = lp;
  unbounded.row_upper = {kInf, kInf};
  unbounded.col_upper[0] = kInf;
  r = solver.solve(unbounded);
  CHECK(r.status == samaya::Status::kUnbounded);
  CHECK(r.verified);
  CHECK(samaya::verify_unbounded_ray(unbounded, r.unbounded_ray).ok);
  CHECK(std::isnan(r.objective));

  params.verify = false;
  r = samaya::Solver(params).solve(lp);
  CHECK(r.status == samaya::Status::kOptimal);
  CHECK(!r.verified);
}
