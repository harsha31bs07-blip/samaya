#include <vector>

#include "samaya/verify.hpp"
#include "test_framework.hpp"

using samaya::kInf;
using samaya::Model;
using samaya::VerifyReport;

namespace {

// max 3x + 2y  s.t.  x + y <= 4,  x + 3y <= 6,  0 <= x <= 3,  y >= 0.  Optimum (3, 1).
Model tiny_lp() {
  Model m;
  m.sense = samaya::ObjSense::kMaximize;
  m.obj = {3, 2};
  m.col_lower = {0, 0};
  m.col_upper = {3, kInf};
  m.col_type = {samaya::VarType::kContinuous, samaya::VarType::kContinuous};
  m.row_lower = {-kInf, -kInf};
  m.row_upper = {4, 6};
  m.A = samaya::SparseMatrix::from_triplets(2, 2, {{0, 0, 1}, {0, 1, 1}, {1, 0, 1}, {1, 1, 3}});
  return m;
}

}  // namespace

TEST(verify_primal_feasible_point_and_objective) {
  const Model m = tiny_lp();
  const VerifyReport r = samaya::verify_primal(m, std::vector<double>{3, 1});
  CHECK(r.ok);
  CHECK_EQ(r.objective, 11.0);
  CHECK_EQ(r.max_row_violation, 0.0);
}

TEST(verify_primal_detects_violations) {
  const Model m = tiny_lp();
  CHECK(!samaya::verify_primal(m, std::vector<double>{3.1, 0}).ok);   // Column bound.
  CHECK(!samaya::verify_primal(m, std::vector<double>{3, 1.5}).ok);   // Row.
  CHECK(!samaya::verify_primal(m, std::vector<double>{3}).ok);        // Length.
  CHECK(samaya::verify_primal(m, std::vector<double>{3, 1 + 1e-9}).ok);  // Within tolerance.

  Model mip = m;
  mip.col_type[1] = samaya::VarType::kInteger;
  CHECK(samaya::verify_primal(mip, std::vector<double>{3, 1}).ok);
  const VerifyReport frac = samaya::verify_primal(mip, std::vector<double>{2.5, 0.5});
  CHECK(!frac.ok);
  CHECK_EQ(frac.max_integrality_violation, 0.5);
}

TEST(verify_primal_quadratic_objective) {
  Model m = tiny_lp();
  // Lower triangle of Q = [[2, 1], [1, 4]]: 1/2 x'Qx = x^2 + x y + 2 y^2.
  m.Q = samaya::SparseMatrix::from_triplets(2, 2, {{0, 0, 2}, {1, 0, 1}, {1, 1, 4}});
  const VerifyReport r = samaya::verify_primal(m, std::vector<double>{3, 1});
  CHECK_EQ(r.objective, 11.0 + 9.0 + 3.0 + 2.0);
}

TEST(verify_lp_optimality_accepts_valid_duals_only) {
  const Model m = tiny_lp();
  const std::vector<double> x{3, 1};
  // The optimum is degenerate; y = (2, 0) and y = (0, 2/3) are both valid.
  CHECK(samaya::verify_lp_optimality(m, x, std::vector<double>{2, 0}).ok);
  CHECK(samaya::verify_lp_optimality(m, x, std::vector<double>{0, 2.0 / 3.0}).ok);
  // Wrong sign on a row dual.
  CHECK(!samaya::verify_lp_optimality(m, x, std::vector<double>{-1, 0}).ok);
  // Feasible but not optimal point with its natural duals.
  const VerifyReport sub = samaya::verify_lp_optimality(m, std::vector<double>{0, 0},
                                                        std::vector<double>{0, 0});
  CHECK(!sub.ok);
  CHECK(sub.max_dual_violation > 0.5);

  Model min_model = m;
  min_model.sense = samaya::ObjSense::kMinimize;
  // For "min 3x + 2y" the optimum is (0, 0) with zero duals; (3, 1) is not optimal.
  CHECK(samaya::verify_lp_optimality(min_model, std::vector<double>{0, 0},
                                     std::vector<double>{0, 0})
            .ok);
  CHECK(!samaya::verify_lp_optimality(min_model, x, std::vector<double>{2, 0}).ok);
}

TEST(verify_infeasibility_certificate) {
  // x + y >= 4 and x + y <= 2 with x, y >= 0.
  Model m;
  m.obj = {1, 1};
  m.col_lower = {0, 0};
  m.col_upper = {kInf, kInf};
  m.col_type.assign(2, samaya::VarType::kContinuous);
  m.row_lower = {4, -kInf};
  m.row_upper = {kInf, 2};
  m.A = samaya::SparseMatrix::from_triplets(2, 2, {{0, 0, 1}, {0, 1, 1}, {1, 0, 1}, {1, 1, 1}});
  CHECK(samaya::verify_infeasibility(m, std::vector<double>{1, -1}).ok);
  CHECK(samaya::verify_infeasibility(m, std::vector<double>{-3, 3}).ok);
  CHECK(!samaya::verify_infeasibility(m, std::vector<double>{1, 0}).ok);
  CHECK(!samaya::verify_infeasibility(m, std::vector<double>{0, 0}).ok);
  m.row_lower[0] = 2;  // Now feasible: no certificate can exist.
  CHECK(!samaya::verify_infeasibility(m, std::vector<double>{1, -1}).ok);
}

TEST(verify_unbounded_ray) {
  // min -x - y  s.t.  x - y <= 1,  x, y >= 0.
  Model m;
  m.obj = {-1, -1};
  m.col_lower = {0, 0};
  m.col_upper = {kInf, kInf};
  m.col_type.assign(2, samaya::VarType::kContinuous);
  m.row_lower = {-kInf};
  m.row_upper = {1};
  m.A = samaya::SparseMatrix::from_triplets(1, 2, {{0, 0, 1}, {0, 1, -1}});
  CHECK(samaya::verify_unbounded_ray(m, std::vector<double>{1, 1}).ok);
  CHECK(samaya::verify_unbounded_ray(m, std::vector<double>{0, 2}).ok);
  CHECK(!samaya::verify_unbounded_ray(m, std::vector<double>{1, 0}).ok);   // Leaves the row.
  CHECK(!samaya::verify_unbounded_ray(m, std::vector<double>{-1, 0}).ok);  // Leaves a bound.
  m.sense = samaya::ObjSense::kMaximize;
  CHECK(!samaya::verify_unbounded_ray(m, std::vector<double>{1, 1}).ok);  // Not improving.
}

namespace {

Model continuous(int n) {
  Model m;
  m.col_lower.assign(static_cast<std::size_t>(n), 0.0);
  m.col_upper.assign(static_cast<std::size_t>(n), 10.0);
  m.col_type.assign(static_cast<std::size_t>(n), samaya::VarType::kContinuous);
  return m;
}

// min -x - y  s.t.  x y <= 4 (a bilinear row with no linear part), 0 <= x, y <= 10.
// KKT point (2, 2): the gradient (-1, -1) equals lambda (y, x) = lambda (2, 2) with
// lambda = -1/2, nonpositive as a row at its upper bound requires.
Model bilinear_cap() {
  Model m = continuous(2);
  m.obj = {-1, -1};
  m.row_lower = {-kInf};
  m.row_upper = {4};
  m.A = samaya::SparseMatrix(1, 2);
  m.row_quadratic = {{0, 0, 1, 1.0}};
  return m;
}

}  // namespace

TEST(verify_primal_includes_quadratic_rows) {
  Model m = continuous(2);
  m.obj = {0, 0};
  // 2u^2 + 6uv + 5v^2 + u (QCMATRIX diagonal 2 and 5, symmetric off-diagonal 3) <= 91.
  m.row_lower = {-kInf};
  m.row_upper = {91};
  m.A = samaya::SparseMatrix::from_triplets(1, 2, {{0, 0, 1}});
  m.row_quadratic = {{0, 0, 0, 2.0}, {0, 0, 1, 6.0}, {0, 1, 1, 5.0}};
  CHECK(m.validate().empty());
  CHECK(m.problem_class() == samaya::ProblemClass::kQCP);
  CHECK(samaya::verify_primal(m, std::vector<double>{2, 3}).ok);   // 89 + 2 = 91.
  m.row_upper = {90};
  const VerifyReport r = samaya::verify_primal(m, std::vector<double>{2, 3});
  CHECK(!r.ok);  // Violated only through the quadratic terms.
  CHECK(r.max_row_violation > 0.0);
}

TEST(verify_qcp_kkt_accepts_the_multipliers_of_a_kkt_point_only) {
  const Model m = bilinear_cap();
  const std::vector<double> x{2, 2};
  CHECK(samaya::verify_qcp_kkt(m, x, std::vector<double>{-0.5}).ok);
  CHECK(!samaya::verify_qcp_kkt(m, x, std::vector<double>{-0.4}).ok);  // Wrong multiplier.
  CHECK(!samaya::verify_qcp_kkt(m, x, std::vector<double>{0.5}).ok);   // Wrong sign.
  // Feasible but not stationary: at (1, 3) the row is not active, so no multiplier can help.
  CHECK(!samaya::verify_qcp_kkt(m, std::vector<double>{1, 3}, std::vector<double>{-0.5}).ok);
  CHECK(!samaya::verify_qcp_kkt(m, std::vector<double>{1, 3}, std::vector<double>{0.0}).ok);
  // A Jacobian read as linear (dh/dx = 1) would accept lambda = -1; the right one does not.
  CHECK(!samaya::verify_qcp_kkt(m, x, std::vector<double>{-1.0}).ok);
}

TEST(verify_qcp_kkt_diagonal_term_and_maximization) {
  // max x  s.t.  x^2 <= 4: KKT at x = 2 with y = 1/4 (d = 1 - 2x y = 0; for a maximization a
  // row at its upper bound needs y >= 0). A diagonal term's derivative is 2 v x.
  Model m = continuous(1);
  m.sense = samaya::ObjSense::kMaximize;
  m.obj = {1};
  m.row_lower = {-kInf};
  m.row_upper = {4};
  m.A = samaya::SparseMatrix(1, 1);
  m.row_quadratic = {{0, 0, 0, 1.0}};
  CHECK(samaya::verify_qcp_kkt(m, std::vector<double>{2}, std::vector<double>{0.25}).ok);
  CHECK(!samaya::verify_qcp_kkt(m, std::vector<double>{2}, std::vector<double>{0.5}).ok);
}

TEST(verify_qcp_kkt_certifies_stationarity_not_local_optimality) {
  // min z  s.t.  z - u v = 0,  -1 <= u, v <= 1: the origin is a KKT point (lambda = 1), yet
  // (t, -t, -t^2) is feasible and better arbitrarily close to it. The check must accept it:
  // it certifies first-order stationarity only, which is all the SLP status claims.
  Model m;
  m.obj = {0, 0, 1};
  m.col_lower = {-1, -1, -kInf};
  m.col_upper = {1, 1, kInf};
  m.col_type.assign(3, samaya::VarType::kContinuous);
  m.row_lower = {0};
  m.row_upper = {0};
  m.A = samaya::SparseMatrix::from_triplets(1, 3, {{0, 2, 1}});
  m.row_quadratic = {{0, 0, 1, -1.0}};
  CHECK(samaya::verify_qcp_kkt(m, std::vector<double>{0, 0, 0}, std::vector<double>{1}).ok);
  CHECK(samaya::verify_primal(m, std::vector<double>{0.1, -0.1, -0.01}).ok);
}