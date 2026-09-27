#pragma once

#include <vector>

#include "core/log.hpp"
#include "samaya/model.hpp"
#include "samaya/types.hpp"

namespace samaya {

enum class IpmStatus : std::uint8_t {
  kOptimal,         // Residuals and gap within the tolerance (the verifier still decides).
  kNotConvex,       // Q is not positive semidefinite (up to IpmOptions::convexity_tol).
  kIterationLimit,  // Also returned when the iterates diverge (infeasible or unbounded models).
  kTimeLimit,
  kNumericalError,
};

const char* to_string(IpmStatus status);

struct IpmOptions {
  // Relative primal residual, dual residual and gap. The target is tight because the strict
  // verifier needs each complementary pair to be within about 1e-6 of its bound or of zero
  // (products near 1e-12); if the iterates stall before it, the best iterate within
  // acceptable_tol is returned as optimal and the verifier decides.
  double tol = 1e-11;
  double acceptable_tol = 1e-8;
  int stall_iterations = 50;  // Iterations without a 10x better best merit that count as a stall.
  int max_iterations = 300;
  double time_limit = kInf;
  // Q counts as convex if Q + convexity_tol * max|q_jj| I has an LDL^T with positive pivots.
  double convexity_tol = 1e-7;
};

struct IpmResult {
  IpmStatus status = IpmStatus::kNumericalError;
  std::vector<double> x;  // Column values.
  std::vector<double> y;  // Row duals, in the convention of verify_lp_optimality.
  // The polished point (empty if polishing was not possible): the bounds the final iterate
  // treats as active are fixed and the rest solved with zero reduced costs, one KKT solve (the
  // "solution polishing" of OSQP, Stellato et al. 2020). Interior-point iterates stay a little
  // off their bounds, which a strict optimality check rejects; the caller verifies both points.
  std::vector<double> polished_x;
  std::vector<double> polished_y;
  int iterations = 0;
  double primal_residual = kInf;  // Relative, at the last iterate (scaled model).
  double dual_residual = kInf;
  double gap = kInf;
  int regularized_pivots = 0;  // Summed over all factorizations.
};

// Primal-dual interior-point method for convex QP (LP when Q is empty):
//   min c'x + 1/2 x'Qx  s.t.  row_lower <= A x <= row_upper,  col_lower <= x <= col_upper.
//
// Rows with two different bounds get a slack w = a_i'x bounded like the row; fixed columns are
// substituted out; everything else keeps its bounds (free columns stay free). Each iteration
// solves the regularized quasi-definite system
//   [ -(Q + Theta^-1 + rho I)  A_hat' ] [dx]   [r1]
//   [        A_hat            delta I ] [dy] = [r2]
// with LdlFactor (the ordering is computed once), refined against the unregularized matrix.
// Mehrotra predictor-corrector (Mehrotra, SIAM J. Optim. 1992) with the fraction-to-boundary
// rule; separate primal and dual steps for LP, one common step for QP. Ruiz equilibration of
// [Q A'; A 0] first (Ruiz 2001).
//
// Integer columns are treated as continuous. The model's objective sense is respected.
IpmResult solve_ipm(const Model& model, const IpmOptions& options, const Logger& log);

}  // namespace samaya
