#pragma once

#include <span>
#include <string>

#include "samaya/model.hpp"

namespace samaya {

// Independent solution checks. They use only the original model and the reported vectors,
// recompute every derived quantity (row activities, reduced costs) in extended precision, and
// share no code with the solvers.
//
// Violations are relative: a bound or row violation is divided by 1 + the magnitude of the
// quantities involved (the bound and the largest term |a_ij x_j|), and a dual violation by 1 +
// the magnitude of the terms forming the reduced cost.
struct VerifyTolerances {
  double primal = 1e-6;
  double dual = 1e-6;
  double integrality = 1e-6;
};

struct VerifyReport {
  bool ok = true;
  double max_bound_violation = 0.0;
  double max_row_violation = 0.0;
  double max_integrality_violation = 0.0;
  double max_dual_violation = 0.0;
  double objective = 0.0;  // c'x + 1/2 x'Qx + offset, summed in extended precision.
  std::string message;     // Describes the worst failure when !ok.
};

// Checks bounds, rows and integrality of x and recomputes the objective.
VerifyReport verify_primal(const Model& model, std::span<const double> x,
                           const VerifyTolerances& tol = {});

// Checks x for primal feasibility and the row duals y for LP optimality: with d = c - A'y
// (recomputed), every variable and row that is not at its lower (upper) bound must have a
// reduced cost that does not allow improvement, respecting the objective sense.
VerifyReport verify_lp_optimality(const Model& model, std::span<const double> x,
                                  std::span<const double> y, const VerifyTolerances& tol = {});

// The same for a convex QP: the reduced costs are d = c + Qx - A'y, with Q the symmetric matrix
// whose lower triangle is model.Q (objective c'x + 1/2 x'Qx). Convexity itself is not checked.
VerifyReport verify_qp_optimality(const Model& model, std::span<const double> x,
                                  std::span<const double> y, const VerifyTolerances& tol = {});

// Checks a Farkas certificate y (one multiplier per row): the value of y'(A x) - y'r over all x
// within the column bounds and r within the row bounds must exclude zero, which proves that
// A x = r has no solution in the bounds.
VerifyReport verify_infeasibility(const Model& model, std::span<const double> y,
                                  const VerifyTolerances& tol = {});

// Checks an unbounded ray v: it improves the objective, and moving along it from any feasible
// point never violates a bound or a row (it only moves in directions where they are infinite).
VerifyReport verify_unbounded_ray(const Model& model, std::span<const double> ray,
                                  const VerifyTolerances& tol = {});

}  // namespace samaya
