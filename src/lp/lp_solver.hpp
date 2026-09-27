#pragma once

#include <vector>

#include "core/log.hpp"
#include "linalg/scaling.hpp"
#include "lp/simplex.hpp"
#include "samaya/model.hpp"
#include "samaya/params.hpp"

namespace samaya {

// Solution of an LP in the space of the original model (unscaled, original objective sense).
struct LpResult {
  SimplexStatus status = SimplexStatus::kNumericalError;
  double objective = 0.0;              // c'x + offset, valid when optimal.
  std::vector<double> col_value;       // x (n).
  std::vector<double> row_activity;    // A x (m).
  std::vector<double> row_dual;        // y with c - A'y = reduced costs (m).
  std::vector<double> col_dual;        // Reduced costs (n).
  std::vector<double> dual_ray;        // Farkas certificate when infeasible (m).
  std::vector<double> primal_ray;      // Improving direction when unbounded (n).
  long long iterations = 0;
  SimplexStats stats;
};

struct LpSolveOptions {
  bool scale = true;
  SimplexOptions simplex;
};

// The model in computational form (minimization), scaled by `scaling`:
// min cost'v s.t. [A -I] v = 0, lower <= v <= upper. Used by solve_lp and by PDLP.
LpProblem make_problem(const Model& model, const Scaling& scaling);

// Solves the LP relaxation of `model` (integrality is ignored; Q must be empty) with the dual
// simplex method.
LpResult solve_lp(const Model& model, const LpSolveOptions& options, const Logger& log);

// Simplex options derived from user-facing parameters.
LpSolveOptions lp_options_from_params(const Params& params);

}  // namespace samaya
