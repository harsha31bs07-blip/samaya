#pragma once

#include <vector>

#include "core/log.hpp"
#include "samaya/model.hpp"
#include "samaya/params.hpp"
#include "samaya/solver.hpp"

namespace samaya {

// Solves a model with quadratic rows (Model::row_quadratic; continuous columns, linear
// objective) by trust-region penalty SLP (docs/SLP_PLAN.md): every LP is solved by the dual
// simplex and verified, and a point is reported as kFirstOrderStationary only after
// verify_primal and verify_qcp_kkt pass on the original model. That certifies a feasible KKT
// point, not a local or global optimum. Params::mip_start is the starting point (NaN: unknown).
void solve_qcp(const Model& model, const Params& params, const Logger& log, Result& result);

// For tests: the LP the loop solves at `center` with trust radius `radius` and penalty 1 (its
// columns are the model's, then two elastics per quadratic row), with the scales computed at
// `center`. Lets a test check the trust bounds directly (docs/SLP_PLAN.md section 5).
Model slp_linearization_for_test(const Model& model, const std::vector<double>& center,
                                 double radius);

}  // namespace samaya
