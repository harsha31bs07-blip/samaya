#pragma once

#include <vector>

#include "core/log.hpp"
#include "lp/lp_solver.hpp"
#include "samaya/params.hpp"
#include "samaya/solver.hpp"

namespace samaya {

// Runs the independent verifier on an LP outcome. Returns true when the claim is proven (or is a
// limit status, which claims nothing).
bool verify_lp(const Model& model, const LpResult& lp, Result& result);

// Solves an LP with the dual simplex. If the outcome fails verification the solve is repeated
// with tighter tolerances and then without scaling. Returns the last outcome; `verified` says
// whether it passed (limit statuses claim nothing and count as verified, and nothing is checked
// when Params::verify is off: callers that need a proof must also check the status).
LpResult solve_lp_verified(const Model& model, const Params& params, const Logger& log,
                           Result& result, bool& verified,
                           const std::vector<VarStatus>& start = {},
                           long long max_iterations = -1);

}  // namespace samaya
