#pragma once

#include <vector>

#include "core/log.hpp"
#include "samaya/model.hpp"
#include "samaya/status.hpp"

namespace samaya {

struct MiqpOptions {
  double time_limit = kInf;
  long long node_limit = -1;  // -1 = unlimited.
  double rel_gap = 1e-4;
  double abs_gap = 1e-6;
  double integrality_tol = 1e-6;
};

struct MiqpOutcome {
  Status status = Status::kNotSolved;
  double objective = kInf;  // Of x, in the model's sense (+-inf without a solution).
  double bound = -kInf;     // Best proven bound, in the model's sense.
  std::vector<double> x;    // Best integer solution (empty if none).
  long long nodes = 0;
  long long ipm_iterations = 0;
};

// Mixed-integer convex QP by branch and bound over the QP relaxations, each solved by the
// interior-point method (solve_ipm). Best-first node selection, branching on the most
// fractional integer column; at the root and every kRoundingFrequency nodes the relaxation is
// rounded, its integers fixed and the remaining continuous QP solved (a rounding heuristic).
//
// A node is pruned as infeasible only when the dual simplex, on the node's constraints alone,
// returns a verified Farkas certificate; a node whose relaxation cannot be solved otherwise makes
// the search incomplete (kNumericalError unless a limit is hit first). Every incumbent passes
// verify_primal (bounds, rows, integrality) before it is accepted. Relaxation values are the
// objective of the interior-point solution, converged to a relative gap of about 1e-9, less a
// relative margin of kBoundMargin; the optimality proof rests on that accuracy.
MiqpOutcome solve_miqp(const Model& model, const MiqpOptions& options, const Logger& log);

}  // namespace samaya
