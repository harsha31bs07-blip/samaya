#pragma once

#include <vector>

#include "core/log.hpp"
#include "lp/simplex.hpp"

namespace samaya {

// Options for the PDLP solver.
struct PdlpOptions {
  double tol = 1e-8;              // Target relative tolerance for primal/dual residual and gap.
  long long max_iterations = -1;  // -1: automatic, proportional to problem size.
  double time_limit = kInf;       // Wall-clock time limit in seconds.
  bool use_gpu = false;           // Use GPU kernels (WP2); ignored if SAMAYA_CUDA is off.
};

// Result of a PDLP solve.
struct PdlpResult {
  SimplexStatus status = SimplexStatus::kNumericalError;
  std::vector<double> x;           // Primal solution (n + m).
  std::vector<double> y;           // Dual solution (m).
  std::vector<double> dual_ray;    // Farkas certificate when infeasible (m).
  std::vector<double> primal_ray;  // Improving direction when unbounded (n).
  long long iterations = 0;
  double primal_residual = 0.0;  // ||A x - r|| / (1 + ||A|| ||x|| + ||r||)
  double dual_residual = 0.0;    // ||A^T y + d - c|| / (1 + ||A|| ||y|| + ||c||)
  double gap = 0.0;              // |c^T x - b^T y| / (1 + |c^T x| + |b^T y|)
  SimplexStats stats;
};

// Solves the scaled LP problem using PDLP.
// Works on the LpProblem in computational form: min cost'v s.t. [A -I] v = 0, lower <= v <= upper.
// The solution vectors (x, y) are in the scaled space.
PdlpResult solve_pdlp(const LpProblem& lp, const PdlpOptions& options, const Logger& log);

}  // namespace samaya