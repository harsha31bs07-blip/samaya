#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "samaya/types.hpp"

namespace samaya {

enum class LpMethod : std::uint8_t {
  kAuto,        // Pick per model; concurrent when threads allow.
  kDualSimplex,
  kPrimalSimplex,
  kBarrier,     // Interior point + crossover.
  kPdlp,        // First-order primal-dual method, GPU when available.
  kConcurrent,  // Race dual simplex, barrier and PDLP.
};

struct Params {
  // Limits.
  double time_limit = kInf;  // Seconds of wall-clock time.
  long long node_limit = -1;  // -1 = unlimited.
  int threads = 0;            // 0 = hardware concurrency.

  // Tolerances.
  double primal_feasibility_tol = 1e-7;
  double dual_feasibility_tol = 1e-7;
  double integrality_tol = 1e-6;
  double mip_rel_gap = 1e-4;
  double mip_abs_gap = 1e-6;

  // Algorithm selection.
  LpMethod lp_method = LpMethod::kAuto;
  bool presolve = true;
  bool use_gpu = false;
  // After barrier or PDLP: turn the answer into a basic (vertex) solution by crossover and the
  // simplex, and report that. Off: the interior point is reported when it verifies. A point that
  // fails verification goes to the simplex from the crossover basis either way.
  bool crossover = false;
  // Simplex iterations the crossover may take before the dual simplex starts over from the slack
  // basis: a guess that needs more was poor (dfl001 from a 1e-4 PDLP point: 289,256 iterations
  // and a numerical error, where a fresh start needs about 23,000). -1: automatic,
  // 2 (rows + columns) and at least 10,000.
  long long crossover_max_iterations = -1;

  // Independently re-check every solution before reporting it.
  bool verify = true;

  // MILP start, one value per column of the model (NaN = unknown), e.g. the previous plan when
  // re-planning. Used as the incumbent if feasible, else completed or repaired (see
  // MipOptions::start).
  std::vector<double> mip_start;

  // 0 = silent, 1 = summary, 2 = iteration log, 3 = debug.
  int log_level = 1;
};

}  // namespace samaya
