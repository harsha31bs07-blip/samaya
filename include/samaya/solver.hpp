#pragma once

#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "samaya/model.hpp"
#include "samaya/params.hpp"
#include "samaya/status.hpp"

namespace samaya {

struct Result {
  Status status = Status::kNotSolved;
  std::string message;

  // Objective value of the returned solution and the best proven bound (equal for LP/QP optima).
  double objective = std::numeric_limits<double>::quiet_NaN();
  double dual_bound = std::numeric_limits<double>::quiet_NaN();

  std::vector<double> col_value;     // x
  std::vector<double> row_activity;  // A x
  std::vector<double> row_dual;      // y with c - A'y = reduced costs (LP/QP only)
  std::vector<double> col_dual;      // reduced costs (LP/QP only)

  // Certificates: a Farkas multiplier per row when infeasible, an improving direction when
  // unbounded.
  std::vector<double> infeasibility_certificate;
  std::vector<double> unbounded_ray;

  // Whether the status and solution passed the independent verifier (Params::verify), and the
  // largest relative violations it measured.
  bool verified = false;
  double max_primal_violation = 0.0;
  double max_dual_violation = 0.0;

  double solve_seconds = 0.0;
  long long simplex_iterations = 0;
  long long barrier_iterations = 0;
  long long nodes = 0;
};

class Solver {
 public:
  explicit Solver(Params params = {}) : params_(std::move(params)) {}

  const Params& params() const { return params_; }
  Params& params() { return params_; }

  Result solve(const Model& model) const;

 private:
  Params params_;
};

// Semantic version of the library, e.g. "0.1.0".
const char* version();

// The GPU samaya can use for PDLP (Params::use_gpu): whether this build has the CUDA kernels, and
// the name of the NVIDIA GPU it would run on (empty when there is no usable GPU or driver).
struct GpuInfo {
  bool built_with_cuda = false;
  std::string device;
};
GpuInfo gpu_info();

}  // namespace samaya
