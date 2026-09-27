#include "samaya/status.hpp"

namespace samaya {

const char* to_string(Status status) {
  switch (status) {
    case Status::kNotSolved: return "not_solved";
    case Status::kOptimal: return "optimal";
    case Status::kInfeasible: return "infeasible";
    case Status::kUnbounded: return "unbounded";
    case Status::kInfeasibleOrUnbounded: return "infeasible_or_unbounded";
    case Status::kTimeLimit: return "time_limit";
    case Status::kIterationLimit: return "iteration_limit";
    case Status::kNodeLimit: return "node_limit";
    case Status::kNumericalError: return "numerical_error";
    case Status::kInvalidModel: return "invalid_model";
    case Status::kNotImplemented: return "not_implemented";
    case Status::kNotConvex: return "not_convex";
  }
  return "unknown";
}

}  // namespace samaya
