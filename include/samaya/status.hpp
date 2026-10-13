#pragma once

#include <cstdint>

namespace samaya {

enum class Status : std::uint8_t {
  kNotSolved,
  kOptimal,
  kInfeasible,
  kUnbounded,
  kInfeasibleOrUnbounded,
  kTimeLimit,
  kIterationLimit,
  kNodeLimit,
  kNumericalError,
  kInvalidModel,
  kNotImplemented,
  kNotConvex,  // A QP whose Q is not positive semidefinite (only convex QP is solved).
  // A QCP (quadratic rows, docs/SLP_PLAN.md): a feasible point where the KKT conditions hold.
  // Local optimality is not certified, and the point need not be a global optimum.
  kFirstOrderStationary,
};

const char* to_string(Status status);

}  // namespace samaya
