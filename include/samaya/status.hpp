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
};

const char* to_string(Status status);

}  // namespace samaya
