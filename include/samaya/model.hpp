#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "samaya/sparse_matrix.hpp"
#include "samaya/types.hpp"

namespace samaya {

enum class ObjSense : std::int8_t { kMinimize = 1, kMaximize = -1 };

enum class VarType : std::uint8_t { kContinuous, kInteger };

// kQCP: quadratic terms in constraint rows (docs/SLP_PLAN.md), solved to a first-order
// stationary point by SLP.
enum class ProblemClass : std::uint8_t { kLP, kQP, kMILP, kMIQP, kQCP };

const char* to_string(ProblemClass problem_class);

// A quadratic term value * x_col1 * x_col2 in the activity of `row`, with col1 <= col2 and one
// term per (row, col1, col2). x'Qx in a row (the CPLEX/Gurobi QCMATRIX convention) gives
// value = Q_ii on the diagonal and Q_ij + Q_ji off it.
struct RowQuadraticTerm {
  Index row = 0;
  Index col1 = 0;
  Index col2 = 0;
  double value = 0.0;
};

// An optimization model in the form
//
//   min / max   c'x + 1/2 x'Qx + obj_offset
//   subject to  row_lower <= A x + q(x) <= row_upper
//               col_lower <=  x  <= col_upper
//               x_j integer for col_type[j] == kInteger
//
// Infinite bounds are represented by +/- kInf. Q is stored as its lower triangle (row >= col)
// and is empty for linear objectives. q(x) holds the row_quadratic terms (bilinear blend
// qualities, for instance) and is empty for linear rows.
struct Model {
  std::string name;
  ObjSense sense = ObjSense::kMinimize;
  double obj_offset = 0.0;

  std::vector<double> obj;
  std::vector<double> col_lower;
  std::vector<double> col_upper;
  std::vector<VarType> col_type;
  std::vector<std::string> col_names;

  std::vector<double> row_lower;
  std::vector<double> row_upper;
  std::vector<std::string> row_names;

  SparseMatrix A;
  SparseMatrix Q;
  std::vector<RowQuadraticTerm> row_quadratic;

  Index num_rows() const { return static_cast<Index>(row_lower.size()); }
  Index num_cols() const { return static_cast<Index>(col_lower.size()); }
  Index num_integers() const;
  bool is_mip() const { return num_integers() > 0; }
  bool is_qp() const { return !Q.empty(); }
  bool has_quadratic_rows() const { return !row_quadratic.empty(); }
  ProblemClass problem_class() const;

  // Returns an empty string when the model is consistent, otherwise a description of the first
  // problem found (mismatched sizes, NaNs, infinite bounds on the wrong side, ...). Crossed
  // bounds (lower > upper) are valid input: they make the model infeasible.
  std::string validate() const;

  // Describes the first column or row whose lower bound exceeds its upper bound, if any.
  std::string crossed_bounds() const;
};

}  // namespace samaya
