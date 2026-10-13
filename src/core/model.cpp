#include "samaya/model.hpp"

#include <algorithm>
#include <cmath>
#include <string>

namespace samaya {

const char* to_string(ProblemClass problem_class) {
  switch (problem_class) {
    case ProblemClass::kLP: return "LP";
    case ProblemClass::kQP: return "QP";
    case ProblemClass::kMILP: return "MILP";
    case ProblemClass::kMIQP: return "MIQP";
    case ProblemClass::kQCP: return "QCP";
  }
  return "unknown";
}

Index Model::num_integers() const {
  return static_cast<Index>(std::count(col_type.begin(), col_type.end(), VarType::kInteger));
}

ProblemClass Model::problem_class() const {
  // Quadratic rows win: the solver refuses them together with integers or a quadratic
  // objective, rather than letting the MILP or QP code ignore them.
  if (has_quadratic_rows()) return ProblemClass::kQCP;
  if (is_mip()) return is_qp() ? ProblemClass::kMIQP : ProblemClass::kMILP;
  return is_qp() ? ProblemClass::kQP : ProblemClass::kLP;
}

namespace {

std::string check_bounds(const char* what, const std::vector<double>& lower,
                         const std::vector<double>& upper) {
  for (std::size_t i = 0; i < lower.size(); ++i) {
    if (std::isnan(lower[i]) || std::isnan(upper[i])) {
      return std::string(what) + " " + std::to_string(i) + " has a NaN bound";
    }
    if (lower[i] == kInf || upper[i] == -kInf) {
      return std::string(what) + " " + std::to_string(i) + " has an infinite bound on the wrong side";
    }
  }
  return {};
}

std::string first_crossed(const char* what, const std::vector<double>& lower,
                          const std::vector<double>& upper) {
  for (std::size_t i = 0; i < lower.size(); ++i) {
    if (lower[i] > upper[i]) {
      return std::string(what) + " " + std::to_string(i) + " has lower bound > upper bound";
    }
  }
  return {};
}

}  // namespace

std::string Model::validate() const {
  const auto n = static_cast<std::size_t>(num_cols());
  const auto m = static_cast<std::size_t>(num_rows());
  if (obj.size() != n || col_upper.size() != n || col_type.size() != n) {
    return "column arrays have inconsistent sizes";
  }
  if (row_upper.size() != m) return "row arrays have inconsistent sizes";
  if (!col_names.empty() && col_names.size() != n) return "col_names has the wrong size";
  if (!row_names.empty() && row_names.size() != m) return "row_names has the wrong size";
  if (A.rows() != num_rows() || A.cols() != num_cols()) return "constraint matrix has wrong shape";
  if (!Q.empty() && (Q.rows() != num_cols() || Q.cols() != num_cols())) {
    return "Q matrix has wrong shape";
  }
  if (!std::isfinite(obj_offset)) return "objective offset is not finite";

  for (std::size_t j = 0; j < n; ++j) {
    if (!std::isfinite(obj[j])) return "objective coefficient " + std::to_string(j) + " is not finite";
  }
  for (double v : A.values()) {
    if (!std::isfinite(v)) return "constraint matrix contains a non-finite value";
  }

  const auto q_start = Q.col_start();
  const auto q_index = Q.row_index();
  for (Index j = 0; j < Q.cols(); ++j) {
    for (NnzIndex p = q_start[j]; p < q_start[j + 1]; ++p) {
      if (q_index[p] < j) return "Q must be stored as its lower triangle";
      if (!std::isfinite(Q.values()[p])) return "Q contains a non-finite value";
    }
  }

  for (std::size_t k = 0; k < row_quadratic.size(); ++k) {
    const RowQuadraticTerm& t = row_quadratic[k];
    if (t.row < 0 || t.row >= num_rows() || t.col1 < 0 || t.col2 >= num_cols() ||
        t.col1 > t.col2) {
      return "quadratic row term " + std::to_string(k) + " has an invalid row or column";
    }
    if (!std::isfinite(t.value)) {
      return "quadratic row term " + std::to_string(k) + " is not finite";
    }
  }

  if (std::string e = check_bounds("column", col_lower, col_upper); !e.empty()) return e;
  if (std::string e = check_bounds("row", row_lower, row_upper); !e.empty()) return e;
  return {};
}

std::string Model::crossed_bounds() const {
  if (std::string e = first_crossed("column", col_lower, col_upper); !e.empty()) return e;
  return first_crossed("row", row_lower, row_upper);
}

}  // namespace samaya
