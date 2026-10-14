#include "lp/crossover.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>

#include "linalg/basis_factor.hpp"

namespace samaya {

namespace {

// Match Simplex::rebuild's cap so crossover cannot spend indefinitely on numerical repairs.
constexpr int kMaxFactorizations = 4;

// Distances and reduced costs are compared relative to the size of the value and of the cost,
// so that the score does not depend on the units of a column or a row.
double relative(double amount, double size) { return amount / (1.0 + std::fabs(size)); }

// Nonbasic status at the bound nearest to value (lower on a tie, which also covers fixed
// variables); free variables sit at zero.
VarStatus nearest_bound(double value, double lower, double upper) {
  const bool has_lower = lower > -kInf;
  const bool has_upper = upper < kInf;
  if (has_lower && has_upper) {
    return value - lower <= upper - value ? VarStatus::kAtLower : VarStatus::kAtUpper;
  }
  if (has_lower) return VarStatus::kAtLower;
  if (has_upper) return VarStatus::kAtUpper;
  return VarStatus::kAtZero;
}

}  // namespace

std::vector<VarStatus> crossover_basis(const Model& model, std::span<const double> x,
                                       std::span<const double> y, CrossoverStats* stats) {
  const Index m = model.num_rows();
  const Index n = model.num_cols();
  assert(x.size() == static_cast<std::size_t>(n) && y.size() == static_cast<std::size_t>(m));
  const Index nt = n + m;

  // The values, bounds and reduced costs of the n + m variables of the computational form. The
  // logical of row i carries the row activity a_i'x and the row bounds; its reduced cost is the
  // row dual (its column in [A -I] is -e_i and its cost is zero).
  std::vector<double> activity(static_cast<std::size_t>(m), 0.0);
  model.A.multiply(x, activity);
  std::vector<double> aty(static_cast<std::size_t>(n), 0.0);
  model.A.multiply_transpose(y, aty);

  std::vector<double> value(static_cast<std::size_t>(nt));
  std::vector<double> lower(static_cast<std::size_t>(nt));
  std::vector<double> upper(static_cast<std::size_t>(nt));
  // Relative distance to the nearest bound.
  std::vector<double> slack(static_cast<std::size_t>(nt));
  std::vector<double> score(static_cast<std::size_t>(nt));
  for (Index j = 0; j < nt; ++j) {
    double reduced_cost = 0.0;
    double cost = 0.0;
    if (j < n) {
      value[j] = x[j];
      lower[j] = model.col_lower[j];
      upper[j] = model.col_upper[j];
      cost = model.obj[j];
      reduced_cost = cost - aty[j];
    } else {
      const Index i = j - n;
      value[j] = activity[i];
      lower[j] = model.row_lower[i];
      upper[j] = model.row_upper[i];
      reduced_cost = y[i];
    }
    // A point slightly outside a bound counts as on it.
    const double to_lower = lower[j] > -kInf ? std::max(0.0, value[j] - lower[j]) : kInf;
    const double to_upper = upper[j] < kInf ? std::max(0.0, upper[j] - value[j]) : kInf;
    slack[j] = relative(std::min(to_lower, to_upper), value[j]);
    const double dual = relative(std::fabs(reduced_cost), cost);
    // In [0, 1]: 1 for a variable strictly between its bounds with a zero reduced cost (basic),
    // 0 for one on a bound with a nonzero reduced cost (nonbasic). Free variables score 1.
    score[j] = slack[j] == kInf ? 1.0 : slack[j] + dual > 0.0 ? slack[j] / (slack[j] + dual) : 0.0;
  }

  // The m best scores are basic; ties go to the variable farther from its bounds.
  std::vector<Index> order(static_cast<std::size_t>(nt));
  std::iota(order.begin(), order.end(), Index{0});
  std::stable_sort(order.begin(), order.end(), [&](Index a, Index b) {
    if (score[a] != score[b]) return score[a] > score[b];
    return slack[a] > slack[b];
  });
  std::vector<VarStatus> status(static_cast<std::size_t>(nt));
  for (Index k = 0; k < nt; ++k) {
    const Index j = order[k];
    status[j] = k < m ? VarStatus::kBasic : nearest_bound(value[j], lower[j], upper[j]);
  }

  // Use the same variable order as the simplex warm start, and check independence before the
  // point values are lost. A repair's leaving variable must use this point, not a fresh zero.
  std::vector<Index> basic;
  basic.reserve(static_cast<std::size_t>(m));
  for (Index j = 0; j < nt; ++j) {
    if (status[j] == VarStatus::kBasic) basic.push_back(j);
  }
  CrossoverStats measurements;
  BasisFactor::Options options;
  options.collect_stats = stats != nullptr;
  BasisFactor factor(model.A, options);
  for (int attempt = 0; attempt < kMaxFactorizations; ++attempt) {
    ++measurements.factorizations;
    measurements.rank_deficiency = factor.factorize(basic);
    if (measurements.rank_deficiency == 0 || attempt + 1 == kMaxFactorizations) break;
    const auto& positions = factor.singular_positions();
    const auto& rows = factor.unpivoted_rows();
    bool stopped = false;
    for (std::size_t k = 0; k < positions.size(); ++k) {
      const Index entering = n + rows[k];
      if (status[entering] == VarStatus::kBasic) {
        stopped = true;
        break;
      }
      const Index leaving = basic[positions[k]];
      status[leaving] = nearest_bound(value[leaving], lower[leaving], upper[leaving]);
      basic[positions[k]] = entering;
      status[entering] = VarStatus::kBasic;
      ++measurements.repairs;
    }
    if (stopped) break;
  }
  if (stats != nullptr) {
    for (const Index j : basic) {
      if (j < n) {
        ++measurements.basic_structurals;
        if (lower[j] == upper[j]) ++measurements.basic_fixed_structurals;
      } else if (lower[j] == upper[j]) {
        ++measurements.basic_equality_logicals;
      } else {
        ++measurements.basic_other_logicals;
      }
    }
    measurements.factorization_seconds = factor.stats().rebuild_seconds;
    *stats = measurements;
  }
  return status;
}

}  // namespace samaya
