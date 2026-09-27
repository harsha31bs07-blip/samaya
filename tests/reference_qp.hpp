#pragma once

// Dense reference for small convex QPs: enumerates the active sets and returns the first KKT
// point. For a convex QP any point satisfying the KKT conditions is optimal, so no search order
// matters. Obviously correct rather than fast: only for a handful of columns and rows.

#include <cmath>
#include <vector>

#include "dense_reference.hpp"
#include "samaya/model.hpp"

namespace samaya::test {

struct ReferenceQpResult {
  enum class Status { kOptimal, kNoKktPoint } status = Status::kNoKktPoint;
  double objective = 0.0;  // In the model's sense, offset included.
  std::vector<double> x;
};

// Requires finite column bounds on every column (a bounded feasible set, so an optimum exists
// whenever the model is feasible) and a convex objective in the model's sense.
inline ReferenceQpResult reference_qp(const Model& model) {
  const int n = model.num_cols();
  const int m = model.num_rows();
  const double sense = model.sense == ObjSense::kMaximize ? -1.0 : 1.0;
  DenseMatrix Q(static_cast<std::size_t>(n), std::vector<double>(static_cast<std::size_t>(n), 0.0));
  {
    const auto s = model.Q.col_start();
    const auto r = model.Q.row_index();
    const auto v = model.Q.values();
    for (int j = 0; j < model.Q.cols(); ++j) {
      for (NnzIndex p = s[j]; p < s[j + 1]; ++p) {
        Q[r[p]][j] = sense * v[p];
        Q[j][r[p]] = sense * v[p];
      }
    }
  }
  std::vector<double> c(static_cast<std::size_t>(n));
  for (int j = 0; j < n; ++j) c[j] = sense * model.obj[j];
  // Constraints: the columns (unit normals), then the rows.
  DenseMatrix g;
  std::vector<double> lo;
  std::vector<double> up;
  for (int j = 0; j < n; ++j) {
    std::vector<double> e(static_cast<std::size_t>(n), 0.0);
    e[j] = 1.0;
    g.push_back(e);
    lo.push_back(model.col_lower[j]);
    up.push_back(model.col_upper[j]);
  }
  {
    DenseMatrix rows(static_cast<std::size_t>(m),
                     std::vector<double>(static_cast<std::size_t>(n), 0.0));
    const auto s = model.A.col_start();
    const auto r = model.A.row_index();
    const auto v = model.A.values();
    for (int j = 0; j < n; ++j) {
      for (NnzIndex p = s[j]; p < s[j + 1]; ++p) rows[r[p]][j] = v[p];
    }
    for (int i = 0; i < m; ++i) {
      bool empty = true;
      for (const double coef : rows[i]) empty = empty && coef == 0.0;
      if (empty) {  // An empty row holds for every x iff 0 is within its bounds.
        if (model.row_lower[i] > 0.0 || model.row_upper[i] < 0.0) return {};
        continue;
      }
      g.push_back(rows[i]);
      lo.push_back(model.row_lower[i]);
      up.push_back(model.row_upper[i]);
    }
  }
  const int k_total = static_cast<int>(g.size());
  // state[k]: 0 inactive, 1 at its lower side, 2 at its upper side (equalities always 1).
  std::vector<int> state(static_cast<std::size_t>(k_total), 0);
  const auto next = [&]() {
    for (int k = 0; k < k_total; ++k) {
      if (lo[k] == up[k]) continue;  // Always active.
      int s = state[k] + 1;
      if (s == 1 && !std::isfinite(lo[k])) s = 2;
      if (s == 2 && !std::isfinite(up[k])) s = 3;
      if (s < 3) {
        state[k] = s;
        return true;
      }
      state[k] = 0;
    }
    return false;
  };
  for (int k = 0; k < k_total; ++k) {
    if (lo[k] == up[k]) state[k] = 1;
  }
  ReferenceQpResult best;
  do {
    std::vector<int> active;
    for (int k = 0; k < k_total; ++k) {
      if (state[k] != 0) active.push_back(k);
    }
    if (static_cast<int>(active.size()) > n) continue;
    const int a = static_cast<int>(active.size());
    DenseMatrix kkt(static_cast<std::size_t>(n + a),
                    std::vector<double>(static_cast<std::size_t>(n + a), 0.0));
    std::vector<double> rhs(static_cast<std::size_t>(n + a), 0.0);
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) kkt[i][j] = Q[i][j];
      rhs[i] = -c[i];
    }
    for (int t = 0; t < a; ++t) {
      const int k = active[t];
      for (int j = 0; j < n; ++j) {
        kkt[n + t][j] = g[k][j];
        kkt[j][n + t] = -g[k][j];  // Q x + c - G' lambda = 0.
      }
      rhs[n + t] = state[k] == 1 ? lo[k] : up[k];
    }
    std::vector<double> z;
    if (!dense_solve(kkt, rhs, z)) continue;
    bool ok = true;
    for (int k = 0; k < k_total && ok; ++k) {
      double act = 0.0;
      double scale = 1.0;
      for (int j = 0; j < n; ++j) {
        act += g[k][j] * z[j];
        scale += std::fabs(g[k][j] * z[j]);
      }
      const double tol = 1e-9 * (scale + std::fabs(std::isfinite(lo[k]) ? lo[k] : 0.0) +
                                 std::fabs(std::isfinite(up[k]) ? up[k] : 0.0));
      if (act < lo[k] - tol || act > up[k] + tol) ok = false;
    }
    // Multipliers: a lower side (g'x >= lo) needs lambda >= 0, an upper side lambda <= 0.
    for (int t = 0; t < a && ok; ++t) {
      const int k = active[t];
      const double lambda = z[static_cast<std::size_t>(n + t)];
      if (lo[k] == up[k]) continue;
      if (state[k] == 1 && lambda < -1e-9) ok = false;
      if (state[k] == 2 && lambda > 1e-9) ok = false;
    }
    if (!ok) continue;
    best.status = ReferenceQpResult::Status::kOptimal;
    best.x.assign(z.begin(), z.begin() + n);
    double obj = 0.0;
    for (int i = 0; i < n; ++i) {
      obj += c[i] * z[i];
      for (int j = 0; j < n; ++j) obj += 0.5 * z[i] * Q[i][j] * z[j];
    }
    best.objective = sense * obj + model.obj_offset;
    return best;
  } while (next());
  return best;
}

}  // namespace samaya::test
