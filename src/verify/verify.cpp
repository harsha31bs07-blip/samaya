#include "samaya/verify.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

#include "core/xreal.hpp"

namespace samaya {
namespace {

// Extended precision for sums and products: a double-double, the same bits on every platform.
using Real = Xreal;

void fail(VerifyReport& report, double violation, double tolerance, const char* what, Index index) {
  if (violation <= tolerance) return;
  if (report.ok || report.message.empty()) {
    char buf[160];
    std::snprintf(buf, sizeof buf, "%s %d violated by %.3e (tolerance %.1e)", what, index,
                  violation, tolerance);
    report.message = buf;
  }
  report.ok = false;
}

// Relative violation of lo <= v <= up given the magnitude of the terms forming v.
double bound_violation(double v, double lo, double up, double scale) {
  double violation = 0.0;
  if (v < lo) violation = (lo - v) / (1.0 + std::max(std::fabs(lo), scale));
  if (v > up) violation = (v - up) / (1.0 + std::max(std::fabs(up), scale));
  return violation;
}

}  // namespace

VerifyReport verify_primal(const Model& model, std::span<const double> x,
                           const VerifyTolerances& tol) {
  VerifyReport report;
  const Index n = model.num_cols();
  const Index m = model.num_rows();
  if (x.size() != static_cast<std::size_t>(n)) {
    report.ok = false;
    report.message = "solution vector has the wrong length";
    return report;
  }

  for (Index j = 0; j < n; ++j) {
    if (!std::isfinite(x[j])) {
      report.ok = false;
      report.message = "solution contains a non-finite value";
      return report;
    }
    const double v = bound_violation(x[j], model.col_lower[j], model.col_upper[j], 0.0);
    report.max_bound_violation = std::max(report.max_bound_violation, v);
    fail(report, v, tol.primal, "column bound", j);
    if (model.col_type[j] == VarType::kInteger) {
      const double frac = std::fabs(x[j] - std::round(x[j]));
      report.max_integrality_violation = std::max(report.max_integrality_violation, frac);
      fail(report, frac, tol.integrality, "integrality of column", j);
    }
  }

  std::vector<Real> activity(static_cast<std::size_t>(m), Real());
  std::vector<double> magnitude(static_cast<std::size_t>(m), 0.0);
  const auto start = model.A.col_start();
  const auto index = model.A.row_index();
  const auto value = model.A.values();
  for (Index j = 0; j < n; ++j) {
    for (NnzIndex p = start[j]; p < start[j + 1]; ++p) {
      const Real term = static_cast<Real>(value[p]) * x[j];
      activity[index[p]] += term;
      magnitude[index[p]] = std::max(magnitude[index[p]], static_cast<double>(fabs(term)));
    }
  }
  for (Index i = 0; i < m; ++i) {
    const double v = bound_violation(static_cast<double>(activity[i]), model.row_lower[i],
                                     model.row_upper[i], magnitude[i]);
    report.max_row_violation = std::max(report.max_row_violation, v);
    fail(report, v, tol.primal, "row", i);
  }

  Real objective = model.obj_offset;
  for (Index j = 0; j < n; ++j) objective += static_cast<Real>(model.obj[j]) * x[j];
  const auto q_start = model.Q.col_start();
  const auto q_index = model.Q.row_index();
  const auto q_value = model.Q.values();
  for (Index j = 0; j < model.Q.cols(); ++j) {
    for (NnzIndex p = q_start[j]; p < q_start[j + 1]; ++p) {
      // Lower triangle: off-diagonal entries appear once but stand for two symmetric terms.
      const Real term = static_cast<Real>(q_value[p]) * x[q_index[p]] * x[j];
      objective += q_index[p] == j ? 0.5 * term : term;
    }
  }
  report.objective = static_cast<double>(objective);
  return report;
}

namespace {

// The optimality conditions of LP and convex QP: with d = c + Qx - A'y (Q empty for LP), every
// column and row not at its lower (upper) bound must have a reduced cost that allows no
// improvement. `quadratic` says whether Q enters d.
VerifyReport verify_kkt(const Model& model, std::span<const double> x, std::span<const double> y,
                        const VerifyTolerances& tol, bool quadratic) {
  VerifyReport report = verify_primal(model, x, tol);
  const Index n = model.num_cols();
  const Index m = model.num_rows();
  if (y.size() != static_cast<std::size_t>(m)) {
    report.ok = false;
    report.message = "dual vector has the wrong length";
    return report;
  }
  // For maximization the optimality conditions hold with the reduced costs negated.
  const double sense = model.sense == ObjSense::kMaximize ? -1.0 : 1.0;
  const auto at_bound = [&](double v, double bound, double scale) {
    return std::isfinite(bound) &&
           std::fabs(v - bound) <= tol.primal * (1.0 + std::max(std::fabs(bound), scale));
  };

  const auto start = model.A.col_start();
  const auto index = model.A.row_index();
  const auto value = model.A.values();
  std::vector<Real> activity(static_cast<std::size_t>(m), Real());
  std::vector<double> magnitude(static_cast<std::size_t>(m), 0.0);
  // (Qx)_j from the lower triangle, in extended precision.
  std::vector<Real> qx(static_cast<std::size_t>(n), Real());
  std::vector<double> q_scale(static_cast<std::size_t>(n), 0.0);
  if (quadratic) {
    const auto q_start = model.Q.col_start();
    const auto q_index = model.Q.row_index();
    const auto q_value = model.Q.values();
    for (Index j = 0; j < model.Q.cols(); ++j) {
      for (NnzIndex p = q_start[j]; p < q_start[j + 1]; ++p) {
        const Index i = q_index[p];
        const Real a = static_cast<Real>(q_value[p]) * x[j];
        qx[i] += a;
        q_scale[i] += static_cast<double>(fabs(a));
        if (i != j) {
          const Real b = static_cast<Real>(q_value[p]) * x[i];
          qx[j] += b;
          q_scale[j] += static_cast<double>(fabs(b));
        }
      }
    }
  }
  for (Index j = 0; j < n; ++j) {
    Real d = model.obj[j] + qx[j];
    double scale = std::fabs(model.obj[j]) + q_scale[j];
    for (NnzIndex p = start[j]; p < start[j + 1]; ++p) {
      const Real term = static_cast<Real>(value[p]) * y[index[p]];
      d -= term;
      scale += static_cast<double>(fabs(term));
      activity[index[p]] += static_cast<Real>(value[p]) * x[j];
      magnitude[index[p]] =
          std::max(magnitude[index[p]], std::fabs(value[p] * x[j]));
    }
    const double dj = sense * static_cast<double>(d);
    const bool lower = at_bound(x[j], model.col_lower[j], 0.0);
    const bool upper = at_bound(x[j], model.col_upper[j], 0.0);
    double violation = 0.0;
    if (!lower) violation = std::max(violation, dj);
    if (!upper) violation = std::max(violation, -dj);
    violation /= 1.0 + scale;
    report.max_dual_violation = std::max(report.max_dual_violation, violation);
    fail(report, violation, tol.dual, "reduced cost of column", j);
  }
  for (Index i = 0; i < m; ++i) {
    if (!std::isfinite(y[i])) {
      report.ok = false;
      report.message = "dual vector contains a non-finite value";
      return report;
    }
    // A row is a variable r_i = a_i'x with reduced cost y_i.
    const double yi = sense * y[i];
    const double r = static_cast<double>(activity[i]);
    const bool lower = at_bound(r, model.row_lower[i], magnitude[i]);
    const bool upper = at_bound(r, model.row_upper[i], magnitude[i]);
    double violation = 0.0;
    if (!lower) violation = std::max(violation, yi);
    if (!upper) violation = std::max(violation, -yi);
    violation /= 1.0 + std::fabs(yi);
    report.max_dual_violation = std::max(report.max_dual_violation, violation);
    fail(report, violation, tol.dual, "dual of row", i);
  }
  return report;
}

}  // namespace

VerifyReport verify_lp_optimality(const Model& model, std::span<const double> x,
                                  std::span<const double> y, const VerifyTolerances& tol) {
  return verify_kkt(model, x, y, tol, false);
}

VerifyReport verify_qp_optimality(const Model& model, std::span<const double> x,
                                  std::span<const double> y, const VerifyTolerances& tol) {
  return verify_kkt(model, x, y, tol, true);
}

VerifyReport verify_infeasibility(const Model& model, std::span<const double> y,
                                  const VerifyTolerances& tol) {
  VerifyReport report;
  const Index n = model.num_cols();
  const Index m = model.num_rows();
  if (y.size() != static_cast<std::size_t>(m)) {
    report.ok = false;
    report.message = "certificate has the wrong length";
    return report;
  }
  double y_max = 0.0;
  for (const double v : y) {
    if (!std::isfinite(v)) {
      report.ok = false;
      report.message = "certificate contains a non-finite value";
      return report;
    }
    y_max = std::max(y_max, std::fabs(v));
  }
  if (y_max == 0.0) {
    report.ok = false;
    report.message = "certificate is zero";
    return report;
  }

  // f(x, r) = w'x - y'r with w = A'y (y normalized to max |y_i| = 1). Bound f over the box.
  Real lo = 0.0;
  Real hi = 0.0;
  double scale = 0.0;
  const auto add_term = [&](Real coef, double lower, double upper) {
    if (coef > 0) {
      lo += lower > -kInf ? coef * lower : -INFINITY;
      hi += upper < kInf ? coef * upper : INFINITY;
    } else if (coef < 0) {
      lo += upper < kInf ? coef * upper : -INFINITY;
      hi += lower > -kInf ? coef * lower : INFINITY;
    }
  };
  const auto start = model.A.col_start();
  const auto index = model.A.row_index();
  const auto value = model.A.values();
  for (Index j = 0; j < n; ++j) {
    Real w = 0.0;
    double w_scale = 0.0;
    for (NnzIndex p = start[j]; p < start[j + 1]; ++p) {
      const Real term = static_cast<Real>(value[p]) * (y[index[p]] / y_max);
      w += term;
      w_scale += static_cast<double>(fabs(term));
    }
    // Cancellation noise must not turn an infinite bound into an infinite range.
    if (std::fabs(static_cast<double>(w)) <= 1e-9 * w_scale) w = 0.0;
    add_term(w, model.col_lower[j], model.col_upper[j]);
    if (w != 0) {
      scale = std::max({scale, std::fabs(static_cast<double>(w) * model.col_lower[j]),
                        std::fabs(static_cast<double>(w) * model.col_upper[j])});
    }
  }
  for (Index i = 0; i < m; ++i) {
    const Real coef = -static_cast<Real>(y[i] / y_max);
    add_term(coef, model.row_lower[i], model.row_upper[i]);
    if (coef != 0) {
      scale = std::max({scale, std::fabs(static_cast<double>(coef) * model.row_lower[i]),
                        std::fabs(static_cast<double>(coef) * model.row_upper[i])});
    }
  }
  if (!std::isfinite(scale)) scale = 0.0;
  const double margin = tol.primal * (1.0 + scale);
  if (lo > margin || hi < -margin) return report;
  report.ok = false;
  char buf[160];
  std::snprintf(buf, sizeof buf, "certificate range [%.3e, %.3e] does not exclude zero",
                static_cast<double>(lo), static_cast<double>(hi));
  report.message = buf;
  return report;
}

VerifyReport verify_unbounded_ray(const Model& model, std::span<const double> ray,
                                  const VerifyTolerances& tol) {
  VerifyReport report;
  const Index n = model.num_cols();
  const Index m = model.num_rows();
  if (ray.size() != static_cast<std::size_t>(n)) {
    report.ok = false;
    report.message = "ray has the wrong length";
    return report;
  }
  double v_max = 0.0;
  for (const double v : ray) v_max = std::max(v_max, std::fabs(v));
  if (!(v_max > 0.0) || !std::isfinite(v_max)) {
    report.ok = false;
    report.message = "ray is zero or not finite";
    return report;
  }

  const double sense = model.sense == ObjSense::kMaximize ? -1.0 : 1.0;
  Real slope = 0.0;
  double c_max = 0.0;
  for (Index j = 0; j < n; ++j) {
    const double v = ray[j] / v_max;
    slope += static_cast<Real>(sense * model.obj[j]) * v;
    c_max = std::max(c_max, std::fabs(model.obj[j]));
    double violation = 0.0;
    if (v < 0.0 && model.col_lower[j] > -kInf) violation = -v;
    if (v > 0.0 && model.col_upper[j] < kInf) violation = v;
    report.max_bound_violation = std::max(report.max_bound_violation, violation);
    fail(report, violation, tol.primal, "ray leaves the bounds of column", j);
  }
  std::vector<Real> av(static_cast<std::size_t>(m), Real());
  std::vector<double> magnitude(static_cast<std::size_t>(m), 0.0);
  const auto start = model.A.col_start();
  const auto index = model.A.row_index();
  const auto value = model.A.values();
  for (Index j = 0; j < n; ++j) {
    for (NnzIndex p = start[j]; p < start[j + 1]; ++p) {
      const double term = value[p] * ray[j] / v_max;
      av[index[p]] += term;
      magnitude[index[p]] = std::max(magnitude[index[p]], std::fabs(term));
    }
  }
  for (Index i = 0; i < m; ++i) {
    const double a = static_cast<double>(av[i]);
    double violation = 0.0;
    if (a < 0.0 && model.row_lower[i] > -kInf) violation = -a / (1.0 + magnitude[i]);
    if (a > 0.0 && model.row_upper[i] < kInf) violation = a / (1.0 + magnitude[i]);
    report.max_row_violation = std::max(report.max_row_violation, violation);
    fail(report, violation, tol.primal, "ray leaves the bounds of row", i);
  }
  if (!(slope < -tol.dual * (1.0 + c_max))) {
    report.ok = false;
    report.message = "ray does not improve the objective";
  }
  return report;
}

}  // namespace samaya
