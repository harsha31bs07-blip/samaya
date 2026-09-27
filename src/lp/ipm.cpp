#include "lp/ipm.hpp"

#include <algorithm>
#include <cmath>
#include <utility>

#include "linalg/ldl.hpp"

namespace samaya {

namespace {

// Fraction of the step to the boundary taken (keeps the iterates strictly interior).
constexpr double kStepFraction = 0.995;
// Static regularization of the KKT system (scaled model): -rho on the primal block, +delta on the
// dual block. Small enough that iterative refinement against the unregularized system removes
// their effect, large enough to keep the matrix quasi-definite with free columns or dependent
// rows (Friedlander and Orban, Math. Program. Comput. 2012).
constexpr double kPrimalRegularization = 1e-8;
constexpr double kDualRegularization = 1e-8;
// Pivots of the wrong sign or smaller than kPivotTol are replaced by +-kPivotReplacement
// (Clarabel's defaults: 1e-13 and 2e-7).
constexpr double kPivotTol = 1e-13;
constexpr double kPivotReplacement = 2e-7;
constexpr int kRefinementSteps = 3;
constexpr int kRuizPasses = 10;
// Iterates larger than this (relative to the start) mean the model is infeasible or unbounded.
constexpr double kDivergence = 1e12;
// Polishing pins an active bound through a diagonal of this size (scaled model): the solve then
// puts the column within about 1/kPolishWeight of the bound, and it is snapped exactly after.
constexpr double kPolishWeight = 1e10;
// Active-set rounds of polishing, and its tolerance for a bound or reduced-cost violation.
constexpr int kPolishRounds = 10;
constexpr double kPolishTol = 1e-9;
// A step shorter than kShortStep makes the next iteration center with sigma >= kShortStepSigma.
constexpr double kShortStep = 0.1;
constexpr double kShortStepSigma = 0.5;

// Largest magnitude; NaN if any entry is NaN (so a broken iterate never looks converged).
double inf_norm(const std::vector<double>& v) {
  double m = 0.0;
  for (const double x : v) {
    if (!(std::fabs(x) <= m)) m = std::fabs(x);
  }
  return m;
}

double dot(const std::vector<double>& a, const std::vector<double>& b) {
  double s = 0.0;
  for (std::size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
  return s;
}

// The problem in interior-point form: min c'z + 1/2 z'Qz, A z = b, l <= z <= u.
struct Form {
  Index n = 0;  // Columns z (kept model columns, then row slacks).
  Index m = 0;  // Equality constraints.
  std::vector<Triplet> a;  // A (m x n).
  std::vector<Triplet> q;  // Q (n x n), both triangles.
  std::vector<double> b, c, l, u;
};

}  // namespace

const char* to_string(IpmStatus status) {
  switch (status) {
    case IpmStatus::kOptimal: return "optimal";
    case IpmStatus::kNotConvex: return "not convex";
    case IpmStatus::kIterationLimit: return "iteration limit";
    case IpmStatus::kTimeLimit: return "time limit";
    case IpmStatus::kNumericalError: return "numerical error";
  }
  return "unknown";
}

IpmResult solve_ipm(const Model& model, const IpmOptions& options, const Logger& log) {
  const Timer timer;
  IpmResult result;
  const Index n0 = model.num_cols();
  const Index m0 = model.num_rows();
  const double sense = model.sense == ObjSense::kMaximize ? -1.0 : 1.0;

  // ---- Interior-point form: fixed columns substituted, ranged and one-sided rows get slacks.
  Form f;
  std::vector<Index> col_pos(static_cast<std::size_t>(n0), -1);
  for (Index j = 0; j < n0; ++j) {
    if (model.col_lower[j] != model.col_upper[j]) col_pos[j] = f.n++;
  }
  const Index kept = f.n;
  std::vector<double> shift(static_cast<std::size_t>(m0), 0.0);  // Fixed columns' activity.
  const auto start = model.A.col_start();
  const auto index = model.A.row_index();
  const auto value = model.A.values();
  for (Index j = 0; j < n0; ++j) {
    if (col_pos[j] >= 0) continue;
    for (NnzIndex p = start[j]; p < start[j + 1]; ++p) {
      shift[index[p]] += value[p] * model.col_lower[j];
    }
  }
  std::vector<Index> row_pos(static_cast<std::size_t>(m0), -1);
  std::vector<Index> slack_of(static_cast<std::size_t>(m0), -1);
  for (Index i = 0; i < m0; ++i) {
    const double rl = model.row_lower[i];
    const double ru = model.row_upper[i];
    if (rl == -kInf && ru == kInf) continue;  // A free row constrains nothing.
    row_pos[i] = f.m++;
    if (rl == ru) {
      f.b.push_back(rl - shift[i]);
    } else {
      f.b.push_back(0.0);
      slack_of[i] = f.n++;
    }
  }
  f.c.assign(static_cast<std::size_t>(f.n), 0.0);
  f.l.assign(static_cast<std::size_t>(f.n), 0.0);
  f.u.assign(static_cast<std::size_t>(f.n), 0.0);
  for (Index j = 0; j < n0; ++j) {
    if (col_pos[j] < 0) continue;
    const auto k = static_cast<std::size_t>(col_pos[j]);
    f.c[k] = sense * model.obj[j];
    f.l[k] = model.col_lower[j];
    f.u[k] = model.col_upper[j];
    for (NnzIndex p = start[j]; p < start[j + 1]; ++p) {
      if (row_pos[index[p]] >= 0) f.a.push_back({row_pos[index[p]], col_pos[j], value[p]});
    }
  }
  for (Index i = 0; i < m0; ++i) {
    if (slack_of[i] < 0) continue;
    const auto s = static_cast<std::size_t>(slack_of[i]);
    f.a.push_back({row_pos[i], slack_of[i], -1.0});
    f.l[s] = model.row_lower[i] - shift[i];
    f.u[s] = model.row_upper[i] - shift[i];
  }
  // Q: lower triangle in the model; both triangles here. Entries with a fixed column become
  // linear terms (q_jk x_k with x_k fixed) or drop into the constant.
  double q_diag_max = 0.0;
  {
    const auto qs = model.Q.col_start();
    const auto qi = model.Q.row_index();
    const auto qv = model.Q.values();
    for (Index j = 0; j < model.Q.cols(); ++j) {
      for (NnzIndex p = qs[j]; p < qs[j + 1]; ++p) {
        const Index i = qi[p];
        const double v = sense * qv[p];
        const Index pi = col_pos[i];
        const Index pj = col_pos[j];
        if (pi >= 0 && pj >= 0) {
          f.q.push_back({pi, pj, v});
          if (pi != pj) f.q.push_back({pj, pi, v});
          if (pi == pj) q_diag_max = std::max(q_diag_max, std::fabs(v));
        } else if (pi >= 0) {
          f.c[static_cast<std::size_t>(pi)] += v * model.col_lower[j];
        } else if (pj >= 0) {
          f.c[static_cast<std::size_t>(pj)] += v * model.col_lower[i];
        }
      }
    }
  }
  const bool has_q = !f.q.empty();
  const Index n = f.n;
  const Index m = f.m;
  const auto nn = static_cast<std::size_t>(n);
  const auto mm = static_cast<std::size_t>(m);

  // ---- Convexity: Q + tol * max|q_jj| I must have an LDL^T with positive pivots (the inertia of
  // an LDL^T is that of the matrix, Sylvester's law).
  if (has_q) {
    std::vector<Triplet> t;
    for (const Triplet& e : f.q) {
      if (e.row <= e.col) t.push_back(e);
    }
    const double shift_q = options.convexity_tol * std::max(q_diag_max, 1.0);
    for (Index j = 0; j < kept; ++j) t.push_back({j, j, shift_q});
    const SparseMatrix upper = SparseMatrix::from_triplets(kept, kept, std::move(t));
    LdlFactor qf;
    qf.analyze(upper);
    const std::vector<int> positive(static_cast<std::size_t>(kept), 1);
    if (qf.factorize(upper, positive, 0.0, 1.0) > 0) {
      result.status = IpmStatus::kNotConvex;
      log.log(1, "IPM: Q is not positive semidefinite");
      return result;
    }
  }

  // ---- Ruiz equilibration of [Q A'; A 0]: column scales dc (n), row scales dr (m).
  std::vector<double> dc(nn, 1.0);
  std::vector<double> dr(mm, 1.0);
  for (int pass = 0; pass < kRuizPasses; ++pass) {
    std::vector<double> col_norm(nn, 0.0);
    std::vector<double> row_norm(mm, 0.0);
    for (const Triplet& e : f.a) {
      const auto i = static_cast<std::size_t>(e.row);
      const auto j = static_cast<std::size_t>(e.col);
      const double v = std::fabs(e.value * dr[i] * dc[j]);
      col_norm[j] = std::max(col_norm[j], v);
      row_norm[i] = std::max(row_norm[i], v);
    }
    for (const Triplet& e : f.q) {
      const auto i = static_cast<std::size_t>(e.row);
      const auto j = static_cast<std::size_t>(e.col);
      col_norm[j] = std::max(col_norm[j], std::fabs(e.value * dc[i] * dc[j]));
    }
    for (std::size_t j = 0; j < nn; ++j) {
      if (col_norm[j] > 0.0) dc[j] /= std::sqrt(col_norm[j]);
    }
    for (std::size_t i = 0; i < mm; ++i) {
      if (row_norm[i] > 0.0) dr[i] /= std::sqrt(row_norm[i]);
    }
  }
  for (Triplet& e : f.a) {
    e.value *= dr[static_cast<std::size_t>(e.row)] * dc[static_cast<std::size_t>(e.col)];
  }
  for (Triplet& e : f.q) {
    e.value *= dc[static_cast<std::size_t>(e.row)] * dc[static_cast<std::size_t>(e.col)];
  }
  for (std::size_t j = 0; j < nn; ++j) {
    f.c[j] *= dc[j];
    f.l[j] /= dc[j];
    f.u[j] /= dc[j];
  }
  for (std::size_t i = 0; i < mm; ++i) f.b[i] *= dr[i];
  // Objective scaling: costs and Q to a largest entry of about 1.
  double cost_scale = inf_norm(f.c);
  for (const Triplet& e : f.q) cost_scale = std::max(cost_scale, std::fabs(e.value));
  cost_scale = cost_scale > 0.0 ? 1.0 / cost_scale : 1.0;
  for (double& v : f.c) v *= cost_scale;
  for (Triplet& e : f.q) e.value *= cost_scale;

  const SparseMatrix A = SparseMatrix::from_triplets(m, n, f.a);
  const SparseMatrix Q = SparseMatrix::from_triplets(n, n, f.q);
  std::vector<double> q_diag(nn, 0.0);
  for (const Triplet& e : f.q) {
    if (e.row == e.col) q_diag[static_cast<std::size_t>(e.row)] += e.value;
  }
  std::vector<char> has_l(nn), has_u(nn);
  int bounded = 0;
  for (std::size_t j = 0; j < nn; ++j) {
    has_l[j] = std::isfinite(f.l[j]) ? 1 : 0;
    has_u[j] = std::isfinite(f.u[j]) ? 1 : 0;
    bounded += has_l[j] + has_u[j];
  }
  const double b_norm = inf_norm(f.b);
  const double c_norm = inf_norm(f.c);

  // ---- KKT pattern: upper triangle of [-(Q + Theta^-1 + rho I)  A'; A  delta I].
  std::vector<Triplet> kkt_fixed;  // Everything but the primal diagonal.
  for (const Triplet& e : f.q) {
    if (e.row < e.col) kkt_fixed.push_back({e.row, e.col, -e.value});
  }
  for (const Triplet& e : f.a) kkt_fixed.push_back({e.col, n + e.row, e.value});
  for (Index i = 0; i < m; ++i) kkt_fixed.push_back({n + i, n + i, kDualRegularization});
  std::vector<int> sign(nn + mm, 1);
  std::fill(sign.begin(), sign.begin() + static_cast<std::ptrdiff_t>(nn), -1);
  const auto kkt_matrix = [&](const std::vector<double>& theta_inv) {
    std::vector<Triplet> t = kkt_fixed;
    for (Index j = 0; j < n; ++j) {
      const auto uj = static_cast<std::size_t>(j);
      t.push_back({j, j, -(q_diag[uj] + theta_inv[uj] + kPrimalRegularization)});
    }
    return SparseMatrix::from_triplets(n + m, n + m, std::move(t));
  };
  LdlFactor ldl;

  std::vector<double> z(nn, 0.0), y(mm, 0.0), zl(nn, 0.0), zu(nn, 0.0);
  // Bound slacks z - l and u - z, kept as their own vectors: recomputing them by subtraction
  // loses them to round-off once they fall below eps |z| (then Theta^-1 divides by zero).
  std::vector<double> sl(nn, 0.0), su(nn, 0.0);
  std::vector<double> rp(mm), rd(nn), qz(nn), aty(nn), theta_inv(nn, 1.0);
  std::vector<double> dz(nn), dy(mm), dzl(nn), dzu(nn);
  std::vector<double> dz_aff(nn), dzl_aff(nn), dzu_aff(nn);
  std::vector<double> rl_c(nn), ru_c(nn);
  std::vector<double> rhs(nn + mm), sol(nn + mm), res(nn + mm), tmp_n(nn), tmp_m(mm);

  // Solves the Newton system for (dz, dy) given the complementarity right-hand sides.
  const auto newton = [&]() {
    for (std::size_t j = 0; j < nn; ++j) {
      double r = rd[j];
      if (has_l[j]) r -= rl_c[j] / sl[j];
      if (has_u[j]) r += ru_c[j] / su[j];
      rhs[j] = r;
    }
    for (std::size_t i = 0; i < mm; ++i) rhs[nn + i] = rp[i];
    sol = rhs;
    ldl.solve(sol);
    // Iterative refinement against the unregularized matrix.
    for (int step = 0; step < kRefinementSteps; ++step) {
      std::vector<double> x1(sol.begin(), sol.begin() + static_cast<std::ptrdiff_t>(nn));
      std::vector<double> x2(sol.begin() + static_cast<std::ptrdiff_t>(nn), sol.end());
      Q.multiply(x1, tmp_n);
      std::vector<double> atx2(nn, 0.0);
      A.multiply_transpose(x2, atx2);
      A.multiply(x1, tmp_m);
      double worst = 0.0;
      for (std::size_t j = 0; j < nn; ++j) {
        res[j] = rhs[j] - (-(tmp_n[j] + theta_inv[j] * x1[j]) + atx2[j]);
        worst = std::max(worst, std::fabs(res[j]));
      }
      for (std::size_t i = 0; i < mm; ++i) {
        res[nn + i] = rhs[nn + i] - tmp_m[i];
        worst = std::max(worst, std::fabs(res[nn + i]));
      }
      if (worst <= 1e-14 * (1.0 + inf_norm(rhs))) break;
      ldl.solve(res);
      for (std::size_t k = 0; k < nn + mm; ++k) sol[k] += res[k];
    }
    if (log.enabled(3)) {
      std::vector<double> x1(sol.begin(), sol.begin() + static_cast<std::ptrdiff_t>(nn));
      std::vector<double> x2(sol.begin() + static_cast<std::ptrdiff_t>(nn), sol.end());
      Q.multiply(x1, tmp_n);
      std::vector<double> atx2(nn, 0.0);
      A.multiply_transpose(x2, atx2);
      A.multiply(x1, tmp_m);
      double worst = 0.0;
      for (std::size_t j = 0; j < nn; ++j) {
        const double lhs = -(tmp_n[j] + theta_inv[j] * x1[j]) + atx2[j];
        worst = std::max(worst, std::fabs(rhs[j] - lhs));
      }
      for (std::size_t i = 0; i < mm; ++i) {
        worst = std::max(worst, std::fabs(rhs[nn + i] - tmp_m[i]));
      }
      log.log(3, "IPM   newton residual %.2e (rhs %.2e, sol %.2e)", worst, inf_norm(rhs),
              inf_norm(sol));
    }
    for (std::size_t j = 0; j < nn; ++j) {
      dz[j] = sol[j];
      dzl[j] = has_l[j] ? (rl_c[j] - zl[j] * dz[j]) / sl[j] : 0.0;
      dzu[j] = has_u[j] ? (ru_c[j] + zu[j] * dz[j]) / su[j] : 0.0;
    }
    for (std::size_t i = 0; i < mm; ++i) dy[i] = sol[nn + i];
  };
  const auto max_steps = [&](double& alpha_p, double& alpha_d) {
    alpha_p = 1.0;
    alpha_d = 1.0;
    for (std::size_t j = 0; j < nn; ++j) {
      if (has_l[j] && dz[j] < 0.0) alpha_p = std::min(alpha_p, -sl[j] / dz[j]);
      if (has_u[j] && dz[j] > 0.0) alpha_p = std::min(alpha_p, su[j] / dz[j]);
      if (has_l[j] && dzl[j] < 0.0) alpha_d = std::min(alpha_d, -zl[j] / dzl[j]);
      if (has_u[j] && dzu[j] < 0.0) alpha_d = std::min(alpha_d, -zu[j] / dzu[j]);
    }
  };

  // ---- Starting point (Mehrotra, SIAM J. Optim. 1992, adapted to bounds): the least-norm z with
  // A z = b and the least-squares y for A'y ~ c + Qz, from the KKT system with Theta^-1 = I;
  // then the bound slacks and duals are shifted into the interior, first so the most negative
  // becomes positive, then so the complementarity products are balanced.
  {
    ldl.analyze(kkt_matrix(theta_inv));
    result.regularized_pivots += ldl.factorize(kkt_matrix(theta_inv), sign, kPivotTol,
                                               kPivotReplacement);
    std::vector<double> v(nn + mm, 0.0);
    for (std::size_t i = 0; i < mm; ++i) v[nn + i] = f.b[i];
    ldl.solve(v);  // [-(Q+I) A'; A 0][z; y] = [0; b].
    for (std::size_t j = 0; j < nn; ++j) z[j] = v[j];
    std::fill(v.begin(), v.end(), 0.0);
    for (std::size_t j = 0; j < nn; ++j) v[j] = f.c[j];
    ldl.solve(v);  // [-(Q+I) A'; A 0][r; y] = [c; 0]: A'y approximates c.
    for (std::size_t i = 0; i < mm; ++i) y[i] = v[nn + i];
    Q.multiply(z, qz);
    std::fill(aty.begin(), aty.end(), 0.0);
    A.multiply_transpose(y, aty);
    // Bound slacks of z and the reduced costs split into bound duals.
    double min_slack = kInf;
    double min_dual = kInf;
    for (std::size_t j = 0; j < nn; ++j) {
      const double d = f.c[j] + qz[j] - aty[j];
      if (has_l[j] && has_u[j]) {
        zl[j] = std::max(d, 0.0);
        zu[j] = std::max(-d, 0.0);
      } else if (has_l[j]) {
        zl[j] = d;
      } else if (has_u[j]) {
        zu[j] = -d;
      }
      if (has_l[j]) {
        min_slack = std::min(min_slack, z[j] - f.l[j]);
        min_dual = std::min(min_dual, zl[j]);
      }
      if (has_u[j]) {
        min_slack = std::min(min_slack, f.u[j] - z[j]);
        min_dual = std::min(min_dual, zu[j]);
      }
    }
    if (bounded > 0) {
      const double dp = std::max(-1.5 * min_slack, 0.0);
      const double dd = std::max(-1.5 * min_dual, 0.0);
      double prod = 0.0;
      double sum_slack = 0.0;
      double sum_dual = 0.0;
      for (std::size_t j = 0; j < nn; ++j) {
        if (has_l[j]) {
          const double slack_l = z[j] - f.l[j] + dp;
          prod += slack_l * (zl[j] + dd);
          sum_slack += slack_l;
          sum_dual += zl[j] + dd;
        }
        if (has_u[j]) {
          const double slack_u = f.u[j] - z[j] + dp;
          prod += slack_u * (zu[j] + dd);
          sum_slack += slack_u;
          sum_dual += zu[j] + dd;
        }
      }
      // Keep both shifts positive even when the least-squares point is already interior.
      const double shift_p = dp + (sum_dual > 0.0 ? 0.5 * prod / sum_dual : 0.0) + 1e-2;
      const double shift_d = dd + (sum_slack > 0.0 ? 0.5 * prod / sum_slack : 0.0) + 1e-2;
      for (std::size_t j = 0; j < nn; ++j) {
        if (has_l[j] && has_u[j]) {
          // A box: move inside by the shift, at most to the midpoint.
          const double margin = std::min(shift_p, 0.5 * (f.u[j] - f.l[j]));
          z[j] = std::clamp(z[j], f.l[j] + margin, f.u[j] - margin);
        } else if (has_l[j]) {
          z[j] = std::max(z[j], f.l[j]) + shift_p;
        } else if (has_u[j]) {
          z[j] = std::min(z[j], f.u[j]) - shift_p;
        }
        if (has_l[j]) zl[j] = std::max(zl[j], 0.0) + shift_d;
        if (has_u[j]) zu[j] = std::max(zu[j], 0.0) + shift_d;
      }
    }
  }
  for (std::size_t j = 0; j < nn; ++j) {
    if (has_l[j]) sl[j] = z[j] - f.l[j];
    if (has_u[j]) su[j] = f.u[j] - z[j];
  }
  const double z0_norm = inf_norm(z);

  result.status = IpmStatus::kIterationLimit;
  // The best iterate by max(primal residual, dual residual, gap).
  double best_merit = kInf;
  int best_iter = 0;
  double stall_merit = kInf;
  int stall_start = 0;
  std::vector<double> best_z, best_y, best_zl, best_zu, best_sl, best_su;
  double last_step = 1.0;
  for (int iter = 0; iter <= options.max_iterations; ++iter) {
    // Residuals and complementarity.
    A.multiply(z, tmp_m);
    for (std::size_t i = 0; i < mm; ++i) rp[i] = f.b[i] - tmp_m[i];
    Q.multiply(z, qz);
    std::fill(aty.begin(), aty.end(), 0.0);
    A.multiply_transpose(y, aty);
    double comp = 0.0;
    double lz = 0.0;  // l'zl - u'zu, finite bounds only.
    for (std::size_t j = 0; j < nn; ++j) {
      rd[j] = f.c[j] + qz[j] - aty[j] - zl[j] + zu[j];
      if (has_l[j]) {
        comp += sl[j] * zl[j];
        lz += f.l[j] * zl[j];
      }
      if (has_u[j]) {
        comp += su[j] * zu[j];
        lz -= f.u[j] * zu[j];
      }
    }
    const double mu = bounded > 0 ? comp / bounded : 0.0;
    const double zqz = dot(z, qz);
    const double pobj = dot(f.c, z) + 0.5 * zqz;
    const double dobj = dot(f.b, y) - 0.5 * zqz + lz;
    result.primal_residual = inf_norm(rp) / (1.0 + b_norm);
    result.dual_residual = inf_norm(rd) / (1.0 + c_norm);
    result.gap = std::fabs(pobj - dobj) / (1.0 + std::fabs(pobj));
    result.iterations = iter;
    log.log(2, "IPM %3d  pobj %+.10e  dobj %+.10e  pres %.2e  dres %.2e  gap %.2e  mu %.2e", iter,
            pobj / cost_scale, dobj / cost_scale, result.primal_residual, result.dual_residual,
            result.gap, mu);
    if (result.primal_residual <= options.tol && result.dual_residual <= options.tol &&
        result.gap <= options.tol) {
      result.status = IpmStatus::kOptimal;
      break;
    }
    const double merit = std::max({result.primal_residual, result.dual_residual, result.gap});
    if (!std::isfinite(merit) || !std::isfinite(mu)) {
      log.log(1, "IPM: iterate is not finite");
      result.status = IpmStatus::kNumericalError;
      if (best_merit <= options.acceptable_tol) {
        z = best_z;
        y = best_y;
        zl = best_zl;
        zu = best_zu;
        sl = best_sl;
        su = best_su;
        result.status = IpmStatus::kOptimal;
      }
      break;
    }
    if (merit < best_merit) {
      best_merit = merit;
      best_iter = iter;
      best_z = z;
      best_y = y;
      best_zl = zl;
      best_zu = zu;
      best_sl = sl;
      best_su = su;
    }
    if (merit < 0.1 * stall_merit) {
      stall_merit = merit;
      stall_start = iter;
    }
    const bool stalled = iter - stall_start >= options.stall_iterations;
    if (stalled || iter == options.max_iterations) {
      if (best_merit <= options.acceptable_tol) {
        // Stalled short of the target but within the acceptable tolerance: return the best
        // iterate; the caller's verifier decides whether it is optimal.
        z = best_z;
        y = best_y;
        zl = best_zl;
        zu = best_zu;
        sl = best_sl;
        su = best_su;
        result.status = IpmStatus::kOptimal;
        log.log(1, "IPM: stalled at %.1e; using iterate %d (merit %.1e)", merit, best_iter,
                best_merit);
      }
      break;
    }
    if (timer.seconds() > options.time_limit) {
      result.status = IpmStatus::kTimeLimit;
      break;
    }
    if (inf_norm(z) > kDivergence * (1.0 + z0_norm) || inf_norm(y) > kDivergence) {
      log.log(1, "IPM: iterates diverge (infeasible or unbounded model)");
      break;
    }

    // Factorize with the current Theta^-1.
    for (std::size_t j = 0; j < nn; ++j) {
      double t = 0.0;
      if (has_l[j]) t += zl[j] / sl[j];
      if (has_u[j]) t += zu[j] / su[j];
      theta_inv[j] = t;
    }
    const SparseMatrix kkt = kkt_matrix(theta_inv);
    result.regularized_pivots += ldl.factorize(kkt, sign, kPivotTol, kPivotReplacement);

    // Predictor (affine scaling): complementarity driven to zero.
    for (std::size_t j = 0; j < nn; ++j) {
      rl_c[j] = has_l[j] ? -sl[j] * zl[j] : 0.0;
      ru_c[j] = has_u[j] ? -su[j] * zu[j] : 0.0;
    }
    newton();
    double ap = 0.0;
    double ad = 0.0;
    max_steps(ap, ad);
    if (has_q) ap = ad = std::min(ap, ad);
    double comp_aff = 0.0;
    for (std::size_t j = 0; j < nn; ++j) {
      if (has_l[j]) comp_aff += (sl[j] + ap * dz[j]) * (zl[j] + ad * dzl[j]);
      if (has_u[j]) comp_aff += (su[j] - ap * dz[j]) * (zu[j] + ad * dzu[j]);
    }
    const double mu_aff = bounded > 0 ? comp_aff / bounded : 0.0;
    // Mehrotra's centering, but at least kShortStepSigma after a short step: iterates pressed
    // against the boundary need centering before they can move (otherwise complementarity
    // collapses to zero while the residuals stall).
    double sigma = mu > 0.0 ? std::pow(mu_aff / mu, 3.0) : 0.0;
    if (last_step < kShortStep) sigma = std::max(sigma, kShortStepSigma);
    dz_aff = dz;
    dzl_aff = dzl;
    dzu_aff = dzu;

    // Corrector: centering plus the second-order term of the predictor.
    for (std::size_t j = 0; j < nn; ++j) {
      rl_c[j] = has_l[j] ? sigma * mu - sl[j] * zl[j] - dz_aff[j] * dzl_aff[j] : 0.0;
      ru_c[j] = has_u[j] ? sigma * mu - su[j] * zu[j] + dz_aff[j] * dzu_aff[j] : 0.0;
    }
    newton();
    max_steps(ap, ad);
    if (has_q) ap = ad = std::min(ap, ad);

    ap = std::min(1.0, kStepFraction * ap);
    ad = std::min(1.0, kStepFraction * ad);
    for (std::size_t j = 0; j < nn; ++j) {
      z[j] += ap * dz[j];
      sl[j] += ap * dz[j];
      su[j] -= ap * dz[j];
      zl[j] += ad * dzl[j];
      zu[j] += ad * dzu[j];
    }
    for (std::size_t i = 0; i < mm; ++i) y[i] += ad * dy[i];
    last_step = std::min(ap, ad);
    if (log.enabled(3)) {
      double min_sl = kInf;
      double min_zl = kInf;
      for (std::size_t j = 0; j < nn; ++j) {
        if (has_l[j]) min_sl = std::min(min_sl, sl[j]);
        if (has_u[j]) min_sl = std::min(min_sl, su[j]);
        if (has_l[j]) min_zl = std::min(min_zl, zl[j]);
        if (has_u[j]) min_zl = std::min(min_zl, zu[j]);
      }
      log.log(3, "IPM   steps %.2e %.2e sigma %.2e |dz| %.2e |dy| %.2e min slack %.2e "
              "min dual %.2e", ap, ad, sigma, inf_norm(dz), inf_norm(dy), min_sl, min_zl);
    }
  }

  // ---- Polishing, a primal-dual active-set loop (Hintermueller, Ito and Kunisch, SIAM J. Optim.
  // 2002) started from the bounds the iterate treats as active (slack below its dual). Each round
  // solves  -(Q + W) z + A'y = c - W bound,  A z = b  with W = kPolishWeight on the pinned
  // columns and 0 elsewhere (so the free columns get zero reduced costs), then pins free columns
  // that left their bounds and releases pinned ones whose reduced cost has the wrong sign, until
  // nothing changes. Interior-point iterates stay a little off their bounds, which a strict
  // optimality check rejects; the caller verifies both this point and the iterate.
  std::vector<double> pz;
  std::vector<double> py;
  if (result.status == IpmStatus::kOptimal) {
    std::vector<double> weight(nn, 0.0);
    std::vector<double> target(nn, 0.0);
    std::vector<signed char> active(nn, 0);  // -1 at lower, +1 at upper.
    const auto pin = [&](std::size_t j, signed char side) {
      active[j] = side;
      target[j] = side < 0 ? f.l[j] : f.u[j];
      weight[j] = kPolishWeight;
    };
    for (std::size_t j = 0; j < nn; ++j) {
      const double s_low = has_l[j] ? sl[j] : kInf;
      const double s_up = has_u[j] ? su[j] : kInf;
      if (has_l[j] && s_low < zl[j] && s_low <= s_up) {
        pin(j, -1);
      } else if (has_u[j] && s_up < zu[j]) {
        pin(j, 1);
      }
    }
    std::vector<double> r0(nn + mm);
    std::vector<double> v;
    bool settled = false;
    for (int round = 0; round < kPolishRounds && !settled; ++round) {
      ldl.factorize(kkt_matrix(weight), sign, kPivotTol, kPivotReplacement);
      for (std::size_t k = 0; k < nn; ++k) r0[k] = f.c[k] - weight[k] * target[k];
      for (std::size_t i = 0; i < mm; ++i) r0[nn + i] = f.b[i];
      v = r0;
      ldl.solve(v);
      for (int step = 0; step < kRefinementSteps; ++step) {
        std::vector<double> x1(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(nn));
        std::vector<double> x2(v.begin() + static_cast<std::ptrdiff_t>(nn), v.end());
        Q.multiply(x1, tmp_n);
        std::vector<double> atx2(nn, 0.0);
        A.multiply_transpose(x2, atx2);
        A.multiply(x1, tmp_m);
        for (std::size_t k = 0; k < nn; ++k) {
          res[k] = r0[k] - (-(tmp_n[k] + weight[k] * x1[k]) + atx2[k]);
        }
        for (std::size_t i = 0; i < mm; ++i) res[nn + i] = r0[nn + i] - tmp_m[i];
        ldl.solve(res);
        for (std::size_t k = 0; k < nn + mm; ++k) v[k] += res[k];
      }
      // Reduced costs d = c + Qz - A'y of the round's point.
      std::vector<double> x1(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(nn));
      std::vector<double> x2(v.begin() + static_cast<std::ptrdiff_t>(nn), v.end());
      Q.multiply(x1, tmp_n);
      std::vector<double> atx2(nn, 0.0);
      A.multiply_transpose(x2, atx2);
      settled = true;
      for (std::size_t j = 0; j < nn; ++j) {
        if (active[j] == 0) {
          if (has_l[j] && x1[j] < f.l[j] - kPolishTol * (1.0 + std::fabs(f.l[j]))) {
            pin(j, -1);
            settled = false;
          } else if (has_u[j] && x1[j] > f.u[j] + kPolishTol * (1.0 + std::fabs(f.u[j]))) {
            pin(j, 1);
            settled = false;
          }
        } else {
          const double d = f.c[j] + tmp_n[j] - atx2[j];
          if ((active[j] < 0 && d < -kPolishTol) || (active[j] > 0 && d > kPolishTol)) {
            active[j] = 0;
            weight[j] = 0.0;
            settled = false;
          }
        }
      }
    }
    log.log(2, "IPM: polishing %s", settled ? "settled" : "did not settle");
    pz.assign(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(nn));
    py.assign(v.begin() + static_cast<std::ptrdiff_t>(nn), v.end());
    for (std::size_t j = 0; j < nn; ++j) {
      if (active[j] != 0) pz[j] = target[j];
      if (has_l[j]) pz[j] = std::max(pz[j], f.l[j]);
      if (has_u[j]) pz[j] = std::min(pz[j], f.u[j]);
    }
  }

  // ---- Back to the model: unscale, restore fixed columns, map the row duals.
  result.x.assign(static_cast<std::size_t>(n0), 0.0);
  for (Index j = 0; j < n0; ++j) {
    const Index k = col_pos[j];
    result.x[j] = k >= 0 ? z[static_cast<std::size_t>(k)] * dc[static_cast<std::size_t>(k)]
                         : model.col_lower[j];
  }
  result.y.assign(static_cast<std::size_t>(m0), 0.0);
  for (Index i = 0; i < m0; ++i) {
    const Index r = row_pos[i];
    if (r < 0) continue;
    result.y[i] = sense * y[static_cast<std::size_t>(r)] * dr[static_cast<std::size_t>(r)] /
                  cost_scale;
  }
  if (!pz.empty()) {
    result.polished_x.assign(static_cast<std::size_t>(n0), 0.0);
    result.polished_y.assign(static_cast<std::size_t>(m0), 0.0);
    for (Index j = 0; j < n0; ++j) {
      const Index k = col_pos[j];
      result.polished_x[j] = k >= 0 ? pz[static_cast<std::size_t>(k)] *
                                          dc[static_cast<std::size_t>(k)]
                                    : model.col_lower[j];
    }
    for (Index i = 0; i < m0; ++i) {
      const Index r = row_pos[i];
      if (r < 0) continue;
      result.polished_y[i] = sense * py[static_cast<std::size_t>(r)] *
                             dr[static_cast<std::size_t>(r)] / cost_scale;
    }
  }
  log.log(1, "IPM: %s after %d iterations, %.2f s (%d regularized pivots, factor %lld nonzeros)",
          to_string(result.status), result.iterations, timer.seconds(), result.regularized_pivots,
          static_cast<long long>(ldl.factor_nnz()));
  return result;
}

}  // namespace samaya
