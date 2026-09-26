#include "mip/cuts.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <set>

namespace samaya {

namespace {

constexpr double kMinFractionality = 0.01;  // Gomory: basic value at least this far from integer.
constexpr double kTinyCoefficient = 1e-12;
constexpr double kMaxDynamism = 1e6;    // Largest ratio of cut coefficient magnitudes.
constexpr double kRhsRelaxation = 1e-9;  // Relative safety margin on the right-hand side.
constexpr double kMinViolation = 1e-6;
constexpr int kMaxMirDivisors = 8;
constexpr NnzIndex kMaxRowLength = 1000;  // MIR and cover rows longer than this are skipped.

bool is_integer_col(const Model& model, Index j) {
  return model.col_type[j] == VarType::kInteger;
}

// Accumulates a cut as dense column coefficients plus a right-hand side.
class CutBuilder {
 public:
  explicit CutBuilder(const CutContext& ctx)
      : ctx_(ctx),
        dense_(static_cast<std::size_t>(ctx.model.num_cols()), 0.0),
        mark_(static_cast<std::size_t>(ctx.model.num_cols()), 0) {}

  void add(Index j, double v) {
    if (!mark_[j]) {
      mark_[j] = 1;
      touched_.push_back(j);
    }
    dense_[j] += v;
  }
  // Adds v times row i's activity.
  void add_row(Index i, double v) {
    const auto start = ctx_.At.col_start();
    const auto index = ctx_.At.row_index();
    const auto value = ctx_.At.values();
    for (NnzIndex p = start[i]; p < start[i + 1]; ++p) add(index[p], v * value[p]);
  }
  // Returns  sum dense_j x_j  (sign) rhs  as a >= cut: sign +1 keeps it, -1 negates a <= cut.
  Cut take(double rhs, double sign) {
    Cut cut;
    for (const Index j : touched_) {
      if (dense_[j] != 0.0) {
        cut.index.push_back(j);
        cut.value.push_back(sign * dense_[j]);
      }
      dense_[j] = 0.0;
      mark_[j] = 0;
    }
    touched_.clear();
    cut.lower = sign * rhs;
    return cut;
  }
  void clear() { take(0.0, 1.0); }

 private:
  const CutContext& ctx_;
  std::vector<double> dense_;
  std::vector<char> mark_;
  std::vector<Index> touched_;
};

double violation_of(const CutContext& ctx, const Cut& cut) {
  double activity = 0.0;
  for (std::size_t k = 0; k < cut.index.size(); ++k) activity += cut.value[k] * ctx.x[cut.index[k]];
  return cut.lower - activity;
}

double norm_of(const Cut& cut) {
  double s = 0.0;
  for (const double v : cut.value) s += v * v;
  return std::sqrt(s);
}

}  // namespace

bool finalize_cut(const CutContext& ctx, Cut& cut) {
  if (cut.index.empty() || !std::isfinite(cut.lower)) return false;
  double max_abs = 0.0;
  for (const double v : cut.value) {
    if (!std::isfinite(v)) return false;
    max_abs = std::max(max_abs, std::fabs(v));
  }
  if (max_abs == 0.0) return false;
  // Drop tiny coefficients: pi_j x_j <= max(pi_j l_j, pi_j u_j), so the rest must cover the
  // right-hand side minus that amount.
  std::size_t kept = 0;
  double min_abs = kInf;
  for (std::size_t k = 0; k < cut.index.size(); ++k) {
    const Index j = cut.index[k];
    const double v = cut.value[k];
    if (std::fabs(v) < kTinyCoefficient * max_abs) {
      const double largest = std::max(v * ctx.lower[j], v * ctx.upper[j]);
      if (!std::isfinite(largest)) return false;
      cut.lower -= largest;
      continue;
    }
    min_abs = std::min(min_abs, std::fabs(v));
    cut.index[kept] = j;
    cut.value[kept] = v;
    ++kept;
  }
  cut.index.resize(kept);
  cut.value.resize(kept);
  if (kept == 0 || max_abs / min_abs > kMaxDynamism) return false;
  cut.lower -= kRhsRelaxation * std::max(1.0, std::fabs(cut.lower));
  const double violation = violation_of(ctx, cut);
  if (violation <= kMinViolation * (1.0 + std::fabs(cut.lower))) return false;
  cut.efficacy = violation / norm_of(cut);
  return true;
}

bool gomory_mixed_integer_cut(const CutContext& ctx, Index k, const std::vector<double>& row,
                              const std::vector<VarStatus>& status, Cut& cut) {
  const Model& model = ctx.model;
  const Index n = model.num_cols();
  const Index m = model.num_rows();
  // Tableau row: x_k + sum_N a_j v_j = 0. With v_j = bound_j +- t_j (t_j >= 0) it becomes
  // x_k + sum a'_j t_j = b'.
  struct Term {
    Index var;
    double a;      // a'_j
    double bound;  // The bound t_j is measured from.
    bool at_upper;
  };
  std::vector<Term> terms;
  double b = 0.0;
  double row_max = 0.0;
  for (Index j = 0; j < n + m; ++j) {
    if (j != k) row_max = std::max(row_max, std::fabs(row[j]));
  }
  for (Index j = 0; j < n + m; ++j) {
    if (j == k || status[j] == VarStatus::kBasic) continue;
    const double a = row[j];
    if (std::fabs(a) <= kTinyCoefficient * std::max(1.0, row_max)) continue;
    const double lo = j < n ? ctx.lower[j] : model.row_lower[j - n];
    const double up = j < n ? ctx.upper[j] : model.row_upper[j - n];
    if (status[j] == VarStatus::kAtZero) return false;  // Free nonbasic: no t_j >= 0 form.
    const bool at_upper = status[j] == VarStatus::kAtUpper;
    const double bound = at_upper ? up : lo;
    if (!std::isfinite(bound)) return false;
    b -= a * bound;
    if (lo == up) continue;  // Fixed: t_j = 0.
    terms.push_back({j, at_upper ? -a : a, bound, at_upper});
  }
  const double xk = ctx.x[k];
  if (std::fabs(b - xk) > 1e-6 * (1.0 + std::fabs(xk))) return false;  // Inaccurate row.
  const double f0 = b - std::floor(b);
  if (f0 < kMinFractionality || f0 > 1.0 - kMinFractionality) return false;

  // GMI:  sum_int min(f_j/f0, (1-f_j)/(1-f0)) t_j
  //     + sum_cont (a_j > 0 ? a_j/f0 : -a_j/(1-f0)) t_j >= 1.
  CutBuilder builder(ctx);
  double rhs = 1.0;
  for (const Term& t : terms) {
    const bool integral = t.var < n && is_integer_col(model, t.var) &&
                          t.bound == std::round(t.bound);
    double coef;
    if (integral) {
      const double fj = t.a - std::floor(t.a);
      coef = fj <= f0 ? fj / f0 : (1.0 - fj) / (1.0 - f0);
    } else {
      coef = t.a >= 0.0 ? t.a / f0 : -t.a / (1.0 - f0);
    }
    if (coef == 0.0) continue;
    // t = v - bound (at lower) or bound - v (at upper).
    const double v_coef = t.at_upper ? -coef : coef;
    rhs += v_coef * t.bound;
    if (t.var < n) {
      builder.add(t.var, v_coef);
    } else {
      builder.add_row(t.var - n, v_coef);
    }
  }
  cut = builder.take(rhs, 1.0);
  return true;
}

void separate_mir(const CutContext& ctx, std::vector<Cut>& cuts) {
  const Model& model = ctx.model;
  const auto start = ctx.At.col_start();
  const auto index = ctx.At.row_index();
  const auto value = ctx.At.values();
  CutBuilder builder(ctx);
  struct IntTerm {
    Index col;
    double c;  // Coefficient of t_j >= 0.
    double t;  // t_j at the LP point.
    double bound;
    bool at_upper;
  };
  struct ContTerm {
    Index col;
    double h;
    double bound;
    bool at_upper;
  };
  std::vector<IntTerm> ints;
  std::vector<ContTerm> conts;
  for (Index i = 0; i < ctx.original_rows; ++i) {
    if (start[i + 1] - start[i] > kMaxRowLength || start[i + 1] == start[i]) continue;
    for (int side = 0; side < 2; ++side) {
      // side 0: a x <= ru;  side 1: -a x <= -rl.
      const double bound = side == 0 ? model.row_upper[i] : model.row_lower[i];
      if (!std::isfinite(bound)) continue;
      const double sign = side == 0 ? 1.0 : -1.0;
      double beta = sign * bound;
      ints.clear();
      conts.clear();
      bool ok = true;
      for (NnzIndex p = start[i]; p < start[i + 1] && ok; ++p) {
        const Index j = index[p];
        const double c = sign * value[p];
        const double lo = ctx.lower[j];
        const double up = ctx.upper[j];
        const double xj = ctx.x[j];
        const bool use_lower =
            std::isfinite(lo) && (!std::isfinite(up) || xj - lo <= up - xj);
        if (!use_lower && !std::isfinite(up)) {
          ok = false;
          break;
        }
        const double b = use_lower ? lo : up;
        beta -= c * b;
        if (lo == up) continue;
        const double coef = use_lower ? c : -c;
        if (is_integer_col(model, j)) {
          ints.push_back({j, coef, use_lower ? xj - lo : up - xj, b, !use_lower});
        } else {
          conts.push_back({j, coef, b, !use_lower});
        }
      }
      if (!ok || ints.empty()) continue;

      // Divisors: coefficients of integer columns away from their bound.
      std::vector<double> divisors;
      for (const IntTerm& t : ints) {
        if (t.t <= 1e-6) continue;
        const double d = std::fabs(t.c);
        if (d < 1e-6) continue;
        if (std::none_of(divisors.begin(), divisors.end(),
                         [&](double e) { return std::fabs(e - d) <= 1e-9 * d; })) {
          divisors.push_back(d);
        }
        if (static_cast<int>(divisors.size()) >= kMaxMirDivisors) break;
      }
      const auto try_divisor = [&](double delta, Cut& out) {
        const double scaled = beta / delta;
        const double f = scaled - std::floor(scaled);
        if (f < 0.05 || f > 0.95) return false;
        // sum F(c_j/delta) t_j + sum_{h<0} h_j s_j / (delta (1-f)) <= floor(beta/delta).
        double rhs = std::floor(scaled);
        for (const IntTerm& t : ints) {
          const double a = t.c / delta;
          const double fa = a - std::floor(a);
          const double coef = std::floor(a) + std::max(0.0, fa - f) / (1.0 - f);
          if (coef == 0.0) continue;
          const double x_coef = t.at_upper ? -coef : coef;
          rhs += x_coef * t.bound;  // coef * t with t = x - l or u - x.
          builder.add(t.col, x_coef);
        }
        for (const ContTerm& s : conts) {
          if (s.h >= 0.0) continue;
          const double coef = s.h / (delta * (1.0 - f));
          const double x_coef = s.at_upper ? -coef : coef;
          rhs += x_coef * s.bound;
          builder.add(s.col, x_coef);
        }
        out = builder.take(rhs, -1.0);
        return finalize_cut(ctx, out);
      };
      Cut best;
      double best_delta = 0.0;
      for (const double delta : divisors) {
        Cut c;
        if (try_divisor(delta, c) && c.efficacy > best.efficacy) {
          best = std::move(c);
          best_delta = delta;
        }
      }
      if (best_delta > 0.0) {
        for (const double factor : {0.5, 0.25, 0.125}) {
          Cut c;
          if (try_divisor(best_delta * factor, c) && c.efficacy > best.efficacy) {
            best = std::move(c);
          }
        }
        cuts.push_back(std::move(best));
      }
    }
  }
}

namespace {

// Aggregated c-MIR limits: starting rows per call, rows added to one aggregation, and the
// smallest distance of a continuous column from its substituted bound worth eliminating.
constexpr int kMaxAggregationStarts = 300;
constexpr int kMaxAggregations = 5;
constexpr double kMinBoundDistance = 1e-6;

// A variable bound x_j <= coef * y (upper) or x_j >= coef * y (lower) with y binary, from a
// two-entry row with right-hand side 0.
struct VariableBound {
  Index y = -1;
  double coef = 0.0;
};

void find_variable_bounds(const CutContext& ctx, std::vector<VariableBound>& upper,
                          std::vector<VariableBound>& lower) {
  const Model& model = ctx.model;
  const auto start = ctx.At.col_start();
  const auto index = ctx.At.row_index();
  const auto value = ctx.At.values();
  upper.assign(static_cast<std::size_t>(model.num_cols()), {});
  lower.assign(static_cast<std::size_t>(model.num_cols()), {});
  for (Index i = 0; i < ctx.original_rows; ++i) {
    if (start[i + 1] - start[i] != 2) continue;
    for (int side = 0; side < 2; ++side) {
      if ((side == 0 ? model.row_upper[i] : model.row_lower[i]) != 0.0) continue;
      const double sign = side == 0 ? 1.0 : -1.0;  // sign * (a x + b y) <= 0.
      for (int k = 0; k < 2; ++k) {
        const Index j = index[start[i] + k];
        const Index y = index[start[i] + 1 - k];
        const double a = sign * value[start[i] + k];
        const double b = sign * value[start[i] + 1 - k];
        if (is_integer_col(model, j) || !is_integer_col(model, y) || ctx.lower[y] != 0.0 ||
            ctx.upper[y] != 1.0 || a == 0.0) {
          continue;
        }
        const double coef = -b / a;
        if (a > 0.0) {  // x <= coef y
          if (upper[j].y < 0 || coef < upper[j].coef) upper[j] = {y, coef};
        } else {  // x >= coef y
          if (lower[j].y < 0 || coef > lower[j].coef) lower[j] = {y, coef};
        }
      }
    }
  }
}

// sum_j coef_j x_j <= beta, built from original rows.
struct Aggregate {
  std::vector<double> coef;
  std::vector<char> mark;
  std::vector<Index> cols;
  double beta = 0.0;

  explicit Aggregate(Index n)
      : coef(static_cast<std::size_t>(n), 0.0), mark(static_cast<std::size_t>(n), 0) {}
  void clear() {
    for (const Index j : cols) {
      coef[j] = 0.0;
      mark[j] = 0;
    }
    cols.clear();
    beta = 0.0;
  }
  // Adds lambda times row i, using the side that keeps it a valid <= inequality.
  bool add_row(const CutContext& ctx, Index i, double lambda) {
    const double side = lambda > 0.0 ? ctx.model.row_upper[i] : ctx.model.row_lower[i];
    if (!std::isfinite(side)) return false;
    const auto start = ctx.At.col_start();
    const auto index = ctx.At.row_index();
    const auto value = ctx.At.values();
    for (NnzIndex p = start[i]; p < start[i + 1]; ++p) {
      const Index j = index[p];
      if (!mark[j]) {
        mark[j] = 1;
        cols.push_back(j);
      }
      coef[j] += lambda * value[p];
    }
    beta += lambda * side;
    return true;
  }
};

// How a continuous column is replaced by a slack s >= 0 in the MIR: s = x - l, s = u - x,
// s = c y - x (variable upper bound) or s = x - c y (variable lower bound).
enum class Substitution : std::uint8_t { kLower, kUpper, kVub, kVlb };

struct ContinuousChoice {
  Substitution kind = Substitution::kLower;
  double distance = kInf;  // The slack at the LP point.
};

ContinuousChoice choose_substitution(const CutContext& ctx, Index j,
                                     const std::vector<VariableBound>& vub,
                                     const std::vector<VariableBound>& vlb) {
  const double xj = ctx.x[j];
  ContinuousChoice best;
  const auto consider = [&](Substitution kind, double distance) {
    // Variable bounds win ties: they carry the binary into the cut (flow-cover strength).
    const bool variable = kind == Substitution::kVub || kind == Substitution::kVlb;
    if (distance < best.distance - 1e-12 || (variable && distance <= best.distance + 1e-12)) {
      best = {kind, distance};
    }
  };
  if (std::isfinite(ctx.lower[j])) consider(Substitution::kLower, xj - ctx.lower[j]);
  if (std::isfinite(ctx.upper[j])) consider(Substitution::kUpper, ctx.upper[j] - xj);
  if (vub[j].y >= 0) consider(Substitution::kVub, vub[j].coef * ctx.x[vub[j].y] - xj);
  if (vlb[j].y >= 0) consider(Substitution::kVlb, xj - vlb[j].coef * ctx.x[vlb[j].y]);
  return best;
}

// c-MIR on an aggregated row with bound substitution; the most efficacious cut over the
// divisors, if any is violated.
bool mir_on_aggregate(const CutContext& ctx, const Aggregate& agg,
                      const std::vector<VariableBound>& vub, const std::vector<VariableBound>& vlb,
                      CutBuilder& builder, Aggregate& ints_work, Cut& best) {
  const Model& model = ctx.model;
  struct ContTerm {
    Index col;
    double h;
    double bound;
    Substitution kind;
  };
  struct IntTerm {
    Index col;
    double c;
    double t;
    double bound;
    bool at_upper;
  };
  std::vector<ContTerm> conts;
  std::vector<IntTerm> ints;
  ints_work.clear();
  double beta = agg.beta;
  for (const Index j : agg.cols) {
    const double c = agg.coef[j];
    if (c == 0.0) continue;
    if (is_integer_col(model, j)) {
      ints_work.coef[j] += c;
      if (!ints_work.mark[j]) {
        ints_work.mark[j] = 1;
        ints_work.cols.push_back(j);
      }
      continue;
    }
    if (ctx.lower[j] == ctx.upper[j]) {
      beta -= c * ctx.lower[j];
      continue;
    }
    const ContinuousChoice choice = choose_substitution(ctx, j, vub, vlb);
    if (!std::isfinite(choice.distance)) return false;
    const auto add_int = [&](Index y, double v) {
      ints_work.coef[y] += v;
      if (!ints_work.mark[y]) {
        ints_work.mark[y] = 1;
        ints_work.cols.push_back(y);
      }
    };
    switch (choice.kind) {
      case Substitution::kLower:  // c x = c s + c l
        beta -= c * ctx.lower[j];
        conts.push_back({j, c, ctx.lower[j], choice.kind});
        break;
      case Substitution::kUpper:  // c x = -c s + c u
        beta -= c * ctx.upper[j];
        conts.push_back({j, -c, ctx.upper[j], choice.kind});
        break;
      case Substitution::kVub:  // c x = c v y - c s
        add_int(vub[j].y, c * vub[j].coef);
        conts.push_back({j, -c, 0.0, choice.kind});
        break;
      case Substitution::kVlb:  // c x = c v y + c s
        add_int(vlb[j].y, c * vlb[j].coef);
        conts.push_back({j, c, 0.0, choice.kind});
        break;
    }
  }
  for (const Index j : ints_work.cols) {
    const double c = ints_work.coef[j];
    const double lo = ctx.lower[j];
    const double up = ctx.upper[j];
    const double xj = ctx.x[j];
    const bool use_lower = std::isfinite(lo) && (!std::isfinite(up) || xj - lo <= up - xj);
    if (!use_lower && !std::isfinite(up)) return false;
    const double b = use_lower ? lo : up;
    beta -= c * b;
    if (lo == up || c == 0.0) continue;
    ints.push_back({j, use_lower ? c : -c, use_lower ? xj - lo : up - xj, b, !use_lower});
  }
  if (ints.empty()) return false;

  std::vector<double> divisors;
  for (const IntTerm& t : ints) {
    if (t.t <= 1e-6) continue;
    const double d = std::fabs(t.c);
    if (d < 1e-6) continue;
    if (std::none_of(divisors.begin(), divisors.end(),
                     [&](double e) { return std::fabs(e - d) <= 1e-9 * d; })) {
      divisors.push_back(d);
    }
    if (static_cast<int>(divisors.size()) >= kMaxMirDivisors) break;
  }
  const auto try_divisor = [&](double delta, Cut& out) {
    const double scaled = beta / delta;
    const double f = scaled - std::floor(scaled);
    if (f < 0.05 || f > 0.95) return false;
    double rhs = std::floor(scaled);
    for (const IntTerm& t : ints) {
      const double a = t.c / delta;
      const double fa = a - std::floor(a);
      const double coef = std::floor(a) + std::max(0.0, fa - f) / (1.0 - f);
      if (coef == 0.0) continue;
      const double x_coef = t.at_upper ? -coef : coef;
      rhs += x_coef * t.bound;
      builder.add(t.col, x_coef);
    }
    for (const ContTerm& s : conts) {
      if (s.h >= 0.0) continue;
      const double coef = s.h / (delta * (1.0 - f));
      switch (s.kind) {
        case Substitution::kLower:  // coef (x - l)
          rhs += coef * s.bound;
          builder.add(s.col, coef);
          break;
        case Substitution::kUpper:  // coef (u - x)
          rhs -= coef * s.bound;
          builder.add(s.col, -coef);
          break;
        case Substitution::kVub:  // coef (v y - x)
          builder.add(vub[s.col].y, coef * vub[s.col].coef);
          builder.add(s.col, -coef);
          break;
        case Substitution::kVlb:  // coef (x - v y)
          builder.add(s.col, coef);
          builder.add(vlb[s.col].y, -coef * vlb[s.col].coef);
          break;
      }
    }
    out = builder.take(rhs, -1.0);
    return finalize_cut(ctx, out);
  };
  bool found = false;
  double best_delta = 0.0;
  for (const double delta : divisors) {
    Cut c;
    if (try_divisor(delta, c) && (!found || c.efficacy > best.efficacy)) {
      best = std::move(c);
      best_delta = delta;
      found = true;
    }
  }
  if (found) {
    for (const double factor : {0.5, 0.25, 0.125}) {
      Cut c;
      if (try_divisor(best_delta * factor, c) && c.efficacy > best.efficacy) best = std::move(c);
    }
  }
  return found;
}

// Simple generalized flow cover inequality (Van Roy and Wolsey) on the single-node flow set of an
// aggregated row sum_j a_j x_j <= beta. Each term becomes a flow z_j = |a_j| x_j with capacity
// z_j <= u_j y_j: a binary in the row is its own flow (u = |a|, y = itself), a column with a
// variable upper bound x <= c y uses it, any other column with a finite upper bound has y fixed
// at 1. With N+ the terms of positive coefficient, N- the others, a cover C+ of N+ with
// lambda = sum_{C+} u_j - beta > 0 gives
//   sum_{C+} z_j + sum_{C+} (u_j - lambda)^+ (1 - y_j) - sum_{L-} lambda y_j - sum_{N- \ L-} z_j
//     <= beta,
// with L- the N- terms where lambda y_j < z_j at the LP point (constant y_j = 1 counts as the
// constant lambda). Columns with nonzero lower bounds are shifted; rows with general integer
// columns unbounded above are skipped.
bool flow_cover_on_aggregate(const CutContext& ctx, const Aggregate& agg,
                             const std::vector<VariableBound>& vub, CutBuilder& builder,
                             Cut& out) {
  const Model& model = ctx.model;
  struct Flow {
    Index x;        // The column of the flow.
    double a;       // |a_j|: z = a x (after the shift by the lower bound).
    double u;       // Capacity of z.
    Index y;        // Its binary, or -1 when y is the constant 1.
    double z;       // z at the LP point.
    double ylp;     // y at the LP point.
    bool positive;  // In N+.
  };
  std::vector<Flow> flows;
  double beta = agg.beta;
  for (const Index j : agg.cols) {
    const double c = agg.coef[j];
    if (c == 0.0) continue;
    const double lo = ctx.lower[j];
    const double up = ctx.upper[j];
    if (!std::isfinite(lo)) return false;
    const bool binary = is_integer_col(model, j) && lo == 0.0 && up == 1.0;
    Flow f{j, std::fabs(c), kInf, -1, 0.0, 1.0, c > 0.0};
    if (binary) {
      f.u = f.a;
      f.y = j;
      f.ylp = ctx.x[j];
      f.z = f.a * ctx.x[j];
    } else {
      beta -= c * lo;  // Shift to x' = x - lo >= 0.
      f.z = f.a * (ctx.x[j] - lo);
      if (lo == 0.0 && vub[j].y >= 0 && !is_integer_col(model, j)) {
        f.u = f.a * vub[j].coef;
        f.y = vub[j].y;
        f.ylp = ctx.x[f.y];
      } else if (std::isfinite(up)) {
        f.u = f.a * (up - lo);
      }
    }
    flows.push_back(f);
  }
  // Cover: N+ flows with a finite capacity, those with y near 1 first, until they exceed beta.
  std::vector<std::size_t> order;
  for (std::size_t k = 0; k < flows.size(); ++k) {
    if (flows[k].positive && std::isfinite(flows[k].u)) order.push_back(k);
  }
  std::sort(order.begin(), order.end(), [&](std::size_t p, std::size_t q) {
    return 1.0 - flows[p].ylp < 1.0 - flows[q].ylp ||
           (1.0 - flows[p].ylp == 1.0 - flows[q].ylp && flows[p].u > flows[q].u);
  });
  std::vector<char> in_cover(flows.size(), 0);
  double capacity = 0.0;
  const double eps = 1e-9 * (1.0 + std::fabs(beta));
  for (const std::size_t k : order) {
    if (capacity > beta + eps) break;
    in_cover[k] = 1;
    capacity += flows[k].u;
  }
  const double lambda = capacity - beta;
  if (!(lambda > eps)) return false;

  // sum_{C+} a x + sum_{C++} (u - lambda)(1 - y) - sum_{L-} lambda y - sum_{rest of N-} a x
  //   <= beta   (x shifted by its lower bound; constants move to the right-hand side).
  double rhs = beta;
  for (std::size_t k = 0; k < flows.size(); ++k) {
    const Flow& f = flows[k];
    const double lo = f.y == f.x ? 0.0 : ctx.lower[f.x];
    if (f.positive) {
      if (!in_cover[k]) continue;
      if (f.y == f.x) {
        // A binary flow: a y + (a - lambda)^+ (1 - y).
        const double extra = std::max(0.0, f.u - lambda);
        builder.add(f.x, f.a - extra);
        rhs -= extra;
      } else {
        builder.add(f.x, f.a);
        rhs += f.a * lo;
        if (f.y >= 0 && f.u > lambda) {
          builder.add(f.y, -(f.u - lambda));
          rhs -= f.u - lambda;
        }
      }
    } else {
      // On a tie (x = y = 0 at the LP point) a variable-bound flow enters as lambda y, which
      // stays strong wherever the LP reroutes flow through it (x = M y there); a binary's
      // own flow keeps x.
      const bool use_y = f.y == -1        ? lambda < f.z
                         : f.y == f.x    ? lambda * f.ylp < f.z
                                         : lambda * f.ylp <= f.z;
      if (use_y) {
        if (f.y >= 0) {
          builder.add(f.y, -lambda);
        } else {
          rhs += lambda;  // -lambda * 1 on the left.
        }
      } else if (f.y == f.x) {
        builder.add(f.x, -f.a);
      } else {
        builder.add(f.x, -f.a);
        rhs -= f.a * lo;
      }
    }
  }
  out = builder.take(rhs, -1.0);
  return finalize_cut(ctx, out);
}

// Distance of a continuous column from its nearer simple bound: which columns are worth
// eliminating by aggregation (a column on a variable bound still counts as interior here, as in
// SCIP's aggregation heuristic; the substitution then uses the variable bound).
double simple_bound_distance(const CutContext& ctx, Index j) {
  const double xj = ctx.x[j];
  return std::min(std::isfinite(ctx.lower[j]) ? xj - ctx.lower[j] : kInf,
                  std::isfinite(ctx.upper[j]) ? ctx.upper[j] - xj : kInf);
}

}  // namespace

Index separate_aggregated_mir(const CutContext& ctx, std::vector<Cut>& cuts, Index first_row) {
  const Model& model = ctx.model;
  const Index n = model.num_cols();
  const auto rstart = ctx.At.col_start();
  const auto cstart = model.A.col_start();
  const auto cindex = model.A.row_index();
  const auto cvalue = model.A.values();
  std::vector<VariableBound> vub;
  std::vector<VariableBound> vlb;
  find_variable_bounds(ctx, vub, vlb);
  std::vector<double> activity(static_cast<std::size_t>(model.num_rows()), 0.0);
  model.A.multiply(ctx.x, activity);

  // Starting rows: those with a continuous column away from its substituted bound, from
  // first_row on (wrapping around) until the cap; the next call continues where this one stopped.
  std::vector<Index> starts;
  const Index rows = ctx.original_rows;
  Index next_row = 0;
  if (first_row < 0 || first_row >= rows) first_row = 0;
  for (Index k = 0; k < rows; ++k) {
    const Index i = first_row + k < rows ? first_row + k : first_row + k - rows;
    if (static_cast<int>(starts.size()) >= kMaxAggregationStarts) {
      next_row = i;
      break;
    }
    if (rstart[i + 1] - rstart[i] > kMaxRowLength) continue;
    const auto index = ctx.At.row_index();
    for (NnzIndex p = rstart[i]; p < rstart[i + 1]; ++p) {
      const Index j = index[p];
      if (!is_integer_col(model, j) && simple_bound_distance(ctx, j) > kMinBoundDistance) {
        starts.push_back(i);
        break;
      }
    }
  }
  CutBuilder builder(ctx);
  Aggregate agg(n);
  Aggregate ints_work(n);
  std::vector<Index> used;

  for (const Index i0 : starts) {
    // Both finite sides, for an equality too: c-MIR and flow covers give different cuts for
    // a'x <= b and -a'x <= -b. On a demand node's balance row (inflow - outflow = d) the <= side
    // gives x_a <= d y_a + (M - d) (outflows' y), the >= side the cut-set inequality
    // d (inflows' y) + (unused inflows' x) >= d.
    const double ru = model.row_upper[i0];
    const double rl = model.row_lower[i0];
    for (int side = 0; side < 2; ++side) {
      const double lambda = side == 0 ? 1.0 : -1.0;
      if (!std::isfinite(side == 0 ? ru : rl)) continue;
      agg.clear();
      agg.add_row(ctx, i0, lambda);
      used.assign(1, i0);
      for (int step = 0; step <= kMaxAggregations; ++step) {
        Cut cut;
        if (mir_on_aggregate(ctx, agg, vub, vlb, builder, ints_work, cut) ||
            flow_cover_on_aggregate(ctx, agg, vub, builder, cut)) {
          cuts.push_back(std::move(cut));
          break;
        }
        if (step == kMaxAggregations) break;
        // Eliminate the continuous column farthest from its bound.
        Index pick = -1;
        double pick_distance = kMinBoundDistance;
        for (const Index j : agg.cols) {
          if (agg.coef[j] == 0.0 || is_integer_col(model, j)) continue;
          const double d = simple_bound_distance(ctx, j);
          if (d > pick_distance) {
            pick_distance = d;
            pick = j;
          }
        }
        if (pick < 0) break;
        // With the tightest row containing it that is not used yet.
        Index row = -1;
        double row_lambda = 0.0;
        double row_slack = kInf;
        for (NnzIndex p = cstart[pick]; p < cstart[pick + 1]; ++p) {
          const Index i = cindex[p];
          if (i >= ctx.original_rows || cvalue[p] == 0.0 ||
              std::find(used.begin(), used.end(), i) != used.end() ||
              rstart[i + 1] - rstart[i] > kMaxRowLength) {
            continue;
          }
          const double l = -agg.coef[pick] / cvalue[p];
          const double bound = l > 0.0 ? model.row_upper[i] : model.row_lower[i];
          if (!std::isfinite(bound)) continue;
          const double slack = std::fabs(bound - activity[i]);
          if (slack < row_slack) {
            row_slack = slack;
            row = i;
            row_lambda = l;
          }
        }
        if (row < 0) break;
        agg.add_row(ctx, row, row_lambda);
        agg.coef[pick] = 0.0;  // Eliminated exactly.
        used.push_back(row);
      }
    }
  }
  return next_row;
}

void separate_knapsack_covers(const CutContext& ctx, std::vector<Cut>& cuts) {
  const Model& model = ctx.model;
  const auto start = ctx.At.col_start();
  const auto index = ctx.At.row_index();
  const auto value = ctx.At.values();
  CutBuilder builder(ctx);
  struct Item {
    Index col;
    double w;
    double y;  // LP value of the (possibly complemented) binary.
    bool complemented;
  };
  std::vector<Item> items;
  for (Index i = 0; i < ctx.original_rows; ++i) {
    if (start[i + 1] - start[i] > kMaxRowLength || start[i + 1] - start[i] < 2) continue;
    for (int side = 0; side < 2; ++side) {
      const double bound = side == 0 ? model.row_upper[i] : model.row_lower[i];
      if (!std::isfinite(bound)) continue;
      const double sign = side == 0 ? 1.0 : -1.0;
      double capacity = sign * bound;
      items.clear();
      bool ok = true;
      for (NnzIndex p = start[i]; p < start[i + 1]; ++p) {
        const Index j = index[p];
        const double w = sign * value[p];
        const bool binary = is_integer_col(model, j) && ctx.lower[j] == 0.0 && ctx.upper[j] == 1.0;
        if (!binary) {
          // Relax to the smallest contribution.
          const double least = w > 0.0 ? w * ctx.lower[j] : w * ctx.upper[j];
          if (!std::isfinite(least)) {
            ok = false;
            break;
          }
          capacity -= least;
          continue;
        }
        if (w > 0.0) {
          items.push_back({j, w, ctx.x[j], false});
        } else {
          capacity -= w;  // w x = w - w (1 - x).
          items.push_back({j, -w, 1.0 - ctx.x[j], true});
        }
      }
      if (!ok || items.size() < 2 || capacity < 0.0) continue;
      const double total = std::accumulate(items.begin(), items.end(), 0.0,
                                           [](double s, const Item& it) { return s + it.w; });
      const double eps = 1e-9 * (1.0 + std::fabs(capacity));
      if (total <= capacity + eps) continue;

      // Greedy cover: prefer items with large LP value per unit weight.
      std::sort(items.begin(), items.end(), [](const Item& a, const Item& b) {
        return (1.0 - a.y) * b.w < (1.0 - b.y) * a.w;
      });
      std::vector<std::size_t> cover;
      double weight = 0.0;
      for (std::size_t t = 0; t < items.size() && weight <= capacity + eps; ++t) {
        cover.push_back(t);
        weight += items[t].w;
      }
      if (weight <= capacity + eps) continue;
      // Make it minimal, removing items with the smallest LP value first.
      std::vector<std::size_t> order = cover;
      std::sort(order.begin(), order.end(),
                [&](std::size_t a, std::size_t b) { return items[a].y < items[b].y; });
      for (const std::size_t t : order) {
        if (weight - items[t].w > capacity + eps) {
          weight -= items[t].w;
          cover.erase(std::find(cover.begin(), cover.end(), t));
        }
      }
      double lhs = 0.0;
      double max_w = 0.0;
      for (const std::size_t t : cover) {
        lhs += items[t].y;
        max_w = std::max(max_w, items[t].w);
      }
      const double size = static_cast<double>(cover.size());
      if (lhs <= size - 1.0 + 1e-6) continue;
      // Extended cover: every item at least as heavy as the heaviest cover item joins.
      std::vector<char> in_cover(items.size(), 0);
      for (const std::size_t t : cover) in_cover[t] = 1;
      double rhs = size - 1.0;  // sum_E y <= |C| - 1, with y = x or 1 - x.
      for (std::size_t t = 0; t < items.size(); ++t) {
        if (!in_cover[t] && items[t].w < max_w) continue;
        if (items[t].complemented) {
          builder.add(items[t].col, -1.0);
          rhs -= 1.0;
        } else {
          builder.add(items[t].col, 1.0);
        }
      }
      Cut cut = builder.take(rhs, -1.0);
      if (finalize_cut(ctx, cut)) cuts.push_back(std::move(cut));
    }
  }
}

namespace {

// A literal is 2 j (x_j = 1) or 2 j + 1 (x_j = 0).
constexpr int kMaxCliqueSeeds = 500;        // Fractional literals tried per call.
constexpr long long kMaxCliqueWork = 2000000;  // Neighbour visits per call (dense graphs).

}  // namespace

void separate_cliques(const CutContext& ctx, std::vector<Cut>& cuts) {
  const Model& model = ctx.model;
  const Index n = model.num_cols();
  const auto start = ctx.At.col_start();
  const auto index = ctx.At.row_index();
  const auto value = ctx.At.values();
  const auto is_binary = [&](Index j) {
    return is_integer_col(model, j) && ctx.lower[j] == 0.0 && ctx.upper[j] == 1.0;
  };

  // Clique table: from each row side, its heaviest literals that pairwise exceed the room.
  std::vector<std::vector<Index>> cliques;
  struct Item {
    Index lit;
    double w;
  };
  std::vector<Item> items;
  for (Index i = 0; i < ctx.original_rows; ++i) {
    if (start[i + 1] - start[i] > kMaxRowLength || start[i + 1] - start[i] < 2) continue;
    for (int side = 0; side < 2; ++side) {
      const double bound = side == 0 ? model.row_upper[i] : model.row_lower[i];
      if (!std::isfinite(bound)) continue;
      const double sign = side == 0 ? 1.0 : -1.0;
      double room = sign * bound;  // For  sum sign a_j x_j <= sign bound.
      items.clear();
      bool ok = true;
      for (NnzIndex p = start[i]; p < start[i + 1]; ++p) {
        const Index j = index[p];
        const double w = sign * value[p];
        if (w == 0.0) continue;
        if (!is_binary(j)) {
          const double least = w > 0.0 ? w * ctx.lower[j] : w * ctx.upper[j];
          if (!std::isfinite(least)) {
            ok = false;
            break;
          }
          room -= least;
        } else if (w > 0.0) {
          items.push_back({2 * j, w});  // x_j = 1 uses w of the room.
        } else {
          room -= w;                         // w x_j = w - w (1 - x_j):
          items.push_back({2 * j + 1, -w});  // x_j = 0 uses -w of it.
        }
      }
      if (!ok || items.size() < 2 || room < 0.0) continue;
      std::sort(items.begin(), items.end(), [](const Item& a, const Item& b) { return a.w > b.w; });
      // Round-off must never create a conflict: two literals conflict only if they exceed the
      // room by a margin.
      const double eps = 1e-9 * (1.0 + std::fabs(room));
      std::size_t k = 1;
      while (k < items.size() && items[k - 1].w + items[k].w > room + eps) ++k;
      if (k < 2) continue;
      std::vector<Index> clique;
      clique.reserve(k);
      for (std::size_t t = 0; t < k; ++t) clique.push_back(items[t].lit);
      cliques.push_back(std::move(clique));
    }
  }
  if (cliques.empty()) return;
  std::vector<std::vector<int>> lit_cliques(static_cast<std::size_t>(2 * n));
  for (std::size_t c = 0; c < cliques.size(); ++c) {
    for (const Index l : cliques[c]) lit_cliques[l].push_back(static_cast<int>(c));
  }
  const auto lit_value = [&](Index l) { return l % 2 == 0 ? ctx.x[l / 2] : 1.0 - ctx.x[l / 2]; };

  // Seeds: fractional literals with a conflict, largest LP value first.
  std::vector<Index> seeds;
  for (Index l = 0; l < 2 * n; ++l) {
    const double v = lit_value(l);
    if (!lit_cliques[l].empty() && v > kMinViolation && v < 1.0 - kMinViolation) seeds.push_back(l);
  }
  std::sort(seeds.begin(), seeds.end(),
            [&](Index a, Index b) { return lit_value(a) > lit_value(b); });
  if (static_cast<int>(seeds.size()) > kMaxCliqueSeeds) seeds.resize(kMaxCliqueSeeds);

  // count[l]: members of the current clique that l conflicts with; stamp[l] marks the literals
  // already counted for the member being added.
  std::vector<int> count(static_cast<std::size_t>(2 * n), 0);
  std::vector<int> stamp(static_cast<std::size_t>(2 * n), -1);
  std::vector<Index> touched;
  std::set<std::vector<Index>> found;
  long long work = 0;
  CutBuilder builder(ctx);
  int stamp_id = 0;
  for (const Index seed : seeds) {
    if (work > kMaxCliqueWork) break;
    std::vector<Index> members;
    // Adds l to the clique and counts it for every literal it conflicts with.
    const auto add_member = [&](Index l) {
      members.push_back(l);
      ++stamp_id;
      for (const int c : lit_cliques[l]) {
        for (const Index o : cliques[static_cast<std::size_t>(c)]) {
          ++work;
          if (o == l || stamp[o] == stamp_id) continue;
          stamp[o] = stamp_id;
          if (count[o] == 0) touched.push_back(o);
          ++count[o];
        }
      }
    };
    add_member(seed);
    double weight = lit_value(seed);
    for (;;) {
      // The candidate conflicting with every member that has the largest LP value (a zero one
      // still makes the clique, and so the cut, stronger elsewhere).
      Index best = -1;
      double best_value = -1.0;
      const auto size = static_cast<int>(members.size());
      work += static_cast<long long>(touched.size());
      for (const Index o : touched) {
        if (count[o] != size) continue;
        if (std::find(members.begin(), members.end(), o ^ 1) != members.end()) continue;
        const double v = lit_value(o);
        if (v > best_value) {
          best_value = v;
          best = o;
        }
      }
      if (best < 0) break;
      add_member(best);
      weight += std::max(0.0, best_value);
    }
    for (const Index o : touched) count[o] = 0;
    touched.clear();
    if (members.size() < 2 || weight <= 1.0 + kMinViolation) continue;
    std::vector<Index> key = members;
    std::sort(key.begin(), key.end());
    if (!found.insert(key).second) continue;
    // sum over x_j literals of x_j + sum over complemented ones of (1 - x_j) <= 1.
    double rhs = 1.0;
    for (const Index l : members) {
      if (l % 2 == 0) {
        builder.add(l / 2, 1.0);
      } else {
        builder.add(l / 2, -1.0);
        rhs -= 1.0;
      }
    }
    Cut cut = builder.take(rhs, -1.0);
    if (finalize_cut(ctx, cut)) cuts.push_back(std::move(cut));
  }
}

}  // namespace samaya
