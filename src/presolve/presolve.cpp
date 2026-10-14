#include "presolve/presolve.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <unordered_map>

namespace samaya {

namespace {

// A free column singleton is substituted out only if its coefficient is not tiny relative to the
// rest of its row (as in threshold pivoting), to keep the substitution stable.
constexpr double kSubstitutionPivotRatio = 0.01;
// Prove redundant bounds with an outward allowance, rather than accepting a near miss.
constexpr double kImpliedFreeErrorTol = 1e-12;
// Reject ill-conditioned retained residuals; never approximate them by zero. Applied both
// against the row's largest coefficient and against the operands of the update itself
// (cancelled()): when old - term cancels to a few ulps of the operands, the computed result is
// rounding noise, and a row left with only noise coefficients would turn into a wrong bound.
constexpr double kSubstitutionResidualRatio = 1e-9;
constexpr Index kSubstitutionRowLength = 10;  // Bound the first version's equation work.
constexpr Index kSubstitutionColLength = 10;  // Bound the number of other rows updated.
constexpr Index kSubstitutionFillLimit = 10;  // Follows HiGHS's default; not tuned.
constexpr double kIntegerDataLimit = 2147483648.0;  // Exact integer arithmetic fits int64_t.
constexpr double kBoundRoundoff = 4.0 * 0x1p-52;

// True when the updated coefficient old - term is nonzero but tiny against its operands.
bool cancelled(double old, double term, double result) {
  return result != 0.0 &&
         std::fabs(result) < kSubstitutionResidualRatio * (std::fabs(old) + std::fabs(term));
}

}  // namespace

Presolve::Presolve(const Model& model, const PresolveOptions& options)
    : model_(model),
      options_(options),
      m_(model.num_rows()),
      n_(model.num_cols()),
      sense_(model.sense == ObjSense::kMaximize ? -1.0 : 1.0) {
  cost_.resize(static_cast<std::size_t>(n_));
  for (Index j = 0; j < n_; ++j) cost_[j] = sense_ * model.obj[j];
  offset_ = sense_ * model.obj_offset;
  lower_ = model.col_lower;
  upper_ = model.col_upper;
  row_lower_ = model.row_lower;
  row_upper_ = model.row_upper;
  row_active_.assign(static_cast<std::size_t>(m_), 1);
  col_active_.assign(static_cast<std::size_t>(n_), 1);
  row_count_.assign(static_cast<std::size_t>(m_), 0);
  col_count_.assign(static_cast<std::size_t>(n_), 0);
  const auto start = model.A.col_start();
  const auto index = model.A.row_index();
  const auto val = model.A.values();
  rows_.resize(static_cast<std::size_t>(m_));
  cols_.resize(static_cast<std::size_t>(n_));
  for (Index j = 0; j < n_; ++j) {
    col_count_[j] = static_cast<Index>(start[j + 1] - start[j]);
    cols_[j].resize(static_cast<std::size_t>(col_count_[j]));
    for (NnzIndex p = start[j]; p < start[j + 1]; ++p) ++row_count_[index[p]];
  }
  for (Index i = 0; i < m_; ++i) rows_[i].resize(static_cast<std::size_t>(row_count_[i]));
  std::vector<Index> next(static_cast<std::size_t>(m_), 0);
  // Scattering CSC columns in increasing order exactly reproduces transpose()'s row order.
  for (Index j = 0; j < n_; ++j) {
    auto& col = cols_[j];
    for (NnzIndex p = start[j]; p < start[j + 1]; ++p) {
      const Index i = index[p];
      rows_[i][next[i]++] = {j, val[p]};
      col[static_cast<std::size_t>(p - start[j])] = {i, val[p]};
    }
  }
}

double Presolve::tol(double v) const {
  return options_.feasibility_tol * (1.0 + (std::isfinite(v) ? std::fabs(v) : 0.0));
}

bool Presolve::add_entry(Index i, Index j, double a) {
  if (!row_active_[i] || !col_active_[j] || a == 0.0) return false;
  const auto by_index = [](const Entry& e, Index k) { return e.index < k; };
  auto& row = rows_[i];
  auto& col = cols_[j];
  const auto r = std::lower_bound(row.begin(), row.end(), j, by_index);
  if (r != row.end() && r->index == j) return false;
  const auto c = std::lower_bound(col.begin(), col.end(), i, by_index);
  row.insert(r, {j, a});
  col.insert(c, {i, a});
  ++row_count_[i];
  ++col_count_[j];
  return true;
}

bool Presolve::change_coefficient(Index i, Index j, double a) {
  if (!row_active_[i] || !col_active_[j]) return false;
  if (a == 0.0) return remove_entry(i, j);
  const auto by_index = [](const Entry& e, Index k) { return e.index < k; };
  auto& row = rows_[i];
  auto& col = cols_[j];
  const auto r = std::lower_bound(row.begin(), row.end(), j, by_index);
  if (r == row.end() || r->index != j) return false;
  const auto c = std::lower_bound(col.begin(), col.end(), i, by_index);
  r->a = a;
  c->a = a;
  return true;
}

bool Presolve::remove_entry(Index i, Index j) {
  const auto by_index = [](const Entry& e, Index k) { return e.index < k; };
  auto& row = rows_[i];
  auto& col = cols_[j];
  const auto r = std::lower_bound(row.begin(), row.end(), j, by_index);
  if (r == row.end() || r->index != j) return false;
  const auto c = std::lower_bound(col.begin(), col.end(), i, by_index);
  row.erase(r);
  col.erase(c);
  if (row_active_[i] && col_active_[j]) {
    --row_count_[i];
    --col_count_[j];
  }
  return true;
}

void Presolve::snapshot_column(Index j, std::size_t& first, std::size_t& last) {
  first = column_terms_.size();
  if (col_active_[j]) {
    for (const auto& e : cols_[j]) {
      if (row_active_[e.index]) column_terms_.push_back(e);
    }
  }
  last = column_terms_.size();
}

void Presolve::remove_row(Index i) {
  if (!row_active_[i]) return;
  row_active_[i] = 0;
  row_count_[i] = 0;
  ++stats_.rows_removed;
  for (const auto& e : rows_[i]) {
    if (col_active_[e.index]) --col_count_[e.index];
  }
}

void Presolve::remove_col(Index j) {
  if (!col_active_[j]) return;
  col_active_[j] = 0;
  col_count_[j] = 0;
  ++stats_.cols_removed;
  for (const auto& e : cols_[j]) {
    if (row_active_[e.index]) --row_count_[e.index];
  }
}

void Presolve::fix_col(Index j, double value) {
  for (const auto& e : cols_[j]) {
    const Index i = e.index;
    if (!row_active_[i]) continue;
    if (row_lower_[i] > -kInf) row_lower_[i] -= e.a * value;
    if (row_upper_[i] < kInf) row_upper_[i] -= e.a * value;
  }
  offset_ += cost_[j] * value;
  remove_col(j);
}

bool Presolve::set_bounds(Index j, double lower, double upper) {
  if (is_int(j)) {
    if (lower > -kInf) lower = std::ceil(lower - options_.integrality_tol);
    if (upper < kInf) upper = std::floor(upper + options_.integrality_tol);
  }
  lower = std::max(lower, lower_[j]);
  upper = std::min(upper, upper_[j]);
  if (lower > upper) {
    if (lower > upper + tol(upper)) return false;
    // Crossed within the tolerance: keep whichever bound was already there.
    if (lower == lower_[j]) {
      upper = lower;
    } else {
      lower = upper;
    }
  }
  if (lower != lower_[j] || upper != upper_[j]) ++stats_.bounds_tightened;
  lower_[j] = lower;
  upper_[j] = upper;
  return true;
}

void Presolve::activity_bounds(Index i, double& min_act, double& max_act) const {
  min_act = 0.0;
  max_act = 0.0;
  for (const auto& e : rows_[i]) {
    const Index j = e.index;
    if (!col_active_[j]) continue;
    const double a = e.a;
    if (a > 0.0) {
      min_act += lower_[j] > -kInf ? a * lower_[j] : -kInf;
      max_act += upper_[j] < kInf ? a * upper_[j] : kInf;
    } else {
      min_act += upper_[j] < kInf ? a * upper_[j] : -kInf;
      max_act += lower_[j] > -kInf ? a * lower_[j] : kInf;
    }
  }
}

void Presolve::other_activity_bounds(Index i, Index j, ActivitySide& lower,
                                     ActivitySide& upper) const {
  lower = {};
  upper = {};
  const auto add = [](ActivitySide& side, double term) {
    if (std::isinf(term)) {
      ++side.infinite;
    } else {
      side.sum += term;
      side.magnitude += std::fabs(term);
    }
  };
  // Never subtract the singleton's contribution from an activity that may be infinite.
  for (const auto& e : rows_[i]) {
    const Index k = e.index;
    if (k == j || !col_active_[k]) continue;
    add(lower, e.a * (e.a > 0.0 ? lower_[k] : upper_[k]));
    add(upper, e.a * (e.a > 0.0 ? upper_[k] : lower_[k]));
  }
}

bool Presolve::implied_free(Index i, Index j, double a, bool& tolerance_only) const {
  tolerance_only = false;
  ActivitySide lower, upper;
  other_activity_bounds(i, j, lower, upper);
  const double rl = row_lower_[i];
  const double ru = row_upper_[i];
  if (!std::isfinite(a) || a == 0.0 || std::isnan(rl) || std::isnan(ru) ||
      !std::isfinite(lower.sum) || !std::isfinite(upper.sum) ||
      !std::isfinite(lower.magnitude) || !std::isfinite(upper.magnitude)) return false;

  double lo = (std::isfinite(rl) && upper.infinite == 0) ? (rl - upper.sum) / a
                                                       : -kInf / a;
  double up = (std::isfinite(ru) && lower.infinite == 0) ? (ru - lower.sum) / a
                                                       : kInf / a;
  if (a < 0.0) std::swap(lo, up);
  if (std::isnan(lo) || std::isnan(up) || lo > up) return false;
  const bool prove_lower = lower_[j] != -kInf;
  const bool prove_upper = upper_[j] != kInf;
  if ((prove_lower && (!std::isfinite(lo) || !(lo >= lower_[j]))) ||
      (prove_upper && (!std::isfinite(up) || !(up <= upper_[j])))) return false;

  // Include the row bound in the subtraction's magnitude and the column bound in the
  // comparison's scale. Round the shrunken bounds inward so addition cannot erase the margin.
  const double abs_a = std::fabs(a);
  const double lo_scale = a > 0.0 ? std::fabs(rl) + upper.magnitude
                                 : std::fabs(ru) + lower.magnitude;
  const double up_scale = a > 0.0 ? std::fabs(ru) + lower.magnitude
                                 : std::fabs(rl) + upper.magnitude;
  const double lo_error = kImpliedFreeErrorTol * (lo_scale / abs_a + std::fabs(lower_[j]));
  const double up_error = kImpliedFreeErrorTol * (up_scale / abs_a + std::fabs(upper_[j]));
  const bool lower_ok = !prove_lower ||
      (std::isfinite(lo_error) && lo >= std::nextafter(lower_[j] + lo_error, kInf));
  const bool upper_ok = !prove_upper ||
      (std::isfinite(up_error) && up <= std::nextafter(upper_[j] - up_error, -kInf));
  tolerance_only = !lower_ok || !upper_ok;
  return !tolerance_only;
}

bool Presolve::process_rows(bool& changed) {
  for (Index i = 0; i < m_; ++i) {
    if (!row_active_[i]) continue;
    const double rl = row_lower_[i];
    const double ru = row_upper_[i];
    if (row_count_[i] == 0) {
      if (rl > tol(rl) || ru < -tol(ru)) return false;
      records_.push_back({.kind = Kind::kRemoveRow, .row = i});
      remove_row(i);
      ++stats_.empty_rows;
      changed = true;
      continue;
    }
    if (row_count_[i] == 1) {
      const auto entry = std::find_if(rows_[i].begin(), rows_[i].end(),
                                      [&](const Entry& e) { return col_active_[e.index]; });
      const Index j = entry->index;
      const double a = entry->a;
      double lo = a > 0.0 ? rl / a : ru / a;
      double up = a > 0.0 ? ru / a : rl / a;
      if (std::isnan(lo)) lo = -kInf;
      if (std::isnan(up)) up = kInf;
      records_.push_back({.kind = Kind::kSingletonRow,
                          .row = i,
                          .col = j,
                          .a = a,
                          .cost = cost_[j],
                          .old_lower = lower_[j],
                          .old_upper = upper_[j]});
      snapshot_column(j, records_.back().first, records_.back().last);
      remove_row(i);
      if (!set_bounds(j, lo, up)) return false;
      ++stats_.singleton_rows;
      changed = true;
      continue;
    }

    double min_act;
    double max_act;
    activity_bounds(i, min_act, max_act);
    if (min_act > ru + tol(ru) || max_act < rl - tol(rl)) return false;
    const bool lower_redundant = rl == -kInf || min_act >= rl - tol(rl);
    const bool upper_redundant = ru == kInf || max_act <= ru + tol(ru);
    if (lower_redundant && upper_redundant) {
      records_.push_back({.kind = Kind::kRemoveRow, .row = i});
      remove_row(i);
      ++stats_.redundant_rows;
      changed = true;
      continue;
    }
    const bool force_min = min_act > -kInf && ru < kInf && min_act >= ru - tol(ru);
    const bool force_max = max_act < kInf && rl > -kInf && max_act <= rl + tol(rl);
    if (force_min || force_max) {
      // Every column sits at the bound that minimizes (maximizes) the activity.
      Record record{.kind = Kind::kForcingRow, .row = i, .flag1 = force_min};
      record.first = terms_.size();
      for (const auto& e : rows_[i]) {
        const Index j = e.index;
        if (!col_active_[j]) continue;
        const double a = e.a;
        const double v = (a > 0.0) == force_min ? lower_[j] : upper_[j];
        Term term{j, a, cost_[j], v, lower_[j] < upper_[j]};
        snapshot_column(j, term.column_first, term.column_last);
        terms_.push_back(term);
      }
      record.last = terms_.size();
      records_.push_back(record);
      remove_row(i);
      for (std::size_t t = record.first; t < record.last; ++t) {
        fix_col(terms_[t].col, terms_[t].value);
      }
      ++stats_.forcing_rows;
      changed = true;
    }
  }
  return true;
}

bool Presolve::process_cols(bool& changed) {
  for (Index j = 0; j < n_; ++j) {
    if (!col_active_[j]) continue;
    if (lower_[j] == upper_[j]) {
      records_.push_back({.kind = Kind::kFixCol, .col = j, .value = lower_[j]});
      fix_col(j, lower_[j]);
      ++stats_.fixed_cols;
      changed = true;
      continue;
    }
    const double c = cost_[j];
    // Whether decreasing (increasing) x_j can never make an active row infeasible.
    bool down_ok = true;
    bool up_ok = true;
    for (const auto& e : cols_[j]) {
      if (!down_ok && !up_ok) break;
      const Index i = e.index;
      if (!row_active_[i]) continue;
      const bool has_lower = row_lower_[i] > -kInf;
      const bool has_upper = row_upper_[i] < kInf;
      if (e.a > 0.0) {
        down_ok = down_ok && !has_lower;
        up_ok = up_ok && !has_upper;
      } else {
        down_ok = down_ok && !has_upper;
        up_ok = up_ok && !has_lower;
      }
    }
    double value;
    if (c >= 0.0 && down_ok && lower_[j] > -kInf) {
      value = lower_[j];
    } else if (c <= 0.0 && up_ok && upper_[j] < kInf) {
      value = upper_[j];
    } else if (c == 0.0 && down_ok && up_ok) {
      value = std::clamp(0.0, lower_[j], upper_[j]);  // Free and cost-free: any value will do.
    } else {
      continue;  // Includes unbounded directions, which are left to the simplex.
    }
    records_.push_back({.kind = Kind::kFixCol, .col = j, .value = value});
    const bool empty = col_count_[j] == 0;
    fix_col(j, value);
    if (empty) {
      ++stats_.empty_cols;
    } else {
      ++stats_.dominated_cols;
    }
    changed = true;
  }
  return true;
}

bool Presolve::process_singleton_cols(bool& changed) {
  for (Index j = 0; j < n_; ++j) {
    if (!col_active_[j] || col_count_[j] != 1) continue;
    if (options_.mip && is_int(j)) continue;
    const auto entry = std::find_if(cols_[j].begin(), cols_[j].end(),
                                    [&](const Entry& e) { return row_active_[e.index]; });
    const Index i = entry->index;
    const double a = entry->a;
    const bool equality = row_lower_[i] == row_upper_[i];
    if (!equality && cost_[j] != 0.0) continue;

    double row_max = 0.0;
    for (const auto& e : rows_[i]) {
      if (col_active_[e.index]) row_max = std::max(row_max, std::fabs(e.a));
    }
    if (std::fabs(a) < kSubstitutionPivotRatio * row_max) continue;

    const bool free = lower_[j] == -kInf && upper_[j] == kInf;
    if (!free) {
      bool tolerance_only;
      // Bounds can change between reductions: prove this immediately before substitution.
      if (!implied_free(i, j, a, tolerance_only)) {
        if (tolerance_only) ++stats_.implied_free_rejected_by_tolerance;
        continue;
      }
    }

    Record record{.kind = equality ? Kind::kSubstituteCol : Kind::kSlackCol,
                  .row = i,
                  .col = j,
                  .a = a,
                  .cost = cost_[j],
                  .rhs_lower = row_lower_[i],
                  .rhs_upper = row_upper_[i]};
    record.first = terms_.size();
    for (const auto& e : rows_[i]) {
      const Index k = e.index;
      if (k == j || !col_active_[k]) continue;
      terms_.push_back({k, e.a, 0.0, 0.0, false});
    }
    record.last = terms_.size();
    if (equality && cost_[j] != 0.0) {
      // c_j x_j = c_j (b - sum_k a_k x_k) / a moves into the other costs and the offset.
      const double ratio = cost_[j] / a;
      for (std::size_t t = record.first; t < record.last; ++t) {
        cost_[terms_[t].col] -= ratio * terms_[t].a;
      }
      offset_ += ratio * row_lower_[i];
    }
    records_.push_back(record);
    remove_col(j);
    remove_row(i);
    if (free) {
      ++stats_.free_col_singletons;
    } else {
      ++stats_.implied_free_col_singletons;
    }
    changed = true;
  }
  return true;
}

bool Presolve::process_duplicate_rows(bool& changed) {
  // Hash each row's pattern and values normalized by its first active coefficient.
  const auto first_active = [&](Index i) {
    std::size_t p = 0;
    while (p < rows_[i].size() && !col_active_[rows_[i][p].index]) ++p;
    return p;
  };
  std::unordered_map<std::uint64_t, std::vector<Index>> buckets;
  for (Index i = 0; i < m_; ++i) {
    if (!row_active_[i] || row_count_[i] < 2) continue;
    const auto& row = rows_[i];
    const std::size_t p0 = first_active(i);
    const double scale = 1.0 / row[p0].a;
    std::uint64_t h = 1469598103934665603ULL;
    for (std::size_t p = p0; p < row.size(); ++p) {
      if (!col_active_[row[p].index]) continue;
      const auto q = static_cast<std::int64_t>(std::llround(row[p].a * scale * 1e6));
      h = (h ^ static_cast<std::uint64_t>(row[p].index)) * 1099511628211ULL;
      h = (h ^ static_cast<std::uint64_t>(q)) * 1099511628211ULL;
    }
    buckets[h].push_back(i);
  }

  // factor with a_k = factor * a_i, or 0 if the active parts are not parallel.
  const auto parallel = [&](Index i, Index k) {
    if (row_count_[i] != row_count_[k]) return 0.0;
    const auto& row_i = rows_[i];
    const auto& row_k = rows_[k];
    std::size_t p = first_active(i);
    std::size_t q = first_active(k);
    const double factor = row_k[q].a / row_i[p].a;
    if (!std::isfinite(factor) || factor == 0.0) return 0.0;
    while (p < row_i.size() || q < row_k.size()) {
      while (p < row_i.size() && !col_active_[row_i[p].index]) ++p;
      while (q < row_k.size() && !col_active_[row_k[q].index]) ++q;
      const bool end_i = p >= row_i.size();
      const bool end_k = q >= row_k.size();
      if (end_i || end_k) return end_i && end_k ? factor : 0.0;
      if (row_i[p].index != row_k[q].index) return 0.0;
      // Require the same computed ratio and an exact, bitwise product back-check. A tiny
      // residual can distinguish feasible equations with very large solutions.
      if (row_k[q].a / row_i[p].a != factor ||
          std::bit_cast<std::uint64_t>(row_k[q].a) !=
              std::bit_cast<std::uint64_t>(factor * row_i[p].a)) {
        return 0.0;
      }
      ++p;
      ++q;
    }
    return factor;
  };

  for (auto& [hash, rows] : buckets) {
    if (rows.size() < 2) continue;
    for (std::size_t s = 0; s < rows.size(); ++s) {
      const Index i = rows[s];
      if (!row_active_[i]) continue;
      for (std::size_t t = s + 1; t < rows.size(); ++t) {
        const Index k = rows[t];
        if (!row_active_[k]) continue;
        const double factor = parallel(i, k);
        if (factor == 0.0) continue;
        // Row k in terms of row i's activity.
        double lo = factor > 0.0 ? row_lower_[k] / factor : row_upper_[k] / factor;
        double up = factor > 0.0 ? row_upper_[k] / factor : row_lower_[k] / factor;
        if (std::isnan(lo)) lo = -kInf;
        if (std::isnan(up)) up = kInf;
        Record record{.kind = Kind::kDuplicateRow, .row = i, .row2 = k, .value = factor};
        record.flag1 = lo > row_lower_[i];
        record.flag2 = up < row_upper_[i];
        double new_lower = std::max(lo, row_lower_[i]);
        double new_upper = std::min(up, row_upper_[i]);
        if (new_lower > new_upper) {
          if (new_lower > new_upper + tol(new_upper)) return false;
          // Crossed within the tolerance: make the row an equality at the upper bound.
          new_lower = new_upper;
          record.flag1 = record.flag2;
        }
        row_lower_[i] = new_lower;
        row_upper_[i] = new_upper;
        records_.push_back(record);
        remove_row(k);
        ++stats_.duplicate_rows;
        changed = true;
      }
    }
  }
  return true;
}

bool Presolve::process_doubleton_rows(bool& changed) {
  if (!options_.mip) return true;  // Dual restoration for transferred bounds is a separate rule.
  struct Update {
    Index row;
    double coefficient;
    double lower;
    double upper;
    bool had_survivor;
  };
  for (Index i = 0; i < m_; ++i) {
    if (!row_active_[i] || row_count_[i] != 2 || row_lower_[i] != row_upper_[i]) continue;
    Entry pivot{}, survivor{};
    bool first = true;
    for (const auto& e : rows_[i]) {
      if (!col_active_[e.index]) continue;
      if (first) pivot = e;
      else survivor = e;
      first = false;
    }
    const auto prefer = [&](const Entry& e, const Entry& f) {
      if (is_int(e.index) != is_int(f.index)) return !is_int(e.index);
      if (col_count_[e.index] != col_count_[f.index]) {
        return col_count_[e.index] < col_count_[f.index];
      }
      return std::fabs(e.a) > std::fabs(f.a);
    };
    if (prefer(survivor, pivot)) std::swap(pivot, survivor);
    const Index x = pivot.index, y = survivor.index;
    const double a = pivot.a, b = survivor.a, c = row_lower_[i];
    if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(c)) {
      ++stats_.doubleton_rejected_nonfinite;
      continue;
    }
    if (std::fabs(a) < kSubstitutionPivotRatio * std::max(std::fabs(a), std::fabs(b))) {
      ++stats_.doubleton_rejected_pivot;
      continue;
    }
    if (is_int(x)) {
      const auto integral = [](double v) {
        return std::fabs(v) <= kIntegerDataLimit && v == std::trunc(v);
      };
      if (!is_int(y) || !integral(a) || !integral(b) || !integral(c)) {
        ++stats_.doubleton_rejected_integer_data;
        continue;
      }
      const auto ia = static_cast<std::int64_t>(a);
      if (static_cast<std::int64_t>(b) % ia != 0 || static_cast<std::int64_t>(c) % ia != 0) {
        ++stats_.doubleton_rejected_divisibility;
        continue;
      }
    }

    // The lower/upper y endpoints come from these x bounds; signs determine their order.
    const bool decreasing = (a > 0.0) == (b > 0.0);
    const double xlo = decreasing ? upper_[x] : lower_[x];
    const double xup = decreasing ? lower_[x] : upper_[x];
    double lo = std::isfinite(xlo) ? (c - a * xlo) / b : -kInf;
    double up = std::isfinite(xup) ? (c - a * xup) / b : kInf;
    const auto endpoint_error = [&](double bound, double endpoint) {
      if (!std::isfinite(bound)) return 0.0;
      return kBoundRoundoff * ((std::fabs(c) + std::fabs(a * bound)) / std::fabs(b) +
                              std::fabs(endpoint));
    };
    const double lo_error = endpoint_error(xlo, lo), up_error = endpoint_error(xup, up);
    if ((std::isfinite(xlo) && !std::isfinite(lo)) ||
        (std::isfinite(xup) && !std::isfinite(up)) ||
        !std::isfinite(lo_error) || !std::isfinite(up_error)) {
      ++stats_.doubleton_rejected_nonfinite;
      continue;
    }
    // Outward endpoints: the computed ones widened by their rounding error, so they never cut
    // off a y that is feasible in exact arithmetic.
    double outer_lo = std::isfinite(lo) ? std::nextafter(lo - lo_error, -kInf) : lo;
    double outer_up = std::isfinite(up) ? std::nextafter(up + up_error, kInf) : up;
    // Infeasibility is claimed only when the interval stays empty after also allowing x's bound
    // and the row their feasibility tolerances (as set_bounds does); a crossing inside them is
    // skipped below, never collapsed.
    const auto slack = [&](double bound) {
      return (tol(bound) * std::fabs(a) + tol(c)) / std::fabs(b);
    };
    double proof_lo = std::isfinite(lo) ? std::nextafter(outer_lo - slack(xlo), -kInf) : lo;
    double proof_up = std::isfinite(up) ? std::nextafter(outer_up + slack(xup), kInf) : up;
    if (is_int(y)) {
      outer_lo = std::ceil(outer_lo);
      outer_up = std::floor(outer_up);
      proof_lo = std::ceil(proof_lo);
      proof_up = std::floor(proof_up);
    }
    if (std::max(lower_[y], proof_lo) > std::min(upper_[y], proof_up)) return false;
    if (is_int(y)) {
      // Start from the outward integers and move inward: rounding the computed endpoint
      // directly can step past an exactly feasible integer (a true endpoint of 3 computed as
      // 3 + 1 ulp would start at 4).
      lo = std::max(lower_[y], outer_lo);
      up = std::min(upper_[y], outer_up);
      bool ambiguous = false, nonfinite = false;
      const auto check_endpoint = [&](double& v, bool lower) {
        if (!std::isfinite(v)) return;
        // Validate both x bounds; this also covers endpoints retained from y's old bounds.
        for (int attempt = 0; attempt < 2; ++attempt) {
          const double product = b * v;
          const double value = (c - product) / a;
          const double error = kBoundRoundoff * (std::fabs(c) + std::fabs(product)) / std::fabs(a);
          if (!std::isfinite(product) || !std::isfinite(value) || !std::isfinite(error)) {
            nonfinite = true;
            return;
          }
          // When b * v, c - b * v and the division are all exact (zero FMA and two-sum
          // residuals, as with integral data), x(v) is compared without an allowance: otherwise
          // an x landing exactly on its bound, the common integral case, would be ambiguous.
          const double diff = c - product;
          const double minus_part = diff - c;          // Knuth's two-sum of c and -product.
          const double c_part = diff - minus_part;
          const bool exact = std::fma(b, v, -product) == 0.0 &&
                             (c - c_part) + (-product - minus_part) == 0.0 &&
                             std::fma(a, value, -diff) == 0.0;
          // Outside means outside by more than x's feasibility tolerance too: the solver would
          // accept a nearer point, so cutting it off is not proven; that is ambiguous.
          const double allowance = exact ? 0.0 : error;
          const bool inside = exact ? value >= lower_[x] && value <= upper_[x]
                                    : value - lower_[x] > error && upper_[x] - value > error;
          if (inside) return;
          const bool outside = lower_[x] - value > allowance + tol(lower_[x]) ||
                               value - upper_[x] > allowance + tol(upper_[x]);
          if (!outside || attempt == 1) {
            ambiguous = true;
            return;
          }
          const double next = v + (lower ? 1.0 : -1.0);
          if (next == v || !std::isfinite(next)) {
            ambiguous = true;
            return;
          }
          v = next;
        }
      };
      check_endpoint(lo, true);
      check_endpoint(up, false);
      if (nonfinite) {
        ++stats_.doubleton_rejected_nonfinite;
        continue;
      }
      if (ambiguous) {
        ++stats_.doubleton_rejected_integer_bound;
        continue;
      }
    } else {
      lo = std::max(lower_[y], lo);
      up = std::min(upper_[y], up);
    }
    if (lo > up) {
      ++stats_.doubleton_rejected_crossed_bounds;
      continue;
    }

    std::vector<Update> updates;
    bool nonfinite = false, residual = false;
    for (const auto& e : cols_[x]) {
      const Index r = e.index;
      if (r == i || !row_active_[r]) continue;
      double old = 0.0, largest = 0.0;
      bool had_survivor = false;
      for (const auto& f : rows_[r]) {
        if (!col_active_[f.index] || f.index == x) continue;
        if (f.index == y) {
          old = f.a;
          had_survivor = true;
        } else {
          largest = std::max(largest, std::fabs(f.a));
        }
      }
      const double term = e.a * b / a;
      const double coefficient = old - term;
      const double shift = e.a * c / a;
      if (cancelled(old, term, coefficient)) residual = true;
      const double rl = std::isfinite(row_lower_[r]) ? row_lower_[r] - shift : row_lower_[r];
      const double ru = std::isfinite(row_upper_[r]) ? row_upper_[r] - shift : row_upper_[r];
      nonfinite = nonfinite || !std::isfinite(coefficient) || !std::isfinite(shift) ||
          (std::isfinite(row_lower_[r]) && !std::isfinite(rl)) ||
          (std::isfinite(row_upper_[r]) && !std::isfinite(ru));
      largest = std::max(largest, std::fabs(coefficient));
      if (coefficient != 0.0 && std::fabs(coefficient) < kSubstitutionResidualRatio * largest) {
        residual = true;
      }
      // Check all retained nonzeros, not only the changed survivor coefficient.
      for (const auto& f : rows_[r]) {
        if (col_active_[f.index] && f.index != x && f.index != y &&
            std::fabs(f.a) < kSubstitutionResidualRatio * largest) residual = true;
      }
      updates.push_back({r, coefficient, rl, ru, had_survivor});
    }
    const double new_cost = cost_[y] - cost_[x] * b / a;
    const double new_offset = offset_ + cost_[x] * c / a;
    if (nonfinite || !std::isfinite(new_cost) || !std::isfinite(new_offset)) {
      ++stats_.doubleton_rejected_nonfinite;
      continue;
    }
    if (residual) {
      ++stats_.doubleton_rejected_residual;
      continue;
    }

    Record record{.kind = Kind::kDoubleton, .row = i, .col = x, .a = a,
                  .cost = cost_[x], .old_lower = lower_[x], .old_upper = upper_[x],
                  .rhs_lower = c};
    record.survivor = y;
    record.b = b;
    record.survivor_lower = lower_[y];
    record.survivor_upper = upper_[y];
    record.survivor_cost = new_cost;
    snapshot_column(x, record.first, record.last);
    // Stage the post-substitution y snapshot too, before modifying either adjacency list.
    std::vector<Entry> new_column;
    for (const auto& e : cols_[y]) {
      if (e.index != i && row_active_[e.index]) new_column.push_back(e);
    }
    for (const auto& u : updates) {
      const auto it = std::lower_bound(new_column.begin(), new_column.end(), u.row,
          [](const Entry& e, Index r) { return e.index < r; });
      if (it != new_column.end() && it->index == u.row) {
        if (u.coefficient == 0.0) new_column.erase(it);
        else it->a = u.coefficient;
      } else if (u.coefficient != 0.0) {
        new_column.insert(it, {u.row, u.coefficient});
      }
    }
    record.survivor_first = column_terms_.size();
    column_terms_.insert(column_terms_.end(), new_column.begin(), new_column.end());
    record.survivor_last = column_terms_.size();
    records_.push_back(record);

    for (const auto& u : updates) {
      if (u.had_survivor) change_coefficient(u.row, y, u.coefficient);
      else if (u.coefficient != 0.0) add_entry(u.row, y, u.coefficient);
      remove_entry(u.row, x);
      row_lower_[u.row] = u.lower;
      row_upper_[u.row] = u.upper;
    }
    if (lo != lower_[y] || up != upper_[y]) ++stats_.bounds_tightened;
    lower_[y] = lo;
    upper_[y] = up;
    cost_[y] = new_cost;
    offset_ = new_offset;
    remove_row(i);
    remove_col(x);
    if (is_int(x)) ++stats_.doubleton_integer_cols;
    else ++stats_.doubleton_continuous_cols;
    changed = true;
  }
  return true;
}

bool Presolve::process_substitutions(bool& changed) {
  if (!options_.substitute) return true;
  struct Update {
    Index row;
    std::vector<Entry> entries;
    double lower;
    double upper;
  };
  for (Index i = 0; i < m_; ++i) {
    if (!row_active_[i] || row_count_[i] < 2 || row_lower_[i] != row_upper_[i]) continue;
    if (row_count_[i] > kSubstitutionRowLength) {
      ++stats_.substitution_rejected_row_length;
      continue;
    }
    const double rhs = row_lower_[i];
    std::vector<Entry> equation;
    double largest = 0.0;
    for (const auto& e : rows_[i]) {
      if (!col_active_[e.index]) continue;
      equation.push_back(e);
      largest = std::max(largest, std::fabs(e.a));
    }
    Entry pivot{-1, 0.0};
    const auto integer = [&](Index j) { return options_.mip && is_int(j); };
    const auto prefer = [&](const Entry& e, const Entry& f) {
      if (integer(e.index) != integer(f.index)) return !integer(e.index);
      // Row length is shared, so comparing col lengths compares the Markowitz counts.
      if (col_count_[e.index] != col_count_[f.index]) {
        return col_count_[e.index] < col_count_[f.index];
      }
      return std::fabs(e.a) > std::fabs(f.a);
    };
    for (const auto& e : equation) {
      if (col_count_[e.index] < 2) continue;  // C1 owns singletons.
      if (col_count_[e.index] - 1 > kSubstitutionColLength) {
        ++stats_.substitution_rejected_col_length;
        continue;
      }
      if (!std::isfinite(e.a) || !std::isfinite(largest) || !std::isfinite(rhs)) {
        ++stats_.substitution_rejected_nonfinite;
        continue;
      }
      if (std::fabs(e.a) < kSubstitutionPivotRatio * largest) {
        ++stats_.substitution_rejected_pivot;
        continue;
      }
      if (integer(e.index)) {
        const auto integral = [](double v) {
          return std::fabs(v) <= kIntegerDataLimit && v == std::trunc(v);
        };
        bool valid = integral(e.a) && integral(rhs);
        for (const auto& f : equation) {
          if (f.index != e.index) valid = valid && is_int(f.index) && integral(f.a);
        }
        if (valid) {
          const auto a = static_cast<std::int64_t>(e.a);
          valid = static_cast<std::int64_t>(rhs) % a == 0;
          for (const auto& f : equation) {
            if (f.index != e.index) valid = valid && static_cast<std::int64_t>(f.a) % a == 0;
          }
        }
        if (!valid) {
          ++stats_.substitution_rejected_integer;
          continue;
        }
      }
      // Bounds can change between equations and passes: prove this candidate afresh.
      bool tolerance_only = false;
      if (!implied_free(i, e.index, e.a, tolerance_only)) {
        if (tolerance_only) ++stats_.substitution_rejected_by_tolerance;
        else ++stats_.substitution_rejected_not_implied_free;
        continue;
      }
      if (pivot.index < 0 || prefer(e, pivot)) pivot = e;
    }
    if (pivot.index < 0) continue;
    const Index x = pivot.index;
    const double a = pivot.a;
    std::vector<Update> updates;
    Index fill = 0;
    bool nonfinite = false, residual = false;
    for (const auto& e : cols_[x]) {
      const Index r = e.index;
      if (r == i || !row_active_[r]) continue;
      Update update{r, {}, row_lower_[r], row_upper_[r]};
      for (const auto& f : rows_[r]) {
        if (col_active_[f.index] && f.index != x) update.entries.push_back(f);
      }
      const double factor = e.a / a;
      const double shift = factor * rhs;
      nonfinite = nonfinite || !std::isfinite(factor) || !std::isfinite(shift);
      if (std::isfinite(update.lower)) update.lower -= shift;
      if (std::isfinite(update.upper)) update.upper -= shift;
      nonfinite = nonfinite || (std::isfinite(row_lower_[r]) && !std::isfinite(update.lower)) ||
          (std::isfinite(row_upper_[r]) && !std::isfinite(update.upper));
      for (const auto& f : equation) {
        if (f.index == x) continue;
        const auto it = std::lower_bound(update.entries.begin(), update.entries.end(), f.index,
            [](const Entry& g, Index j) { return g.index < j; });
        const bool had_entry = it != update.entries.end() && it->index == f.index;
        const double old = had_entry ? it->a : 0.0;
        const double term = factor * f.a;
        const double coefficient = old - term;
        nonfinite = nonfinite || !std::isfinite(coefficient);
        if (cancelled(old, term, coefficient)) residual = true;
        if (had_entry) {
          if (coefficient == 0.0) {
            update.entries.erase(it);
            --fill;
          } else {
            it->a = coefficient;
          }
        } else if (coefficient != 0.0) {
          update.entries.insert(it, {f.index, coefficient});
          ++fill;
        }
      }
      double row_largest = 0.0;
      for (const auto& f : update.entries) row_largest = std::max(row_largest, std::fabs(f.a));
      for (const auto& f : update.entries) {
        if (std::fabs(f.a) < kSubstitutionResidualRatio * row_largest) residual = true;
      }
      updates.push_back(std::move(update));
    }
    std::vector<Entry> costs;
    const double cost_factor = cost_[x] / a;
    const double new_offset = offset_ + cost_factor * rhs;
    nonfinite = nonfinite || !std::isfinite(cost_factor) || !std::isfinite(new_offset);
    for (const auto& e : equation) {
      if (e.index == x) continue;
      const double cost = cost_[e.index] - cost_factor * e.a;
      costs.push_back({e.index, cost});
      nonfinite = nonfinite || !std::isfinite(cost);
    }
    if (nonfinite) {
      ++stats_.substitution_rejected_nonfinite;
      continue;
    }
    if (residual) {
      ++stats_.substitution_rejected_residual;
      continue;
    }
    // Net growth excludes the removed pivot entries and the removed equation.
    if (fill > kSubstitutionFillLimit) {
      ++stats_.substitution_rejected_fill;
      continue;
    }
    Record record{.kind = Kind::kSubstitution, .row = i, .col = x, .a = a,
                  .cost = cost_[x], .rhs_lower = rhs};
    record.first = terms_.size();
    for (const auto& e : equation) {
      if (e.index != x) terms_.push_back({e.index, e.a, cost_[e.index], 0.0, false});
    }
    record.last = terms_.size();
    record.column_first = column_terms_.size();
    for (const auto& e : cols_[x]) {
      if (e.index != i && row_active_[e.index]) column_terms_.push_back(e);
    }
    record.column_last = column_terms_.size();
    records_.push_back(record);  // Save both original snapshots before any mutation.
    for (const auto& u : updates) {
      for (const auto& e : equation) {
        if (e.index == x) continue;
        const auto it = std::lower_bound(u.entries.begin(), u.entries.end(), e.index,
            [](const Entry& f, Index j) { return f.index < j; });
        const double coefficient = it != u.entries.end() && it->index == e.index ? it->a : 0.0;
        if (!change_coefficient(u.row, e.index, coefficient) && coefficient != 0.0) {
          add_entry(u.row, e.index, coefficient);
        }
      }
      remove_entry(u.row, x);
      row_lower_[u.row] = u.lower;
      row_upper_[u.row] = u.upper;
    }
    for (const auto& e : costs) cost_[e.index] = e.a;
    offset_ = new_offset;
    remove_row(i);
    remove_col(x);
    if (integer(x)) ++stats_.substitution_integer_cols;
    else ++stats_.substitution_continuous_cols;
    changed = true;
  }
  return true;
}

PresolveStatus Presolve::run() {
  for (Index j = 0; j < n_; ++j) {
    if (is_int(j) && !set_bounds(j, lower_[j], upper_[j])) return PresolveStatus::kInfeasible;
  }
  for (int pass = 0; pass < options_.max_passes; ++pass) {
    bool changed = false;
    if (!process_rows(changed) || !process_cols(changed) || !process_singleton_cols(changed) ||
        !process_duplicate_rows(changed) || !process_doubleton_rows(changed) ||
        !process_substitutions(changed)) {
      return PresolveStatus::kInfeasible;
    }
    ++stats_.passes;
    if (!changed) break;
  }
  build_reduced();
  return PresolveStatus::kReduced;
}

void Presolve::build_reduced() {
  std::vector<Index> new_col(static_cast<std::size_t>(n_), -1);
  std::vector<Index> new_row(static_cast<std::size_t>(m_), -1);
  col_map_.clear();
  row_map_.clear();
  for (Index j = 0; j < n_; ++j) {
    if (!col_active_[j]) continue;
    new_col[j] = static_cast<Index>(col_map_.size());
    col_map_.push_back(j);
  }
  for (Index i = 0; i < m_; ++i) {
    if (!row_active_[i]) continue;
    new_row[i] = static_cast<Index>(row_map_.size());
    row_map_.push_back(i);
  }
  reduced_ = Model{};
  reduced_.name = model_.name;
  reduced_.sense = model_.sense;
  reduced_.obj_offset = sense_ * offset_;
  for (const Index j : col_map_) {
    reduced_.obj.push_back(sense_ * cost_[j]);
    reduced_.col_lower.push_back(lower_[j]);
    reduced_.col_upper.push_back(upper_[j]);
    reduced_.col_type.push_back(model_.col_type[j]);
    if (!model_.col_names.empty()) reduced_.col_names.push_back(model_.col_names[j]);
  }
  for (const Index i : row_map_) {
    reduced_.row_lower.push_back(row_lower_[i]);
    reduced_.row_upper.push_back(row_upper_[i]);
    if (!model_.row_names.empty()) reduced_.row_names.push_back(model_.row_names[i]);
  }
  std::vector<Triplet> triplets;
  std::size_t nnz = 0;
  for (const Index j : col_map_) nnz += static_cast<std::size_t>(col_count_[j]);
  triplets.reserve(nnz);
  for (Index j = 0; j < n_; ++j) {
    if (new_col[j] < 0) continue;
    for (const auto& e : cols_[j]) {
      if (new_row[e.index] >= 0) triplets.push_back({new_row[e.index], new_col[j], e.a});
    }
  }
  reduced_.A = SparseMatrix::from_triplets(static_cast<Index>(row_map_.size()),
                                           static_cast<Index>(col_map_.size()),
                                           std::move(triplets));
}

void Presolve::postsolve(const std::vector<double>& reduced_x,
                         const std::vector<double>& reduced_y, std::vector<double>& x,
                         std::vector<double>& y) const {
  const bool duals = reduced_y.size() == row_map_.size() && !options_.mip;
  x.assign(static_cast<std::size_t>(n_), 0.0);
  y.assign(static_cast<std::size_t>(m_), 0.0);
  for (std::size_t k = 0; k < col_map_.size(); ++k) x[col_map_[k]] = reduced_x[k];
  if (duals) {
    for (std::size_t k = 0; k < row_map_.size(); ++k) y[row_map_[k]] = sense_ * reduced_y[k];
  }

  // Column snapshots preserve the coefficients and summation order at each reduction. Rows
  // removed before that reduction still have zero duals, so their entries were omitted.
  const auto reduced_cost = [&](std::size_t first, std::size_t last, double cost) {
    double d = cost;
    for (std::size_t t = first; t < last; ++t) {
      const Entry& e = column_terms_[t];
      d -= e.a * y[e.index];
    }
    return d;
  };
  const auto at = [&](double v, double bound) {
    return std::isfinite(bound) && std::fabs(v - bound) <= tol(bound);
  };

  for (auto it = records_.rbegin(); it != records_.rend(); ++it) {
    const Record& r = *it;
    switch (r.kind) {
      case Kind::kFixCol:
        x[r.col] = r.value;
        break;
      case Kind::kRemoveRow:
        y[r.row] = 0.0;
        break;
      case Kind::kSingletonRow: {
        if (!duals) break;
        // If the column rests on a bound that came from this row, the row takes its reduced
        // cost.
        const double d = reduced_cost(r.first, r.last, r.cost);
        const double v = x[r.col];
        if ((d > 0.0 && !at(v, r.old_lower)) || (d < 0.0 && !at(v, r.old_upper))) {
          y[r.row] = d / r.a;
        }
        break;
      }
      case Kind::kForcingRow: {
        for (std::size_t t = r.first; t < r.last; ++t) x[terms_[t].col] = terms_[t].value;
        if (!duals) break;
        // The row is at its upper bound (minimum activity) or lower bound (maximum activity);
        // choose y of that sign so every column's reduced cost fits the bound it sits on.
        double yi = 0.0;
        for (std::size_t t = r.first; t < r.last; ++t) {
          const Term& term = terms_[t];
          if (!term.boxed) continue;
          const double ratio = reduced_cost(term.column_first, term.column_last, term.cost) / term.a;
          yi = r.flag1 ? std::min(yi, ratio) : std::max(yi, ratio);
        }
        y[r.row] = yi;
        break;
      }
      case Kind::kSubstituteCol: {
        double s = 0.0;
        for (std::size_t t = r.first; t < r.last; ++t) s += terms_[t].a * x[terms_[t].col];
        x[r.col] = (r.rhs_lower - s) / r.a;
        if (duals) y[r.row] = r.cost / r.a;
        break;
      }
      case Kind::kDoubleton:
        x[r.col] = (r.rhs_lower - r.b * x[r.survivor]) / r.a;
        break;
      case Kind::kSubstitution: {
        double s = 0.0;
        for (std::size_t t = r.first; t < r.last; ++t) s += terms_[t].a * x[terms_[t].col];
        x[r.col] = (r.rhs_lower - s) / r.a;
        // y_e = (c_x - sum a_rx y_r) / a makes d_x = 0; transferred costs and row
        // coefficients cancel in d_k = c_k - sum a_rk y_r - a_ek y_e for every other k.
        if (duals) y[r.row] = reduced_cost(r.column_first, r.column_last, r.cost) / r.a;
        break;
      }
      case Kind::kSlackCol: {
        double s = 0.0;
        for (std::size_t t = r.first; t < r.last; ++t) s += terms_[t].a * x[terms_[t].col];
        const double target = std::clamp(s, r.rhs_lower, r.rhs_upper);
        x[r.col] = (target - s) / r.a;
        if (duals) y[r.row] = 0.0;
        break;
      }
      case Kind::kDuplicateRow: {
        if (!duals) break;
        const double yi = y[r.row];
        if ((yi > 0.0 && r.flag1) || (yi < 0.0 && r.flag2)) {
          y[r.row2] = yi / r.value;
          y[r.row] = 0.0;
        }
        break;
      }
    }
  }
  if (duals) {
    for (double& v : y) v *= sense_;
  }
}

}  // namespace samaya
