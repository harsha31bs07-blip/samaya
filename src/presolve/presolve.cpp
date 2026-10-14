#include "presolve/presolve.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <unordered_map>

namespace samaya {

namespace {

// A free column singleton is substituted out only if its coefficient is not tiny relative to the
// rest of its row (as in threshold pivoting), to keep the substitution stable.
constexpr double kSubstitutionPivotRatio = 0.01;
// Relative tolerance for two rows being parallel.
constexpr double kParallelTol = 1e-12;
// Prove redundant bounds with an outward allowance, rather than accepting a near miss.
constexpr double kImpliedFreeErrorTol = 1e-12;

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
    while (p < row_i.size() || q < row_k.size()) {
      while (p < row_i.size() && !col_active_[row_i[p].index]) ++p;
      while (q < row_k.size() && !col_active_[row_k[q].index]) ++q;
      const bool end_i = p >= row_i.size();
      const bool end_k = q >= row_k.size();
      if (end_i || end_k) return end_i && end_k ? factor : 0.0;
      if (row_i[p].index != row_k[q].index) return 0.0;
      if (std::fabs(row_k[q].a - factor * row_i[p].a) > kParallelTol * std::fabs(row_k[q].a)) {
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

PresolveStatus Presolve::run() {
  for (Index j = 0; j < n_; ++j) {
    if (is_int(j) && !set_bounds(j, lower_[j], upper_[j])) return PresolveStatus::kInfeasible;
  }
  for (int pass = 0; pass < options_.max_passes; ++pass) {
    bool changed = false;
    if (!process_rows(changed) || !process_cols(changed) || !process_singleton_cols(changed) ||
        !process_duplicate_rows(changed)) {
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
