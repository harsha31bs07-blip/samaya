#include "linalg/basis_factor.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>

namespace samaya {
namespace {

// Doubly linked lists of items bucketed by their current nonzero count.
class CountLists {
 public:
  void init(Index items, Index max_count) {
    head_.assign(static_cast<std::size_t>(max_count) + 1, -1);
    next_.assign(static_cast<std::size_t>(items), -1);
    prev_.assign(static_cast<std::size_t>(items), -1);
    count_.assign(static_cast<std::size_t>(items), -1);
  }
  Index head(Index count) const { return head_[count]; }
  Index next(Index item) const { return next_[item]; }
  bool contains(Index item) const { return count_[item] >= 0; }

  void insert(Index item, Index count) {
    count_[item] = count;
    prev_[item] = -1;
    next_[item] = head_[count];
    if (head_[count] != -1) prev_[head_[count]] = item;
    head_[count] = item;
  }
  void remove(Index item) {
    const Index c = count_[item];
    if (c < 0) return;
    if (prev_[item] != -1) {
      next_[prev_[item]] = next_[item];
    } else {
      head_[c] = next_[item];
    }
    if (next_[item] != -1) prev_[next_[item]] = prev_[item];
    count_[item] = -1;
  }
  void update(Index item, Index count) {
    if (count_[item] == count) return;
    remove(item);
    insert(item, count);
  }

 private:
  std::vector<Index> head_, next_, prev_, count_;
};

void erase_index(std::vector<Index>& list, Index value) {
  const auto it = std::find(list.begin(), list.end(), value);
  assert(it != list.end());
  *it = list.back();
  list.pop_back();
}

}  // namespace

BasisFactor::BasisFactor(const SparseMatrix& A, Options options)
    : A_(A), options_(options), m_(A.rows()), n_(A.cols()) {
  work_.assign(static_cast<std::size_t>(m_), 0.0);
  update_work_.assign(static_cast<std::size_t>(m_), 0.0);
}

void BasisFactor::clear() {
  l_pivot_row_.clear();
  l_start_.assign(1, 0);
  l_index_.clear();
  l_value_.clear();
  r_pivot_row_.clear();
  r_start_.assign(1, 0);
  r_index_.clear();
  r_value_.clear();
  const auto m = static_cast<std::size_t>(m_);
  prow_.assign(m, -1);
  pos_of_row_.assign(m, -1);
  diag_.assign(m, 0.0);
  ucol_.assign(m, {});
  urow_.assign(m, {});
  order_.clear();
  order_slot_.assign(m, -1);
  u_nnz_ = 0;
  num_updates_ = 0;
  extra_solve_work_ = 0.0;
  valid_ = false;
  singular_positions_.clear();
  unpivoted_rows_.clear();
}

void BasisFactor::erase_position(std::vector<Entry>& list, Index position) {
  for (std::size_t k = 0; k < list.size(); ++k) {
    if (list[k].first == position) {
      list[k] = list.back();
      list.pop_back();
      return;
    }
  }
  assert(false && "entry not found");
}

BasisFactor::RebuildMeasurement::RebuildMeasurement(BasisFactor& f, RebuildReason reason)
    : factor(f) {
  if (!factor.options_.collect_stats) return;
  start = std::chrono::steady_clock::now();
  ++factor.stats_.rebuilds;
  ++factor.stats_.reasons[static_cast<std::size_t>(reason)];
  factor.stats_.old_l_nnz.add(static_cast<long long>(factor.l_index_.size()));
  factor.stats_.old_u_nnz.add(factor.u_nnz_ + static_cast<long long>(
      std::count_if(factor.diag_.begin(), factor.diag_.end(), [](double v) { return v != 0.0; })));
  factor.stats_.old_r_nnz.add(static_cast<long long>(factor.r_index_.size()));
}

BasisFactor::RebuildMeasurement::~RebuildMeasurement() {
  if (!factor.options_.collect_stats) return;
  const double seconds = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - start).count();
  factor.stats_.rebuild_seconds += seconds;
  factor.stats_.max_rebuild_seconds = std::max(factor.stats_.max_rebuild_seconds, seconds);
  factor.stats_.l_nnz.add(static_cast<long long>(factor.l_index_.size()));
  factor.stats_.u_nnz.add(factor.u_nnz_ + static_cast<long long>(factor.order_.size()));
  factor.stats_.r_nnz.add(static_cast<long long>(factor.r_index_.size()));
}

Index BasisFactor::factorize(std::span<const Index> basic, RebuildReason reason) {
  assert(basic.size() == static_cast<std::size_t>(m_));
  const RebuildMeasurement measurement(*this, reason);
  clear();
  const Index m = m_;
  if (m == 0) {
    valid_ = true;
    return 0;
  }

  // Active submatrix, by column (with values) and by row (pattern only).
  std::vector<std::vector<Entry>> col(static_cast<std::size_t>(m));
  std::vector<std::vector<Index>> row(static_cast<std::size_t>(m));
  const auto a_start = A_.col_start();
  const auto a_index = A_.row_index();
  const auto a_value = A_.values();
  for (Index k = 0; k < m; ++k) {
    const Index j = basic[k];
    if (j < n_) {
      for (NnzIndex p = a_start[j]; p < a_start[j + 1]; ++p) {
        col[k].push_back({a_index[p], a_value[p]});
        row[a_index[p]].push_back(k);
      }
    } else {
      col[k].push_back({j - n_, -1.0});
      row[j - n_].push_back(k);
    }
  }

  CountLists col_lists;
  CountLists row_lists;
  col_lists.init(m, m);
  row_lists.init(m, m);
  for (Index k = 0; k < m; ++k) col_lists.insert(k, static_cast<Index>(col[k].size()));
  for (Index i = 0; i < m; ++i) row_lists.insert(i, static_cast<Index>(row[i].size()));

  std::vector<int> mark(static_cast<std::size_t>(m), 0);
  std::vector<Index> where(static_cast<std::size_t>(m), 0);
  int stamp = 0;

  const auto column_max = [&](Index k) {
    double mx = 0.0;
    for (const Entry& e : col[k]) mx = std::max(mx, std::fabs(e.second));
    return mx;
  };
  const auto value_at = [&](Index k, Index i) {
    for (const Entry& e : col[k]) {
      if (e.first == i) return e.second;
    }
    return 0.0;
  };

  Index pivots = 0;
  while (pivots + static_cast<Index>(singular_positions_.size()) < m) {
    // Columns without entries are linearly dependent on the pivoted ones.
    for (Index k = col_lists.head(0); k != -1; k = col_lists.head(0)) {
      col_lists.remove(k);
      singular_positions_.push_back(k);
    }
    if (pivots + static_cast<Index>(singular_positions_.size()) == m) break;

    // Markowitz search, lowest counts first; singletons have cost 0 and are taken at once.
    Index best_row = -1;
    Index best_col = -1;
    double best_value = 0.0;
    long long best_cost = std::numeric_limits<long long>::max();
    int searched = 0;
    const auto consider = [&](Index i, Index k, double v, double colmax) {
      if (std::fabs(v) <= options_.abs_pivot_tol) return;
      if (std::fabs(v) < options_.pivot_threshold * colmax) return;
      const long long cost = static_cast<long long>(col[k].size() - 1) *
                             static_cast<long long>(row[i].size() - 1);
      if (cost < best_cost || (cost == best_cost && std::fabs(v) > std::fabs(best_value))) {
        best_cost = cost;
        best_row = i;
        best_col = k;
        best_value = v;
      }
    };
    for (Index c = 1; c <= m; ++c) {
      bool done = false;
      for (Index k = col_lists.head(c); k != -1 && !done; k = col_lists.next(k)) {
        const double colmax = column_max(k);
        for (const Entry& e : col[k]) consider(e.first, k, e.second, colmax);
        if (best_col != -1 && (++searched >= options_.search_limit || best_cost == 0)) done = true;
      }
      for (Index i = row_lists.head(c); i != -1 && !done; i = row_lists.next(i)) {
        for (Index k : row[i]) consider(i, k, value_at(k, i), column_max(k));
        if (best_col != -1 && (++searched >= options_.search_limit || best_cost == 0)) done = true;
      }
      if (done) break;
    }

    if (best_col == -1) {
      // No acceptable pivot anywhere: the remaining columns are numerically dependent.
      for (Index k = 0; k < m; ++k) {
        if (col_lists.contains(k)) {
          col_lists.remove(k);
          singular_positions_.push_back(k);
        }
      }
      break;
    }

    // Eliminate with pivot (p, q).
    const Index p = best_row;
    const Index q = best_col;
    const double pivot = best_value;
    col_lists.remove(q);
    row_lists.remove(p);

    const NnzIndex l_begin = static_cast<NnzIndex>(l_index_.size());
    for (const Entry& e : col[q]) {
      if (e.first == p) continue;
      l_index_.push_back(e.first);
      l_value_.push_back(e.second / pivot);
      erase_index(row[e.first], q);
    }
    const NnzIndex l_end = static_cast<NnzIndex>(l_index_.size());
    if (l_end > l_begin) {
      l_pivot_row_.push_back(p);
      l_start_.push_back(l_end);
    }

    prow_[q] = p;
    diag_[q] = pivot;
    order_slot_[q] = static_cast<Index>(order_.size());
    order_.push_back(q);

    for (Index k : row[p]) {
      if (k == q) continue;
      // Move a_pk from the active column into U.
      double a_pk = 0.0;
      for (std::size_t t = 0; t < col[k].size(); ++t) {
        if (col[k][t].first == p) {
          a_pk = col[k][t].second;
          col[k][t] = col[k].back();
          col[k].pop_back();
          break;
        }
      }
      ucol_[k].push_back({q, a_pk});
      urow_[q].push_back({k, a_pk});
      ++u_nnz_;
      if (l_end == l_begin) {
        col_lists.update(k, static_cast<Index>(col[k].size()));
        continue;
      }

      // col_k -= a_pk * l, with fill-in.
      ++stamp;
      for (std::size_t t = 0; t < col[k].size(); ++t) {
        mark[col[k][t].first] = stamp;
        where[col[k][t].first] = static_cast<Index>(t);
      }
      for (NnzIndex t = l_begin; t < l_end; ++t) {
        const Index i = l_index_[t];
        const double delta = -l_value_[t] * a_pk;
        if (mark[i] == stamp) {
          col[k][where[i]].second += delta;
        } else {
          mark[i] = stamp;
          where[i] = static_cast<Index>(col[k].size());
          col[k].push_back({i, delta});
          row[i].push_back(k);
        }
      }
      // Drop entries that cancelled out.
      for (NnzIndex t = l_begin; t < l_end; ++t) {
        const Index i = l_index_[t];
        const Index w = where[i];
        if (std::fabs(col[k][w].second) > options_.drop_tol) continue;
        const Index last = static_cast<Index>(col[k].size()) - 1;
        if (w != last) {
          col[k][w] = col[k][last];
          where[col[k][w].first] = w;
        }
        col[k].pop_back();
        erase_index(row[i], k);
      }
      col_lists.update(k, static_cast<Index>(col[k].size()));
    }
    for (NnzIndex t = l_begin; t < l_end; ++t) {
      const Index i = l_index_[t];
      row_lists.update(i, static_cast<Index>(row[i].size()));
    }
    row[p].clear();
    col[q].clear();
    ++pivots;
  }

  if (!singular_positions_.empty()) {
    for (Index i = 0; i < m; ++i) {
      if (row_lists.contains(i)) unpivoted_rows_.push_back(i);
    }
    assert(unpivoted_rows_.size() == singular_positions_.size());
    return static_cast<Index>(singular_positions_.size());
  }

  for (Index k = 0; k < m; ++k) pos_of_row_[prow_[k]] = k;
  initial_nnz_ = factor_nnz();
  valid_ = true;
  return 0;
}

NnzIndex BasisFactor::factor_nnz() const {
  return static_cast<NnzIndex>(l_index_.size() + r_index_.size()) + u_nnz_ + m_;
}

bool BasisFactor::should_refactor() const {
  return num_updates_ >= options_.max_updates || factor_nnz() > 3 * initial_nnz_ + m_ ||
         update_work_exceeds_rebuild();
}

bool BasisFactor::update_work_exceeds_rebuild() const {
  return options_.rebuild_cost > 0.0 &&
         extra_solve_work_ >= options_.rebuild_cost * static_cast<double>(initial_nnz_);
}

BasisFactor::RebuildReason BasisFactor::refactor_reason() const {
  // Match the short-circuit priority of should_refactor() when several rules fire.
  if (num_updates_ >= options_.max_updates) return RebuildReason::kInterval;
  if (factor_nnz() > 3 * initial_nnz_ + m_) return RebuildReason::kFillLimit;
  return RebuildReason::kUpdateWork;
}

void BasisFactor::ftran(std::vector<double>& x, std::vector<double>* spike,
                        MeasurementCount* nonzeros) const {
  assert(valid_);
  long long l_entries = 0;
  for (std::size_t e = 0; e < l_pivot_row_.size(); ++e) {
    const double xp = x[l_pivot_row_[e]];
    if (xp == 0.0) continue;
    l_entries += l_start_[e + 1] - l_start_[e];
    for (NnzIndex t = l_start_[e]; t < l_start_[e + 1]; ++t) x[l_index_[t]] -= l_value_[t] * xp;
  }
  for (std::size_t e = 0; e < r_pivot_row_.size(); ++e) {
    double sum = 0.0;
    for (NnzIndex t = r_start_[e]; t < r_start_[e + 1]; ++t) sum += r_value_[t] * x[r_index_[t]];
    x[r_pivot_row_[e]] -= sum;
  }
  if (spike != nullptr) *spike = x;

  std::vector<double>& y = work_;
  long long nnz = 0;
  long long u_entries = 0;
  for (std::size_t slot = order_.size(); slot-- > 0;) {
    const Index r = order_[slot];
    if (r < 0) continue;
    const double xp = x[prow_[r]];
    if (xp == 0.0) {
      y[r] = 0.0;
      continue;
    }
    const double yr = xp / diag_[r];
    y[r] = yr;
    if (nonzeros != nullptr && yr != 0.0) ++nnz;
    u_entries += static_cast<long long>(ucol_[r].size());
    for (const Entry& e : ucol_[r]) x[prow_[e.first]] -= e.second * yr;
  }
  x.swap(y);
  if (nonzeros != nullptr) nonzeros->add(nnz);
  if (options_.collect_stats) {
    ++stats_.ftran_calls;
    stats_.ftran_l_etas += static_cast<long long>(l_pivot_row_.size());
    stats_.ftran_l_entries += l_entries;
    stats_.ftran_r_entries += static_cast<long long>(r_index_.size());
    stats_.ftran_u_slots += static_cast<long long>(order_.size());
    stats_.ftran_u_entries += u_entries;
  }
}

void BasisFactor::btran(std::vector<double>& c) const {
  assert(valid_);
  std::vector<double>& w = work_;
  long long u_entries = 0;
  long long r_entries = 0;
  for (const Index r : order_) {
    if (r < 0) continue;
    const double val = c[r];
    if (val == 0.0) {
      w[prow_[r]] = 0.0;
      continue;
    }
    const double z = val / diag_[r];
    w[prow_[r]] = z;
    u_entries += static_cast<long long>(urow_[r].size());
    for (const Entry& e : urow_[r]) c[e.first] -= e.second * z;
  }
  for (std::size_t e = r_pivot_row_.size(); e-- > 0;) {
    const double wp = w[r_pivot_row_[e]];
    if (wp == 0.0) continue;
    r_entries += r_start_[e + 1] - r_start_[e];
    for (NnzIndex t = r_start_[e]; t < r_start_[e + 1]; ++t) w[r_index_[t]] -= r_value_[t] * wp;
  }
  if (options_.collect_stats) {
    ++stats_.btran_calls;
    stats_.btran_u_slots += static_cast<long long>(order_.size());
    stats_.btran_u_entries += u_entries;
    stats_.btran_r_entries += r_entries;
    stats_.btran_l_entries += static_cast<long long>(l_index_.size());
  }
  for (std::size_t e = l_pivot_row_.size(); e-- > 0;) {
    double sum = 0.0;
    for (NnzIndex t = l_start_[e]; t < l_start_[e + 1]; ++t) sum += l_value_[t] * w[l_index_[t]];
    w[l_pivot_row_[e]] -= sum;
  }
  c.swap(w);
}

bool BasisFactor::update(Index r, const std::vector<double>& spike) {
  assert(valid_);
  valid_ = false;  // Until the update completes successfully.

  // Remove the outgoing column from U.
  for (const Entry& e : ucol_[r]) erase_position(urow_[e.first], r);
  u_nnz_ -= static_cast<NnzIndex>(ucol_[r].size());
  ucol_[r].clear();

  // Row r of U now has off-diagonal entries in columns that follow it in the pivot order.
  // Moving pivot r to the end requires eliminating them with a row eta.
  std::vector<double>& w = update_work_;
  for (const Entry& e : urow_[r]) {
    w[e.first] = e.second;
    erase_position(ucol_[e.first], r);
  }
  u_nnz_ -= static_cast<NnzIndex>(urow_[r].size());
  urow_[r].clear();

  const Index p = prow_[r];
  const auto slot_r = static_cast<std::size_t>(order_slot_[r]);
  const auto eta_begin = static_cast<NnzIndex>(r_index_.size());
  for (std::size_t slot = slot_r + 1; slot < order_.size(); ++slot) {
    const Index k = order_[slot];
    if (k < 0) continue;
    const double wk = w[k];
    if (wk == 0.0) continue;
    w[k] = 0.0;
    const double mult = wk / diag_[k];
    r_index_.push_back(prow_[k]);
    r_value_.push_back(mult);
    for (const Entry& e : urow_[k]) w[e.first] -= mult * e.second;
  }
  const auto eta_end = static_cast<NnzIndex>(r_index_.size());

  // Apply the new row eta to the spike to get the replacement column of U.
  double diag = spike[p];
  for (NnzIndex t = eta_begin; t < eta_end; ++t) diag -= r_value_[t] * spike[r_index_[t]];
  double spike_max = 0.0;
  for (const double v : spike) spike_max = std::max(spike_max, std::fabs(v));
  if (std::fabs(diag) <= options_.abs_pivot_tol || std::fabs(diag) < 1e-14 * spike_max) {
    return false;
  }
  if (eta_end > eta_begin) {
    r_pivot_row_.push_back(p);
    r_start_.push_back(eta_end);
  }

  for (Index i = 0; i < m_; ++i) {
    if (i == p) continue;
    const double v = spike[i];
    if (std::fabs(v) <= options_.drop_tol) continue;
    const Index k = pos_of_row_[i];
    ucol_[r].push_back({k, v});
    urow_[k].push_back({r, v});
    ++u_nnz_;
  }
  diag_[r] = diag;

  order_[slot_r] = -1;
  order_slot_[r] = static_cast<Index>(order_.size());
  order_.push_back(r);
  ++num_updates_;
  // Every later iteration's solves (BTRAN, entering FTRAN, DSE FTRAN) also pass the entries
  // that the updates so far added to the factor.
  constexpr double kSolvesPerIteration = 3.0;
  const NnzIndex added = factor_nnz() - initial_nnz_;
  if (added > 0) extra_solve_work_ += kSolvesPerIteration * static_cast<double>(added);
  valid_ = true;
  return true;
}

}  // namespace samaya
