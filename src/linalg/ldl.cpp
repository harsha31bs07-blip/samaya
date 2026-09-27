#include "linalg/ldl.hpp"

#include <algorithm>
#include <cassert>
#include <utility>

#include "linalg/ordering.hpp"

namespace samaya {

void LdlFactor::analyze(const SparseMatrix& upper, std::vector<Index> perm) {
  n_ = upper.cols();
  const auto n = static_cast<std::size_t>(n_);
  perm_ = perm.empty() ? minimum_degree_order(upper) : std::move(perm);
  inverse_ = inverse_permutation(perm_);

  // Permuted upper triangle: entry (i, j) of the input goes to (min, max) of the new positions.
  // A missing diagonal gets a structural entry (source -1, value 0).
  struct Entry {
    Index row;
    NnzIndex source;
  };
  std::vector<std::vector<Entry>> cols(n);
  std::vector<char> has_diagonal(n, 0);
  const auto start = upper.col_start();
  const auto index = upper.row_index();
  for (Index j = 0; j < n_; ++j) {
    for (NnzIndex p = start[j]; p < start[j + 1]; ++p) {
      const Index i = index[p];
      if (i > j) continue;  // Only the upper triangle is read.
      const Index pi = inverse_[static_cast<std::size_t>(i)];
      const Index pj = inverse_[static_cast<std::size_t>(j)];
      const Index r = std::min(pi, pj);
      const Index c = std::max(pi, pj);
      cols[static_cast<std::size_t>(c)].push_back({r, p});
      if (r == c) has_diagonal[static_cast<std::size_t>(c)] = 1;
    }
  }
  pu_start_.assign(n + 1, 0);
  pu_index_.clear();
  pu_source_.clear();
  for (std::size_t c = 0; c < n; ++c) {
    if (!has_diagonal[c]) cols[c].push_back({static_cast<Index>(c), -1});
    std::sort(cols[c].begin(), cols[c].end(),
              [](const Entry& a, const Entry& b) { return a.row < b.row; });
    for (const Entry& e : cols[c]) {
      pu_index_.push_back(e.row);
      pu_source_.push_back(e.source);
    }
    pu_start_[c + 1] = static_cast<NnzIndex>(pu_index_.size());
  }

  // Elimination tree and column counts of L (QDLDL's etree).
  etree_.assign(n, -1);
  std::vector<NnzIndex> count(n, 0);
  std::vector<Index> flag(n, -1);
  for (Index j = 0; j < n_; ++j) {
    flag[static_cast<std::size_t>(j)] = j;
    for (NnzIndex p = pu_start_[static_cast<std::size_t>(j)];
         p < pu_start_[static_cast<std::size_t>(j) + 1]; ++p) {
      Index i = pu_index_[static_cast<std::size_t>(p)];
      while (i != j && flag[static_cast<std::size_t>(i)] != j) {
        if (etree_[static_cast<std::size_t>(i)] == -1) etree_[static_cast<std::size_t>(i)] = j;
        ++count[static_cast<std::size_t>(i)];
        flag[static_cast<std::size_t>(i)] = j;
        i = etree_[static_cast<std::size_t>(i)];
      }
    }
  }
  l_start_.assign(n + 1, 0);
  for (std::size_t i = 0; i < n; ++i) l_start_[i + 1] = l_start_[i] + count[i];
  const auto nnz = static_cast<std::size_t>(l_start_[n]);
  l_index_.assign(nnz, 0);
  l_value_.assign(nnz, 0.0);
  d_.assign(n, 0.0);
  y_.assign(n, 0.0);
  y_mark_.assign(n, 0);
  y_index_.assign(n, 0);
  stack_.assign(n, 0);
  next_.assign(n, 0);
  work_.assign(n, 0.0);
}

Index LdlFactor::factorize(const SparseMatrix& upper, std::span<const int> sign, double pivot_tol,
                           double replacement) {
  const auto n = static_cast<std::size_t>(n_);
  const auto values = upper.values();
  Index replaced = 0;
  for (std::size_t i = 0; i < n; ++i) next_[i] = l_start_[i];
  for (Index k = 0; k < n_; ++k) {
    const auto uk = static_cast<std::size_t>(k);
    // Scatter column k of the permuted upper triangle into y, and collect the pattern of row k of
    // L: the nodes reached from its entries by walking up the elimination tree, in topological
    // order (each walk is pushed on the stack and then reversed onto y_index_).
    Index nnz_y = 0;
    double dk = 0.0;
    for (NnzIndex p = pu_start_[uk]; p < pu_start_[uk + 1]; ++p) {
      const Index i = pu_index_[static_cast<std::size_t>(p)];
      const NnzIndex source = pu_source_[static_cast<std::size_t>(p)];
      const double v = source >= 0 ? values[static_cast<std::size_t>(source)] : 0.0;
      if (i == k) {
        dk += v;
        continue;
      }
      y_[static_cast<std::size_t>(i)] += v;
      Index node = i;
      Index top = 0;
      while (node != -1 && node < k && !y_mark_[static_cast<std::size_t>(node)]) {
        y_mark_[static_cast<std::size_t>(node)] = 1;
        stack_[static_cast<std::size_t>(top++)] = node;
        node = etree_[static_cast<std::size_t>(node)];
      }
      while (top > 0) y_index_[static_cast<std::size_t>(nnz_y++)] = stack_[--top];
    }
    // Eliminate in reverse order of discovery, which is a topological order of the tree.
    for (Index t = nnz_y; t-- > 0;) {
      const Index c = y_index_[static_cast<std::size_t>(t)];
      const auto uc = static_cast<std::size_t>(c);
      const double yc = y_[uc];
      const NnzIndex end = next_[uc];
      for (NnzIndex p = l_start_[uc]; p < end; ++p) {
        y_[static_cast<std::size_t>(l_index_[static_cast<std::size_t>(p)])] -=
            l_value_[static_cast<std::size_t>(p)] * yc;
      }
      const double lkc = yc / d_[uc];
      l_index_[static_cast<std::size_t>(end)] = k;
      l_value_[static_cast<std::size_t>(end)] = lkc;
      next_[uc] = end + 1;
      dk -= yc * lkc;
      y_[uc] = 0.0;
      y_mark_[uc] = 0;
    }
    const double s = sign[static_cast<std::size_t>(perm_[uk])] > 0 ? 1.0 : -1.0;
    if (!(s * dk >= pivot_tol)) {
      dk = s * replacement;
      ++replaced;
    }
    d_[uk] = dk;
  }
  return replaced;
}

void LdlFactor::solve(std::vector<double>& x) const {
  const auto n = static_cast<std::size_t>(n_);
  assert(x.size() == n);
  std::vector<double>& w = work_;
  for (std::size_t k = 0; k < n; ++k) w[k] = x[static_cast<std::size_t>(perm_[k])];
  for (std::size_t i = 0; i < n; ++i) {
    const double wi = w[i];
    if (wi == 0.0) continue;
    for (NnzIndex p = l_start_[i]; p < l_start_[i + 1]; ++p) {
      w[static_cast<std::size_t>(l_index_[static_cast<std::size_t>(p)])] -=
          l_value_[static_cast<std::size_t>(p)] * wi;
    }
  }
  for (std::size_t i = 0; i < n; ++i) w[i] /= d_[i];
  for (std::size_t i = n; i-- > 0;) {
    double sum = w[i];
    for (NnzIndex p = l_start_[i]; p < l_start_[i + 1]; ++p) {
      const auto up = static_cast<std::size_t>(p);
      sum -= l_value_[up] * w[static_cast<std::size_t>(l_index_[up])];
    }
    w[i] = sum;
  }
  for (std::size_t k = 0; k < n; ++k) x[static_cast<std::size_t>(perm_[k])] = w[k];
}

}  // namespace samaya
