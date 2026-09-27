#include "linalg/ordering.hpp"

#include <algorithm>
#include <cmath>
#include <set>
#include <utility>

namespace samaya {

namespace {

// Nodes denser than this are ordered last, outside the graph (as AMD does): a dense row of a KKT
// system would otherwise make every degree update touch it. Threshold as in AMD's default,
// max(16, 10 sqrt(n)).
Index dense_threshold(Index n) {
  return std::max<Index>(16, static_cast<Index>(10.0 * std::sqrt(static_cast<double>(n))));
}

}  // namespace

std::vector<Index> minimum_degree_order(const SparseMatrix& pattern) {
  const Index n = pattern.cols();
  const auto nn = static_cast<std::size_t>(n);
  std::vector<std::vector<Index>> vars(nn);  // Variable neighbours (uneliminated).
  {
    const auto start = pattern.col_start();
    const auto index = pattern.row_index();
    for (Index j = 0; j < n; ++j) {
      for (NnzIndex p = start[j]; p < start[j + 1]; ++p) {
        const Index i = index[p];
        if (i == j) continue;
        vars[static_cast<std::size_t>(i)].push_back(j);
        vars[static_cast<std::size_t>(j)].push_back(i);
      }
    }
    for (auto& v : vars) {
      std::sort(v.begin(), v.end());
      v.erase(std::unique(v.begin(), v.end()), v.end());
    }
  }

  std::vector<char> eliminated(nn, 0);  // Eliminated, or set aside as dense.
  std::vector<Index> dense;
  const Index threshold = dense_threshold(n);
  for (Index i = 0; i < n; ++i) {
    if (static_cast<Index>(vars[static_cast<std::size_t>(i)].size()) > threshold) {
      eliminated[static_cast<std::size_t>(i)] = 1;
      dense.push_back(i);
    }
  }

  std::vector<std::vector<Index>> elems(nn);  // Elements adjacent to each variable.
  std::vector<std::vector<Index>> evars(nn);  // Variables of each element (an eliminated pivot).
  std::vector<char> absorbed(nn, 0);
  std::vector<Index> mark(nn, -1);
  std::vector<Index> mark2(nn, -1);
  Index stamp = 0;
  Index stamp2 = 0;
  std::vector<Index> degree(nn, 0);
  std::set<std::pair<Index, Index>> queue;  // (degree, variable)

  const auto compute_degree = [&](Index i) {
    ++stamp2;
    mark2[static_cast<std::size_t>(i)] = stamp2;
    Index count = 0;
    for (const Index v : vars[static_cast<std::size_t>(i)]) {
      if (mark2[static_cast<std::size_t>(v)] != stamp2) {
        mark2[static_cast<std::size_t>(v)] = stamp2;
        ++count;
      }
    }
    for (const Index e : elems[static_cast<std::size_t>(i)]) {
      auto& list = evars[static_cast<std::size_t>(e)];
      std::size_t keep = 0;
      for (const Index v : list) {
        if (eliminated[static_cast<std::size_t>(v)]) continue;
        list[keep++] = v;
        if (mark2[static_cast<std::size_t>(v)] != stamp2) {
          mark2[static_cast<std::size_t>(v)] = stamp2;
          ++count;
        }
      }
      list.resize(keep);
    }
    return count;
  };

  for (Index i = 0; i < n; ++i) {
    if (eliminated[static_cast<std::size_t>(i)]) continue;
    auto& v = vars[static_cast<std::size_t>(i)];
    v.erase(std::remove_if(v.begin(), v.end(),
                           [&](Index u) { return eliminated[static_cast<std::size_t>(u)] != 0; }),
            v.end());
    degree[static_cast<std::size_t>(i)] = static_cast<Index>(v.size());
    queue.insert({degree[static_cast<std::size_t>(i)], i});
  }

  std::vector<Index> perm;
  perm.reserve(nn);
  std::vector<Index> lp;
  while (!queue.empty()) {
    const Index p = queue.begin()->second;
    queue.erase(queue.begin());
    const auto up = static_cast<std::size_t>(p);
    eliminated[up] = 1;
    perm.push_back(p);

    // The new element: p's variables and the variables of its elements, which it absorbs.
    ++stamp;
    mark[up] = stamp;
    lp.clear();
    for (const Index v : vars[up]) {
      const auto uv = static_cast<std::size_t>(v);
      if (!eliminated[uv] && mark[uv] != stamp) {
        mark[uv] = stamp;
        lp.push_back(v);
      }
    }
    for (const Index e : elems[up]) {
      const auto ue = static_cast<std::size_t>(e);
      if (absorbed[ue]) continue;
      for (const Index v : evars[ue]) {
        const auto uv = static_cast<std::size_t>(v);
        if (!eliminated[uv] && mark[uv] != stamp) {
          mark[uv] = stamp;
          lp.push_back(v);
        }
      }
      absorbed[ue] = 1;
      std::vector<Index>().swap(evars[ue]);
    }
    std::vector<Index>().swap(vars[up]);
    std::vector<Index>().swap(elems[up]);
    evars[up] = lp;

    for (const Index i : lp) {
      const auto ui = static_cast<std::size_t>(i);
      auto& el = elems[ui];
      el.erase(std::remove_if(el.begin(), el.end(),
                              [&](Index e) { return absorbed[static_cast<std::size_t>(e)] != 0; }),
               el.end());
      el.push_back(p);
      // Variables also in the new element are reachable through it.
      auto& vl = vars[ui];
      vl.erase(std::remove_if(vl.begin(), vl.end(),
                              [&](Index v) {
                                const auto uv = static_cast<std::size_t>(v);
                                return eliminated[uv] != 0 || mark[uv] == stamp;
                              }),
               vl.end());
    }
    for (const Index i : lp) {
      const auto ui = static_cast<std::size_t>(i);
      queue.erase({degree[ui], i});
      degree[ui] = compute_degree(i);
      queue.insert({degree[ui], i});
    }
  }
  perm.insert(perm.end(), dense.begin(), dense.end());
  return perm;
}

std::vector<Index> inverse_permutation(const std::vector<Index>& perm) {
  std::vector<Index> inverse(perm.size());
  for (std::size_t k = 0; k < perm.size(); ++k) {
    inverse[static_cast<std::size_t>(perm[k])] = static_cast<Index>(k);
  }
  return inverse;
}

}  // namespace samaya
