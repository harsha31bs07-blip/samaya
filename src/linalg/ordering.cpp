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

std::vector<Index> minimum_degree_order(const SparseMatrix& pattern,
                                        std::size_t exact_degree_work) {
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
  // w[e] = |Le \ Lp| for the elements next to the new element's variables (AMD's approximate
  // degree), valid where wmark[e] == the current pivot's stamp.
  std::vector<Index> w(nn, 0);
  std::vector<Index> wmark(nn, -1);
  std::vector<Index> mark2(nn, -1);  // For the exact degree.
  Index stamp2 = 0;
  Index stamp = 0;
  Index remaining = 0;  // Variables not yet eliminated (dense ones excluded).
  std::vector<Index> degree(nn, 0);
  std::set<std::pair<Index, Index>> queue;  // (degree, variable)

  for (Index i = 0; i < n; ++i) {
    if (eliminated[static_cast<std::size_t>(i)]) continue;
    auto& v = vars[static_cast<std::size_t>(i)];
    v.erase(std::remove_if(v.begin(), v.end(),
                           [&](Index u) { return eliminated[static_cast<std::size_t>(u)] != 0; }),
            v.end());
    degree[static_cast<std::size_t>(i)] = static_cast<Index>(v.size());
    queue.insert({degree[static_cast<std::size_t>(i)], i});
    ++remaining;
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
    --remaining;

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
    // Approximate degrees (Amestoy, Davis and Duff 1996): w[e] = |Le \ Lp| for every element e
    // next to Lp, from one pass over their element lists; the exact degree, a union over all of
    // a variable's elements, cost about 16 s on a 400 x 400 transportation KKT system.
    for (const Index i : lp) {
      for (const Index e : elems[static_cast<std::size_t>(i)]) {
        const auto ue = static_cast<std::size_t>(e);
        if (e == p) continue;
        if (wmark[ue] != stamp) {
          wmark[ue] = stamp;
          w[ue] = static_cast<Index>(evars[ue].size());
        }
        --w[ue];
      }
    }
    const Index lp_size = static_cast<Index>(lp.size());
    for (const Index i : lp) {
      const auto ui = static_cast<std::size_t>(i);
      Index external = 0;  // Sum of |Le \ Lp| over i's other elements.
      auto& el = elems[ui];
      std::size_t keep = 0;
      for (const Index e : el) {
        const auto ue = static_cast<std::size_t>(e);
        if (absorbed[ue]) continue;
        if (e != p && w[ue] == 0) {
          // Aggressive absorption: Le is inside Lp, so the new element covers it.
          absorbed[ue] = 1;
          std::vector<Index>().swap(evars[ue]);
          continue;
        }
        if (e != p) external += w[ue];
        el[keep++] = e;
      }
      el.resize(keep);
      const Index own = static_cast<Index>(vars[ui].size());
      std::size_t work = vars[ui].size();
      for (const Index e : el) work += evars[static_cast<std::size_t>(e)].size();
      Index d = 0;
      // Exact while cheap (kExactDegreeWork): small models keep the exact minimum-degree order,
      // which the barrier was developed on (with AMD's bound everywhere, beaconfd and forplan
      // stopped verifying from the barrier alone); a 400 x 400 transportation KKT system orders
      // in about 2 s instead of 16 s.
      if (work <= exact_degree_work) {
        // Exact: i's own neighbours and the variables of its elements (p's among them), counted
        // once each.
        ++stamp2;
        mark2[ui] = stamp2;
        for (const Index v : vars[ui]) {
          if (mark2[static_cast<std::size_t>(v)] != stamp2) {
            mark2[static_cast<std::size_t>(v)] = stamp2;
            ++d;
          }
        }
        for (const Index e : el) {
          for (const Index v : evars[static_cast<std::size_t>(e)]) {
            if (mark2[static_cast<std::size_t>(v)] != stamp2) {
              mark2[static_cast<std::size_t>(v)] = stamp2;
              ++d;
            }
          }
        }
      } else {
        // AMD's bound: no more than the variables left, than the old degree plus the new
        // element, or than i's own neighbours plus the new element plus the other elements
        // outside it.
        d = std::min({remaining - 1, degree[ui] + lp_size - 1, own + lp_size - 1 + external});
      }
      queue.erase({degree[ui], i});
      degree[ui] = std::max<Index>(d, 0);
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
