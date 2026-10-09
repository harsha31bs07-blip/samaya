#pragma once

#include <cstddef>
#include <vector>

#include "samaya/sparse_matrix.hpp"

namespace samaya {

// Fill-reducing symmetric ordering: minimum degree on the quotient graph (George and Liu, "The
// evolution of the minimum degree ordering algorithm", SIAM Review 1989), with exact external
// degrees and element absorption. `pattern` is a square matrix whose nonzero pattern, taken
// symmetrically (entries (i, j) and (j, i) alike, the diagonal ignored), is the graph to order.
// Returns perm with perm[k] = the original index eliminated k-th.
//
// A degree is updated exactly when that touches at most exact_degree_work entries, and with AMD's
// approximate bound (Amestoy, Davis and Duff 1996) above it: small models keep the exact order,
// large ones stay fast. 0 approximates every update; SIZE_MAX makes every update exact.
inline constexpr std::size_t kExactDegreeWork = 4096;
std::vector<Index> minimum_degree_order(const SparseMatrix& pattern,
                                        std::size_t exact_degree_work = kExactDegreeWork);

// Inverse of a permutation: inverse[perm[k]] = k.
std::vector<Index> inverse_permutation(const std::vector<Index>& perm);

}  // namespace samaya
