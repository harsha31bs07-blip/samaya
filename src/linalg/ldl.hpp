#pragma once

#include <span>
#include <vector>

#include "samaya/sparse_matrix.hpp"

namespace samaya {

// Sparse LDL^T factorization of a symmetric quasi-definite matrix P K P^T = L D L^T (L unit lower
// triangular, D diagonal), for the interior-point KKT systems.
//
// A quasi-definite matrix [-E  A^T; A  F] with E, F positive definite has such a factorization
// for every symmetric permutation (Vanderbei, SIAM J. Optim. 1995), so the fill-reducing ordering
// is chosen once from the pattern and no pivoting is needed. Each pivot must have the sign given
// for its row; a pivot of the wrong sign or smaller than the regularization is replaced by
// sign * regularization (static regularization, as in QDLDL-based and Clarabel-type solvers) and
// counted, and the solve then refines against the unregularized matrix.
//
// The numeric phase is the up-looking algorithm of QDLDL (Stellato et al., OSQP, 2020): row k of
// L is found by walking the elimination tree from the entries of column k of the upper triangle.
class LdlFactor {
 public:
  // Orders and analyzes the pattern of `upper`, the upper triangle (row <= col) of an n x n
  // symmetric matrix in CSC form (entries below the diagonal are ignored). The pattern must stay
  // the same for factorize(); values may change. `perm` fixes the ordering (perm[k] = original
  // index eliminated k-th); empty means minimum degree.
  void analyze(const SparseMatrix& upper, std::vector<Index> perm = {});

  // Factorizes the values of `upper` (same pattern as analyze()). sign[i] is the required sign of
  // the pivot of original row i (+1 or -1). A pivot d with sign * d < pivot_tol is replaced by
  // sign * replacement. Returns the number of pivots replaced.
  Index factorize(const SparseMatrix& upper, std::span<const int> sign, double pivot_tol,
                  double replacement);

  // Solves (P^T L D L^T P) x = b in place: the factorized, possibly regularized matrix.
  void solve(std::vector<double>& x) const;

  Index size() const { return n_; }
  NnzIndex factor_nnz() const { return static_cast<NnzIndex>(l_index_.size()); }
  const std::vector<Index>& permutation() const { return perm_; }

 private:
  Index n_ = 0;
  std::vector<Index> perm_;      // perm_[k] = original index at position k.
  std::vector<Index> inverse_;   // inverse_[i] = position of original index i.
  // The permuted upper triangle: its pattern, and for each of its entries the position of the
  // entry of the input matrix it takes its value from.
  std::vector<NnzIndex> pu_start_;
  std::vector<Index> pu_index_;
  std::vector<NnzIndex> pu_source_;
  std::vector<Index> etree_;     // Parent in the elimination tree, -1 for roots.
  std::vector<NnzIndex> l_start_;
  std::vector<Index> l_index_;
  std::vector<double> l_value_;
  std::vector<double> d_;
  // Scratch for factorize() and solve().
  std::vector<double> y_;
  std::vector<char> y_mark_;
  std::vector<Index> y_index_;
  std::vector<Index> stack_;
  std::vector<NnzIndex> next_;
  mutable std::vector<double> work_;
};

}  // namespace samaya
