#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include <vector>

#include "linalg/ldl.hpp"
#include "linalg/ordering.hpp"
#include "test_framework.hpp"

using samaya::Index;
using samaya::LdlFactor;
using samaya::SparseMatrix;
using samaya::Triplet;

namespace {

// A random quasi-definite matrix K = [-E  A^T; A  F]: E and F sparse, symmetric and strictly
// diagonally dominant (so positive definite), A sparse. Returned as its full symmetric triplets.
std::vector<Triplet> random_quasi_definite(int n1, int n2, double density, std::mt19937& rng) {
  std::uniform_real_distribution<double> value(-1.0, 1.0);
  std::bernoulli_distribution keep(density);
  const int n = n1 + n2;
  std::vector<double> row_sum(static_cast<std::size_t>(n), 0.0);
  std::vector<Triplet> t;
  const auto add_sym = [&](int i, int j, double v) {
    t.push_back({i, j, v});
    t.push_back({j, i, v});
  };
  // Off-diagonal entries of E (block 1) and F (block 2), and of A (coupling).
  for (int j = 0; j < n; ++j) {
    for (int i = j + 1; i < n; ++i) {
      if (!keep(rng)) continue;
      const double v = value(rng);
      const bool same_block = (i < n1) == (j < n1);
      add_sym(i, j, same_block && j < n1 ? -v : v);
      if (same_block) {
        row_sum[static_cast<std::size_t>(i)] += std::fabs(v);
        row_sum[static_cast<std::size_t>(j)] += std::fabs(v);
      }
    }
  }
  for (int i = 0; i < n; ++i) {
    const double d = row_sum[static_cast<std::size_t>(i)] + 0.1 + std::fabs(value(rng));
    t.push_back({i, i, i < n1 ? -d : d});
  }
  return t;
}

SparseMatrix upper_triangle(int n, const std::vector<Triplet>& full) {
  std::vector<Triplet> t;
  for (const Triplet& e : full) {
    if (e.row <= e.col) t.push_back(e);
  }
  return SparseMatrix::from_triplets(n, n, std::move(t));
}

double residual(int n, const std::vector<Triplet>& full, const std::vector<double>& x,
                const std::vector<double>& b) {
  std::vector<double> r = b;
  for (const Triplet& e : full) {
    r[static_cast<std::size_t>(e.row)] -= e.value * x[static_cast<std::size_t>(e.col)];
  }
  double rmax = 0.0;
  double bmax = 0.0;
  for (int i = 0; i < n; ++i) {
    rmax = std::max(rmax, std::fabs(r[static_cast<std::size_t>(i)]));
    bmax = std::max(bmax, std::fabs(b[static_cast<std::size_t>(i)]));
  }
  return rmax / (1.0 + bmax);
}

// Nonzeros of L for the pattern of `upper` under `perm` (natural order if empty).
samaya::NnzIndex fill(const SparseMatrix& upper, std::vector<Index> perm) {
  LdlFactor f;
  f.analyze(upper, std::move(perm));
  return f.factor_nnz();
}

}  // namespace

TEST(ldl_solves_random_quasi_definite_systems) {
  // Random quasi-definite matrices of 2..80 rows with the pivot signs of the blocks: every solve
  // must reproduce b to 1e-10, and no pivot may need regularization.
  std::mt19937 rng(11);
  std::uniform_real_distribution<double> value(-1.0, 1.0);
  int systems = 0;
  double worst = 0.0;
  Index replaced = 0;
  for (int k = 0; k < 300; ++k) {
    const int n1 = 1 + static_cast<int>(rng() % 50);
    const int n2 = 1 + static_cast<int>(rng() % 30);
    const int n = n1 + n2;
    const double density = std::uniform_real_distribution<double>(0.02, 0.3)(rng);
    const std::vector<Triplet> full = random_quasi_definite(n1, n2, density, rng);
    const SparseMatrix upper = upper_triangle(n, full);
    std::vector<int> sign(static_cast<std::size_t>(n), 1);
    std::fill(sign.begin(), sign.begin() + n1, -1);
    LdlFactor f;
    f.analyze(upper);
    replaced += f.factorize(upper, sign, 1e-13, 1e-8);
    for (int rhs = 0; rhs < 2; ++rhs) {
      std::vector<double> b(static_cast<std::size_t>(n));
      for (double& v : b) v = value(rng);
      std::vector<double> x = b;
      f.solve(x);
      worst = std::max(worst, residual(n, full, x, b));
    }
    ++systems;
  }
  std::printf("  %d systems, worst relative residual %.2e, %d pivots regularized\n", systems,
              worst, replaced);
  CHECK(worst <= 1e-10);
  CHECK_EQ(replaced, 0);
}

TEST(ldl_values_can_change_between_factorizations) {
  // The IPM refactorizes the same pattern with new values each iteration.
  std::mt19937 rng(12);
  const int n1 = 30;
  const int n2 = 20;
  const int n = n1 + n2;
  std::vector<Triplet> full = random_quasi_definite(n1, n2, 0.1, rng);
  LdlFactor f;
  f.analyze(upper_triangle(n, full));
  std::vector<int> sign(static_cast<std::size_t>(n), 1);
  std::fill(sign.begin(), sign.begin() + n1, -1);
  for (int round = 0; round < 5; ++round) {
    for (Triplet& e : full) {
      if (e.row == e.col) e.value *= 1.0 + 0.5 * round;  // Diagonal changes, as with Theta.
    }
    const SparseMatrix upper = upper_triangle(n, full);
    CHECK_EQ(f.factorize(upper, sign, 1e-13, 1e-8), 0);
    std::vector<double> b(static_cast<std::size_t>(n), 1.0);
    std::vector<double> x = b;
    f.solve(x);
    CHECK(residual(n, full, x, b) <= 1e-10);
  }
}

TEST(ldl_regularizes_wrong_sign_pivots) {
  // A zero (2,2) block: the pivots of block 2 are zero before elimination; with a nonsingular
  // coupling the elimination makes them positive, but a zero coupling row leaves a zero pivot,
  // which must be replaced (not divided by) and counted.
  const int n = 3;
  const std::vector<Triplet> full = {{0, 0, -2.0}, {1, 0, 1.0}, {0, 1, 1.0}, {1, 1, 0.0},
                                     {2, 2, 0.0}};
  std::vector<Triplet> t;
  for (const Triplet& e : full) {
    if (e.row <= e.col && e.value != 0.0) t.push_back(e);
  }
  t.push_back({2, 2, 0.0});
  const SparseMatrix upper = SparseMatrix::from_triplets(n, n, t);
  LdlFactor f;
  f.analyze(upper, {0, 1, 2});  // Natural order: pivots -2, then 0.5, then 0.
  const std::vector<int> sign = {-1, 1, 1};
  CHECK_EQ(f.factorize(upper, sign, 1e-13, 1e-8), 1);  // Row 2 is empty.
  std::vector<double> x = {1.0, 1.0, 0.0};
  f.solve(x);
  for (const double v : x) CHECK(std::isfinite(v));
  // Rows 0-1 are exact: -2 x0 + x1 = 1, x0 = 1 (x1 has a 0.5 pivot after elimination).
  CHECK_NEAR(x[0], 1.0, 1e-12);
  CHECK_NEAR(x[1], 3.0, 1e-12);
}

TEST(minimum_degree_reduces_fill) {
  // A 2-D grid Laplacian (the classic fill test) and an arrow matrix whose dense row comes first:
  // minimum degree must beat the natural order on both, and put the dense row last (fill n - 1).
  const int g = 20;
  const int n = g * g;
  std::vector<Triplet> t;
  for (int r = 0; r < g; ++r) {
    for (int c = 0; c < g; ++c) {
      const int i = r * g + c;
      t.push_back({i, i, 4.0});
      if (c + 1 < g) t.push_back({i, i + 1, -1.0});
      if (r + 1 < g) t.push_back({i, i + g, -1.0});
    }
  }
  const SparseMatrix grid = SparseMatrix::from_triplets(n, n, t);
  std::vector<Index> natural(static_cast<std::size_t>(n));
  std::iota(natural.begin(), natural.end(), Index{0});
  const samaya::NnzIndex md = fill(grid, {});
  const samaya::NnzIndex nat = fill(grid, natural);
  std::printf("  grid %dx%d: fill %lld (minimum degree) vs %lld (natural)\n", g, g,
              static_cast<long long>(md), static_cast<long long>(nat));
  CHECK(md < nat);

  const int m = 200;
  std::vector<Triplet> a;
  for (int i = 0; i < m; ++i) {
    a.push_back({i, i, 1.0});
    if (i > 0) a.push_back({0, i, 1.0});
  }
  const SparseMatrix arrow = SparseMatrix::from_triplets(m, m, a);
  std::vector<Index> nat_arrow(static_cast<std::size_t>(m));
  std::iota(nat_arrow.begin(), nat_arrow.end(), Index{0});
  CHECK_EQ(fill(arrow, {}), m - 1);
  CHECK_EQ(fill(arrow, nat_arrow), static_cast<samaya::NnzIndex>(m) * (m - 1) / 2);

  const std::vector<Index> perm = samaya::minimum_degree_order(grid);
  std::vector<Index> sorted = perm;
  std::sort(sorted.begin(), sorted.end());
  CHECK(sorted == natural);  // A permutation.
}
