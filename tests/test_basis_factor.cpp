#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

#include "dense_reference.hpp"
#include "linalg/basis_factor.hpp"
#include "test_framework.hpp"

using samaya::BasisFactor;
using samaya::Index;
using samaya::SparseMatrix;
using samaya::Triplet;
using samaya::test::DenseMatrix;

namespace {

SparseMatrix random_matrix(int m, int n, double density, std::mt19937& rng) {
  std::uniform_real_distribution<double> value(-1.0, 1.0);
  // Callers pass 3.0 / m, which exceeds 1 for tiny m; a probability above 1 is undefined (MSVC's
  // debug library asserts). Capped, it means "keep every entry", as it did with libstdc++.
  std::bernoulli_distribution keep(std::min(1.0, density));
  std::vector<Triplet> t;
  for (int j = 0; j < n; ++j) {
    t.push_back({static_cast<Index>(rng() % m), j, value(rng) + 2.0});
    for (int i = 0; i < m; ++i) {
      if (keep(rng)) t.push_back({i, j, value(rng)});
    }
  }
  return SparseMatrix::from_triplets(m, n, std::move(t));
}

// Dense B for the given basis of [A -I].
DenseMatrix dense_basis(const SparseMatrix& a, const std::vector<Index>& basic) {
  const Index m = a.rows();
  DenseMatrix b(m, std::vector<double>(m, 0.0));
  for (Index k = 0; k < m; ++k) {
    const Index j = basic[k];
    if (j < a.cols()) {
      for (auto p = a.col_start()[j]; p < a.col_start()[j + 1]; ++p) {
        b[a.row_index()[p]][k] = a.values()[p];
      }
    } else {
      b[j - a.cols()][k] = -1.0;
    }
  }
  return b;
}

std::vector<double> column(const SparseMatrix& a, Index j) {
  std::vector<double> c(a.rows(), 0.0);
  if (j < a.cols()) {
    for (auto p = a.col_start()[j]; p < a.col_start()[j + 1]; ++p) c[a.row_index()[p]] = a.values()[p];
  } else {
    c[j - a.cols()] = -1.0;
  }
  return c;
}

double max_abs(const std::vector<double>& v) {
  double m = 0.0;
  for (double x : v) m = std::max(m, std::fabs(x));
  return m;
}

// Checks ftran and btran against the dense basis on random right-hand sides.
void check_solves(const BasisFactor& f, const DenseMatrix& b, std::mt19937& rng) {
  std::uniform_real_distribution<double> value(-1.0, 1.0);
  const std::size_t m = b.size();
  for (int trial = 0; trial < 3; ++trial) {
    std::vector<double> rhs(m);
    for (double& v : rhs) v = (rng() % 3 == 0) ? 0.0 : value(rng);
    std::vector<double> z = rhs;
    f.ftran(z);
    const std::vector<double> bz = samaya::test::multiply(b, z);
    double err = 0.0;
    for (std::size_t i = 0; i < m; ++i) err = std::max(err, std::fabs(bz[i] - rhs[i]));
    CHECK(err <= 1e-9 * (1.0 + max_abs(z)));

    std::vector<double> w = rhs;
    f.btran(w);
    const std::vector<double> btw = samaya::test::multiply(samaya::test::transpose(b), w);
    err = 0.0;
    for (std::size_t i = 0; i < m; ++i) err = std::max(err, std::fabs(btw[i] - rhs[i]));
    CHECK(err <= 1e-9 * (1.0 + max_abs(w)));
  }
}

// Random basis mixing structural and logical columns; repaired until nonsingular.
std::vector<Index> random_basis(const SparseMatrix& a, BasisFactor& f, std::mt19937& rng) {
  const Index m = a.rows();
  const Index n = a.cols();
  std::vector<Index> all(static_cast<std::size_t>(n + m));
  for (Index j = 0; j < n + m; ++j) all[j] = j;
  std::shuffle(all.begin(), all.end(), rng);
  std::vector<Index> basic(all.begin(), all.begin() + m);
  for (int attempt = 0; attempt < 5; ++attempt) {
    if (f.factorize(basic) == 0) return basic;
    for (std::size_t k = 0; k < f.singular_positions().size(); ++k) {
      basic[f.singular_positions()[k]] = n + f.unpivoted_rows()[k];
    }
  }
  CHECK(false);
  return basic;
}

}  // namespace

TEST(factor_slack_basis) {
  std::mt19937 rng(1);
  const SparseMatrix a = random_matrix(6, 4, 0.3, rng);
  BasisFactor f(a);
  std::vector<Index> basic;
  for (Index i = 0; i < 6; ++i) basic.push_back(4 + i);
  REQUIRE(f.factorize(basic) == 0);
  check_solves(f, dense_basis(a, basic), rng);
}

TEST(factor_matches_dense_on_random_bases) {
  std::mt19937 rng(42);
  for (int trial = 0; trial < 60; ++trial) {
    const int m = 1 + static_cast<int>(rng() % 40);
    const int n = m + static_cast<int>(rng() % 40);
    const double density = (trial % 3 == 0) ? 0.5 : 3.0 / m;
    const SparseMatrix a = random_matrix(m, n, density, rng);
    BasisFactor f(a);
    const std::vector<Index> basic = random_basis(a, f, rng);
    check_solves(f, dense_basis(a, basic), rng);
  }
}

TEST(factor_forrest_tomlin_updates_match_dense) {
  std::mt19937 rng(7);
  int updates_done = 0;
  int refactors = 0;
  for (int trial = 0; trial < 20; ++trial) {
    const int m = 5 + static_cast<int>(rng() % 60);
    const int n = 2 * m;
    const SparseMatrix a = random_matrix(m, n, 4.0 / m, rng);
    BasisFactor::Options options;
    options.max_updates = 25;
    BasisFactor f(a, options);
    std::vector<Index> basic = random_basis(a, f, rng);
    for (int it = 0; it < 80; ++it) {
      // Pick an entering column not in the basis with a safe pivot element.
      const Index j = static_cast<Index>(rng() % (n + m));
      if (std::find(basic.begin(), basic.end(), j) != basic.end()) continue;
      std::vector<double> alpha = column(a, j);
      std::vector<double> spike;
      f.ftran(alpha, &spike);
      const Index r = static_cast<Index>(rng() % m);
      if (std::fabs(alpha[r]) < 1e-2 * max_abs(alpha)) continue;
      basic[r] = j;
      if (f.should_refactor() || !f.update(r, spike)) {
        REQUIRE(f.factorize(basic) == 0);
        ++refactors;
      } else {
        ++updates_done;
      }
      check_solves(f, dense_basis(a, basic), rng);
    }
  }
  CHECK(updates_done > 500);
  CHECK(refactors > 0);
}

TEST(factor_detects_rank_deficiency) {
  // Columns 0 and 1 are identical; column 2 is their negation plus nothing new.
  const SparseMatrix a = SparseMatrix::from_triplets(
      3, 3, {{0, 0, 1.0}, {1, 0, 2.0}, {0, 1, 1.0}, {1, 1, 2.0}, {2, 2, 5.0}});
  BasisFactor f(a);
  std::vector<Index> basic{0, 1, 2};
  CHECK_EQ(f.factorize(basic), 1);
  REQUIRE(f.singular_positions().size() == 1);
  const Index pos = f.singular_positions()[0];
  CHECK(pos == 0 || pos == 1);
  basic[pos] = 3 + f.unpivoted_rows()[0];
  REQUIRE(f.factorize(basic) == 0);
  std::mt19937 rng(3);
  check_solves(f, dense_basis(a, basic), rng);
}

TEST(factor_detects_numerically_tiny_column) {
  const SparseMatrix a = SparseMatrix::from_triplets(2, 2, {{0, 0, 1.0}, {1, 1, 1e-13}});
  BasisFactor f(a);
  CHECK_EQ(f.factorize(std::vector<Index>{0, 1}), 1);
  CHECK_EQ(f.singular_positions()[0], 1);
  CHECK_EQ(f.unpivoted_rows()[0], 1);
}

TEST(factor_larger_sparse_basis_limits_fill) {
  std::mt19937 rng(11);
  const int m = 400;
  const SparseMatrix a = random_matrix(m, 2 * m, 3.0 / m, rng);
  BasisFactor f(a);
  std::vector<Index> basic(m);
  for (int k = 0; k < m; ++k) basic[k] = k;  // All structural: a genuine nucleus.
  for (int attempt = 0; attempt < 5 && f.factorize(basic) != 0; ++attempt) {
    for (std::size_t k = 0; k < f.singular_positions().size(); ++k) {
      basic[f.singular_positions()[k]] = 2 * m + f.unpivoted_rows()[k];
    }
  }
  REQUIRE(f.factorize(basic) == 0);
  check_solves(f, dense_basis(a, basic), rng);
  // Markowitz ordering keeps the factor far sparser than the dense m*m/2.
  CHECK(f.factor_nnz() < 20 * a.nnz());
}
