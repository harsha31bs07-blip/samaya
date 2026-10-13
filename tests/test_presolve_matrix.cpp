#include <algorithm>
#include <random>
#include <vector>

#include "presolve/presolve.hpp"
#include "test_framework.hpp"

namespace samaya {

// Exercise internal matrix primitives without exposing them in the Presolve API.
struct PresolveTestAccess {
  static bool add(Presolve& p, Index i, Index j, double a) { return p.add_entry(i, j, a); }
  static bool change(Presolve& p, Index i, Index j, double a) {
    return p.change_coefficient(i, j, a);
  }
  static bool remove(Presolve& p, Index i, Index j) { return p.remove_entry(i, j); }
  static void deactivate_row(Presolve& p, Index i) { p.remove_row(i); }
  static void deactivate_col(Presolve& p, Index j) { p.remove_col(j); }

  static void check(const Presolve& p, const std::vector<std::vector<double>>& dense) {
    std::vector<std::vector<double>> rows(p.m_, std::vector<double>(p.n_, 0.0));
    std::vector<std::vector<double>> cols = rows;
    std::vector<Index> row_counts(p.m_, 0), col_counts(p.n_, 0);
    for (Index i = 0; i < p.m_; ++i) {
      Index previous = -1;
      for (const auto& e : p.rows_[i]) {
        CHECK(e.index > previous);
        CHECK(e.a != 0.0);
        previous = e.index;
        rows[i][e.index] = e.a;
        if (p.row_active_[i] && p.col_active_[e.index]) ++row_counts[i];
      }
    }
    for (Index j = 0; j < p.n_; ++j) {
      Index previous = -1;
      for (const auto& e : p.cols_[j]) {
        CHECK(e.index > previous);
        CHECK(e.a != 0.0);
        previous = e.index;
        cols[e.index][j] = e.a;
        if (p.col_active_[j] && p.row_active_[e.index]) ++col_counts[j];
      }
    }
    CHECK(rows == dense);
    CHECK(cols == dense);
    CHECK(p.row_count_ == row_counts);
    CHECK(p.col_count_ == col_counts);
  }
};

}  // namespace samaya

TEST(presolve_mutable_matrix_random_operations) {
  constexpr samaya::Index kRows = 11;
  constexpr samaya::Index kCols = 13;
  samaya::Model model;
  model.obj.assign(kCols, 0.0);
  model.col_lower.assign(kCols, 0.0);
  model.col_upper.assign(kCols, samaya::kInf);
  model.col_type.assign(kCols, samaya::VarType::kContinuous);
  model.row_lower.assign(kRows, -samaya::kInf);
  model.row_upper.assign(kRows, samaya::kInf);
  model.A = samaya::SparseMatrix::from_triplets(kRows, kCols, {{3, 7, 2.0}});
  samaya::Presolve p(model);
  using Access = samaya::PresolveTestAccess;
  std::vector<std::vector<double>> dense(kRows, std::vector<double>(kCols, 0.0));
  dense[3][7] = 2.0;
  Access::check(p, dense);
  // Both the row and column lose their last entry.
  CHECK(Access::remove(p, 3, 7));
  dense[3][7] = 0.0;
  Access::check(p, dense);
  CHECK(!Access::remove(p, 3, 7));
  CHECK(!Access::change(p, 3, 7, 1.0));
  CHECK(!Access::add(p, 3, 7, 0.0));

  std::mt19937 rng(5140);
  for (int operation = 0; operation < 5000; ++operation) {
    const auto i = static_cast<samaya::Index>(rng() % kRows);
    const auto j = static_cast<samaya::Index>(rng() % kCols);
    double a = static_cast<double>(static_cast<int>(rng() % 17) - 8);
    const auto kind = rng() % 3;
    const bool present = dense[i][j] != 0.0;
    if (kind == 0) {
      if (a == 0.0) a = 1.0;
      CHECK_EQ(Access::add(p, i, j, a), !present);
      if (!present) dense[i][j] = a;
    } else if (kind == 1) {
      CHECK_EQ(Access::change(p, i, j, a), present);
      if (present) dense[i][j] = a;
    } else {
      CHECK_EQ(Access::remove(p, i, j), present);
      dense[i][j] = 0.0;
    }
    Access::check(p, dense);
  }

  // Deactivation keeps stored coefficients in both directions but counts active entries only.
  Access::deactivate_row(p, 2);
  Access::check(p, dense);
  CHECK(!Access::add(p, 2, 4, 5.0));
  CHECK(!Access::change(p, 2, 4, 5.0));
  Access::deactivate_col(p, 5);
  Access::check(p, dense);
  CHECK(!Access::add(p, 4, 5, 5.0));
  CHECK(!Access::change(p, 4, 5, 5.0));
  // Physical removal of inactive entries must not decrement counts a second time.
  for (samaya::Index j = 0; j < kCols; ++j) {
    CHECK_EQ(Access::remove(p, 2, j), dense[2][j] != 0.0);
    dense[2][j] = 0.0;
    Access::check(p, dense);
  }
  for (samaya::Index i = 0; i < kRows; ++i) {
    CHECK_EQ(Access::remove(p, i, 5), dense[i][5] != 0.0);
    dense[i][5] = 0.0;
    Access::check(p, dense);
  }
}
