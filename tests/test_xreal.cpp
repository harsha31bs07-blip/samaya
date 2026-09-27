// Xreal (src/core/xreal.hpp): the double-double type that replaces long double, so that the
// verifier computes the same bits on every platform. The checks are exact, not approximate:
// integer arithmetic where the result fits in 64 bits, std::fma (correctly rounded, hence exact
// for a product's rounding error) elsewhere. None relies on long double, which is plain double
// with MSVC.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <utility>
#include <vector>

#include "core/xreal.hpp"
#include "samaya/verify.hpp"
#include "test_framework.hpp"

using samaya::Xreal;

TEST(xreal_two_sum_is_exact_on_integers) {
  // a and b are integer-valued doubles below 2^60, so their exact values, a + b and its rounding
  // error are integers that fit in int64 and s + e == a + b can be checked exactly. (The exact
  // values come from the doubles themselves: MSVC at /O2 folds a round-trip check such as
  // (int64)(double)i == i to true even when i is not representable.)
  std::mt19937_64 rng(7);
  for (int k = 0; k < 20000; ++k) {
    const int shift = static_cast<int>(rng() % 60);
    const std::int64_t ia = static_cast<std::int64_t>(rng() >> (4 + shift)) * (k % 2 ? -1 : 1);
    const std::int64_t ib = static_cast<std::int64_t>(rng() >> (4 + (k % 60)));
    const double a = static_cast<double>(ia);  // rounded to an integer when |ia| > 2^53
    const double b = static_cast<double>(ib);
    double s = 0.0, e = 0.0;
    Xreal::two_sum(a, b, s, e);
    CHECK_EQ(s, a + b);
    const std::int64_t exact = static_cast<std::int64_t>(a) + static_cast<std::int64_t>(b);
    CHECK_EQ(static_cast<std::int64_t>(s) + static_cast<std::int64_t>(e), exact);
  }
}

TEST(xreal_two_product_matches_fma) {
  // fma(a, b, -p) is a*b - p rounded once; that difference is representable, so fma returns it
  // exactly and must equal Dekker's error term.
  std::mt19937_64 rng(11);
  std::uniform_real_distribution<double> mant(-1.0, 1.0);
  std::uniform_int_distribution<int> expo(-300, 300);
  int checked = 0;
  for (int k = 0; k < 20000; ++k) {
    const double a = std::ldexp(mant(rng), expo(rng));
    const double b = std::ldexp(mant(rng), expo(rng) / 2);
    double p = 0.0, e = 0.0;
    Xreal::two_product(a, b, p, e);
    CHECK_EQ(p, a * b);
    // Where a*b underflows toward subnormals the error is not representable; skip those.
    if (std::fabs(p) < 1e-290) continue;
    CHECK_EQ(e, std::fma(a, b, -p));
    ++checked;
  }
  CHECK(checked > 15000);
  // Magnitudes above 2^996 take the scaled split.
  for (const double a : {1.7e308, -1.5e307, 8.98846567431158e307}) {
    const double b = 0.4999999999999999;
    double p = 0.0, e = 0.0;
    Xreal::two_product(a, b, p, e);
    CHECK_EQ(p, a * b);
    CHECK_EQ(e, std::fma(a, b, -p));
  }
}

TEST(xreal_recovers_cancelled_sums_that_double_loses) {
  Xreal x = 1e16;
  x += 1.0;
  x += -1e16;
  CHECK_EQ(static_cast<double>(x), 1.0);
  // Products that cancel in pairs, interleaved with small integers: the exact total is the
  // sum of the small integers. The products are integers near 2^80 (60 significant bits, so
  // they need Dekker's exact product) and every partial sum stays within 106 bits, so the
  // double-double must be exact where plain double loses the small terms.
  std::mt19937_64 rng(3);
  std::vector<std::pair<double, double>> terms;
  double expected = 0.0;
  for (int k = 0; k < 200; ++k) {
    const double a = std::ldexp(static_cast<double>((rng() >> 34) | (1ULL << 29)), 20);
    const double b = static_cast<double>((rng() >> 34) | (1ULL << 29));
    terms.push_back({a, b});
    terms.push_back({-a, b});
    const double small = static_cast<double>(static_cast<int>(rng() % 7) - 3);
    terms.push_back({small, 1.0});
    expected += small;
  }
  std::shuffle(terms.begin(), terms.end(), rng);
  Xreal dot;
  double plain = 0.0;
  for (const auto& [a, b] : terms) {
    dot += Xreal(a) * b;
    plain += a * b;
  }
  CHECK_EQ(static_cast<double>(dot), expected);
  CHECK(plain != expected);  // the case really needs the extra precision
}

TEST(xreal_infinities_and_comparisons) {
  Xreal lo = 1.0;
  lo += -std::numeric_limits<double>::infinity();
  CHECK(std::isinf(static_cast<double>(lo)) && static_cast<double>(lo) < 0);
  lo += 5.0;
  CHECK(std::isinf(static_cast<double>(lo)));
  Xreal hi = std::numeric_limits<double>::infinity();
  CHECK(hi > 1e308);
  CHECK(!(hi < 0.0));
  Xreal t = 1e16;
  t += 1.0;  // hi = 1e16, lo = 1: greater than 1e16 as a double
  CHECK(t > 1e16);
  CHECK(fabs(-t) == t);
  CHECK(Xreal(0.0) == 0.0);
  CHECK(-Xreal(2.0) < 0);
}

TEST(xreal_verifier_objective_keeps_cancelled_terms) {
  // obj = 1e16 x0 + x1 - 1e16 x2 at x = (1, 1, 1): exactly 1. With 53-bit sums (long double on
  // MSVC) it would come out as 0 or 2.
  samaya::Model m;
  m.obj = {1e16, 1.0, -1e16};
  m.col_lower = {0, 0, 0};
  m.col_upper = {2, 2, 2};
  m.col_type = {samaya::VarType::kContinuous, samaya::VarType::kContinuous,
                samaya::VarType::kContinuous};
  m.A = samaya::SparseMatrix::from_triplets(0, 3, {});
  const samaya::VerifyReport r = samaya::verify_primal(m, std::vector<double>{1, 1, 1});
  CHECK(r.ok);
  CHECK_EQ(r.objective, 1.0);
}
