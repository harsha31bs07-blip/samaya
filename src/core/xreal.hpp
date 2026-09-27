#pragma once

// Extended precision that is the same on every platform and compiler: a double-double
// (hi + lo, about 106 significant bits), built from IEEE double +, - and * only. Sums use
// Knuth's TwoSum and products Dekker's TwoProduct with Veltkamp splitting, so no FMA is needed
// and the result does not depend on the CPU. It replaces `long double`, which is 80-bit (64-bit
// significand) with GCC on x86-64 but plain double with MSVC.
//
// The error terms are exact only under IEEE arithmetic: these routines must not be compiled
// with -ffast-math or FP contraction (CMakeLists.txt sets -ffp-contract=off and /fp:precise).
// Infinities pass through (hi = +-inf, lo = 0); NaN propagates in hi.

#include <cmath>

namespace samaya {

class Xreal {
 public:
  constexpr Xreal() = default;
  constexpr Xreal(double v) : hi_(v), lo_(0.0) {}  // NOLINT: implicit, like a wider float type

  double hi() const { return hi_; }
  double lo() const { return lo_; }
  explicit operator double() const { return hi_ + lo_; }

  // s + e == a + b exactly (Knuth), for any finite a and b.
  static void two_sum(double a, double b, double& s, double& e) {
    s = a + b;
    const double bb = s - a;
    e = (a - (s - bb)) + (b - bb);
  }

  // p + e == a * b exactly (Dekker), for finite a, b whose product neither overflows nor
  // underflows.
  static void two_product(double a, double b, double& p, double& e) {
    p = a * b;
    if (!std::isfinite(p)) {
      e = 0.0;
      return;
    }
    double ah = 0.0, al = 0.0, bh = 0.0, bl = 0.0;
    split(a, ah, al);
    split(b, bh, bl);
    e = ((ah * bh - p) + ah * bl + al * bh) + al * bl;
  }

  Xreal& operator+=(const Xreal& o) {
    if (!std::isfinite(hi_) || !std::isfinite(o.hi_)) return set_nonfinite(hi_ + o.hi_);
    double s = 0.0, e = 0.0, t = 0.0, f = 0.0;
    two_sum(hi_, o.hi_, s, e);
    two_sum(lo_, o.lo_, t, f);
    e += t;
    quick_two_sum(s, e, s, e);
    e += f;
    quick_two_sum(s, e, hi_, lo_);
    if (!std::isfinite(hi_)) lo_ = 0.0;
    return *this;
  }
  Xreal& operator-=(const Xreal& o) { return *this += -o; }

  // Multiplication by a double: exact when this is a double (lo = 0), since then it is Dekker's
  // product; otherwise the usual double-double product (error about 2^-104 relative).
  Xreal& operator*=(double b) {
    if (!std::isfinite(hi_) || !std::isfinite(b)) return set_nonfinite(hi_ * b);
    double p = 0.0, e = 0.0;
    two_product(hi_, b, p, e);
    if (!std::isfinite(p)) return set_nonfinite(p);
    e += lo_ * b;
    quick_two_sum(p, e, hi_, lo_);
    return *this;
  }

  Xreal operator-() const {
    Xreal r;
    r.hi_ = -hi_;
    r.lo_ = -lo_;
    return r;
  }

  friend Xreal operator+(Xreal a, const Xreal& b) { return a += b; }
  friend Xreal operator-(Xreal a, const Xreal& b) { return a -= b; }
  friend Xreal operator*(Xreal a, double b) { return a *= b; }
  friend Xreal operator*(double a, Xreal b) { return b *= a; }

  // Comparisons on the exact value hi + lo (hi carries the rounded value, lo the rest).
  friend bool operator==(const Xreal& a, const Xreal& b) {
    return a.hi_ == b.hi_ && a.lo_ == b.lo_;
  }
  friend bool operator!=(const Xreal& a, const Xreal& b) { return !(a == b); }
  friend bool operator<(const Xreal& a, const Xreal& b) {
    return a.hi_ < b.hi_ || (a.hi_ == b.hi_ && a.lo_ < b.lo_);
  }
  friend bool operator>(const Xreal& a, const Xreal& b) { return b < a; }
  friend bool operator<=(const Xreal& a, const Xreal& b) { return !(b < a); }
  friend bool operator>=(const Xreal& a, const Xreal& b) { return !(a < b); }

  friend Xreal fabs(const Xreal& a) {
    return a.hi_ < 0.0 || (a.hi_ == 0.0 && a.lo_ < 0.0) ? -a : a;
  }

 private:
  // Veltkamp's split: a == hi + lo with hi holding the top 26 bits. Magnitudes above 2^996
  // are scaled down by an exact power of two first, or 2^27 * a would overflow.
  static void split(double a, double& hi, double& lo) {
    constexpr double kSplitter = 134217729.0;              // 2^27 + 1
    constexpr double kSplitLimit = 6.69692879491417e+299;  // 2^996
    constexpr double kDown = 3.7252902984619140625e-09;    // 2^-28
    constexpr double kUp = 268435456.0;                    // 2^28
    if (a > kSplitLimit || a < -kSplitLimit) {
      a *= kDown;
      const double t = kSplitter * a;
      hi = t - (t - a);
      lo = a - hi;
      hi *= kUp;
      lo *= kUp;
      return;
    }
    const double t = kSplitter * a;
    hi = t - (t - a);
    lo = a - hi;
  }

  // s + e == a + b exactly, given |a| >= |b| (or a == 0).
  static void quick_two_sum(double a, double b, double& s, double& e) {
    s = a + b;
    e = b - (s - a);
  }

  Xreal& set_nonfinite(double v) {
    hi_ = v;
    lo_ = 0.0;
    return *this;
  }

  double hi_ = 0.0;
  double lo_ = 0.0;
};

}  // namespace samaya
