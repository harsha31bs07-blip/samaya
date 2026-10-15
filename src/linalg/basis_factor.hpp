#pragma once

#include <array>
#include <chrono>
#include <cstdint>
#include <span>
#include <utility>
#include <vector>

#include "samaya/sparse_matrix.hpp"

namespace samaya {

// Integer observations (including zero), used for solve densities and factor fill.
struct MeasurementCount {
  long long count = 0;
  long long sum = 0;
  long long maximum = 0;
  void add(long long value) {
    ++count;
    sum += value;
    if (value > maximum) maximum = value;
  }
};

// Sparse LU factorization of a simplex basis with Forrest–Tomlin updates.
//
// The basis matrix B consists of m columns of [A  -I]: basis entry j < n is column j of A and
// j >= n is the logical column -e_(j-n). "Positions" 0..m-1 index the columns of B; "rows" index
// the rows of A.
//
// Factorization: Markowitz pivot selection with threshold partial pivoting. Column and row
// singletons have Markowitz cost zero and are always taken first, so the triangular part of the
// basis is eliminated before the nucleus without fill-in.
//
// Representation: B^-1 = U^-1 R L, where L is a product of column etas from elimination, R a
// product of row etas from Forrest–Tomlin updates, and U is upper triangular up to a symmetric
// permutation that changes as updates move pivots to the end of the pivot order.
class BasisFactor {
 public:
  enum class RebuildReason : std::uint8_t {
    kInitial, kPhase, kInterval, kFillLimit, kUnsafeUpdate, kMismatch,
    kNumerical, kConfirmation, kRankRepair, kUpdateWork, kCount
  };
  struct Stats {
    long long rebuilds = 0;  // Factorization attempts, including rank-repair retries.
    std::array<long long, static_cast<std::size_t>(RebuildReason::kCount)> reasons{};
    MeasurementCount l_nnz, u_nnz, r_nnz;  // After factorization; U includes its diagonal.
    MeasurementCount old_l_nnz, old_u_nnz, old_r_nnz;  // Immediately before clearing updates.
    double rebuild_seconds = 0.0;
    double max_rebuild_seconds = 0.0;
    // Solve work, in loop iterations and entries touched per part (deterministic, so it shows
    // where FTRAN and BTRAN spend their time without timers). L etas and U slots count every
    // one visited, including the ones skipped for a zero pivot entry.
    long long ftran_calls = 0, ftran_l_etas = 0, ftran_l_entries = 0, ftran_r_entries = 0;
    long long ftran_u_slots = 0, ftran_u_entries = 0;
    long long btran_calls = 0, btran_u_slots = 0, btran_u_entries = 0, btran_r_entries = 0;
    long long btran_l_entries = 0;
  };
  struct Options {
    double pivot_threshold = 0.1;  // Relative threshold against the column max.
    double abs_pivot_tol = 1e-11;  // Smaller pivots count as zero (rank deficiency).
    double drop_tol = 1e-14;       // Fill-in below this magnitude is dropped.
    int search_limit = 4;          // Markowitz search stops after this many candidate lines.
    int max_updates = 100;         // Refactorization interval (an upper bound with rebuild_cost).
    // Amortized refactorization: rebuild once the extra solve work caused by the updates since
    // the last factorization (3 solves per iteration over the entries added to R and U) exceeds
    // the estimated cost of a factorization, rebuild_cost x its factor entries in units of one
    // solve's cost per entry. Step-0 profiles (SIMPLEX_SPEED_PLAN.md) measured that ratio at
    // about 70 on most Netlib models and 150-340 on the large ones (pilot, pilot87, fit2p, the
    // 467k-row LP); the cost near the optimal interval is flat (the interval goes with the square
    // root of the ratio), so one constant serves. 0 disables the rule.
    double rebuild_cost = 200.0;
    bool collect_stats = false;
  };

  explicit BasisFactor(const SparseMatrix& A) : BasisFactor(A, Options{}) {}
  BasisFactor(const SparseMatrix& A, Options options);

  // Factorizes the basis. Returns the rank deficiency: 0 on success. Otherwise the factor is not
  // usable and singular_positions() / unpivoted_rows() (equal length) say which basis columns
  // could not be pivoted and which rows were left; replacing basis position
  // singular_positions()[k] with the logical of unpivoted_rows()[k] gives a nonsingular basis.
  Index factorize(std::span<const Index> basic, RebuildReason reason = RebuildReason::kInitial);

  const std::vector<Index>& singular_positions() const { return singular_positions_; }
  const std::vector<Index>& unpivoted_rows() const { return unpivoted_rows_; }

  // Solves B z = x. Input is indexed by row, output by basis position (in place). If spike is
  // non-null it receives the partially transformed vector R L x needed by update().
  // Optional density observation is counted as output entries are written, with no extra scan.
  void ftran(std::vector<double>& x, std::vector<double>* spike = nullptr,
             MeasurementCount* nonzeros = nullptr) const;

  // Solves B^T z = c. Input is indexed by basis position, output by row (in place).
  void btran(std::vector<double>& c) const;

  // Replaces the column at basis position `position` by the column whose ftran spike is given.
  // Returns false when the update is numerically unsafe; the factor is then invalid and
  // factorize() must be called with the new basis.
  bool update(Index position, const std::vector<double>& spike);

  int num_updates() const { return num_updates_; }
  bool should_refactor() const;
  bool update_work_exceeds_rebuild() const;
  Index rows() const { return m_; }
  NnzIndex factor_nnz() const;
  RebuildReason refactor_reason() const;
  const Stats& stats() const { return stats_; }
  void reset_stats() { stats_ = {}; }

 private:
  // Records every return path, including a rank-deficient factorization.
  struct RebuildMeasurement {
    BasisFactor& factor;
    std::chrono::steady_clock::time_point start{};
    RebuildMeasurement(BasisFactor& factor, RebuildReason reason);
    ~RebuildMeasurement();
  };
  using Entry = std::pair<Index, double>;

  void clear();
  static void erase_position(std::vector<Entry>& list, Index position);

  const SparseMatrix& A_;
  Options options_;
  mutable Stats stats_;  // The solves are const but count their work.
  Index m_;
  Index n_;

  // L column etas: eta e has pivot row l_pivot_row_[e] and entries l_start_[e]..l_start_[e+1].
  std::vector<Index> l_pivot_row_;
  std::vector<NnzIndex> l_start_;
  std::vector<Index> l_index_;
  std::vector<double> l_value_;

  // Forrest–Tomlin row etas: x[pivot] -= sum value * x[index].
  std::vector<Index> r_pivot_row_;
  std::vector<NnzIndex> r_start_;
  std::vector<Index> r_index_;
  std::vector<double> r_value_;

  // U, indexed by basis position on both sides. Entry (k, u) in ucol_[r] means row prow_[k] has
  // coefficient u in column r; urow_ holds the same entries by row.
  std::vector<Index> prow_;
  std::vector<Index> pos_of_row_;
  std::vector<double> diag_;
  std::vector<std::vector<Entry>> ucol_;
  std::vector<std::vector<Entry>> urow_;
  std::vector<Index> order_;       // Pivot order; -1 marks a pivot moved to the end.
  std::vector<Index> order_slot_;  // Slot of each position in order_.
  NnzIndex u_nnz_ = 0;
  NnzIndex initial_nnz_ = 0;

  int num_updates_ = 0;
  double extra_solve_work_ = 0.0;  // Accumulated solve entries added by updates (see Options).
  bool valid_ = false;
  std::vector<Index> singular_positions_;
  std::vector<Index> unpivoted_rows_;

  mutable std::vector<double> work_;  // Scratch for ftran/btran; contents undefined.
  std::vector<double> update_work_;   // All zero between calls to update().
};

}  // namespace samaya
