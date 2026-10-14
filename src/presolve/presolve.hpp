#pragma once

#include <cstdint>
#include <vector>

#include "samaya/model.hpp"

namespace samaya {

enum class PresolveStatus : std::uint8_t {
  kReduced,     // reduced() holds an equivalent (possibly empty) model.
  kInfeasible,  // A reduction proved the model infeasible.
};

struct PresolveOptions {
  double feasibility_tol = 1e-9;
  double integrality_tol = 1e-6;
  // Keep reductions valid for integer columns (exact divisibility for substitution, integral
  // bounds). Postsolve then restores primal values only.
  bool mip = false;
  // Implied-free column substitution through equations (C3). Off only for A/B measurements
  // and for tests that isolate an earlier rule.
  bool substitute = true;
  int max_passes = 50;
};

struct PresolveStats {
  int passes = 0;
  Index rows_removed = 0;
  Index cols_removed = 0;
  Index empty_rows = 0;
  Index singleton_rows = 0;
  Index forcing_rows = 0;
  Index redundant_rows = 0;
  Index duplicate_rows = 0;
  Index fixed_cols = 0;
  Index empty_cols = 0;
  Index dominated_cols = 0;
  Index free_col_singletons = 0;
  Index implied_free_col_singletons = 0;
  Index implied_free_rejected_by_tolerance = 0;
  Index bounds_tightened = 0;
  Index doubleton_continuous_cols = 0;
  Index doubleton_integer_cols = 0;
  Index doubleton_rejected_pivot = 0;
  Index doubleton_rejected_integer_data = 0;
  Index doubleton_rejected_divisibility = 0;
  Index doubleton_rejected_integer_bound = 0;
  Index doubleton_rejected_crossed_bounds = 0;
  Index doubleton_rejected_residual = 0;
  Index doubleton_rejected_nonfinite = 0;
  Index substitution_continuous_cols = 0;
  Index substitution_integer_cols = 0;
  Index substitution_rejected_row_length = 0;
  Index substitution_rejected_col_length = 0;
  Index substitution_rejected_not_implied_free = 0;
  Index substitution_rejected_by_tolerance = 0;
  Index substitution_rejected_pivot = 0;
  Index substitution_rejected_integer = 0;
  Index substitution_rejected_nonfinite = 0;
  Index substitution_rejected_residual = 0;
  Index substitution_rejected_fill = 0;
};

// LP / MILP presolve: removes rows and columns, then maps a solution of
// the reduced model back to the original.
//
// Reductions: empty rows and columns, fixed columns, singleton rows (turned into column bounds),
// forcing and redundant rows (from row activity bounds), dominated columns (a column whose cost
// and coefficients all push it to one finite bound is fixed there), free or implied-free column
// singletons (substituted out of an equality row, or absorbing an inequality row when the column
// has no cost), duplicate (parallel) rows, doubleton equations (MIP only), and implied-free
// column substitution through equations with limited fill.
//
// Postsolve undoes the reductions in reverse order and restores primal values and, in LP mode,
// row duals: each step keeps the solution optimal for the model as it was before that reduction,
// so the result satisfies the optimality conditions of the original model. Reduced costs follow
// from the original model as c - A'y.
class Presolve {
 public:
  Presolve(const Model& model, const PresolveOptions& options = {});

  PresolveStatus run();

  const Model& reduced() const { return reduced_; }
  const PresolveStats& stats() const { return stats_; }
  // Original index of each reduced column / row.
  const std::vector<Index>& col_map() const { return col_map_; }
  const std::vector<Index>& row_map() const { return row_map_; }

  // Maps a reduced solution (x: reduced columns, y: reduced row duals in the convention
  // c - A'y of the reduced model; may be empty for MIP) to the original model.
  void postsolve(const std::vector<double>& reduced_x, const std::vector<double>& reduced_y,
                 std::vector<double>& x, std::vector<double>& y) const;

 private:
  friend struct PresolveTestAccess;
  friend struct PresolveSubstitutionTestAccess;
  struct Entry {
    Index index;  // Column in rows_, row in cols_.
    double a;
  };
  enum class Kind : std::uint8_t {
    kFixCol,          // Column fixed at `value` (fixed, empty or dominated).
    kRemoveRow,       // Empty or redundant row: y = 0.
    kSingletonRow,    // Row turned into bounds on `col`.
    kForcingRow,      // Row forcing its columns to the bounds in `terms` (lower: activity min).
    kSubstituteCol,   // Free singleton column solved from equality `row`.
    kSlackCol,        // Free cost-free singleton column absorbing inequality `row`.
    kDuplicateRow,    // `row2` merged into `row` with a_row2 = factor * a_row.
    kDoubleton,       // MIP only: a * col + b * survivor = rhs_lower.
    kSubstitution,    // Implied-free col solved from row and substituted through other rows.
  };
  struct Term {
    Index col;
    double a;
    double cost;   // Column cost when the reduction was made.
    double value;  // Forcing rows: the bound the column is fixed at.
    bool boxed;    // Forcing rows: the column had distinct finite-or-not bounds (lower < upper).
    std::size_t column_first = 0;  // Column snapshot in column_terms_.
    std::size_t column_last = 0;
  };
  struct Record {
    Kind kind;
    Index row = -1;
    Index col = -1;
    Index row2 = -1;
    double a = 0.0;
    double value = 0.0;
    double cost = 0.0;
    double old_lower = 0.0;
    double old_upper = 0.0;
    double rhs_lower = 0.0;
    double rhs_upper = 0.0;
    bool flag1 = false;  // Forcing: at activity minimum. Duplicate: lower bound from row2.
    bool flag2 = false;  // Duplicate: upper bound from row2.
    // [first, last) in column_terms_ for singleton rows / doubletons, otherwise terms_.
    std::size_t first = 0;
    std::size_t last = 0;
    std::size_t column_first = 0;  // Substitution: col snapshot excluding row.
    std::size_t column_last = 0;
    Index survivor = -1;
    double b = 0.0;
    // Doubletons retain pre-transfer bounds and the staged post-substitution column and cost
    // so the later LP dual implementation can use the same saved reduction state.
    double survivor_lower = 0.0;
    double survivor_upper = 0.0;
    double survivor_cost = 0.0;
    std::size_t survivor_first = 0;
    std::size_t survivor_last = 0;
  };

  bool is_int(Index j) const { return model_.col_type[j] == VarType::kInteger; }
  double tol(double v) const;

  // Sorted, unique, nonzero coefficients in both directions. Insertion/change is refused for
  // inactive rows or columns; changing to zero removes the entry. Deactivation is logical:
  // both lists retain the same stored entries, while counts include only active endpoints.
  bool add_entry(Index i, Index j, double a);
  bool change_coefficient(Index i, Index j, double a);
  bool remove_entry(Index i, Index j);
  void snapshot_column(Index j, std::size_t& first, std::size_t& last);

  void remove_row(Index i);
  void remove_col(Index j);
  void fix_col(Index j, double value);
  bool set_bounds(Index j, double lower, double upper);  // False if crossed (infeasible).
  void activity_bounds(Index i, double& min_act, double& max_act) const;
  struct ActivitySide {
    double sum = 0.0;
    double magnitude = 0.0;
    Index infinite = 0;
  };
  void other_activity_bounds(Index i, Index j, ActivitySide& lower, ActivitySide& upper) const;
  bool implied_free(Index i, Index j, double a, bool& tolerance_only) const;

  bool process_cols(bool& changed);
  bool process_rows(bool& changed);
  bool process_singleton_cols(bool& changed);
  bool process_duplicate_rows(bool& changed);
  bool process_doubleton_rows(bool& changed);
  bool process_substitutions(bool& changed);
  void build_reduced();

  const Model& model_;
  PresolveOptions options_;
  Index m_;
  Index n_;
  double sense_;
  std::vector<std::vector<Entry>> rows_;
  std::vector<std::vector<Entry>> cols_;

  std::vector<double> cost_;  // Minimization costs.
  std::vector<double> lower_;
  std::vector<double> upper_;
  std::vector<double> row_lower_;
  std::vector<double> row_upper_;
  double offset_ = 0.0;
  std::vector<char> row_active_;
  std::vector<char> col_active_;
  std::vector<Index> row_count_;
  std::vector<Index> col_count_;

  std::vector<Record> records_;
  std::vector<Term> terms_;
  std::vector<Entry> column_terms_;  // Active column coefficients at the time of a reduction.

  Model reduced_;
  std::vector<Index> col_map_;
  std::vector<Index> row_map_;
  PresolveStats stats_;
};

}  // namespace samaya
