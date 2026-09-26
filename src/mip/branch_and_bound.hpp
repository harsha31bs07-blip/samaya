#pragma once

#include <map>
#include <memory>
#include <optional>
#include <random>
#include <vector>

#include "core/log.hpp"
#include "linalg/scaling.hpp"
#include "lp/simplex.hpp"
#include "mip/cuts.hpp"
#include "samaya/model.hpp"
#include "samaya/status.hpp"

namespace samaya {

struct MipOptions {
  double time_limit = kInf;
  long long node_limit = -1;  // -1: unlimited.
  double rel_gap = 1e-4;
  double abs_gap = 1e-6;
  double integrality_tol = 1e-6;
  double feasibility_tol = 1e-6;  // Row and bound tolerance for accepting a solution.
  // Pseudocosts count as reliable after this many observations in each direction; before that a
  // candidate is strong-branched.
  int reliability = 4;
  int max_strong_branching = 8;       // Strong-branching candidates per node.
  long long strong_iterations = 200;  // Dual simplex iterations per strong-branching child.
  // Beyond the soft limit of open nodes the search switches to depth first (a LIFO stack), which
  // stops the open list from growing; beyond the hard limit it stops with kNodeLimit instead of
  // exhausting memory. An open node costs about 400 bytes, so the defaults stay near 2 GB.
  std::size_t max_open_nodes_soft = 2500000;
  std::size_t max_open_nodes = 5000000;
  // Root cutting planes (Gomory mixed-integer, c-MIR, lifted knapsack covers).
  bool cuts = true;
  int max_cut_rounds = 20;
  int max_cuts_per_round = 200;
  // Cuts in the tree: pool cuts violated at a node, and fresh c-MIR and cover cuts separated with
  // the root bounds (so valid everywhere) at every tree_separation_frequency-th depth; a node's
  // cuts stay only if they close enough of its gap. Sequential search only: node bases then keep
  // one row layout. SCIP's default frequency for these separators is 10.
  bool tree_cuts = true;
  int tree_separation_frequency = 10;
  // Conflict analysis: proofs from infeasible and cut-off node LPs are stored and propagated. A
  // proof may involve at most conflict_density x columns (+ 10): a dense one rarely prunes and
  // costs every node a pass over its columns (SCIP's conflict/maxvarsfac is 0.15).
  bool conflicts = true;
  double conflict_density = 0.15;
  // Kept tree cuts stay in every later node LP, so they may grow the LP by at most this share of
  // its rows: a node LP's cost grows about linearly with the rows, and on small models
  // (markshare_4_0: 6 rows) node throughput is what solves them.
  double tree_cut_row_fraction = 0.25;
  // Primal heuristics: feasibility pump, diving, and the RENS/RINS sub-MIPs (solved by a nested
  // search on a presolved copy; a nested search never starts sub-MIPs itself).
  bool heuristics = true;
  // Tree search threads. The search starts sequentially; once it has processed
  // parallel_start_nodes nodes and the tree is still open, the open nodes go to a shared pool
  // served by this many threads, each with its own copy of the LP (the root, its cuts and
  // bounds are shared). The run is then not deterministic. 1 keeps the sequential search.
  int threads = 1;
  long long parallel_start_nodes = 200;
  // Root reductions for integer columns: coefficient tightening and probing on binaries.
  bool probing = true;
  // Restart once after the root when it fixed a large share of the integer columns.
  bool restart = true;
  bool sub_mip_heuristics = true;
  // A known point in the model's columns, e.g. the previous plan when re-planning; NaN marks an
  // unknown value. It becomes the incumbent if it is feasible; otherwise its integer values are
  // fixed and the continuous columns re-solved; otherwise it is Feasibility Jump's starting point,
  // which repairs it.
  std::vector<double> start;
  // Only solutions strictly better than this objective (in the model's sense) are accepted; the
  // search prunes against it as if it were an incumbent. Used by the sub-MIPs.
  std::optional<double> objective_cutoff;
  // Tests: a known feasible (e.g. optimal) solution. Every cut is checked against it and a cut
  // that separates it is counted in MipOutcome::debug_cut_violations. Setting it also checks the
  // Farkas ray of every infeasible node LP, kept as a proof or not, with verify_infeasibility; a
  // missing ray counts too (debug_farkas_failures).
  std::vector<double> debug_solution;
};

struct MipOutcome {
  // kOptimal (gap closed), kInfeasible, kInfeasibleOrUnbounded, kUnbounded, kTimeLimit,
  // kNodeLimit or kNumericalError (a node LP could not be solved, so the search is incomplete).
  Status status = Status::kNotSolved;
  double objective = kInf;  // Of `x`, in the model's sense; +-inf without a solution.
  double bound = -kInf;     // Best proven bound, in the model's sense.
  std::vector<double> x;    // Best integer solution found (empty if none).
  long long nodes = 0;
  long long lp_iterations = 0;
  long long strong_branching_iterations = 0;
  int heuristic_solutions = 0;
  int coefficients_tightened = 0;
  int probing_fixed = 0;
  int probing_tightened = 0;
  long long heuristic_lp_iterations = 0;  // Diving and feasibility-pump LPs (part of lp_iterations).
  int threads_used = 1;
  long long reduced_cost_fixings = 0;  // Column bounds tightened by reduced costs.
  bool restarted = false;              // The search restarted after the root.
  int cut_rounds = 0;
  int cuts_added = 0;             // Cuts in the LP after the root (non-binding ones removed).
  int tree_cuts = 0;              // Tree cuts kept in the LP (pool and fresh).
  long long tree_cuts_separated = 0;  // Tree cuts appended, including those removed again.
  long long conflicts_found = 0;      // Proofs kept in the conflict pool.
  long long conflict_prunes = 0;      // Nodes pruned by a stored proof before their LP.
  long long conflict_infeasible_lps = 0;  // Infeasible node LPs whose Farkas ray was analyzed.
  long long farkas_proofs = 0;            // Of conflicts_found: proofs from Farkas rays.
  long long debug_conflict_violations = 0;  // Proofs that exclude debug_solution (must be 0).
  // Farkas rays of infeasible node LPs (kept as a proof or not) that, mapped to the original rows,
  // verify_infeasibility rejects on their node's bounds, or missing rays (must be 0; checked only
  // with debug_solution).
  long long debug_farkas_failures = 0;
  double root_bound = -kInf;      // Root LP bound before and after cuts, in the model's sense.
  double root_bound_cuts = -kInf;
  long long debug_cut_violations = 0;
  // Root bounds or tightened rows that exclude debug_solution (must stay 0).
  long long debug_reduction_violations = 0;
};

// LP-based branch-and-bound for mixed-integer linear programs.
//
// Node relaxations are solved by the dual simplex, warm-started from the parent's optimal basis,
// on one scaled copy of the LP whose column bounds change from node to node. Each node first runs
// bound propagation on the integer columns. Branching uses pseudocosts with reliability
// initialization by strong branching (the product score); strong branching that proves a child
// infeasible or cut off tightens the node instead. Nodes are selected by best bound, with plunging
// (a depth-first dive into a child) while the dive stays promising. Heuristics: lock-based simple
// rounding at every node, and rounding followed by an LP over the continuous columns at the root
// and periodically. A candidate solution is accepted only if it satisfies the rows, bounds and
// integrality within the tolerances.
class BranchAndBound {
 public:
  BranchAndBound(const Model& model, const MipOptions& options, const Logger& log);
  ~BranchAndBound();

  MipOutcome solve();

 private:
  struct BoundChange {
    Index col;
    double lower;
    double upper;
  };
  // Bound changes are stored as a tree shared by the nodes: a node's bounds are the root bounds
  // plus the changes of every segment from the root down, then its own `local` changes. Memory
  // is proportional to the changes, not to nodes x depth.
  struct PathSegment {
    std::shared_ptr<const PathSegment> parent;
    std::vector<BoundChange> changes;
  };
  struct Node {
    double bound = -kInf;  // Lower bound (minimization) inherited from the parent.
    int depth = 0;
    std::shared_ptr<const PathSegment> path;
    std::vector<BoundChange> local;  // The branching, then propagation at this node.
    std::shared_ptr<const std::vector<VarStatus>> basis;
    Index branch_col = -1;  // The branching that created this node, for pseudocosts.
    bool branch_up = false;
    double branch_distance = 0.0;
  };
  enum class NodeResult : std::uint8_t { kPruned, kBranched, kStopped, kFailed, kUnbounded };

  // Relaxation.
  void build_relaxation();
  void set_bound(Index j, double lower, double upper);
  void apply_node_bounds(const Node& node);
  SimplexStatus solve_relaxation(const std::vector<VarStatus>* start, long long iteration_limit);
  double relaxation_objective() const;
  std::vector<VarStatus> current_basis() const;

  MipOutcome restart();

  // Parallel tree search (parallel.cpp).
  struct Shared;
  void run_parallel(bool& unbounded, bool& stopped, bool& gap_closed);
  void adopt(const BranchAndBound& master, int id);
  void worker_loop(Shared& shared);
  void publish_incumbent();
  void pull_incumbent();

  // Root reductions (probing.cpp).
  int tighten_coefficients();
  bool probe();

  // Heuristics (heuristics.cpp).
  void run_feasibility_jump();
  enum class DiveRule : std::uint8_t { kFractional, kCoefficient, kPseudocost, kGuided };
  void run_heuristics(const Node& node, const std::vector<double>& x,
                      const std::vector<VarStatus>& basis);
  void dive(DiveRule rule, const std::vector<double>& x, const std::vector<VarStatus>& basis,
            long long budget);
  void feasibility_pump(const std::vector<double>& x, const std::vector<VarStatus>& basis,
                        long long budget);
  void rens(const std::vector<double>& x);
  void rins(const std::vector<double>& x);
  void sub_mip(const std::vector<double>& lower, const std::vector<double>& upper,
               const char* name);
  void undo_bounds(std::size_t mark);

  // Cuts.
  void root_cut_loop();
  // Separates cuts at the LP point x_ and appends the selected ones as rows. At the root (tree
  // false) with the current bounds, Gomory cuts included. In the tree with the root bounds,
  // without Gomory cuts (their tableau rows are only locally valid), plus violated pool cuts,
  // at most tree_max_cuts; fresh separation only if `fresh`.
  int separate_and_add_cuts(bool tree = false, bool fresh = true, int tree_max_cuts = 0);
  // Gomory mixed-integer cuts from the current LP's tableau rows (valid under ctx's bounds).
  void separate_gomory(const CutContext& ctx, std::vector<Cut>& candidates);
  // Adds tree cuts at the current node when allowed; returns the number of rows appended.
  int tree_cut_round(int depth);
  void pool_rows(const std::vector<Index>& rows);
  void remove_rows(const std::vector<Index>& rows);

  // Conflict analysis (conflicts.cpp). A node LP that is infeasible or exceeds the cutoff yields a
  // proof: a combination of the rows (Farkas ray or LP duals) that no point within the node's
  // bounds satisfies. The proof holds everywhere, so it is stored and propagated at later nodes.
  struct Conflict {
    std::vector<Index> index;
    std::vector<double> value;
    double lower = -kInf;  // lower <= a'x <= upper, where the upper side of a cutoff proof is
    double upper = kInf;   // upper + the current cutoff (it tightens as the incumbent improves).
    bool cutoff_relative = false;
    long long last_used = 0;  // Node count when it last pruned or tightened.
  };
  void analyze_infeasible_lp();
  void analyze_cutoff_lp();
  // Keeps a proof  lower <= g'x <= upper  (g over the columns, dense) if it is sparse enough;
  // returns whether it was kept.
  bool add_conflict(const std::vector<double>& g, double lower, double upper, bool relative);
  bool propagate_conflicts(std::vector<BoundChange>* record);

  // Search.
  NodeResult process_node(Node& node, std::vector<Node>& children);
  bool propagate(std::vector<Index> changed, std::vector<BoundChange>* record);
  void reduced_cost_fixing(Node& node, double objective, const std::vector<double>& reduced,
                           const std::vector<VarStatus>& basis);
  Index select_branching(const std::vector<Index>& fractional, const std::vector<double>& x,
                         double objective, const std::vector<VarStatus>& basis,
                         bool& node_infeasible, BoundChange& tighten, bool& has_tighten);
  bool strong_budget_left() const;
  std::size_t plunge_child(const std::vector<Node>& children) const;
  double pseudocost(Index j, bool up) const;
  void record_pseudocost(Index j, bool up, double gain, double distance);
  double cutoff() const;
  double effective_bound(double objective) const;
  double best_bound(double extra) const;
  double remaining_time() const;  // Seconds left, minus the time to release the stored nodes.
  bool time_up() const;

  // Solutions.
  bool is_fractional(double v) const;
  bool improves(double value) const;  // Better than the incumbent (minimization).
  bool try_solution(std::vector<double> x);
  void simple_rounding(const std::vector<double>& x);
  // Fixes the integers at their rounded values and re-solves the continuous columns; true if
  // that gave a new incumbent.
  bool round_and_solve(const std::vector<double>& x, const std::vector<VarStatus>& basis,
                       long long iteration_limit);
  // Uses options_.start (see MipOptions::start); returns true if it gave an incumbent.
  bool use_start();

  Model model_;  // Own copy: cuts are appended as rows.
  MipOptions options_;
  const Logger& log_;
  Timer timer_;
  Index m_;
  Index n_;
  Index original_rows_;
  double sense_;
  std::vector<double> cost_;  // Minimization costs.
  double offset_;
  bool integral_objective_ = false;
  std::vector<Index> integers_;
  SparseMatrix At_;
  std::vector<int> down_locks_;
  std::vector<int> up_locks_;

  Logger lp_log_{0};
  Scaling scaling_;
  LpProblem lp_;
  std::unique_ptr<Simplex> simplex_;
  std::vector<double> x_;  // Relaxation solution in the original space.
  SimplexStatus lp_status_ = SimplexStatus::kNumericalError;

  std::vector<double> root_lower_;
  std::vector<double> root_upper_;
  std::vector<double> lower_;
  std::vector<double> upper_;
  std::vector<Index> touched_;
  std::vector<char> is_touched_;
  std::vector<VarStatus> root_basis_;
  // Globally valid cuts not in the LP: root cuts that ended up non-binding. Tree nodes add the ones
  // their LP point violates.
  std::vector<Cut> cut_pool_;
  std::vector<Conflict> conflicts_;
  double tree_cut_seconds_ = 0.0;
  // Tree cut rounds so far, and those whose cuts stayed in the LP.
  long long tree_cut_rounds_ = 0;
  long long tree_cut_kept_rounds_ = 0;
  // Rounds whose cuts were removed back off: the next round waits this many nodes, doubling after
  // each failure, back to 1 after a success. Fresh separation that finds nothing backs off the
  // same way on its own interval.
  long long next_tree_cut_node_ = 0;
  long long tree_cut_interval_ = 1;
  long long next_fresh_node_ = 0;
  long long fresh_interval_ = 1;

  std::vector<char> row_mark_;  // Propagation queue membership.

  std::vector<double> pc_sum_[2];  // [0] down, [1] up: summed gain per unit change.
  std::vector<int> pc_count_[2];
  double pc_total_sum_[2] = {0.0, 0.0};
  long long pc_total_count_[2] = {0, 0};

  std::multimap<double, Node> open_;
  std::vector<Node> dive_stack_;  // Depth-first mode once open_ reaches the soft limit.
  double incumbent_value_ = kInf;
  std::vector<double> incumbent_;
  double pruned_bound_ = kInf;  // Smallest bound among nodes pruned by the cutoff.
  MipOutcome outcome_;
  bool incomplete_ = false;
  bool open_limit_logged_ = false;

  // Heuristic state. While logging_undo_ is set, set_bound records the previous bounds so a dive
  // can backtrack and restore the node.
  std::vector<BoundChange> undo_log_;
  bool logging_undo_ = false;
  int next_dive_rule_ = 0;
  long long next_dive_node_ = 0;
  long long dive_interval_ = 10;
  long long next_rins_node_ = 0;
  double rins_incumbent_ = kInf;
  // The incumbent value Feasibility Jump found (kInf if none): its point ignores the objective,
  // so the root still runs the objective feasibility pump while it is the incumbent.
  double jump_incumbent_value_ = kInf;
  double sub_mip_seconds_ = 0.0;
  std::mt19937 rng_{12345};

  // Parallel search: the shared pool (null when sequential), this worker's index and the pool
  // size last seen (for the time reserve).
  Shared* shared_ = nullptr;
  int worker_id_ = 0;
  std::size_t shared_stored_ = 0;
};

}  // namespace samaya
