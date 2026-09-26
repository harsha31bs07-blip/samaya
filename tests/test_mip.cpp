// Branch-and-bound against the reference MILP solver (depth-first branch-and-bound on the dense
// reference simplex) on random bounded MILPs, plus limits and known models.

#include <chrono>
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

#include "core/log.hpp"
#include "mip/branch_and_bound.hpp"
#include "mip/cuts.hpp"
#include "mip/feasibility_jump.hpp"
#include "reference_milp.hpp"
#include "samaya/io.hpp"
#include "samaya/solver.hpp"
#include "samaya/verify.hpp"
#include "test_framework.hpp"

using samaya::Index;
using samaya::kInf;
using samaya::Model;
using samaya::Status;
using samaya::test::ReferenceMilpResult;

namespace {

enum class MilpFamily { kMixed, kPureInteger, kKnapsack, kEquality };

// Bounded random MILP. Rows are built around a point that is integral in the integer columns, so
// most models are feasible; kEquality adds equality rows over integer columns (often infeasible).
Model random_milp(MilpFamily family, std::mt19937& rng) {
  const auto uniform_int = [&](int lo, int hi) {
    return std::uniform_int_distribution<int>(lo, hi)(rng);
  };
  const auto chance = [&](double p) { return std::bernoulli_distribution(p)(rng); };
  const auto real = [&](double lo, double hi) {
    return std::uniform_real_distribution<double>(lo, hi)(rng);
  };
  Model model;
  model.sense = chance(0.4) ? samaya::ObjSense::kMaximize : samaya::ObjSense::kMinimize;
  const int n = family == MilpFamily::kKnapsack ? uniform_int(4, 14) : uniform_int(2, 10);
  std::vector<double> point(static_cast<std::size_t>(n));
  for (int j = 0; j < n; ++j) {
    const bool integer = family != MilpFamily::kMixed || chance(0.6);
    const int lo = family == MilpFamily::kKnapsack ? 0 : uniform_int(-3, 2);
    const int span = family == MilpFamily::kKnapsack ? uniform_int(1, 3) : uniform_int(1, 6);
    model.col_lower.push_back(lo);
    model.col_upper.push_back(lo + span);
    model.col_type.push_back(integer ? samaya::VarType::kInteger : samaya::VarType::kContinuous);
    point[j] = integer ? uniform_int(lo, lo + span) : real(lo, lo + span);
    model.obj.push_back(family == MilpFamily::kKnapsack ? uniform_int(1, 20)
                                                        : uniform_int(-9, 9) * real(0.5, 1.5));
  }
  if (family == MilpFamily::kKnapsack) model.sense = samaya::ObjSense::kMaximize;

  std::vector<samaya::Triplet> t;
  const int m = family == MilpFamily::kKnapsack ? uniform_int(1, 3) : uniform_int(1, 8);
  for (int i = 0; i < m; ++i) {
    double act = 0.0;
    double weight_sum = 0.0;
    for (int j = 0; j < n; ++j) {
      if (family != MilpFamily::kKnapsack && !chance(0.6)) continue;
      double a = family == MilpFamily::kKnapsack ? uniform_int(1, 25) : uniform_int(-6, 6);
      if (a == 0.0) a = 3.0;
      if (family == MilpFamily::kMixed && chance(0.3)) a *= real(0.5, 1.5);
      t.push_back({i, j, a});
      act += a * point[j];
      weight_sum += a * model.col_upper[j];
    }
    double lo = -kInf;
    double up = kInf;
    if (family == MilpFamily::kKnapsack) {
      up = std::floor(weight_sum * real(0.3, 0.7));
    } else if (family == MilpFamily::kEquality && chance(0.5)) {
      // Shifted right-hand sides may leave no integer solution even when the LP is feasible.
      const double r = real(0, 1);
      lo = up = act + (r < 0.3 ? uniform_int(-2, 2) : r < 0.5 ? real(0.1, 0.9) : 0.0);
    } else {
      const double slack = chance(0.5) ? 0.0 : real(0, 4);
      const double r = real(0, 1);
      if (r < 0.4) {
        up = act + slack;
      } else if (r < 0.8) {
        lo = act - slack;
      } else {
        lo = act - slack;
        up = act + slack + 1;
      }
      // Fractional right-hand sides make the LP bound weaker than the integer optimum.
      if (chance(0.5) && up < kInf) up += real(0, 0.9);
    }
    model.row_lower.push_back(lo);
    model.row_upper.push_back(up);
  }
  model.A = samaya::SparseMatrix::from_triplets(m, n, std::move(t));
  return model;
}

const char* name(MilpFamily f) {
  switch (f) {
    case MilpFamily::kMixed: return "mixed";
    case MilpFamily::kPureInteger: return "pure";
    case MilpFamily::kKnapsack: return "knapsack";
    case MilpFamily::kEquality: return "equality";
  }
  return "?";
}

struct MilpTally {
  int optimal = 0;
  int infeasible = 0;
  int failures = 0;
  long long nodes = 0;
};

MilpTally cross_check(MilpFamily family, int count, unsigned seed, bool presolve) {
  std::mt19937 rng(seed);
  MilpTally tally;
  samaya::Params params;
  params.log_level = 0;
  params.presolve = presolve;
  params.mip_rel_gap = 0.0;
  params.mip_abs_gap = 1e-9;
  for (int k = 0; k < count; ++k) {
    const Model model = random_milp(family, rng);
    const ReferenceMilpResult ref = samaya::test::reference_milp(model);
    if (ref.status == ReferenceMilpResult::Status::kNodeLimit) continue;
    const samaya::Result got = samaya::Solver(params).solve(model);
    bool ok = false;
    if (ref.status == ReferenceMilpResult::Status::kOptimal) {
      ++tally.optimal;
      ok = got.status == Status::kOptimal && got.verified &&
           std::fabs(got.objective - ref.objective) <= 1e-6 * (1.0 + std::fabs(ref.objective)) &&
           std::fabs(got.dual_bound - got.objective) <= 1e-6 * (1.0 + std::fabs(ref.objective));
    } else if (ref.status == ReferenceMilpResult::Status::kInfeasible) {
      ++tally.infeasible;
      ok = got.status == Status::kInfeasible;
    }
    tally.nodes += got.nodes;
    if (!ok) {
      ++tally.failures;
      std::fprintf(stderr, "  %s #%d (%dx%d): reference %d %.9g, solver %s %.9g bound %.9g\n",
                   name(family), k, model.num_rows(), model.num_cols(),
                   static_cast<int>(ref.status), ref.objective, samaya::to_string(got.status),
                   got.objective, got.dual_bound);
    }
  }
  std::printf("  %-9s presolve %d: %d optimal, %d infeasible, %d failures, %lld nodes\n",
              name(family), presolve, tally.optimal, tally.infeasible, tally.failures,
              tally.nodes);
  return tally;
}

// Set partitioning with a planted partition: every row is covered exactly once. Rounding the LP
// point rarely gives a partition, so these need the pump, diving or a sub-MIP.
Model set_partitioning(std::mt19937& rng) {
  const auto uniform_int = [&](int lo, int hi) {
    return std::uniform_int_distribution<int>(lo, hi)(rng);
  };
  const int m = uniform_int(15, 30);
  std::vector<std::vector<int>> columns;
  // The planted partition: consecutive blocks of 1-3 rows.
  for (int i = 0; i < m;) {
    const int size = std::min(m - i, uniform_int(1, 3));
    std::vector<int> col;
    for (int k = 0; k < size; ++k) col.push_back(i + k);
    columns.push_back(col);
    i += size;
  }
  const int extra = uniform_int(2 * m, 4 * m);
  for (int k = 0; k < extra; ++k) {
    std::vector<int> col;
    const int size = uniform_int(2, 4);
    while (static_cast<int>(col.size()) < size) {
      const int r = uniform_int(0, m - 1);
      if (std::find(col.begin(), col.end(), r) == col.end()) col.push_back(r);
    }
    columns.push_back(col);
  }
  std::shuffle(columns.begin(), columns.end(), rng);
  Model model;
  std::vector<samaya::Triplet> t;
  for (std::size_t j = 0; j < columns.size(); ++j) {
    // Costs slightly below the size favour the random (non-planted) columns in the LP.
    model.obj.push_back(static_cast<double>(columns[j].size()) *
                        (0.8 + 0.4 * (uniform_int(0, 99) / 100.0)));
    model.col_lower.push_back(0);
    model.col_upper.push_back(1);
    model.col_type.push_back(samaya::VarType::kInteger);
    for (const int r : columns[j]) t.push_back({r, static_cast<Index>(j), 1.0});
  }
  model.row_lower.assign(static_cast<std::size_t>(m), 1.0);
  model.row_upper.assign(static_cast<std::size_t>(m), 1.0);
  model.A = samaya::SparseMatrix::from_triplets(m, static_cast<Index>(columns.size()),
                                                std::move(t));
  return model;
}

// Fixed-charge facility location: open y_i (binary, fixed cost), ship x_ij <= d_j y_i
// (variable upper bounds), meet each demand, and respect capacities sum_j x_ij <= cap_i y_i.
Model fixed_charge(std::mt19937& rng) {
  const auto uniform_int = [&](int lo, int hi) {
    return std::uniform_int_distribution<int>(lo, hi)(rng);
  };
  const int facilities = uniform_int(2, 4);
  const int customers = uniform_int(2, 4);
  Model model;
  std::vector<samaya::Triplet> t;
  std::vector<double> demand;
  double total = 0.0;
  for (int j = 0; j < customers; ++j) {
    demand.push_back(uniform_int(3, 12) + 0.5 * uniform_int(0, 1));
    total += demand.back();
  }
  const auto y = [&](int i) { return static_cast<Index>(i); };
  const auto x = [&](int i, int j) { return static_cast<Index>(facilities + i * customers + j); };
  for (int i = 0; i < facilities; ++i) {
    model.obj.push_back(uniform_int(10, 40));
    model.col_lower.push_back(0);
    model.col_upper.push_back(1);
    model.col_type.push_back(samaya::VarType::kInteger);
  }
  for (int i = 0; i < facilities; ++i) {
    for (int j = 0; j < customers; ++j) {
      model.obj.push_back(uniform_int(1, 9));
      model.col_lower.push_back(0);
      model.col_upper.push_back(kInf);
      model.col_type.push_back(samaya::VarType::kContinuous);
    }
  }
  Index row = 0;
  for (int j = 0; j < customers; ++j, ++row) {
    for (int i = 0; i < facilities; ++i) t.push_back({row, x(i, j), 1.0});
    model.row_lower.push_back(demand[j]);
    model.row_upper.push_back(demand[j]);
  }
  for (int i = 0; i < facilities; ++i, ++row) {
    const double cap = std::floor(total * (0.4 + 0.1 * uniform_int(0, 6)));
    for (int j = 0; j < customers; ++j) t.push_back({row, x(i, j), 1.0});
    t.push_back({row, y(i), -cap});
    model.row_lower.push_back(-kInf);
    model.row_upper.push_back(0.0);
  }
  for (int i = 0; i < facilities; ++i) {
    for (int j = 0; j < customers; ++j, ++row) {
      t.push_back({row, x(i, j), 1.0});
      t.push_back({row, y(i), -demand[j]});
      model.row_lower.push_back(-kInf);
      model.row_upper.push_back(0.0);
    }
  }
  model.A = samaya::SparseMatrix::from_triplets(row, static_cast<Index>(model.obj.size()),
                                                std::move(t));
  return model;
}

// Multi-dimensional knapsack over 16-22 binaries and 2-4 rows with correlated weights: small
// enough for the reference search, deep enough that the search branches past the root.
Model multi_knapsack(std::mt19937& rng) {
  const auto uniform_int = [&](int lo, int hi) {
    return std::uniform_int_distribution<int>(lo, hi)(rng);
  };
  const int n = uniform_int(16, 22);
  const int m = uniform_int(2, 4);
  Model model;
  model.sense = samaya::ObjSense::kMaximize;
  std::vector<samaya::Triplet> t;
  std::vector<double> total(static_cast<std::size_t>(m), 0.0);
  for (int j = 0; j < n; ++j) {
    double value = 0.0;
    for (int i = 0; i < m; ++i) {
      const int w = uniform_int(10, 60);
      t.push_back({i, j, static_cast<double>(w)});
      total[static_cast<std::size_t>(i)] += w;
      value += w;
    }
    model.obj.push_back(value / m + uniform_int(-5, 5));
    model.col_lower.push_back(0);
    model.col_upper.push_back(1);
    model.col_type.push_back(samaya::VarType::kInteger);
  }
  for (int i = 0; i < m; ++i) {
    model.row_lower.push_back(-kInf);
    model.row_upper.push_back(std::floor(0.5 * total[static_cast<std::size_t>(i)]));
  }
  model.A = samaya::SparseMatrix::from_triplets(m, n, std::move(t));
  return model;
}

// Transshipment network with big-M arcs in both directions between linked nodes, as in
// p200x1188c and mc11 (every link there is an arc pair). A demand node then also has outflows,
// so the flow cover on the <= side of its balance row, x_a <= d y_a + (M - d) (outflows' y), is
// weak; the >= side gives the cut-set inequality: some inflow arc must be open. Node 0 supplies
// everything; `arcs` returns each arc's (tail, head), flows first, then their binaries.
Model transshipment_network(std::mt19937& rng, std::vector<std::pair<int, int>>& arcs) {
  // A big M far above every demand (d / M < 0.05, as on mc11: 1 / 212): c-MIR rejects such a
  // fraction, so the flow covers must produce the cut-set inequalities.
  constexpr double kBigM = 1000.0;
  const auto uniform_int = [&](int lo, int hi) {
    return std::uniform_int_distribution<int>(lo, hi)(rng);
  };
  const int nodes = uniform_int(4, 6);
  std::vector<double> demand(static_cast<std::size_t>(nodes), 0.0);
  double supply = 0.0;
  for (int v = 1; v < nodes; ++v) {
    demand[v] = uniform_int(1, 6);
    supply += demand[v];
  }
  demand[0] = -supply;
  arcs.clear();
  const auto link = [&](int a, int b) {
    arcs.push_back({a, b});
    arcs.push_back({b, a});
  };
  for (int v = 1; v < nodes; ++v) link(uniform_int(0, v - 1), v);  // Connected.
  const int extra = uniform_int(1, 2);
  for (int k = 0; k < extra; ++k) {
    const int a = uniform_int(0, nodes - 1);
    const int b = uniform_int(0, nodes - 1);
    if (a != b) link(a, b);
  }
  Model model;
  std::vector<samaya::Triplet> t;
  const auto na = static_cast<Index>(arcs.size());
  for (Index a = 0; a < na; ++a) {
    model.obj.push_back(uniform_int(1, 5));
    model.col_lower.push_back(0);
    model.col_upper.push_back(kInf);
    model.col_type.push_back(samaya::VarType::kContinuous);
  }
  for (Index a = 0; a < na; ++a) {
    model.obj.push_back(uniform_int(10, 60));
    model.col_lower.push_back(0);
    model.col_upper.push_back(1);
    model.col_type.push_back(samaya::VarType::kInteger);
  }
  Index row = 0;
  for (int v = 0; v < nodes; ++v, ++row) {  // inflow - outflow = demand.
    for (Index a = 0; a < na; ++a) {
      if (arcs[a].second == v) t.push_back({row, a, 1.0});
      if (arcs[a].first == v) t.push_back({row, a, -1.0});
    }
    model.row_lower.push_back(demand[v]);
    model.row_upper.push_back(demand[v]);
  }
  for (Index a = 0; a < na; ++a, ++row) {
    t.push_back({row, a, 1.0});
    t.push_back({row, na + a, -kBigM});
    model.row_lower.push_back(-kInf);
    model.row_upper.push_back(0.0);
  }
  model.A = samaya::SparseMatrix::from_triplets(row, 2 * na, std::move(t));
  return model;
}

// Fixed-charge network flow with big-M arcs (x_a <= M y_a, M = total supply), as in p200x1188c
// and mc11: flow conservation rows over continuous flows only, so c-MIR finds nothing and the
// root gap needs flow covers.
Model big_m_network(std::mt19937& rng) {
  const auto uniform_int = [&](int lo, int hi) {
    return std::uniform_int_distribution<int>(lo, hi)(rng);
  };
  const int nodes = uniform_int(3, 5);
  std::vector<double> demand(static_cast<std::size_t>(nodes), 0.0);
  double supply = 0.0;
  for (int v = 1; v < nodes; ++v) {
    demand[v] = uniform_int(1, 9);
    supply += demand[v];
  }
  demand[0] = -supply;  // Node 0 supplies everything.
  std::vector<std::pair<int, int>> arcs;
  for (int v = 1; v < nodes; ++v) arcs.push_back({uniform_int(0, v - 1), v});  // Connected.
  const int extra = uniform_int(2, 5);
  for (int k = 0; k < extra; ++k) {
    const int a = uniform_int(0, nodes - 1);
    const int b = uniform_int(0, nodes - 1);
    if (a != b) arcs.push_back({a, b});
  }
  Model model;
  std::vector<samaya::Triplet> t;
  const auto na = static_cast<Index>(arcs.size());
  for (Index a = 0; a < na; ++a) {  // Flows, then their binaries.
    model.obj.push_back(uniform_int(1, 5));
    model.col_lower.push_back(0);
    model.col_upper.push_back(kInf);
    model.col_type.push_back(samaya::VarType::kContinuous);
  }
  for (Index a = 0; a < na; ++a) {
    model.obj.push_back(uniform_int(10, 60));
    model.col_lower.push_back(0);
    model.col_upper.push_back(1);
    model.col_type.push_back(samaya::VarType::kInteger);
  }
  Index row = 0;
  for (int v = 0; v < nodes; ++v, ++row) {  // inflow - outflow = demand.
    for (Index a = 0; a < na; ++a) {
      if (arcs[a].second == v) t.push_back({row, a, 1.0});
      if (arcs[a].first == v) t.push_back({row, a, -1.0});
    }
    model.row_lower.push_back(demand[v]);
    model.row_upper.push_back(demand[v]);
  }
  for (Index a = 0; a < na; ++a, ++row) {
    t.push_back({row, a, 1.0});
    t.push_back({row, na + a, -supply});
    model.row_lower.push_back(-kInf);
    model.row_upper.push_back(0.0);
  }
  model.A = samaya::SparseMatrix::from_triplets(row, 2 * na, std::move(t));
  return model;
}

// Each row times a factor between 1e-3 and 1e3: the same model, but its LP is scaled, so anything
// derived from the scaled LP (a proof, say) must map the scaled multipliers back to the rows.
void scale_rows_randomly(Model& model, std::mt19937& rng) {
  std::vector<double> factor(static_cast<std::size_t>(model.num_rows()));
  for (double& f : factor) f = std::pow(10.0, std::uniform_real_distribution<double>(-3, 3)(rng));
  std::vector<samaya::Triplet> t;
  const auto start = model.A.col_start();
  const auto index = model.A.row_index();
  const auto value = model.A.values();
  for (Index j = 0; j < model.num_cols(); ++j) {
    for (samaya::NnzIndex p = start[j]; p < start[j + 1]; ++p) {
      t.push_back({index[p], j, value[p] * factor[index[p]]});
    }
  }
  for (Index i = 0; i < model.num_rows(); ++i) {
    model.row_lower[i] *= factor[i];
    model.row_upper[i] *= factor[i];
  }
  model.A = samaya::SparseMatrix::from_triplets(model.num_rows(), model.num_cols(), std::move(t));
}

// Capacitated lot sizing in batches over 5-8 periods: production x_p (continuous) comes in
// batches of `size` units (x_p <= size z_p, z_p integer in [0, 3]), and the stock s_p carried out
// of period p is capped below one batch:
//   s_{p-1} + x_p - s_p = d_p  (no stock before the first period),   x_p - size z_p <= 0.
// Propagation tightens only integer columns and checks one row at a time, and the balance rows
// hold no integer column, so a node with too few batches for a demand is found infeasible only by
// its LP, through a balance row and a linking row together. Every column is bounded and the data
// are integers, so each such LP misses its demand by at least one unit. Demands of at most three
// batches keep every model feasible (make d_p in period p), and with 4^8 integer points the
// reference search has at most 2 * 4^8 nodes, below its limit of 200000.
Model lot_sizing(std::mt19937& rng) {
  const auto uniform_int = [&](int lo, int hi) {
    return std::uniform_int_distribution<int>(lo, hi)(rng);
  };
  constexpr int kBatches = 3;
  const int periods = uniform_int(5, 8);
  const int size = uniform_int(4, 12);
  // Below one batch: a missing batch can be made up from stock only in part.
  const int stock = uniform_int(1, size - 1);
  const auto x = [&](int p) { return static_cast<Index>(p); };
  const auto s = [&](int p) { return static_cast<Index>(periods + p); };
  const auto z = [&](int p) { return static_cast<Index>(2 * periods + p); };
  Model model;
  for (int p = 0; p < periods; ++p) {  // Production: a cost per unit.
    model.obj.push_back(uniform_int(1, 5));
    model.col_lower.push_back(0);
    model.col_upper.push_back(kBatches * size);
    model.col_type.push_back(samaya::VarType::kContinuous);
  }
  for (int p = 0; p < periods; ++p) {  // Stock: a holding cost per unit.
    model.obj.push_back(uniform_int(1, 3));
    model.col_lower.push_back(0);
    model.col_upper.push_back(stock);
    model.col_type.push_back(samaya::VarType::kContinuous);
  }
  for (int p = 0; p < periods; ++p) {  // Batches: a fixed cost each.
    model.obj.push_back(uniform_int(10, 40));
    model.col_lower.push_back(0);
    model.col_upper.push_back(kBatches);
    model.col_type.push_back(samaya::VarType::kInteger);
  }
  std::vector<samaya::Triplet> t;
  Index row = 0;
  for (int p = 0; p < periods; ++p, ++row) {  // s_{p-1} + x_p - s_p = d_p.
    if (p > 0) t.push_back({row, s(p - 1), 1.0});
    t.push_back({row, x(p), 1.0});
    t.push_back({row, s(p), -1.0});
    const double demand = uniform_int(size, kBatches * size);
    model.row_lower.push_back(demand);
    model.row_upper.push_back(demand);
  }
  for (int p = 0; p < periods; ++p, ++row) {  // x_p - size z_p <= 0.
    t.push_back({row, x(p), 1.0});
    t.push_back({row, z(p), -static_cast<double>(size)});
    model.row_lower.push_back(-kInf);
    model.row_upper.push_back(0.0);
  }
  model.A = samaya::SparseMatrix::from_triplets(row, static_cast<Index>(3 * periods), std::move(t));
  return model;
}

}  // namespace

TEST(mip_heuristics_find_solutions_at_the_root) {
  // The root alone (one node, no cuts): simple rounding and round-and-solve against the pump,
  // diving and RENS. Every solution must pass the verifier on the model.
  const samaya::Logger quiet(0);
  int found[2] = {0, 0};
  int models = 0;
  std::mt19937 rng(31);
  for (int k = 0; k < 40; ++k) {
    const Model model = set_partitioning(rng);
    ++models;
    for (int h = 0; h < 2; ++h) {
      samaya::MipOptions options;
      options.cuts = false;
      options.node_limit = 1;
      options.heuristics = h == 1;
      const samaya::MipOutcome out = samaya::BranchAndBound(model, options, quiet).solve();
      if (out.x.empty()) continue;
      ++found[h];
      const samaya::VerifyReport report = samaya::verify_primal(model, out.x, {});
      CHECK(report.ok);
      CHECK(out.bound <= out.objective + 1e-9);
    }
  }
  std::printf("  %d set-partitioning models: solution at the root %d without heuristics, %d "
              "with\n", models, found[0], found[1]);
  CHECK(found[1] >= 36);
  CHECK(found[1] > found[0]);
}

TEST(mip_objective_cutoff_accepts_only_better_solutions) {
  const samaya::Logger quiet(0);
  std::mt19937 rng(41);
  int checked = 0;
  for (int k = 0; k < 120; ++k) {
    const Model model = random_milp(k % 2 ? MilpFamily::kMixed : MilpFamily::kKnapsack, rng);
    const ReferenceMilpResult ref = samaya::test::reference_milp(model);
    if (ref.status != ReferenceMilpResult::Status::kOptimal) continue;
    ++checked;
    const double worse = model.sense == samaya::ObjSense::kMaximize ? -1.0 : 1.0;
    samaya::MipOptions options;
    options.rel_gap = 0.0;
    options.abs_gap = 1e-9;
    // Nothing is better than the optimum (the margin covers round-off between the reference's
    // objective and ours).
    options.objective_cutoff = ref.objective - worse * 1e-6 * (1.0 + std::fabs(ref.objective));
    CHECK(samaya::BranchAndBound(model, options, quiet).solve().x.empty());
    // A looser cutoff still finds the optimum.
    options.objective_cutoff = ref.objective + worse * (1.0 + std::fabs(ref.objective));
    const samaya::MipOutcome out = samaya::BranchAndBound(model, options, quiet).solve();
    CHECK(out.status == Status::kOptimal);
    CHECK(std::fabs(out.objective - ref.objective) <= 1e-6 * (1.0 + std::fabs(ref.objective)));
  }
  CHECK(checked > 60);
}

TEST(mip_random_models_match_reference) {
  for (const bool presolve : {true, false}) {
    const MilpTally mixed = cross_check(MilpFamily::kMixed, 250, 11, presolve);
    const MilpTally pure = cross_check(MilpFamily::kPureInteger, 250, 12, presolve);
    const MilpTally knapsack = cross_check(MilpFamily::kKnapsack, 200, 13, presolve);
    const MilpTally equality = cross_check(MilpFamily::kEquality, 250, 14, presolve);
    CHECK_EQ(mixed.failures + pure.failures + knapsack.failures + equality.failures, 0);
    CHECK(mixed.optimal > 150);
    CHECK(pure.optimal > 150);
    CHECK(knapsack.optimal > 150);
    CHECK(equality.infeasible > 20);
    CHECK(knapsack.nodes > 200);  // The search actually branches.
  }
}

TEST(mip_cuts_never_separate_an_optimal_solution) {
  // Root cuts are checked against a known optimal solution (the reference's); a valid cut can
  // never separate it. The search runs without presolve so the solution applies directly.
  const samaya::Logger quiet(0);
  long long violations = 0;
  long long cuts = 0;
  int models = 0;
  int tightened = 0;
  for (const MilpFamily family :
       {MilpFamily::kMixed, MilpFamily::kPureInteger, MilpFamily::kKnapsack}) {
    std::mt19937 rng(100 + static_cast<unsigned>(family));
    for (int k = 0; k < 200; ++k) {
      const Model model = random_milp(family, rng);
      const ReferenceMilpResult ref = samaya::test::reference_milp(model);
      if (ref.status != ReferenceMilpResult::Status::kOptimal) continue;
      samaya::MipOptions options;
      options.rel_gap = 0.0;
      options.abs_gap = 1e-9;
      options.debug_solution = ref.x;
      samaya::BranchAndBound search(model, options, quiet);
      const samaya::MipOutcome out = search.solve();
      ++models;
      violations += out.debug_cut_violations;
      cuts += out.cut_rounds > 0 ? out.cuts_added : 0;
      const double sense = model.sense == samaya::ObjSense::kMaximize ? -1.0 : 1.0;
      if (sense * out.root_bound_cuts > sense * out.root_bound + 1e-9) ++tightened;
      CHECK(out.status == Status::kOptimal);
      CHECK(std::fabs(out.objective - ref.objective) <= 1e-6 * (1 + std::fabs(ref.objective)));
    }
  }
  std::printf("  %d models, %lld cuts kept, %d root bounds tightened, %lld violations\n",
              models, cuts, tightened, violations);
  CHECK_EQ(violations, 0);
  CHECK(models > 400);
  CHECK(tightened > 100);
}

TEST(mip_root_reductions_keep_an_optimal_solution) {
  // Coefficient tightening and probing change the root bounds and rows; a known optimal
  // solution must stay feasible for them, and the search must still find the optimum.
  const samaya::Logger quiet(0);
  long long violations = 0;
  long long reductions = 0;
  int models = 0;
  for (const MilpFamily family :
       {MilpFamily::kMixed, MilpFamily::kPureInteger, MilpFamily::kKnapsack}) {
    std::mt19937 rng(200 + static_cast<unsigned>(family));
    for (int k = 0; k < 200; ++k) {
      Model model = random_milp(family, rng);
      // Binaries make the reductions apply more often.
      if (k % 2 == 0) {
        for (Index j = 0; j < model.num_cols(); ++j) {
          if (model.col_type[j] == samaya::VarType::kInteger) {
            model.col_lower[j] = 0;
            model.col_upper[j] = 1;
          }
        }
      }
      const ReferenceMilpResult ref = samaya::test::reference_milp(model);
      if (ref.status != ReferenceMilpResult::Status::kOptimal) continue;
      samaya::MipOptions options;
      options.rel_gap = 0.0;
      options.abs_gap = 1e-9;
      options.cuts = false;
      options.debug_solution = ref.x;
      const samaya::MipOutcome out = samaya::BranchAndBound(model, options, quiet).solve();
      ++models;
      violations += out.debug_reduction_violations;
      reductions += out.coefficients_tightened + out.probing_fixed + out.probing_tightened;
      CHECK(out.status == Status::kOptimal);
      CHECK(std::fabs(out.objective - ref.objective) <= 1e-6 * (1 + std::fabs(ref.objective)));
    }
  }
  // 5 x1 + x2 + x3 <= 6 over binaries is x1 + x2 + x3 <= 2 on integer points.
  Model knap;
  knap.sense = samaya::ObjSense::kMaximize;
  knap.obj = {5, 3, 3};  // LP bound 10 on the original row (x2 = x3 = 1, x1 = 0.8), 8 after.
  knap.col_lower = {0, 0, 0};
  knap.col_upper = {1, 1, 1};
  knap.col_type.assign(3, samaya::VarType::kInteger);
  knap.row_lower = {-kInf};
  knap.row_upper = {6};
  knap.A = samaya::SparseMatrix::from_triplets(1, 3, {{0, 0, 5}, {0, 1, 1}, {0, 2, 1}});
  samaya::MipOptions options;
  options.debug_solution = {1, 1, 0};
  const samaya::MipOutcome out = samaya::BranchAndBound(knap, options, quiet).solve();
  CHECK(out.coefficients_tightened >= 1);
  CHECK_EQ(out.debug_reduction_violations, 0);
  CHECK_NEAR(out.objective, 8.0, 1e-9);
  CHECK_NEAR(out.root_bound, 8.0, 1e-9);  // The tightened row makes the LP integral.

  std::printf("  %d models, %lld reductions, %lld violations\n", models, reductions, violations);
  CHECK_EQ(violations, 0);
  CHECK(models > 400);
  CHECK(reductions > 50);
}

TEST(mip_cuts_on_fixed_charge_models) {
  // Fixed-charge facility location (big-M rows x <= u y): no cut may separate a known optimal
  // solution, and the cuts must tighten the root.
  const samaya::Logger quiet(0);
  std::mt19937 rng(71);
  long long violations = 0;
  int models = 0;
  int tightened = 0;
  double closed = 0.0;
  for (int k = 0; k < 150; ++k) {
    const Model model = fixed_charge(rng);
    const ReferenceMilpResult ref = samaya::test::reference_milp(model);
    if (ref.status != ReferenceMilpResult::Status::kOptimal) continue;
    samaya::MipOptions options;
    options.rel_gap = 0.0;
    options.abs_gap = 1e-9;
    options.probing = false;
    options.debug_solution = ref.x;
    const samaya::MipOutcome out = samaya::BranchAndBound(model, options, quiet).solve();
    ++models;
    violations += out.debug_cut_violations;
    CHECK(out.status == Status::kOptimal);
    CHECK(std::fabs(out.objective - ref.objective) <= 1e-6 * (1 + std::fabs(ref.objective)));
    const double gap = ref.objective - out.root_bound;
    if (gap > 1e-6 && out.root_bound_cuts > out.root_bound + 1e-6) {
      ++tightened;
      closed += (out.root_bound_cuts - out.root_bound) / gap;
    }
  }
  std::printf("  %d fixed-charge models, %d root bounds tightened (%.0f%% of the gap closed on "
              "average), %lld violations\n", models, tightened,
              tightened > 0 ? 100.0 * closed / tightened : 0.0, violations);
  CHECK_EQ(violations, 0);
  CHECK(models > 100);
  CHECK(tightened > models / 2);
}

TEST(mip_tree_cuts_never_separate_an_optimal_solution) {
  // Tree cuts are separated with the root bounds so they hold at every node; one derived with a
  // node's own bounds could cut off the optimum elsewhere. Separating at every node (with room for
  // as many cuts as rows) on the random and fixed-charge families, no cut may separate a known
  // optimal solution, whether it stays in the LP or is removed again, and the optimum must match
  // the reference.
  const samaya::Logger quiet(0);
  long long violations = 0;
  long long separated = 0;
  long long kept = 0;
  int models = 0;
  int with_tree_cuts = 0;
  std::mt19937 rng(131);
  for (int k = 0; k < 200; ++k) {
    Model model;
    switch (k % 4) {
      case 0: model = random_milp(MilpFamily::kMixed, rng); break;
      case 1: model = multi_knapsack(rng); break;
      case 2: model = fixed_charge(rng); break;
      default: model = big_m_network(rng); break;
    }
    const ReferenceMilpResult ref = samaya::test::reference_milp(model);
    if (ref.status != ReferenceMilpResult::Status::kOptimal) continue;
    samaya::MipOptions options;
    options.rel_gap = 0.0;
    options.abs_gap = 1e-9;
    options.probing = false;
    options.restart = false;
    options.tree_separation_frequency = 1;
    options.tree_cut_row_fraction = 1.0;
    options.debug_solution = ref.x;
    const samaya::MipOutcome out = samaya::BranchAndBound(model, options, quiet).solve();
    ++models;
    violations += out.debug_cut_violations;
    separated += out.tree_cuts_separated;
    kept += out.tree_cuts;
    with_tree_cuts += out.tree_cuts_separated > 0;
    CHECK(out.status == Status::kOptimal);
    CHECK(std::fabs(out.objective - ref.objective) <= 1e-6 * (1 + std::fabs(ref.objective)));
  }
  std::printf("  %d models, %d with tree cuts (%lld separated, %lld kept), %lld violations\n",
              models, with_tree_cuts, separated, kept, violations);
  CHECK_EQ(violations, 0);
  CHECK(models > 150);
  CHECK(with_tree_cuts > 30);
}

TEST(mip_start_becomes_the_incumbent_or_is_completed) {
  // With heuristics off and no nodes, a solution can only come from the start. A feasible start
  // (the reference optimum) must be taken exactly as it is; with the continuous values removed
  // (NaN) the integer values must be kept and the continuous columns re-solved, which gives the
  // optimum again. A start outside the bounds must be ignored safely.
  const samaya::Logger quiet(0);
  std::mt19937 rng(151);
  int models = 0;
  int taken = 0;
  int completed = 0;
  for (int k = 0; k < 200; ++k) {
    const Model model = k % 2 == 0 ? random_milp(MilpFamily::kMixed, rng) : fixed_charge(rng);
    const ReferenceMilpResult ref = samaya::test::reference_milp(model);
    if (ref.status != ReferenceMilpResult::Status::kOptimal) continue;
    ++models;
    const double tol = 1e-6 * (1.0 + std::fabs(ref.objective));
    samaya::MipOptions options;
    options.heuristics = false;
    options.node_limit = 0;
    options.start = ref.x;
    samaya::MipOutcome out = samaya::BranchAndBound(model, options, quiet).solve();
    // Taken as it is (the plan does not move), not merely matched in objective by a re-solve.
    bool same = !out.x.empty();
    for (Index j = 0; same && j < model.num_cols(); ++j) {
      const bool integer = model.col_type[j] == samaya::VarType::kInteger;
      same = out.x[j] == (integer ? std::round(ref.x[j]) : ref.x[j]);
    }
    if (same) ++taken;

    for (Index j = 0; j < model.num_cols(); ++j) {
      if (model.col_type[j] == samaya::VarType::kContinuous) options.start[j] = std::nan("");
    }
    out = samaya::BranchAndBound(model, options, quiet).solve();
    // (Without continuous columns this is the same start as above.)
    if (!out.x.empty() && std::fabs(out.objective - ref.objective) <= tol) ++completed;

    options.start.assign(static_cast<std::size_t>(model.num_cols()), 1e9);
    options.node_limit = -1;
    options.heuristics = true;
    out = samaya::BranchAndBound(model, options, quiet).solve();
    CHECK(out.status == Status::kOptimal);
    CHECK(std::fabs(out.objective - ref.objective) <= tol);
  }
  std::printf("  %d models: start taken %d, completed from its integers %d\n", models, taken,
              completed);
  CHECK(models > 100);
  CHECK_EQ(taken, models);
  CHECK_EQ(completed, models);
}

TEST(mip_conflicts_never_exclude_an_optimal_solution) {
  // Proofs from infeasible node LPs (Farkas rays) hold at every node, and proofs from cut-off node
  // LPs (duals) hold for every improving solution. Checked against a known optimum on families
  // with many infeasible and cut-off nodes, the optimum must match the reference and no proof
  // may exclude the optimum while it would still improve on the incumbent.
  const samaya::Logger quiet(0);
  long long violations = 0;
  long long found = 0;
  long long prunes = 0;
  long long infeasible_lps = 0;
  long long farkas = 0;
  long long farkas_failures = 0;
  int models = 0;
  std::mt19937 rng(171);
  for (int k = 0; k < 120; ++k) {
    Model model;
    switch (k % 3) {
      case 0: model = multi_knapsack(rng); break;
      case 1: model = random_milp(MilpFamily::kEquality, rng); break;
      default: model = fixed_charge(rng); break;
    }
    const ReferenceMilpResult ref = samaya::test::reference_milp(model);
    if (ref.status != ReferenceMilpResult::Status::kOptimal) continue;
    scale_rows_randomly(model, rng);
    samaya::MipOptions options;
    options.rel_gap = 0.0;
    options.abs_gap = 1e-9;
    // Without root reductions, cuts and heuristics the trees are deep enough to prune by proofs.
    options.probing = false;
    options.restart = false;
    options.cuts = false;
    options.heuristics = false;
    options.conflict_density = 1.0;  // Keep every proof: validity does not depend on density.
    options.debug_solution = ref.x;
    const samaya::MipOutcome out = samaya::BranchAndBound(model, options, quiet).solve();
    ++models;
    violations += out.debug_conflict_violations;
    found += out.conflicts_found;
    prunes += out.conflict_prunes;
    infeasible_lps += out.conflict_infeasible_lps;
    farkas += out.farkas_proofs;
    farkas_failures += out.debug_farkas_failures;
    CHECK(out.status == Status::kOptimal);
    CHECK(std::fabs(out.objective - ref.objective) <= 1e-6 * (1 + std::fabs(ref.objective)));
  }
  std::printf("  %d models, %lld proofs stored (%lld from %lld infeasible LPs), %lld nodes pruned "
              "by them, %lld violations, %lld rays rejected\n",
              models, found, farkas, infeasible_lps, prunes, violations, farkas_failures);
  CHECK_EQ(violations, 0);
  CHECK_EQ(farkas_failures, 0);
  CHECK(models > 80);
  CHECK(found > 200);
  CHECK(prunes > 20);
}

TEST(mip_farkas_rays_prove_their_nodes_infeasible) {
  // Lot sizing with its rows scaled as above: an infeasible node there is found only by its LP, and
  // its Farkas ray combines a balance row with a linking row (see lot_sizing). Any multipliers give
  // a valid proof, so a wrong mapping of the ray to the original rows would only lose proofs; it
  // is caught by checking every ray with verify_infeasibility on its node's bounds. Strong
  // branching is off: it proves infeasible children before they get a node LP, and those are not
  // analyzed (the test above keeps it on).
  const samaya::Logger quiet(0);
  constexpr int kModels = 60;
  long long violations = 0;
  long long found = 0;
  long long prunes = 0;
  long long infeasible_lps = 0;
  long long farkas = 0;
  long long farkas_failures = 0;
  long long reference_nodes = 0;
  int models = 0;
  std::mt19937 rng(172);
  for (int k = 0; k < kModels; ++k) {
    Model model = lot_sizing(rng);
    const ReferenceMilpResult ref = samaya::test::reference_milp(model);
    REQUIRE(ref.status == ReferenceMilpResult::Status::kOptimal);
    reference_nodes += ref.nodes;
    scale_rows_randomly(model, rng);
    samaya::MipOptions options;
    options.rel_gap = 0.0;
    options.abs_gap = 1e-9;
    options.probing = false;
    options.restart = false;
    options.cuts = false;
    options.heuristics = false;
    options.max_strong_branching = 0;
    options.conflict_density = 1.0;
    options.debug_solution = ref.x;
    const samaya::MipOutcome out = samaya::BranchAndBound(model, options, quiet).solve();
    ++models;
    violations += out.debug_conflict_violations;
    found += out.conflicts_found;
    prunes += out.conflict_prunes;
    infeasible_lps += out.conflict_infeasible_lps;
    farkas += out.farkas_proofs;
    farkas_failures += out.debug_farkas_failures;
    CHECK(out.status == Status::kOptimal);
    CHECK(std::fabs(out.objective - ref.objective) <= 1e-6 * (1 + std::fabs(ref.objective)));
  }
  std::printf("  %d lot-sizing models (%lld reference nodes): %lld proofs stored (%lld from %lld "
              "infeasible LPs), %lld nodes pruned by them, %lld violations, %lld rays rejected\n",
              models, reference_nodes, found, farkas, infeasible_lps, prunes, violations,
              farkas_failures);
  CHECK_EQ(models, kModels);
  CHECK_EQ(violations, 0);
  CHECK_EQ(farkas_failures, 0);
  // Floors against a vacuous run, about a third of the counts at this seed (265 infeasible LPs,
  // all 265 rays kept as proofs, 63 prunes).
  CHECK(infeasible_lps > 80);
  CHECK(farkas > 80);
  CHECK(prunes > 20);
}

TEST(mip_conflicts_from_strong_branching_children_hold) {
  // Lot sizing again, now with strong branching: a child it proves infeasible or cut off yields a
  // proof too, and each Farkas ray is checked on that child's bounds (the strong-branched column
  // is still tightened when the ray is analyzed).
  const samaya::Logger quiet(0);
  constexpr int kModels = 60;
  long long violations = 0;
  long long found = 0;
  long long from_strong_branching = 0;
  long long prunes = 0;
  long long infeasible_lps = 0;
  long long farkas = 0;
  long long farkas_failures = 0;
  int models = 0;
  std::mt19937 rng(173);
  for (int k = 0; k < kModels; ++k) {
    Model model = lot_sizing(rng);
    const ReferenceMilpResult ref = samaya::test::reference_milp(model);
    REQUIRE(ref.status == ReferenceMilpResult::Status::kOptimal);
    scale_rows_randomly(model, rng);
    samaya::MipOptions options;
    options.rel_gap = 0.0;
    options.abs_gap = 1e-9;
    options.probing = false;
    options.restart = false;
    options.cuts = false;
    options.heuristics = false;
    options.conflict_density = 1.0;
    options.debug_solution = ref.x;
    const samaya::MipOutcome out = samaya::BranchAndBound(model, options, quiet).solve();
    ++models;
    violations += out.debug_conflict_violations;
    found += out.conflicts_found;
    from_strong_branching += out.strong_branching_conflicts;
    prunes += out.conflict_prunes;
    infeasible_lps += out.conflict_infeasible_lps;
    farkas += out.farkas_proofs;
    farkas_failures += out.debug_farkas_failures;
    CHECK(out.status == Status::kOptimal);
    CHECK(std::fabs(out.objective - ref.objective) <= 1e-6 * (1 + std::fabs(ref.objective)));
  }
  std::printf("  %d lot-sizing models: %lld proofs stored (%lld from strong branching, %lld from "
              "%lld infeasible LPs), %lld nodes pruned by them, %lld violations, %lld rays "
              "rejected\n",
              models, found, from_strong_branching, farkas, infeasible_lps, prunes, violations,
              farkas_failures);
  CHECK_EQ(models, kModels);
  CHECK_EQ(violations, 0);
  CHECK_EQ(farkas_failures, 0);
  // Floors against a vacuous run, about a third of the counts at this seed (387 proofs from
  // strong branching, 269 Farkas proofs).
  CHECK(from_strong_branching > 120);
  CHECK(farkas > 80);
}

TEST(mip_conflicts_with_an_integral_objective_hold) {
  // Multi-dimensional knapsacks with integer values: an improving solution's objective is then an
  // integer below the cutoff, and cut-off proofs are stated against that integer. They must
  // still never exclude an improving optimum.
  const samaya::Logger quiet(0);
  constexpr int kModels = 40;
  long long violations = 0;
  long long found = 0;
  long long cutoff_proofs = 0;
  long long prunes = 0;
  int models = 0;
  std::mt19937 rng(174);
  for (int k = 0; k < kModels; ++k) {
    Model model = multi_knapsack(rng);
    for (double& c : model.obj) c = std::round(c);
    const ReferenceMilpResult ref = samaya::test::reference_milp(model);
    REQUIRE(ref.status == ReferenceMilpResult::Status::kOptimal);
    scale_rows_randomly(model, rng);
    samaya::MipOptions options;
    options.rel_gap = 0.0;
    options.abs_gap = 1e-9;
    options.probing = false;
    options.restart = false;
    options.cuts = false;
    options.heuristics = false;
    options.conflict_density = 1.0;
    options.debug_solution = ref.x;
    const samaya::MipOutcome out = samaya::BranchAndBound(model, options, quiet).solve();
    ++models;
    violations += out.debug_conflict_violations;
    found += out.conflicts_found;
    cutoff_proofs += out.conflicts_found - out.farkas_proofs;
    prunes += out.conflict_prunes;
    CHECK(out.status == Status::kOptimal);
    CHECK(std::fabs(out.objective - ref.objective) <= 1e-6 * (1 + std::fabs(ref.objective)));
  }
  std::printf("  %d knapsacks with integer values: %lld proofs stored (%lld cut-off proofs), %lld "
              "nodes pruned by them, %lld violations\n",
              models, found, cutoff_proofs, prunes, violations);
  CHECK_EQ(models, kModels);
  CHECK_EQ(violations, 0);
  // Floors against a vacuous run, about a third of the counts at this seed (539 cut-off proofs,
  // 2062 prunes).
  CHECK(cutoff_proofs > 170);
  CHECK(prunes > 680);
}

TEST(mip_parallel_search_matches_reference) {
  // Four threads from the first node (no sequential ramp-up) on the random families: the
  // parallel search must reach the same optimum and status as the reference.
  const samaya::Logger quiet(0);
  int failures = 0;
  int models = 0;
  int parallel = 0;
  for (const MilpFamily family : {MilpFamily::kMixed, MilpFamily::kPureInteger,
                                  MilpFamily::kKnapsack, MilpFamily::kEquality}) {
    std::mt19937 rng(300 + static_cast<unsigned>(family));
    for (int k = 0; k < 120; ++k) {
      const Model model = random_milp(family, rng);
      const ReferenceMilpResult ref = samaya::test::reference_milp(model);
      if (ref.status == ReferenceMilpResult::Status::kNodeLimit) continue;
      samaya::MipOptions options;
      options.rel_gap = 0.0;
      options.abs_gap = 1e-9;
      options.threads = 4;
      options.parallel_start_nodes = 0;
      options.cuts = k % 2 == 0;
      const samaya::MipOutcome out = samaya::BranchAndBound(model, options, quiet).solve();
      ++models;
      parallel += out.threads_used > 1;
      bool ok = false;
      if (ref.status == ReferenceMilpResult::Status::kOptimal) {
        ok = out.status == Status::kOptimal &&
             std::fabs(out.objective - ref.objective) <= 1e-6 * (1 + std::fabs(ref.objective)) &&
             std::fabs(out.bound - out.objective) <= 1e-6 * (1 + std::fabs(ref.objective));
      } else {
        ok = out.status == Status::kInfeasible;
      }
      if (!ok) {
        ++failures;
        std::fprintf(stderr, "  %s #%d: reference %d %.9g, parallel %s %.9g bound %.9g\n",
                     name(family), k, static_cast<int>(ref.status), ref.objective,
                     samaya::to_string(out.status), out.objective, out.bound);
      }
    }
  }
  // Harder models, so that the threads actually share a pool: set partitioning with a node limit
  // high enough to finish.
  std::mt19937 rng(77);
  int agree = 0;
  int compared = 0;
  for (int k = 0; k < 6; ++k) {
    const Model model = set_partitioning(rng);
    samaya::MipOptions options;
    options.rel_gap = 0.0;
    options.abs_gap = 1e-9;
    options.parallel_start_nodes = 0;
    const samaya::MipOutcome one = samaya::BranchAndBound(model, options, quiet).solve();
    options.threads = 4;
    const samaya::MipOutcome four = samaya::BranchAndBound(model, options, quiet).solve();
    if (one.status != Status::kOptimal) continue;
    ++compared;
    agree += four.status == Status::kOptimal &&
             std::fabs(four.objective - one.objective) <= 1e-6 * (1 + std::fabs(one.objective));
  }
  std::printf("  %d random models (%d searched in parallel), %d failures; set partitioning %d/%d "
              "agree\n", models, parallel, failures, agree, compared);
  CHECK_EQ(failures, 0);
  CHECK(models > 400);
  CHECK(parallel > 50);
  CHECK_EQ(agree, compared);
}

TEST(mip_flow_covers_on_big_m_networks) {
  // Big-M fixed-charge networks: root cuts (flow covers on the conservation rows and their
  // aggregations) must never separate a known optimum and must close a real part of the gap.
  const samaya::Logger quiet(0);
  std::mt19937 rng(91);
  long long violations = 0;
  int models = 0;
  double closed = 0.0;
  int gaps = 0;
  for (int k = 0; k < 200; ++k) {
    const Model model = big_m_network(rng);
    const ReferenceMilpResult ref = samaya::test::reference_milp(model);
    if (ref.status != ReferenceMilpResult::Status::kOptimal) continue;
    samaya::MipOptions options;
    options.rel_gap = 0.0;
    options.abs_gap = 1e-9;
    options.probing = false;
    options.heuristics = false;
    options.debug_solution = ref.x;
    const samaya::MipOutcome out = samaya::BranchAndBound(model, options, quiet).solve();
    ++models;
    violations += out.debug_cut_violations;
    CHECK(out.status == Status::kOptimal);
    CHECK(std::fabs(out.objective - ref.objective) <= 1e-6 * (1 + std::fabs(ref.objective)));
    const double gap = ref.objective - out.root_bound;
    if (gap > 1e-6) {
      ++gaps;
      closed += (out.root_bound_cuts - out.root_bound) / gap;
    }
  }
  const double average = gaps > 0 ? closed / gaps : 0.0;
  std::printf("  %d big-M network models, %.0f%% of the root gap closed on average, %lld "
              "violations\n", models, 100.0 * average, violations);
  CHECK_EQ(violations, 0);
  CHECK(models > 150);
  CHECK(average > 0.3);
}

TEST(mip_root_cuts_are_valid_on_transshipment_networks) {
  // Root cuts on transshipment networks (every link an arc pair, as on p200x1188c and mc11): no
  // cut may separate the optimum. The root bound is reported against the LP with every
  // single-node cut-set inequality added (each demand node has an open inflow arc); on models
  // this small the other cuts reach it too, so the separator tests above check the mechanism.
  const samaya::Logger quiet(0);
  std::mt19937 rng(92);
  long long violations = 0;
  int models = 0;
  int gaps = 0;
  int reached = 0;
  double closed = 0.0;
  std::vector<std::pair<int, int>> arcs;
  for (int k = 0; k < 150; ++k) {
    const Model model = transshipment_network(rng, arcs);
    const ReferenceMilpResult ref = samaya::test::reference_milp(model);
    if (ref.status != ReferenceMilpResult::Status::kOptimal) continue;
    // The LP with the cut-set rows: sum of y_a over the arcs into v >= 1 for each demand node v.
    Model with_cut_sets = model;
    {
      const auto na = static_cast<Index>(arcs.size());
      std::vector<samaya::Triplet> t;
      const auto start = model.A.col_start();
      const auto index = model.A.row_index();
      const auto value = model.A.values();
      for (Index j = 0; j < model.num_cols(); ++j) {
        for (samaya::NnzIndex p = start[j]; p < start[j + 1]; ++p) {
          t.push_back({index[p], j, value[p]});
        }
      }
      Index row = model.num_rows();
      int max_node = 0;
      for (const auto& [a, b] : arcs) max_node = std::max({max_node, a, b});
      for (int v = 1; v <= max_node; ++v, ++row) {
        for (Index a = 0; a < na; ++a) {
          if (arcs[a].second == v) t.push_back({row, na + a, 1.0});
        }
        with_cut_sets.row_lower.push_back(1.0);
        with_cut_sets.row_upper.push_back(kInf);
      }
      with_cut_sets.A = samaya::SparseMatrix::from_triplets(row, model.num_cols(), std::move(t));
    }
    const samaya::test::ReferenceResult cut_set_lp =
        samaya::test::ReferenceLp(with_cut_sets).solve();
    REQUIRE(cut_set_lp.status == samaya::test::ReferenceResult::Status::kOptimal);
    samaya::MipOptions options;
    options.rel_gap = 0.0;
    options.abs_gap = 1e-9;
    options.probing = false;
    options.heuristics = false;
    options.debug_solution = ref.x;
    const samaya::MipOutcome out = samaya::BranchAndBound(model, options, quiet).solve();
    ++models;
    violations += out.debug_cut_violations;
    CHECK(out.status == Status::kOptimal);
    CHECK(std::fabs(out.objective - ref.objective) <= 1e-6 * (1 + std::fabs(ref.objective)));
    const double target = cut_set_lp.objective;
    const double gap = target - out.root_bound;
    if (gap > 1e-6 * (1.0 + std::fabs(target))) {
      ++gaps;
      closed += std::min(1.0, (out.root_bound_cuts - out.root_bound) / gap);
      if (out.root_bound_cuts >= target - 1e-6 * (1.0 + std::fabs(target))) ++reached;
    }
  }
  const double average = gaps > 0 ? closed / gaps : 0.0;
  std::printf("  %d transshipment networks, %d with a cut-set gap: %.0f%% of it closed on "
              "average, bound reached on %d, %lld violations\n",
              models, gaps, 100.0 * average, reached, violations);
  CHECK_EQ(violations, 0);
  CHECK(models > 100);
  CHECK(gaps > 0);
}

TEST(mip_flow_cover_separates_the_cut_set_side_of_a_transit_node) {
  // A demand node that also forwards flow: x1 + x2 - x3 = d (inflows 1, 2; outflow 3), with
  // x_a <= M y_a and M >> d. At the LP point x1 = d, y1 = d / M, and the outflow arc is open
  // (y3 = 1) without flow. The <= side's flow cover, x_a <= d y_a + (M - d) y3, holds there;
  // only the >= side gives a violated cut: an inflow arc must be open (d y1 + d y2 >= d, or with
  // x2 for y2). Every cut must hold at every integer-feasible point.
  constexpr double kDemand = 7.0;
  constexpr double kBigM = 1000.0;  // d / M < 0.05: c-MIR cannot cut here.
  Model model;
  for (int a = 0; a < 3; ++a) {
    model.obj.push_back(1.0);
    model.col_lower.push_back(0);
    model.col_upper.push_back(kInf);
    model.col_type.push_back(samaya::VarType::kContinuous);
  }
  for (int a = 0; a < 3; ++a) {
    model.obj.push_back(20.0);
    model.col_lower.push_back(0);
    model.col_upper.push_back(1);
    model.col_type.push_back(samaya::VarType::kInteger);
  }
  std::vector<samaya::Triplet> t = {{0, 0, 1.0}, {0, 1, 1.0}, {0, 2, -1.0}};
  model.row_lower = {kDemand};
  model.row_upper = {kDemand};
  for (Index a = 0; a < 3; ++a) {
    t.push_back({1 + a, a, 1.0});
    t.push_back({1 + a, 3 + a, -kBigM});
    model.row_lower.push_back(-kInf);
    model.row_upper.push_back(0.0);
  }
  model.A = samaya::SparseMatrix::from_triplets(4, 6, std::move(t));
  const samaya::SparseMatrix at = model.A.transpose();
  const std::vector<double> x = {kDemand, 0, 0, kDemand / kBigM, 0, 1};
  const samaya::CutContext ctx{model, at, model.col_lower, model.col_upper, x, 4};
  std::vector<samaya::Cut> cuts;
  samaya::separate_aggregated_mir(ctx, cuts);
  double best_violation = 0.0;
  int invalid = 0;
  for (const samaya::Cut& c : cuts) {
    const auto activity = [&](const std::vector<double>& p) {
      double sum = 0.0;
      for (std::size_t k = 0; k < c.index.size(); ++k) sum += c.value[k] * p[c.index[k]];
      return sum;
    };
    best_violation = std::max(best_violation, c.lower - activity(x));
    // Integer-feasible points: any open arcs, the inflow d + x3 on one open inflow arc, and
    // x3 in {0, 5} when the outflow arc is open.
    for (int open = 0; open < 8; ++open) {
      for (int carrier = 0; carrier < 2; ++carrier) {
        if (!(open >> carrier & 1)) continue;
        for (const double out : {0.0, 5.0}) {
          if (out > 0.0 && !(open >> 2 & 1)) continue;
          std::vector<double> p(6, 0.0);
          p[carrier] = kDemand + out;
          p[2] = out;
          for (int a = 0; a < 3; ++a) p[3 + a] = open >> a & 1;
          if (activity(p) < c.lower - 1e-9) ++invalid;
        }
      }
    }
  }
  // The same LP point with the flow moved to inflow arc 2: a cut that wrote the unused inflow as
  // x2 (not d y2) no longer cuts it off, and the LP would just reroute round after round.
  const std::vector<double> rerouted = {0, kDemand, 0, 0, kDemand / kBigM, 1};
  double rerouted_violation = 0.0;
  for (const samaya::Cut& c : cuts) {
    double sum = 0.0;
    for (std::size_t k = 0; k < c.index.size(); ++k) sum += c.value[k] * rerouted[c.index[k]];
    rerouted_violation = std::max(rerouted_violation, c.lower - sum);
  }
  std::printf("  %zu cuts, largest violation of the LP point %.3f (rerouted %.3f), %d invalid\n",
              cuts.size(), best_violation, rerouted_violation, invalid);
  CHECK_EQ(invalid, 0);
  CHECK(best_violation > 0.5);  // d y1 + d y2 >= d is violated by d (1 - d / M) = 6.95.
  CHECK(rerouted_violation > 0.5);
}

TEST(mip_flow_cover_separates_the_cut_set_side_of_a_balance_row) {
  // x0 + x1 + x2 = d: an uncapacitated backup supply x0 (no binary) and two arcs x_a <= M y_a,
  // arc 1 also capped at d / 2. At the LP point x0 = x1 = d / 2, y1 = x1 / M: every flow of an
  // arc sits at a simple bound, so only the balance row is a start row (as on mc11, where the
  // start cap is used up by balance rows). Its >= side gives d y1 + d y2 + x0 >= d, violated by
  // d / 2 - d y1; its <= side cuts nothing here.
  constexpr double kDemand = 7.0;
  constexpr double kBigM = 1000.0;
  Model model;
  model.obj = {50.0, 1.0, 1.0, 20.0, 20.0};
  model.col_lower = {0, 0, 0, 0, 0};
  model.col_upper = {kInf, kDemand / 2, kInf, 1, 1};
  model.col_type = {samaya::VarType::kContinuous, samaya::VarType::kContinuous,
                    samaya::VarType::kContinuous, samaya::VarType::kInteger,
                    samaya::VarType::kInteger};
  std::vector<samaya::Triplet> t = {{0, 0, 1.0}, {0, 1, 1.0}, {0, 2, 1.0}, {1, 1, 1.0},
                                    {1, 3, -kBigM}, {2, 2, 1.0}, {2, 4, -kBigM}};
  model.row_lower = {kDemand, -kInf, -kInf};
  model.row_upper = {kDemand, 0.0, 0.0};
  model.A = samaya::SparseMatrix::from_triplets(3, 5, std::move(t));
  const samaya::SparseMatrix at = model.A.transpose();
  const std::vector<double> x = {kDemand / 2, kDemand / 2, 0, kDemand / 2 / kBigM, 0};
  const samaya::CutContext ctx{model, at, model.col_lower, model.col_upper, x, 3};
  std::vector<samaya::Cut> cuts;
  samaya::separate_aggregated_mir(ctx, cuts);
  double best_violation = 0.0;
  int invalid = 0;
  for (const samaya::Cut& c : cuts) {
    const auto activity = [&](const std::vector<double>& p) {
      double sum = 0.0;
      for (std::size_t k = 0; k < c.index.size(); ++k) sum += c.value[k] * p[c.index[k]];
      return sum;
    };
    best_violation = std::max(best_violation, c.lower - activity(x));
    // Integer-feasible points: any open arcs; the demand from one open arc (arc 1 up to d / 2),
    // the rest from the backup.
    for (int open = 0; open < 4; ++open) {
      for (int carrier = -1; carrier < 2; ++carrier) {
        if (carrier >= 0 && !(open >> carrier & 1)) continue;
        std::vector<double> p(5, 0.0);
        if (carrier >= 0) p[1 + carrier] = carrier == 0 ? kDemand / 2 : kDemand;
        p[0] = kDemand - p[1] - p[2];
        p[3] = open & 1;
        p[4] = open >> 1 & 1;
        if (activity(p) < c.lower - 1e-9) ++invalid;
      }
    }
  }
  // The same point with arc 2 carrying the flow instead: the cut-set inequality covers both
  // arcs, so it cuts this point off too.
  const std::vector<double> rerouted = {kDemand / 2, 0, kDemand / 2, 0, kDemand / 2 / kBigM};
  double rerouted_violation = 0.0;
  for (const samaya::Cut& c : cuts) {
    double sum = 0.0;
    for (std::size_t k = 0; k < c.index.size(); ++k) sum += c.value[k] * rerouted[c.index[k]];
    rerouted_violation = std::max(rerouted_violation, c.lower - sum);
  }
  std::printf("  %zu cuts, largest violation of the LP point %.3f (rerouted %.3f), %d invalid\n",
              cuts.size(), best_violation, rerouted_violation, invalid);
  CHECK_EQ(invalid, 0);
  CHECK(best_violation > 0.5);  // d y1 + d y2 + x0 >= d: violated by d / 2 - d y1 = 3.475.
  CHECK(rerouted_violation > 0.5);
}

TEST(mip_flow_cover_separates_single_node_big_m) {
  // One demand node: x1 + x2 + x3 = d, x_a <= M y_a with M >> d. The LP point sends d on arc 1
  // with y1 = d / M. The separator must cut it off (e.g. y1 + y2 + y3 >= 1), and every cut must
  // hold at every integer-feasible vertex (each nonempty set of open arcs, all flow on one of
  // them).
  constexpr double kDemand = 7.0;
  constexpr double kBigM = 1000.0;  // d / M < 0.05: c-MIR cannot cut here, only a flow cover.
  Model model;
  for (int a = 0; a < 3; ++a) {
    model.obj.push_back(1.0);
    model.col_lower.push_back(0);
    model.col_upper.push_back(kInf);
    model.col_type.push_back(samaya::VarType::kContinuous);
  }
  for (int a = 0; a < 3; ++a) {
    model.obj.push_back(20.0);
    model.col_lower.push_back(0);
    model.col_upper.push_back(1);
    model.col_type.push_back(samaya::VarType::kInteger);
  }
  std::vector<samaya::Triplet> t = {{0, 0, 1.0}, {0, 1, 1.0}, {0, 2, 1.0}};
  model.row_lower = {kDemand};
  model.row_upper = {kDemand};
  for (Index a = 0; a < 3; ++a) {
    t.push_back({1 + a, a, 1.0});
    t.push_back({1 + a, 3 + a, -kBigM});
    model.row_lower.push_back(-kInf);
    model.row_upper.push_back(0.0);
  }
  model.A = samaya::SparseMatrix::from_triplets(4, 6, std::move(t));
  const samaya::SparseMatrix at = model.A.transpose();
  const std::vector<double> x = {kDemand, 0, 0, kDemand / kBigM, 0, 0};
  const samaya::CutContext ctx{model, at, model.col_lower, model.col_upper, x, 4};
  std::vector<samaya::Cut> cuts;
  samaya::separate_aggregated_mir(ctx, cuts);
  REQUIRE(!cuts.empty());
  double best_violation = 0.0;
  int invalid = 0;
  for (const samaya::Cut& c : cuts) {
    const auto activity = [&](const std::vector<double>& p) {
      double s = 0.0;
      for (std::size_t k = 0; k < c.index.size(); ++k) s += c.value[k] * p[c.index[k]];
      return s;
    };
    best_violation = std::max(best_violation, c.lower - activity(x));
    for (int open = 1; open < 8; ++open) {
      for (int carrier = 0; carrier < 3; ++carrier) {
        if (!(open >> carrier & 1)) continue;
        std::vector<double> p(6, 0.0);
        p[carrier] = kDemand;
        for (int a = 0; a < 3; ++a) p[3 + a] = open >> a & 1;
        if (activity(p) < c.lower - 1e-9) ++invalid;
      }
    }
  }
  std::printf("  %zu cuts, largest violation of the LP point %.3f, %d invalid\n", cuts.size(),
              best_violation, invalid);
  CHECK_EQ(invalid, 0);
  CHECK(best_violation > 0.5);  // y1 + y2 + y3 >= 1 is violated by 1 - 7/1000.
}

TEST(mip_restart_after_root_matches_reference) {
  // Knapsacks with 18-22 items: once the root heuristics find a good solution, reduced-cost
  // fixing fixes many items at the root and the search restarts on the presolved rest. The
  // result must match the reference, whether or not it restarted.
  const samaya::Logger quiet(0);
  std::mt19937 rng(61);
  int restarted = 0;
  int failures = 0;
  int models = 0;
  for (int k = 0; k < 60; ++k) {
    const auto uniform_int = [&](int lo, int hi) {
      return std::uniform_int_distribution<int>(lo, hi)(rng);
    };
    Model model;
    model.sense = samaya::ObjSense::kMaximize;
    const int n = uniform_int(18, 22);
    std::vector<samaya::Triplet> t;
    double total = 0.0;
    for (int j = 0; j < n; ++j) {
      const int w = uniform_int(5, 40);
      model.obj.push_back(w * uniform_int(8, 12) / 10.0 + uniform_int(0, 3));
      model.col_lower.push_back(0);
      model.col_upper.push_back(1);
      model.col_type.push_back(samaya::VarType::kInteger);
      t.push_back({0, j, static_cast<double>(w)});
      total += w;
    }
    model.row_lower = {-kInf};
    model.row_upper = {std::floor(total * 0.45)};
    model.A = samaya::SparseMatrix::from_triplets(1, n, std::move(t));
    const ReferenceMilpResult ref = samaya::test::reference_milp(model);
    if (ref.status != ReferenceMilpResult::Status::kOptimal) continue;
    ++models;
    samaya::MipOptions options;
    options.rel_gap = 0.0;
    options.abs_gap = 1e-9;
    const samaya::MipOutcome out = samaya::BranchAndBound(model, options, quiet).solve();
    restarted += out.restarted;
    if (out.status != Status::kOptimal ||
        std::fabs(out.objective - ref.objective) > 1e-6 * (1 + std::fabs(ref.objective)) ||
        std::fabs(out.bound - out.objective) > 1e-6 * (1 + std::fabs(ref.objective))) {
      ++failures;
    }
  }
  std::printf("  %d knapsacks, %d restarted, %d failures\n", models, restarted, failures);
  CHECK_EQ(failures, 0);
  CHECK(models > 40);
  CHECK(restarted > 5);
}

TEST(mip_feasibility_jump_finds_verified_points) {
  // Feasibility Jump alone (no LP) on feasible random models and set partitioning: every point it
  // returns must satisfy all rows, bounds and integrality, and it must find most of them.
  std::mt19937 rng(81);
  int found = 0;
  int models = 0;
  int invalid = 0;
  for (int k = 0; k < 160; ++k) {
    const Model model = k % 2 == 0 ? random_milp(MilpFamily::kMixed, rng) : set_partitioning(rng);
    const ReferenceMilpResult ref = samaya::test::reference_milp(model);
    if (ref.status != ReferenceMilpResult::Status::kOptimal) continue;
    ++models;
    const samaya::FeasibilityJumpResult jump = samaya::feasibility_jump(
        model, model.num_rows(), model.col_lower, model.col_upper, 5.0, 2000000, 7u + k);
    if (!jump.found) continue;
    ++found;
    samaya::VerifyTolerances tol;
    tol.primal = 1e-7;
    if (!samaya::verify_primal(model, jump.x, tol).ok) ++invalid;
  }
  std::printf("  %d feasible models, feasibility jump found %d points, %d invalid\n", models,
              found, invalid);
  CHECK_EQ(invalid, 0);
  CHECK(models > 100);
  CHECK(found * 10 >= models * 8);
}

TEST(mip_random_models_match_reference_without_cuts) {
  std::mt19937 rng(21);
  samaya::Params params;
  params.log_level = 0;
  int failures = 0;
  for (int k = 0; k < 150; ++k) {
    const Model model = random_milp(k % 2 ? MilpFamily::kMixed : MilpFamily::kKnapsack, rng);
    const ReferenceMilpResult ref = samaya::test::reference_milp(model);
    if (ref.status != ReferenceMilpResult::Status::kOptimal) continue;
    samaya::MipOptions options;
    options.cuts = false;
    options.rel_gap = 0.0;
    options.abs_gap = 1e-9;
    const samaya::Logger quiet(0);
    const samaya::MipOutcome out = samaya::BranchAndBound(model, options, quiet).solve();
    if (out.status != Status::kOptimal ||
        std::fabs(out.objective - ref.objective) > 1e-6 * (1 + std::fabs(ref.objective))) {
      ++failures;
    }
  }
  CHECK_EQ(failures, 0);
}

TEST(mip_node_limit_reports_limit_with_valid_bound) {
  // A knapsack whose LP bound is loose needs more than one node.
  Model model;
  model.sense = samaya::ObjSense::kMaximize;
  const double w[] = {12, 7, 11, 8, 9, 6, 5, 13, 10, 4};
  const double v[] = {24, 13, 23, 15, 16, 11, 10, 25, 19, 7};
  std::vector<samaya::Triplet> t;
  for (int j = 0; j < 10; ++j) {
    model.obj.push_back(v[j]);
    model.col_lower.push_back(0);
    model.col_upper.push_back(1);
    model.col_type.push_back(samaya::VarType::kInteger);
    t.push_back({0, j, w[j]});
  }
  model.row_lower.push_back(-kInf);
  model.row_upper.push_back(26);
  model.A = samaya::SparseMatrix::from_triplets(1, 10, std::move(t));
  const ReferenceMilpResult ref = samaya::test::reference_milp(model);
  REQUIRE(ref.status == ReferenceMilpResult::Status::kOptimal);

  samaya::Params params;
  params.log_level = 0;
  params.presolve = false;
  params.mip_rel_gap = 0.0;
  params.node_limit = 1;
  const samaya::Result limited = samaya::Solver(params).solve(model);
  CHECK(limited.status == Status::kNodeLimit || limited.status == Status::kOptimal);
  CHECK(limited.dual_bound >= ref.objective - 1e-9);  // Maximization: bound from above.

  // An open-node cap stops the search cleanly with a valid bound instead of exhausting memory;
  // the soft cap switches to depth-first search first.
  const samaya::Logger quiet(0);
  samaya::MipOptions capped;
  capped.rel_gap = 0.0;
  capped.cuts = false;
  capped.max_open_nodes_soft = 2;
  capped.max_open_nodes = 3;
  const samaya::MipOutcome out = samaya::BranchAndBound(model, capped, quiet).solve();
  CHECK(out.status == Status::kNodeLimit || out.status == Status::kOptimal);
  CHECK(out.bound >= ref.objective - 1e-9);
  capped.max_open_nodes = 1000;
  const samaya::MipOutcome dfs = samaya::BranchAndBound(model, capped, quiet).solve();
  CHECK(dfs.status == Status::kOptimal);
  CHECK_NEAR(dfs.objective, ref.objective, 1e-9);

  params.node_limit = -1;
  const samaya::Result full = samaya::Solver(params).solve(model);
  CHECK(full.status == Status::kOptimal);
  CHECK(full.verified);
  CHECK_NEAR(full.objective, ref.objective, 1e-9);
}

TEST(mip_time_limit_is_respected) {
  // A market-share model (equality rows over binaries with slack penalties) is far too hard to
  // finish in the limit; the search must stop within it and report a valid bound.
  std::mt19937 rng(7);
  constexpr int kRows = 4;
  constexpr int kBinaries = 40;
  Model model;
  std::vector<samaya::Triplet> t;
  for (int j = 0; j < kBinaries; ++j) {
    model.obj.push_back(0.0);
    model.col_lower.push_back(0);
    model.col_upper.push_back(1);
    model.col_type.push_back(samaya::VarType::kInteger);
  }
  for (int i = 0; i < kRows; ++i) {
    double sum = 0.0;
    for (int j = 0; j < kBinaries; ++j) {
      const double a = std::uniform_int_distribution<int>(0, 99)(rng);
      t.push_back({i, j, a});
      sum += a;
    }
    // sum_j a_ij x_j + s_i^- - s_i^+ = floor(sum / 2), minimizing the slacks.
    for (int k = 0; k < 2; ++k) {
      t.push_back({i, static_cast<Index>(model.obj.size()), k == 0 ? 1.0 : -1.0});
      model.obj.push_back(1.0);
      model.col_lower.push_back(0);
      model.col_upper.push_back(kInf);
      model.col_type.push_back(samaya::VarType::kContinuous);
    }
    model.row_lower.push_back(std::floor(sum / 2));
    model.row_upper.push_back(std::floor(sum / 2));
  }
  model.A = samaya::SparseMatrix::from_triplets(kRows, static_cast<Index>(model.obj.size()),
                                                std::move(t));
  samaya::Params params;
  params.log_level = 0;
  params.time_limit = 0.5;
  params.node_limit = 5000000;  // A backstop so a search that ignores the clock still ends.
  const auto start = std::chrono::steady_clock::now();
  const samaya::Result result = samaya::Solver(params).solve(model);
  const double elapsed =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  CHECK(result.status == Status::kTimeLimit || result.status == Status::kOptimal);
  // Generous slack for sanitizer builds; the overshoot this guards against was seconds.
  CHECK(elapsed <= params.time_limit + 0.5);
  if (result.status == Status::kTimeLimit && std::isfinite(result.objective)) {
    CHECK(result.verified);
    CHECK(result.dual_bound <= result.objective + 1e-9);
  }
}

TEST(mip_nearly_integral_lp_solution_is_repaired_not_pruned) {
  // Regression (cases/mrpl.py crude, small): the root LP after cuts is integral within the
  // tolerance, but rounding a cargo binary (coefficient 130 in a tank balance) breaks the row by
  // more than the row tolerance. The node must be repaired by re-solving the continuous columns,
  // not pruned as infeasible. HiGHS: optimal 708749.0576.
  const Model model =
      samaya::read_mps(std::string(SAMAYA_TEST_DATA_DIR) + "/mrpl_crude_small.mps");
  samaya::Params params;
  params.log_level = 0;
  const samaya::Result result = samaya::Solver(params).solve(model);
  CHECK(result.status == Status::kOptimal);
  CHECK(result.verified);
  CHECK(std::fabs(result.objective - 708749.0576) <= 1e-4 * 708749.0576);
}

TEST(mip_dives_stop_on_columns_already_at_an_integral_bound) {
  // Regression (MIPLIB 2017 neos-2657525-crna, CC BY 4.0): a node dive met a general integer at
  // 1907.99999 with lower bound 1908. It counts as fractional, but rounding it changes no bound,
  // so the dive repeated the step without an LP iteration until the time limit (300k steps). The
  // search must reach its node limit well within the time limit. Known optimum 1.810748.
  const Model model =
      samaya::read_mps(std::string(SAMAYA_TEST_DATA_DIR) + "/neos-2657525-crna.mps");
  samaya::Params params;
  params.log_level = 0;
  params.threads = 1;
  params.node_limit = 5000;
  params.time_limit = 600.0;
  const samaya::Result result = samaya::Solver(params).solve(model);
  CHECK(result.status == Status::kNodeLimit);
  if (std::isfinite(result.objective)) {
    CHECK(result.verified);
    CHECK(result.objective >= 1.810748 - 1e-6);
  }
}

TEST(mip_infeasible_and_unbounded_models) {
  // 2x = 1 with x integer: LP feasible, integer infeasible.
  Model infeasible = samaya::read_mps_from_string(
      "NAME X\nROWS\n N obj\n E r\nCOLUMNS\n MARKER 'MARKER' 'INTORG'\n x obj 1 r 2\n"
      " MARKER 'MARKER' 'INTEND'\nRHS\n rhs r 1\nBOUNDS\n UP b x 10\nENDATA\n");
  samaya::Params params;
  params.log_level = 0;
  CHECK(samaya::Solver(params).solve(infeasible).status == Status::kInfeasible);
  params.presolve = false;
  CHECK(samaya::Solver(params).solve(infeasible).status == Status::kInfeasible);

  // max x + y with x integer free: unbounded relaxation.
  Model unbounded = samaya::read_mps_from_string(
      "NAME X\nOBJSENSE\n MAX\nROWS\n N obj\n L r\nCOLUMNS\n MARKER 'MARKER' 'INTORG'\n"
      " x obj 1 r 1\n MARKER 'MARKER' 'INTEND'\n y obj 1 r -1\nRHS\n rhs r 3\nBOUNDS\n"
      " PL b x\nENDATA\n");
  const samaya::Result r = samaya::Solver(params).solve(unbounded);
  CHECK(r.status == Status::kUnbounded || r.status == Status::kInfeasibleOrUnbounded);
}
