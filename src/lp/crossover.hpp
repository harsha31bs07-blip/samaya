#pragma once

#include <span>
#include <vector>

#include "lp/simplex.hpp"
#include "samaya/model.hpp"

namespace samaya {

struct CrossoverStats {
  Index basic_structurals = 0;
  Index basic_equality_logicals = 0;
  Index basic_other_logicals = 0;
  Index basic_fixed_structurals = 0;  // Subset of basic_structurals.
  Index repairs = 0;  // Positions replaced across all attempts.
  int factorizations = 0;
  Index rank_deficiency = 0;  // From the last factorization; zero on success.
  double factorization_seconds = 0.0;
  // What snapping the nonbasic variables to their bounds does to the point (measured on the
  // final basis when it factorized). Superbasics: nonbasics whose point value is off the chosen
  // bound by more than 1e-9, 1e-6 and 1e-3 relative.
  Index superbasic_e9 = 0;
  Index superbasic_e6 = 0;
  Index superbasic_e3 = 0;
  Index infeasible_basics_before = 0;  // Basic values of the point outside their bounds.
  Index infeasible_basics_after = 0;   // Basic values after the nonbasics are snapped.
  double max_basic_change = 0.0;       // Largest change of a basic value from the snap.
  double max_violation_after = 0.0;    // Largest bound violation of a basic value after it.
};

// Crossover: a simplex basis guessed from a primal-dual point of an LP, e.g. the answer of the
// interior-point method or of PDLP (ARCHITECTURE.md WP3 M5).
//
// Returns n + m statuses in the computational form of make_problem (the structural columns, then
// the row activities), exactly m of them basic, for Simplex::solve(statuses). The selected basis
// is repaired with logicals to factorize without rank deficiency in the original model's space.
// Repairs stop after four factorizations or if an entering logical is already basic, leaving any
// remaining deficiency to the simplex. The simplex factorizes a scaled copy, so it can still
// occasionally repair the basis; its iterations take the basis to a verified vertex.
//
// x are the column values and y the row duals in the convention of verify_lp_optimality
// (d = c - A'y). Each variable gets a score for being basic from its distance to its nearest
// bound against the size of its reduced cost: near an optimum, complementary slackness makes one
// of the two close to zero. The m best are initially basic; the others, including columns leaving
// during rank repair, sit at the bound nearest to their value in the supplied point.
std::vector<VarStatus> crossover_basis(const Model& model, std::span<const double> x,
                                       std::span<const double> y,
                                       CrossoverStats* stats = nullptr);

}  // namespace samaya
