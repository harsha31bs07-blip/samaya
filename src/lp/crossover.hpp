#pragma once

#include <span>
#include <vector>

#include "lp/simplex.hpp"
#include "samaya/model.hpp"

namespace samaya {

// Crossover: a simplex basis guessed from a primal-dual point of an LP, e.g. the answer of the
// interior-point method or of PDLP (ARCHITECTURE.md WP3 M5).
//
// Returns n + m statuses in the computational form of make_problem (the structural columns, then
// the row activities), exactly m of them basic, for Simplex::solve(statuses). The simplex then
// finishes the job: it repairs a rank-deficient guess with logicals, and its iterations take the
// basis to a verified vertex. A good guess makes that a short cleanup instead of a full solve.
//
// x are the column values and y the row duals in the convention of verify_lp_optimality
// (d = c - A'y). Each variable gets a score for being basic from its distance to its nearest
// bound against the size of its reduced cost: near an optimum, complementary slackness makes one
// of the two close to zero. The m best are basic; the others sit at their nearest bound.
std::vector<VarStatus> crossover_basis(const Model& model, std::span<const double> x,
                                       std::span<const double> y);

}  // namespace samaya
