#pragma once

#include "core/log.hpp"
#include "samaya/solver.hpp"

namespace samaya {

// Ranging for an optimal, verified LP answer `result` of `model` (Result::status kOptimal,
// col_value and row_dual filled). Returns available = false with a message for anything else.
//
// The ranges belong to a basis, and the reported answer (from presolve, barrier or PDLP) has
// none on the original model, so one is rebuilt from it. With dual degeneracy that basis can
// have other duals or binding bounds than the reported answer. On success the answer in `result`
// is therefore replaced by the basis's own solution, after it passes the verifier with the same
// objective, so that the reported duals, binding bounds and ranges always belong together.
Ranging compute_ranging(const Model& model, Result& result, const Params& params,
                        const Logger& log);

}  // namespace samaya
