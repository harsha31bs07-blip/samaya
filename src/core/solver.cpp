#include "samaya/solver.hpp"

#include <algorithm>
#include <string>
#include <thread>

#include "core/log.hpp"
#include "core/xreal.hpp"
#include "lp/ipm.hpp"
#include "lp/lp_solver.hpp"
#include "lp/pdlp.hpp"
#include "mip/branch_and_bound.hpp"
#include "presolve/presolve.hpp"
#include "qp/miqp.hpp"
#include "samaya/verify.hpp"

#ifndef SAMAYA_VERSION_STRING
#define SAMAYA_VERSION_STRING "0.0.0"
#endif

namespace samaya {

const char* version() { return SAMAYA_VERSION_STRING; }

namespace {

constexpr double kPdlpTolerance = 1e-8;

Status to_status(SimplexStatus s) {
  switch (s) {
    case SimplexStatus::kOptimal: return Status::kOptimal;
    case SimplexStatus::kInfeasible: return Status::kInfeasible;
    case SimplexStatus::kUnbounded: return Status::kUnbounded;
    case SimplexStatus::kIterationLimit: return Status::kIterationLimit;
    case SimplexStatus::kTimeLimit: return Status::kTimeLimit;
    case SimplexStatus::kNumericalError: return Status::kNumericalError;
  }
  return Status::kNumericalError;
}

// Runs the independent verifier on an LP outcome. Returns true when the claim is proven (or is a
// limit status, which claims nothing).
bool verify_lp(const Model& model, const LpResult& lp, Result& result) {
  VerifyReport report;
  switch (lp.status) {
    case SimplexStatus::kOptimal:
      report = verify_lp_optimality(model, lp.col_value, lp.row_dual);
      break;
    case SimplexStatus::kInfeasible:
      report = verify_infeasibility(model, lp.dual_ray);
      break;
    case SimplexStatus::kUnbounded:
      report = verify_unbounded_ray(model, lp.primal_ray);
      break;
    default:
      return true;
  }
  result.max_primal_violation = std::max(report.max_bound_violation, report.max_row_violation);
  result.max_dual_violation = report.max_dual_violation;
  result.message = report.message;
  return report.ok;
}

// Solves an LP with the dual simplex. If the outcome fails verification the solve is repeated
// with tighter tolerances and then without scaling. Returns the last outcome; `verified` says
// whether it passed (limit statuses claim nothing and count as verified).
LpResult solve_lp_verified(const Model& model, const Params& params, const Logger& log,
                           Result& result, bool& verified) {
  LpSolveOptions attempts[3];
  attempts[0] = lp_options_from_params(params);
  attempts[1] = attempts[0];
  attempts[1].simplex.primal_tol = std::min(attempts[1].simplex.primal_tol, 1e-9);
  attempts[1].simplex.dual_tol = std::min(attempts[1].simplex.dual_tol, 1e-9);
  attempts[2] = attempts[1];
  attempts[2].scale = false;
  const int num_attempts = params.verify ? 3 : 1;

  LpResult lp;
  verified = false;
  long long iterations = 0;
  for (int a = 0; a < num_attempts; ++a) {
    lp = solve_lp(model, attempts[a], log);
    iterations += lp.iterations;
    if (!params.verify) break;
    verified = verify_lp(model, lp, result);
    if (verified) break;
    log.log(1, "Verification of %s result failed (%s); retrying with stricter settings",
            to_string(lp.status), result.message.c_str());
  }
  lp.iterations = iterations;
  return lp;
}

// Presolves the LP, solves the reduced model and maps an optimal solution back. Returns false
// (and leaves `lp` unspecified) when the caller should solve the original model instead: the
// reduced model was not solved to a verified optimum, or presolve found infeasibility (the
// original solve then produces a certificate for the original model).
bool solve_presolved_lp(const Model& model, const Params& params, const Logger& log,
                        Result& result, LpResult& lp, long long& iterations) {
  Presolve presolve(model);
  if (presolve.run() != PresolveStatus::kReduced) {
    log.log(1, "Presolve: infeasibility detected; solving the original model for a certificate");
    return false;
  }
  const Model& reduced = presolve.reduced();
  const PresolveStats& st = presolve.stats();
  log.log(1, "Presolve: %d rows, %d cols -> %d rows, %d cols, %lld nonzeros (%d passes)",
          model.num_rows(), model.num_cols(), reduced.num_rows(), reduced.num_cols(),
          static_cast<long long>(reduced.A.nnz()), st.passes);
  log.log(2,
          "Presolve: rows: %d empty, %d singleton, %d forcing, %d redundant, %d duplicate; "
          "cols: %d fixed, %d empty, %d dominated, %d free singleton; %d bounds tightened",
          st.empty_rows, st.singleton_rows, st.forcing_rows, st.redundant_rows,
          st.duplicate_rows, st.fixed_cols, st.empty_cols, st.dominated_cols,
          st.free_col_singletons, st.bounds_tightened);

  LpResult reduced_lp;
  if (reduced.num_cols() > 0) {
    bool verified = false;
    reduced_lp = solve_lp_verified(reduced, params, log, result, verified);
    iterations += reduced_lp.iterations;
    if (reduced_lp.status != SimplexStatus::kOptimal || !verified) {
      log.log(1, "Presolved model: %s; solving the original model",
              to_string(reduced_lp.status));
      return false;
    }
  }

  lp = LpResult{};
  lp.status = SimplexStatus::kOptimal;
  presolve.postsolve(reduced_lp.col_value, reduced_lp.row_dual, lp.col_value, lp.row_dual);
  const Index n = model.num_cols();
  const Index m = model.num_rows();
  lp.row_activity.assign(static_cast<std::size_t>(m), 0.0);
  model.A.multiply(lp.col_value, lp.row_activity);
  lp.col_dual.assign(static_cast<std::size_t>(n), 0.0);
  model.A.multiply_transpose(lp.row_dual, lp.col_dual);
  Xreal objective = model.obj_offset;  // double-double: the same bits on every platform
  for (Index j = 0; j < n; ++j) {
    lp.col_dual[j] = model.obj[j] - lp.col_dual[j];
    objective += Xreal(model.obj[j]) * lp.col_value[j];
  }
  lp.objective = static_cast<double>(objective);
  if (params.verify && !verify_lp(model, lp, result)) {
    log.log(1, "Postsolved solution failed verification (%s); solving the original model",
            result.message.c_str());
    return false;
  }
  return true;
}

// Fills an optimal LP/QP result from x and y (row duals): activities, reduced costs
// d = c + Qx - A'y and the verified objective.
void fill_continuous_optimum(const Model& model, std::vector<double> x, std::vector<double> y,
                             const VerifyReport& report, Result& result) {
  const Index n = model.num_cols();
  result.objective = report.objective;
  result.dual_bound = report.objective;
  result.max_primal_violation = std::max(report.max_bound_violation, report.max_row_violation);
  result.max_dual_violation = report.max_dual_violation;
  result.row_activity.assign(static_cast<std::size_t>(model.num_rows()), 0.0);
  model.A.multiply(x, result.row_activity);
  std::vector<double> aty(static_cast<std::size_t>(n), 0.0);
  model.A.multiply_transpose(y, aty);
  result.col_dual.assign(static_cast<std::size_t>(n), 0.0);
  for (Index j = 0; j < n; ++j) result.col_dual[j] = model.obj[j] - aty[j];
  const auto q_start = model.Q.col_start();
  const auto q_index = model.Q.row_index();
  const auto q_value = model.Q.values();
  for (Index j = 0; j < model.Q.cols(); ++j) {
    for (NnzIndex p = q_start[j]; p < q_start[j + 1]; ++p) {
      const Index i = q_index[p];
      result.col_dual[i] += q_value[p] * x[j];
      if (i != j) result.col_dual[j] += q_value[p] * x[i];
    }
  }
  result.col_value = std::move(x);
  result.row_dual = std::move(y);
}

// LP by the interior-point method (--lp-method barrier). Returns false when the caller should
// fall back to the dual simplex: the IPM did not converge, or its point failed the strict
// verifier (there is no crossover yet, so the point is not polished).
bool solve_barrier_lp(const Model& model, const Params& params, const Logger& log,
                      Result& result) {
  IpmOptions options;
  options.time_limit = params.time_limit;
  const IpmResult ipm = solve_ipm(model, options, log);
  result.barrier_iterations = ipm.iterations;
  if (ipm.status == IpmStatus::kTimeLimit) {
    result.status = Status::kTimeLimit;
    return false;
  }
  if (ipm.status != IpmStatus::kOptimal) {
    log.log(1, "Barrier: %s; using the dual simplex", to_string(ipm.status));
    return false;
  }
  // The last iterate first; the polished point if the iterate does not pass.
  VerifyReport report = verify_lp_optimality(model, ipm.x, ipm.y);
  if (report.ok || !params.verify) {
    result.status = Status::kOptimal;
    result.verified = params.verify;
    fill_continuous_optimum(model, ipm.x, ipm.y, report, result);
    return true;
  }
  if (!ipm.polished_x.empty()) {
    const VerifyReport polished = verify_lp_optimality(model, ipm.polished_x, ipm.polished_y);
    if (!polished.ok) log.log(2, "Barrier: polished point rejected (%s)", polished.message.c_str());
    if (polished.ok) {
      result.status = Status::kOptimal;
      result.verified = true;
      fill_continuous_optimum(model, ipm.polished_x, ipm.polished_y, polished, result);
      return true;
    }
  }
  log.log(1, "Barrier: point failed verification (%s); using the dual simplex",
          report.message.c_str());
  return false;
}

void solve_lp_model(const Model& model, const Params& params, const Logger& log, Result& result);

// Convex QP by the interior-point method. The answer is reported optimal only if
// verify_qp_optimality accepts it. If the IPM does not converge, the constraints alone are
// solved as an LP by the dual simplex: an infeasible QP then gets a verified Farkas certificate.
void solve_qp_model(const Model& model, const Params& params, const Logger& log,
                    Result& result) {
  const Timer timer;
  IpmOptions options;
  options.time_limit = params.time_limit;
  const IpmResult ipm = solve_ipm(model, options, log);
  result.barrier_iterations = ipm.iterations;
  if (ipm.status == IpmStatus::kNotConvex) {
    result.status = Status::kNotConvex;
    result.message = "Q is not positive semidefinite; only convex QP is solved";
    return;
  }
  if (ipm.status == IpmStatus::kOptimal) {
    // The last iterate first; the polished point if the iterate does not pass.
    const VerifyReport report = verify_qp_optimality(model, ipm.x, ipm.y);
    if (!params.verify || report.ok) {
      result.status = Status::kOptimal;
      result.verified = params.verify;
      fill_continuous_optimum(model, ipm.x, ipm.y, report, result);
      return;
    }
    if (!ipm.polished_x.empty()) {
      const VerifyReport polished = verify_qp_optimality(model, ipm.polished_x, ipm.polished_y);
      if (!polished.ok) log.log(2, "QP: polished point rejected (%s)", polished.message.c_str());
      if (polished.ok) {
        result.status = Status::kOptimal;
        result.verified = true;
        fill_continuous_optimum(model, ipm.polished_x, ipm.polished_y, polished, result);
        return;
      }
    }
    result.status = Status::kNumericalError;
    result.message = "interior-point solution failed verification: " + report.message;
    result.max_primal_violation = std::max(report.max_bound_violation, report.max_row_violation);
    result.max_dual_violation = report.max_dual_violation;
    return;
  }
  if (ipm.status == IpmStatus::kTimeLimit) {
    result.status = Status::kTimeLimit;
    return;
  }
  // Not converged: decide feasibility with the simplex on the constraints alone.
  Model feasibility = model;
  feasibility.Q = SparseMatrix();
  std::fill(feasibility.obj.begin(), feasibility.obj.end(), 0.0);
  Params lp_params = params;
  lp_params.time_limit = params.time_limit - timer.seconds();
  lp_params.lp_method = LpMethod::kDualSimplex;
  Result lp;
  solve_lp_model(feasibility, lp_params, log, lp);
  result.simplex_iterations = lp.simplex_iterations;
  if (lp.status == Status::kInfeasible) {
    result.status = Status::kInfeasible;
    result.verified = lp.verified;
    result.infeasibility_certificate = std::move(lp.infeasibility_certificate);
    return;
  }
  result.status = ipm.status == IpmStatus::kIterationLimit ? Status::kIterationLimit
                                                             : Status::kNumericalError;
  result.message = std::string("interior point: ") + to_string(ipm.status) +
                   (lp.status == Status::kOptimal ? " (the constraints are feasible)" : "");
}

// LP by PDLP (--lp-method pdlp; on the GPU with --gpu when built with SAMAYA_CUDA). Returns false
// when the caller should fall back to the dual simplex: PDLP stopped without an answer, or its
// answer failed the strict verifier (first-order solutions are less accurate and there is no
// crossover yet; ARCHITECTURE.md WP1 M4).
bool solve_pdlp_lp(const Model& model, const Params& params, const Logger& log, Result& result) {
  const Scaling scaled = compute_scaling(model.A);
  const LpProblem scaled_lp = make_problem(model, scaled);
  PdlpOptions pdlp_opts;
  // PDLP's strict tolerance (ARCHITECTURE.md WP1 M5): at the simplex's 1e-7 a first-order
  // method leaves objective errors near 1e-5; at 1e-8 its answers agree with the other engines.
  pdlp_opts.tol = std::min(params.dual_feasibility_tol, kPdlpTolerance);
  pdlp_opts.max_iterations = -1;  // Automatic, proportional to the problem size.
  pdlp_opts.time_limit = params.time_limit;
  pdlp_opts.use_gpu = params.use_gpu;
  const PdlpResult pdlp_result = solve_pdlp(scaled_lp, pdlp_opts, log);
  result.barrier_iterations = pdlp_result.iterations;

  // Map the PDLP result (scaled, minimization) to the model.
  LpResult lp;
  lp.status = pdlp_result.status;
  lp.iterations = pdlp_result.iterations;
  lp.stats = pdlp_result.stats;
  const Index n = model.num_cols();
  const double sense = model.sense == ObjSense::kMaximize ? -1.0 : 1.0;
  if (pdlp_result.x.size() == static_cast<std::size_t>(n + model.num_rows())) {
    lp.col_value.assign(pdlp_result.x.begin(), pdlp_result.x.begin() + n);
    lp.row_activity.assign(pdlp_result.x.begin() + n, pdlp_result.x.end());
    lp.row_dual = pdlp_result.y;
    scaled.unscale_cols(lp.col_value);
    scaled.unscale_row_duals(lp.row_dual);
    // PDLP's Lagrangian is c'x + y'(Ax - s), so its reduced costs are c + A'y: negate y for the
    // convention d = c - A'y, and apply the objective sense.
    for (double& y : lp.row_dual) y *= -sense;
    lp.col_dual.assign(static_cast<std::size_t>(n), 0.0);
    model.A.multiply_transpose(lp.row_dual, lp.col_dual);
    Xreal objective = model.obj_offset;  // double-double: the same bits on every platform
    for (Index j = 0; j < n; ++j) {
      lp.col_dual[j] = model.obj[j] - lp.col_dual[j];
      objective += Xreal(model.obj[j]) * lp.col_value[j];
    }
    lp.objective = static_cast<double>(objective);
  }
  if (lp.status == SimplexStatus::kInfeasible && !pdlp_result.dual_ray.empty()) {
    lp.dual_ray = pdlp_result.dual_ray;
    scaled.unscale_row_duals(lp.dual_ray);
    for (double& y : lp.dual_ray) y *= sense;
  } else if (lp.status == SimplexStatus::kUnbounded && !pdlp_result.primal_ray.empty()) {
    lp.primal_ray = pdlp_result.primal_ray;
    scaled.unscale_cols(lp.primal_ray);  // In the scaled space of x, like the solution.
  }
  if (lp.status == SimplexStatus::kTimeLimit) {
    result.status = Status::kTimeLimit;
    return false;
  }
  const bool answered = lp.status == SimplexStatus::kOptimal ||
                        lp.status == SimplexStatus::kInfeasible ||
                        lp.status == SimplexStatus::kUnbounded;
  const bool verified = answered && (!params.verify || verify_lp(model, lp, result));
  if (!verified) {
    log.log(1, "PDLP: %s%s; using the dual simplex", to_string(lp.status),
            answered ? " failed verification" : "");
    result.message.clear();
    return false;
  }
  result.status = to_status(lp.status);
  result.verified = params.verify;
  result.message.clear();
  result.col_value = std::move(lp.col_value);
  result.row_activity = std::move(lp.row_activity);
  if (lp.status == SimplexStatus::kOptimal) {
    result.objective = lp.objective;
    result.dual_bound = lp.objective;
    result.row_dual = std::move(lp.row_dual);
    result.col_dual = std::move(lp.col_dual);
  }
  result.infeasibility_certificate = std::move(lp.dual_ray);
  result.unbounded_ray = std::move(lp.primal_ray);
  return true;
}

// Mixed-integer convex QP by branch and bound over the QP relaxations. The returned solution is
// checked against the original model (bounds, rows, integrality); the optimality of the bound
// rests on the search, as for MILP.
void solve_miqp_model(const Model& model, const Params& params, const Logger& log,
                      Result& result) {
  MiqpOptions options;
  options.time_limit = params.time_limit;
  options.node_limit = params.node_limit;
  options.rel_gap = params.mip_rel_gap;
  options.abs_gap = params.mip_abs_gap;
  options.integrality_tol = params.integrality_tol;
  const MiqpOutcome out = solve_miqp(model, options, log);
  result.status = out.status;
  result.nodes = out.nodes;
  result.barrier_iterations = out.ipm_iterations;
  result.dual_bound = out.bound;
  if (out.status == Status::kNotConvex) {
    result.message = "Q is not positive semidefinite; only convex MIQP is solved";
  }
  if (out.x.empty()) return;
  VerifyTolerances tol;
  tol.integrality = params.integrality_tol;
  const VerifyReport report = verify_primal(model, out.x, tol);
  result.max_primal_violation = std::max(report.max_bound_violation, report.max_row_violation);
  result.objective = report.objective;
  result.col_value = out.x;
  result.row_activity.assign(static_cast<std::size_t>(model.num_rows()), 0.0);
  model.A.multiply(result.col_value, result.row_activity);
  if (params.verify) {
    result.verified = report.ok;
    if (!report.ok) {
      result.message = "solution failed verification: " + report.message;
      result.status = Status::kNumericalError;
    }
  }
  if (result.status == Status::kOptimal) {
    result.dual_bound = model.sense == ObjSense::kMinimize
                            ? std::min(result.dual_bound, result.objective)
                            : std::max(result.dual_bound, result.objective);
  }
}

void solve_lp_model(const Model& model, const Params& requested, const Logger& log,
                    Result& result) {
  const Timer timer;
  Params params = requested;
  if (requested.lp_method == LpMethod::kBarrier || requested.lp_method == LpMethod::kPdlp) {
    const bool solved = requested.lp_method == LpMethod::kBarrier
                            ? solve_barrier_lp(model, requested, log, result)
                            : solve_pdlp_lp(model, requested, log, result);
    // A time limit in the first method is final; after any other failure the dual simplex gets
    // what is left of the limit.
    if (solved || result.status == Status::kTimeLimit) return;
    params.time_limit = requested.time_limit - timer.seconds();
  }
  LpResult lp;
  long long iterations = 0;
  bool verified = false;
  if (params.presolve && solve_presolved_lp(model, params, log, result, lp, iterations)) {
    verified = params.verify;
  } else {
    lp = solve_lp_verified(model, params, log, result, verified);
    iterations += lp.iterations;
  }

  result.status = to_status(lp.status);
  result.simplex_iterations = iterations;
  result.verified = verified && (lp.status == SimplexStatus::kOptimal ||
                                 lp.status == SimplexStatus::kInfeasible ||
                                 lp.status == SimplexStatus::kUnbounded);
  if (params.verify && !verified) {
    result.message = std::string("result failed verification: ") + result.message;
    result.status = Status::kNumericalError;
  } else if (verified) {
    result.message.clear();
  }

  result.col_value = std::move(lp.col_value);
  result.row_activity = std::move(lp.row_activity);
  if (lp.status == SimplexStatus::kOptimal) {
    result.objective = lp.objective;
    result.dual_bound = lp.objective;
    result.row_dual = std::move(lp.row_dual);
    result.col_dual = std::move(lp.col_dual);
  }
  result.infeasibility_certificate = std::move(lp.dual_ray);
  result.unbounded_ray = std::move(lp.primal_ray);
}

// Presolves the MILP, runs branch-and-bound on the reduced model and maps the best solution
// back. The returned solution is checked against the original model (bounds, rows and
// integrality); the optimality of the bound rests on the search itself.
void solve_mip_model(const Model& model, const Params& params, const Logger& log,
                     Result& result) {
  const Timer timer;
  PresolveOptions presolve_options;
  presolve_options.mip = true;
  presolve_options.integrality_tol = params.integrality_tol;
  Presolve presolve(model, presolve_options);
  if (params.presolve) {
    if (presolve.run() == PresolveStatus::kInfeasible) {
      result.status = Status::kInfeasible;
      result.message = "presolve proved infeasibility";
      log.log(1, "Presolve: infeasible");
      return;
    }
    const Model& reduced = presolve.reduced();
    log.log(1, "Presolve: %d rows, %d cols (%d integer) -> %d rows, %d cols (%d integer)",
            model.num_rows(), model.num_cols(), model.num_integers(), reduced.num_rows(),
            reduced.num_cols(), reduced.num_integers());
  }
  const Model& work = params.presolve ? presolve.reduced() : model;

  MipOutcome outcome;
  if (work.num_cols() == 0) {
    outcome.status = Status::kOptimal;
    outcome.objective = outcome.bound = work.obj_offset;
  } else {
    MipOptions options;
    options.time_limit = params.time_limit - timer.seconds();
    options.node_limit = params.node_limit;
    options.rel_gap = params.mip_rel_gap;
    options.abs_gap = params.mip_abs_gap;
    options.integrality_tol = params.integrality_tol;
    options.threads = params.threads > 0
                          ? params.threads
                          : std::max(1, static_cast<int>(std::thread::hardware_concurrency()));
    if (params.mip_start.size() == static_cast<std::size_t>(model.num_cols())) {
      if (!params.presolve) {
        options.start = params.mip_start;
      } else {
        for (const Index j : presolve.col_map()) options.start.push_back(params.mip_start[j]);
      }
    }
    BranchAndBound search(work, options, log);
    outcome = search.solve();
  }
  result.status = outcome.status;
  result.nodes = outcome.nodes;
  result.simplex_iterations = outcome.lp_iterations;
  result.dual_bound = outcome.bound;
  log.log(1, "MIP: %lld nodes, %lld LP iterations (%lld in strong branching), %d heuristic "
          "solutions", outcome.nodes, outcome.lp_iterations, outcome.strong_branching_iterations,
          outcome.heuristic_solutions);

  const bool has_solution = !outcome.x.empty() || (work.num_cols() == 0 &&
                                                   outcome.status == Status::kOptimal);
  if (!has_solution) return;
  std::vector<double> x;
  if (params.presolve) {
    std::vector<double> unused_y;
    presolve.postsolve(outcome.x, {}, x, unused_y);
  } else {
    x = std::move(outcome.x);
  }
  VerifyTolerances tol;
  tol.integrality = params.integrality_tol;
  const VerifyReport report = verify_primal(model, x, tol);
  result.max_primal_violation = std::max(report.max_bound_violation, report.max_row_violation);
  result.objective = report.objective;
  result.col_value = std::move(x);
  result.row_activity.assign(static_cast<std::size_t>(model.num_rows()), 0.0);
  model.A.multiply(result.col_value, result.row_activity);
  if (params.verify) {
    result.verified = report.ok;
    if (!report.ok) {
      result.message = "solution failed verification: " + report.message;
      result.status = Status::kNumericalError;
    }
  }
  if (result.status == Status::kOptimal) {
    // The bound cannot pass the solution's objective.
    result.dual_bound = model.sense == ObjSense::kMinimize
                            ? std::min(result.dual_bound, result.objective)
                            : std::max(result.dual_bound, result.objective);
  }
}

}  // namespace

Result Solver::solve(const Model& model) const {
  const Timer timer;
  const Logger log(params_.log_level);
  Result result;

  if (std::string error = model.validate(); !error.empty()) {
    result.status = Status::kInvalidModel;
    result.message = std::move(error);
    log.log(1, "Invalid model: %s", result.message.c_str());
    result.solve_seconds = timer.seconds();
    return result;
  }

  log.log(1, "samaya %s: %s with %d rows, %d cols, %lld nonzeros", version(),
          to_string(model.problem_class()), model.num_rows(), model.num_cols(),
          static_cast<long long>(model.A.nnz()));

  if (std::string crossed = model.crossed_bounds(); !crossed.empty()) {
    // Infeasible by inspection; the crossed bound itself is the proof.
    result.status = Status::kInfeasible;
    result.verified = true;
    result.message = std::move(crossed);
    log.log(1, "Infeasible: %s", result.message.c_str());
    result.solve_seconds = timer.seconds();
    return result;
  }

  switch (model.problem_class()) {
    case ProblemClass::kLP:
      solve_lp_model(model, params_, log, result);
      break;
    case ProblemClass::kMILP:
      solve_mip_model(model, params_, log, result);
      break;
    case ProblemClass::kQP:
      solve_qp_model(model, params_, log, result);
      break;
    case ProblemClass::kMIQP:
      solve_miqp_model(model, params_, log, result);
      break;
  }

  result.solve_seconds = timer.seconds();
  if (result.status == Status::kOptimal) {
    log.log(1, "Optimal objective %.12g after %lld simplex iterations, %.3f s%s", result.objective,
            result.simplex_iterations, result.solve_seconds,
            result.verified ? " (verified)" : "");
  } else {
    log.log(1, "Status %s after %lld simplex iterations, %.3f s%s%s", to_string(result.status),
            result.simplex_iterations, result.solve_seconds, result.message.empty() ? "" : ": ",
            result.message.c_str());
  }
  return result;
}

}  // namespace samaya
