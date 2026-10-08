#pragma once

#include <memory>
#include <vector>

#include "samaya/sparse_matrix.hpp"

namespace samaya {
namespace gpu {

// The per-iteration arithmetic of PDHG (solve_pdlp) on the GPU. The host keeps every decision
// of the algorithm (termination, restarts, infeasibility and unboundedness detection), so CPU and
// GPU runs follow the same logic and differ only in the order of floating-point sums. The
// residuals and distances those decisions need are summed on the device (evaluate,
// restart_distance), so the iterates are read back only for an answer or a detection check.
//
// Problem: min cost'v s.t. [A -I] v = 0, lower <= v <= upper, v = (x, r) with n + m entries,
// y the m row duals. One step, from v, v_prev and y:
//   v_bar = v + theta (v - v_prev)
//   y_new = y + sigma (A v_bar_x - v_bar_r)
//   v_new = proj_[lower, upper](v - tau (cost + [A -I]' y_new))
//   v_erg, y_erg: running averages of v_new, y_new with weight 1 / (iter + 1).
// The header is free of CUDA types so that C++ sources can include it.
class PdhgDevice {
 public:
  // Uploads A (m x n, CSC), its transpose and the vectors. Returns null without a usable GPU.
  static std::unique_ptr<PdhgDevice> create(const SparseMatrix& A, const SparseMatrix& At,
                                            const std::vector<double>& cost,
                                            const std::vector<double>& lower,
                                            const std::vector<double>& upper);
  ~PdhgDevice();
  PdhgDevice(const PdhgDevice&) = delete;
  PdhgDevice& operator=(const PdhgDevice&) = delete;

  // Computes v_new, y_new and the averages for iteration `iter` (asynchronous).
  void step(double tau, double sigma, double theta, long long iter);
  // v_prev <- v, v <- v_new, y <- y_new (pointer swaps).
  void advance();
  // v_prev <- v: drops the extrapolation momentum (a restart).
  void restart();
  // Copies the iterates to the host (synchronizes). Any vector argument may be null.
  void download(std::vector<double>* v, std::vector<double>* v_new, std::vector<double>* y,
                std::vector<double>* y_new, std::vector<double>* v_erg,
                std::vector<double>* y_erg) const;

  // The data evaluate() maps the iterates back with: the equilibration's column and row scales
  // (x = col_scale v_x, r = v_r / row_scale, y_original = row_scale y) and the cost and bounds
  // before equilibration (n + m each). Call once after create().
  void set_original(const std::vector<double>& col_scale, const std::vector<double>& row_scale,
                    const std::vector<double>& cost, const std::vector<double>& lower,
                    const std::vector<double>& upper);

  enum class Iterate { kCurrent, kAverage };  // (v_new, y_new) or (v_erg, y_erg) after step().

  // Sums over an iterate in the original space, from which the host computes PDLP's relative
  // residuals and gap: ||A x - r||^2, the squared dual residual, the primal and dual objectives,
  // ||x||^2, ||r||^2 and ||y||^2. Summed in a fixed order, so a run is repeatable. Synchronizes.
  struct Sums {
    double primal_sq = 0.0;
    double dual_sq = 0.0;
    double primal_obj = 0.0;
    double dual_obj = 0.0;
    double x_sq = 0.0;
    double r_sq = 0.0;
    double y_sq = 0.0;
  };
  Sums evaluate(Iterate which);
  // Squared distances ||v - v_restart||^2 and ||y - y_restart||^2 of an iterate from the last
  // restart point (scaled space; the primal weight update). Synchronizes.
  void restart_distance(Iterate which, double& dv_sq, double& dy_sq);
  // Restarts at an iterate: it becomes the current point, the next one (so that advance() leaves
  // v_prev = v) and the restart point.
  void restart_to(Iterate which);

 private:
  struct Impl;
  explicit PdhgDevice(std::unique_ptr<Impl> impl);
  std::unique_ptr<Impl> impl_;
};

}  // namespace gpu
}  // namespace samaya
