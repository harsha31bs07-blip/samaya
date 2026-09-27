#pragma once

#include <memory>
#include <vector>

#include "samaya/sparse_matrix.hpp"

namespace samaya {
namespace gpu {

// The per-iteration arithmetic of PDHG (solve_pdlp) on the GPU. The host keeps every decision
// of the algorithm (termination, restarts, infeasibility and unboundedness detection) and reads
// the iterates back only when it evaluates them, so CPU and GPU runs follow the same logic and
// differ only in the order of floating-point sums.
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

 private:
  struct Impl;
  explicit PdhgDevice(std::unique_ptr<Impl> impl);
  std::unique_ptr<Impl> impl_;
};

}  // namespace gpu
}  // namespace samaya
