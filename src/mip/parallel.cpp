// Parallel tree search. After the sequential ramp-up the open nodes move to a pool shared by
// `threads` workers: the calling search (worker 0, which keeps the sub-MIP heuristics) and
// copies of it that own their LP and simplex. A node's bound changes are a shared, immutable
// path, so any worker can process any node. Workers share the pool, the incumbent and the
// global bound through one mutex; pseudocosts start from the ramp-up's and are then per worker.

#include <algorithm>
#include <atomic>
#include <cmath>
#include <condition_variable>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <thread>
#include <vector>

#include "mip/branch_and_bound.hpp"
#include "mip/search_constants.hpp"

namespace samaya {

struct BranchAndBound::Shared {
  std::mutex mutex;
  std::condition_variable wake;
  std::multimap<double, Node> pool;
  // Per worker: the smallest bound among the nodes it holds (kInf when it holds none). With the
  // pool this gives the global bound.
  std::vector<double> held;
  int workers = 1;
  int idle = 0;
  bool stop = false;
  bool unbounded = false;
  bool stopped = false;
  bool gap_closed = false;
  std::atomic<long long> nodes{0};
  double incumbent_value = kInf;
  std::vector<double> incumbent;
  std::atomic<double> incumbent_hint{kInf};  // Read without the lock to skip needless syncs.
};

void BranchAndBound::adopt(const BranchAndBound& master, int id) {
  worker_id_ = id;
  original_rows_ = master.original_rows_;
  timer_ = master.timer_;
  root_lower_ = master.root_lower_;
  root_upper_ = master.root_upper_;
  for (Index j = 0; j < n_; ++j) {
    lower_[j] = root_lower_[j];
    upper_[j] = root_upper_[j];
    lp_.lower[j] = lower_[j] / scaling_.col[j];
    lp_.upper[j] = upper_[j] / scaling_.col[j];
  }
  root_basis_ = master.root_basis_;
  down_locks_ = master.down_locks_;
  up_locks_ = master.up_locks_;
  for (int d = 0; d < 2; ++d) {
    pc_sum_[d] = master.pc_sum_[d];
    pc_count_[d] = master.pc_count_[d];
    pc_total_sum_[d] = master.pc_total_sum_[d];
    pc_total_count_[d] = master.pc_total_count_[d];
  }
  incumbent_value_ = master.incumbent_value_;
  incumbent_ = master.incumbent_;
  // The master's proofs hold for every worker: they use only the rows and the root bounds.
  conflicts_ = master.conflicts_;
  for (Conflict& c : conflicts_) c.last_used = 0;
  rng_.seed(12345u + static_cast<unsigned>(id));
  next_dive_rule_ = id;  // Workers start with different diving rules.
}

void BranchAndBound::publish_incumbent() {
  std::lock_guard<std::mutex> lock(shared_->mutex);
  if (incumbent_value_ < shared_->incumbent_value) {
    shared_->incumbent_value = incumbent_value_;
    shared_->incumbent = incumbent_;
    shared_->incumbent_hint.store(incumbent_value_, std::memory_order_relaxed);
  }
}

void BranchAndBound::pull_incumbent() {
  if (shared_->incumbent_hint.load(std::memory_order_relaxed) >= incumbent_value_) return;
  std::lock_guard<std::mutex> lock(shared_->mutex);
  if (shared_->incumbent_value < incumbent_value_) {
    incumbent_value_ = shared_->incumbent_value;
    incumbent_ = shared_->incumbent;
  }
}

void BranchAndBound::worker_loop(Shared& sh) {
  const auto id = static_cast<std::size_t>(worker_id_);
  std::optional<Node> current;
  double last_log = timer_.seconds();
  const auto held_bound = [&]() {
    double b = current ? current->bound : kInf;
    for (const Node& node : dive_stack_) b = std::min(b, node.bound);
    return b;
  };
  // Under the lock: return every held node to the pool.
  const auto give_back = [&]() {
    if (current) {
      sh.pool.emplace(current->bound, std::move(*current));
      current.reset();
    }
    for (Node& node : dive_stack_) sh.pool.emplace(node.bound, std::move(node));
    dive_stack_.clear();
    sh.held[id] = kInf;
  };
  const auto finish = [&](bool* flag) {
    std::lock_guard<std::mutex> lock(sh.mutex);
    if (flag != nullptr) *flag = true;
    sh.stop = true;
    give_back();
    sh.wake.notify_all();
  };

  for (;;) {
    pull_incumbent();
    if (!current && !dive_stack_.empty()) {
      current = std::move(dive_stack_.back());
      dive_stack_.pop_back();
    }
    if (!current) {
      std::unique_lock<std::mutex> lock(sh.mutex);
      sh.held[id] = kInf;
      for (;;) {
        if (sh.stop) return;
        if (!sh.pool.empty()) {
          auto it = sh.pool.begin();
          current = std::move(it->second);
          sh.pool.erase(it);
          sh.held[id] = current->bound;
          shared_stored_ = sh.pool.size();
          break;
        }
        if (sh.idle + 1 == sh.workers) {
          // Every other worker is waiting and the pool is empty: the tree is exhausted.
          sh.stop = true;
          sh.wake.notify_all();
          return;
        }
        ++sh.idle;
        sh.wake.wait(lock);
        --sh.idle;
      }
    }
    const double node_bound = effective_bound(current->bound);
    if (node_bound >= cutoff()) {
      pruned_bound_ = std::min(pruned_bound_, node_bound);
      current.reset();
      continue;
    }
    {
      std::lock_guard<std::mutex> lock(sh.mutex);
      if (sh.stop) {
        give_back();
        return;
      }
      sh.held[id] = held_bound();
      double best = sh.pool.empty() ? kInf : sh.pool.begin()->first;
      for (const double b : sh.held) best = std::min(best, b);
      best = effective_bound(best);
      const double incumbent = sh.incumbent_value;
      if (incumbent < kInf &&
          incumbent - best <= std::max(options_.abs_gap, options_.rel_gap * std::fabs(incumbent))) {
        sh.gap_closed = true;
        sh.stop = true;
        give_back();
        sh.wake.notify_all();
        return;
      }
      if (time_up() || sh.pool.size() >= options_.max_open_nodes ||
          (options_.node_limit >= 0 && sh.nodes.load() >= options_.node_limit)) {
        sh.stopped = true;
        sh.stop = true;
        give_back();
        sh.wake.notify_all();
        return;
      }
      shared_stored_ = sh.pool.size();
    }

    ++sh.nodes;
    ++outcome_.nodes;
    std::vector<Node> children;
    const NodeResult result = process_node(*current, children);
    if (id == 0 && log_.level() >= 1 && timer_.seconds() - last_log >= 5.0) {
      last_log = timer_.seconds();
      std::lock_guard<std::mutex> lock(sh.mutex);
      double best = sh.pool.empty() ? kInf : sh.pool.begin()->first;
      for (const double b : sh.held) best = std::min(best, b);
      log_.log(1, "MIP nodes %lld, open %zu, bound %.10g, incumbent %.10g, %d threads, %.1f s",
               sh.nodes.load(), sh.pool.size(), sense_ * effective_bound(best),
               sense_ * sh.incumbent_value, sh.workers, timer_.seconds());
    }
    if (result == NodeResult::kUnbounded) {
      finish(&sh.unbounded);
      return;
    }
    if (result == NodeResult::kStopped) {
      finish(&sh.stopped);
      return;
    }
    if (result == NodeResult::kFailed) incomplete_ = true;
    if (result != NodeResult::kBranched) {
      current.reset();
      continue;
    }
    const std::size_t dive = plunge_child(children);
    Node& other = children[1 - dive];
    Node& next = children[dive];
    if (!dive_stack_.empty() || shared_stored_ >= options_.max_open_nodes_soft) {
      // Memory mode, as in the sequential search: depth first on a local stack.
      dive_stack_.push_back(std::move(other));
      current = std::move(next);
      continue;
    }
    bool keep = false;
    {
      std::lock_guard<std::mutex> lock(sh.mutex);
      sh.pool.emplace(other.bound, std::move(other));
      if (incumbent_value_ == kInf) {
        keep = next.depth < search::kMaxPlungeDepth;
      } else {
        const double best = sh.pool.begin()->first;
        keep = effective_bound(next.bound) <=
               best + search::kPlungeGapFraction * (incumbent_value_ - best);
      }
      if (keep) {
        sh.held[id] = next.bound;
      } else {
        sh.pool.emplace(next.bound, std::move(next));
        sh.held[id] = kInf;
      }
      shared_stored_ = sh.pool.size();
    }
    if (keep) {
      current = std::move(next);
      sh.wake.notify_one();
    } else {
      current.reset();
      sh.wake.notify_all();
    }
  }
}

void BranchAndBound::run_parallel(bool& unbounded, bool& stopped, bool& gap_closed) {
  const int threads = options_.threads;
  Shared sh;
  sh.workers = threads;
  sh.held.assign(static_cast<std::size_t>(threads), kInf);
  for (auto& [bound, node] : open_) sh.pool.emplace(bound, std::move(node));
  open_.clear();
  for (Node& node : dive_stack_) sh.pool.emplace(node.bound, std::move(node));
  dive_stack_.clear();
  sh.incumbent_value = incumbent_value_;
  sh.incumbent = incumbent_;
  sh.incumbent_hint.store(incumbent_value_);
  sh.nodes.store(outcome_.nodes);
  log_.log(1, "MIP: %d threads from node %lld, %zu open nodes", threads, outcome_.nodes,
           sh.pool.size());

  // Build every worker before any starts: adopt() reads this search's state.
  const Logger quiet(0);
  MipOptions worker_options = options_;
  worker_options.threads = 1;
  worker_options.sub_mip_heuristics = false;
  worker_options.debug_solution.clear();
  worker_options.start.clear();
  std::vector<std::unique_ptr<BranchAndBound>> workers;
  for (int k = 1; k < threads; ++k) {
    workers.push_back(std::make_unique<BranchAndBound>(model_, worker_options, quiet));
    workers.back()->adopt(*this, k);
    workers.back()->shared_ = &sh;
  }
  shared_ = &sh;
  worker_id_ = 0;
  std::vector<std::thread> running;
  running.reserve(workers.size());
  for (auto& worker : workers) {
    running.emplace_back([&sh, w = worker.get()] { w->worker_loop(sh); });
  }
  worker_loop(sh);
  for (std::thread& t : running) t.join();
  shared_ = nullptr;
  shared_stored_ = 0;

  for (const auto& worker : workers) {
    outcome_.lp_iterations += worker->outcome_.lp_iterations;
    outcome_.strong_branching_iterations += worker->outcome_.strong_branching_iterations;
    outcome_.heuristic_solutions += worker->outcome_.heuristic_solutions;
    outcome_.heuristic_lp_iterations += worker->outcome_.heuristic_lp_iterations;
    outcome_.reduced_cost_fixings += worker->outcome_.reduced_cost_fixings;
    outcome_.conflicts_found += worker->outcome_.conflicts_found;
    outcome_.conflict_prunes += worker->outcome_.conflict_prunes;
    outcome_.conflict_infeasible_lps += worker->outcome_.conflict_infeasible_lps;
    outcome_.farkas_proofs += worker->outcome_.farkas_proofs;
    outcome_.strong_branching_conflicts += worker->outcome_.strong_branching_conflicts;
    pruned_bound_ = std::min(pruned_bound_, worker->pruned_bound_);
    incomplete_ = incomplete_ || worker->incomplete_;
  }
  outcome_.nodes = sh.nodes.load();
  outcome_.threads_used = threads;
  incumbent_value_ = sh.incumbent_value;
  incumbent_ = std::move(sh.incumbent);
  for (auto& [bound, node] : sh.pool) open_.emplace(bound, std::move(node));
  unbounded = sh.unbounded;
  stopped = sh.stopped;
  gap_closed = sh.gap_closed;
}

}  // namespace samaya
