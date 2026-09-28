// samaya command-line interface.
//
//   samaya [options] model.mps
//
// With --json, a single JSON object summarizing the run is printed as the last line of stdout;
// bench/harness.py relies on that.

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <shellapi.h>
#endif

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include "samaya.hpp"

namespace {

// Paths are UTF-8 strings throughout. Opening them through std::filesystem::path makes MSVC use
// the wide-character file APIs, so names outside the ANSI code page (e.g. Devanagari) work.
std::filesystem::path utf8_path(const std::string& s) {
  return std::filesystem::path(std::u8string(s.begin(), s.end()));
}

#ifdef _WIN32
// The command line as UTF-8. Windows passes argv in the ANSI code page, which loses characters.
std::vector<std::string> utf8_arguments() {
  int n = 0;
  LPWSTR* wide = CommandLineToArgvW(GetCommandLineW(), &n);
  std::vector<std::string> out;
  for (int i = 0; i < n; ++i) {
    const int len = WideCharToMultiByte(CP_UTF8, 0, wide[i], -1, nullptr, 0, nullptr, nullptr);
    std::string arg(static_cast<std::size_t>(len > 1 ? len - 1 : 0), '\0');
    if (len > 1) WideCharToMultiByte(CP_UTF8, 0, wide[i], -1, arg.data(), len, nullptr, nullptr);
    out.push_back(std::move(arg));
  }
  LocalFree(wide);
  return out;
}
#endif

constexpr int kExitOk = 0;
constexpr int kExitReadError = 1;
constexpr int kExitUsage = 2;

void print_usage(std::FILE* out) {
  std::fprintf(out,
               "Usage: samaya [options] <model.mps>\n"
               "\n"
               "Options:\n"
               "  --stats               Print model statistics and exit without solving\n"
               "  --json                Print a JSON run summary as the last line of output\n"
               "  --solution <file>     Write status, objective, primal and dual values to a file\n"
               "  --mip-start <file>    Start a MILP from a solution file written by --solution\n"
               "                        (columns matched by name; e.g. yesterday's plan)\n"
               "  --time-limit <sec>    Wall-clock time limit\n"
               "  --threads <n>         Worker threads (0 = all cores)\n"
               "  --mip-gap <gap>       Relative MIP gap tolerance\n"
               "  --lp-method <m>       auto | dual | primal | barrier | pdlp | concurrent\n"
               "  --no-presolve         Disable presolve\n"
               "  --gpu                 Run PDLP's iterations on an NVIDIA GPU when available\n"
               "                        (with --lp-method pdlp; the CPU otherwise)\n"
               "  --gpu-info            Print, as JSON, whether this build has CUDA and which GPU\n"
               "                        it would use, and exit\n"
               "  --log-level <0-3>     Verbosity (default 1)\n"
               "  --version             Print the version and exit\n"
               "  --help                Print this message and exit\n");
}

bool parse_lp_method(std::string_view s, samaya::LpMethod& out) {
  using samaya::LpMethod;
  if (s == "auto") out = LpMethod::kAuto;
  else if (s == "dual") out = LpMethod::kDualSimplex;
  else if (s == "primal") out = LpMethod::kPrimalSimplex;
  else if (s == "barrier") out = LpMethod::kBarrier;
  else if (s == "pdlp") out = LpMethod::kPdlp;
  else if (s == "concurrent") out = LpMethod::kConcurrent;
  else return false;
  return true;
}

std::string json_string(std::string_view s) {
  std::string out = "\"";
  for (const char c : s) {
    switch (c) {
      case '"': out += "\\\""; break;
      case '\\': out += "\\\\"; break;
      case '\n': out += "\\n"; break;
      case '\t': out += "\\t"; break;
      default:
        if (static_cast<unsigned char>(c) < 0x20) {
          char buf[8];
          std::snprintf(buf, sizeof buf, "\\u%04x", c);
          out += buf;
        } else {
          out += c;
        }
    }
  }
  return out + "\"";
}

std::string json_number(double v) {
  if (!std::isfinite(v)) return "null";
  char buf[32];
  std::snprintf(buf, sizeof buf, "%.17g", v);
  return buf;
}

// Plain-text solution file: a header, then one line per column and per row.
bool write_solution(const std::string& file, const samaya::Model& model,
                    const samaya::Result& result) {
  std::ofstream out(utf8_path(file));
  if (!out) return false;
  out.precision(17);
  out << "status " << samaya::to_string(result.status) << "\n";
  out << "verified " << (result.verified ? "yes" : "no") << "\n";
  if (std::isfinite(result.objective)) out << "objective " << result.objective << "\n";
  const bool has_duals = !result.row_dual.empty();
  if (!result.col_value.empty()) {
    out << "\ncolumns " << model.num_cols() << "\n# name value" << (has_duals ? " reduced_cost" : "")
        << "\n";
    for (samaya::Index j = 0; j < model.num_cols(); ++j) {
      out << (model.col_names.empty() ? "c" + std::to_string(j) : model.col_names[j]) << " "
          << result.col_value[j];
      if (has_duals) out << " " << result.col_dual[j];
      out << "\n";
    }
    out << "\nrows " << model.num_rows() << "\n# name activity" << (has_duals ? " dual" : "")
        << "\n";
    for (samaya::Index i = 0; i < model.num_rows(); ++i) {
      out << (model.row_names.empty() ? "r" + std::to_string(i) : model.row_names[i]) << " "
          << result.row_activity[i];
      if (has_duals) out << " " << result.row_dual[i];
      out << "\n";
    }
  }
  if (!result.infeasibility_certificate.empty()) {
    out << "\ninfeasibility_certificate " << model.num_rows() << "\n";
    for (samaya::Index i = 0; i < model.num_rows(); ++i) {
      out << (model.row_names.empty() ? "r" + std::to_string(i) : model.row_names[i]) << " "
          << result.infeasibility_certificate[i] << "\n";
    }
  }
  if (!result.unbounded_ray.empty()) {
    out << "\nunbounded_ray " << model.num_cols() << "\n";
    for (samaya::Index j = 0; j < model.num_cols(); ++j) {
      out << (model.col_names.empty() ? "c" + std::to_string(j) : model.col_names[j]) << " "
          << result.unbounded_ray[j] << "\n";
    }
  }
  return static_cast<bool>(out);
}

// Reads the column values of a --solution file into one value per model column, matched by name
// (NaN where the file has no value). Returns the number of columns matched, or -1 if the file
// cannot be read.
int read_mip_start(const std::string& file, const samaya::Model& model,
                   std::vector<double>& start) {
  std::ifstream in(utf8_path(file));
  if (!in) return -1;
  std::unordered_map<std::string, samaya::Index> index;
  for (samaya::Index j = 0; j < model.num_cols(); ++j) {
    index.emplace(model.col_names.empty() ? "c" + std::to_string(j) : model.col_names[j], j);
  }
  start.assign(static_cast<std::size_t>(model.num_cols()), std::nan(""));
  int matched = 0;
  std::string line;
  bool in_columns = false;
  while (std::getline(in, line)) {
    if (line.rfind("columns ", 0) == 0) {
      in_columns = true;
      continue;
    }
    if (!in_columns || line.empty() || line.front() == '#') {
      if (in_columns && line.empty()) break;
      continue;
    }
    std::istringstream fields(line);
    std::string name;
    double value = 0.0;
    if (!(fields >> name >> value)) continue;
    const auto it = index.find(name);
    if (it == index.end()) continue;
    start[static_cast<std::size_t>(it->second)] = value;
    ++matched;
  }
  return matched;
}

}  // namespace

int main(int argc, char** argv) {
#ifdef _WIN32
  SetConsoleOutputCP(CP_UTF8);
  std::vector<std::string> utf8_args = utf8_arguments();
  std::vector<char*> utf8_argv;
  for (std::string& arg : utf8_args) utf8_argv.push_back(arg.data());
  utf8_argv.push_back(nullptr);
  argc = static_cast<int>(utf8_args.size());
  argv = utf8_argv.data();
#endif
  samaya::Params params;
  bool stats_only = false;
  bool json = false;
  std::string path;
  std::string solution_path;
  std::string start_path;

  for (int i = 1; i < argc; ++i) {
    const std::string_view arg = argv[i];
    const auto next_value = [&](const char* name) -> const char* {
      if (i + 1 >= argc) {
        std::fprintf(stderr, "samaya: %s requires a value\n", name);
        std::exit(kExitUsage);
      }
      return argv[++i];
    };
    const auto next_number = [&](const char* name) -> double {
      const char* text = next_value(name);
      char* end = nullptr;
      const double v = std::strtod(text, &end);
      if (end == text || *end != '\0') {
        std::fprintf(stderr, "samaya: invalid value '%s' for %s\n", text, name);
        std::exit(kExitUsage);
      }
      return v;
    };

    if (arg == "--help" || arg == "-h") {
      print_usage(stdout);
      return kExitOk;
    } else if (arg == "--gpu-info") {
      const samaya::GpuInfo gpu = samaya::gpu_info();
      std::printf("{\"cuda\":%s,\"gpu\":%s}\n", gpu.built_with_cuda ? "true" : "false",
                  gpu.device.empty() ? "null" : json_string(gpu.device).c_str());
      return kExitOk;
    } else if (arg == "--version") {
      std::printf("samaya %s\n", samaya::version());
      return kExitOk;
    } else if (arg == "--stats") {
      stats_only = true;
    } else if (arg == "--json") {
      json = true;
    } else if (arg == "--solution") {
      solution_path = next_value("--solution");
    } else if (arg == "--mip-start") {
      start_path = next_value("--mip-start");
    } else if (arg == "--time-limit") {
      params.time_limit = next_number("--time-limit");
    } else if (arg == "--threads") {
      params.threads = static_cast<int>(next_number("--threads"));
    } else if (arg == "--mip-gap") {
      params.mip_rel_gap = next_number("--mip-gap");
    } else if (arg == "--lp-method") {
      const char* method = next_value("--lp-method");
      if (!parse_lp_method(method, params.lp_method)) {
        std::fprintf(stderr, "samaya: unknown LP method '%s'\n", method);
        return kExitUsage;
      }
    } else if (arg == "--no-presolve") {
      params.presolve = false;
    } else if (arg == "--gpu") {
      params.use_gpu = true;
    } else if (arg == "--log-level") {
      params.log_level = static_cast<int>(next_number("--log-level"));
    } else if (!arg.empty() && arg.front() == '-') {
      std::fprintf(stderr, "samaya: unknown option '%s'\n", argv[i]);
      print_usage(stderr);
      return kExitUsage;
    } else if (path.empty()) {
      path = arg;
    } else {
      std::fprintf(stderr, "samaya: more than one model file given\n");
      return kExitUsage;
    }
  }
  if (path.empty()) {
    print_usage(stderr);
    return kExitUsage;
  }

  const auto read_start = std::chrono::steady_clock::now();
  samaya::Model model;
  try {
    model = samaya::read_mps(path);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "samaya: error reading '%s': %s\n", path.c_str(), e.what());
    if (json) {
      std::printf("{\"instance\":%s,\"status\":\"read_error\",\"message\":%s}\n",
                  json_string(path).c_str(), json_string(e.what()).c_str());
    }
    return kExitReadError;
  }
  const double read_seconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - read_start).count();

  const samaya::ModelStats stats = samaya::compute_stats(model);
  if (params.log_level > 0 || stats_only) {
    samaya::print_stats(std::cout, model, stats);
    std::cout << "Read time " << read_seconds << " s\n";
    std::cout.flush();
  }

  if (!start_path.empty()) {
    const int matched = read_mip_start(start_path, model, params.mip_start);
    if (matched < 0) {
      std::fprintf(stderr, "samaya: cannot read start '%s'\n", start_path.c_str());
      return kExitReadError;
    }
    if (params.log_level > 0) {
      std::printf("MIP start: %d of %d columns from %s\n", matched, model.num_cols(),
                  start_path.c_str());
    }
  }

  samaya::Result result;
  if (!stats_only) result = samaya::Solver(params).solve(model);

  if (!solution_path.empty() && !stats_only && !write_solution(solution_path, model, result)) {
    std::fprintf(stderr, "samaya: cannot write solution to '%s'\n", solution_path.c_str());
  }

  if (json) {
    std::printf(
        "{\"instance\":%s,\"name\":%s,\"class\":\"%s\",\"rows\":%d,\"cols\":%d,\"nnz\":%lld,"
        "\"integers\":%d,\"status\":\"%s\",\"verified\":%s,\"objective\":%s,"
        "\"dual_bound\":%s,\"max_primal_violation\":%s,\"max_dual_violation\":%s,"
        "\"simplex_iterations\":%lld,\"barrier_iterations\":%lld,\"read_seconds\":%s,"
        "\"solve_seconds\":%s,\"nodes\":%lld,"
        "\"message\":%s}\n",
        json_string(path).c_str(), json_string(model.name).c_str(),
        samaya::to_string(model.problem_class()), stats.rows, stats.cols,
        static_cast<long long>(stats.nnz), stats.integers, samaya::to_string(result.status),
        result.verified ? "true" : "false", json_number(result.objective).c_str(),
        json_number(result.dual_bound).c_str(), json_number(result.max_primal_violation).c_str(),
        json_number(result.max_dual_violation).c_str(), result.simplex_iterations,
        result.barrier_iterations,
        json_number(read_seconds).c_str(), json_number(result.solve_seconds).c_str(),
        result.nodes, json_string(result.message).c_str());
  }
  return kExitOk;
}
