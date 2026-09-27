#pragma once

// The portable part of samaya Studio (no Windows APIs, unit-tested on Linux too): a small JSON
// reader and writer for the page's messages, the reader of samaya's solution file, and the
// independent re-check. The re-check parses the MPS file with its own code and recomputes rows,
// bounds, integrality and the objective from the solution alone; it shares no code with samaya's
// reader (src/io) or verifier (src/verify), so it checks them rather than repeating them.

#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace studio {

// JSON --------------------------------------------------------------------------------------------
struct Json {
  enum Kind { kNull, kBool, kNumber, kString, kArray, kObject };
  Kind kind = kNull;
  bool boolean = false;
  double number = 0.0;
  std::string string;
  std::vector<Json> array;
  std::vector<std::pair<std::string, Json>> object;

  const Json* get(std::string_view key) const;
  std::string str(std::string_view key, const std::string& fallback = "") const;
  double num(std::string_view key, double fallback) const;
  bool flag(std::string_view key, bool fallback) const;
};

// Parses a JSON text; throws std::runtime_error on malformed input.
Json parse_json(std::string_view text);
// A JSON string literal (quotes and escapes); the text is UTF-8 and stays UTF-8.
std::string json_string(std::string_view text);
// A JSON number; NaN and infinities become null (JSON has no such numbers).
std::string json_number(double v);

// samaya's solution file (apps/cli, --solution) -------------------------------------------------
struct SolutionFile {
  std::string status;
  bool verified = false;
  double objective = 0.0;
  bool has_objective = false;
  std::vector<std::pair<std::string, double>> columns;  // name, value
  struct Row {
    std::string name;
    double activity = 0.0;
    double dual = 0.0;
  };
  std::vector<Row> rows;
};
SolutionFile read_solution(const std::string& text);

// Independent re-check ------------------------------------------------------------------------
struct Recheck {
  bool available = false;  // false when there is no solution or the model is not free MPS
  std::string reason;      // why it is not available
  bool maximize = false;
  bool quadratic = false;  // QPS model: the objective check is left out (linear rows only)
  long long rows = 0, cols = 0, integers = 0, missing = 0;
  double row_violation = 0.0;  // worst relative row violation
  std::string row_name;
  double bound_violation = 0.0;
  double integrality = 0.0;
  double objective = 0.0;       // recomputed c'x + constant
  double objective_diff = 0.0;  // |recomputed - reported| / (1 + |reported|)
};
// model_text: the MPS file; the solution's objective is compared with the recomputed one.
Recheck recheck(const std::string& model_text, const SolutionFile& solution);

// Decodes standard base64 (the page sends files it builds, e.g. workbooks, this way).
std::string base64_decode(std::string_view text);

}  // namespace studio
