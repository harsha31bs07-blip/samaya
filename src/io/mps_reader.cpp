// MPS / QPS reader.
//
// Supports the sections NAME, OBJSENSE, OBJNAME, ROWS, COLUMNS (with INTORG/INTEND markers),
// RHS, RANGES, BOUNDS (UP LO FX FR MI PL BV LI UI), QUADOBJ, QMATRIX, QSECTION (objective only)
// and ENDATA. Section headers start in column 1; data lines start with whitespace.
//
// The file is first read as free MPS (whitespace-separated fields). If that fails, it is read
// again as fixed MPS, where data fields sit in columns 2-3, 5-12, 15-22, 25-36, 40-47 and 50-61
// and names may contain spaces (Netlib's FORPLAN, for example).
//
// Conventions (matching common practice in CPLEX / Gurobi / HiGHS readers):
//   * The first N row is the objective unless OBJNAME names another; other N rows are dropped.
//   * An RHS entry on the objective row sets the objective offset to minus that value.
//   * Columns default to [0, +inf), including integer columns declared with markers.
//   * UP with a negative value on a column whose lower bound was not set makes the lower -inf.
//   * |value| >= 1e30 is treated as infinite.

#include <algorithm>
#include <charconv>
#include <filesystem>
#include <cstdint>
#include <fstream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "samaya/io.hpp"

namespace samaya {
namespace {

enum class Section : std::uint8_t {
  kNone,
  kObjSense,
  kObjName,
  kRows,
  kColumns,
  kRhs,
  kRanges,
  kBounds,
  kQuadObj,   // One triangle of Q.
  kQMatrix,   // Full symmetric Q.
  kEnd,
};

enum class RowKind : std::uint8_t { kE, kL, kG };

bool iequals(std::string_view a, std::string_view b) {
  if (a.size() != b.size()) return false;
  for (std::size_t i = 0; i < a.size(); ++i) {
    const auto lower = [](char c) { return (c >= 'A' && c <= 'Z') ? static_cast<char>(c + 32) : c; };
    if (lower(a[i]) != lower(b[i])) return false;
  }
  return true;
}

class MpsParser {
 public:
  explicit MpsParser(bool fixed_format) : fixed_(fixed_format) {}

  Model parse(std::string_view text) {
    std::size_t pos = 0;
    while (pos <= text.size() && section_ != Section::kEnd) {
      std::size_t end = text.find('\n', pos);
      if (end == std::string_view::npos) end = text.size();
      std::string_view line = text.substr(pos, end - pos);
      if (!line.empty() && line.back() == '\r') line.remove_suffix(1);
      ++line_no_;
      parse_line(line);
      if (end == text.size()) break;
      pos = end + 1;
    }
    return finish();
  }

 private:
  [[noreturn]] void fail(const std::string& message) const { throw ParseError(message, line_no_); }

  void tokenize(std::string_view line) {
    tokens_.clear();
    std::size_t i = 0;
    while (i < line.size()) {
      while (i < line.size() && (line[i] == ' ' || line[i] == '\t')) ++i;
      const std::size_t start = i;
      while (i < line.size() && line[i] != ' ' && line[i] != '\t') ++i;
      if (i > start) tokens_.push_back(line.substr(start, i - start));
    }
  }

  // Fixed MPS: the type field (columns 2-3) and five data fields, trimmed; empty fields are
  // skipped, which leaves the same token layout as free MPS.
  void tokenize_fixed(std::string_view line) {
    static constexpr std::size_t kStart[] = {1, 4, 14, 24, 39, 49};
    static constexpr std::size_t kEnd[] = {3, 12, 22, 36, 47, 61};
    tokens_.clear();
    for (std::size_t f = 0; f < std::size(kStart); ++f) {
      if (kStart[f] >= line.size()) break;
      std::string_view field = line.substr(kStart[f], std::min(kEnd[f], line.size()) - kStart[f]);
      while (!field.empty() && (field.front() == ' ' || field.front() == '\t')) field.remove_prefix(1);
      while (!field.empty() && (field.back() == ' ' || field.back() == '\t')) field.remove_suffix(1);
      if (!field.empty()) tokens_.push_back(field);
    }
  }

  double number(std::string_view s) const {
    if (!s.empty() && s.front() == '+') s.remove_prefix(1);
    double v = 0.0;
    const auto [ptr, ec] = std::from_chars(s.data(), s.data() + s.size(), v);
    if (ec != std::errc() || ptr != s.data() + s.size()) {
      fail("invalid number '" + std::string(s) + "'");
    }
    if (v >= kInfiniteBound) return kInf;
    if (v <= -kInfiniteBound) return -kInf;
    return v;
  }

  void parse_line(std::string_view line) {
    if (line.empty() || line.front() == '*') return;
    const bool header = line.front() != ' ' && line.front() != '\t';
    if (fixed_ && !header) {
      tokenize_fixed(line);
    } else {
      tokenize(line);
    }
    if (tokens_.empty()) return;
    if (header) {
      parse_header();
    } else {
      parse_data();
    }
  }

  void parse_header() {
    const std::string_view key = tokens_[0];
    if (key == "NAME") {
      if (tokens_.size() > 1) model_.name = std::string(tokens_[1]);
      section_ = Section::kNone;
    } else if (key == "OBJSENSE") {
      section_ = Section::kObjSense;
      if (tokens_.size() > 1) set_sense(tokens_[1]);
    } else if (key == "OBJNAME") {
      section_ = Section::kObjName;
      if (tokens_.size() > 1) set_objective_name(tokens_[1]);
    } else if (key == "ROWS") {
      section_ = Section::kRows;
    } else if (key == "COLUMNS") {
      section_ = Section::kColumns;
    } else if (key == "RHS") {
      section_ = Section::kRhs;
    } else if (key == "RANGES") {
      section_ = Section::kRanges;
    } else if (key == "BOUNDS") {
      section_ = Section::kBounds;
    } else if (key == "QUADOBJ") {
      section_ = Section::kQuadObj;
    } else if (key == "QMATRIX") {
      section_ = Section::kQMatrix;
    } else if (key == "QSECTION") {
      if (tokens_.size() < 2 || tokens_[1] != objective_name_) {
        fail("QSECTION for constraints (quadratic constraints) is not supported");
      }
      section_ = Section::kQMatrix;
    } else if (key == "ENDATA") {
      section_ = Section::kEnd;
    } else if (key == "SOS" || key == "QCMATRIX" || key == "INDICATORS" || key == "CSECTION") {
      fail("section " + std::string(key) + " is not supported");
    } else {
      fail("unknown section '" + std::string(key) + "'");
    }
  }

  void parse_data() {
    switch (section_) {
      case Section::kObjSense: set_sense(tokens_[0]); break;
      case Section::kObjName: set_objective_name(tokens_[0]); break;
      case Section::kRows: parse_row(); break;
      case Section::kColumns: parse_column(); break;
      case Section::kRhs: parse_rhs_or_range(/*is_range=*/false); break;
      case Section::kRanges: parse_rhs_or_range(/*is_range=*/true); break;
      case Section::kBounds: parse_bound(); break;
      case Section::kQuadObj: parse_quadratic(/*full_matrix=*/false); break;
      case Section::kQMatrix: parse_quadratic(/*full_matrix=*/true); break;
      case Section::kNone:
      case Section::kEnd: fail("data line outside of a section");
    }
  }

  void set_sense(std::string_view s) {
    if (iequals(s, "MAX") || iequals(s, "MAXIMIZE")) {
      model_.sense = ObjSense::kMaximize;
    } else if (iequals(s, "MIN") || iequals(s, "MINIMIZE")) {
      model_.sense = ObjSense::kMinimize;
    } else {
      fail("unknown objective sense '" + std::string(s) + "'");
    }
  }

  void set_objective_name(std::string_view name) {
    if (!rows_seen_.empty()) fail("OBJNAME must precede ROWS");
    objective_name_ = std::string(name);
  }

  void parse_row() {
    if (tokens_.size() < 2) fail("ROWS entry needs a type and a name");
    const std::string_view type = tokens_[0];
    std::string name(tokens_[1]);
    if (!rows_seen_.insert(name).second) fail("duplicate row '" + name + "'");

    if (iequals(type, "N")) {
      if (objective_name_.empty()) objective_name_ = name;
      if (name != objective_name_) free_rows_.insert(std::move(name));
      return;
    }
    RowKind kind;
    if (iequals(type, "E")) {
      kind = RowKind::kE;
    } else if (iequals(type, "L")) {
      kind = RowKind::kL;
    } else if (iequals(type, "G")) {
      kind = RowKind::kG;
    } else {
      fail("unknown row type '" + std::string(type) + "'");
    }
    row_index_.emplace(name, static_cast<Index>(row_kind_.size()));
    row_kind_.push_back(kind);
    model_.row_names.push_back(std::move(name));
  }

  // Returns the row index, or -1 for the objective and -2 for dropped free rows.
  Index find_row(std::string_view name) const {
    const std::string key(name);
    if (key == objective_name_) return -1;
    if (const auto it = row_index_.find(key); it != row_index_.end()) return it->second;
    if (free_rows_.count(key) != 0) return -2;
    fail("unknown row '" + key + "'");
  }

  Index find_col(std::string_view name) const {
    const auto it = col_index_.find(std::string(name));
    if (it == col_index_.end()) fail("unknown column '" + std::string(name) + "'");
    return it->second;
  }

  Index add_or_get_col(std::string_view name) {
    const auto [it, inserted] = col_index_.emplace(std::string(name), model_.num_cols());
    if (inserted) {
      model_.col_names.emplace_back(name);
      model_.obj.push_back(0.0);
      model_.col_lower.push_back(0.0);
      model_.col_upper.push_back(kInf);
      model_.col_type.push_back(in_integer_block_ ? VarType::kInteger : VarType::kContinuous);
      lower_set_.push_back(0);
    }
    return it->second;
  }

  void parse_column() {
    if (tokens_.size() >= 3 && tokens_[1] == "'MARKER'") {
      if (tokens_[2] == "'INTORG'") {
        in_integer_block_ = true;
      } else if (tokens_[2] == "'INTEND'") {
        in_integer_block_ = false;
      } else {
        fail("unknown marker " + std::string(tokens_[2]));
      }
      return;
    }
    if (tokens_.size() != 3 && tokens_.size() != 5) fail("COLUMNS entry needs 3 or 5 fields");
    const Index col = add_or_get_col(tokens_[0]);
    for (std::size_t k = 1; k + 1 < tokens_.size(); k += 2) {
      const Index row = find_row(tokens_[k]);
      const double value = number(tokens_[k + 1]);
      if (row == -1) {
        model_.obj[col] += value;
      } else if (row >= 0) {
        triplets_.push_back({row, col, value});
      }
    }
  }

  void parse_rhs_or_range(bool is_range) {
    // Optional leading set name makes the field count odd.
    const std::size_t first = tokens_.size() % 2 == 1 ? 1 : 0;
    if (tokens_.size() - first != 2 && tokens_.size() - first != 4) {
      fail(is_range ? "RANGES entry has the wrong number of fields"
                    : "RHS entry has the wrong number of fields");
    }
    if (rhs_.size() < row_kind_.size()) {
      rhs_.resize(row_kind_.size(), 0.0);
      range_.resize(row_kind_.size(), 0.0);
      has_range_.resize(row_kind_.size(), 0);
    }
    for (std::size_t k = first; k + 1 < tokens_.size(); k += 2) {
      const Index row = find_row(tokens_[k]);
      const double value = number(tokens_[k + 1]);
      if (is_range) {
        if (row < 0) fail("RANGES entry on an objective or free row");
        range_[row] = value;
        has_range_[row] = 1;
      } else if (row == -1) {
        model_.obj_offset = -value;
      } else if (row >= 0) {
        rhs_[row] = value;
      }
    }
  }

  void parse_bound() {
    if (tokens_.size() < 2) fail("BOUNDS entry is too short");
    const std::string_view type = tokens_[0];
    const bool needs_value = !(iequals(type, "FR") || iequals(type, "MI") || iequals(type, "PL") ||
                               iequals(type, "BV"));

    // Layout is TYPE [set] COLUMN [value]; detect whether the optional set name is present.
    std::size_t col_pos = 1;
    if (needs_value) {
      if (tokens_.size() == 4) col_pos = 2;
      else if (tokens_.size() != 3) fail("BOUNDS entry has the wrong number of fields");
    } else if (tokens_.size() == 4) {
      col_pos = 2;
    } else if (tokens_.size() == 3 && col_index_.count(std::string(tokens_[2])) != 0) {
      col_pos = 2;
    } else if (tokens_.size() > 4) {
      fail("BOUNDS entry has the wrong number of fields");
    }
    const Index col = find_col(tokens_[col_pos]);
    const double value = needs_value ? number(tokens_[col_pos + 1]) : 0.0;
    double& lower = model_.col_lower[col];
    double& upper = model_.col_upper[col];

    if (iequals(type, "UP") || iequals(type, "UI")) {
      upper = value;
      if (value < 0.0 && lower == 0.0 && !lower_set_[col]) lower = -kInf;
      if (iequals(type, "UI")) model_.col_type[col] = VarType::kInteger;
    } else if (iequals(type, "LO") || iequals(type, "LI")) {
      lower = value;
      lower_set_[col] = 1;
      if (iequals(type, "LI")) model_.col_type[col] = VarType::kInteger;
    } else if (iequals(type, "FX")) {
      lower = upper = value;
      lower_set_[col] = 1;
    } else if (iequals(type, "FR")) {
      lower = -kInf;
      upper = kInf;
      lower_set_[col] = 1;
    } else if (iequals(type, "MI")) {
      lower = -kInf;
      lower_set_[col] = 1;
    } else if (iequals(type, "PL")) {
      upper = kInf;
    } else if (iequals(type, "BV")) {
      lower = 0.0;
      upper = 1.0;
      lower_set_[col] = 1;
      model_.col_type[col] = VarType::kInteger;
    } else if (iequals(type, "SC")) {
      fail("semi-continuous bounds (SC) are not supported");
    } else {
      fail("unknown bound type '" + std::string(type) + "'");
    }
  }

  void parse_quadratic(bool full_matrix) {
    if (tokens_.size() != 3) fail("quadratic entry needs 3 fields");
    const Index a = find_col(tokens_[0]);
    const Index b = find_col(tokens_[1]);
    const double value = number(tokens_[2]);
    if (full_matrix) {
      // Both triangles are listed; keep the lower one.
      if (a >= b) q_triplets_.push_back({a, b, value});
    } else {
      q_triplets_.push_back({std::max(a, b), std::min(a, b), value});
    }
  }

  Model finish() {
    if (objective_name_.empty() && row_kind_.empty() && model_.num_cols() == 0) {
      fail("no ROWS or COLUMNS found");
    }
    const std::size_t m = row_kind_.size();
    rhs_.resize(m, 0.0);
    range_.resize(m, 0.0);
    has_range_.resize(m, 0);
    model_.row_lower.resize(m);
    model_.row_upper.resize(m);
    for (std::size_t i = 0; i < m; ++i) {
      const double b = rhs_[i];
      double lo = -kInf;
      double up = kInf;
      switch (row_kind_[i]) {
        case RowKind::kE: lo = up = b; break;
        case RowKind::kL: up = b; break;
        case RowKind::kG: lo = b; break;
      }
      if (has_range_[i]) {
        const double r = range_[i];
        const double ar = r < 0.0 ? -r : r;
        switch (row_kind_[i]) {
          case RowKind::kE:
            if (r >= 0.0) up = b + ar;
            else lo = b - ar;
            break;
          case RowKind::kL: lo = b - ar; break;
          case RowKind::kG: up = b + ar; break;
        }
      }
      model_.row_lower[i] = lo;
      model_.row_upper[i] = up;
    }

    const Index rows = static_cast<Index>(m);
    const Index cols = model_.num_cols();
    model_.A = SparseMatrix::from_triplets(rows, cols, std::move(triplets_));
    if (!q_triplets_.empty()) {
      model_.Q = SparseMatrix::from_triplets(cols, cols, std::move(q_triplets_));
    }
    line_no_ = 0;
    if (const std::string error = model_.validate(); !error.empty()) fail(error);
    return std::move(model_);
  }

  bool fixed_;
  Model model_;
  Section section_ = Section::kNone;
  long long line_no_ = 0;
  std::vector<std::string_view> tokens_;

  std::string objective_name_;
  std::unordered_set<std::string> rows_seen_;
  std::unordered_set<std::string> free_rows_;
  std::unordered_map<std::string, Index> row_index_;
  std::unordered_map<std::string, Index> col_index_;
  std::vector<RowKind> row_kind_;
  std::vector<double> rhs_;
  std::vector<double> range_;
  std::vector<std::uint8_t> has_range_;
  std::vector<std::uint8_t> lower_set_;
  bool in_integer_block_ = false;

  std::vector<Triplet> triplets_;
  std::vector<Triplet> q_triplets_;
};

}  // namespace

Model read_mps_from_string(std::string_view text) {
  try {
    return MpsParser(/*fixed_format=*/false).parse(text);
  } catch (const ParseError& free_error) {
    try {
      return MpsParser(/*fixed_format=*/true).parse(text);
    } catch (const ParseError&) {
      throw free_error;  // The free-format message is usually the more useful one.
    }
  }
}

Model read_mps(const std::string& path) {
  // The path is UTF-8; through std::filesystem::path MSVC opens it with the wide-character API.
  const std::filesystem::path file(std::u8string(path.begin(), path.end()));
  std::ifstream in(file, std::ios::binary);
  if (!in) throw std::runtime_error("cannot open '" + path + "'");
  std::ostringstream buffer;
  buffer << in.rdbuf();
  const std::string text = std::move(buffer).str();
  return read_mps_from_string(text);
}

}  // namespace samaya
