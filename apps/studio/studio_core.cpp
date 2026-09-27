#include "studio_core.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace studio {

// JSON ------------------------------------------------------------------------------------------

const Json* Json::get(std::string_view key) const {
  for (const auto& [k, v] : object) {
    if (k == key) return &v;
  }
  return nullptr;
}
std::string Json::str(std::string_view key, const std::string& fallback) const {
  const Json* v = get(key);
  return v && v->kind == kString ? v->string : fallback;
}
double Json::num(std::string_view key, double fallback) const {
  const Json* v = get(key);
  return v && v->kind == kNumber ? v->number : fallback;
}
bool Json::flag(std::string_view key, bool fallback) const {
  const Json* v = get(key);
  return v && v->kind == kBool ? v->boolean : fallback;
}

namespace {

class JsonParser {
 public:
  explicit JsonParser(std::string_view text) : s_(text) {}

  Json parse() {
    Json v = value();
    skip();
    if (i_ != s_.size()) fail("trailing characters");
    return v;
  }

 private:
  [[noreturn]] void fail(const char* what) const {
    throw std::runtime_error(std::string("JSON: ") + what + " at " + std::to_string(i_));
  }
  void skip() {
    while (i_ < s_.size() && (s_[i_] == ' ' || s_[i_] == '\t' || s_[i_] == '\n' || s_[i_] == '\r'))
      ++i_;
  }
  bool eat(char c) {
    skip();
    if (i_ < s_.size() && s_[i_] == c) {
      ++i_;
      return true;
    }
    return false;
  }
  Json value() {
    skip();
    if (i_ >= s_.size()) fail("unexpected end");
    const char c = s_[i_];
    Json v;
    if (c == '{') {
      ++i_;
      v.kind = Json::kObject;
      if (eat('}')) return v;
      do {
        skip();
        std::string key = string_literal();
        if (!eat(':')) fail("expected ':'");
        v.object.emplace_back(std::move(key), value());
      } while (eat(','));
      if (!eat('}')) fail("expected '}'");
    } else if (c == '[') {
      ++i_;
      v.kind = Json::kArray;
      if (eat(']')) return v;
      do {
        v.array.push_back(value());
      } while (eat(','));
      if (!eat(']')) fail("expected ']'");
    } else if (c == '"') {
      v.kind = Json::kString;
      v.string = string_literal();
    } else if (s_.substr(i_, 4) == "true") {
      i_ += 4;
      v.kind = Json::kBool;
      v.boolean = true;
    } else if (s_.substr(i_, 5) == "false") {
      i_ += 5;
      v.kind = Json::kBool;
    } else if (s_.substr(i_, 4) == "null") {
      i_ += 4;
    } else {
      const std::string tail(s_.substr(i_, std::min<std::size_t>(64, s_.size() - i_)));
      char* end = nullptr;
      v.number = std::strtod(tail.c_str(), &end);
      if (end == tail.c_str()) fail("bad value");
      v.kind = Json::kNumber;
      i_ += static_cast<std::size_t>(end - tail.c_str());
    }
    return v;
  }
  static void append_utf8(std::string& out, unsigned cp) {
    if (cp < 0x80) {
      out += static_cast<char>(cp);
    } else if (cp < 0x800) {
      out += static_cast<char>(0xC0 | (cp >> 6));
      out += static_cast<char>(0x80 | (cp & 0x3F));
    } else if (cp < 0x10000) {
      out += static_cast<char>(0xE0 | (cp >> 12));
      out += static_cast<char>(0x80 | ((cp >> 6) & 0x3F));
      out += static_cast<char>(0x80 | (cp & 0x3F));
    } else {
      out += static_cast<char>(0xF0 | (cp >> 18));
      out += static_cast<char>(0x80 | ((cp >> 12) & 0x3F));
      out += static_cast<char>(0x80 | ((cp >> 6) & 0x3F));
      out += static_cast<char>(0x80 | (cp & 0x3F));
    }
  }
  unsigned hex4() {
    if (i_ + 4 > s_.size()) fail("short \\u escape");
    unsigned v = 0;
    for (int k = 0; k < 4; ++k) {
      const char c = s_[i_++];
      v <<= 4;
      if (c >= '0' && c <= '9')
        v |= static_cast<unsigned>(c - '0');
      else if (c >= 'a' && c <= 'f')
        v |= static_cast<unsigned>(c - 'a' + 10);
      else if (c >= 'A' && c <= 'F')
        v |= static_cast<unsigned>(c - 'A' + 10);
      else
        fail("bad \\u escape");
    }
    return v;
  }
  std::string string_literal() {
    if (i_ >= s_.size() || s_[i_] != '"') fail("expected string");
    ++i_;
    std::string out;
    while (true) {
      if (i_ >= s_.size()) fail("unterminated string");
      const char c = s_[i_++];
      if (c == '"') break;
      if (c != '\\') {
        out += c;
        continue;
      }
      if (i_ >= s_.size()) fail("bad escape");
      const char e = s_[i_++];
      switch (e) {
        case '"':
          out += '"';
          break;
        case '\\':
          out += '\\';
          break;
        case '/':
          out += '/';
          break;
        case 'b':
          out += '\b';
          break;
        case 'f':
          out += '\f';
          break;
        case 'n':
          out += '\n';
          break;
        case 'r':
          out += '\r';
          break;
        case 't':
          out += '\t';
          break;
        case 'u': {
          unsigned cp = hex4();
          if (cp >= 0xD800 && cp < 0xDC00 && i_ + 6 <= s_.size() && s_[i_] == '\\' &&
              s_[i_ + 1] == 'u') {
            i_ += 2;
            const unsigned lo = hex4();
            cp = 0x10000 + ((cp - 0xD800) << 10) + (lo - 0xDC00);
          }
          append_utf8(out, cp);
          break;
        }
        default:
          fail("bad escape");
      }
    }
    return out;
  }

  std::string_view s_;
  std::size_t i_ = 0;
};

}  // namespace

Json parse_json(std::string_view text) {
  return JsonParser(text).parse();
}

std::string json_string(std::string_view text) {
  std::string out = "\"";
  for (const char c : text) {
    switch (c) {
      case '"':
        out += "\\\"";
        break;
      case '\\':
        out += "\\\\";
        break;
      case '\n':
        out += "\\n";
        break;
      case '\r':
        out += "\\r";
        break;
      case '\t':
        out += "\\t";
        break;
      default:
        if (static_cast<unsigned char>(c) < 0x20) {
          char buf[8];
          std::snprintf(buf, sizeof buf, "\\u%04x",
                        static_cast<unsigned>(static_cast<unsigned char>(c)));
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

// Solution file ---------------------------------------------------------------------------------

namespace {

std::vector<std::string_view> tokens(std::string_view line) {
  std::vector<std::string_view> out;
  std::size_t i = 0;
  while (i < line.size()) {
    while (i < line.size() && (line[i] == ' ' || line[i] == '\t' || line[i] == '\r')) ++i;
    const std::size_t start = i;
    while (i < line.size() && line[i] != ' ' && line[i] != '\t' && line[i] != '\r') ++i;
    if (i > start) out.push_back(line.substr(start, i - start));
  }
  return out;
}

double to_double(std::string_view s) {
  const std::string t(s);
  char* end = nullptr;
  const double v = std::strtod(t.c_str(), &end);
  if (end == t.c_str()) throw std::runtime_error("not a number: '" + t + "'");
  return v;
}

template <typename F>
void for_each_line(const std::string& text, F&& f) {
  std::size_t start = 0;
  while (start <= text.size()) {
    std::size_t end = text.find('\n', start);
    if (end == std::string::npos) end = text.size();
    f(std::string_view(text).substr(start, end - start));
    start = end + 1;
  }
}

}  // namespace

SolutionFile read_solution(const std::string& text) {
  SolutionFile sol;
  std::string section;
  for_each_line(text, [&](std::string_view line) {
    const auto t = tokens(line);
    if (t.empty() || t[0][0] == '#') return;
    if (t[0] == "columns" || t[0] == "rows" || t[0] == "infeasibility_certificate" ||
        t[0] == "unbounded_ray") {
      section = std::string(t[0]);
      return;
    }
    if (section.empty()) {
      if (t[0] == "status" && t.size() > 1) sol.status = std::string(t[1]);
      if (t[0] == "verified" && t.size() > 1) sol.verified = t[1] == "yes";
      if (t[0] == "objective" && t.size() > 1) {
        sol.objective = to_double(t[1]);
        sol.has_objective = true;
      }
    } else if (section == "columns" && t.size() >= 2) {
      sol.columns.emplace_back(std::string(t[0]), to_double(t[1]));
    } else if (section == "rows" && t.size() >= 2) {
      sol.rows.push_back(
          {std::string(t[0]), to_double(t[1]), t.size() > 2 ? to_double(t[2]) : 0.0});
    }
  });
  return sol;
}

// Independent re-check ------------------------------------------------------------------------

namespace {

constexpr double kInf = std::numeric_limits<double>::infinity();

// Neumaier's compensated sum: a running sum and its rounding error, so long sums of mixed
// magnitudes keep their small terms (plenty against the verifier's 1e-6 tolerances).
struct Sum {
  double s = 0.0;
  double c = 0.0;
  void add(double v) {
    const double t = s + v;
    c += std::fabs(s) >= std::fabs(v) ? (s - t) + v : (v - t) + s;
    s = t;
  }
  double value() const { return s + c; }
};

struct RowInfo {
  char type = 'N';
  double rhs = 0.0;
  double range = 0.0;
  bool has_range = false;
};

}  // namespace

Recheck recheck(const std::string& model_text, const SolutionFile& solution) {
  Recheck out;
  if (solution.columns.empty()) {
    out.reason = "no solution to re-check";
    return out;
  }
  std::unordered_map<std::string, std::size_t> row_index;
  std::vector<std::string> row_names;
  std::vector<RowInfo> rows;
  std::string objective_row;
  struct Entry {
    std::size_t row;
    double value;
  };
  std::unordered_map<std::string, std::size_t> col_index;
  std::vector<std::string> col_names;
  std::vector<std::vector<Entry>> col_entries;
  std::vector<double> lower, upper;
  std::vector<bool> is_integer;
  std::string section;
  bool integer_block = false;
  bool ok = true;

  const auto column = [&](std::string_view name) -> std::size_t {
    const std::string key(name);
    auto it = col_index.find(key);
    if (it != col_index.end()) return it->second;
    const std::size_t j = col_names.size();
    col_index.emplace(key, j);
    col_names.push_back(key);
    col_entries.emplace_back();
    lower.push_back(0.0);
    upper.push_back(kInf);
    is_integer.push_back(integer_block);
    return j;
  };

  for_each_line(model_text, [&](std::string_view line) {
    if (!ok || line.empty() || line[0] == '*') return;
    const auto t = tokens(line);
    if (t.empty()) return;
    const bool header = line[0] != ' ' && line[0] != '\t';
    if (header) {
      section = std::string(t[0]);
      if (section == "OBJSENSE" && t.size() > 1) out.maximize = t[1] == "MAX" || t[1] == "MAXIMIZE";
      if (section == "QUADOBJ" || section == "QMATRIX" || section == "QSECTION")
        out.quadratic = true;
      return;
    }
    try {
      if (section == "OBJSENSE") {
        out.maximize = t[0] == "MAX" || t[0] == "MAXIMIZE";
      } else if (section == "ROWS" && t.size() >= 2) {
        RowInfo r;
        r.type = t[0][0];
        if (r.type == 'N' && objective_row.empty()) objective_row = std::string(t[1]);
        row_index.emplace(std::string(t[1]), rows.size());
        row_names.emplace_back(t[1]);
        rows.push_back(r);
      } else if (section == "COLUMNS") {
        if (t.size() >= 3 && t[1] == "'MARKER'") {
          integer_block = t[2] == "'INTORG'";
          return;
        }
        const std::size_t j = column(t[0]);
        for (std::size_t k = 1; k + 1 < t.size(); k += 2) {
          auto it = row_index.find(std::string(t[k]));
          if (it == row_index.end()) throw std::runtime_error("unknown row");
          col_entries[j].push_back({it->second, to_double(t[k + 1])});
        }
      } else if (section == "RHS" || section == "RANGES") {
        // An optional set name comes first: pairs start at 1 when the count is odd.
        for (std::size_t k = t.size() % 2; k + 1 < t.size(); k += 2) {
          auto it = row_index.find(std::string(t[k]));
          if (it == row_index.end()) throw std::runtime_error("unknown row");
          if (section == "RHS") {
            rows[it->second].rhs = to_double(t[k + 1]);
          } else {
            rows[it->second].range = to_double(t[k + 1]);
            rows[it->second].has_range = true;
          }
        }
      } else if (section == "BOUNDS" && t.size() >= 2) {
        // "TYPE [set] column [value]": the set name is optional, so count the tokens.
        const std::string kind(t[0]);
        const bool valued = kind != "FR" && kind != "MI" && kind != "PL" && kind != "BV";
        const bool named = t.size() >= (valued ? 4u : 3u);
        const std::size_t j = column(t[named ? 2 : 1]);
        const double v = valued ? to_double(t[named ? 3 : 2]) : 0.0;
        if (kind == "UP" || kind == "UI") upper[j] = v;
        if (kind == "LO" || kind == "LI") lower[j] = v;
        if (kind == "FX") lower[j] = upper[j] = v;
        if (kind == "FR") lower[j] = -kInf, upper[j] = kInf;
        if (kind == "MI") lower[j] = -kInf;
        if (kind == "PL") upper[j] = kInf;
        if (kind == "BV") lower[j] = 0.0, upper[j] = 1.0, is_integer[j] = true;
        if (kind == "UI" || kind == "LI") is_integer[j] = true;
      }
    } catch (const std::exception&) {
      ok = false;
    }
  });
  if (!ok || rows.empty()) {
    out.reason = "the model is not in free MPS format";
    return out;
  }

  std::unordered_map<std::string, double> x;
  x.reserve(solution.columns.size());
  for (const auto& [name, value] : solution.columns) x.emplace(name, value);

  std::vector<Sum> activity(rows.size());
  std::vector<double> magnitude(rows.size(), 0.0);
  for (std::size_t j = 0; j < col_names.size(); ++j) {
    auto it = x.find(col_names[j]);
    if (it == x.end()) {
      ++out.missing;
      continue;
    }
    const double v = it->second;
    for (const Entry& e : col_entries[j]) {
      activity[e.row].add(e.value * v);
      magnitude[e.row] = std::max(magnitude[e.row], std::fabs(e.value * v));
    }
    const double lo = lower[j], up = upper[j];
    if (v < lo)
      out.bound_violation = std::max(out.bound_violation, (lo - v) / (1.0 + std::fabs(lo)));
    if (v > up)
      out.bound_violation = std::max(out.bound_violation, (v - up) / (1.0 + std::fabs(up)));
    if (is_integer[j]) {
      ++out.integers;
      out.integrality = std::max(out.integrality, std::fabs(v - std::round(v)));
    }
  }
  for (std::size_t i = 0; i < rows.size(); ++i) {
    const RowInfo& r = rows[i];
    if (r.type == 'N') continue;
    double lo = -kInf, up = kInf;
    const double b = r.rhs, R = std::fabs(r.range);
    if (r.type == 'E') {
      lo = up = b;
      if (r.has_range) (r.range > 0 ? up : lo) = r.range > 0 ? b + R : b - R;
    } else if (r.type == 'L') {
      up = b;
      if (r.has_range) lo = b - R;
    } else if (r.type == 'G') {
      lo = b;
      if (r.has_range) up = b + R;
    }
    const double a = activity[i].value();
    const double viol = std::max({lo - a, a - up, 0.0});
    const double scale = 1.0 + std::max({std::isfinite(lo) ? std::fabs(lo) : 0.0,
                                         std::isfinite(up) ? std::fabs(up) : 0.0, magnitude[i]});
    if (viol / scale > out.row_violation) {
      out.row_violation = viol / scale;
      out.row_name = row_names[i];
    }
  }
  if (!objective_row.empty()) {
    const std::size_t o = row_index[objective_row];
    // The RHS of the objective row is minus its constant term (MPS convention).
    out.objective = activity[o].value() - rows[o].rhs;
    if (solution.has_objective && !out.quadratic) {
      out.objective_diff =
          std::fabs(out.objective - solution.objective) / (1.0 + std::fabs(solution.objective));
    }
  }
  out.rows = static_cast<long long>(
      std::count_if(rows.begin(), rows.end(), [](const RowInfo& r) { return r.type != 'N'; }));
  out.cols = static_cast<long long>(col_names.size());
  out.available = true;
  return out;
}

std::string base64_decode(std::string_view text) {
  auto value = [](char c) -> int {
    if (c >= 'A' && c <= 'Z') return c - 'A';
    if (c >= 'a' && c <= 'z') return c - 'a' + 26;
    if (c >= '0' && c <= '9') return c - '0' + 52;
    if (c == '+') return 62;
    if (c == '/') return 63;
    return -1;
  };
  std::string out;
  out.reserve(text.size() * 3 / 4);
  unsigned buffer = 0;
  int bits = 0;
  for (const char c : text) {
    const int v = value(c);
    if (v < 0) continue;  // padding and whitespace
    buffer = (buffer << 6) | static_cast<unsigned>(v);
    bits += 6;
    if (bits >= 8) {
      bits -= 8;
      out += static_cast<char>((buffer >> bits) & 0xFF);
    }
  }
  return out;
}

}  // namespace studio
