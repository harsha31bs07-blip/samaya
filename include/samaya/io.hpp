#pragma once

#include <stdexcept>
#include <string>
#include <string_view>

#include "samaya/model.hpp"

namespace samaya {

class ParseError : public std::runtime_error {
 public:
  ParseError(const std::string& message, long long line)
      : std::runtime_error(line > 0 ? "line " + std::to_string(line) + ": " + message : message),
        line_(line) {}

  long long line() const { return line_; }

 private:
  long long line_;
};

// Reads a model in (free or fixed) MPS format, including the QPS extensions QUADOBJ / QMATRIX.
// Names may contain spaces only in fixed MPS, which is tried when free MPS fails. Throws
// ParseError on malformed input and std::runtime_error if the file cannot be opened. The path is
// UTF-8 on every platform (on Windows it is opened with the wide-character API).
Model read_mps(const std::string& path);
Model read_mps_from_string(std::string_view text);

}  // namespace samaya
