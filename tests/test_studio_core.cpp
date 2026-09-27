// samaya Studio's portable core (apps/studio/studio_core.*): the JSON used by the page, the
// solution-file reader and the independent re-check.

#include <cmath>
#include <string>

#include "studio_core.hpp"
#include "test_framework.hpp"

namespace {

// max 3x + 2y + 0.5  s.t.  x + y <= 4 (c1),  x + 3y >= 2 (c2),  x <= 3,  y integer.
const char* kModel =
    "NAME TINY\n"
    "OBJSENSE\n    MAX\n"
    "ROWS\n N obj\n L c1\n G c2\n"
    "COLUMNS\n x obj 3 c1 1\n x c2 1\n"
    " M1 'MARKER' 'INTORG'\n y obj 2 c1 1\n y c2 3\n M2 'MARKER' 'INTEND'\n"
    "RHS\n rhs c1 4 c2 2\n rhs obj -0.5\n"
    "BOUNDS\n UP BND x 3\n UP y 10\n"
    "ENDATA\n";

std::string solution_text(double x, double y, double objective) {
  return "status optimal\nverified yes\nobjective " + std::to_string(objective) +
         "\n\ncolumns 2\n# name value\nx " + std::to_string(x) + "\ny " + std::to_string(y) +
         "\n\nrows 2\n# name activity\nc1 " + std::to_string(x + y) + "\nc2 " +
         std::to_string(x + 3 * y) + "\n";
}

}  // namespace

TEST(studio_json_round_trip_and_escapes) {
  const studio::Json j = studio::parse_json(
      R"({"type":"solve","id":7,"path":"C:\\Users\\H B\\m\u00e9.mps","options":{"presolve":true,"gap":1e-4},"list":[1,2]})");
  CHECK_EQ(j.str("type"), std::string("solve"));
  CHECK_EQ(j.num("id", 0), 7.0);
  CHECK_EQ(j.str("path"), std::string("C:\\Users\\H B\\m\xC3\xA9.mps"));
  REQUIRE(j.get("options") != nullptr);
  CHECK(j.get("options")->flag("presolve", false));
  CHECK_EQ(j.get("options")->num("gap", 0), 1e-4);
  CHECK_EQ(j.get("list")->array.size(), std::size_t{2});
  CHECK_EQ(studio::json_string("a\"b\\c\nd"), std::string("\"a\\\"b\\\\c\\nd\""));
  CHECK_EQ(studio::json_number(INFINITY), std::string("null"));
  CHECK_THROWS(studio::parse_json("{\"a\":}"), std::runtime_error);
}

TEST(studio_reads_solution_files) {
  const studio::SolutionFile s = studio::read_solution(
      "status optimal\nverified yes\nobjective 12.5\n\ncolumns 2\n# name value reduced_cost\n"
      "x 3 0\ny 1 -0.5\n\nrows 1\n# name activity dual\nc1 4 2.5\n");
  CHECK_EQ(s.status, std::string("optimal"));
  CHECK(s.verified);
  CHECK_EQ(s.objective, 12.5);
  REQUIRE(s.columns.size() == 2);
  CHECK_EQ(s.columns[1].first, std::string("y"));
  REQUIRE(s.rows.size() == 1);
  CHECK_EQ(s.rows[0].dual, 2.5);
}

TEST(studio_recheck_accepts_the_optimum_and_catches_violations) {
  // The optimum: x = 3, y = 1, objective 3*3 + 2*1 + 0.5 = 11.5.
  const studio::Recheck good =
      studio::recheck(kModel, studio::read_solution(solution_text(3, 1, 11.5)));
  REQUIRE(good.available);
  CHECK(good.maximize);
  CHECK_EQ(good.rows, 2);
  CHECK_EQ(good.integers, 1);
  CHECK_EQ(good.row_violation, 0.0);
  CHECK_EQ(good.bound_violation, 0.0);
  CHECK_EQ(good.integrality, 0.0);
  CHECK_NEAR(good.objective, 11.5, 1e-12);
  CHECK(good.objective_diff < 1e-12);
  // x = 3.5 breaks its bound and row c1 (3.5 + 1 > 4); y = 1.25 is fractional; the objective
  // given does not match the point.
  const studio::Recheck bad =
      studio::recheck(kModel, studio::read_solution(solution_text(3.5, 1.25, 11.5)));
  REQUIRE(bad.available);
  CHECK(bad.bound_violation > 0.1);
  CHECK(bad.row_violation > 0.1);
  CHECK_EQ(bad.row_name, std::string("c1"));
  CHECK_NEAR(bad.integrality, 0.25, 1e-12);
  CHECK(bad.objective_diff > 0.1);
  // No solution (an infeasible model): nothing to re-check.
  CHECK(!studio::recheck(kModel, studio::read_solution("status infeasible\nverified yes\n"))
             .available);
}

TEST(studio_base64_decodes_binary) {
  CHECK_EQ(studio::base64_decode("UEsDBA=="), std::string("PK\x03\x04", 4));
  CHECK_EQ(studio::base64_decode("aGVsbG8gd29ybGQ="), std::string("hello world"));
}
