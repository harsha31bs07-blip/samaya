#include <filesystem>
#include <fstream>
#include <string>

#include "samaya/io.hpp"
#include "test_framework.hpp"

using samaya::kInf;
using samaya::Model;
using samaya::ParseError;
using samaya::read_mps_from_string;
using samaya::VarType;

TEST(mps_reads_tiny_lp_file) {
  const Model m = samaya::read_mps(std::string(SAMAYA_TEST_DATA_DIR) + "/tiny_lp.mps");
  CHECK_EQ(m.name, "TINY_LP");
  CHECK(m.sense == samaya::ObjSense::kMaximize);
  REQUIRE(m.num_rows() == 2);
  REQUIRE(m.num_cols() == 2);
  CHECK_EQ(m.A.nnz(), 4);
  CHECK_EQ(m.obj[0], 3.0);
  CHECK_EQ(m.obj[1], 2.0);
  CHECK_EQ(m.row_lower[0], -kInf);
  CHECK_EQ(m.row_upper[0], 4.0);
  CHECK_EQ(m.row_upper[1], 6.0);
  CHECK_EQ(m.col_upper[0], 3.0);
  CHECK_EQ(m.col_upper[1], kInf);
  CHECK(m.problem_class() == samaya::ProblemClass::kLP);
  CHECK(m.validate().empty());
}

TEST(mps_missing_file_throws) {
  CHECK_THROWS(samaya::read_mps("/nonexistent/file.mps"), std::runtime_error);
}

TEST(mps_reads_utf8_path_with_spaces) {
  // A file name with Devanagari and a space. The path is passed as UTF-8; on Windows it must be
  // opened with the wide-character API (sources compile with /utf-8).
  const std::string name = "samaya \u092E\u0949\u0921\u0932 test.mps";
  const std::filesystem::path file = std::filesystem::temp_directory_path() /
                                     std::filesystem::path(std::u8string(name.begin(), name.end()));
  {
    std::ofstream out(file);
    out << "NAME UNI\nROWS\n N obj\n L c1\nCOLUMNS\n x obj 1 c1 1\n y obj 2 c1 1\nRHS\n rhs c1 4\n"
           "ENDATA\n";
  }
  const std::u8string u8 = file.u8string();
  const Model m = samaya::read_mps(std::string(u8.begin(), u8.end()));
  CHECK_EQ(m.num_rows(), 1);
  CHECK_EQ(m.num_cols(), 2);
  std::filesystem::remove(file);
}

TEST(mps_row_types_ranges_and_offset) {
  const Model m = read_mps_from_string(R"(NAME RANGES
ROWS
 N obj
 E e_pos
 E e_neg
 L l
 G g
 N unused
COLUMNS
 x obj 1 e_pos 1
 x e_neg 1 l 1
 x g 1 unused 7
RHS
 obj -2.5
 e_pos 10 e_neg 10
 l 10 g 10
RANGES
 rng e_pos 4 e_neg -4
 rng l 3 g -3
ENDATA
)");
  REQUIRE(m.num_rows() == 4);
  CHECK_EQ(m.obj_offset, 2.5);
  CHECK_EQ(m.row_lower[0], 10.0);  // E, R > 0: [rhs, rhs + |R|]
  CHECK_EQ(m.row_upper[0], 14.0);
  CHECK_EQ(m.row_lower[1], 6.0);   // E, R < 0: [rhs - |R|, rhs]
  CHECK_EQ(m.row_upper[1], 10.0);
  CHECK_EQ(m.row_lower[2], 7.0);   // L: [rhs - |R|, rhs]
  CHECK_EQ(m.row_upper[2], 10.0);
  CHECK_EQ(m.row_lower[3], 10.0);  // G: [rhs, rhs + |R|]
  CHECK_EQ(m.row_upper[3], 13.0);
  CHECK_EQ(m.A.nnz(), 4);          // Entry on the extra N row is dropped.
}

TEST(mps_bound_types) {
  const Model m = read_mps_from_string(R"(NAME BOUNDS
ROWS
 N obj
 L c
COLUMNS
 a c 1
 b c 1
 c1 c 1
 d c 1
 e c 1
 f c 1
 g c 1
 h c 1
 i c 1
RHS
 c 1
BOUNDS
 UP BND a 5
 LO BND b -2
 FX BND c1 3.5
 FR BND d
 MI BND e
 UP BND f -1
 BV BND g
 LI BND h 2
 UI BND h 9
 LO BND i 1
 UP BND i 1e30
ENDATA
)");
  REQUIRE(m.num_cols() == 9);
  CHECK_EQ(m.col_lower[0], 0.0);
  CHECK_EQ(m.col_upper[0], 5.0);
  CHECK_EQ(m.col_lower[1], -2.0);
  CHECK_EQ(m.col_lower[2], 3.5);
  CHECK_EQ(m.col_upper[2], 3.5);
  CHECK_EQ(m.col_lower[3], -kInf);
  CHECK_EQ(m.col_upper[3], kInf);
  CHECK_EQ(m.col_lower[4], -kInf);
  CHECK_EQ(m.col_upper[4], kInf);
  CHECK_EQ(m.col_lower[5], -kInf);  // Negative UP with default lower bound.
  CHECK_EQ(m.col_upper[5], -1.0);
  CHECK(m.col_type[6] == VarType::kInteger);
  CHECK_EQ(m.col_upper[6], 1.0);
  CHECK(m.col_type[7] == VarType::kInteger);
  CHECK_EQ(m.col_lower[7], 2.0);
  CHECK_EQ(m.col_upper[7], 9.0);
  CHECK_EQ(m.col_lower[8], 1.0);
  CHECK_EQ(m.col_upper[8], kInf);
  CHECK(m.col_type[0] == VarType::kContinuous);
}

TEST(mps_integer_markers) {
  const Model m = read_mps_from_string(R"(NAME MIP
ROWS
 N obj
 G c
COLUMNS
 x obj 1 c 1
 MARKER 'MARKER' 'INTORG'
 y obj 2 c 1
 z obj 3 c 1
 MARKER 'MARKER' 'INTEND'
 w obj 4 c 1
RHS
 RHS c 1
ENDATA
)");
  REQUIRE(m.num_cols() == 4);
  CHECK(m.col_type[0] == VarType::kContinuous);
  CHECK(m.col_type[1] == VarType::kInteger);
  CHECK(m.col_type[2] == VarType::kInteger);
  CHECK(m.col_type[3] == VarType::kContinuous);
  CHECK_EQ(m.col_upper[1], kInf);
  CHECK_EQ(m.num_integers(), 2);
  CHECK(m.problem_class() == samaya::ProblemClass::kMILP);
}

TEST(mps_quadobj_and_qmatrix) {
  const char* header = R"(NAME QP
ROWS
 N obj
 E c
COLUMNS
 x obj 1 c 1
 y obj 1 c 1
RHS
 c 1
)";
  const Model a = read_mps_from_string(std::string(header) + R"(QUADOBJ
 x x 2
 x y 1
 y y 4
ENDATA
)");
  const Model b = read_mps_from_string(std::string(header) + R"(QMATRIX
 x x 2
 x y 1
 y x 1
 y y 4
ENDATA
)");
  for (const Model* m : {&a, &b}) {
    REQUIRE(m->is_qp());
    CHECK(m->problem_class() == samaya::ProblemClass::kQP);
    CHECK_EQ(m->Q.nnz(), 3);
    // Lower triangle: column 0 holds (0,0)=2 and (1,0)=1, column 1 holds (1,1)=4.
    CHECK_EQ(m->Q.row_index()[1], 1);
    CHECK_EQ(m->Q.values()[1], 1.0);
    CHECK_EQ(m->Q.values()[2], 4.0);
    CHECK(m->validate().empty());
  }
}

TEST(mps_objsense_inline_and_crlf) {
  const Model m = read_mps_from_string(
      "NAME X\r\nOBJSENSE MAXIMIZE\r\nROWS\r\n N obj\r\nCOLUMNS\r\n x obj 1\r\nENDATA\r\n");
  CHECK(m.sense == samaya::ObjSense::kMaximize);
  CHECK_EQ(m.num_cols(), 1);
  CHECK_EQ(m.num_rows(), 0);
}

// A fixed-MPS data line: fields start in columns 2, 5, 15, 25, 40 and 50.
std::string fixed_line(const std::string& type, const std::string& f1, const std::string& f2 = "",
                       const std::string& f3 = "", const std::string& f4 = "",
                       const std::string& f5 = "") {
  std::string line(61, ' ');
  const std::string* fields[] = {&type, &f1, &f2, &f3, &f4, &f5};
  const std::size_t start[] = {1, 4, 14, 24, 39, 49};
  for (std::size_t f = 0; f < 6; ++f) line.replace(start[f], fields[f]->size(), *fields[f]);
  line.erase(line.find_last_not_of(' ') + 1);
  return line + "\n";
}

TEST(mps_fixed_format_names_with_spaces) {
  const std::string text =
      "NAME          FIXED\nROWS\n" + fixed_line("N", "COST") + fixed_line("L", "LIM 1") +
      fixed_line("G", "MY ROW") + "COLUMNS\n" + fixed_line("", "MARKER", "'MARKER'", "", "'INTORG'") +
      fixed_line("", "X ONE", "COST", "1.0", "LIM 1", "1.0") +
      fixed_line("", "X ONE", "MY ROW", "1.0") +
      fixed_line("", "MARKER", "'MARKER'", "", "'INTEND'") +
      fixed_line("", "X TWO", "COST", "2.", "LIM 1", "1.") +
      fixed_line("", "X TWO", "MY ROW", "-1.") + "RHS\n" +
      fixed_line("", "RHS 1", "LIM 1", "4.", "MY ROW", "-2.") + "RANGES\n" +
      fixed_line("", "RNG 1", "LIM 1", "3.") + "BOUNDS\n" + fixed_line("UP", "BND 1", "X TWO", "3.") +
      "ENDATA\n";
  const Model m = read_mps_from_string(text);
  CHECK_EQ(m.num_rows(), 2);
  CHECK_EQ(m.num_cols(), 2);
  CHECK(m.col_names[0] == "X ONE");
  CHECK(m.row_names[1] == "MY ROW");
  CHECK(m.col_type[0] == samaya::VarType::kInteger);
  CHECK(m.col_type[1] == samaya::VarType::kContinuous);
  CHECK_EQ(m.obj[1], 2.0);
  CHECK_EQ(m.row_lower[0], 1.0);
  CHECK_EQ(m.row_upper[0], 4.0);
  CHECK_EQ(m.row_lower[1], -2.0);
  CHECK_EQ(m.col_upper[1], 3.0);
  CHECK_EQ(m.A.nnz(), 4);
}

TEST(mps_errors_report_line_numbers) {
  try {
    read_mps_from_string("NAME X\nROWS\n N obj\nCOLUMNS\n x nosuchrow 1\nENDATA\n");
    CHECK(false);
  } catch (const ParseError& e) {
    CHECK_EQ(e.line(), 5);
    CHECK(std::string(e.what()).find("nosuchrow") != std::string::npos);
  }
  CHECK_THROWS(read_mps_from_string("NAME X\nROWS\n Q r\nENDATA\n"), ParseError);
  CHECK_THROWS(read_mps_from_string("NAME X\nROWS\n N obj\n L r\n L r\nENDATA\n"), ParseError);
  CHECK_THROWS(read_mps_from_string("NAME X\nROWS\n N o\nCOLUMNS\n x o abc\nENDATA\n"), ParseError);
  CHECK_THROWS(read_mps_from_string("NAME X\nBOGUS\nENDATA\n"), ParseError);
  CHECK_THROWS(read_mps_from_string(
                   "NAME X\nROWS\n N o\nCOLUMNS\n x o 1\nBOUNDS\n SC B x 1\nENDATA\n"),
               ParseError);
  // Crossed bounds make a model infeasible, not malformed, so they parse.
  const Model crossed = read_mps_from_string(
      "NAME X\nROWS\n N o\nCOLUMNS\n x o 1\nBOUNDS\n LO B x 2\n UP B x 1\nENDATA\n");
  CHECK(!crossed.crossed_bounds().empty());
  // Infinite bounds on the wrong side are caught by model validation.
  CHECK_THROWS(read_mps_from_string(
                   "NAME X\nROWS\n N o\nCOLUMNS\n x o 1\nBOUNDS\n LO B x 1e30\nENDATA\n"),
               ParseError);
}
