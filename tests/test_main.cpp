#include <cstring>
#include <exception>
#include <iostream>

#include "test_framework.hpp"

#if defined(_MSC_VER) && defined(_DEBUG)
#include <crtdbg.h>
#include <cstdlib>
#endif

int main(int argc, char** argv) {
#if defined(_MSC_VER) && defined(_DEBUG)
  // MSVC's debug library reports failed assertions in a dialog that waits for a click; send the
  // reports to stderr instead, so a failure fails the run (unattended, as in CI).
  for (const int kind : {_CRT_WARN, _CRT_ERROR, _CRT_ASSERT}) {
    _CrtSetReportMode(kind, _CRTDBG_MODE_FILE);
    _CrtSetReportFile(kind, _CRTDBG_FILE_STDERR);
  }
  _set_abort_behavior(0, _WRITE_ABORT_MSG | _CALL_REPORTFAULT);
#endif
  using namespace samaya::test;
  const char* filter = argc > 1 ? argv[1] : nullptr;
  int ran = 0;
  int failed_cases = 0;
  for (const Case& c : registry()) {
    if (filter != nullptr && std::strstr(c.name, filter) == nullptr) continue;
    const int before = failures();
    try {
      c.body();
    } catch (const AbortCase&) {
    } catch (const std::exception& e) {
      report(c.name, 0, std::string("unexpected exception: ") + e.what());
    }
    ++ran;
    if (failures() != before) {
      ++failed_cases;
      std::cerr << "[FAIL] " << c.name << "\n";
    } else {
      // Flushed, so that a crash in the next case (an abort in a debug library) names the right one.
      std::cout << "[ OK ] " << c.name << std::endl;
    }
  }
  std::cout << ran - failed_cases << "/" << ran << " test cases passed\n";
  return failed_cases == 0 && ran > 0 ? 0 : 1;
}
