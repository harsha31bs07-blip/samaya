add_executable(samaya_tests
  test_basis_factor.cpp
  test_ldl.cpp
  test_lp_known.cpp
  test_lp_random.cpp
  test_mip.cpp
  test_main.cpp
  test_model_and_solver.cpp
  test_mps_reader.cpp
  test_pdlp.cpp
  test_presolve.cpp
  test_qp.cpp
  test_scaling.cpp
  test_sparse_matrix.cpp
  test_verify.cpp
  test_xreal.cpp
  test_studio_core.cpp
  ${PROJECT_SOURCE_DIR}/apps/studio/studio_core.cpp
)
target_link_libraries(samaya_tests PRIVATE samaya::samaya samaya_options)
# Tests exercise internal components (factorization, simplex) directly.
target_include_directories(samaya_tests PRIVATE ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/apps/studio)
target_compile_definitions(samaya_tests PRIVATE
  SAMAYA_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/instances")

if(SAMAYA_CUDA)
  # GPU tests need the CUDA runtime.
  target_sources(samaya_tests PRIVATE test_gpu.cpp)
  if(WIN32)
    target_link_libraries(samaya_tests PRIVATE CUDA::cudart_static)
  else()
    target_link_libraries(samaya_tests PRIVATE CUDA::cudart)
  endif()
  target_compile_definitions(samaya_tests PRIVATE SAMAYA_HAVE_CUDA
                             SAMAYA_CUDA_MIN_ARCH=${SAMAYA_CUDA_MIN_ARCH})
endif()

add_test(NAME unit_tests COMMAND samaya_tests)

if(SAMAYA_BUILD_CLI)
  add_test(NAME cli_stats COMMAND samaya_cli --stats ${CMAKE_CURRENT_SOURCE_DIR}/instances/tiny_lp.mps)
  set_tests_properties(cli_stats PROPERTIES PASS_REGULAR_EXPRESSION "rows 2 .*cols 2")
  add_test(NAME cli_json COMMAND samaya_cli --json --log-level 0
    ${CMAKE_CURRENT_SOURCE_DIR}/instances/tiny_lp.mps)
  set_tests_properties(cli_json PROPERTIES PASS_REGULAR_EXPRESSION
    "\"status\":\"optimal\",\"verified\":true,\"objective\":11[,.]")
  add_test(NAME cli_gpu_info COMMAND samaya_cli --gpu-info)
  if(SAMAYA_CUDA)
    set_tests_properties(cli_gpu_info PROPERTIES PASS_REGULAR_EXPRESSION "^\\{\"cuda\":true,\"gpu\":")
  else()
    set_tests_properties(cli_gpu_info PROPERTIES PASS_REGULAR_EXPRESSION "^\\{\"cuda\":false,\"gpu\":null\\}")
  endif()
  add_test(NAME cli_missing_file COMMAND samaya_cli /nonexistent.mps)
  set_tests_properties(cli_missing_file PROPERTIES WILL_FAIL TRUE)
  # The run report (cases/report.py): both formats on a generated MRPL case.
  find_package(Python3 COMPONENTS Interpreter)
  if(Python3_Interpreter_FOUND)
    add_test(NAME report_formats COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/report_smoke.py
      $<TARGET_FILE:samaya_cli> ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/report_smoke)
    set_tests_properties(report_formats PROPERTIES PASS_REGULAR_EXPRESSION "report smoke test passed")
  endif()
endif()
