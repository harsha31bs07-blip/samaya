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
)
target_link_libraries(samaya_tests PRIVATE samaya::samaya samaya_options)
# Tests exercise internal components (factorization, simplex) directly.
target_include_directories(samaya_tests PRIVATE ${PROJECT_SOURCE_DIR}/src)
target_compile_definitions(samaya_tests PRIVATE
  SAMAYA_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/instances")

if(SAMAYA_CUDA)
  # GPU tests need the CUDA runtime.
  target_sources(samaya_tests PRIVATE test_gpu.cpp)
  target_link_libraries(samaya_tests PRIVATE CUDA::cudart)
  target_compile_definitions(samaya_tests PRIVATE SAMAYA_HAVE_CUDA)
endif()

add_test(NAME unit_tests COMMAND samaya_tests)

if(SAMAYA_BUILD_CLI)
  add_test(NAME cli_stats COMMAND samaya_cli --stats ${CMAKE_CURRENT_SOURCE_DIR}/instances/tiny_lp.mps)
  set_tests_properties(cli_stats PROPERTIES PASS_REGULAR_EXPRESSION "rows 2 .*cols 2")
  add_test(NAME cli_json COMMAND samaya_cli --json --log-level 0
    ${CMAKE_CURRENT_SOURCE_DIR}/instances/tiny_lp.mps)
  set_tests_properties(cli_json PROPERTIES PASS_REGULAR_EXPRESSION
    "\"status\":\"optimal\",\"verified\":true,\"objective\":11[,.]")
  add_test(NAME cli_missing_file COMMAND samaya_cli /nonexistent.mps)
  set_tests_properties(cli_missing_file PROPERTIES WILL_FAIL TRUE)
endif()
