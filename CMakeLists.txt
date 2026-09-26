cmake_minimum_required(VERSION 3.22)

project(samaya
  VERSION 0.1.0
  DESCRIPTION "Sovereign LP / MILP / QP optimization solver core"
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(GNUInstallDirs)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(SAMAYA_BUILD_TESTS "Build unit tests" ON)
option(SAMAYA_BUILD_CLI "Build the samaya command-line tool" ON)
option(SAMAYA_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)
option(SAMAYA_SANITIZE "Build with AddressSanitizer + UndefinedBehaviorSanitizer" OFF)
option(SAMAYA_CUDA "Enable CUDA kernels (PDLP, SpMV); lands in Phase 3" OFF)

if(SAMAYA_CUDA)
  enable_language(CUDA)
  find_package(CUDAToolkit REQUIRED)
endif()

# ---------------------------------------------------------------------------
# Compiler flags shared by every target.
# ---------------------------------------------------------------------------
add_library(samaya_options INTERFACE)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(samaya_options INTERFACE -Wall -Wextra -Wpedantic -Wshadow)
  if(SAMAYA_WARNINGS_AS_ERRORS)
    target_compile_options(samaya_options INTERFACE -Werror)
  endif()
  if(SAMAYA_SANITIZE)
    target_compile_options(samaya_options INTERFACE
      -fsanitize=address,undefined -fno-omit-frame-pointer -fno-sanitize-recover=all)
    target_link_options(samaya_options INTERFACE -fsanitize=address,undefined)
  endif()
elseif(MSVC)
  target_compile_options(samaya_options INTERFACE /W4 /permissive-)
  if(SAMAYA_WARNINGS_AS_ERRORS)
    target_compile_options(samaya_options INTERFACE /WX)
  endif()
endif()

# ---------------------------------------------------------------------------
# Core library.
# ---------------------------------------------------------------------------
add_library(samaya
  src/core/log.cpp
  src/core/model.cpp
  src/core/model_stats.cpp
  src/core/solver.cpp
  src/core/status.cpp
  src/core/c_api.cpp
  src/io/mps_reader.cpp
  src/linalg/basis_factor.cpp
  src/lp/lp_solver.cpp
  src/lp/simplex.cpp
  src/mip/branch_and_bound.cpp
  src/mip/cuts.cpp
  src/mip/heuristics.cpp
  src/mip/probing.cpp
  src/mip/parallel.cpp
  src/mip/conflicts.cpp
  src/mip/feasibility_jump.cpp
  src/presolve/presolve.cpp
  src/linalg/scaling.cpp
  src/linalg/sparse_matrix.cpp
  src/verify/verify.cpp
)
add_library(samaya::samaya ALIAS samaya)
target_include_directories(samaya
  PUBLIC  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
          $<INSTALL_INTERFACE:include>
  PRIVATE ${PROJECT_SOURCE_DIR}/src)
target_link_libraries(samaya PRIVATE samaya_options)
target_compile_definitions(samaya PRIVATE SAMAYA_VERSION_STRING="${PROJECT_VERSION}")

find_package(Threads REQUIRED)
target_link_libraries(samaya PUBLIC Threads::Threads)

if(SAMAYA_BUILD_CLI)
  add_subdirectory(apps/cli)
endif()

if(SAMAYA_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

install(TARGETS samaya EXPORT samayaTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
