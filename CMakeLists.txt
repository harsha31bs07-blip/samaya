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
option(SAMAYA_BUILD_STUDIO "Build samaya Studio, the Windows desktop app (WebView2)" ${WIN32})
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
  # C++ sources only: nvcc passes host flags to its own generated code, which -Wpedantic rejects.
  target_compile_options(samaya_options INTERFACE
    "$<$<COMPILE_LANGUAGE:CXX>:-Wall;-Wextra;-Wpedantic;-Wshadow>")
  if(SAMAYA_WARNINGS_AS_ERRORS)
    target_compile_options(samaya_options INTERFACE "$<$<COMPILE_LANGUAGE:CXX>:-Werror>")
  endif()
  if(SAMAYA_SANITIZE)
    target_compile_options(samaya_options INTERFACE
      "$<$<COMPILE_LANGUAGE:CXX>:-fsanitize=address,undefined;-fno-omit-frame-pointer;-fno-sanitize-recover=all>")
    target_link_options(samaya_options INTERFACE -fsanitize=address,undefined)
  endif()
  # Exact IEEE arithmetic: Xreal (src/core/xreal.hpp) needs a * b + c rounded twice, never fused
  # into an FMA, or its error terms are wrong. Off everywhere, so builds with -march=native agree.
  target_compile_options(samaya_options INTERFACE "$<$<COMPILE_LANGUAGE:CXX>:-ffp-contract=off>")
elseif(MSVC)
  # /fp:precise: no FMA contraction and no reassociation, as Xreal needs (see above).
  target_compile_options(samaya_options INTERFACE
    "$<$<COMPILE_LANGUAGE:CXX>:/W4;/permissive-;/utf-8;/fp:precise;/Zc:__cplusplus>")
  if(SAMAYA_WARNINGS_AS_ERRORS)
    target_compile_options(samaya_options INTERFACE "$<$<COMPILE_LANGUAGE:CXX>:/WX>")
  endif()
  if(SAMAYA_SANITIZE)
    # MSVC has AddressSanitizer but no UBSan (UBSan runs in the Linux asan preset).
    target_compile_options(samaya_options INTERFACE "$<$<COMPILE_LANGUAGE:CXX>:/fsanitize=address;/Zi>")
    target_link_options(samaya_options INTERFACE /INCREMENTAL:NO /DEBUG)
    # The /RTC run-time checks cannot be combined with ASan.
    string(REGEX REPLACE "/RTC[1csu]+" "" CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
  endif()
  # Windows gives a program 1 MB of stack (Linux 8 MB); std::thread workers inherit this size.
  target_link_options(samaya_options INTERFACE /STACK:8388608)
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
  src/linalg/ldl.cpp
  src/linalg/ordering.cpp
  src/lp/ipm.cpp
  src/lp/lp_solver.cpp
  src/lp/pdlp.cpp
  src/lp/simplex.cpp
  src/mip/branch_and_bound.cpp
  src/mip/cuts.cpp
  src/mip/heuristics.cpp
  src/mip/probing.cpp
  src/mip/parallel.cpp
  src/mip/conflicts.cpp
  src/mip/feasibility_jump.cpp
  src/presolve/presolve.cpp
  src/qp/miqp.cpp
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

if(SAMAYA_CUDA)
  # GPU kernels for PDLP (WP2): our own SpMV and vector kernels, CUDA runtime only.
  target_sources(samaya PRIVATE
    src/gpu/device.cu
    src/gpu/spmv.cu
    src/gpu/vector_ops.cu
    src/gpu/pdhg_core.cu
  )
  target_link_libraries(samaya PRIVATE CUDA::cudart)
  target_compile_definitions(samaya PRIVATE SAMAYA_HAVE_CUDA)
endif()

if(SAMAYA_BUILD_CLI)
  add_subdirectory(apps/cli)
  if(WIN32 AND SAMAYA_BUILD_STUDIO)
    add_subdirectory(apps/studio)
  endif()
endif()

if(SAMAYA_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

install(TARGETS samaya EXPORT samayaTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
