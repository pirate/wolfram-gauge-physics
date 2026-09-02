cmake_minimum_required(VERSION 3.20)
project(wolfram_gauge_physics LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(WGPHYSICS_BUILD_REWRITE_HARNESS
  "Build the executable backed by WolframInstitute/HypergraphRewritingEngine" ON)

add_library(wgphysics INTERFACE)
target_include_directories(wgphysics INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/src")

add_executable(fiber_demo examples/fiber_rewrite.cpp)
target_link_libraries(fiber_demo PRIVATE wgphysics)
target_compile_options(fiber_demo PRIVATE -Wall -Wextra -Wpedantic)

add_executable(research_demo examples/research_pipeline.cpp)
target_link_libraries(research_demo PRIVATE wgphysics)
target_compile_options(research_demo PRIVATE -Wall -Wextra -Wpedantic)

add_executable(wgphysics_census src/census.cpp)
target_link_libraries(wgphysics_census PRIVATE wgphysics)
target_compile_options(wgphysics_census PRIVATE -Wall -Wextra -Wpedantic)

enable_testing()
add_executable(infragauge_test tests/infragauge.cpp)
target_link_libraries(infragauge_test PRIVATE wgphysics)
target_compile_options(infragauge_test PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME infragauge COMMAND infragauge_test)
set_tests_properties(infragauge PROPERTIES LABELS wgphysics)

add_executable(research_test tests/research.cpp)
target_link_libraries(research_test PRIVATE wgphysics)
target_compile_options(research_test PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME research COMMAND research_test)
set_tests_properties(research PROPERTIES LABELS wgphysics)

if(WGPHYSICS_BUILD_REWRITE_HARNESS)
  include(FetchContent)
  if(POLICY CMP0169)
    cmake_policy(SET CMP0169 OLD)
  endif()
  if(POLICY CMP0097)
    cmake_policy(SET CMP0097 NEW)
  endif()

  set(BUILD_TESTS OFF CACHE BOOL "" FORCE)
  set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
  set(BUILD_BENCHMARKS OFF CACHE BOOL "" FORCE)
  set(BUILD_GPU OFF CACHE BOOL "" FORCE)
  set(BUILD_VISUALIZATION OFF CACHE BOOL "" FORCE)
  set(BUILD_WOLFRAM_LANGUAGE_PACLET OFF CACHE BOOL "" FORCE)
  set(HG_ENGINE_STATS OFF CACHE BOOL "" FORCE)

  FetchContent_Declare(
    hypergraph_rewriting_engine
    GIT_REPOSITORY https://github.com/WolframInstitute/HypergraphRewritingEngine.git
    GIT_TAG 03fe60ddf338983060b6bb4b23e8b4b5d7ae7337
    GIT_SHALLOW TRUE
    # The engine sources do not require its documentation-tool submodules.
    # Skipping them also avoids nested example remotes that use SSH-only URLs.
    GIT_SUBMODULES ""
    GIT_SUBMODULES_RECURSE FALSE
  )

  # The pinned upstream top-level file currently assumes CMAKE_SOURCE_DIR. Patch
  # only FetchContent's disposable copy so it can be embedded as a subproject.
  FetchContent_GetProperties(hypergraph_rewriting_engine)
  if(NOT hypergraph_rewriting_engine_POPULATED)
    FetchContent_Populate(hypergraph_rewriting_engine)
    set(_engine_cmake "${hypergraph_rewriting_engine_SOURCE_DIR}/CMakeLists.txt")
    file(READ "${_engine_cmake}" _engine_contents)
    string(REPLACE "\${CMAKE_SOURCE_DIR}" "\${CMAKE_CURRENT_SOURCE_DIR}"
      _engine_contents "${_engine_contents}")
    string(REPLACE "\${CMAKE_BINARY_DIR}" "\${CMAKE_CURRENT_BINARY_DIR}"
      _engine_contents "${_engine_contents}")
    file(WRITE "${_engine_cmake}" "${_engine_contents}")
    add_subdirectory(
      "${hypergraph_rewriting_engine_SOURCE_DIR}"
      "${hypergraph_rewriting_engine_BINARY_DIR}"
      EXCLUDE_FROM_ALL
    )
  endif()

  add_executable(wgphysics_evolve src/main.cpp)
  target_link_libraries(wgphysics_evolve PRIVATE hypergraph)
  target_include_directories(wgphysics_evolve PRIVATE
    "${hypergraph_rewriting_engine_SOURCE_DIR}/common/include"
    "${hypergraph_rewriting_engine_BINARY_DIR}/generated"
  )
  target_compile_options(wgphysics_evolve PRIVATE -Wall -Wextra -Wpedantic)
endif()
