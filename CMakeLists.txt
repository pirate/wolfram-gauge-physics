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

add_executable(cell_dynamics_demo examples/cell_dynamics.cpp)
target_link_libraries(cell_dynamics_demo PRIVATE wgphysics)
target_compile_options(cell_dynamics_demo PRIVATE -Wall -Wextra -Wpedantic)

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

add_executable(observables_test tests/observables.cpp)
target_link_libraries(observables_test PRIVATE wgphysics)
target_compile_options(observables_test PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME observables COMMAND observables_test)
set_tests_properties(observables PROPERTIES LABELS wgphysics)

add_executable(gauge_interactions_test tests/gauge_interactions.cpp)
target_link_libraries(gauge_interactions_test PRIVATE wgphysics)
target_compile_options(gauge_interactions_test PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME gauge_interactions COMMAND gauge_interactions_test)
set_tests_properties(gauge_interactions PROPERTIES LABELS wgphysics)

add_executable(cell_chain_test tests/cell_chain.cpp)
target_link_libraries(cell_chain_test PRIVATE wgphysics)
target_compile_options(cell_chain_test PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME cell_chain COMMAND cell_chain_test)
set_tests_properties(cell_chain PROPERTIES LABELS wgphysics)

add_executable(mesh_dynamics_test tests/mesh_dynamics.cpp)
target_link_libraries(mesh_dynamics_test PRIVATE wgphysics)
target_compile_options(mesh_dynamics_test PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME mesh_dynamics COMMAND mesh_dynamics_test)
set_tests_properties(mesh_dynamics PROPERTIES LABELS wgphysics)

add_executable(shared_edge_test tests/shared_edge.cpp)
target_link_libraries(shared_edge_test PRIVATE wgphysics)
target_compile_options(shared_edge_test PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME shared_edge COMMAND shared_edge_test)
set_tests_properties(shared_edge PROPERTIES LABELS wgphysics)

add_executable(three_face_test tests/three_face.cpp)
target_link_libraries(three_face_test PRIVATE wgphysics)
target_compile_options(three_face_test PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME three_face COMMAND three_face_test)
set_tests_properties(three_face PROPERTIES LABELS wgphysics)

add_executable(mixed_bank_test tests/mixed_bank.cpp)
target_link_libraries(mixed_bank_test PRIVATE wgphysics)
target_compile_options(mixed_bank_test PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME mixed_bank COMMAND mixed_bank_test)
set_tests_properties(mixed_bank PROPERTIES LABELS wgphysics)

add_executable(equivariant_pair_search_test tests/equivariant_pair_search.cpp)
target_link_libraries(equivariant_pair_search_test PRIVATE wgphysics)
target_compile_options(equivariant_pair_search_test PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME equivariant_pair_search COMMAND equivariant_pair_search_test)
set_tests_properties(equivariant_pair_search PROPERTIES LABELS wgphysics)

add_executable(wgphysics_involution_census src/involution_census.cpp)
target_link_libraries(wgphysics_involution_census PRIVATE wgphysics)
target_compile_options(wgphysics_involution_census PRIVATE -Wall -Wextra -Wpedantic)

add_executable(wgphysics_involution_braid_search src/involution_braid_search.cpp)
target_link_libraries(wgphysics_involution_braid_search PRIVATE wgphysics)
target_compile_options(wgphysics_involution_braid_search PRIVATE -Wall -Wextra -Wpedantic)

add_executable(wgphysics_diffusion src/diffusion.cpp)
target_link_libraries(wgphysics_diffusion PRIVATE wgphysics)
target_compile_options(wgphysics_diffusion PRIVATE -Wall -Wextra -Wpedantic)

add_executable(wgphysics_mesh_experiments src/mesh_experiments.cpp)
target_link_libraries(wgphysics_mesh_experiments PRIVATE wgphysics)
target_compile_options(wgphysics_mesh_experiments PRIVATE -Wall -Wextra -Wpedantic)

add_executable(wgphysics_mesh_orbits src/mesh_orbits.cpp)
target_link_libraries(wgphysics_mesh_orbits PRIVATE wgphysics)
target_compile_options(wgphysics_mesh_orbits PRIVATE -Wall -Wextra -Wpedantic)

add_executable(wgphysics_triple_census src/triple_census.cpp)
target_link_libraries(wgphysics_triple_census PRIVATE wgphysics)
target_compile_options(wgphysics_triple_census PRIVATE -Wall -Wextra -Wpedantic)

add_executable(wgphysics_three_face_experiments src/three_face_experiments.cpp)
target_link_libraries(wgphysics_three_face_experiments PRIVATE wgphysics)
target_compile_options(wgphysics_three_face_experiments PRIVATE -Wall -Wextra -Wpedantic)

add_executable(wgphysics_three_face_replay src/three_face_replay.cpp)
target_link_libraries(wgphysics_three_face_replay PRIVATE wgphysics)
target_compile_options(wgphysics_three_face_replay PRIVATE -Wall -Wextra -Wpedantic)

add_executable(wgphysics_mixed_bank_experiments src/mixed_bank_experiments.cpp)
target_link_libraries(wgphysics_mixed_bank_experiments PRIVATE wgphysics)
target_compile_options(wgphysics_mixed_bank_experiments PRIVATE -Wall -Wextra -Wpedantic)

add_executable(wgphysics_chain_experiments src/chain_experiments.cpp)
target_link_libraries(wgphysics_chain_experiments PRIVATE wgphysics)
target_compile_options(wgphysics_chain_experiments PRIVATE -Wall -Wextra -Wpedantic)

add_executable(wgphysics_pair_census src/pair_census.cpp)
target_link_libraries(wgphysics_pair_census PRIVATE wgphysics)
target_compile_options(wgphysics_pair_census PRIVATE -Wall -Wextra -Wpedantic)

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
  target_link_libraries(wgphysics_evolve PRIVATE hypergraph wgphysics)
  target_include_directories(wgphysics_evolve PRIVATE
    "${hypergraph_rewriting_engine_SOURCE_DIR}/common/include"
    "${hypergraph_rewriting_engine_BINARY_DIR}/generated"
  )
  target_compile_options(wgphysics_evolve PRIVATE -Wall -Wextra -Wpedantic)

  add_executable(wgphysics_product_evolve src/product_main.cpp)
  target_link_libraries(wgphysics_product_evolve PRIVATE hypergraph wgphysics)
  target_include_directories(wgphysics_product_evolve PRIVATE
    "${hypergraph_rewriting_engine_SOURCE_DIR}/common/include"
    "${hypergraph_rewriting_engine_BINARY_DIR}/generated"
  )
  target_compile_options(wgphysics_product_evolve PRIVATE -Wall -Wextra -Wpedantic)

  add_executable(product_evolution_test tests/product_evolution.cpp)
  target_link_libraries(product_evolution_test PRIVATE hypergraph wgphysics)
  target_include_directories(product_evolution_test PRIVATE
    "${hypergraph_rewriting_engine_SOURCE_DIR}/common/include"
    "${hypergraph_rewriting_engine_BINARY_DIR}/generated"
  )
  target_compile_options(product_evolution_test PRIVATE -Wall -Wextra -Wpedantic)
  add_test(NAME product_evolution COMMAND product_evolution_test)
  set_tests_properties(product_evolution PROPERTIES LABELS wgphysics)

  add_executable(wgphysics_product_census src/product_census.cpp)
  target_link_libraries(wgphysics_product_census PRIVATE hypergraph wgphysics)
  target_include_directories(wgphysics_product_census PRIVATE
    "${hypergraph_rewriting_engine_SOURCE_DIR}/common/include"
    "${hypergraph_rewriting_engine_BINARY_DIR}/generated"
  )
  target_compile_options(wgphysics_product_census PRIVATE -Wall -Wextra -Wpedantic)

  add_executable(wgphysics_dynamics_census src/dynamics_census.cpp)
  target_link_libraries(wgphysics_dynamics_census PRIVATE hypergraph wgphysics)
  target_include_directories(wgphysics_dynamics_census PRIVATE
    "${hypergraph_rewriting_engine_SOURCE_DIR}/common/include"
    "${hypergraph_rewriting_engine_BINARY_DIR}/generated"
  )
  target_compile_options(wgphysics_dynamics_census PRIVATE -Wall -Wextra -Wpedantic)
endif()
