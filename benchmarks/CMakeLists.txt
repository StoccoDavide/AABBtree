# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                        #
#                                                                                                 #
# The AABBtree project is distributed under the BSD 2-Clause License.                             #
#                                                                                                 #
# Davide Stocco                                                                 Enrico Bertolazzi #
# University of Trento                                                       University of Trento #
# davide.stocco@unitn.it                                               enrico.bertolazzi@unitn.it #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

include(AABBtreeMatplot)

file(GLOB_RECURSE BENCHMARK_FUNCTIONS "${CMAKE_CURRENT_SOURCE_DIR}/benchmark_functions.cc")
add_executable(benchmark_functions ${BENCHMARK_FUNCTIONS})
target_link_libraries(benchmark_functions PRIVATE AABBtree Matplot++::matplot)

file(GLOB_RECURSE BENCHMARK_ASTEROIDS "${CMAKE_CURRENT_SOURCE_DIR}/benchmark_asteroids.cc")
add_executable(benchmark_asteroids ${BENCHMARK_ASTEROIDS})
target_link_libraries(benchmark_asteroids PRIVATE AABBtree Matplot++::matplot)
