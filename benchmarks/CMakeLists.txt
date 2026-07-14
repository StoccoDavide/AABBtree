# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                  #
#                                                                           #
# The AABBtree project is distributed under the BSD 2-Clause License.       #
#                                                                           #
# Davide Stocco Enrico Bertolazzi                                           #
# University of Trento University of Trento                                 #
# davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

include(AABBtreeMatplot)

file(GLOB_RECURSE BENCHMARK_FUNCTIONS "${CMAKE_CURRENT_SOURCE_DIR}/benchmark_functions.cc")
add_executable(benchmark_functions ${BENCHMARK_FUNCTIONS})
target_link_libraries(benchmark_functions PRIVATE AABBtree Matplot++::matplot)

file(GLOB_RECURSE BENCHMARK_ASTEROIDS_CLUSTERING_CUSTOM "${CMAKE_CURRENT_SOURCE_DIR}/benchmark_asteroids_clustering_custom.cc")
add_executable(benchmark_asteroids_clustering_custom ${BENCHMARK_ASTEROIDS_CLUSTERING_CUSTOM})
target_link_libraries(benchmark_asteroids_clustering_custom PRIVATE AABBtree Matplot++::matplot)

file(GLOB_RECURSE BENCHMARK_ASTEROIDS_CLUSTERING_HAUSDORFF "${CMAKE_CURRENT_SOURCE_DIR}/benchmark_asteroids_clustering_hausdorff.cc")
add_executable(benchmark_asteroids_clustering_hausdorff ${BENCHMARK_ASTEROIDS_CLUSTERING_HAUSDORFF})
target_link_libraries(benchmark_asteroids_clustering_hausdorff PRIVATE AABBtree Matplot++::matplot)

file(GLOB_RECURSE BENCHMARK_ASTEROIDS_CLUSTERING_FRECHET "${CMAKE_CURRENT_SOURCE_DIR}/benchmark_asteroids_clustering_frechet.cc")
add_executable(benchmark_asteroids_clustering_frechet ${BENCHMARK_ASTEROIDS_CLUSTERING_FRECHET})
target_link_libraries(benchmark_asteroids_clustering_frechet PRIVATE AABBtree Matplot++::matplot)

file(GLOB_RECURSE BENCHMARK_ASTEROIDS_ORBIT "${CMAKE_CURRENT_SOURCE_DIR}/benchmark_asteroids_orbit.cc")
add_executable(benchmark_asteroids_orbit ${BENCHMARK_ASTEROIDS_ORBIT})
target_link_libraries(benchmark_asteroids_orbit PRIVATE AABBtree Matplot++::matplot)

file(GLOB_RECURSE BENCHMARK_ASTEROIDS_COMPLEXITY "${CMAKE_CURRENT_SOURCE_DIR}/benchmark_asteroids_complexity.cc")
add_executable(benchmark_asteroids_complexity ${BENCHMARK_ASTEROIDS_COMPLEXITY})
target_link_libraries(benchmark_asteroids_complexity PRIVATE AABBtree Matplot++::matplot)
