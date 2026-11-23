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

set(AABBTREE_BENCHMARKS
  benchmark_functions
  benchmark_asteroids_clustering_custom
  benchmark_asteroids_clustering_hausdorff
  benchmark_asteroids_clustering_frechet
  benchmark_asteroids_orbit
  benchmark_asteroids_complexity
)

foreach(BENCHMARK_NAME IN LISTS AABBTREE_BENCHMARKS)
  add_executable(${BENCHMARK_NAME} ${BENCHMARK_NAME}.cc)
  target_link_libraries(${BENCHMARK_NAME} PRIVATE
    AABBtree::AABBtree aabbtree_project_options Matplot++::matplot)
endforeach()
