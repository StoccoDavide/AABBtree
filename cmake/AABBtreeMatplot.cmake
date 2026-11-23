# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                  #
#                                                                           #
# The AABBtree project is distributed under the BSD 2-Clause License.       #
#                                                                           #
# Davide Stocco Enrico Bertolazzi                                           #
# University of Trento University of Trento                                 #
# davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Matplot++ is used only by the tests and the benchmarks to plot the results.
# When plotting is disabled an empty target is defined, so that the targets can
# always link to Matplot++::matplot.

if(AABBTREE_ENABLE_PLOTTING)
  set(MATPLOT_REQUIRED_VERSION 1.2.0)
  list(APPEND CMAKE_PREFIX_PATH "${AABBTREE_THIRD_PARTY_DIR}")
  find_package(Matplot++ ${MATPLOT_REQUIRED_VERSION} NO_MODULE QUIET)
  if(NOT TARGET Matplot++::matplot)
    message(FATAL_ERROR "AABBtree: plotting is enabled but Matplot++ ${MATPLOT_REQUIRED_VERSION} "
      "was not found. Install it, or configure with -DAABBTREE_ENABLE_PLOTTING=OFF.")
  endif()
  get_target_property(MATPLOT_INCLUDE_DIRS Matplot++::matplot INTERFACE_INCLUDE_DIRECTORIES)
  message(STATUS "AABBtree: found Matplot++ in ${MATPLOT_INCLUDE_DIRS}")
elseif(NOT TARGET Matplot++::matplot)
  add_library(Matplot++::matplot INTERFACE IMPORTED)
  message(STATUS "AABBtree: Matplot++ plotting disabled")
endif()
