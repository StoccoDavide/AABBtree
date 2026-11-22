# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                  #
#                                                                           #
# The AABBtree project is distributed under the BSD 2-Clause License.       #
#                                                                           #
# Davide Stocco Enrico Bertolazzi                                           #
# University of Trento University of Trento                                 #
# davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Look for OpenMP. Apple Clang does not ship the OpenMP runtime: in that case
# look for a libomp installed with Homebrew or MacPorts and pass the hints to
# FindOpenMP. If OpenMP is not found the library is used serially.

find_package(OpenMP QUIET COMPONENTS CXX)

if(NOT OpenMP_CXX_FOUND AND APPLE AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
  find_path(AABBTREE_LIBOMP_INCLUDE_DIR omp.h
    HINTS
      /opt/homebrew/opt/libomp/include
      /usr/local/opt/libomp/include
      /opt/local/include/libomp
  )
  find_library(AABBTREE_LIBOMP_LIBRARY NAMES omp
    HINTS
      /opt/homebrew/opt/libomp/lib
      /usr/local/opt/libomp/lib
      /opt/local/lib/libomp
  )
  if(AABBTREE_LIBOMP_INCLUDE_DIR AND AABBTREE_LIBOMP_LIBRARY)
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I${AABBTREE_LIBOMP_INCLUDE_DIR}" CACHE STRING "" FORCE)
    set(OpenMP_CXX_LIB_NAMES "omp" CACHE STRING "" FORCE)
    set(OpenMP_omp_LIBRARY "${AABBTREE_LIBOMP_LIBRARY}" CACHE FILEPATH "" FORCE)
    unset(OpenMP_CXX_FOUND)
    unset(OpenMP_CXX_FOUND CACHE)
    find_package(OpenMP QUIET COMPONENTS CXX)
  endif()
endif()

if(OpenMP_CXX_FOUND)
  message(STATUS "AABBtree: OpenMP ${OpenMP_CXX_VERSION} enabled")
else()
  message(STATUS "AABBtree: OpenMP not found, batch queries will run serially")
endif()
