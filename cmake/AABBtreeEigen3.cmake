# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                  #
#                                                                           #
# The AABBtree project is distributed under the BSD 2-Clause License.       #
#                                                                           #
# Davide Stocco Enrico Bertolazzi                                           #
# University of Trento University of Trento                                 #
# davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

set(EIGEN_REQUIRED_VERSION 3.4.0)

list(APPEND CMAKE_PREFIX_PATH "${AABBTREE_THIRD_PARTY_DIR}")
find_package(Eigen3 ${EIGEN_REQUIRED_VERSION} NO_MODULE QUIET)

if(TARGET Eigen3::Eigen)
  get_target_property(EIGEN_INCLUDE_DIRS Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
  message(STATUS "AABBtree: found Eigen3 in ${EIGEN_INCLUDE_DIRS}")
else()
  # Eigen is header-only: download it and define the target directly. As an
  # imported target its headers are treated as system headers, as those of an
  # installed Eigen.
  include(AABBtreeFetch)
  aabbtree_fetch(Eigen3
    "https://gitlab.com/libeigen/eigen/-/archive/${EIGEN_REQUIRED_VERSION}/eigen-${EIGEN_REQUIRED_VERSION}.tar.gz"
    HEADERS_ONLY)
  add_library(Eigen3::Eigen INTERFACE IMPORTED GLOBAL)
  set_target_properties(Eigen3::Eigen PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${eigen3_SOURCE_DIR}")
endif()
