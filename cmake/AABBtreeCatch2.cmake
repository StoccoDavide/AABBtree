# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                  #
#                                                                           #
# The AABBtree project is distributed under the BSD 2-Clause License.       #
#                                                                           #
# Davide Stocco Enrico Bertolazzi                                           #
# University of Trento University of Trento                                 #
# davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

set(CATCH2_REQUIRED_VERSION 3.7.1)

list(APPEND CMAKE_PREFIX_PATH "${AABBTREE_THIRD_PARTY_DIR}")
find_package(Catch2 ${CATCH2_REQUIRED_VERSION} NO_MODULE QUIET)

if(TARGET Catch2::Catch2WithMain)
  message(STATUS "AABBtree: found Catch2 ${Catch2_VERSION}")
else()
  include(AABBtreeFetch)
  aabbtree_fetch(Catch2
    "https://github.com/catchorg/Catch2/archive/refs/tags/v${CATCH2_REQUIRED_VERSION}.tar.gz")
  if(NOT TARGET Catch2::Catch2WithMain)
    message(FATAL_ERROR "AABBtree: Catch2::Catch2WithMain target not found after FetchContent")
  endif()
endif()
