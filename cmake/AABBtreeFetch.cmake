# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                  #
#                                                                           #
# The AABBtree project is distributed under the BSD 2-Clause License.       #
#                                                                           #
# Davide Stocco Enrico Bertolazzi                                           #
# University of Trento University of Trento                                 #
# davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Make a dependency available with FetchContent.
#
# The sources are downloaded once into ${AABBTREE_THIRD_PARTY_DIR}/<name>-src
# and shared by all the build directories, while the dependency is built inside
# each build directory (so builds with different compilers or options do not
# overwrite each other).
#
#   aabbtree_fetch(<name> <url> [HEADERS_ONLY])
#
# With HEADERS_ONLY the sources are only downloaded: the CMake project of the
# dependency is not configured (useful for header-only libraries such as Eigen,
# whose project would also register its tests and build directory).

include(FetchContent)

# Extract the archives with the current time as timestamp (CMake >= 3.24)
if(POLICY CMP0135)
  cmake_policy(SET CMP0135 NEW)
endif()

macro(aabbtree_fetch NAME URL)
  cmake_parse_arguments(_aabbtree_fetch "HEADERS_ONLY" "" "" ${ARGN})
  string(TOLOWER "${NAME}" _aabbtree_name_lower)
  string(TOUPPER "${NAME}" _aabbtree_name_upper)
  set(_aabbtree_source_dir "${AABBTREE_THIRD_PARTY_DIR}/${_aabbtree_name_lower}-src")
  if(EXISTS "${_aabbtree_source_dir}/CMakeLists.txt")
    # Already downloaded: use the sources without downloading them again
    set(FETCHCONTENT_SOURCE_DIR_${_aabbtree_name_upper} "${_aabbtree_source_dir}")
    message(STATUS "AABBtree: using ${NAME} sources in ${_aabbtree_source_dir}")
  else()
    message(STATUS "AABBtree: downloading ${NAME} to ${_aabbtree_source_dir}")
  endif()
  if(_aabbtree_fetch_HEADERS_ONLY)
    # A SOURCE_SUBDIR without a CMakeLists.txt: the sources are downloaded but
    # not added to the build
    fetchcontent_declare(${NAME} URL "${URL}" SOURCE_DIR "${_aabbtree_source_dir}"
      SOURCE_SUBDIR "aabbtree-headers-only")
  else()
    fetchcontent_declare(${NAME} URL "${URL}" SOURCE_DIR "${_aabbtree_source_dir}")
  endif()
  fetchcontent_makeavailable(${NAME})
endmacro()
