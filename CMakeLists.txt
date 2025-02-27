# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Copyright (c) 2025, Davide Stocco, Mattia Piazza and Enrico Bertolazzi.                         #
#                                                                                                 #
# The AABBtree project is distributed under the BSD 2-Clause License.                             #
#                                                                                                 #
# Davide Stocco                                                                 Enrico Bertolazzi #
# University of Trento                                                       University of Trento #
# davide.stocco@unitn.it                                               enrico.bertolazzi@unitn.it #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

cmake_minimum_required(VERSION 3.14)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(AABBtreeVersion)
get_version_from_git()

project(AABBtree
  VERSION ${AABBTREE_VERSION}
  DESCRIPTION "A C++ non-recursive ND AABB tree"
  HOMEPAGE_URL "https://github.com/StoccoDavide/AABBtree"
  LANGUAGES CXX
)

if(NOT AABBTREE_VERSION_PRERELEASE)
  set(AABBTREE_FULL_VERSION "${AABBTREE_VERSION}")
  write_version_file()
else()
  set(AABBTREE_FULL_VERSION "${AABBTREE_VERSION}-${AABBTREE_VERSION_PRERELEASE}")
endif()
message(STATUS "AABBtree: version ${AABBTREE_FULL_VERSION}")

set(AABBTREE_THIRD_PARTY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/third_party")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

if(NOT CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/lib")
  message(STATUS "AABBtree: CMAKE_INSTALL_PREFIX not specified, defaulting to ${CMAKE_CURRENT_SOURCE_DIR}")
endif()

# Dependencies
include(AABBtreeEigen3)

# Library definition
add_library(AABBtree INTERFACE)
add_library(AABBtree::AABBtree ALIAS AABBtree)

target_link_libraries(AABBtree INTERFACE Eigen3::Eigen)

target_include_directories(AABBtree INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

# INFO: Ensure this options are only available to developers. In theory we could leave them available,
# but the call to `target_sources` breaks `fetchcontent` compatibility due to absolute paths being
# added to `INTERFACE_SOURCES`. I tried solving it, but it seems to be poorly documented, supported, etc.
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug")
    message(STATUS "AABBtree: CMAKE_BUILD_TYPE not specified, defaulting to ${CMAKE_BUILD_TYPE}")
  endif()

  option(BUILD_TESTS "Build tests" OFF)
  option(BUILD_BENCHMARKS "Build benchmarks" OFF)
  option(BUILD_EXAMPLES "Build examples" OFF)

  file(GLOB_RECURSE HEADER_FILES_HH "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hh")
  file(GLOB_RECURSE HEADER_FILES_HXX "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hxx")
  message(STATUS "AABBtree: found ${HEADER_FILES_HXX} header files")
  foreach(HEADER_FILE IN LISTS HEADER_FILES_HH HEADER_FILES_HXX)
    file(
      RELATIVE_PATH
      RELATIVE_PATH_HEADER_FILE
      ${CMAKE_CURRENT_SOURCE_DIR} ${HEADER_FILE}
    )
    list(APPEND RELATIVE_HEADER_FILES ${RELATIVE_PATH_HEADER_FILE})
  endforeach()
  target_sources(AABBtree INTERFACE $<BUILD_INTERFACE:${RELATIVE_HEADER_FILES}>)

  if(BUILD_TESTS)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/tests)
  endif()

  if(BUILD_BENCHMARKS)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/benchmarks)
  endif()

  if(BUILD_EXAMPLES)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/examples)
  endif()
endif()

if(MSVC)
  target_compile_options(AABBtree INTERFACE /W4 /WX)
else()
  target_compile_options(AABBtree INTERFACE -Wall -Wextra -pedantic -Werror)
endif()

# Installation
install(TARGETS AABBtree EXPORT AABBtreeTargets)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION include)

include(CMakePackageConfigHelpers)

export(EXPORT AABBtreeTargets
  FILE ${CMAKE_CURRENT_BINARY_DIR}/AABBtreeTargets.cmake
  NAMESPACE AABBtree::
)

install(EXPORT AABBtreeTargets
  FILE AABBtreeTargets.cmake
  NAMESPACE AABBtree::
  DESTINATION lib/cmake/AABBtree
)

configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/AABBtreeConfig.cmake
  INSTALL_DESTINATION lib/cmake/AABBtree
)

write_basic_package_version_file(
  AABBtreeConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/AABBtreeConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/AABBtreeConfigVersion.cmake
  DESTINATION lib/cmake/AABBtree
)
