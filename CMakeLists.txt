# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                  #
#                                                                           #
# The AABBtree project is distributed under the BSD 2-Clause License.       #
#                                                                           #
# Davide Stocco Enrico Bertolazzi                                           #
# University of Trento University of Trento                                 #
# davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

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

# True when AABBtree is the main project, false when it is included by another
# project (add_subdirectory or FetchContent)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(AABBTREE_IS_TOP_LEVEL ON)
else()
  set(AABBTREE_IS_TOP_LEVEL OFF)
endif()

if(NOT AABBTREE_VERSION_PRERELEASE)
  set(AABBTREE_FULL_VERSION "${AABBTREE_VERSION}")
  if(AABBTREE_IS_TOP_LEVEL)
    write_version_file()
  endif()
else()
  set(AABBTREE_FULL_VERSION "${AABBTREE_VERSION}-${AABBTREE_VERSION_PRERELEASE}")
endif()
message(STATUS "AABBtree: version ${AABBTREE_FULL_VERSION}")

set(AABBTREE_THIRD_PARTY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/third_party")

# Settings that only make sense when AABBtree is the main project (they must
# not change the configuration of a parent project)
if(AABBTREE_IS_TOP_LEVEL)
  if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "AABBtree: CMAKE_BUILD_TYPE not specified, defaulting to Release")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
  endif()
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
  if(CMAKE_GENERATOR STREQUAL "Xcode")
    set(CMAKE_OSX_DEPLOYMENT_TARGET "15.0" CACHE STRING "Minimum macOS version")
    set(CMAKE_XCODE_SCHEME_WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
    set(CMAKE_XCODE_GENERATE_SCHEME TRUE)
  endif()
endif()

# ------------------------------------------------------------------------------
# Library (header-only)
# ------------------------------------------------------------------------------

include(GNUInstallDirs)
include(AABBtreeEigen3)

add_library(AABBtree INTERFACE)
add_library(AABBtree::AABBtree ALIAS AABBtree)

# Only the usage requirements needed to compile the headers are propagated to
# the users of the library: C++20, the include directory and the dependencies
target_compile_features(AABBtree INTERFACE cxx_std_20)
target_link_libraries(AABBtree INTERFACE Eigen3::Eigen)
target_include_directories(AABBtree INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# Optional OpenMP support (parallel build, batch and tree-tree queries)
option(AABBTREE_ENABLE_OPENMP "Enable OpenMP parallelism in AABBtree" ON)
set(AABBTREE_USE_OPENMP OFF)
if(AABBTREE_ENABLE_OPENMP)
  include(AABBtreeOpenMP)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(AABBtree INTERFACE OpenMP::OpenMP_CXX)
    set(AABBTREE_USE_OPENMP ON)
  endif()
endif()

# ------------------------------------------------------------------------------
# Tests, benchmarks and examples (only when AABBtree is the main project)
# ------------------------------------------------------------------------------

if(AABBTREE_IS_TOP_LEVEL)
  option(BUILD_TESTS "Build tests" ON)
  option(BUILD_BENCHMARKS "Build benchmarks" ON)
  option(BUILD_EXAMPLES "Build examples" ON)
  option(AABBTREE_ENABLE_NATIVE_OPTIMIZATIONS
    "Compile tests, benchmarks and examples with -march=native -ffast-math" ON)
  if(CMAKE_GENERATOR STREQUAL "Xcode")
    option(AABBTREE_ENABLE_PLOTTING "Enable plotting with Matplot++" ON)
  else()
    option(AABBTREE_ENABLE_PLOTTING "Enable plotting with Matplot++" OFF)
  endif()

  # Show the headers in the IDEs. A separate target is used, so that the
  # exported AABBtree target does not carry any source path.
  file(GLOB_RECURSE AABBTREE_HEADERS CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hh")
  add_custom_target(AABBtree_headers SOURCES ${AABBTREE_HEADERS})

  # Compiler options for the code of this project: strict warnings, and
  # optionally aggressive optimizations. They are not propagated to the users
  # of the library.
  add_library(aabbtree_project_options INTERFACE)
  if(MSVC)
    target_compile_options(aabbtree_project_options INTERFACE /W4 /WX)
  else()
    target_compile_options(aabbtree_project_options INTERFACE
      -Wall -Wextra -pedantic -Werror)
    if(AABBTREE_ENABLE_NATIVE_OPTIMIZATIONS)
      set(AABBTREE_NATIVE_FLAGS -march=native -ffast-math -fno-finite-math-only -funroll-loops)
      target_compile_options(aabbtree_project_options INTERFACE
        "$<$<OR:$<CONFIG:Release>,$<CONFIG:RelWithDebInfo>>:${AABBTREE_NATIVE_FLAGS}>")
    endif()
  endif()
  if(AABBTREE_ENABLE_PLOTTING)
    target_compile_definitions(aabbtree_project_options INTERFACE AABBTREE_ENABLE_PLOTTING)
  endif()

  if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
  endif()
  if(BUILD_BENCHMARKS)
    add_subdirectory(benchmarks)
  endif()
  if(BUILD_EXAMPLES)
    add_subdirectory(examples)
  endif()
endif()

# ------------------------------------------------------------------------------
# Installation and package configuration
# ------------------------------------------------------------------------------

include(CMakePackageConfigHelpers)

set(AABBTREE_INSTALL_CMAKEDIR "${CMAKE_INSTALL_DATADIR}/cmake/AABBtree")

install(TARGETS AABBtree EXPORT AABBtreeTargets)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

export(EXPORT AABBtreeTargets
  FILE ${CMAKE_CURRENT_BINARY_DIR}/AABBtreeTargets.cmake
  NAMESPACE AABBtree::
)

install(EXPORT AABBtreeTargets
  FILE AABBtreeTargets.cmake
  NAMESPACE AABBtree::
  DESTINATION ${AABBTREE_INSTALL_CMAKEDIR}
)

configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/AABBtreeConfig.cmake
  INSTALL_DESTINATION ${AABBTREE_INSTALL_CMAKEDIR}
)

# The library is header-only: the package fits any architecture
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/AABBtreeConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
  ARCH_INDEPENDENT
)

# The OpenMP module is shipped with the package, so that the users find the
# OpenMP runtime the same way (e.g. libomp with Apple Clang)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/AABBtreeOpenMP.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/AABBtreeOpenMP.cmake COPYONLY)

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/AABBtreeConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/AABBtreeConfigVersion.cmake
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/AABBtreeOpenMP.cmake
  DESTINATION ${AABBTREE_INSTALL_CMAKEDIR}
)
