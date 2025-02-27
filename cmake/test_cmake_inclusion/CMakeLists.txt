cmake_minimum_required(VERSION 3.14)

project(test)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

option(ADD_SUBDIRECTORY_METHOD OFF)
option(FETCH_CONTENT_METHOD ON)
option(FIND_PACKAGE_METHOD OFF)

set(CMAKE_BUILD_TYPE Release)

add_executable(test test.cc)

if(ADD_SUBDIRECTORY_METHOD)
  add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/third_party/AABBtree")
endif()

if(FIND_PACKAGE_METHOD)
  list(APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/third_party")
  find_package(
    AABBtree
    0.0.0
    NO_MODULE
  )
endif()

if(FETCH_CONTENT_METHOD)
  include(FetchContent)

  set(FETCHCONTENT_BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/third_party")
  fetchcontent_declare(
    AABBtree
    GIT_REPOSITORY "https://github.com/StoccoDavide/AABBtree.git"
    GIT_TAG        main
  )

  fetchcontent_makeavailable(AABBtree)
endif()

target_link_libraries(test PRIVATE AABBtree::AABBtree)
