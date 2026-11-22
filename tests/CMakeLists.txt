# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                  #
#                                                                           #
# The AABBtree project is distributed under the BSD 2-Clause License.       #
#                                                                           #
# Davide Stocco Enrico Bertolazzi                                           #
# University of Trento University of Trento                                 #
# davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

include(AABBtreeCatch2)
include(AABBtreeMatplot)

file(GLOB_RECURSE TEST_BOX "${CMAKE_CURRENT_SOURCE_DIR}/test_box.cc")
add_executable(test_box ${TEST_BOX})
target_link_libraries(test_box PRIVATE AABBtree Catch2::Catch2WithMain Matplot++::matplot)

file(GLOB_RECURSE TEST_RAY "${CMAKE_CURRENT_SOURCE_DIR}/test_ray.cc")
add_executable(test_ray ${TEST_RAY})
target_link_libraries(test_ray PRIVATE AABBtree Catch2::Catch2WithMain Matplot++::matplot)

file(GLOB_RECURSE TEST_TREE_OBJS "${CMAKE_CURRENT_SOURCE_DIR}/test_tree_objs.cc")
add_executable(test_tree_objs ${TEST_TREE_OBJS})
target_link_libraries(test_tree_objs PRIVATE AABBtree Catch2::Catch2WithMain Matplot++::matplot)

file(GLOB_RECURSE TEST_TREE_TREE "${CMAKE_CURRENT_SOURCE_DIR}/test_tree_tree.cc")
add_executable(test_tree_tree ${TEST_TREE_TREE})
target_link_libraries(test_tree_tree PRIVATE AABBtree Catch2::Catch2WithMain Matplot++::matplot)

file(GLOB_RECURSE TEST_TREE_DISTANCE "${CMAKE_CURRENT_SOURCE_DIR}/test_tree_distance.cc")
add_executable(test_tree_distance ${TEST_TREE_DISTANCE})
target_link_libraries(test_tree_distance PRIVATE AABBtree Catch2::Catch2WithMain Matplot++::matplot)

file(GLOB_RECURSE TEST_TREE_INTERSECTION "${CMAKE_CURRENT_SOURCE_DIR}/test_tree_intersection.cc")
add_executable(test_tree_intersection ${TEST_TREE_INTERSECTION})
target_link_libraries(test_tree_intersection PRIVATE AABBtree Catch2::Catch2WithMain Matplot++::matplot)

file(GLOB_RECURSE TEST_TREE_BRUTEFORCE "${CMAKE_CURRENT_SOURCE_DIR}/test_tree_bruteforce.cc")
add_executable(test_tree_bruteforce ${TEST_TREE_BRUTEFORCE})
target_link_libraries(test_tree_bruteforce PRIVATE AABBtree Catch2::Catch2WithMain Matplot++::matplot)

# Register the tests with CTest
enable_testing()
foreach(TEST_NAME test_box test_ray test_tree_objs test_tree_tree test_tree_distance
                  test_tree_intersection test_tree_bruteforce)
  add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endforeach()
