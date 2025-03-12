# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                        #
#                                                                                                 #
# The AABBtree project is distributed under the BSD 2-Clause License.                             #
#                                                                                                 #
# Davide Stocco                                                                 Enrico Bertolazzi #
# University of Trento                                                       University of Trento #
# davide.stocco@unitn.it                                               enrico.bertolazzi@unitn.it #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

include(AABBtreeCatch2)
include(AABBtreeMatplot)

file(GLOB_RECURSE TEST_BOX "${CMAKE_CURRENT_SOURCE_DIR}/test_box.cc")
add_executable(test_box ${TEST_BOX})
target_link_libraries(test_box PRIVATE AABBtree Catch2::Catch2WithMain Matplot++::matplot)

file(GLOB_RECURSE TEST_RAY "${CMAKE_CURRENT_SOURCE_DIR}/test_ray.cc")
add_executable(test_ray ${TEST_RAY})
target_link_libraries(test_ray PRIVATE AABBtree Catch2::Catch2WithMain Matplot++::matplot)

file(GLOB_RECURSE TEST_TREE "${CMAKE_CURRENT_SOURCE_DIR}/test_tree.cc")
add_executable(test_tree ${TEST_TREE})
target_link_libraries(test_tree PRIVATE AABBtree Catch2::Catch2WithMain Matplot++::matplot)
