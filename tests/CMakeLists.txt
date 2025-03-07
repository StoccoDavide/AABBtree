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

#file(GLOB_RECURSE TEST_RECURSIVE "${CMAKE_CURRENT_SOURCE_DIR}/test_recursive.cc")
#add_executable(test_recursive ${TEST_RECURSIVE})
#target_link_libraries(test_recursive PRIVATE AABBtree Catch2::Catch2WithMain Matplot++::matplot)

file(GLOB_RECURSE TEST_NONRECURSIVE "${CMAKE_CURRENT_SOURCE_DIR}/test_nonrecursive.cc")
add_executable(test_nonrecursive ${TEST_NONRECURSIVE})
target_link_libraries(test_nonrecursive PRIVATE AABBtree Catch2::Catch2WithMain Matplot++::matplot)
