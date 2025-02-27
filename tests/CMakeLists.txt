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

file(GLOB_RECURSE TEST_ALIGNEDBOX "${CMAKE_CURRENT_SOURCE_DIR}/test_alignedbox.cc")
add_executable(test_alignedbox ${TEST_ALIGNEDBOX})
target_link_libraries(test_alignedbox PRIVATE AABBtree Catch2::Catch2WithMain)

file(GLOB_RECURSE TEST_AABBOX "${CMAKE_CURRENT_SOURCE_DIR}/test_aabbox.cc")
add_executable(test_aabbox ${TEST_AABBOX})
target_link_libraries(test_aabbox PRIVATE AABBtree Catch2::Catch2WithMain)
