# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                        #
#                                                                                                 #
# The AABBtree project is distributed under the BSD 2-Clause License.                             #
#                                                                                                 #
# Davide Stocco                                                                 Enrico Bertolazzi #
# University of Trento                                                       University of Trento #
# davide.stocco@unitn.it                                               enrico.bertolazzi@unitn.it #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

file(GLOB_RECURSE EXAMPLE_BOX "${CMAKE_CURRENT_SOURCE_DIR}/example_box.cc")
add_executable(example_box ${EXAMPLE_BOX})
target_link_libraries(example_box PRIVATE AABBtree)

file(GLOB_RECURSE EXAMPLE_RAY "${CMAKE_CURRENT_SOURCE_DIR}/example_ray.cc")
add_executable(example_ray ${EXAMPLE_RAY})
target_link_libraries(example_ray PRIVATE AABBtree)

file(GLOB_RECURSE EXAMPLE_TREE "${CMAKE_CURRENT_SOURCE_DIR}/example_tree.cc")
add_executable(example_tree ${EXAMPLE_TREE})
target_link_libraries(example_tree PRIVATE AABBtree)
