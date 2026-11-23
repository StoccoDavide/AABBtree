# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                  #
#                                                                           #
# The AABBtree project is distributed under the BSD 2-Clause License.       #
#                                                                           #
# Davide Stocco Enrico Bertolazzi                                           #
# University of Trento University of Trento                                 #
# davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

set(AABBTREE_EXAMPLES
  example_box
  example_ray
  example_tree
)

foreach(EXAMPLE_NAME IN LISTS AABBTREE_EXAMPLES)
  add_executable(${EXAMPLE_NAME} ${EXAMPLE_NAME}.cc)
  target_link_libraries(${EXAMPLE_NAME} PRIVATE AABBtree::AABBtree aabbtree_project_options)
endforeach()
