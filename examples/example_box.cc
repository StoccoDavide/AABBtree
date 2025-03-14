/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The AABBtree project is distributed under the BSD 2-Clause License.                           *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// A simple example of how to use the AABBtree library's Box class: for more details, see the documentation.

#include <iostream>
#include "AABBtree.hh"
using namespace AABBtree;

int main() {

  using Real = double;
  using Point = Eigen::Matrix<Real, 2, 1>;
  using Box = AABBtree::Box<Real, 2>;

  // Define min and max corners for a 2D box
  Point min_corner(0.0, 0.0);
  Point max_corner(5.0, 3.0);

  // Create a Box instance
  Box box(min_corner, max_corner);

  // Check if the box is empty
  std::cout << "Is box empty? " << (box.is_empty() ? "Yes" : "No") << std::endl;

  // Check if the box contains a point
  Point point(2.0, 1.0);
  std::cout << "Does the box contain the point? " << (box.contains(point) ? "Yes" : "No") << std::endl;

  // Check if the box intersects another box
  Box other_box(Point(1.0, 1.0), Point(3.0, 2.0));
  std::cout << "Does the box intersect the other box? " << (box.intersects(other_box) ? "Yes" : "No") << std::endl;

  // Compute the intersection of two boxes
  Box intersection_box;
  box.intersect(other_box, intersection_box);

  // Find the closest point and its distance in the box to a given point
  Point closest_point;
  Real min_distance{box.interior_distance(point, closest_point)};
  std::cout << "Minimum distance: " << min_distance << std::endl;

  // Find the farthest point and its distance in the box to a given point
  Point farthest_point;
  Real max_distance{box.exterior_distance(point, farthest_point)};
  std::cout << "Maximum distance: " << max_distance << std::endl;

  // Find the closest points between two boxes
  Point closest_point_box;
  Point closest_point_other_box;
  Real min_distance_boxes{box.interior_distance(other_box, closest_point_box, closest_point_other_box)};
  std::cout << "Minimum distance between boxes: " << min_distance_boxes << std::endl;

  // Find the farthest points between two boxes
  Point farthest_point_box;
  Point farthest_point_other_box;
  Real max_distance_boxes{box.exterior_distance(other_box, farthest_point_box, farthest_point_other_box)};
  std::cout << "Maximum distance between boxes: " << max_distance_boxes << std::endl;

  return 0;
}
