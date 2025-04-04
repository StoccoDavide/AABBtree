/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The AABBtree project is distributed under the BSD 2-Clause License.                           *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// A simple example of how to use the AABBtree library's Ray class: for more details, see the documentation.

#include <iostream>
#include "AABBtree.hh"
using namespace AABBtree;

// Main function
int main() {

  using Real = double;
  using Point = Eigen::Matrix<Real, 2, 1>;
  using Vector = Eigen::Matrix<Real, 2, 1>;
  using Box = AABBtree::Box<Real, 2>;
  using Ray = AABBtree::Ray<Real, 2>;

  // Define origin and direction for a ray
  Point origin(0.0, 0.0);
  Vector direction(1.0, 1.0);

  // Create a Ray instance
  Ray ray(origin, direction);

  // Define min and max corners for a 2D box
  Point min_corner(0.0, 0.0);
  Point max_corner(5.0, 3.0);

  // Create a Box instance
  Box box(min_corner, max_corner);

  // Check if the ray intersects the box
  std::cout << "Does the ray intersect the box? " << (ray.intersects(box) ? "Yes\n" : "No\n");

  // Find the intersection point if the ray intersects the box
  Point closest_intersection_point, farthest_intersection_point;
  ray.intersect(box, closest_intersection_point, farthest_intersection_point);

  // Find the closest point and its distance in the box to the ray
  Point closest_point_on_ray, closest_point_on_box;
  Real min_distance{ray.interior_distance(box, closest_point_on_ray, closest_point_on_box)};
  std::cout << "Minimum distance: " << min_distance << '\n';

  // Find the farthest point and its distance in the box to the ray
  Point farthest_point_on_ray, farthest_point_on_box;
  Real max_distance{ray.exterior_distance(box, farthest_point_on_ray, farthest_point_on_box)};
  std::cout << "Maximum distance: " << max_distance << '\n';

  return 0;
}
