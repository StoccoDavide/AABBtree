/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The AABBtree project is distributed under the BSD 2-Clause License.                           *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// A simple example of how to use the AABBtree library's Tree class: for more details, see the documentation.

#include <vector>
#include <iostream>
#include "AABBtree.hh"
using namespace AABBtree;

int main() {

  using Real = double;
  using Point = Eigen::Matrix<Real, 2, 1>;
  using Box = AABBtree::Box<Real, 2>;
  using Ray = AABBtree::Ray<Real, 2>;

  // Define origin and direction for a ray
  Point origin(0.0, 0.0);
  Eigen::Matrix<Real, 2, 1> direction(1.0, 1.0);

  // Create a ray instance
  Ray ray(origin, direction);

  // Define min and max corners for 2D boxes
  Point min_corner1(0.0, 0.0), max_corner1(5.0, 3.0);
  Point min_corner2(3.0, 1.0), max_corner2(8.0, 5.0);
  Point min_corner3(7.0, 2.0), max_corner3(10.0, 6.0);

  // Create some box instances
  Box box1(min_corner1, max_corner1);
  Box box2(min_corner2, max_corner2);
  Box box3(min_corner3, max_corner3);

  // Create a list of boxes to build the tree
  std::unique_ptr<AABBtree::BoxUniquePtrList<Real, 2>> boxes =
    std::make_unique<AABBtree::BoxUniquePtrList<Real, 2>>();

  boxes->push_back(std::make_unique<Box>(box1));
  boxes->push_back(std::make_unique<Box>(box2));
  boxes->push_back(std::make_unique<Box>(box3));

  // Build the AABB tree
  AABBtree::Tree<Real, 2> tree;
  tree.build(std::move(boxes));
  tree.print(std::cout);

  // Intersect ray with the tree
  IndexSet candidates;
  if (tree.intersect(ray, candidates)) {
    std::cout << "Ray intersects with the following boxes:\n";
    for (const auto& idx : candidates) {
      std::cout << "Box " << idx << '\n';
    }
  } else {
    std::cout << "Ray does not intersect any boxes.\n";
  }

  // Check self-intersection of the tree
  candidates.clear();
  if (tree.self_intersect(candidates)) {
    std::cout << "Self-intersection detected between the following boxes:\n";
    for (const auto& idx : candidates) {
      std::cout << "Box " << idx << '\n';
    }
  } else {
    std::cout << "No self-intersections detected in the tree.\n";
  }

  // Compute the minimum distance between the tree and a point
  Point point(2.0, 1.0);
  IndexSet candidates_distance;
  Real min_distance{tree.distance(point, candidates_distance)};
  std::cout << "Minimum distance between the tree and the point: " << min_distance << '\n';

  // Let us build a tree with a single box
  Box box4(Point(0.0, 0.0), Point(1.0, 1.0));
  std::unique_ptr<AABBtree::BoxUniquePtrList<Real, 2>> boxes_single =
    std::make_unique<AABBtree::BoxUniquePtrList<Real, 2>>();
  boxes_single->push_back(std::make_unique<Box>(box4));
  AABBtree::Tree<Real, 2> tree_single;
  tree_single.build(std::move(boxes_single));
  tree_single.print(std::cout);

  // Intersect two trees
  IndexMap candidates_map;
  if (tree.intersect(tree_single, candidates_map)) {
    std::cout << "Intersection detected between the following boxes:\n";
    for (const auto& [key, values] : candidates_map) {
      for (const auto& val : values) {
        std::cout << "Box " << key << " and Box " << val << '\n';
      }
    }
  } else {
    std::cout << "No intersections detected between the trees.\n";
  }

  return 0;
}
