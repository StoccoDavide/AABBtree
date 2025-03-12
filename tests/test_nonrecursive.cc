/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The AABBtree project is distributed under the BSD 2-Clause License.                           *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


// C++17 standard libraries
#include <vector>
#include <iostream>

// AABBtree library
#include "AABBtree.hh"
using namespace AABBtree;

// Catch2 library
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
using namespace Catch::Matchers;

// Matplot++ library
#ifdef AABBTREE_ENABLE_PLOTTING
#include <matplot/matplot.h>
using namespace matplot;
static auto fig{figure(true)};
static axes_handle ax{fig->current_axes()};
#endif

// Test utilities
#include "test_utilities.hh"

#ifdef AABBTREE_ENABLE_PLOTTING
#ifndef SET_PLOT
#define SET_PLOT \
xlim(ax, {-3.0, 3.0}); xlabel(ax, "x"); \
ylim(ax, {-3.0, 3.0}); ylabel(ax, "y"); \
grid(ax, true);
#endif
#endif

TEMPLATE_TEST_CASE("NonRecursive", "[template]", float, double) {

  std::vector<std::string> colors = {"r", "g", "b", "c", "m", "y", "k", "w"};

  using BoxUniquePtrList = AABBtree::BoxUniquePtrList<TestType, 2>;
  using Vector = AABBtree::Vector<TestType, 2>;
  using Box = AABBtree::Box<TestType, 2>;
  using Ray = AABBtree::Ray<TestType, 2>;

  #ifdef AABBTREE_ENABLE_PLOTTING
  SET_PLOT
  title(ax, "Build");
  #endif

  // Build segments
  Integer const n{100};
  TestType const scale{2.5};
  TestType const length{0.2};
  std::vector<Segment<TestType>> segments(n);
  #ifdef AABBTREE_ENABLE_PLOTTING
  ax->hold(true);
  #endif
  Eigen::Matrix<TestType, 2, 2> R;
  for (Integer i{0}; i < n; ++i) {
    segments[i].p_1() = scale*Vector::Random();
    R = Eigen::Matrix<TestType, 2, 2>::Random();
    R.col(0).normalize(); R.col(1).normalize();
    segments[i].p_2() = segments[i].p_1() + R*Vector(length, 0.0);
    #ifdef AABBTREE_ENABLE_PLOTTING
    plot_segment<TestType>(segments[i].p_1(), segments[i].p_2(), colors[6], 1.0);
    #endif
  }

  // Intersect segments
  std::vector<Vector> segment_points;
  for (Integer i{0}; i < n; ++i) {
    for (Integer j{i+1}; j < n; ++j) {
      Vector point;
      if (segments[i].intersect(segments[j], point)) {
        segment_points.push_back(point);
        #ifdef AABBTREE_ENABLE_PLOTTING
        plot_point<TestType, 2>(point, colors[0], 5.0);
        #endif
      }
  }}

  // Segments boxes
  std::unique_ptr<BoxUniquePtrList> boxes = std::make_unique<BoxUniquePtrList>();
  for (Integer i{0}; i < n; ++i) {
    Box const box{segments[i].box()};
    boxes->push_back(std::make_unique<Box>(box));
    #ifdef AABBTREE_ENABLE_PLOTTING
    plot_box<TestType, 2>(box, colors[3], 0.25);
    #endif
  }

  // Build tree
  AABBtree::NonRecursive<TestType, 2> tree;
  tree.build(std::move(boxes));
  tree.print(std::cout);
  #ifdef AABBTREE_ENABLE_PLOTTING
  plot_tree<TestType, 2>(tree, colors[2], 0.5);
  show(fig);
  #endif

  // Intersect tree
  IndexSet candidates;
  std::vector<Vector> tree_points;
  for (Integer i{0}; i < n; ++i) {
    Vector point;
    if (tree.intersect(segments[i].box(), candidates)) {
      std::vector<Integer> candidates_vec(candidates.begin(), candidates.end());
      for (Integer i{0}; i < static_cast<Integer>(candidates_vec.size()); ++i) {
        for (Integer j{i+1}; j < static_cast<Integer>(candidates_vec.size()); ++j) {
          std::cout << "Segment " << i << " intersects with segment " << j << std::endl;
          if (segments[candidates_vec[i]].intersect(segments[candidates_vec[j]], point)) {
            tree_points.push_back(point);
            #ifdef AABBTREE_ENABLE_PLOTTING
            plot_point<TestType, 2>(point, colors[1], 5.0);
            #endif
          }
        }
      }
    }
  }
  #ifdef AABBTREE_ENABLE_PLOTTING
  show(fig);
  #endif

  // Intersect ray
  Ray ray(0.0, 0.0, 1.0, 1.0);
  for (Integer i{0}; i < n; ++i) {
    if (tree.intersect(ray, candidates)) {
      for (const auto& i : candidates) {
        #ifdef AABBTREE_ENABLE_PLOTTING
        plot_box<TestType, 2>(*tree.box(i), colors[0], 2.0);
        #endif
      }
    }
  }
  #ifdef AABBTREE_ENABLE_PLOTTING
  plot_ray<TestType, 2>(ray, colors[0], 1.0);
  show(fig);
  #endif

  // Box distance
  Vector const point(0.0, 0.0);
  #ifdef AABBTREE_ENABLE_PLOTTING
  plot_point<TestType, 2>(point, colors[0], 5.0);
  #endif
  TestType const distance{tree.distance(point, candidates)};
  std::cout << "Distance: " << distance << std::endl;
  for (const auto& i : candidates) {
    #ifdef AABBTREE_ENABLE_PLOTTING
    plot_box<TestType, 2>(*tree.box(i), colors[0], 2.0);
    #endif
  }
  #ifdef AABBTREE_ENABLE_PLOTTING
  show(fig);
  ax->clear();
  #endif




  // Check results
  //REQUIRE(segment_points.size() == tree_points.size());
  //for (Integer i{0}; i < static_cast<Integer>(segment_points.size()); ++i) {
  //  REQUIRE_THAT(segment_points[i], IsApprox(tree_points[i]));
  //}
}
