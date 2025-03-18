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
#include "TestUtilities.hh"
using namespace TestUtilities;

#ifdef AABBTREE_ENABLE_PLOTTING
#ifndef SET_PLOT
#define SET_PLOT \
xlim(ax, {-3.0, 3.0}); xlabel(ax, "x"); \
ylim(ax, {-3.0, 3.0}); ylabel(ax, "y"); \
grid(ax, true);
#endif
#endif

TEMPLATE_TEST_CASE("Tree-Tree", "[template]", float, double) {

  std::vector<std::string> colors = {"r", "g", "b", "c", "m", "y", "k", "w"};

  using BoxUniquePtrList = AABBtree::BoxUniquePtrList<TestType, 2>;
  using Vector = AABBtree::Vector<TestType, 2>;
  using Box = AABBtree::Box<TestType, 2>;

  #ifdef AABBTREE_ENABLE_PLOTTING
  SET_PLOT
  title(ax, "Build");
  #endif

  // Build segments batch 1
  Integer const n_1{50};
  TestType const scale_1{2.0};
  TestType const length_1{0.125};
  std::vector<Segment<TestType>> segments_1(n_1);
  #ifdef AABBTREE_ENABLE_PLOTTING
  ax->hold(true);
  #endif
  for (Integer i{0}; i < n_1; ++i) {
    Eigen::Matrix<TestType, 2, 2> R;
    segments_1[i].point(0) = scale_1*Vector::Random();
    R = Eigen::Matrix<TestType, 2, 2>::Random();
    R.col(0).normalize(); R.col(1).normalize();
    segments_1[i].point(1) = segments_1[i].point(0) + R*Vector(length_1, 0.0);
    #ifdef AABBTREE_ENABLE_PLOTTING
    plot_segment<TestType>(segments_1[i].point(0), segments_1[i].point(1), colors[6], 1.0);
    #endif
  }

  // Build segments batch 2
  Integer const n_2{75};
  TestType const scale_2{2.0};
  TestType const length_2{0.25};
  std::vector<Segment<TestType>> segments_2(n_2);
  #ifdef AABBTREE_ENABLE_PLOTTING
  ax->hold(true);
  #endif
  for (Integer i{0}; i < n_2; ++i) {
    Eigen::Matrix<TestType, 2, 2> R;
    segments_2[i].point(0) = scale_2*Vector::Random();
    R = Eigen::Matrix<TestType, 2, 2>::Random();
    R.col(0).normalize(); R.col(1).normalize();
    segments_2[i].point(1) = segments_2[i].point(0) + R*Vector(length_2, 0.0);
    #ifdef AABBTREE_ENABLE_PLOTTING
    plot_segment<TestType>(segments_2[i].point(0), segments_2[i].point(1), colors[2], 1.0);
    #endif
  }

  // Intersect segments of batch 1 with batch 2
  std::vector<Vector> segment_points;
  for (Integer i{0}; i < n_1; ++i) {
    for (Integer j{0}; j < n_2; ++j) {
      Vector point;
      if (segments_1[i].intersect(segments_2[j], point)) {
        segment_points.push_back(point);
        #ifdef AABBTREE_ENABLE_PLOTTING
        plot_point<TestType, 2>(point, colors[0], 5.0);
        #endif
  }}}

  // Segments 1 boxes
  std::unique_ptr<BoxUniquePtrList> boxes_1 = std::make_unique<BoxUniquePtrList>();
  for (Integer i{0}; i < n_1; ++i) {
    Box const box{segments_1[i].bounding_box()};
    boxes_1->push_back(std::make_unique<Box>(box));
    #ifdef AABBTREE_ENABLE_PLOTTING
    plot_box<TestType, 2>(box, colors[3], 0.25);
    #endif
  }

  // Segments 2 boxes
  std::unique_ptr<BoxUniquePtrList> boxes_2 = std::make_unique<BoxUniquePtrList>();
  for (Integer i{0}; i < n_2; ++i) {
    Box const box{segments_2[i].bounding_box()};
    boxes_2->push_back(std::make_unique<Box>(box));
    #ifdef AABBTREE_ENABLE_PLOTTING
    plot_box<TestType, 2>(box, colors[4], 0.25);
    #endif
  }

  // Build tree 1
  AABBtree::Tree<TestType, 2> tree_1;
  tree_1.build(std::move(boxes_1));
  tree_1.print(std::cout);
  #ifdef AABBTREE_ENABLE_PLOTTING
  plot_tree<TestType, 2>(tree_1, colors[2], 0.5);
  #endif

  // Build tree 2
  AABBtree::Tree<TestType, 2> tree_2;
  tree_2.build(std::move(boxes_2));
  tree_2.print(std::cout);
  #ifdef AABBTREE_ENABLE_PLOTTING
  plot_tree<TestType, 2>(tree_2, colors[6], 0.5);
  show(fig);
  #endif

  // Intersect trees
  IndexMap candidates;
  std::vector<Vector> tree_points;
  Vector point;
  if (tree_1.intersect(tree_2, candidates)) {
    for (auto const & [key, value] : candidates) {
      for (auto const & val : value) {
        if (segments_1[key].intersect(segments_2[val], point)) {
          tree_points.push_back(point);
          #ifdef AABBTREE_ENABLE_PLOTTING
          plot_point<TestType, 2>(point, colors[1], 5.0);
          #endif
  }}}}
  #ifdef AABBTREE_ENABLE_PLOTTING
  show(fig);
  #endif

  // Trees distance
  for (auto const & b_1 : tree_1.boxes()) {
    for (auto const & b_2 : tree_2.boxes()) {
      if (b_1->interior_distance(*b_2) == 0.0) {
        #ifdef AABBTREE_ENABLE_PLOTTING
        plot_box<TestType, 2>(*b_1, colors[0], 2.0);
        plot_box<TestType, 2>(*b_2, colors[0], 2.0);
        #endif
  }}}
  TestType const distance{tree_1.distance(tree_2, candidates)};
  if (distance >= 0.0) {
    #ifdef AABBTREE_ENABLE_PLOTTING
    for (auto const & [key, value] : candidates) {
      plot_box<TestType, 2>(*tree_1.box(key), colors[1], 1.0);
      for (auto const & val : value) {
        plot_box<TestType, 2>(*tree_2.box(val), colors[1], 1.0);
    }}
    #endif
  }
  #ifdef AABBTREE_ENABLE_PLOTTING
  show(fig);
  ax->clear();
  #endif
}
