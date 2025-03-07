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

  //using IndexSet = AABBtree::IndexSet;

  using Vector = AABBtree::Vector<TestType, 2>;
  using Box = AABBtree::Box<TestType, 2>;
  using BoxUniquePtrList = AABBtree::BoxUniquePtrList<TestType, 2>;

  SET_PLOT
  title(ax, "Build");

  // Build segments
  Integer const n{10};
  TestType const scale{2.0};
  TestType const length{0.5};
  std::vector<Segment<TestType>> segments(n);
  ax->hold(true);
  Eigen::Matrix<TestType, 2, 2> R;
  for (Integer i{0}; i < n; ++i) {
    segments[i].p_1() = scale*Vector::Random();
    R = Eigen::Matrix<TestType, 2, 2>::Random();
    R.col(0).normalize(); R.col(1).normalize();
    segments[i].p_2() = segments[i].p_1() + R*Vector(length, 0.0);
    plot_segment<TestType>(segments[i].p_1(), segments[i].p_2(), colors[6], 1.0);
  }
  show(fig);

  // Intersect segments
  std::vector<Vector> points;
  for (Integer i{0}; i < n; ++i) {
    for (Integer j{i+1}; j < n; ++j) {
      Vector point;
      if (segments[i].intersect(segments[j], point)) {
        points.push_back(point);
        plot_point<TestType, 2>(point, colors[0], 5.0);
      }
  }}
  show(fig);

  // Segments boxes
  std::unique_ptr<BoxUniquePtrList> boxes = std::make_unique<BoxUniquePtrList>();
  for (Integer i{0}; i < n; ++i) {
    Box const box{segments[i].box()};
    boxes->push_back(std::make_unique<Box>(box));
    plot_box<TestType, 2>(box, colors[3], 0.25);
  }
  show(fig);

  // Build tree
  AABBtree::NonRecursive<TestType, 2> tree;
  tree.build(std::move(boxes));
  plot_tree<TestType, 2>(tree, colors[2], 0.5);
  show(fig);

  // Intersect tree
  //IndexSet candidates;
  //for (Integer i{0}; i < n; ++i) {
  //  for (Integer j{i+1}; j < n; ++j) {
  //    Vector point;
  //    if (tree.intersect(segments[i])) {
  //      points.push_back(point);
  //      plot_point<TestType, 2>(point, colors[3], 0.5);
  //    }
  //}}
  //plot_tree<TestType, 2>(tree, colors[1], 1.0); show(fig);

  ax->hold(false);
  ax->clear();
}
