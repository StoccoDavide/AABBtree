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

TEMPLATE_TEST_CASE("Tree-Point-Ray)", "[template]", float, double) {

  std::vector<std::string> colors = {"r", "g", "b", "c", "m", "y", "k", "w"};

  using BoxUniquePtrList = AABBtree::BoxUniquePtrList<TestType, 2>;
  using Vector = AABBtree::Vector<TestType, 2>;
  using Point = AABBtree::Point<TestType, 2>;
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
    segments[i].point(0) = scale*Vector::Random();
    R = Eigen::Matrix<TestType, 2, 2>::Random();
    R.col(0).normalize(); R.col(1).normalize();
    segments[i].point(1) = segments[i].point(0) + R*Vector(length, 0.0);
    #ifdef AABBTREE_ENABLE_PLOTTING
    plot_segment<TestType>(segments[i].point(0), segments[i].point(1), colors[6], 1.0);
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
  }}}

  // Segments boxes
  std::unique_ptr<BoxUniquePtrList> boxes = std::make_unique<BoxUniquePtrList>();
  for (Integer i{0}; i < n; ++i) {
    Box const box{segments[i].bounding_box()};
    boxes->push_back(std::make_unique<Box>(box));
    #ifdef AABBTREE_ENABLE_PLOTTING
    plot_box<TestType, 2>(box, colors[3], 0.25);
    #endif
  }

  // Build tree
  AABBtree::Tree<TestType, 2> tree;
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
    if (tree.intersect(segments[i].bounding_box(), candidates)) {
      std::vector<Integer> candidates_vec(candidates.begin(), candidates.end());
      for (Integer i{0}; i < static_cast<Integer>(candidates_vec.size()); ++i) {
        for (Integer j{i+1}; j < static_cast<Integer>(candidates_vec.size()); ++j) {
          if (segments[candidates_vec[i]].intersect(segments[candidates_vec[j]], point)) {
            tree_points.push_back(point);
            #ifdef AABBTREE_ENABLE_PLOTTING
            plot_point<TestType, 2>(point, colors[1], 5.0);
            #endif
  }}}}}
  #ifdef AABBTREE_ENABLE_PLOTTING
  show(fig);
  #endif

  // Self intersect tree
  //candidates.clear();
  //tree_points.clear();
  //for (Integer i{0}; i < n; ++i) {
  //  Vector point;
  //  if (tree.self_intersect(candidates)) {
  //    std::vector<Integer> candidates_vec(candidates.begin(), candidates.end());
  //    for (Integer i{0}; i < static_cast<Integer>(candidates_vec.size()); ++i) {
  //      for (Integer j{i+1}; j < static_cast<Integer>(candidates_vec.size()); ++j) {
  //        if (segments[candidates_vec[i]].intersect(segments[candidates_vec[j]], point)) {
  //          tree_points.push_back(point);
  //          #ifdef AABBTREE_ENABLE_PLOTTING
  //          plot_point<TestType, 2>(point, colors[2], 5.0);
  //          #endif
  //}}}}}
  //#ifdef AABBTREE_ENABLE_PLOTTING
  //show(fig);
  //#endif

  // Intersect ray
  Ray ray(0.0, 0.0, 1.0, 1.0);
  for (Integer i{0}; i < n; ++i) {
    if (tree.intersect(ray, candidates)) {
      #ifdef AABBTREE_ENABLE_PLOTTING
      for (const auto & i : candidates) {
        plot_box<TestType, 2>(*tree.box(i), colors[0], 2.0);
      }
      #endif
  }}
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
  if (distance > 0.0) {
    #ifdef AABBTREE_ENABLE_PLOTTING
    for (const auto & i : candidates) {
      plot_box<TestType, 2>(*tree.box(i), colors[0], 2.0);
    }
    #endif
  }
  #ifdef AABBTREE_ENABLE_PLOTTING
  show(fig);
  #endif

  // Within distance
  TestType const max_distance{1.0};
  auto distance_func = [](Point const & p, Box const & b) -> TestType {return b.interior_distance(p);};
  #ifdef AABBTREE_ENABLE_PLOTTING
  for (auto const & b : tree.boxes()) {
    if (distance_func(point, *b) <= max_distance) {
      plot_box<TestType, 2>(*b, colors[0], 5.0);
    }
  }
  #endif
  TestType const distance_within{tree.within_distance(point, max_distance, candidates, distance_func)};
  if (distance_within > 0.0) {
    #ifdef AABBTREE_ENABLE_PLOTTING
    plot_circle<TestType, 2>(point, max_distance, colors[0], 2.0);
    for (const auto & i : candidates) {
      plot_box<TestType, 2>(*tree.box(i), colors[1], 2.5);
    }
    #endif
  }
  #ifdef AABBTREE_ENABLE_PLOTTING
  show(fig);
  ax->clear();
  #endif
}
