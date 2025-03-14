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
#include <cmath>

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

// Benchmark utilities
#include "BenchmarkUtilities.hh"
using namespace BenchmarkUtilities;

#ifdef AABBTREE_ENABLE_PLOTTING
#ifndef SET_PLOT
#define SET_PLOT \
xlim(ax, {-3.0, 3.0}); xlabel(ax, "x"); \
ylim(ax, {-3.0, 3.0}); ylabel(ax, "y"); \
grid(ax, true);
#endif
#endif

// High-frequency sine function
template <typename Real>
Real sine(Real const & x, Real const & f) {
  return std::sin(f*x*x);
}
template <typename Real>
Eigen::Vector<Real, Eigen::Dynamic> sine(Eigen::Vector<Real, Eigen::Dynamic> const & x, Real const & f) {
  Eigen::Vector<Real, Eigen::Dynamic> y(x.size());
  for (Integer i = 0; i < x.size(); ++i) {
    y(i) = sine(x(i), f);
  }
  return y;
}

// Main function
int main() {

  using Real = float;
  using VectorX = Eigen::Vector<Real, Eigen::Dynamic>;
  using Vector2 = Eigen::Vector<Real, 2>;
  using Box = AABBtree::Box<Real, 2>;
  using BoxUniquePtrList = AABBtree::BoxUniquePtrList<Real, 2>;
  using Segment = BenchmarkUtilities::Segment<Real>;

  //std::vector<std::string> colors = {"r", "g", "b", "c", "m", "y", "k", "w"};

  // Define the sampling points
  Integer const n_points{101};
  Real const bound{2.0};
  VectorX x_f1 = VectorX::LinSpaced(n_points, -bound, bound);
  VectorX y_f2 = VectorX::LinSpaced(n_points, -bound, bound);

  // Compute the function values
  Real const amplitude{2.0};
  VectorX y_f1(n_points); y_f1 = amplitude*sine<Real>(x_f1, 3.0*M_PI);
  VectorX x_f2(n_points); x_f2 = amplitude*sine<Real>(y_f2, 3.0*M_PI);

  // Compute the bounding boxes
  std::vector<Segment> segments_1(n_points-1);
  std::vector<Segment> segments_2(n_points-1);
  for (Integer i{0}; i < n_points; ++i) {
    segments_1[i].point(0) = Vector2(x_f1(i), y_f1(i));
    segments_1[i].point(1) = Vector2(x_f1(i+1), y_f1(i+1));
    segments_2[i].point(0) = Vector2(x_f2(i), y_f2(i));
    segments_2[i].point(1) = Vector2(x_f2(i+1), y_f2(i+1));
  }

  // Build trees
  std::unique_ptr<BoxUniquePtrList> boxes_1 = std::make_unique<BoxUniquePtrList>();
  std::unique_ptr<BoxUniquePtrList> boxes_2 = std::make_unique<BoxUniquePtrList>();
  for (Integer i{0}; i < n_points-1; ++i) {
    Box const box_1{segments_1[i].bounding_box()};
    Box const box_2{segments_2[i].bounding_box()};
    boxes_1->push_back(std::make_unique<Box>(box_1));
    boxes_2->push_back(std::make_unique<Box>(box_2));
  }
  AABBtree::Tree<Real, 2> tree_1;
  AABBtree::Tree<Real, 2> tree_2;
  tree_1.build(std::move(boxes_1));
  tree_2.build(std::move(boxes_2));
  tree_1.print(std::cout);
  tree_2.print(std::cout);

  // Plot the functions
  #ifdef AABBTREE_ENABLE_PLOTTING
  SET_PLOT
  std::vector<Real> x_f1_vec(x_f1.data(), x_f1.data() + x_f1.size());
  std::vector<Real> y_f1_vec(y_f1.data(), y_f1.data() + y_f1.size());
  std::vector<Real> x_f2_vec(x_f2.data(), x_f2.data() + n_points);
  std::vector<Real> y_f2_vec(y_f2.data(), y_f2.data() + n_points);
  ax->plot(x_f1_vec, y_f1_vec); ax->hold(true);
  ax->plot(x_f2_vec, y_f2_vec);
  show(fig);
  #endif

  // Intersect tree 1 with tree 2
  IndexMap candidates;
  if (tree_1.intersect(tree_2, candidates)) {
    std::cout << "Candidates: " << candidates.size() << std::endl;
    std::vector<Vector2> tree_points;
    Vector2 point;
    for (auto const & [key, value] : candidates) {
      for (auto const & val : value) {
        if (segments_1[key].intersect(segments_2[val], point)) {
          tree_points.push_back(point);
          #ifdef AABBTREE_ENABLE_PLOTTING
          plot_box<Real, 2>(*tree_1.box(key), "r", 2.0);
          plot_box<Real, 2>(*tree_2.box(val), "g", 2.0);
          //ax->plot({point.x()}, {point.y()}, "o");
          #endif
  }}}}
  #ifdef AABBTREE_ENABLE_PLOTTING
  show(fig);
  #endif

  return 0;
}