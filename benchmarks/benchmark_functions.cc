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

// Matplot++ library
#ifdef AABBTREE_ENABLE_PLOTTING
#include "Plot2D.hh"
#endif

using namespace AABBtree;

// Benchmark utilities
#include "BenchmarkUtilities.hh"
using namespace BenchmarkUtilities;

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
  using TicToc = BenchmarkUtilities::TicToc<Real>;

  // Start the timer
  TicToc timer;

  // Define the sampling points
  Integer const n_points{1001};
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
    segments_1[i].point(0) << x_f1(i), y_f1(i);
    segments_1[i].point(1) << x_f1(i+1), y_f1(i+1);
    segments_2[i].point(0) << x_f2(i), y_f2(i);
    segments_2[i].point(1) << x_f2(i+1), y_f2(i+1);
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
  timer.tic(); tree_1.build(std::move(boxes_1)); timer.toc();
  std::cout << std::fixed << "Tree 1 built in " << timer.elapsed_us() << " us\n";
  timer.tic(); tree_2.build(std::move(boxes_2)); timer.toc();
  std::cout << "Tree 2 built in " << timer.elapsed_us() << " us\n";
  tree_1.print(std::cout);
  tree_2.print(std::cout);

  // Plot the functions
  #ifdef AABBTREE_ENABLE_PLOTTING
  Plot2D P;
  P.xlim(-3.0, 3.0); P.xlabel("x");
  P.ylim(-3.0, 3.0); P.ylabel("y");
  P.grid(true);
  std::vector<Real> x_f1_vec(x_f1.data(), x_f1.data() + x_f1.size());
  std::vector<Real> y_f1_vec(y_f1.data(), y_f1.data() + y_f1.size());
  std::vector<Real> x_f2_vec(x_f2.data(), x_f2.data() + n_points);
  std::vector<Real> y_f2_vec(y_f2.data(), y_f2.data() + n_points);
  P.plot(x_f1_vec, y_f1_vec);
  P.hold(true);
  P.plot(x_f2_vec, y_f2_vec);
  P.show();
  #endif

  // Create a box
  Box box(-0.05, -0.05, 0.05, 0.05);

  // Intersect box with tree 1
  IndexSet candidates_set;
  bool do_intersect;
  timer.tic(); do_intersect = tree_1.intersect(box, candidates_set); timer.toc();
  std::cout << "Box intersects tree 1 in " << timer.elapsed_us() << " us\n";
  if (do_intersect) {
    std::cout << "Candidates: " << candidates_set.size() << '\n';
    #ifdef AABBTREE_ENABLE_PLOTTING
    for (auto const & i : candidates_set) {
      P.plot_box<Real, 2>( *tree_1.box(i), "r", 0.25);
    }
    #endif
  }

  // Intersect box with tree 2
  timer.tic(); do_intersect = tree_2.intersect(box, candidates_set); timer.toc();
  std::cout << "Box intersects tree 2 in " << timer.elapsed_us() << " us\n";
  if (do_intersect) {
    std::cout << "Candidates: " << candidates_set.size() << '\n';
    #ifdef AABBTREE_ENABLE_PLOTTING
    for (auto const & i : candidates_set) {
      P.plot_box<Real, 2>( *tree_2.box(i), "g", 0.25);
    }
    #endif
  }
  #ifdef AABBTREE_ENABLE_PLOTTING
  P.show();
  #endif

  // Create a ray
  Ray<Real, 2> ray(-2.0, -2.0, 1.0, 1.0);
  ray.direction().normalize();
  #ifdef AABBTREE_ENABLE_PLOTTING
  P.xlim(-3.0, 3.0); P.xlabel("x");
  P.ylim(-3.0, 3.0); P.ylabel("y");
  P.grid(true);
  P.plot_ray<Real, 2>( ray, "b", 1.0);
  #endif

  // Intersect ray with tree 1
  timer.tic(); do_intersect = tree_1.intersect(ray, candidates_set); timer.toc();
  std::cout << "Ray intersects tree 1 in " << timer.elapsed_us() << " us\n";
  if (do_intersect) {
    std::cout << "Candidates: " << candidates_set.size() << '\n';
    #ifdef AABBTREE_ENABLE_PLOTTING
    for (auto const & i : candidates_set) {
      P.plot_box<Real, 2>( *tree_1.box(i), "r", 0.25);
    }
    #endif
  }

  // Intersect ray with tree 2
  timer.tic(); do_intersect = tree_2.intersect(ray, candidates_set); timer.toc();
  std::cout << "Ray intersects tree 2 in " << timer.elapsed_us() << " us\n";
  if (do_intersect) {
    std::cout << "Candidates: " << candidates_set.size() << '\n';
    #ifdef AABBTREE_ENABLE_PLOTTING
    for (auto const & i : candidates_set) {
      P.plot_box<Real, 2>( *tree_1.box(i), "r", 0.25 );
    }
    #endif
  }
  #ifdef AABBTREE_ENABLE_PLOTTING
  P.show();
  #endif

  // Intersect tree 1 with tree 2
  IndexMap candidates_map;
  timer.tic(); do_intersect = tree_1.intersect(tree_2, candidates_map); timer.toc();
  std::cout << "Trees intersect in " << timer.elapsed_us() << " us\n";
  if (do_intersect) {
    std::cout << "Candidates: " << candidates_map.size() << '\n';
    std::vector<Vector2> tree_points;
    Vector2 point;
    for (auto const & [key, value] : candidates_map) {
      for (auto const & val : value) {
        if (segments_1[key].intersect(segments_2[val], point)) {
          tree_points.push_back(point);
          #ifdef AABBTREE_ENABLE_PLOTTING
          P.plot_box<Real, 2>( *tree_1.box(key), "r", 0.25 );
          P.plot_box<Real, 2>( *tree_2.box(val), "g", 0.25 );
          //ax->plot({point.x()}, {point.y()}, "o");
          #endif
  }}}}
  #ifdef AABBTREE_ENABLE_PLOTTING
  P.show();
  #endif

  return 0;
}