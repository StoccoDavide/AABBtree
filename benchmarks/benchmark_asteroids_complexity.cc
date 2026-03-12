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
#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include <sstream>
#include <iomanip>
#include <chrono>

// Matplot++ library
#ifdef AABBTREE_ENABLE_PLOTTING
#include "Plot2D.hh"
#endif

// Benchmark utilities
#include "BenchmarkUtilities.hh"
using namespace BenchmarkUtilities;

using Real = double;
using AABBtree::Integer;

// Main function
std::tuple<Real, Real, Real, Real, Real, Real, Real, Real> ComputeChecks(int n_objects_t1, int n_objects_t2, int n_nodal_objs)
{

  constexpr Real    t_ini{64328.0}; // January 1, 2035
  constexpr Real    t_end{t_ini + 365.0};
  constexpr Integer t_steps{30};
  constexpr Real    tol_trace{1.0e-6};
  constexpr Real    dt{(t_end - t_ini)/t_steps};

  using Vector           = AABBtree::Vector<Real, 3*t_steps>;
  using Box              = AABBtree::Box<Real, 3*t_steps>;
  using Tree             = AABBtree::Tree<Real, 3*t_steps>;
  using Statistics       = Tree::Statistics;
  using BoxUniquePtrList = AABBtree::BoxUniquePtrList<Real, 3*t_steps>;

  // Parse asteroids data
  std::string fname{"./../benchmarks/asteroids.txt"}; // from build directory
  std::vector<Keplerian<Real, Integer>> data_1;
  if (!Parse<Real, Integer>(fname, data_1, n_objects_t1, false)) {
    std::cerr << "Error parsing data 1\n";
    return {};
  }
  std::vector<Keplerian<Real, Integer>> data_2;
  if (!Parse<Real, Integer>(fname, data_2, n_objects_t2, true)) {
    std::cerr << "Error parsing data 2\n";
    return {};
  }

  // Prepare the boxes 1
  std::unique_ptr<BoxUniquePtrList> boxes_1{ std::make_unique<BoxUniquePtrList>()};
  boxes_1->reserve(n_objects_t1);
  for (Integer i{0}; i < n_objects_t1; ++i) {
    Vector x_1, y_1, z_1;
    Keplerian<Real, Integer> data_i1 = data_1[i];
    for (Integer j{0}; j < t_steps; ++j) {
      PropagateOrbit(data_i1, t_ini + j*dt, (j == 0) ? t_ini : t_ini + (j-1)*dt);
      KeplerianToCartesian(data_i1, x_1(j), y_1(j), z_1(j));
    }
    Vector box_min_1, box_max_1;
    for (Integer j{0}; j < t_steps; ++j) {
      box_min_1(3*j+0) = x_1(j) - tol_trace;
      box_min_1(3*j+1) = y_1(j) - tol_trace;
      box_min_1(3*j+2) = z_1(j) - tol_trace;
      box_max_1(3*j+0) = x_1(j) + tol_trace;
      box_max_1(3*j+1) = y_1(j) + tol_trace;
      box_max_1(3*j+2) = z_1(j) + tol_trace;
    }
    boxes_1->push_back(std::make_unique<Box>(box_min_1, box_max_1));
    boxes_1->back()->reorder();
  }

  // Prepare the boxes 2
  std::unique_ptr<BoxUniquePtrList> boxes_2{ std::make_unique<BoxUniquePtrList>()};
  boxes_2->reserve(n_objects_t2);
  for (Integer i{0}; i < n_objects_t2; ++i) {
    Vector x_2, y_2, z_2;
    Keplerian<Real, Integer> data_i2{ data_2[i]};
    for (Integer j{0}; j < t_steps; ++j) {
      PropagateOrbit(data_i2, t_ini + j*dt, (j == 0) ? t_ini : t_ini + (j-1)*dt);
      KeplerianToCartesian(data_i2, x_2(j), y_2(j), z_2(j));
    }

    // swap y and z
    y_2.swap(z_2);

    Vector box_min_2, box_max_2;
    for (Integer j{0}; j < t_steps; ++j) {
      box_min_2(3*j+0) = x_2(j) - tol_trace;
      box_min_2(3*j+1) = y_2(j) - tol_trace;
      box_min_2(3*j+2) = z_2(j) - tol_trace;
      box_max_2(3*j+0) = x_2(j) + tol_trace;
      box_max_2(3*j+1) = y_2(j) + tol_trace;
      box_max_2(3*j+2) = z_2(j) + tol_trace;
    }
    boxes_2->push_back(std::make_unique<Box>(box_min_2, box_max_2));
    boxes_2->back()->reorder();
  }

  // Build the tree
  Tree tree_1, tree_2;
  tree_1.max_nodal_objects(n_nodal_objs);
  tree_2.max_nodal_objects(n_nodal_objs);
  tree_1.build(std::move(boxes_1));
  tree_2.build(std::move(boxes_2));

  // Query tree-tree
  TicToc<Real> timer;
  Statistics stats;
  IndexMap candidates_query_tt;

  Integer checks_i_tt{0};
  timer.tic();
  tree_1.intersect(tree_2, candidates_query_tt);
  tree_1.stats(stats);
  checks_i_tt += stats.check_counter;
  tree_2.intersect(tree_1, candidates_query_tt);
  tree_2.stats(stats);
  checks_i_tt += stats.check_counter;
  timer.toc();
  Real time_i_tt{timer.elapsed_s()};

  Integer checks_d_tt{0};
  timer.tic();
  tree_1.distance(tree_2, candidates_query_tt);
  tree_1.stats(stats);
  checks_d_tt += stats.check_counter;
  tree_2.distance(tree_1, candidates_query_tt);
  tree_2.stats(stats);
  checks_d_tt += stats.check_counter;
  timer.toc();
  Real time_d_tt{timer.elapsed_s()};

  // Query tree-boxes
  IndexSet candidates_query_bb;
  Integer  checks_i_tb{0};
  timer.tic();
  for (Integer i{0}; i < n_objects_t2; ++i) {
    tree_1.intersect(*tree_2.box(i), candidates_query_bb);
    tree_1.stats(stats);
    checks_i_tb += stats.check_counter;
  }
  for (Integer i{0}; i < n_objects_t1; ++i) {
    tree_2.intersect(*tree_1.box(i), candidates_query_bb);
    tree_2.stats(stats);
    checks_i_tb += stats.check_counter;
  }
  timer.toc();
  Real time_i_tb{timer.elapsed_s()};

  Integer checks_d_tb{0};
  timer.tic();
  for (Integer i{0}; i < n_objects_t2; ++i) {
    tree_1.distance(*tree_2.box(i), candidates_query_bb);
    tree_1.stats(stats);
    checks_d_tb += stats.check_counter;
  }
  for (Integer i{0}; i < n_objects_t1; ++i) {
    tree_2.distance(*tree_1.box(i), candidates_query_bb);
    tree_2.stats(stats);
    checks_d_tb += stats.check_counter;
  }
  timer.toc();
  Real time_d_tb{timer.elapsed_s()};

  return std::make_tuple(
    0.5*checks_i_tt, 0.5*checks_i_tb, 0.5*checks_d_tt, 0.5*checks_d_tb,
    0.5*time_i_tt, 0.5*time_i_tb, 0.5*time_d_tt, 0.5*time_d_tb
  );
}

// Main function
int
main() {

  #ifdef AABBTREE_ENABLE_PLOTTING
  Plot2D plot_checks_intersection;
  plot_checks_intersection.grid(true);
  plot_checks_intersection.hold(true);
  plot_checks_intersection.xlabel("Number of objects");
  plot_checks_intersection.ylabel("Number of checks");
  plot_checks_intersection.title("Intersection");
  Plot2D plot_checks_distance;
  plot_checks_distance.grid(true);
  plot_checks_distance.hold(true);
  plot_checks_distance.xlabel("Number of objects");
  plot_checks_distance.ylabel("Number of checks");
  plot_checks_distance.title("Distance");
  Plot2D plot_time_intersection;
  plot_time_intersection.grid(true);
  plot_time_intersection.hold(true);
  plot_time_intersection.xlabel("Number of objects");
  plot_time_intersection.ylabel("Time (s)");
  plot_time_intersection.title("Intersection");
  Plot2D plot_time_distance;
  plot_time_distance.grid(true);
  plot_time_distance.hold(true);
  plot_time_distance.xlabel("Number of objects");
  plot_time_distance.ylabel("Time (s)");
  plot_time_distance.title("Distance");
  #endif

  // Set the number of asteroids for the two trees
  std::vector<std::string> colors = {"r", "g", "b", "c", "m", "y", "k", "r", "g", "b", "c", "m", "y", "k"};
  Integer n_times{5};
  std::vector<Integer>     n_objects = {100, 177, 316, 562, 1000, 1778, 3162, 5623, 10000, 17783, 31623, 56234, 60000};
  std::vector<Integer>     n_nodal_objs = {1, 5, 10, 20, 40};
  std::vector<std::string> L{"Tree-Tree", "Tree-Boxes"};

  // Open files to save the results
  std::ofstream checks_intersection("checks_intersection_aabbtree.txt");
  std::ofstream checks_distance("checks_distance_aabbtree.txt");
  std::ofstream time_intersection("time_intersection_aabbtree.txt");
  std::ofstream time_distance("time_distance_aabbtree.txt");

  std::cout << "----------------------------------------\n";
  std::vector<Real> checks_i_tt(n_objects.size(), 0.0);
  std::vector<Real> checks_i_tb(n_objects.size(), 0.0);
  std::vector<Real> checks_d_tt(n_objects.size(), 0.0);
  std::vector<Real> checks_d_tb(n_objects.size(), 0.0);
  std::vector<Real> time_i_tt(n_objects.size(), 0.0);
  std::vector<Real> time_i_tb(n_objects.size(), 0.0);
  std::vector<Real> time_d_tt(n_objects.size(), 0.0);
  std::vector<Real> time_d_tb(n_objects.size(), 0.0);
  for (Integer j{0}; j < static_cast<Integer>(n_nodal_objs.size()); ++j) {
    for (Integer i{0}; i < static_cast<Integer>(n_objects.size()); ++i) {
      for (Integer k{0}; k < n_times; ++k) {
        auto [c_i_tt_tmp, c_i_tb_tmp, c_d_tt_tmp, c_d_tb_tmp, t_i_tt_tmp, t_i_tb_tmp, t_d_tt_tmp, t_d_tb_tmp] =
          ComputeChecks(n_objects[i], n_objects[i], n_nodal_objs[j]);
          checks_i_tt[i] += c_i_tt_tmp;
          checks_i_tb[i] += c_i_tb_tmp;
          checks_d_tt[i] += c_d_tt_tmp;
          checks_d_tb[i] += c_d_tb_tmp;
          time_i_tt[i]   += t_i_tt_tmp;
          time_i_tb[i]   += t_i_tb_tmp;
          time_d_tt[i]   += t_d_tt_tmp;
          time_d_tb[i]   += t_d_tb_tmp;
      }
      checks_i_tt[i] /= n_times;
      checks_i_tb[i] /= n_times;
      checks_d_tt[i] /= n_times;
      checks_d_tb[i] /= n_times;
      time_i_tt[i]   /= n_times;
      time_i_tb[i]   /= n_times;
      time_d_tt[i]   /= n_times;
      time_d_tb[i]   /= n_times;
      std::cout
        << "Number of asteroids: " << n_objects[i] << std::endl
        << "Number of objects:   " << n_nodal_objs[j] << std::endl
        << "Checks i_tt: " << checks_i_tt[i] << std::endl
        << "Checks i_tb: " << checks_i_tb[i] << std::endl
        << "Checks d_tt: " << checks_d_tt[i] << std::endl
        << "Checks d_tb: " << checks_d_tb[i] << std::endl
        << "Time   i_tt: " << time_i_tt[i] << std::endl
        << "Time   i_tb: " << time_i_tb[i] << std::endl
        << "Time   d_tt: " << time_d_tt[i] << std::endl
        << "Time   d_tb: " << time_d_tb[i] << std::endl
        << "----------------------------------------" << std::endl;
      checks_intersection
        << n_objects[i] << " " << n_nodal_objs[j] << " " << checks_i_tt[i] << " " << checks_i_tb[i] << std::endl;
      checks_distance
        << n_objects[i] << " " << n_nodal_objs[j] << " " << checks_d_tt[i] << " " << checks_d_tb[i] << std::endl;
      time_intersection
        << n_objects[i] << " " << n_nodal_objs[j] << " " << time_i_tt[i] << " " << time_i_tb[i] << std::endl;
      time_distance
        << n_objects[i] << " " << n_nodal_objs[j] << " " << time_d_tt[i] << " " << time_d_tb[i] << std::endl;
    }
    #ifdef AABBTREE_ENABLE_PLOTTING
    plot_checks_intersection.loglog(n_objects, checks_i_tt)->line_width(1.5).color(colors[j]);
    plot_checks_intersection.loglog(n_objects, checks_i_tb)->line_width(1.5).color(colors[j]).line_style("--");
    plot_checks_intersection.legend(L);
    plot_checks_distance.loglog(n_objects, checks_d_tt)->line_width(1.5).color(colors[j]);
    plot_checks_distance.loglog(n_objects, checks_d_tb)->line_width(1.5).color(colors[j]).line_style("--");
    plot_checks_distance.legend(L);
    plot_time_intersection.loglog(n_objects, time_i_tt)->line_width(1.5).color(colors[j]);
    plot_time_intersection.loglog(n_objects, time_i_tb)->line_width(1.5).color(colors[j]).line_style("--");
    plot_time_intersection.legend(L);
    plot_time_distance.loglog(n_objects, time_d_tt)->line_width(1.5).color(colors[j]);
    plot_time_distance.loglog(n_objects, time_d_tb)->line_width(1.5).color(colors[j]).line_style("--");
    plot_time_distance.legend(L);
    #endif
  }
  checks_intersection.close();
  checks_distance.close();
  time_intersection.close();
  time_distance.close();

  #ifdef AABBTREE_ENABLE_PLOTTING
  plot_checks_intersection.show();
  plot_checks_distance.show();
  plot_time_intersection.show();
  plot_time_distance.show();
  #endif

  return 0;
}
