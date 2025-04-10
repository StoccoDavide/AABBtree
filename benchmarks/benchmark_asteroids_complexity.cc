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
std::tuple<Real, Real, Real, Real> ComputeChecks(int n_objects_t1, int n_objects_t2, int n_nodal_objs)
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
  Statistics stats;
  IndexMap candidates_query_tt;

  tree_1.intersect(tree_2, candidates_query_tt);
  tree_1.stats(stats);
  Integer checks_query_it1t2{stats.check_counter};
  tree_2.intersect(tree_1, candidates_query_tt);
  tree_2.stats(stats);
  Integer checks_query_it2t1{stats.check_counter};

  tree_1.distance(tree_2, candidates_query_tt);
  tree_1.stats(stats);
  Integer checks_query_dt1t2{stats.check_counter};
  tree_2.distance(tree_1, candidates_query_tt);
  tree_2.stats(stats);
  Integer checks_query_dt2t1{stats.check_counter};

  // Query tree-boxes
  IndexSet candidates_query_bb;
  Integer  checks_query_it1b2{0};
  for (Integer i{0}; i < n_objects_t2; ++i) {
    tree_1.intersect(*tree_2.box(i), candidates_query_bb);
    tree_1.stats(stats);
    checks_query_it1b2 += stats.check_counter;
  }
  Integer checks_query_it2b1{0};
  for (Integer i{0}; i < n_objects_t1; ++i) {
    tree_2.intersect(*tree_1.box(i), candidates_query_bb);
    tree_2.stats(stats);
    checks_query_it2b1 += stats.check_counter;
  }

  Integer checks_query_dt1b2{0};
  for (Integer i{0}; i < n_objects_t2; ++i) {
    tree_1.distance(*tree_2.box(i), candidates_query_bb);
    tree_1.stats(stats);
    checks_query_dt1b2 += stats.check_counter;
  }
  Integer checks_query_dt2b1{0};
  for (Integer i{0}; i < n_objects_t1; ++i) {
    tree_2.distance(*tree_1.box(i), candidates_query_bb);
    tree_2.stats(stats);
    checks_query_dt2b1 += stats.check_counter;
  }

  return std::make_tuple(
    0.5*(checks_query_it1t2 + checks_query_it2t1),
    0.5*(checks_query_it1b2 + checks_query_it2b1),
    0.5*(checks_query_dt1t2 + checks_query_dt2t1),
    0.5*(checks_query_dt1b2 + checks_query_dt2b1)
  );

}

// Main function
int
main() {

  #ifdef AABBTREE_ENABLE_PLOTTING
  Plot2D DD;
  DD.grid(true);
  DD.hold(true);
  DD.xlabel("Number of objects");
  DD.ylabel("Number of checks");
  DD.title("Distance");
  Plot2D II;
  II.grid(true);
  II.hold(true);
  II.xlabel("Number of objects");
  II.ylabel("Number of checks");
  II.title("Intersection");
  Plot2D RR;
  RR.grid(true);
  RR.hold(true);
  RR.xlabel("Number of objects");
  RR.ylabel("Number of checks");
  RR.title("Intersection");
  #endif

  // Set the number of asteroids for the two trees
  std::vector<std::string> colors = {"r", "g", "b", "c", "m", "y", "k"};
  std::vector<Integer>     n_objects = {50, 100, 200, 400, 800, 1875, 3750, 7500, 15000, 30000, 60000};
  std::vector<Integer>     n_nodal_objs = {1, 5, 10, 20, 40};
  std::vector<std::string> L{"Tree-Tree", "Tree-Boxes"};

  // Open files to save the results
  std::ofstream intersection("intersection.txt");
  std::ofstream distance("distance.txt");
  std::ofstream ratio("ratio.txt");

  std::cout << "----------------------------------------\n";
  std::vector<Real> checks_i_tt(n_objects.size());
  std::vector<Real> checks_i_tb(n_objects.size());
  std::vector<Real> checks_d_tt(n_objects.size());
  std::vector<Real> checks_d_tb(n_objects.size());
  std::vector<Real> checks_r_i(n_objects.size());
  std::vector<Real> checks_r_d(n_objects.size());
  std::vector<Real> n_objects_tt(n_objects.size());
  for (Integer j{0}; j < static_cast<Integer>(n_nodal_objs.size()); ++j) {
    for (Integer i{0}; i < static_cast<Integer>(n_objects.size()); ++i) {
      auto [i_tt, i_tb, d_tt, d_tb] = ComputeChecks(n_objects[i], n_objects[i], n_nodal_objs[j]);
      checks_i_tt[i] = i_tt;
      checks_i_tb[i] = i_tb;
      checks_d_tt[i] = d_tt;
      checks_d_tb[i] = d_tb;
      checks_r_i[i] = i_tb/i_tt;
      checks_r_d[i] = d_tb/d_tt;
      std::cout << "Number of asteroids: " << n_objects[i] << '\n';
      std::cout << "Number of objects:   " << n_nodal_objs[j] << '\n';
      std::cout << "Checks i_tt: " << i_tt << '\n';
      std::cout << "Checks i_tb: " << i_tb << '\n';
      std::cout << "Checks r_i: " << checks_r_i[i] << '\n';
      std::cout << "Checks d_tt: " << d_tt << '\n';
      std::cout << "Checks d_tb: " << d_tb << '\n';
      std::cout << "Checks r_d: " << checks_r_d[i] << '\n';
      std::cout << "----------------------------------------\n";
      intersection  << n_objects[i] << " " << n_nodal_objs[j] << " " << i_tt << " " << i_tb << '\n';
      distance      << n_objects[i] << " " << n_nodal_objs[j] << " " << d_tt << " " << d_tb << '\n';
      ratio         << n_objects[i] << " " << n_nodal_objs[j] << " " << i_tb/i_tt << " " << d_tb/d_tt << '\n';
    }
    #ifdef AABBTREE_ENABLE_PLOTTING
    II.loglog(n_objects, checks_i_tt)->line_width(1.5).color(colors[j]);
    II.loglog(n_objects, checks_i_tb)->line_width(1.5).color(colors[j]).line_style("--");
    II.legend(L);
    DD.loglog(n_objects, checks_d_tt)->line_width(1.5).color(colors[j]);
    DD.loglog(n_objects, checks_d_tb)->line_width(1.5).color(colors[j]).line_style("--");
    DD.legend(L);
    RR.loglog(n_objects, checks_r_i)->line_width(1.5).color(colors[j]);
    RR.loglog(n_objects, checks_r_d)->line_width(1.5).color(colors[j]).line_style("--");
    RR.legend(L);
    #endif
  }
  intersection.close();
  distance.close();
  ratio.close();

  #ifdef AABBTREE_ENABLE_PLOTTING
  II.show();
  DD.show();
  RR.show();
  #endif

  return 0;
}
