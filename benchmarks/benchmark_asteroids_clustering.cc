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
#include <sstream>
#include <iomanip>

// Matplot++ library
#ifdef AABBTREE_ENABLE_PLOTTING
#include "Plot2D.hh"
#endif

// Benchmark utilities
#include "BenchmarkUtilities.hh"
using namespace BenchmarkUtilities;

// Main function
int
main() {

  using Real    = double;
  using Integer = AABBtree::Integer;

  constexpr Integer n_asteroids{60000};
  constexpr Integer n_clusters{60000};
  constexpr Integer n_clusters_to_plot{10};
  constexpr Real    t_ini{64328.0}; // January 1, 2035
  constexpr Integer n_neighbours{10};
  constexpr Real    t_end{t_ini + 30.0};
  constexpr Integer t_steps{4*30};
  constexpr Integer l_trace{10};

  using Vector           = AABBtree::Vector<Real, 3*t_steps>;
  using Box              = AABBtree::Box<Real, 3*t_steps>;
  using Tree             = AABBtree::Tree<Real, 3*t_steps>;
  using BoxUniquePtrList = AABBtree::BoxUniquePtrList<Real, 3*t_steps>;
  using TicToc           = BenchmarkUtilities::TicToc<Real>;

  // Initialize the timer
  TicToc timer;

  // Parse asteroids data
  std::string fname{ "./benchmarks/asteroids.txt" }; // from build directory
  std::vector<Keplerian<Real, Integer>> data;
  timer.tic();
  if (!Parse<Real, Integer>(fname, data, n_asteroids)) {
    std::cerr << "Error parsing data\n";
    return 1;
  }
  timer.toc();
  std::cout << "Data parsed in " << timer.elapsed_us() << " us\n";

  // Set the initial and final times
  Real dt{(t_end - t_ini)/t_steps};

  // Write on a new file
  std::unique_ptr<BoxUniquePtrList> boxes{ std::make_unique<BoxUniquePtrList>() };
  boxes->reserve(n_asteroids);
  std::ofstream asteroids_trace("asteroids_clustering_traces.txt");
  for (Integer i{0}; i < n_asteroids; ++i) {
    asteroids_trace << i;
    Vector x, y, z;
    Keplerian<Real, Integer> data_i = data[i];
    for (Integer j{0}; j < t_steps; ++j) {
      PropagateOrbit(data_i, t_ini + j*dt, (j == 0) ? t_ini : t_ini + (j-1)*dt);
      KeplerianToCartesian(data_i, x(j), y(j), z(j));
      asteroids_trace << " " << x(j) << " " << y(j) << " " << z(j);
    }
    asteroids_trace << '\n';
    Vector box_min, box_max;
    constexpr Real tol_trace{1.0e-6};
    for (Integer j{0}; j < t_steps; ++j) {
      box_min(3*j+0) = x(j) - tol_trace;
      box_min(3*j+1) = y(j) - tol_trace;
      box_min(3*j+2) = z(j) - tol_trace;
      box_max(3*j+0) = x(j) + tol_trace;
      box_max(3*j+1) = y(j) + tol_trace;
      box_max(3*j+2) = z(j) + tol_trace;
    }
    boxes->push_back(std::make_unique<Box>(box_min, box_max));
    boxes->back()->reorder(); // Always wear a helmet
  }
  asteroids_trace.close();
  std::cout << "Boxes prepared in " << timer.elapsed_us() << " us\n";

  // Build the tree
  Tree tree;
  timer.tic(); tree.build(std::move(boxes)); timer.toc();
  std::cout << "Tree built in " << timer.elapsed_us() << " us\n";
  tree.print(std::cout);

  #ifdef AABBTREE_ENABLE_PLOTTING
  //{show(fig);}
  #endif

  // Find the closest clusters
  std::vector<IndexSet> clusters(n_clusters);
  std::vector<Real> clusters_distance(n_clusters);
  timer.tic();
  for (Integer i{0}; i < n_clusters; ++i) {
    if (i % 1000 == 0) {
      std::cout << "Cluster " << i << " of " << n_clusters << '\n';
    }
    // Clustering
    clusters_distance[i] = tree.closest(*tree.box(i), n_neighbours, clusters[i],
    [] (Box const & b1, Box const & b2) {
      Vector v1(b1.baricenter()), v2(b2.baricenter());
      Real distance{0.0}, tmp, dx, dy, dz;
      for (Integer i{0}; i < 3*t_steps; i += 3) {
        tmp = std::numeric_limits<Real>::max();
        for (Integer j{0}; j < 3*t_steps; j += 3) {
          dx = v1(i+0) - v2(j+0);
          dy = v1(i+1) - v2(j+1);
          dz = v1(i+2) - v2(j+2);
          tmp = std::min(tmp, dx*dx + dy*dy + dz*dz);
        }
        distance += std::sqrt(tmp);
      }
      return distance;
    });
  }
  timer.toc();

  // Find the best clusters (minimum distance)
  std::vector<Integer> sorting(n_clusters);
  std::iota(sorting.begin(), sorting.end(), 0);
  std::sort(sorting.begin(), sorting.end(), [&clusters_distance, &clusters](Integer i, Integer j) {
    return clusters_distance[i] < clusters_distance[j] && clusters[i].size() == n_neighbours;
  });

  // Save the clusters
  std::ofstream asteroids_cluster("asteroids_clustering_clusters.txt");
  for ( Integer i{0}; i < n_clusters; ++i ) {
    asteroids_cluster << i << " " << sorting[i] << " " << clusters_distance[sorting[i]];
    for (Integer j : clusters[sorting[i]]) {
      asteroids_cluster << " " << j;
    }
    asteroids_cluster << '\n';
  }
  asteroids_cluster.close();

  // Plot all the asteroids
  #ifdef AABBTREE_ENABLE_PLOTTING
  {
    Plot2D XY;
    Plot2D XZ;
    //plot_set_grid( ax );
    std::vector<Real> x(n_asteroids), y(n_asteroids), z(n_asteroids);
    Real mXY{0.0};
    Real mXZ{0.0};
    for (Integer i{0}; i < n_asteroids; ++i) {
      KeplerianToCartesian(data[i], x[i], y[i], z[i]);
      mXY = std::max(mXY, std::max(std::abs(x[i]), std::abs(y[i])));
      mXZ = std::max(mXZ, std::max(std::abs(x[i]), std::abs(z[i])));
    }
    XY.hold(true);
    XZ.hold(true);
    XY.plot(x, y, ".")->marker_size(0.8).marker_color({0.5, 0.5, 0.5});
    XY.xlim(-mXY, mXY); XY.xlabel("x (AU)");
    XY.ylim(-mXY, mXY); XY.ylabel("y (AU)");
    XZ.plot(x, z, ".")->marker_size(0.8).marker_color({0.5, 0.5, 0.5});
    XZ.xlim(-mXY, mXY); XZ.xlabel("x (AU)");
    XZ.ylim(-mXZ, mXZ); XZ.ylabel("z (AU)");
  }
  #endif

  // Plot the clusters
  #ifdef AABBTREE_ENABLE_PLOTTING
  {
    Plot2D XY;
    Plot2D XZ;
    std::vector<Real> x, y, z;
    for (Integer i{0}; i < n_clusters_to_plot; ++i) {
      x.clear(); y.clear(); z.clear();
      x.resize(n_neighbours); y.resize(n_neighbours); z.resize(n_neighbours);
      Integer j{0};
      for (Integer k : clusters[sorting[i]]) {
        KeplerianToCartesian(data[k], x.at(j), y.at(j), z.at(j)); ++j;
        std::vector<Real> x_trace(l_trace), y_trace(l_trace), z_trace(l_trace);
        Real dt_trace{(t_end - t_ini)/l_trace};
        Keplerian<Real, Integer> data_k = data[k];
        for (Integer l{0}; l < l_trace; ++l) {
          // Backward propagation
          PropagateOrbit(data_k, t_ini + l*dt_trace, (l == 0) ? t_ini : t_ini + (l-1)*dt_trace);
          KeplerianToCartesian(data_k, x_trace[l], y_trace[l], z_trace[l]);
        }
        XY.plot(x_trace, y_trace, "k")->line_width(0.5);
        XZ.plot(x_trace, z_trace, "k")->line_width(0.5);
      }
      XY.plot(x, y, ".")->marker_size(2.5);
      XZ.plot(x, z, ".")->marker_size(2.5);
    }
    XY.show();
    XZ.show();
  }
  #endif

  for (Integer i{0}; i < n_clusters_to_plot; ++i) {
    std::cout << "Cluster " << i << " distance: " << std::scientific <<
    clusters_distance[sorting[i]] << '\n';
  }
  std::cout << "Clusters found in " << timer.elapsed_us() << " us\n";
  std::cout << "Average time per cluster: " << timer.elapsed_us()/n_clusters << " us\n";

  return 0;
}
