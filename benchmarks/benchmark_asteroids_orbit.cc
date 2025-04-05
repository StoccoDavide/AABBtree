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
int main()
{
  using Real = double;
  using Integer = AABBtree::Integer;

  constexpr Integer n_asteroids{60000};
  constexpr Integer n_clusters{60000};
  constexpr Integer n_clusters_to_plot{10};
  constexpr Real    t_ini{64328.0}; // January 1, 2035
  constexpr Integer n_neighbours{2};
  constexpr Real    t_end{t_ini + 10.0*365.0};
  constexpr Integer t_steps{10.0*365};
  #ifdef AABBTREE_ENABLE_PLOTTING
  constexpr Integer l_trace{365};
  #endif

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
  std::ofstream asteroids_trace("asteroids_phasing_traces.txt");
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
    // Phasing
    clusters_distance[i] = tree.closest(*tree.box(i), n_neighbours, clusters[i],
    [] (Box const & b1, Box const & b2) {
      Vector v1(b1.baricenter()), v2(b2.baricenter());
      Real distance{std::numeric_limits<Real>::max()}, dx, dy, dz;
      for (Integer j{0}; j < 3*t_steps; j += 3) {
        dx = v1(j+0) - v2(j+0);
        dy = v1(j+1) - v2(j+1);
        dz = v1(j+2) - v2(j+2);
        distance = std::min(distance, dx*dx + dy*dy + dz*dz);
      }
      return std::sqrt(distance);
    });
  }
  timer.toc();

  // Find the phasing points
  std::vector<Real> phasing_time(n_clusters);
  std::vector<Eigen::Vector3d> phasing_points_1(n_clusters);
  std::vector<Eigen::Vector3d> phasing_points_2(n_clusters);
  for (Integer i{0}; i < n_clusters; ++i) {
    phasing_points_1[i].setZero();
    phasing_points_2[i].setZero();
    phasing_time[i] = 0.0;
    Real dx, dy, dz, tmp, distance{std::numeric_limits<Real>::max()};
    for (Integer j1 : clusters[i]) {
    Vector v1(tree.box(j1)->baricenter());
      for (Integer j2 : clusters[i]) {
        if (j1 == j2) continue;
        Vector v2(tree.box(j2)->baricenter());
        for (Integer k{0}; k < 3*t_steps; k += 3) {
          dx = v1(k+0) - v2(k+0);
          dy = v1(k+1) - v2(k+1);
          dz = v1(k+2) - v2(k+2);
          tmp = dx*dx + dy*dy + dz*dz;
          if (tmp < distance) {
            distance = tmp;
            phasing_points_1[i] << v1(k+0), v1(k+1), v1(k+2);
            phasing_points_2[i] << v2(k+0), v2(k+1), v2(k+2);
            phasing_time[i] = t_ini + (k/3)*dt;
          }
        }
      }
    }
  }

  // Save the phasing points and time
  std::ofstream asteroids_phasing("asteroids_phasing_phasing.txt");
  for (Integer i{0}; i < n_clusters; ++i) {
    asteroids_phasing << i << " " << phasing_time[i] << " " << phasing_points_1[i].transpose() <<
    " " << phasing_points_2[i].transpose() << '\n';
  }
  asteroids_phasing.close();

  // Find the best clusters (minimum distance)
  std::vector<Integer> sorting(n_clusters);
  std::iota(sorting.begin(), sorting.end(), 0);
  std::sort(sorting.begin(), sorting.end(), [&clusters_distance, &clusters](Integer i, Integer j) {
    return clusters_distance[i] < clusters_distance[j] && clusters[i].size() == n_neighbours;
  });

  // Save the clusters
  std::ofstream asteroids_cluster("asteroids_phasing_clusters.txt");
  for (Integer i{0}; i < n_clusters; ++i) {
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
    Plot2D XY, XZ;
    XY.grid(true);
    XZ.grid(true);
      std::vector<Real> x(n_asteroids), y(n_asteroids), z(n_asteroids);
    Real max_xy{0.0};
    Real max_xz{0.0};
    for (Integer i{0}; i < n_asteroids; ++i) {
      KeplerianToCartesian(data[i], x[i], y[i], z[i]);
      max_xy = std::max(max_xy, std::max(std::abs(x[i]), std::abs(y[i])));
      max_xz = std::max(max_xz, std::max(std::abs(x[i]), std::abs(z[i])));
    }
    XY.hold(true);
    XZ.hold(true);
    XY.plot(x, y, ".")->marker_size(0.8).marker_color({0.5, 0.5, 0.5});
    XY.xlim(-max_xy, max_xy); XY.xlabel("x (AU)");
    XY.ylim(-max_xy, max_xy); XY.ylabel("y (AU)");
    XZ.plot(x, z, ".")->marker_size(0.8).marker_color({0.5, 0.5, 0.5});
    XZ.xlim(-max_xy, max_xy); XZ.xlabel("x (AU)");
    XZ.ylim(-max_xz, max_xz); XZ.ylabel("z (AU)");
  }
  #endif

  // Plot the clusters
  #ifdef AABBTREE_ENABLE_PLOTTING
  {
    Plot2D XY, XZ;
    XY.grid(true);
    XZ.grid(true);
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

      // Plot the phasing points
      std::vector<Real> X{phasing_points_1[sorting[i]].x(), phasing_points_2[sorting[i]].x()};
      std::vector<Real> Y{phasing_points_1[sorting[i]].y(), phasing_points_2[sorting[i]].y()};
      std::vector<Real> Z{phasing_points_1[sorting[i]].z(), phasing_points_2[sorting[i]].z()};
      XY.plot( X, Y, "r-" )->line_width(1.0);
      XZ.plot( X, Z, "r-" )->line_width(1.0);
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
