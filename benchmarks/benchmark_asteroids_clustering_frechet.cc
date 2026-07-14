/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The AABBtree project is distributed under the BSD 2-Clause License.       *
 *                                                                           *
 * Davide Stocco Enrico Bertolazzi                                           *
 * University of Trento University of Trento                                 *
 * davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// C++17 standard libraries
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

// Matplot++ library
#ifdef AABBTREE_ENABLE_PLOTTING
#include "Plot2D.hh"
#endif

// Benchmark utilities
#include "BenchmarkUtilities.hh"
using namespace BenchmarkUtilities;

// Main function
int main() {

  using Real = double;
  using Integer = AABBtree::Integer;

  constexpr Integer n_asteroids{60000};
  constexpr Integer n_clusters{60000};
  constexpr Integer n_clusters_to_plot{10};
  constexpr Real t_ini{64328.0}; // January 1, 2035
  constexpr Integer n_neighbours{10};
  constexpr Real t_end{t_ini + 30.0};
  constexpr Integer t_steps{4 * 30};
#ifdef AABBTREE_ENABLE_PLOTTING
  constexpr Integer l_trace{10};
#endif

  using Vector = AABBtree::Vector<Real, 3 * (t_steps - 1)>;
  using Box = AABBtree::Box<Real, 3 * (t_steps - 1)>;
  using Tree = AABBtree::Tree<Real, 3 * (t_steps - 1)>;
  using BoxUniquePtrList = AABBtree::BoxUniquePtrList<Real, 3 * (t_steps - 1)>;
  using TicToc = BenchmarkUtilities::TicToc<Real>;

  // Initialize the timer
  TicToc timer;

  // Parse asteroids data
  std::string fname{"./../benchmarks/asteroids.txt"}; // from build directory
  std::vector<Keplerian<Real, Integer>> data;
  timer.tic();
  if (!Parse<Real, Integer>(fname, data, n_asteroids)) {
    std::cerr << "Error parsing data\n";
    return 1;
  }
  timer.toc();
  std::cout << "Data parsed in " << timer.elapsed_us() << " us\n";

  // Set the initial and final times
  Real dt{(t_end - t_ini) / t_steps};

  // Write on a new file
  std::unique_ptr<BoxUniquePtrList> boxes{std::make_unique<BoxUniquePtrList>()};
  boxes->reserve(n_asteroids - 1);
  std::ofstream asteroids_trace("asteroids_clustering_traces_frechet.txt");
  for (Integer i{0}; i < n_asteroids; ++i) {
    if (i % 1000 == 0) {
      std::cout << "Asteroid " << i << " of " << n_asteroids << std::endl;
    }
    asteroids_trace << i;
    Vector x, y, z;
    Keplerian<Real, Integer> data_i = data[i];
    for (Integer j{0}; j < t_steps; ++j) {
      PropagateOrbit(data_i, t_ini + j * dt,
                     (j == 0) ? t_ini : t_ini + (j - 1) * dt);
      KeplerianToCartesian(data_i, x(j), y(j), z(j));
      asteroids_trace << " " << x(j) << " " << y(j) << " " << z(j);
    }
    asteroids_trace << std::endl;
    Vector box_min, box_max;
    constexpr Real tol_trace{1.0e-6};
    for (Integer j{0}; j < t_steps - 1; ++j) {
      box_min(3 * j + 0) = std::min(x(j), x(j + 1)) - tol_trace;
      box_min(3 * j + 1) = std::min(y(j), y(j + 1)) - tol_trace;
      box_min(3 * j + 2) = std::min(z(j), z(j + 1)) - tol_trace;
      box_max(3 * j + 0) = std::max(x(j), x(j + 1)) + tol_trace;
      box_max(3 * j + 1) = std::max(y(j), y(j + 1)) + tol_trace;
      box_max(3 * j + 2) = std::max(z(j), z(j + 1)) + tol_trace;
    }
    boxes->push_back(std::make_unique<Box>(box_min, box_max));
    boxes->back()->reorder(); // Always wear a helmet
  }
  asteroids_trace.close();
  std::cout << "Boxes prepared in " << timer.elapsed_us() << " us\n";

  // Build the tree
  Tree tree;
  tree.max_dumpings(100);
  tree.max_nodal_objects(50);
  timer.tic();
  tree.build(std::move(boxes));
  timer.toc();
  std::cout << "Tree built in " << timer.elapsed_us() << " us\n";

  // Find the closest clusters
  std::vector<IndexSet> clusters(n_clusters);
  std::vector<Real> clusters_distance(n_clusters);
  // auto frechet_distance_fun = [](Box const &b1, Box const &b2) {
  auto frechet_distance_fun = [](Box const &b1, Box const &b2) -> Real {
    // The input objects b1 and b2 are not single 3D boxes.
    //
    // They are interpreted as complete box trajectories:
    //
    //   b1 = (B1_0, B1_1, ..., B1_{t_steps-1})
    //   b2 = (B2_0, B2_1, ..., B2_{t_steps-1})
    //
    // The k-th box of each trajectory is reconstructed from its minimum
    // and maximum corners:
    //
    //   B1_k = Box(b1.min()(k), b1.max()(k))
    //   B2_k = Box(b2.min()(k), b2.max()(k))

    using Box3 = AABBtree::Box<Real, 3>;

    Vector const &b1_min{b1.min()};
    Vector const &b1_max{b1.max()};
    Vector const &b2_min{b2.min()};
    Vector const &b2_max{b2.max()};

    // Number of boxes in the two trajectories.
    //
    // Here both trajectories are assumed to have the same number of
    // time steps. If later you need two different lengths, replace this
    // by m and n.
    Integer const m{t_steps - 1};
    Integer const n{t_steps - 1};

    // Dynamic programming state:
    //
    // Let D[i][j] be the discrete Fréchet distance between the prefixes
    //
    //   B1_0, ..., B1_i
    //
    // and
    //
    //   B2_0, ..., B2_j.
    //
    // The recurrence is:
    //
    //   D[i][j] =
    //     max(
    //       distance(B1_i, B2_j),
    //       min(D[i-1][j], D[i][j-1], D[i-1][j-1])
    //     )
    //
    // We only store two rows:
    //
    //   previous[j] = D[i-1][j]
    //   current[j]  = D[i][j]
    //
    // Therefore the memory cost is O(t_steps), not O(t_steps^2).
    std::vector<Real> previous(n);
    std::vector<Real> current(n);

    // Reconstruct the first boxes of the two trajectories.
    Box3 const B1_0(b1_min.segment<3>(0), b1_max.segment<3>(0));
    Box3 const B2_0(b2_min.segment<3>(0), b2_max.segment<3>(0));

    // Base case:
    //
    // D[0][0] is simply the distance between the first two boxes.
    previous[0] = B1_0.interior_distance(B2_0);

    // First row initialization:
    //
    // D[0][j] corresponds to matching the first box B1_0 against
    // the prefix:
    //
    //   B2_0, ..., B2_j
    //
    // Since B1 cannot advance, the only admissible path is horizontal:
    //
    //   (0,0) -> (0,1) -> ... -> (0,j)
    //
    // The Fréchet cost is the maximum box-box distance encountered
    // along that path.
    for (Integer j{1}; j < n; ++j) {
      Box3 const B2_j(b2_min.segment<3>(j), b2_max.segment<3>(j));
      Real const d = B1_0.interior_distance(B2_j);
      previous[j] = std::max(previous[j - 1], d);
    }

    // Fill the remaining rows of the DP table.
    for (Integer i{1}; i < m; ++i) {

      // Reconstruct the current box of the first trajectory once per row.
      Box3 const B1_i(b1_min.segment<3>(i), b1_max.segment<3>(i));

      // First column initialization:
      //
      // D[i][0] corresponds to matching the prefix:
      //
      //   B1_0, ..., B1_i
      //
      // against the first box B2_0.
      //
      // Since B2 cannot advance, the only admissible path is vertical:
      //
      //   (0,0) -> (1,0) -> ... -> (i,0)
      current[0] = std::max(previous[0], B1_i.interior_distance(B2_0));

      // Fill the inner cells.
      for (Integer j{1}; j < n; ++j) {

        // Reconstruct the current box of the second trajectory.
        Box3 const B2_j(b2_min.segment<3>(j), b2_max.segment<3>(j));

        // Elementary distance between the current pair of boxes.
        //
        // This is the only place where the actual box-box distance is
        // used. The Fréchet structure is handled by the dynamic
        // programming recurrence below.
        Real const d = B1_i.interior_distance(B2_j);

        // A monotone Fréchet coupling can reach cell (i,j) from only
        // three predecessor cells:
        //
        //   (i-1, j)     : advance in trajectory 1 only
        //   (i,   j-1)   : advance in trajectory 2 only
        //   (i-1, j-1)   : advance in both trajectories
        //
        // Among these three options, choose the one with the smallest
        // previous Fréchet cost.
        Real const best_previous =
            std::min(previous[j], std::min(current[j - 1], previous[j - 1]));

        // Fréchet recurrence:
        //
        // The cost of a coupling is the maximum distance encountered
        // along the whole monotone path. Therefore the value at (i,j)
        // is the maximum between:
        //
        //   1. the current box-box distance d,
        //   2. the best previous path cost.
        current[j] = std::max(d, best_previous);
      }

      // Move the current row into previous.
      //
      // After the swap:
      //
      //   previous[j] contains D[i][j]
      //
      // for the row just computed.
      std::swap(previous, current);
    }

    // The final value is D[t_steps-1][t_steps-1], namely the discrete
    // Fréchet distance between the two complete box trajectories.
    return previous[n - 1];
  };

  timer.tic();
  for (Integer i{0}; i < n_clusters; ++i) {
    clusters_distance[i] = tree.closest(*tree.box(i), n_neighbours, clusters[i],
                                        frechet_distance_fun);
    if (i % 100 == 0) {
      timer.toc();
      std::cout << "Cluster " << i << " of " << n_clusters << std::endl;
      Real estimated_time = timer.elapsed_s() / 100 * (n_clusters - i) / 60;
      std::cout << "Estimated time: " << std::setprecision(2) << std::fixed
                << estimated_time << " min\n";
      timer.tic();
    }
  }
  timer.toc();

  // Find the best clusters (minimum distance)
  std::vector<Integer> sorting(n_clusters);
  std::iota(sorting.begin(), sorting.end(), 0);
  std::sort(sorting.begin(), sorting.end(),
            [&clusters_distance, &clusters](Integer i, Integer j) {
              return clusters_distance[i] < clusters_distance[j] &&
                     clusters[i].size() == n_neighbours;
            });

  // Save the clusters
  std::ofstream asteroids_cluster("asteroids_clustering_clusters_frechet.txt");
  for (Integer i{0}; i < n_clusters; ++i) {
    asteroids_cluster << i << " " << sorting[i] << " "
                      << clusters_distance[sorting[i]];
    for (Integer j : clusters[sorting[i]]) {
      asteroids_cluster << " " << j;
    }
    asteroids_cluster << std::endl;
  }
  asteroids_cluster.close();

// Plot all the asteroids
#ifdef AABBTREE_ENABLE_PLOTTING
  Plot2D XY;
  Plot2D XZ;
  {
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
    XY.xlim(-mXY, mXY);
    XY.xlabel("x (AU)");
    XY.ylim(-mXY, mXY);
    XY.ylabel("y (AU)");
    XZ.plot(x, z, ".")->marker_size(0.8).marker_color({0.5, 0.5, 0.5});
    XZ.xlim(-mXY, mXY);
    XZ.xlabel("x (AU)");
    XZ.ylim(-mXZ, mXZ);
    XZ.ylabel("z (AU)");
  }
#endif

// Plot the clusters
#ifdef AABBTREE_ENABLE_PLOTTING
  {
    std::vector<Real> x, y, z;
    for (Integer i{0}; i < std::min(n_clusters, n_clusters_to_plot); ++i) {
      x.clear();
      y.clear();
      z.clear();
      x.resize(n_neighbours);
      y.resize(n_neighbours);
      z.resize(n_neighbours);
      Integer j{0};
      for (Integer k : clusters[sorting[i]]) {
        KeplerianToCartesian(data[k], x.at(j), y.at(j), z.at(j));
        ++j;
        std::vector<Real> x_trace(l_trace), y_trace(l_trace), z_trace(l_trace);
        Real dt_trace{(t_end - t_ini) / l_trace};
        Keplerian<Real, Integer> data_k = data[k];
        for (Integer l{0}; l < l_trace; ++l) {
          // Backward propagation
          PropagateOrbit(data_k, t_ini + l * dt_trace,
                         (l == 0) ? t_ini : t_ini + (l - 1) * dt_trace);
          KeplerianToCartesian(data_k, x_trace[l], y_trace[l], z_trace[l]);
        }
        XY.plot(x_trace, y_trace, "k")->line_width(0.5);
        XZ.plot(x_trace, z_trace, "k")->line_width(0.5);
      }
      XY.plot(x, y, ".")->marker_size(2.5);
      XZ.plot(x, z, ".")->marker_size(2.5);
    }
  }
  XY.show();
  XZ.show();
#endif

  for (Integer i{0}; i < n_clusters_to_plot; ++i) {
    std::cout << "Cluster " << i << " distance: " << std::scientific
              << clusters_distance[sorting[i]] << std::endl;
  }

  std::cout << "Clusters found in " << timer.elapsed_us() << " us\n";
  std::cout << "Average time per cluster: " << timer.elapsed_us() / n_clusters
            << " us\n";

  return 0;
}
