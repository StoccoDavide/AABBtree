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
#include <matplot/matplot.h>
using namespace matplot;
static auto fig{figure(true)};
static axes_handle ax{fig->current_axes()};
static auto fig_xy{figure(true)};
static axes_handle ax_xy{fig_xy->current_axes()};
static auto fig_xz{figure(true)};
static axes_handle ax_xz{fig_xz->current_axes()};

#ifndef SET_PLOT
#define SET_PLOT \
grid(ax, true); \
grid(ax_xy, true); \
grid(ax_xz, true);
#endif
#endif

// Benchmark utilities
#include "BenchmarkUtilities.hh"
using namespace BenchmarkUtilities;

static const double GM{0.0002959122082855911}; // Gravitational parameter of the Sun in AU^3/day^2

// Keplerian orbital elements
template <typename Real, typename Integer>
struct Keplerian {
  Integer id;   // Asteroid ID
  Real epoch;   // Epoch (Julian date)
  Real a;       // Semi-major axis (AU)
  Real e;       // Eccentricity
  Real i;       // Inclination (rad)
  Real lan;     // Longitude of ascending node (rad)
  Real argperi; // Argument of periapsis (rad)
  Real m;       // Mean anomaly (rad)
};

// Parse orbital data
template <typename Real, typename Integer>
void Parse(std::string const & fname, std::vector<Keplerian<Real, Integer>> & data, Integer const & n) {
  // Clear data
  data.clear();
  data.reserve(n);

  // Open file
  std::ifstream file(fname);
  if (!file) {std::cerr << "Error opening file: " << fname << std::endl; return;}

  // Skip header line
  std::string line;
  std::getline(file, line);

  // Parse data
  Integer count{0};
  while (std::getline(file, line) && count < n) {
    count++;
    std::istringstream iss(line);
    Keplerian<Real, Integer> entry;
    if (!(iss >> entry.id >> entry.epoch >> entry.a >> entry.e >> entry.i >> entry.lan >> entry.argperi >> entry.m))
    {std::cerr << "Error parsing line: " << line << std::endl; continue;}
    entry.i       *= M_PI/180.0;
    entry.lan     *= M_PI/180.0;
    entry.argperi *= M_PI/180.0;
    entry.m       *= M_PI/180.0;
    data.emplace_back(entry);
  }
}

// Solve Kepler's equation to get eccentric anomaly
template <typename Real, typename Integer>
Real SolveKepler(Real M, Real e, Real tol = 1e-6, Integer max_iter = 100) {
  Real E{M};
  for (Integer i{0}; i < max_iter; ++i) {
    Real delta{E - e*std::sin(E) - M};
    if (std::abs(delta) < tol) {break;}
    E -= delta/(1.0 - e*std::cos(E));
    if (i == max_iter - 1) {std::cerr << "Failed to converge" << std::endl;}
  }
  return E;
}

// Convert Keplerian elements to Cartesian coordinates
template <typename Real, typename Integer>
void KeplerianToCartesian(Keplerian<Real, Integer> const & kepl, Real & x, Real & y, Real & z) {
  Real E{SolveKepler<Real, Integer>(kepl.m, kepl.e)};
  Real nu{2.0*std::atan2(std::sqrt(1.0 + kepl.e)*std::sin(E/2), std::sqrt(1.0 - kepl.e)*std::cos(E/2))};
  Real r{kepl.a*(1.0 - kepl.e*std::cos(E))};
  Real cos_lan{std::cos(kepl.lan)}, sin_lan{std::sin(kepl.lan)};
  Real cos_argperi{std::cos(kepl.argperi)}, sin_argperi{std::sin(kepl.argperi)};
  Real cos_i{std::cos(kepl.i)}, sin_i{std::sin(kepl.i)};
  Real cos_nu{std::cos(nu)}, sin_nu{std::sin(nu)};
  Real xp{r*cos_nu}, yp{r*sin_nu};
  x = (cos_lan*cos_argperi - sin_lan*sin_argperi*cos_i)*xp + (-cos_lan*sin_argperi - sin_lan*cos_argperi*cos_i)*yp;
  y = (sin_lan*cos_argperi + cos_lan*sin_argperi*cos_i)*xp + (-sin_lan*sin_argperi + cos_lan*cos_argperi*cos_i)*yp;
  z = (sin_argperi*sin_i)*xp + (cos_argperi*sin_i)*yp;
}

// Propagate orbital elements over time
template <typename Real>
void Propagate(Keplerian<Real, int> & kepl, Real t_end, Real t_ini = 0.0) {
  Real n{std::sqrt(static_cast<Real>(GM)/(kepl.a*kepl.a*kepl.a))}; // Mean motion (rad/day)
  kepl.m = std::fmod(kepl.m + n*(t_ini - t_end), 2.0*M_PI);  // Keep within [0,2π]
  if (kepl.m < 0) {kepl.m += 2.0*M_PI;}
}

// Main function
int main() {

  using Real = double;
  using Integer = AABBtree::Integer;

  // Constants
  constexpr Integer n_asteroids{60000};
  constexpr Integer n_clusters{100};
  constexpr Integer n_neighbours{20};
  constexpr Real t_ini{64328.0}; // January 1, 2035
  constexpr Real t_end{t_ini + 30.0};
  constexpr Integer t_steps{10};

  using Vector = AABBtree::Vector<Real, 3*t_steps>;
  using Box = AABBtree::Box<Real, 3*t_steps>;
  using Tree = AABBtree::Tree<Real, 3*t_steps>;
  using BoxUniquePtrList = AABBtree::BoxUniquePtrList<Real, 3*t_steps>;
  using TicToc = BenchmarkUtilities::TicToc<Real>;

  // Initialize the timer
  TicToc timer;

  // Parse asteroids data
  std::string fname = "./../benchmarks/asteroids.txt"; // from build directory
  std::vector<Keplerian<Real, Integer>> data;
  timer.tic(); Parse<Real, Integer>(fname, data, n_asteroids); timer.toc();
  std::cout << "Data parsed in " << timer.elapsed_us() << " us" << std::endl;

  // Set the initial and final times
  Real dt{(t_end - t_ini)/t_steps};

  // Prepare the boxes
  std::unique_ptr<BoxUniquePtrList> boxes = std::make_unique<BoxUniquePtrList>();
  boxes->reserve(n_asteroids);
  for (Integer i{0}; i < n_asteroids; ++i) {
    Vector x, y, z;
    Keplerian<Real, Integer> data_i = data[i];
    KeplerianToCartesian(data_i, x(0), y(0), z(0));
    for (Integer j{1}; j < t_steps; ++j) {
      Propagate(data_i, t_ini + j*dt, t_ini + (j-1)*dt);
      KeplerianToCartesian(data_i, x(j), y(j), z(j));
    }
    Vector box_min, box_max;
    for (Integer j{0}; j < 3*t_steps; j += 3) {
      box_min[j+0] = std::min({x(j), x(j), x(j)});
      box_min[j+1] = std::min({y(j), y(j), y(j)});
      box_min[j+2] = std::min({z(j), z(j), z(j)});
      box_max[j+0] = std::max({x(j), x(j), x(j)});
      box_max[j+1] = std::max({y(j), y(j), y(j)});
      box_max[j+2] = std::max({z(j), z(j), z(j)});
    }
    boxes->push_back(std::make_unique<Box>(box_min, box_max));
    boxes->back()->reorder(); // Always wear a helmet
  }
  std::cout << "Boxes prepared in " << timer.elapsed_us() << " us" << std::endl;

  // Build the tree
  Tree tree;
  timer.tic(); tree.build(std::move(boxes)); timer.toc();
  std::cout << "Tree built in " << timer.elapsed_us() << " us" << std::endl;
  tree.print(std::cout);

  std::vector<IndexSet> clusters(n_clusters);
  std::vector<Real> clusters_distance(n_clusters);
  timer.tic();
  for (Integer i{0}; i < n_clusters; ++i) {
    clusters_distance[i] = tree.closest(tree.box(i)->baricenter(), n_neighbours, clusters[i]);
      //[] (Vector const & p, Box const & b) {return b.interior_distance(p);});
  }
  timer.toc();

  // Find the best clusters (minimum distance)
  std::vector<Integer> sorting(n_clusters);
  std::iota(sorting.begin(), sorting.end(), 0);
  std::sort(sorting.begin(), sorting.end(), [&clusters_distance](Integer i, Integer j) {
    return clusters_distance[i] < clusters_distance[j];
  });

  // Plot all the asteroids
  #ifdef AABBTREE_ENABLE_PLOTTING
  {
    SET_PLOT
    std::vector<Real> x(n_asteroids), y(n_asteroids), z(n_asteroids);
    Real max_xy{0.0};
    Real max_xz{0.0};
    for (Integer i{0}; i < n_asteroids; ++i) {
      KeplerianToCartesian(data[i], x[i], y[i], z[i]);
      max_xy = std::max(max_xy, std::max(std::abs(x[i]), std::abs(y[i])));
      max_xz = std::max(max_xz, std::max(std::abs(x[i]), std::abs(z[i])));
    }
    ax_xy->hold(true);
    ax_xz->hold(true);
    ax_xy->plot(x, y, ".")->marker_size(0.8).marker_color({0.5, 0.5, 0.5});
    ax_xy->xlim({-max_xy, max_xy}); ax_xy->xlabel("x (AU)");
    ax_xy->ylim({-max_xy, max_xy}); ax_xy->ylabel("y (AU)");
    ax_xz->plot(x, z, ".")->marker_size(0.8).marker_color({0.5, 0.5, 0.5});
    ax_xz->xlim({-max_xy, max_xy}); ax_xz->xlabel("x (AU)");
    ax_xz->ylim({-max_xz, max_xz}); ax_xz->ylabel("z (AU)");
  }
  #endif

  // Plot the clusters
  #ifdef AABBTREE_ENABLE_PLOTTING
  {
    SET_PLOT
    std::vector<Real> x, y, z;
    for (Integer i{0}; i < 10; ++i) {
      x.resize(n_neighbours); y.resize(n_neighbours); z.resize(n_neighbours);
      auto & cluster = clusters[sorting[i]];
      Integer j{0};
      for (auto k : cluster) {
        KeplerianToCartesian(data[k], x[j], y[j], z[j]); ++j;
        Integer l_trace{20};
        std::vector<Real> x_trace(l_trace), y_trace(l_trace), z_trace(l_trace);
        Real dt_trace{(t_end - t_ini)/l_trace};
        Keplerian<Real, Integer> data_k_new = data[k];
        KeplerianToCartesian(data_k_new, x_trace[0], y_trace[0], z_trace[0]);
        for (Integer l{1}; l < l_trace; ++l) {
          Propagate(data_k_new, t_ini + l*dt_trace, t_ini + (l-1)*dt_trace);
          KeplerianToCartesian(data_k_new, x_trace[l], y_trace[l], z_trace[l]);
        }
        ax_xy->plot(x_trace, y_trace, "k")->line_width(0.5);
        ax_xz->plot(x_trace, z_trace, "k")->line_width(0.5);
      }
      ax_xy->plot(x, y, ".")->marker_size(2.5);
      ax_xz->plot(x, z, ".")->marker_size(2.5);
    }
    show(fig_xy);
    show(fig_xz);
  }
  #endif

  for (Integer i{0}; i < 10; ++i) {
    std::cout << "Cluster " << i << " distance: " << std::scientific <<
    clusters_distance[sorting[i]] << std::endl;
  }
  std::cout << "Clusters found in " << timer.elapsed_us() << " us" << std::endl;


  return 0;
}
