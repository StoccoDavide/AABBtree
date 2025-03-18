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
  Integer id; // asteroid ID
  Real epoch; // epoch (Julian date)
  Real a; // semi-major axis (AU)
  Real e; // eccentricity
  Real i; // inclination (rad)
  Real lan; // longitude of ascending node (rad)
  Real argperi; // argument of periapsis (rad)
  Real m; // mean anomaly (rad)
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
    entry.i *= M_PI/180.0; // Convert to radians
    entry.lan *= M_PI/180.0; // Convert to radians
    entry.argperi *= M_PI/180.0; // Convert to radians
    entry.m *= M_PI/180.0; // Convert to radians
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
void Propagate(Keplerian<Real, int> & kepl, Real t, Real t0 = 0.0) {
  Real n{std::sqrt(static_cast<Real>(GM)/(kepl.a*kepl.a*kepl.a))}; // Mean motion (rad/day)
  kepl.m = std::fmod(kepl.m + n*(t - t0), 2.0*M_PI);  // Keep within [0,2π]
  if (kepl.m < 0) {kepl.m += 2.0*M_PI;}
}

// Main function
int main() {
  using Real = double;
  using Integer = AABBtree::Integer;
  using Vector = AABBtree::Vector<Real, 3>;
  using Box = AABBtree::Box<Real, 3>;
  using Tree = AABBtree::Tree<Real, 3>;
  using BoxUniquePtrList = AABBtree::BoxUniquePtrList<Real, 3>;
  using TicToc = BenchmarkUtilities::TicToc<Real>;

  // Initialize the timer
  TicToc timer;

  // Parse asteroids data
  std::string fname = "./../benchmarks/asteroids.txt"; // from build directory
  Integer n{60};
  std::vector<Keplerian<Real, Integer>> data;
  timer.tic(); Parse<Real, Integer>(fname, data, n); timer.toc();
  std::cout << "Data parsed in " << timer.elapsed_us() << " us" << std::endl;

  // Set the initial and final times
  Real t_ini{64328.0}; // January 1, 2035
  Real t_end{t_ini + 10.0};

  // Plot an asteroid and its evolution in time [0, 1month]
  #ifdef AABBTREE_ENABLE_PLOTTING
  {
    SET_PLOT
    Integer id{0}, steps{100};
    std::vector<Real> t(steps), x(steps), y(steps), z(steps);
    Real max_xy{0.0};
    Real t_step{(t_end - t_ini)/(steps - 1)};
    KeplerianToCartesian(data[id], x[0], y[0], z[0]);
    for (Integer i{1}; i < steps; ++i) {
      t[i] = t[i-1] + t_step;
      Propagate(data[id], t[i], t[i-1]);
      KeplerianToCartesian(data[id], x[i], y[i], z[i]);
      max_xy = std::max(max_xy, std::max(std::abs(x[i]), std::abs(y[i])));
    }
    ax->hold(true);
    ax->plot(x, y, "-o");
    ax->xlim({-max_xy, max_xy}); ax->xlabel("x (AU)");
    ax->ylim({-max_xy, max_xy}); ax->ylabel("y (AU)");
    ax->hold(false);
    show(fig);
  }
  #endif

  // Plot the asteroids
  #ifdef AABBTREE_ENABLE_PLOTTING
  {
    SET_PLOT
    std::vector<Real> x(n), y(n), z(n);
    Real max_xy{0.0};
    Real max_xz{0.0};
    for (Integer i{0}; i < n; ++i) {
      KeplerianToCartesian(data[i], x[i], y[i], z[i]);
      max_xy = std::max(max_xy, std::max(std::abs(x[i]), std::abs(y[i])));
      max_xz = std::max(max_xz, std::max(std::abs(x[i]), std::abs(z[i])));
    }
    ax_xy->hold(true);
    ax_xy->plot(x, y, "o");
    ax_xy->xlim({-max_xy, max_xy}); ax_xy->xlabel("x (AU)");
    ax_xy->ylim({-max_xy, max_xy}); ax_xy->ylabel("y (AU)");
    ax_xy->hold(false);
    ax_xz->hold(true);
    ax_xz->plot(x, z, "o");
    ax_xz->xlim({-max_xy, max_xy}); ax_xz->xlabel("x (AU)");
    ax_xz->ylim({-max_xz, max_xz}); ax_xz->ylabel("z (AU)");
    ax_xz->hold(false);
    show(fig_xy);
    show(fig_xz);
  }
  #endif

  // Propagate asteroids over time
  std::vector<Keplerian<Real, Integer>> data_new(data);
  timer.tic();
  for (Integer i{0}; i < n; ++i) {Propagate(data_new[i], t_ini, t_end);}
  timer.toc();

  // Plot the propagated asteroids
  #ifdef AABBTREE_ENABLE_PLOTTING
  {
    std::vector<Real> x_new(n), y_new(n), z_new(n);
    for (Integer i{0}; i < n; ++i) {
      KeplerianToCartesian(data_new[i], x_new[i], y_new[i], z_new[i]);
    }
    ax_xy->hold(true);
    ax_xy->plot(x_new, y_new, "ro");
    ax_xy->hold(false);
    ax_xz->hold(true);
    ax_xz->plot(x_new, z_new, "ro");
    ax_xz->hold(false);
    show(fig_xy);
    show(fig_xz);
  }
  #endif

  // Prepare the boxes
  std::unique_ptr<BoxUniquePtrList> boxes = std::make_unique<BoxUniquePtrList>();
  boxes->reserve(n);
  Real x_min, y_min, z_min, x_max, y_max, z_max;
  //Real a_tol{0.1}, a_cnt{3.0}, a_min{a_cnt-a_tol}, a_max{a_cnt+a_tol};
  //Real e_tol{0.1}, e_cnt{0.1}, e_min{e_cnt-e_tol}, e_max{e_cnt+e_tol};
  //Real i_tol{0.1}, i_cnt{0.1}, i_min{i_cnt-i_tol}, i_max{i_cnt+i_tol};
  //Real lan_tol{0.1}, lan_cnt{0.1}, lan_min{lan_cnt-lan_tol}, lan_max{lan_cnt+lan_tol};
  //Real argperi_tol{0.1}, argperi_cnt{0.1}, argperi_min{argperi_cnt-argperi_tol}, argperi_max{argperi_cnt+argperi_tol};
  //Real m_tol{10e5}, m_cnt{0.1}, m_min{m_cnt-m_tol}, m_max{m_cnt+m_tol};
  for (Integer i{0}; i < n; ++i) {
    KeplerianToCartesian(data[i], x_min, y_min, z_min);
    Propagate(data[i], t_ini, t_end);
    KeplerianToCartesian(data[i], x_max, y_max, z_max);
    boxes->push_back(std::make_unique<Box>(
      Vector(x_min, y_min, z_min), //, a_min, e_min, i_min, lan_min, argperi_min, m_min))
      Vector(x_max, y_max, z_max) //, a_max, e_max, i_max, lan_max, argperi_max, m_max),
    ));
    boxes->back()->reorder();
  }
  std::cout << "Boxes prepared in " << timer.elapsed_us() << " us" << std::endl;

  // Build the tree
  Tree tree;
  timer.tic(); tree.build(std::move(boxes)); timer.toc();
  std::cout << "Tree built in " << timer.elapsed_us() << " us" << std::endl;
  tree.print(std::cout);

  // Find the clusters
  IndexSet cluster;
  Box cluster_box(-5.0, -5.0, -5.0, 5.0, 5.0, 5.0);
  tree.intersect(cluster_box, cluster);
  std::cout << "Cluster size: " << cluster.size() << std::endl;

  return 0;
}
