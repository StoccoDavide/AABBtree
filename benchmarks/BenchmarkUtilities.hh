/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The AABBtree project is distributed under the BSD 2-Clause License.                           *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef INCLUDE_BENCHMARK_UTILITIES_HH
#define INCLUDE_BENCHMARK_UTILITIES_HH

// C++17 standard libraries
#include <vector>
#include <algorithm>

// AABBtree library
#include "AABBtree.hh"

// Matplot++ library
#ifdef AABBTREE_ENABLE_PLOTTING
#include <matplot/matplot.h>
using namespace matplot;
#endif

// Benchmark utilities
namespace BenchmarkUtilities {

  using namespace AABBtree;

  /**
  * \brief Class container for a timer.
  *
  * Class container for a timer that uses the C++ standard library chrono to measure the elapsed time
  * between two points in the code.
  * \tparam Real Scalar number type.
  */
  template <typename Real>
  class TicToc {

    // Basic types definitions
    using Clock = std::chrono::high_resolution_clock; /**< Clock type. */
    using Time  = std::chrono::microseconds; /**< Time type. */

    // Timer variables
    Clock::time_point m_start_time;   /**< Start time.   */
    Clock::time_point m_stop_time;    /**< Stop time.    */
    Time              m_elapsed_time; /**< Elapsed time. */

  public:
    /**
     * Copy constructor for the timer.
     */
    TicToc(const TicToc &) = delete;

    /**
     * Assignment operator for the timer.
     */
    TicToc & operator=(TicToc const &) = delete;

    /**
     * Class constructor for the timer.
     */
    TicToc() : m_elapsed_time(0) {}

    /**
     * Start the timer.
     */
    void tic() {this->m_start_time = Clock::now();}

    /**
     * Stop the timer.
     */
    void
    toc() {
      this->m_stop_time    = Clock::now();
      this->m_elapsed_time = std::chrono::duration_cast<Time>(this->m_stop_time - this->m_start_time);
    }

    /**
     * Get the elapsed time in seconds.
     * \return The elapsed time in seconds.
     */
    Real elapsed_s() const {return static_cast<Real>(1.0e-6)*this->m_elapsed_time.count();}

    /**
     * Get the elapsed time in milliseconds.
     * \return The elapsed time in milliseconds.
     */
    Real elapsed_ms() const {return static_cast<Real>(1.0e-3)*this->m_elapsed_time.count();}

    /**
     * Get the elapsed time in microseconds.
     * \return The elapsed time in microseconds.
     */
    Real elapsed_us() const {return this->m_elapsed_time.count();}

  }; // TicToc

  // Class segment (2D)
  template <typename Real>
  class Segment {
    using Vector2 = AABBtree::Vector<Real, 2>;
    Vector2 m_p_1;
    Vector2 m_p_2;
  public:
    ~Segment() = default;
    Segment() = default;
    Segment(Segment const &) = default;
    Segment(Real t_x_1, Real t_y_1, Real t_x_2, Real t_y_2) : m_p_1(t_x_1, t_y_1), m_p_2(t_x_2, t_y_2) {}
    Segment(Vector2 const & p_1, Vector2 const & p_2) : m_p_1(p_1), m_p_2(p_2) {}

    Vector2 const & p_1() const {return this->m_p_1;}
    Vector2 const & p_2() const {return this->m_p_2;}
    Vector2 & p_1() {return this->m_p_1;}
    Vector2 & p_2() {return this->m_p_2;}

    Vector2 const &
    point(Integer i) const {
      if (i == 0) return this->m_p_1;
      if (i == 1) return this->m_p_2;
      throw std::out_of_range("Segment point index out of range.");
    }

    Vector2 &
    point(Integer i) {
      if (i == 0) return this->m_p_1;
      if (i == 1) return this->m_p_2;
      throw std::out_of_range("Segment point index out of range.");
    }

    Box<Real, 2>
    bounding_box() const {
      return Box<Real, 2>(
        std::min(this->m_p_1[0], this->m_p_2[0]), std::min(this->m_p_1[1], this->m_p_2[1]),
        std::max(this->m_p_1[0], this->m_p_2[0]), std::max(this->m_p_1[1], this->m_p_2[1])
      );
    }

    bool
    intersect( Segment const & segment, Vector2 & point) const {

      Real const & t_x_1{this->m_p_1[0]};   Real const & t_y_1{this->m_p_1[1]};
      Real const & t_x_2{this->m_p_2[0]};   Real const & t_y_2{this->m_p_2[1]};
      Real const & s_x_1{segment.m_p_1[0]}; Real const & s_y_1{segment.m_p_1[1]};
      Real const & s_x_2{segment.m_p_2[0]}; Real const & s_y_2{segment.m_p_2[1]};
      Real const d{(t_x_1 - t_x_2)*(s_y_1 - s_y_2) - (t_y_1 - t_y_2)*(s_x_1 - s_x_2)};

      if (std::abs(d) < std::numeric_limits<Real>::epsilon()) return false;

      Real const t{((t_x_1 - s_x_1)*(s_y_1 - s_y_2) - (t_y_1 - s_y_1)*(s_x_1 - s_x_2))/d};
      Real const u{-((t_x_1 - t_x_2)*(t_y_1 - s_y_1) - (t_y_1 - t_y_2)*(t_x_1 - s_x_1))/d};

      if (t >= 0.0 && t <= 1.0 && u >= 0.0 && u <= 1.0) {
        point.x() = t_x_1 + t*(t_x_2 - t_x_1);
        point.y() = t_y_1 + t*(t_y_2 - t_y_1);
        return true;
      }

      return false;
    }
  };

  } // Asteroids benchmark utilities

  static const double GM{0.0002959122082855911}; // Gravitational parameter of the Sun in AU^3/day^2

  // Keplerian orbital elements
  template <typename Real, typename Integer>
  struct Keplerian {
    Integer id;      // Asteroid ID
    Real    epoch;   // Epoch (Julian date)
    Real    a;       // Semi-major axis (AU)
    Real    e;       // Eccentricity
    Real    i;       // Inclination (rad)
    Real    lan;     // Longitude of ascending node (rad)
    Real    argperi; // Argument of periapsis (rad)
    Real    m;       // Mean anomaly (rad)
  };

  // Parse orbital data
  template <typename Real, typename Integer>
  bool
  Parse(
    std::string                     const & fname,
    std::vector<Keplerian<Real, Integer>> & data,
    Integer                         const & n,
    bool                                    reverse = false
 ) {
    // Clear data
    data.clear();
    data.reserve(n);

    // Open file
    std::ifstream file(fname);
    if (!file) { std::cerr << "Error opening file: " << fname << '\n'; return false;}

    std::vector<std::string> lines;
    std::string line;
    std::getline(file, line); // Skip first line
    while (std::getline(file, line)) {lines.push_back(line);}
    file.close();

    // Reverse lines if needed
    if (reverse) std::reverse(lines.begin(), lines.end());

    // Parse data
    Integer count{0};
    for (const auto & i_line : lines) {
      if (++count > n) {break;}
      std::istringstream iss(i_line);
      Keplerian<Real, Integer> entry;
      iss >> entry.id >> entry.epoch >> entry.a >> entry.e >> entry.i >> entry.lan >> entry.argperi >> entry.m;
      if (!iss) { std::cerr << "Error parsing line: " << i_line << '\n'; continue;}
      entry.i       *= M_PI/180.0;
      entry.lan     *= M_PI/180.0;
      entry.argperi *= M_PI/180.0;
      entry.m       *= M_PI/180.0;
      data.emplace_back(entry);
    }
    return true;
  }

  // Solve Kepler's equation to get eccentric anomaly
  template <typename Real, typename Integer>
  Real SolveKepler(Real M, Real e, Real tol = 1e-8, Integer max_iter = 100) {
    Real E{M};
    for (Integer i{0}; i < max_iter; ++i) {
      Real delta{E - e*std::sin(E) - M};
      if (std::abs(delta) < tol) {break;}
      E -= delta/(1.0 - e*std::cos(E));
      if (i == max_iter - 1) std::cerr << "Failed to converge\n";
    }
    return E;
  }

  // Convert Keplerian elements to Cartesian coordinates
  template <typename Real, typename Integer>
  void
  KeplerianToCartesian( Keplerian<Real, Integer> const & kepl, Real & x, Real & y, Real & z) {
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
  void
  PropagateOrbit( Keplerian<Real, int> & kepl, Real t_end, Real t_ini = 0.0) {
    if (t_end == t_ini) return;
    Real n{std::sqrt(static_cast<Real>(GM)/(kepl.a*kepl.a*kepl.a))}; // Mean motion (rad/day)
    kepl.m = std::fmod(kepl.m + n*(t_end - t_ini), 2.0*M_PI);  // Keep within [0,2π]
    if (kepl.m < 0) kepl.m += 2.0*M_PI;
  }

#endif // INCLUDE_BENCHMARK_UTILITIES_HH
