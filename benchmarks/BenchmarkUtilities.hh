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
    Clock::time_point m_start_time; /**< Start time. */
    Clock::time_point m_stop_time; /**< Stop time. */
    Time m_elapsed_time; /**< Elapsed time. */

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
    void toc()
    {
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

  }; // class TicToc

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
    Vector2 const & point(Integer i) const {
      if (i == 0) {return this->m_p_1;}
      if (i == 1) {return this->m_p_2;}
      throw std::out_of_range("Segment point index out of range.");
    }
    Vector2 & point(Integer i) {
      if (i == 0) {return this->m_p_1;}
      if (i == 1) {return this->m_p_2;}
      throw std::out_of_range("Segment point index out of range.");
    }
    Box<Real, 2> bounding_box() const {
      return Box<Real, 2>(
        std::min(this->m_p_1[0], this->m_p_2[0]), std::min(this->m_p_1[1], this->m_p_2[1]),
        std::max(this->m_p_1[0], this->m_p_2[0]), std::max(this->m_p_1[1], this->m_p_2[1])
      );;
    }
    bool intersect(Segment const & segment, Vector2 & point) const {
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

  #ifdef AABBTREE_ENABLE_PLOTTING

  // Plot a box
  template <typename Real, Integer N>
  void plot_box(Box<Real, N> const & box, std::string const & color, Real const line_width = 1) {
    static_assert(N == 2, "Plotting is only supported for 2D objects.");
    std::vector<Real> min(box.min().data(), box.min().data() + N);
    std::vector<Real> max(box.max().data(), box.max().data() + N);
    if (ax) {
      ax->plot({min[0], max[0], max[0], min[0], min[0]},
        {min[1], min[1], max[1], max[1], min[1]})->line_width(line_width).color(color);
    }
  }

  // Plot a tree
  template <typename Real, Integer N>
  void plot_tree(Tree<Real, N> const & tree, std::string const & color, Real const line_width = 1) {
    static_assert(N == 2, "Plotting is only supported for 2D objects.");
    for (Integer i{0}; i < static_cast<Integer>(tree.size()); ++i) {
      plot_box<Real, N>(tree.node(i).box, color, line_width);
      plot_box<Real, N>(tree.node(i).box_long, color, 2.0*line_width);
    }
  }

  // Plot a point
  template <typename Real, Integer N>
  void plot_point(AABBtree::Vector<Real, N> const & point, std::string const & color,
    Real const marker_size = 0.5) {
    static_assert(N == 2, "Plotting is only supported for 2D objects.");
    std::vector<Real> p(point.data(), point.data() + N);
    if (ax) {
      ax->plot({p[0]}, {p[1]}, "o")->marker_size(marker_size).color(color).marker_face(true);
    }
  }

  // Plot a segment
  template <typename Real, Integer N>
  void plot_segment(AABBtree::Vector<Real, N> const & p_1, AABBtree::Vector<Real, N> const & p_2,
    std::string const & color, Real const line_width = 1) {
    static_assert(N == 2, "Plotting is only supported for 2D objects.");
    std::vector<Real> v_1(p_1.data(), p_1.data() + N);
    std::vector<Real> v_2(p_2.data(), p_2.data() + N);
    if (ax) {
      ax->plot({v_1[0], v_2[0]}, {v_1[1], v_2[1]}, "-o")->color(color).line_width(line_width).
        marker_size(2.0*line_width).marker_face(true);
    }
  }

  // Plot a ray
  template <typename Real, Integer N>
  void plot_ray(AABBtree::Ray<Real, N> const & ray, std::string const & color, Real const line_width = 1) {
    static_assert(N == 2, "Plotting is only supported for 2D objects.");
    Real t{1000.0};
    plot_segment<Real, N>(ray.origin(), ray.origin() + t*ray.direction(), color, line_width);
  }

  // Plot a segment
  template <typename Real>
  void plot_segment(Segment<Real> const & segment, std::string const & color, Real const line_width = 1) {
    plot_segment<Real, 2>(segment.p_1(), segment.p_2(), color, line_width);
  }

  #endif // AABBTREE_ENABLE_PLOTTING

} // namespace BenchmarkUtilities

#endif // INCLUDE_BENCHMARK_UTILITIES_HH
