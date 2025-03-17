/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The AABBtree project is distributed under the BSD 2-Clause License.                           *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef INCLUDE_TEST_UTILITIES_HH
#define INCLUDE_TEST_UTILITIES_HH

// C++17 standard libraries
#include <vector>

// AABBtree library
#include "AABBtree.hh"

// Matplot++ library
#ifdef AABBTREE_ENABLE_PLOTTING
#include <matplot/matplot.h>
using namespace matplot;
#endif

// Test utilities
namespace TestUtilities {

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

  // Plot a circle
  template <typename Real, Integer N>
  void plot_circle(AABBtree::Vector<Real, 2> const & center, Real const radius, std::string const & color,
    Real const line_width = 1, Integer const n_points = 100) {
    std::vector<Real> x, y;
    Real const d_theta{2.0*M_PI/n_points};
    for (Integer i{0}; i < n_points; ++i) {
      Real const theta{d_theta*i};
      x.push_back(center[0] + radius*std::cos(theta));
      y.push_back(center[1] + radius*std::sin(theta));
    }
    if (ax) {
      ax->plot(x, y)->line_width(line_width).color(color);
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

} // namespace TestUtilities

#endif // INCLUDE_TEST_UTILITIES_HH
