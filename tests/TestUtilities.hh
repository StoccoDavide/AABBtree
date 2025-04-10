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
#include "Plot2D.hh"
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
      );
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
  // Plot a segment
  template <typename Real>
  void plot_segment(plot_obj & ax, Segment<Real> const & segment, std::string const & color,
    Real const line_width = 1.0)
  {
    plot_segment<Real, 2>(ax, segment.p_1(), segment.p_2(), color, line_width);
  }
  #endif // AABBTREE_ENABLE_PLOTTING

} // namespace TestUtilities

#endif // INCLUDE_TEST_UTILITIES_HH
