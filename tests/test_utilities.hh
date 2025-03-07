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
  Box<Real, 2> box() const {
    return Box<Real, 2>(
      std::min(this->m_p_1[0], this->m_p_2[0]),
      std::min(this->m_p_1[1], this->m_p_2[1]),
      std::max(this->m_p_1[0], this->m_p_2[0]),
      std::max(this->m_p_1[1], this->m_p_2[1])
    );;
  }
  bool intersect (Segment const & segment, Vector2 & point) const {
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

// Plot a box
template <typename Real, Integer N>
void plot_box(Box<Real, N> const & box, std::string const & color, Real const line_width = 1) {
  static_assert(N == 2, "Plotting is only supported for 2D objects.");
  std::vector<Real> min(box.min().data(), box.min().data() + N);
  std::vector<Real> max(box.max().data(), box.max().data() + N);
  #ifdef AABBTREE_ENABLE_PLOTTING
  if (ax) {
    plot(ax, {min[0], max[0], max[0], min[0], min[0]},
      {min[1], min[1], max[1], max[1], min[1]})->line_width(line_width).color(color);
  }
  #endif
}

// Plot a tree
template <typename Real, Integer N>
void plot_tree(NonRecursive<Real, N> const & tree, std::string const & color, Real const line_width = 1) {
  static_assert(N == 2, "Plotting is only supported for 2D objects.");
  #ifdef AABBTREE_ENABLE_PLOTTING
  for (Integer i{0}; i < static_cast<Integer>(tree.size()); ++i) {
    Box<Real, N> const & box{tree.node(i).box};
    std::cout << "Box " << i << ": " << box.min().transpose() << " " << box.max().transpose() << std::endl;
    plot_box<Real, N>(box, color, line_width);
  }
  #endif
}

// Plot a point
template <typename Real, Integer N>
void plot_point(AABBtree::Vector<Real, N> const & point, std::string const & color,
  Real const marker_size = 0.5) {
  static_assert(N == 2, "Plotting is only supported for 2D objects.");
  std::vector<Real> p(point.data(), point.data() + N);
  #ifdef AABBTREE_ENABLE_PLOTTING
  if (ax) {
    plot(ax, {p[0]}, {p[1]}, "o")->marker_size(marker_size).color(color).marker_face(true);
  }
  #endif
}

// Plot a segment
template <typename Real, Integer N>
void plot_segment(AABBtree::Vector<Real, N> const & p_1, AABBtree::Vector<Real, N> const & p_2,
  std::string const & color, Real const line_width = 1) {
  static_assert(N == 2, "Plotting is only supported for 2D objects.");
  std::vector<Real> v_1(p_1.data(), p_1.data() + N);
  std::vector<Real> v_2(p_2.data(), p_2.data() + N);
  #ifdef AABBTREE_ENABLE_PLOTTING
  if (ax) {
    plot(ax, {v_1[0], v_2[0]}, {v_1[1], v_2[1]}, "-o")->color(color).line_width(line_width).
      marker_size(2.0*line_width).marker_face(true);
  }
  #endif
}

// Plot a segment
template <typename Real>
void plot_segment(Segment<Real> const & segment, std::string const & color, Real const line_width = 1) {
  plot_segment<Real, 2>(segment.p_1(), segment.p_2(), color, line_width);
}

#endif // INCLUDE_TEST_UTILITIES_HH
