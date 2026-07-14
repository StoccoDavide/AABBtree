/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The AABBtree project is distributed under the BSD 2-Clause License.       *
 *                                                                           *
 * Davide Stocco Enrico Bertolazzi                                           *
 * University of Trento University of Trento                                 *
 * davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef INCLUDE_TEST_UTILITIES_HH
#define INCLUDE_TEST_UTILITIES_HH

// C++17 standard libraries
#include <random>
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
template <typename Real> class Segment {
  using Vector2 = AABBtree::Vector<Real, 2>;
  Vector2 m_p_1;
  Vector2 m_p_2;

public:
  ~Segment() = default;
  Segment() = default;
  Segment(Segment const &) = default;
  Segment(Real t_x_1, Real t_y_1, Real t_x_2, Real t_y_2)
      : m_p_1(t_x_1, t_y_1), m_p_2(t_x_2, t_y_2) {}
  Segment(Vector2 const &p_1, Vector2 const &p_2) : m_p_1(p_1), m_p_2(p_2) {}
  Vector2 const &p_1() const { return this->m_p_1; }
  Vector2 const &p_2() const { return this->m_p_2; }
  Vector2 &p_1() { return this->m_p_1; }
  Vector2 &p_2() { return this->m_p_2; }
  Vector2 const &point(Integer i) const {
    if (i == 0) {
      return this->m_p_1;
    }
    if (i == 1) {
      return this->m_p_2;
    }
    throw std::out_of_range("Segment point index out of range.");
  }
  Vector2 &point(Integer i) {
    if (i == 0) {
      return this->m_p_1;
    }
    if (i == 1) {
      return this->m_p_2;
    }
    throw std::out_of_range("Segment point index out of range.");
  }
  Box<Real, 2> bounding_box() const {
    return Box<Real, 2>(std::min(this->m_p_1[0], this->m_p_2[0]),
                        std::min(this->m_p_1[1], this->m_p_2[1]),
                        std::max(this->m_p_1[0], this->m_p_2[0]),
                        std::max(this->m_p_1[1], this->m_p_2[1]));
  }
  bool intersect(Segment const &segment, Vector2 &point) const {
    Real const &t_x_1{this->m_p_1[0]};
    Real const &t_y_1{this->m_p_1[1]};
    Real const &t_x_2{this->m_p_2[0]};
    Real const &t_y_2{this->m_p_2[1]};
    Real const &s_x_1{segment.m_p_1[0]};
    Real const &s_y_1{segment.m_p_1[1]};
    Real const &s_x_2{segment.m_p_2[0]};
    Real const &s_y_2{segment.m_p_2[1]};
    Real const d{(t_x_1 - t_x_2) * (s_y_1 - s_y_2) -
                 (t_y_1 - t_y_2) * (s_x_1 - s_x_2)};
    if (std::abs(d) < std::numeric_limits<Real>::epsilon())
      return false;
    Real const t{((t_x_1 - s_x_1) * (s_y_1 - s_y_2) -
                  (t_y_1 - s_y_1) * (s_x_1 - s_x_2)) /
                 d};
    Real const u{-((t_x_1 - t_x_2) * (t_y_1 - s_y_1) -
                   (t_y_1 - t_y_2) * (t_x_1 - s_x_1)) /
                 d};
    if (t >= 0.0 && t <= 1.0 && u >= 0.0 && u <= 1.0) {
      point.x() = t_x_1 + t * (t_x_2 - t_x_1);
      point.y() = t_y_1 + t * (t_y_2 - t_y_1);
      return true;
    }
    return false;
  }
};

// Random boxes generation utilities (shared by the complexity_randbox_*
// benchmarks)
//
// The parameters below are shared by every complexity_randbox_* benchmark so
// that all of them are tested against the very same (reproducible) random data.

constexpr int RANDBOX_DIMENSION{100};       // Ambient space dimension
constexpr unsigned RANDBOX_SEED_1{123456u}; // Seed for the first set of boxes
constexpr unsigned RANDBOX_SEED_2{123456u}; // Seed for the second set of boxes
constexpr double RANDBOX_DOMAIN_MIN{
    0.0}; // Lower bound of the domain used to place the boxes
constexpr double RANDBOX_DOMAIN_MAX{
    100.0}; // Upper bound of the domain used to place the boxes
constexpr double RANDBOX_SIZE_MIN{0.0}; // Minimum edge length of a random box
constexpr double RANDBOX_SIZE_MAX{1.0}; // Maximum edge length of a random box

// Used as a search radius by point-based data structures (CGAL kd-tree,
// nanoflann) that cannot be queried with an actual box, to approximate a
// box-box intersection query with a fuzzy-sphere/radius query.
constexpr double RANDBOX_QUERY_RADIUS{RANDBOX_SIZE_MAX};

// Axis-aligned box defined through generic (dimension-agnostic) corners
template <typename Real, int N> struct RandomBox {
  std::array<Real, N> min; /**< Minimal corner of the box. */
  std::array<Real, N> max; /**< Maximal corner of the box. */
};

// Generate a reproducible set of random axis-aligned boxes in N dimensions
template <typename Real, int N, typename Integer>
std::vector<RandomBox<Real, N>>
GenerateRandomBoxes(Integer n_boxes, unsigned seed,
                    Real domain_min = static_cast<Real>(RANDBOX_DOMAIN_MIN),
                    Real domain_max = static_cast<Real>(RANDBOX_DOMAIN_MAX),
                    Real size_min = static_cast<Real>(RANDBOX_SIZE_MIN),
                    Real size_max = static_cast<Real>(RANDBOX_SIZE_MAX)) {
  std::vector<RandomBox<Real, N>> boxes;
  boxes.reserve(n_boxes);
  std::mt19937_64 generator(seed);
  std::uniform_real_distribution<Real> center_distribution(domain_min,
                                                           domain_max);
  std::uniform_real_distribution<Real> size_distribution(size_min, size_max);
  for (Integer i{0}; i < n_boxes; ++i) {
    RandomBox<Real, N> box;
    for (int d{0}; d < N; ++d) {
      Real center{center_distribution(generator)};
      Real half_size{static_cast<Real>(0.5) * size_distribution(generator)};
      box.min[d] = center - half_size;
      box.max[d] = center + half_size;
    }
    boxes.push_back(box);
  }
  return boxes;
}

#ifdef AABBTREE_ENABLE_PLOTTING
// Plot a segment
template <typename Real>
void plot_segment(plot_obj &ax, Segment<Real> const &segment,
                  std::string const &color, Real const line_width = 1.0) {
  plot_segment<Real, 2>(ax, segment.p_1(), segment.p_2(), color, line_width);
}
#endif // AABBTREE_ENABLE_PLOTTING

} // namespace TestUtilities

#endif // INCLUDE_TEST_UTILITIES_HH
