/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The AABBtree project is distributed under the BSD 2-Clause License.       *
 *                                                                           *
 * Davide Stocco Enrico Bertolazzi                                           *
 * University of Trento University of Trento                                 *
 * davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// DISCLAIMER: The code in this file is a modified version of the Eigen library.

#pragma once

#ifndef INCLUDE_AABBTREE_RAY_HH
#define INCLUDE_AABBTREE_RAY_HH

#include "AABBtree/Box.hh"

namespace AABBtree {

/**
 * \brief A mathematical ray in N-dimensional space.
 *
 * This class represents a ray defined by an origin point and a direction
 * vector. The ray is the set of points \f$ \mathbf{o} + t\,\mathbf{d} \f$ with
 * \f$ t \geq 0 \f$. The class provides functionality for:
 * - Ray transformations (translation, rotation, scaling);
 * - Distance computations to points and boxes;
 * - Intersection tests with axis-aligned bounding boxes;
 * - Normalization and comparison operations.
 *
 * \tparam Real Floating-point scalar type (must be float or double).
 * \tparam N Dimension of the ambient space.
 * \note The direction vector is not required to be normalized.
 */
template <std::floating_point Real, Integer N> class Ray {
  static_assert(N > 0, "Ray dimension must be positive.");

  using Point = AABBtree::Point<Real, N>;
  using Vector = AABBtree::Vector<Real, N>;
  using Box = AABBtree::Box<Real, N>;

  Point m_origin;     /**< Origin point of the ray. */
  Vector m_direction; /**< Direction vector of the ray. */

public:
  /**
   * Default class destructor.
   */
  ~Ray() = default;

  /**
   * Default class constructor. The origin and the direction are not
   * initialized.
   */
  Ray() = default;

  /**
   * Copy constructor
   * \param[in] r Ray to copy
   */
  Ray(Ray const &r) = default;

  /**
   * Copy assignment operator.
   * \param[in] r Ray to copy
   * \return A reference to the current ray.
   */
  Ray &operator=(Ray const &r) = default;

  /**
   * Class constructor for a ray given an origin and a direction.
   * \param[in] o Origin of the ray.
   * \param[in] d Direction of the ray.
   */
  Ray(Point const &o, Vector const &d) : m_origin(o), m_direction(d) {}

  /**
   * Class constructor for the 1D ray.
   * \param[in] o Origin of the ray.
   * \param[in] d Direction of the ray.
   * \note This constructor is only available for 1D rays.
   */
  Ray(Real const o, Real const d)
    requires(N == 1) : m_origin(o), m_direction(d) {}

  /**
   * Class constructor for the 2D ray.
   * \param[in] o_x Origin \f$ x \f$-axis component.
   * \param[in] o_y Origin \f$ y \f$-axis component.
   * \param[in] d_x Direction \f$ x \f$-axis component.
   * \param[in] d_y Direction \f$ y \f$-axis component.
   * \note This constructor is only available for 2D rays.
   */
  Ray(Real const o_x, Real const o_y, Real const d_x, Real const d_y)
    requires(N == 2)
      : m_origin(o_x, o_y), m_direction(d_x, d_y) {}

  /**
   * Class constructor for the 3D ray.
   * \param[in] o_x Origin \f$ x \f$-axis component.
   * \param[in] o_y Origin \f$ y \f$-axis component.
   * \param[in] o_z Origin \f$ z \f$-axis component.
   * \param[in] d_x Direction \f$ x \f$-axis component.
   * \param[in] d_y Direction \f$ y \f$-axis component.
   * \param[in] d_z Direction \f$ z \f$-axis component.
   * \note This constructor is only available for 3D rays.
   */
  Ray(Real const o_x, Real const o_y, Real const o_z, Real const d_x,
      Real const d_y, Real const d_z)
    requires(N == 3)
      : m_origin(o_x, o_y, o_z), m_direction(d_x, d_y, d_z) {}

  /**
   * Copy constructor for a ray given another ray with a different scalar type.
   * \param[in] r Ray to copy.
   * \tparam OtherReal Type of the scalar coefficients of the ray to copy.
   */
  template <std::floating_point OtherReal>
  explicit Ray(Ray<OtherReal, N> const &r)
      : m_origin(r.origin().template cast<Real>()),
        m_direction(r.direction().template cast<Real>()) {}

  /**
   * Cast the current object to a new scalar type.
   * \tparam NewReal the new scalar type
   * \return A copy of the ray with the new scalar type.
   */
  template <std::floating_point NewReal> Ray<NewReal, N> cast() const {
    return Ray<NewReal, N>(m_origin.template cast<NewReal>(),
                           m_direction.template cast<NewReal>());
  }

  /**
   * Get the reference to the ray origin.
   * \return The reference to the ray origin.
   */
  Point &origin() { return m_origin; }

  /**
   * Get the const reference to the ray origin.
   * \return The const reference to the ray origin.
   */
  Point const &origin() const { return m_origin; }

  /**
   * Get the reference to the ray direction.
   * \return The reference to the ray direction.
   */
  Vector &direction() { return m_direction; }

  /**
   * Get the const reference to the ray direction.
   * \return The const reference to the ray direction.
   */
  Vector const &direction() const { return m_direction; }

  /**
   * Get the point of the ray at a given parameter.
   * \param[in] t Ray parameter (non-negative for points of the ray).
   * \return The point \f$ \mathbf{o} + t\,\mathbf{d} \f$.
   */
  Point point(Real const t) const { return m_origin + t * m_direction; }

  /**
   * Normalize the direction of the ray.
   * \return A reference to the current ray.
   */
  Ray &normalize() {
    m_direction.normalize();
    return *this;
  }

  /**
   * Normalize the direction of the ray.
   * \return A copy of the current normalized ray.
   */
  Ray normalized() const { return Ray(m_origin, m_direction.normalized()); }

  /**
   * Check if the current ray is approximately equal to another ray.
   * \param[in] r Ray to compare with.
   * \param[in] tol Tolerance to use for the comparison.
   */
  bool
  is_approx(Ray const &r,
            Real const tol = Eigen::NumTraits<Real>::dummy_precision()) const {
    return m_origin.isApprox(r.m_origin, tol) &&
           m_direction.isApprox(r.m_direction, tol);
  }

  /**
   * Translate the current ray by a given vector.
   * \param[in] t Vector to translate the ray by.
   * \return A reference to the current ray.
   */
  Ray &translate(Vector const &t) {
    m_origin += t;
    return *this;
  }

  /**
   * Translate the current ray by a given vector.
   * \param[in] t Vector to translate the ray by.
   * \return A copy of the current translated ray.
   */
  Ray translated(Vector const &t) const {
    return Ray(m_origin + t, m_direction);
  }

  /**
   * Transform the current ray.
   * \param[in] t Affine transformation to apply (e.g. \c Eigen::Transform).
   * \return A copy of the current transformed ray.
   * \tparam Transform Type of the transformation: it must support
   * \c t * point and \c t.linear().
   */
  template <typename Transform> Ray transformed(Transform const &t) const {
    return Ray(t * m_origin, t.linear() * m_direction);
  }

  /**
   * Transform the current ray.
   * \param[in] t Affine transformation to apply (e.g. \c Eigen::Transform).
   * \return A reference to the current ray.
   * \tparam Transform Type of the transformation (see \c transformed).
   */
  template <typename Transform> Ray &transform(Transform const &t) {
    *this = this->transformed(t);
    return *this;
  }

  /**
   * Check if the point is on the ray.
   * \param[in] p Point to check.
   * \param[in] tol Angular tolerance (the point may be off the ray by at most
   * \c tol times its distance from the origin).
   * \return True if the point is on the ray, false otherwise.
   */
  bool contains(Point const &p,
                Real tol = Eigen::NumTraits<Real>::dummy_precision()) const {
    Vector const v{p - m_origin};
    Real const v_norm{v.norm()};
    if (v_norm <= tol) {
      return true; // The origin belongs to the ray
    }
    Vector const u{v / v_norm};
    Vector const d{m_direction.normalized()};
    Real const along{u.dot(d)};
    // Component of u orthogonal to d (the sine of the angle between them)
    return along >= -tol && (u - along * d).norm() <= tol;
  }

  /**
   * Check if the current ray intersects a given axis-aligned box.
   * \param[in] b Box to check.
   * \param[in] tol Tolerance to use for the intersection.
   * \return True if the current ray intersects the given box, false otherwise.
   */
  bool intersects(Box const &b,
                  Real tol = Eigen::NumTraits<Real>::dummy_precision()) const {
    Real t_entry, t_exit;
    return this->clip(b, t_entry, t_exit, tol);
  }

  /**
   * Intersect the current ray with a given axis-aligned box.
   * \param[in] b Box to intersect with.
   * \param[out] c Closest intersection point (the origin if it is inside).
   * \param[out] f Farthest intersection point.
   * \param[in] tol Tolerance to use for the intersection.
   * \return True if the current ray intersects the given box, false otherwise
   * (in that case \c c and \c f are not modified).
   */
  bool intersect(Box const &b, Point &c, Point &f,
                 Real tol = Eigen::NumTraits<Real>::dummy_precision()) const {
    Real t_entry, t_exit;
    if (!this->clip(b, t_entry, t_exit, tol)) {
      return false;
    }
    c = this->point(t_entry);
    f = this->point(t_exit);
    return true;
  }

  /**
   * Compute the squared distance between the current ray a given point.
   * \param[in] p Point to compute the squared distance to.
   * \param[in] tol Unused, kept for backward compatibility.
   * \return The squared distance between the ray and the point.
   */
  Real
  squared_distance(Point const &p,
                   Real tol = Eigen::NumTraits<Real>::dummy_precision()) const {
    Point c;
    return this->squared_distance(p, c, tol);
  }

  /**
   * Compute the squared distance between the current ray a given point,
   * returning the closest point on the ray.
   * \param[in] p Point to compute the distance to.
   * \param[out] c Closest point on the ray.
   * \param[in] tol Unused, kept for backward compatibility.
   * \return The squared distance between the ray and the point.
   */
  Real
  squared_distance(Point const &p, Point &c,
                   Real tol = Eigen::NumTraits<Real>::dummy_precision()) const {
    static_cast<void>(tol);
    Real const d2{m_direction.squaredNorm()};
    Real const t{d2 > 0 ? (p - m_origin).dot(m_direction) / d2 : Real(0)};
    c = this->point(std::max(Real(0), t)); // Points behind map to the origin
    return (c - p).squaredNorm();
  }

  /**
   * Compute the distance between the current ray a given point.
   * \param[in] p Point to compute the distance to.
   * \param[in] tol Unused, kept for backward compatibility.
   * \return The distance between the ray and the point.
   */
  Real distance(Point const &p,
                Real tol = Eigen::NumTraits<Real>::dummy_precision()) const {
    return std::sqrt(this->squared_distance(p, tol));
  }

  /**
   * Compute the distance between the current ray a given point, returning the
   * closest point on the ray.
   * \param[in] p Point to compute the distance to.
   * \param[out] c Closest point on the ray.
   * \param[in] tol Unused, kept for backward compatibility.
   * \return The distance between the ray and the point.
   */
  Real distance(Point const &p, Point &c,
                Real tol = Eigen::NumTraits<Real>::dummy_precision()) const {
    return std::sqrt(this->squared_distance(p, c, tol));
  }

  /**
   * Compute the squared \em interior (or \em minimum) distance between the
   * current ray and a given box, i.e. the minimum distance between a point of
   * the ray and a point of the box.
   * \param[in] b Box to compute the squared distance to.
   * \param[in] tol Unused, kept for backward compatibility.
   * \return The squared distance between the ray and the box (zero if they
   * intersect).
   */
  Real squared_interior_distance(
      Box const &b, Real tol = Eigen::NumTraits<Real>::dummy_precision()) const {
    Point p1, p2;
    return this->squared_interior_distance(b, p1, p2, tol);
  }

  /**
   * Compute the squared \em interior (or \em minimum) distance between the
   * current ray and a given box, returning two points at the minimum distance.
   * \param[in] b Box to compute the squared distance to.
   * \param[out] p1 Point at the minimum distance on the current ray.
   * \param[out] p2 Point at the minimum distance on the box.
   * \param[in] tol Unused, kept for backward compatibility.
   * \return The squared distance between the ray and the box (zero if they
   * intersect).
   */
  Real squared_interior_distance(
      Box const &b, Point &p1, Point &p2,
      Real tol = Eigen::NumTraits<Real>::dummy_precision()) const {
    static_cast<void>(tol);
    Real const t{this->minimize_over_ray(b, DistanceKind::INTERIOR)};
    p1 = this->point(t);
    return b.squared_interior_distance(p1, p2);
  }

  /**
   * Compute the \em interior (or \em minimum) distance between the current ray
   * and a given box.
   * \param[in] b Box to compute the distance to.
   * \param[in] tol Unused, kept for backward compatibility.
   * \return The distance between the ray and the box (zero if they intersect).
   */
  Real interior_distance(
      Box const &b, Real tol = Eigen::NumTraits<Real>::dummy_precision()) const {
    return std::sqrt(this->squared_interior_distance(b, tol));
  }

  /**
   * Compute the \em interior (or \em minimum) distance between the current ray
   * and a given box, returning two points at the minimum distance.
   * \param[in] b Box to compute the distance to.
   * \param[out] p1 Point at the minimum distance on the current ray.
   * \param[out] p2 Point at the minimum distance on the box.
   * \param[in] tol Unused, kept for backward compatibility.
   * \return The distance between the ray and the box (zero if they intersect).
   */
  Real interior_distance(
      Box const &b, Point &p1, Point &p2,
      Real tol = Eigen::NumTraits<Real>::dummy_precision()) const {
    return std::sqrt(this->squared_interior_distance(b, p1, p2, tol));
  }

  /**
   * Compute the squared \em exterior distance between the current ray and a
   * given box, i.e. the minimum over the points of the ray of the distance to
   * the farthest point of the box.
   * \param[in] b Box to compute the squared distance to.
   * \param[in] tol Unused, kept for backward compatibility.
   * \return The squared exterior distance between the ray and the box.
   * \note This is an upper bound of the distance between the ray and any point
   * of the box, as the interior distance is a lower bound.
   */
  Real squared_exterior_distance(
      Box const &b, Real tol = Eigen::NumTraits<Real>::dummy_precision()) const {
    Point p1, p2;
    return this->squared_exterior_distance(b, p1, p2, tol);
  }

  /**
   * Compute the squared \em exterior distance between the current ray and a
   * given box (see above), returning the two points at that distance.
   * \param[in] b Box to compute the squared distance to.
   * \param[out] p1 Point of the ray closest to its farthest box point.
   * \param[out] p2 Farthest point of the box from \c p1.
   * \param[in] tol Unused, kept for backward compatibility.
   * \return The squared exterior distance between the ray and the box.
   */
  Real squared_exterior_distance(
      Box const &b, Point &p1, Point &p2,
      Real tol = Eigen::NumTraits<Real>::dummy_precision()) const {
    static_cast<void>(tol);
    Real const t{this->minimize_over_ray(b, DistanceKind::EXTERIOR)};
    p1 = this->point(t);
    return b.squared_exterior_distance(p1, p2);
  }

  /**
   * Compute the \em exterior distance between the current ray and a given box
   * (see \c squared_exterior_distance).
   * \param[in] b Box to compute the distance to.
   * \param[in] tol Unused, kept for backward compatibility.
   * \return The exterior distance between the ray and the box.
   */
  Real exterior_distance(
      Box const &b, Real tol = Eigen::NumTraits<Real>::dummy_precision()) const {
    return std::sqrt(this->squared_exterior_distance(b, tol));
  }

  /**
   * Compute the \em exterior distance between the current ray and a given box
   * (see \c squared_exterior_distance), returning two points at that distance.
   * \param[in] b Box to compute the distance to.
   * \param[out] p1 Point of the ray closest to its farthest box point.
   * \param[out] p2 Farthest point of the box from \c p1.
   * \param[in] tol Unused, kept for backward compatibility.
   * \return The exterior distance between the ray and the box.
   */
  Real exterior_distance(
      Box const &b, Point &p1, Point &p2,
      Real tol = Eigen::NumTraits<Real>::dummy_precision()) const {
    return std::sqrt(this->squared_exterior_distance(b, p1, p2, tol));
  }

  /**
   * Print the ray info to an output stream.
   * \param[in] os Output stream to print the ray info to.
   */
  void print(std::ostream &os) const {
    os << "RAY INFO" << std::endl
       << "\to = " << m_origin.transpose() << std::endl
       << "\td = " << m_direction.transpose() << std::endl;
  }

private:
  /**
   * Clip the ray against a box (slab method).
   * \param[in] b Box to clip against.
   * \param[out] t_entry Parameter where the ray enters the box (0 if the
   * origin is inside).
   * \param[out] t_exit Parameter where the ray exits the box.
   * \param[in] tol Directions components smaller than \c tol are considered
   * zero (the ray is parallel to the corresponding slab).
   * \return True if the ray intersects the box.
   */
  bool clip(Box const &b, Real &t_entry, Real &t_exit, Real const tol) const {
    t_entry = 0; // Only the points with t >= 0 belong to the ray
    t_exit = std::numeric_limits<Real>::infinity();
    for (Integer i{0}; i < N; ++i) {
      Real const o{m_origin[i]};
      Real const d{m_direction[i]};
      if (std::abs(d) <= tol) {
        // Parallel to the slab: either always inside or never
        if (o < b.min()[i] || o > b.max()[i]) {
          return false;
        }
        continue;
      }
      Real t_near{(b.min()[i] - o) / d};
      Real t_far{(b.max()[i] - o) / d};
      if (t_near > t_far) {
        std::swap(t_near, t_far);
      }
      t_entry = std::max(t_entry, t_near);
      t_exit = std::min(t_exit, t_far);
    }
    return t_entry <= t_exit + tol;
  }

  /**
   * Kind of point-box distance minimized along the ray.
   */
  enum class DistanceKind {
    INTERIOR, /**< Distance to the closest point of the box. */
    EXTERIOR  /**< Distance to the farthest point of the box. */
  };

  /**
   * Find the ray parameter \f$ t \geq 0 \f$ that minimizes the squared
   * distance between the ray point \f$ \mathbf{o} + t\,\mathbf{d} \f$ and the
   * box (closest or farthest box point, depending on \c kind).
   *
   * Along each axis the squared distance is a convex, piecewise quadratic
   * function of the coordinate:
   * - interior: \f$ \max(0,\, m - x,\, x - M)^2 \f$, with breakpoints at the
   *   box faces \f$ m \f$ and \f$ M \f$;
   * - exterior: \f$ (|x - c| + h)^2 \f$, with a breakpoint at the box center
   *   \f$ c \f$ (\f$ h \f$ is the half size).
   *
   * Substituting \f$ x = o + t\,d \f$, the sum over the axes is a convex,
   * piecewise quadratic function of \f$ t \f$. Between two consecutive
   * breakpoints it is a single quadratic \f$ A t^2 + B t + C \f$, whose
   * minimum on the interval is found in closed form. The smallest of these
   * minima is the exact global minimum.
   */
  Real minimize_over_ray(Box const &b, DistanceKind const kind) const {
    Point const &b_min{b.min()};
    Point const &b_max{b.max()};
    Real const infinity{std::numeric_limits<Real>::infinity()};

    // Collect the breakpoints (where the ray crosses a face or a center plane)
    std::array<Real, 2 * N + 2> breaks;
    std::size_t n_breaks{0};
    breaks[n_breaks++] = 0;
    for (Integer i{0}; i < N; ++i) {
      Real const d{m_direction[i]};
      if (d == 0) {
        continue;
      }
      if (kind == DistanceKind::INTERIOR) {
        Real const t_1{(b_min[i] - m_origin[i]) / d};
        Real const t_2{(b_max[i] - m_origin[i]) / d};
        if (t_1 > 0) breaks[n_breaks++] = t_1;
        if (t_2 > 0) breaks[n_breaks++] = t_2;
      } else {
        Real const t_c{(Real(0.5) * (b_min[i] + b_max[i]) - m_origin[i]) / d};
        if (t_c > 0) breaks[n_breaks++] = t_c;
      }
    }
    // Insertion sort: there are at most 2N + 1 values (std::sort on such small
    // arrays triggers false -Warray-bounds warnings with GCC)
    for (std::size_t k{1}; k < n_breaks; ++k) {
      for (std::size_t h{k}; h > 0 && breaks[h] < breaks[h - 1]; --h) {
        std::swap(breaks[h], breaks[h - 1]);
      }
    }
    breaks[n_breaks++] = infinity;

    // Minimize the quadratic on each interval [breaks[k], breaks[k+1]]
    Real best_t{0};
    Real best_value{infinity};
    for (std::size_t k{0}; k + 1 < n_breaks; ++k) {
      Real const t_lo{breaks[k]};
      Real const t_hi{breaks[k + 1]};
      if (t_hi <= t_lo) {
        continue; // Repeated breakpoint
      }
      // A point strictly inside the interval selects the quadratic piece
      Real const t_mid{t_hi < infinity ? Real(0.5) * (t_lo + t_hi) : t_lo + 1};

      // Accumulate A t^2 + B t + C, each axis contributes (d t + e)^2
      Real A{0}, B{0}, C{0};
      for (Integer i{0}; i < N; ++i) {
        Real const o{m_origin[i]};
        Real const d{m_direction[i]};
        Real const x{o + t_mid * d};
        Real slope{0}, offset{0};
        if (kind == DistanceKind::INTERIOR) {
          if (x < b_min[i]) { // (m - o - d t)^2
            slope = -d;
            offset = b_min[i] - o;
          } else if (x > b_max[i]) { // (o + d t - M)^2
            slope = d;
            offset = o - b_max[i];
          }
        } else {
          Real const c{Real(0.5) * (b_min[i] + b_max[i])};
          Real const h{Real(0.5) * (b_max[i] - b_min[i])};
          Real const s{x >= c ? Real(1) : Real(-1)}; // (s (o + d t - c) + h)^2
          slope = s * d;
          offset = s * (o - c) + h;
        }
        A += slope * slope;
        B += 2 * slope * offset;
        C += offset * offset;
      }

      // Minimum of the quadratic on the interval
      Real t{t_lo};
      if (A > 0) {
        t = std::clamp(-B / (2 * A), t_lo, t_hi);
      }
      Real const value{(A * t + B) * t + C};
      if (value < best_value) {
        best_value = value;
        best_t = t;
      }
    }
    return best_t;
  }

}; // class Ray

/**
 * Print the ray info to an output stream.
 * \param[in] os Output stream to print the ray info to.
 * \param[in] r Ray to print.
 * \tparam Real Type of the scalar coefficients.
 * \tparam N Dimension of the ambient space.
 */
template <std::floating_point Real, Integer N>
std::ostream &operator<<(std::ostream &os, Ray<Real, N> const &r) {
  r.print(os);
  return os;
}

} // namespace AABBtree

#endif // INCLUDE_AABBTREE_RAY_HH
