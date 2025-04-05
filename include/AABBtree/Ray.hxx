/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                     *
 *                                                                                              *
 * The AABBtree project is distributed under the BSD 2-Clause License.                          *
 *                                                                                              *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// DISCLAIMER: The code in this file is a modified version of the Eigen library.

#pragma once

#ifndef AABBTREE_RAY_HXX
#define AABBTREE_RAY_HXX

#include "AABBtree/Box.hxx"

namespace AABBtree {

  /**
   * \class Ray
   * \brief A mathematical ray in N-dimensional space
   *
   * This class represents a ray defined by an origin point and a direction vector.
   * The ray extends infinitely in the direction of the vector from the origin point.
   *
   * The class provides functionality for:
   * - Ray transformations (translation, rotation, scaling)
   * - Distance computations to points and boxes
   * - Intersection tests with axis-aligned bounding boxes
   * - Normalization and comparison operations
   *
   * \tparam Real Floating-point scalar type (must be float or double)
   * \tparam N Dimension of the space (1, 2, or 3)
   *
   * \note The direction vector is not required to be normalized, but some operations
   *       (like distance calculations) may internally normalize it.
   */
  template <typename Real, Integer N>
  class Ray {
    static_assert(std::is_floating_point<Real>::value, "Ray Real type must be a floating-point type.");
    static_assert(N > 0, "Ray dimension must be positive.");

    // Numerical constants
    constexpr static Real EPS{std::numeric_limits<Real>::epsilon()}; ///< Machine epsilon for floating-point operations
    constexpr static Real MAX{static_cast<Real>(1.0)/EPS};           ///< Large value representing infinity
    constexpr static Real MIN{-MAX};                                 ///< Negative infinity
    constexpr static Real DUMMY_TOL{EPS*static_cast<Real>(100.0)};   ///< Default tolerance for comparisons

    using Point  = Point<Real, N>;
    using Vector = Vector<Real, N>;

    Point  m_origin;    ///< Origin point of the ray
    Vector m_direction; ///< Direction vector of the ray (not necessarily normalized)

  public:
    /**
     * \brief Default destructor
     */
    ~Ray() = default;

    /**
     * \brief Default constructor
     * 
     * Creates a ray with origin at (0,0,...) and direction (1,0,...)
     */
    Ray() = default;

    /**
     * \brief Copy constructor
     * \param[in] r Ray to copy
     */
    Ray(Ray const & r) : m_origin(r.m_origin), m_direction(r.m_direction) {}

    /**
     * \brief Class constructor for a ray given an origin and a direction.
     * \param[in] o Origin of the ray.
     * \param[in] d Direction of the ray.
     */
    Ray(Point const & o, Vector const & d) : m_origin(o), m_direction(d) {}

    /**
     * \brief Class constructor for the 1D ray.
     * \param[in] o Origin of the ray.
     * \param[in] d Direction of the ray.
     * \tparam T Type of the scalar coefficients.
     * \note This constructor is only available for 1D rays.
     */
    template <typename = std::enable_if<N == 1>>
    Ray(Real const o, Real const d) : m_origin(o), m_direction(d) {}

    /**
     * \brief Class constructor for the 2D ray.
     * \param[in] o_x Origin \f$ x \f$-axis component.
     * \param[in] o_y Origin \f$ y \f$-axis component.
     * \param[in] d_x Direction \f$ x \f$-axis component.
     * \param[in] d_y Direction \f$ y \f$-axis component.
     * \tparam T Type of the scalar coefficients.
     * \note This constructor is only available for 2D rays.
     */
    template <typename = std::enable_if<N == 2>>
    Ray(Real const o_x, Real const o_y, Real const d_x, Real const d_y)
    : m_origin(o_x, o_y), m_direction(d_x, d_y) {}

    /**
     * \brief Class constructor for the 3D ray.
     * \param[in] o_x Origin \f$ x \f$-axis component.
     * \param[in] o_y Origin \f$ y \f$-axis component.
     * \param[in] o_z Origin \f$ z \f$-axis component.
     * \param[in] d_x Direction \f$ x \f$-axis component.
     * \param[in] d_y Direction \f$ y \f$-axis component.
     * \param[in] d_z Direction \f$ z \f$-axis component.
     * \tparam T Type of the scalar coefficients.
     * \note This constructor is only available for 3D rays.
     */
    template <typename = std::enable_if<N == 3>>
    Ray(Real const o_x, Real const o_y, Real const o_z, Real const d_x, Real const d_y, Real const d_z)
    : m_origin(o_x, o_y, o_z), m_direction(d_x, d_y, d_z) {}

    /**
     * Copy constructor for a ray given another ray with a different scalar type.
     * \param[in] r Ray to copy.
     * \tparam OtherReal Type of the scalar coefficients of the ray to copy.
     */
    template<typename OtherReal>
    explicit Ray(Ray<OtherReal, N> const & r)
    : m_origin(r.m_origin.template cast<Real>()), m_direction(r.m_direction.template cast<Real>()) {}

    /**
     * \brief Cast the current object to a new scalar type.
     * \tparam NewReal the new scalar type
     * \note If the new real type is equal to the current scalar type currently used, then this
     * function returns a const reference to the current object.
     */
    template<typename NewReal>
    Ray<NewReal, N> cast() const {
      if constexpr (std::is_same<Real, NewReal>::value) return *this;
      return Ray<NewReal, N>(m_origin.template cast<NewReal>(), m_direction.template cast<NewReal>());
    }

    /**
     * \brief Get the reference to the ray origin.
     * \return The reference to the ray origin.
     */
    Point & origin() { return m_origin; }

    /**
     * \brief Get the const reference to the ray origin.
     * \return The const reference to the ray origin.
     */
    Point const & origin() const { return m_origin; }

    /**
     * \brief Get the reference to the ray direction.
     * \return The reference to the ray direction.
     */
    Vector & direction() { return m_direction; }

    /**
     * \brief Get the const reference to the ray direction.
     * \return The const reference to the ray direction.
     */
    Vector const & direction() const { return m_direction; }

    /**
     * \brief Normalize the direction of the ray.
     * \return A reference to the current ray.
     */
    Ray & normalize() { m_direction.normalize(); return *this; }

    /**
     * \brief Normalize the direction of the ray.
     * \return A copy of the current normalized ray.
     */
    Ray normalized() const { return Ray(m_origin, m_direction.normalized()); }

    /**
     * \brief Check if the current ray is approximately equal to another ray.
     * \param[in] r Ray to compare with.
     * \param[in] tol Tolerance to use for the comparison.
     */
    bool is_approx( Ray const & r, Real const tol = DUMMY_TOL ) const
    { return m_origin.isApprox(r.m_origin, tol) && m_direction.isApprox(r.m_direction, tol); }

    /**
     * \brief Translate the current ray by a given vector.
     * \param[in] t Vector to translate the ray by.
     * \return A reference to the current ray.
     */
    Ray & translate(Vector const & t) { m_origin += t; return *this; }

    /**
     * \brief Translate the current ray by a given vector.
     * \param[in] t Vector to translate the ray by.
     * \return A copy of the current translated ray.
     */
    Ray translated(Vector const & t) const { return Ray(m_origin + t, m_direction); }

    /**
     * \brief Transform the current ray by a given vector.
     * \param[in] t Transformation to apply to the ray.
     * \return A copy of the current transformed ray.
     * \tparam Transform Type of the transformation.
     */
    template <typename Transform>
    Ray transformed(Transform const & t) const
    { return Ray(m_origin.transform(t), m_direction.rotate(t)); }

    /**
     * \brief Transform the current ray by a given vector.
     * \param[in] t Transformation to apply to the ray.
     * \return A reference to the current ray.
     * \tparam Transform Type of the transformation.
     */
    template <typename Transform>
    Ray & transform( Transform const & t ) {
      m_origin.transform(t);
      m_direction.rotate(t);
      return *this;
    }

    /**
     * \brief Check if the point is inside the ray.
     * \param[in] p Point to check.
     * \param[in] tol Tolerance to use for the containment.
     * \return True if the point is inside the ray, false otherwise.
     */
    bool contains( Point const & p, Real tol = DUMMY_TOL ) const {
      Vector v((p - m_origin).normalized());
      return std::abs(v.cross(m_direction).norm()) < tol && v.dot(m_direction) >= -tol;
    }

    /**
     * Check if the current ray intersects a given axis-aligned box.
     * \param[in] b Box to check.
     * \param[in] tol Tolerance to use for the intersection.
     * \return True if the current ray intersects the given box, false otherwise.
     */
    bool intersects( Box<Real, N> const & b, Real tol = DUMMY_TOL ) const {
      if (b.contains(m_origin)) return true;
      Point const & b_min{b.min()};
      Point const & b_max{b.max()};
      Vector t_min, t_max; t_min.setConstant(-MAX), t_max.setConstant(MAX);
      for (Integer i{0}; i < N; ++i) {
        if (std::abs(m_direction[i]) > tol) {
          t_min[i] = (b_min[i] - m_origin[i])/m_direction[i];
          t_max[i] = (b_max[i] - m_origin[i])/m_direction[i];
          if (t_min[i] > t_max[i]) {std::swap(t_min[i], t_max[i]);}
        } else if (m_origin[i] < b_min[i] || m_origin[i] > b_max[i]) {
          return false;
        }
      }
      Real t_entry{t_min.maxCoeff()};
      Real t_exit{t_max.minCoeff()};
      return t_entry <= t_exit && t_exit >= -tol;
    }

    /**
     * \brief Check if the current ray intersects a given axis-aligned box.
     * \param[in] b Box to check.
     * \param[out] c Closest intersection point.
     * \param[out] f Farthest intersection point.
     * \param[in] tol Tolerance to use for the intersection.
     * \return True if the current ray intersects the given box, false otherwise.
     */
    bool intersect( Box<Real, N> const & b, Point & c, Point & f, Real tol = DUMMY_TOL ) const {
      if (b.contains(m_origin)) return true;
      Point const & b_min{b.min()};
      Point const & b_max{b.max()};
      Vector t_min, t_max; t_min.setConstant(-MAX), t_max.setConstant(MAX);
      for (Integer i{0}; i < N; ++i) {
        if (std::abs(m_direction[i]) > tol) {
          t_min[i] = (b_min[i] - m_origin[i])/m_direction[i];
          t_max[i] = (b_max[i] - m_origin[i])/m_direction[i];
          if (t_min[i] > t_max[i]) {std::swap(t_min[i], t_max[i]);}
        } else if (m_origin[i] < b_min[i] || m_origin[i] > b_max[i]) {
          return false;
        }
      }
      Real t_entry{t_min.maxCoeff()};
      Real t_exit{t_max.minCoeff()};
      if (t_entry > t_exit && t_exit < -tol) {return false;}
      c = m_origin + t_entry*m_direction;
      f = m_origin + t_exit*m_direction;
      return true;
    }

    /**
     * \brief Compute the squared distance between the current ray a given point.
     * \param[in] p Point to compute the squared distance to.
     * \param[in] tol Tolerance to use for the intersection.
     * \return The squared distance between the ray and the point.
     */
    Real squared_distance( Point const & p, Real tol = DUMMY_TOL ) const {
      Real t{(p - m_origin).dot(m_direction)/m_direction.squaredNorm()};
      return (m_origin + std::max(static_cast<Real>(-tol), t) * m_direction - p).squaredNorm();
    }

    /**
     * \brief Compute the squared distance between the current ray a given point, returning a point at the given distance.
     * \param[in] p Point to compute the distance to.
     * \param[out] c Closest point on the ray.
     * \param[in] tol Tolerance to use for the intersection.
     * \return The squared distance between the ray and the point.
     */
    Real squared_distance( Point const & p, Point & c, Real tol = DUMMY_TOL ) const
    {
      Real t{(p - m_origin).transpose().dot(m_direction)/m_direction.squaredNorm()};
      c = m_origin + std::max(static_cast<Real>(-tol), t) * m_direction;
      return (c - p).squaredNorm();
    }

    /**
     * \brief Compute the distance between the current ray a given point.
     * \param[in] p Point to compute the distance to.
     * \param[in] tol Tolerance to use for the intersection.
     * \return The distance between the ray and the point.
     */
    Real distance( Point const & p, Real tol = DUMMY_TOL ) const
    { return std::sqrt(this->squared_distance(p, tol)); }

    /**
     * \brief Compute the distance between the current ray a given point, returning a point at the given distance.
     * \param[in] p Point to compute the distance to.
     * \param[out] c Closest point on the ray.
     * \param[in] tol Tolerance to use for the intersection.
     * \return The distance between the ray and the point.
     */
    Real distance( Point const & p, Point & c, Real tol = DUMMY_TOL ) const
    { return std::sqrt(this->squared_distance(p, c, tol)); }


    /**
     * \brief Compute the squared \em interior (or \em minimum) distance between the current ray a given box.
     * \param[in] b Box to compute the squared distance to.
     * \param[in] tol Tolerance to use for the distance computation.
     * \return The squared distance between the ray and the box.
     * \note The squared distance is positive if the ray and the box do not intersect, zero otherwise.
     */
    Real squared_interior_distance( Box<Real, N> const & b, Real tol = DUMMY_TOL ) const
    { Point p1, p2; return this->squared_interior_distance(b, p1, p2, tol); }

    /**
     * \brief Compute the squared \em interior (or \em minimum) distance between the current ray a given box,
     *        returning two points at the minimum distance.
     * \param[in] b Box to compute the squared distance to.
     * \param[out] p1 First point at the minimum distance (on the current ray).
     * \param[out] p2 Second point at the minimum distance (on the box).
     * \param[in] tol Tolerance to use for the distance computation.
     * \return The squared distance between the ray and the box.
     * \note The squared distance is positive if the ray and the box do not intersect, zero otherwise.
     */
    Real squared_interior_distance( Box<Real, N> const & b, Point & p1, Point & p2, Real tol = DUMMY_TOL ) const {
      if (b.contains(m_origin)) return 0;
      Point const & b_min{b.min()};
      Point const & b_max{b.max()};

      // Compute intersection parameters
      Vector t_min, t_max; t_min.setConstant(-MAX); t_max.setConstant(MAX);
      for (Integer i{0}; i < N; ++i) {
        if (std::abs(m_direction[i]) > tol) {
          t_min[i] = (b_min[i] - m_origin[i])/m_direction[i];
          t_max[i] = (b_max[i] - m_origin[i])/m_direction[i];
          if (t_min[i] > t_max[i]) {std::swap(t_min[i], t_max[i]);}
        } else if (m_origin[i] < b_min[i] || m_origin[i] > b_max[i]) {
          // Ray is parallel and outside the box and non-intersecting
          b.interior_distance(m_origin, p2);
          Vector sides{b_max - b_min};
          this->distance(p2, p1);
          for (Integer j{0}; j < N; ++j){
            if (j == i) {continue;}
            if (m_direction[j] > 0.0) {p1[j] += 0.5*sides[j]; p2[j] += 0.5*sides[j];}
            else {p1[j] -= 0.5*sides[j]; p2[j] -= 0.5*sides[j];}
          }
          return (p2 - m_origin).norm();
        }
      }

      // Compute global intersection range
      Real t_entry{t_min.maxCoeff()}, t_exit{t_max.minCoeff()};
      if (t_entry <= t_exit && t_exit >= 0) return 0;

      // Compute closest point if no intersection
      for (Integer i{0}; i < N; ++i) {
        if      (m_origin[i] < b_min[i]) { p2[i] = b_min[i]; }
        else if (m_origin[i] > b_max[i]) { p2[i] = b_max[i]; }
        else                             { p2[i] = m_origin[i]; }
      }

      // Project closest point onto the ray
      Vector v(p2 - m_origin);
      Real t_proj{ v.dot(m_direction)/m_direction.squaredNorm() };
      if ( t_proj < 0 ) { p1 = m_origin; return v.norm(); }

      p1 = m_origin + t_proj * m_direction;;
      return (p2 - p1).squaredNorm();
    }

    /**
     * \brief Compute the \em interior (or \em minimum) distance between the current ray a given box.
     * \param[in] b Box to compute the distance to.
     * \param[in] tol Tolerance to use for the distance computation.
     * \return The distance between the ray and the box.
     * \note The distance is positive if the ray and the box do not intersect, zero otherwise.
     */
    Real interior_distance( Box<Real, N> const & b, Real tol = DUMMY_TOL ) const
    { return std::sqrt(this->squared_interior_distance(b, tol)); }

    /**
     * \brief Compute the \em interior (or \em minimum) distance between the current ray a given box,
     * returning two points at the minimum distance.
     * \param[in] b Box to compute the distance to.
     * \param[out] p1 First point at the minimum distance (on the current ray).
     * \param[out] p2 Second point at the minimum distance (on the box).
     * \param[in] tol Tolerance to use for the distance computation.
     * \return The distance between the ray and the box.
     * \note The distance is positive if the ray and the box do not intersect, zero otherwise.
     */
    Real interior_distance( Box<Real, N> const & b, Point & p1, Point & p2, Real tol = DUMMY_TOL ) const
    { return std::sqrt(static_cast<Real>(this->squared_interior_distance(b, p1, p2, tol))); }

    /**
     * \brief Compute the squared \em exterior (or \em maximum) distance between the current ray a given box.
     * \param[in] b Box to compute the squared distance to.
     * \param[in] tol Tolerance to use for the distance computation.
     * \return The squared distance between the ray and the box.
     * \note The squared distance is positive if one box is not contained in the other, zero otherwise.
     */
    Real squared_exterior_distance( Box<Real, N> const & b, Real tol = DUMMY_TOL ) const
    { Point p1, p2; return this->squared_exterior_distance(b, p1, p2, tol); }

    /**
     * \brief Compute the squared \em exterior (or \em maximum) distance between the current ray a given box,
     *        returning two points at the maximum distance.
     * \param[in] b Box to compute the squared distance to.
     * \param[out] p1 First point at the maximum distance (on the current ray).
     * \param[out] p2 Second point at the maximum distance (on the box).
     * \param[in] tol Tolerance to use for the distance computation.
     * \return The squared distance between the ray and the box.
     */
    Real squared_exterior_distance( Box<Real, N> const & b, Point & p1, Point & p2, Real tol = DUMMY_TOL ) const {
      if (b.contains(m_origin)) return 0;
      Point const & b_min{b.min()};
      Point const & b_max{b.max()};

      // Compute intersection parameters
      Vector t_min, t_max; t_min.setConstant(-MAX); t_max.setConstant(MAX);
      for (Integer i{0}; i < N; ++i) {
        if (std::abs(m_direction[i]) > tol) {
          t_min[i] = (b_min[i] - m_origin[i])/m_direction[i];
          t_max[i] = (b_max[i] - m_origin[i])/m_direction[i];
          if (t_min[i] > t_max[i]) {std::swap(t_min[i], t_max[i]);}
        } else if (m_origin[i] < b_min[i] || m_origin[i] > b_max[i]) {
          // Ray is parallel and outside the box and non-intersecting
          b.exterior_distance(m_origin, p2);
          Vector sides{b_max - b_min};
          this->distance(p2, p1);
          for (Integer j{0}; j < N; ++j){
            if (j == i) {continue;}
            if (m_direction[j] > 0.0) {p1[j] -= 0.5*sides[j]; p2[j] -= 0.5*sides[j];}
            else {p1[j] += 0.5*sides[j]; p2[j] += 0.5*sides[j];}
          }
          return (p2 - m_origin).norm();
        }
      }

      // Compute global intersection range
      Real t_entry{t_min.maxCoeff()}, t_exit{t_max.minCoeff()};
      if (t_entry <= t_exit && t_exit >= 0) {return 0.0;}

      // Compute closest point if no intersection
      for (Integer i{0}; i < N; ++i) {
        if      ( m_origin[i] < b_min[i] ) p2[i] = b_max[i];
        else if ( m_origin[i] > b_max[i] ) p2[i] = b_min[i];
        else                               p2[i] = m_origin[i];
      }

      // Project closest point onto the ray
      Vector v(p2 - m_origin);
      Real t_proj{v.dot(m_direction)/m_direction.squaredNorm()};
      if ( t_proj < 0 ) { p1 = m_origin; return v.norm(); }

      p1 = m_origin + t_proj * m_direction;;
      return (p2 - p1).squaredNorm();
    }

    /**
     * \brief Compute the \em exterior (or \em maximum) distance between the current ray a given box.
     * \param[in] b Box to compute the distance to.
     * \param[in] tol Tolerance to use for the distance computation.
     * \return The distance between the ray and the box.
     * \note The distance is positive if the ray and the box do not intersect, zero otherwise.
     */
    Real exterior_distance( Box<Real, N> const & b, Real tol = DUMMY_TOL ) const
    { return std::sqrt(this->squared_exterior_distance(b, tol)); }

    /**
     * \brief Compute the \em exterior (or \em maximum) distance between the current ray a given box,
     *        and return two points at the maximum distance.
     * \param[in] b Box to compute the distance to.
     * \param[out] p1 First point at the maximum distance (on the current ray).
     * \param[out] p2 Second point at the maximum distance (on the box).
     * \param[in] tol Tolerance to use for the distance computation.
     * \return The distance between the ray and the box.
     * \note The distance is positive if the ray and the box do not intersect, zero otherwise.
     */
    Real exterior_distance( Box<Real, N> const & b, Point & p1, Point & p2, Real tol = DUMMY_TOL ) const
    { return std::sqrt(this->squared_exterior_distance(b, p1, p2, tol)); }

    /**
     * \brief Print the ray info to an output stream.
     * \param[in] os Output stream to print the ray info to.
     */
    void print(std::ostream & os) const {
      os <<
        "────────────────────────────────────────────────────────────────────────────────\n" <<
        "RAY INFO\n" <<
        "\to = " << m_origin.transpose()    << '\n' <<
        "\td = " << m_direction.transpose() << '\n' <<
        "────────────────────────────────────────────────────────────────────────────────\n";
    }

  }; // class Ray

  /**
  * Print the ray info to an output stream.
  * \param[in] os Output stream to print the ray info to.
  * \param[in] r Ray to print.
  * \tparam Real Type of the scalar coefficients.
  * \tparam N Dimension of the ambient space.
  */
  template <typename Real, Integer N>
  std::ostream & operator<<(std::ostream & os, Ray<Real, N> const & r) {r.print(os); return os;}

} // namespace AABBtree

#endif // AABBTREE_Ray_HXX
