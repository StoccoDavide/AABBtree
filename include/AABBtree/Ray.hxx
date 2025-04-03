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
  * \brief Class container for a ray.
  *
  * This class represents a ray as a pair of an origin and a direction.
  * \tparam Real Type of the scalar coefficients
  * \tparam N Dimension of the ambient space.
  */
  template <typename Real, Integer N>
  class Ray {
    static_assert(std::is_floating_point<Real>::value, "Ray Real type must be a floating-point type.");
    static_assert(N > 0, "Ray dimension must be positive.");

    constexpr static Real EPS{std::numeric_limits<Real>::epsilon()}; /**> Machine epsilon for the scalar type. */
    constexpr static Real MAX{static_cast<Real>(1.0)/EPS}; /**> Maximum value for the scalar type. */
    constexpr static Real MIN{-MAX}; /**> Minimum value for the scalar type. */
    constexpr static Real DUMMY_TOL{EPS*static_cast<Real>(100.0)}; /**> Dummy tolerance for the scalar type. */

    using Point = Point<Real, N>;
    using Vector = Vector<Real, N>;

    Point m_origin; /**< Origin of the ray. */
    Vector m_direction; /**< Direction of the ray. */

  public:
   /**
    * Class destructor for the ray.
    */
    ~Ray() = default;

    /**
    * Class constructor for a ray.
    */
    Ray() = default;

    /**
    * Copy constructor for a ray.
    * \param[in] r Ray to copy.
    */
    Ray(Ray const & r) : m_origin(r.m_origin), m_direction(r.m_direction) {}

    /**
    * Class constructor for a ray given an origin and a direction.
    * \param[in] o Origin of the ray.
    * \param[in] d Direction of the ray.
    */
    Ray(Point const & o, Vector const & d) : m_origin(o), m_direction(d) {}

    /**
    * Class constructor for the 1D ray.
    * \param[in] o Origin of the ray.
    * \param[in] d Direction of the ray.
    * \tparam T Type of the scalar coefficients.
    * \note This constructor is only available for 1D rays.
    */
    template <typename = std::enable_if<N == 1>>
    Ray(Real const o, Real const d) : m_origin(o), m_direction(d) {}

    /**
    * Class constructor for the 2D ray.
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
    * Class constructor for the 3D ray.
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
    * Cast the current object to a new scalar type.
    * \tparam NewReal the new scalar type
    * \note If the new real type is equal to the current scalar type currently used, then this
    * function returns a const reference to the current object.
    */
    template<typename NewReal>
    Ray<NewReal, N> cast() const
    {
      if constexpr (std::is_same<Real, NewReal>::value) {return *this;}
      return Ray<NewReal, N>(this->m_origin.template cast<NewReal>(), this->m_direction.template cast<NewReal>());
    }

    /**
    * Get the reference to the ray origin.
    * \return The reference to the ray origin.
    */
    Point & origin() {return this->m_origin;}

    /**
    * Get the const reference to the ray origin.
    * \return The const reference to the ray origin.
    */
    Point const & origin() const {return this->m_origin;}

    /**
    * Get the reference to the ray direction.
    * \return The reference to the ray direction.
    */
    Vector & direction() {return this->m_direction;}

    /**
    * Get the const reference to the ray direction.
    * \return The const reference to the ray direction.
    */
    Vector const & direction() const {return this->m_direction;}

    /**
    * Normalize the direction of the ray.
    * \return A reference to the current ray.
    */
    Ray & normalize() {this->m_direction.normalize(); return *this;}

    /**
    * Normalize the direction of the ray.
    * \return A copy of the current normalized ray.
    */
    Ray normalized() const {return Ray(this->m_origin, this->m_direction.normalized());}

    /**
    * Check if the current ray is approximately equal to another ray.
    * \param[in] r Ray to compare with.
    * \param[in] tol Tolerance to use for the comparison.
    */
    bool is_approx(Ray const & r, Real const tol = DUMMY_TOL) const
    {return this->m_origin.isApprox(r.m_origin, tol) && this->m_direction.isApprox(r.m_direction, tol);}

    /**
    * Translate the current ray by a given vector.
    * \param[in] t Vector to translate the ray by.
    * \return A reference to the current ray.
    */
    Ray & translate(Vector const & t) {this->m_origin += t; return *this;}

    /**
    * Translate the current ray by a given vector.
    * \param[in] t Vector to translate the ray by.
    * \return A copy of the current translated ray.
    */
    Ray translated(Vector const & t) const {return Ray(this->m_origin + t, this->m_direction);}

    /**
    * Transform the current ray by a given vector.
    * \param[in] t Transformation to apply to the ray.
    * \return A copy of the current transformed ray.
    * \tparam Transform Type of the transformation.
    */
    template <typename Transform>
    Ray transformed(Transform const & t) const
    {return Ray(this->m_origin.transform(t), this->m_direction.rotate(t));}

    /**
    * Transform the current ray by a given vector.
    * \param[in] t Transformation to apply to the ray.
    * \return A reference to the current ray.
    * \tparam Transform Type of the transformation.
    */
    template <typename Transform>
    Ray & transform(Transform const & t) {
      this->m_origin.transform(t);
      this->m_direction.rotate(t);
      return *this;
    }

    /**
    * Check if the point is inside the ray.
    * \param[in] p Point to check.
    * \param[in] tol Tolerance to use for the containment.
    * \return True if the point is inside the ray, false otherwise.
    */
    bool contains(Point const & p, Real tol = DUMMY_TOL) const
    {
      Vector v((p - this->m_origin).normalized());
      return std::abs(v.cross(this->m_direction).norm()) < tol && v.dot(this->m_direction) >= -tol;
    }

    /**
    * Check if the current ray intersects a given axis-aligned box.
    * \param[in] b Box to check.
    * \param[in] tol Tolerance to use for the intersection.
    * \return True if the current ray intersects the given box, false otherwise.
    */
    bool intersects(Box<Real, N> const & b, Real tol = DUMMY_TOL) const
    {
      if (b.contains(this->m_origin)) {return true;}
      Point const & b_min{b.min()};
      Point const & b_max{b.max()};
      Vector t_min, t_max; t_min.setConstant(-MAX), t_max.setConstant(MAX);
      for (Integer i{0}; i < N; ++i) {
        if (std::abs(this->m_direction[i]) > tol) {
          t_min[i] = (b_min[i] - this->m_origin[i])/this->m_direction[i];
          t_max[i] = (b_max[i] - this->m_origin[i])/this->m_direction[i];
          if (t_min[i] > t_max[i]) {std::swap(t_min[i], t_max[i]);}
        } else if (this->m_origin[i] < b_min[i] || this->m_origin[i] > b_max[i]) {
          return false;
        }
      }
      Real t_entry{t_min.maxCoeff()};
      Real t_exit{t_max.minCoeff()};
      return t_entry <= t_exit && t_exit >= -tol;
    }

    /**
    * Check if the current ray intersects a given axis-aligned box.
    * \param[in] b Box to check.
    * \param[out] c Closest intersection point.
    * \param[out] f Farthest intersection point.
    * \param[in] tol Tolerance to use for the intersection.
    * \return True if the current ray intersects the given box, false otherwise.
    */
    bool intersect(Box<Real, N> const & b, Point & c, Point & f, Real tol = DUMMY_TOL) const
    {
      if (b.contains(this->m_origin)) {return true;}
      Point const & b_min{b.min()};
      Point const & b_max{b.max()};
      Vector t_min, t_max; t_min.setConstant(-MAX), t_max.setConstant(MAX);
      for (Integer i{0}; i < N; ++i) {
        if (std::abs(this->m_direction[i]) > tol) {
          t_min[i] = (b_min[i] - this->m_origin[i])/this->m_direction[i];
          t_max[i] = (b_max[i] - this->m_origin[i])/this->m_direction[i];
          if (t_min[i] > t_max[i]) {std::swap(t_min[i], t_max[i]);}
        } else if (this->m_origin[i] < b_min[i] || this->m_origin[i] > b_max[i]) {
          return false;
        }
      }
      Real t_entry{t_min.maxCoeff()};
      Real t_exit{t_max.minCoeff()};
      if (t_entry > t_exit && t_exit < -tol) {return false;}
      c = this->m_origin + t_entry*this->m_direction;
      f = this->m_origin + t_exit*this->m_direction;
      return true;
    }

    /**
    * Compute the squared distance between the current ray a given point.
    * \param[in] p Point to compute the squared distance to.
    * \param[in] tol Tolerance to use for the intersection.
    * \return The squared distance between the ray and the point.
    */
    Real squared_distance(Point const & p, Real tol = DUMMY_TOL) const
    {
      Real t{(p - this->m_origin).dot(this->m_direction)/this->m_direction.squaredNorm()};
      return (this->m_origin + std::max(static_cast<Real>(-tol), t) * this->m_direction - p).squaredNorm();
    }

    /**
    * Compute the squared distance between the current ray a given point, returning a point at the given distance.
    * \param[in] p Point to compute the distance to.
    * \param[out] c Closest point on the ray.
    * \param[in] tol Tolerance to use for the intersection.
    * \return The squared distance between the ray and the point.
    */
    Real squared_distance(Point const & p, Point & c, Real tol = DUMMY_TOL) const
    {
      Real t{(p - this->m_origin).transpose().dot(this->m_direction)/this->m_direction.squaredNorm()};
      c = this->m_origin + std::max(static_cast<Real>(-tol), t) * this->m_direction;
      return (c - p).squaredNorm();
    }

    /**
    * Compute the distance between the current ray a given point.
    * \param[in] p Point to compute the distance to.
    * \param[in] tol Tolerance to use for the intersection.
    * \return The distance between the ray and the point.
    */
    Real distance(Point const & p, Real tol = DUMMY_TOL) const
    {return std::sqrt(this->squared_distance(p, tol));}

    /**
    * Compute the distance between the current ray a given point, returning a point at the given distance.
    * \param[in] p Point to compute the distance to.
    * \param[out] c Closest point on the ray.
    * \param[in] tol Tolerance to use for the intersection.
    * \return The distance between the ray and the point.
    */
    Real distance(Point const & p, Point & c, Real tol = DUMMY_TOL) const
    {return std::sqrt(this->squared_distance(p, c, tol));}


    /**
    * Compute the squared \em interior (or \em minimum) distance between the current ray a given box.
    * \param[in] b Box to compute the squared distance to.
    * \param[in] tol Tolerance to use for the distance computation.
    * \return The squared distance between the ray and the box.
    * \note The squared distance is positive if the ray and the box do not intersect, zero otherwise.
    */
    Real squared_interior_distance(Box<Real, N> const & b, Real tol = DUMMY_TOL) const
    {Point p1, p2; return this->squared_interior_distance(b, p1, p2, tol);}

    /**
    * Compute the squared \em interior (or \em minimum) distance between the current ray a given box,
    * returning two points at the minimum distance.
    * \param[in] b Box to compute the squared distance to.
    * \param[out] p1 First point at the minimum distance (on the current ray).
    * \param[out] p2 Second point at the minimum distance (on the box).
    * \param[in] tol Tolerance to use for the distance computation.
    * \return The squared distance between the ray and the box.
    * \note The squared distance is positive if the ray and the box do not intersect, zero otherwise.
    */
    Real squared_interior_distance(Box<Real, N> const & b, Point & p1, Point & p2, Real tol = DUMMY_TOL) const
    {
      if (b.contains(this->m_origin)) {return 0.0;}
      Point const & b_min{b.min()};
      Point const & b_max{b.max()};

      // Compute intersection parameters
      Vector t_min, t_max; t_min.setConstant(-MAX); t_max.setConstant(MAX);
      for (Integer i{0}; i < N; ++i) {
        if (std::abs(this->m_direction[i]) > tol) {
          t_min[i] = (b_min[i] - this->m_origin[i])/this->m_direction[i];
          t_max[i] = (b_max[i] - this->m_origin[i])/this->m_direction[i];
          if (t_min[i] > t_max[i]) {std::swap(t_min[i], t_max[i]);}
        } else if (this->m_origin[i] < b_min[i] || this->m_origin[i] > b_max[i]) {
          // Ray is parallel and outside the box and non-intersecting
          b.interior_distance(this->m_origin, p2);
          Vector sides{b_max - b_min};
          this->distance(p2, p1);
          for (Integer j{0}; j < N; ++j){
            if (j == i) {continue;}
            if (this->m_direction[j] > 0.0) {p1[j] += 0.5*sides[j]; p2[j] += 0.5*sides[j];}
            else {p1[j] -= 0.5*sides[j]; p2[j] -= 0.5*sides[j];}
          }
          return (p2 - this->m_origin).norm();
        }
      }

      // Compute global intersection range
      Real t_entry{t_min.maxCoeff()}, t_exit{t_max.minCoeff()};
      if (t_entry <= t_exit && t_exit >= 0) {return 0.0;}

      // Compute closest point if no intersection
      for (Integer i{0}; i < N; ++i) {
        if (this->m_origin[i] < b_min[i]) {p2[i] = b_min[i];}
        else if (this->m_origin[i] > b_max[i]) {p2[i] = b_max[i];}
        else {p2[i] = this->m_origin[i];}
      }

      // Project closest point onto the ray
      Vector v(p2 - this->m_origin);
      Real t_proj{v.dot(this->m_direction)/this->m_direction.squaredNorm()};
      if (t_proj < 0.0) {p1 = this->m_origin; return v.norm();}

      p1 = this->m_origin + t_proj * this->m_direction;;
      return (p2 - p1).squaredNorm();
    }

    /**
    * Compute the \em interior (or \em minimum) distance between the current ray a given box.
    * \param[in] b Box to compute the distance to.
    * \param[in] tol Tolerance to use for the distance computation.
    * \return The distance between the ray and the box.
    * \note The distance is positive if the ray and the box do not intersect, zero otherwise.
    */
    Real interior_distance(Box<Real, N> const & b, Real tol = DUMMY_TOL) const
    {return std::sqrt(this->squared_interior_distance(b, tol));}

    /**
    * Compute the \em interior (or \em minimum) distance between the current ray a given box,
    * returning two points at the minimum distance.
    * \param[in] b Box to compute the distance to.
    * \param[out] p1 First point at the minimum distance (on the current ray).
    * \param[out] p2 Second point at the minimum distance (on the box).
    * \param[in] tol Tolerance to use for the distance computation.
    * \return The distance between the ray and the box.
    * \note The distance is positive if the ray and the box do not intersect, zero otherwise.
    */
    Real interior_distance(Box<Real, N> const & b, Point & p1, Point & p2, Real tol = DUMMY_TOL) const
    {return std::sqrt(static_cast<Real>(this->squared_interior_distance(b, p1, p2, tol)));}

    /**
    * Compute the squared \em exterior (or \em maximum) distance between the current ray a given box.
    * \param[in] b Box to compute the squared distance to.
    * \param[in] tol Tolerance to use for the distance computation.
    * \return The squared distance between the ray and the box.
    * \note The squared distance is positive if one box is not contained in the other, zero otherwise.
    */
    Real squared_exterior_distance(Box<Real, N> const & b, Real tol = DUMMY_TOL) const
    {Point p1, p2; return this->squared_exterior_distance(b, p1, p2, tol);}

    /**
    * Compute the squared \em exterior (or \em maximum) distance between the current ray a given box,
    * returning two points at the maximum distance.
    * \param[in] b Box to compute the squared distance to.
    * \param[out] p1 First point at the maximum distance (on the current ray).
    * \param[out] p2 Second point at the maximum distance (on the box).
    * \param[in] tol Tolerance to use for the distance computation.
    * \return The squared distance between the ray and the box.
    */
    Real squared_exterior_distance(Box<Real, N> const & b, Point & p1, Point & p2, Real tol = DUMMY_TOL) const
    {
      if (b.contains(this->m_origin)) {return 0.0;}
      Point const & b_min{b.min()};
      Point const & b_max{b.max()};

      // Compute intersection parameters
      Vector t_min, t_max; t_min.setConstant(-MAX); t_max.setConstant(MAX);
      for (Integer i{0}; i < N; ++i) {
        if (std::abs(this->m_direction[i]) > tol) {
          t_min[i] = (b_min[i] - this->m_origin[i])/this->m_direction[i];
          t_max[i] = (b_max[i] - this->m_origin[i])/this->m_direction[i];
          if (t_min[i] > t_max[i]) {std::swap(t_min[i], t_max[i]);}
        } else if (this->m_origin[i] < b_min[i] || this->m_origin[i] > b_max[i]) {
          // Ray is parallel and outside the box and non-intersecting
          b.exterior_distance(this->m_origin, p2);
          Vector sides{b_max - b_min};
          this->distance(p2, p1);
          for (Integer j{0}; j < N; ++j){
            if (j == i) {continue;}
            if (this->m_direction[j] > 0.0) {p1[j] -= 0.5*sides[j]; p2[j] -= 0.5*sides[j];}
            else {p1[j] += 0.5*sides[j]; p2[j] += 0.5*sides[j];}
          }
          return (p2 - this->m_origin).norm();
        }
      }

      // Compute global intersection range
      Real t_entry{t_min.maxCoeff()}, t_exit{t_max.minCoeff()};
      if (t_entry <= t_exit && t_exit >= 0) {return 0.0;}

      // Compute closest point if no intersection
      for (Integer i{0}; i < N; ++i) {
        if (this->m_origin[i] < b_min[i]) {p2[i] = b_max[i];}
        else if (this->m_origin[i] > b_max[i]) {p2[i] = b_min[i];}
        else {p2[i] = this->m_origin[i];}
      }

      // Project closest point onto the ray
      Vector v(p2 - this->m_origin);
      Real t_proj{v.dot(this->m_direction)/this->m_direction.squaredNorm()};
      if (t_proj < 0.0) {p1 = this->m_origin; return v.norm();}

      p1 = this->m_origin + t_proj * this->m_direction;;
      return (p2 - p1).squaredNorm();
    }

    /**
    * Compute the \em exterior (or \em maximum) distance between the current ray a given box.
    * \param[in] b Box to compute the distance to.
    * \param[in] tol Tolerance to use for the distance computation.
    * \return The distance between the ray and the box.
    * \note The distance is positive if the ray and the box do not intersect, zero otherwise.
    */
    Real exterior_distance(Box<Real, N> const & b, Real tol = DUMMY_TOL) const
    {return std::sqrt(this->squared_exterior_distance(b, tol));}

    /**
    * Compute the \em exterior (or \em maximum) distance between the current ray a given box,
    * and return two points at the maximum distance.
    * \param[in] b Box to compute the distance to.
    * \param[out] p1 First point at the maximum distance (on the current ray).
    * \param[out] p2 Second point at the maximum distance (on the box).
    * \param[in] tol Tolerance to use for the distance computation.
    * \return The distance between the ray and the box.
    * \note The distance is positive if the ray and the box do not intersect, zero otherwise.
    */
    Real exterior_distance(Box<Real, N> const & b, Point & p1, Point & p2, Real tol = DUMMY_TOL) const
    {return std::sqrt(this->squared_exterior_distance(b, p1, p2, tol));}

    /**
    * Print the ray info to an output stream.
    * \param[in] os Output stream to print the ray info to.
    */
    void print(std::ostream & os) const
    {
      os <<
        "────────────────────────────────────────────────────────────────────────────────" << std::endl <<
        "RAY INFO" << std::endl <<
        "\to = " << this->m_origin.transpose() << std::endl <<
        "\td = " << this->m_direction.transpose()  << std::endl <<
        "────────────────────────────────────────────────────────────────────────────────" << std::endl;
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
