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
    * Copy constructor for a ray given another box.
    * \param[in] r Ray to copy.
    */
    Ray(Ray const & b) : m_origin(b.m_origin), m_direction(b.m_direction) {}

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
    template <typename T = Real>
    Ray<T, 1>(Real const o, Real const d) : m_origin(o), m_direction(d) {}

    /**
    * Class constructor for the 2D ray.
    * \param[in] o_x Origin \f$ x \f$-axis component.
    * \param[in] o_y Origin \f$ y \f$-axis component.
    * \param[in] d_x Direction \f$ x \f$-axis component.
    * \param[in] d_y Direction \f$ y \f$-axis component.
    * \tparam T Type of the scalar coefficients.
    * \note This constructor is only available for 2D rays.
    */
    template <typename T = Real>
    Ray<T, 2>(Real const o_x, Real const o_y, Real const d_x, Real const d_y)
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
    template <typename T = Real>
    Ray<T, 3>(Real const o_x, Real const o_y, Real const o_z, Real const d_x, Real const d_y, Real const d_z)
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
      return Ray<NewReal, N>(this->m_min.template cast<NewReal>(), this->m_max.template cast<NewReal>());
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
      // If origin is inside the box, then the ray intersects the box
      if (b.contains(this->m_origin, tol)) {return true;}

      // For each dimension, the ray enters and exits the bounding box at t_min and t_max
      Point const & b_min{b.min()};
      Point const & b_max{b.max()};
      Vector t_min, t_max;
      for (Integer i{0}; i < N; ++i) {
        if (std::abs(this->m_direction[i]) > tol) {
          t_min[i] = (b_min[i] - this->m_origin[i])/this->m_direction[i];
          t_max[i] = (b_max[i] - this->m_origin[i])/this->m_direction[i];
          if (t_min[i] > t_max[i]) {std::swap(t_min[i], t_max[i]);}
        } else if (this->m_origin[i] < b_min[i] || this->m_origin[i] > b_max[i]) {
          return false;
        }
      }
      return true;
    }

    /**
    * Check if the current ray intersects a given axis-aligned box.
    * \param[in] b Box to check.
    * \param[out] c Closest intersection point.
    * \param[out] f Farthest intersection point.
    * \param[in] tol Tolerance to use for the intersection.
    * \return True if the current ray intersects the given box, false otherwise.
    */
    bool intersection(Box<Real, N> const & b, Point & c, Point & f, Real tol = DUMMY_TOL) const
    {
      // If origin is inside the box, then the ray intersects the box
      if (b.contains(this->m_origin, tol)) {return true;}

      // For each dimension, the ray enters and exits the bounding box at t_min and t_max
      Point const & b_min{b.min()};
      Point const & b_max{b.max()};
      Vector t_min, t_max; t_min.setConstant(MIN), t_max.setConstant(MAX);
      for (Integer i{0}; i < N; ++i) {
        if (std::abs(this->m_direction[i]) > tol) {
          t_min[i] = (b_min[i] - this->m_origin[i])/this->m_direction[i];
          t_max[i] = (b_max[i] - this->m_origin[i])/this->m_direction[i];
          if (t_min[i] > t_max[i]) {std::swap(t_min[i], t_max[i]);}
        } else if (this->m_origin[i] < b_min[i] || this->m_origin[i] > b_max[i]) {
          return false;
        }
      }
      c = this->m_origin + t_min.maxCoeff()*this->m_direction;
      f = this->m_origin + t_max.minCoeff()*this->m_direction;
      return true;
    }

    /**
    * Compute the squared distance between the current ray a given ray.
    * \param[in] r Ray to compute the squared distance to.
    * \param[in] tol Tolerance to use for the distance computation.
    * \return The squared distance between the rays.
    */
    Real squared_distance(Ray const & r, Real tol = DUMMY_TOL) const
    {
      // Alias variables
      Point const & p1{this->m_origin};
      Vector const & d1{this->m_direction};
      Point const & p2{r.m_origin};
      Vector const & d2{r.m_direction};

      // Compute common terms
      Vector p2_p1(p2 - p1);
      Real d1_d2(d1.dot(d2));
      Real p2_p1_d1(p2_p1.dot(d1));
      Real p2_p1_d2(p2_p1.dot(d2));

      // Compute t and s using least-squares solution
      Real t, s, denom{1.0 - d1_d2*d1_d2};
      if (std::abs(denom) < tol) {
        t = 0.0;
        s = std::max(0.0, p2_p1_d1);
      } else {
        t = std::max(0.0, p2_p1_d1 - p2_p1_d2*d1_d2);
        s = std::max(0.0, p2_p1_d2 - p2_p1_d1*d1_d2);
      }

      // Compute distance
      return (this->m_origin + t*d1 - r.m_origin + s*d2).squaredNorm();
    }

    /**
    * Compute the squared distance between the current ray a given ray, returning a point at the given distance.
    * \param[in] r Ray to compute the squared distance to.
    * \param[out] p1 First point at the minimum distance (on the current ray).
    * \param[out] p2 Second point at the minimum distance (on the given ray).
    * \param[in] tol Tolerance to use for the distance computation.
    * \return The squared distance between the rays.
    */
    Real squared_distance(Ray const & r, Point & p1, Point & p2, Real tol = DUMMY_TOL) const
    {
      // Convert to numpy arrays
      Point const & o1{this->m_origin};
      Vector const & d1{this->m_direction};
      Point const & o2{r.m_origin};
      Vector const & d2{r.m_direction};

      // Compute common terms
      Vector o2_o1(o2 - o1);
      Real d1_d2{d1.dot(d2)};
      Real o2_o1_d1{o2_o1.dot(d1)};
      Real o2_o1_d2{o2_o1.dot(d2)};

      // Compute t and s using least-squares solution
      Real t, s, denom{1.0 - d1_d2*d1_d2};
      if (std::abs(denom) < tol) {
        t = 0.0;
        s = std::max(0.0, o2_o1_d1);
      } else {
        t = std::max(0.0, o2_o1_d1 - o2_o1_d2*d1_d2);
        s = std::max(0.0, o2_o1_d2 - o2_o1_d1*d1_d2);
      }

      // Compute closest points
      p1 = o1 + t*d1;
      p2 = o2 + s*d2;

      // Compute distance
      return (p1 - p2).squaredNorm();
    }

    /**
    * Compute the distance between the current ray a given ray, returning a point at the given distance.
    * \param[in] r Ray to compute the squared distance to.
    * \param[out] c Closest point on the ray.
    * \return The distance between the rays.
    */
    Real distance(Ray const & r, Point & c) const {return std::sqrt(this->squared_distance(r, c));}

    /**
    * Compute the squared distance between the current ray a given point.
    * \param[in] p Point to compute the squared distance to.
    * \return The squared distance between the ray and the point.
    */
    Real squared_distance(Point const & p) const
    {
      Real t{std::max(0.0, (p - this->m_origin).dot(this->m_direction)/this->m_direction.squaredNorm())};
      return (this->m_origin + t * this->m_direction - p).squaredNorm();
    }

    /**
    * Compute the squared distance between the current ray a given point, returning a point at the given distance.
    * \param[in] p Point to compute the distance to.
    * \param[out] c Closest point on the ray.
    * \return The squared distance between the ray and the point.
    */
    Real squared_distance(Point const & p, Point & c) const
    {
      Real t{std::max(0.0, (p - this->m_origin).dot(this->m_direction) / this->m_direction.squaredNorm())};
      c = this->m_origin + t * this->m_direction;
      return (c - p).squaredNorm();
    }

    /**
    * Compute the distance between the current ray a given point.
    * \param[in] p Point to compute the distance to.
    * \return The distance between the ray and the point.
    */
    Real distance(Point const & p) const {return std::sqrt(this->squared_distance(p));}

    /**
    * Compute the distance between the current ray a given point, returning a point at the given distance.
    * \param[in] p Point to compute the distance to.
    * \param[out] c Closest point on the ray.
    * \return The distance between the ray and the point.
    */
    Real distance(Point const & p, Point & c) const {return std::sqrt(this->squared_distance(p, c));}

    /**
    * Compute the squared distance between the current ray a given axis-aligned box.
    * \param[in] b Box to compute the squared distance to.
    * \return The squared distance between the ray and the axis-aligned box.
    */
    Real squared_distance(Box<Real, N> const & b, Real tol = DUMMY_TOL) const
    {
      if (b.contains(this->m_origin, tol)) {return -1.0;}

      // For each dimension, the ray enters and exits the bounding box at t_min and t_max
      Point const & b_min{b.min()};
      Point const & b_max{b.max()};
      Vector t_min, t_max; t_min.setConstant(MIN), t_max.setConstant(MAX);
      for (Integer i{0}; i < N; ++i) {
        if (std::abs(this->m_direction[i]) > tol) {
          t_min[i] = (b_min[i] - this->m_origin[i])/this->m_direction[i];
          t_max[i] = (b_max[i] - this->m_origin[i])/this->m_direction[i];
          if (t_min[i] > t_max[i]) {std::swap(t_min[i], t_max[i]);}
        } else if (this->m_origin[i] < b_min[i] || this->m_origin[i] > b_max[i]) {
          return false;
        }
      }
      return (this->m_origin + t_min.maxCoeff()*this->m_direction -
              this->m_origin + t_max.minCoeff()*this->m_direction).norm();
    }

    /**
    * Compute the distance between the current ray a given axis-aligned box, returning a point at the given distance.
    * \param[in] b Box to compute the distance to.
    * \param[out] p1 First point at the minimum distance (on the current ray).
    * \param[out] p2 Second point at the minimum distance (on the given box).
    * \return The distance between the ray and the axis-aligned box. The returned value is positive
    * if the ray origin is outside the box, negative otherwise.
    */
    Real distance(Box<Real, N> const & b, Point & p1, Point & p2, Real tol = DUMMY_TOL) const
    {
      if (b.contains(this->m_origin, tol)) {return -1.0;}

      // For each dimension, the ray enters and exits the bounding box at t_min and t_max
      Point const & b_min{b.min()};
      Point const & b_max{b.max()};
      Vector t_min, t_max; t_min.setConstant(MIN), t_max.setConstant(MAX);
      for (Integer i{0}; i < N; ++i) {
        if (std::abs(this->m_direction[i]) > tol) {
          t_min[i] = (b_min[i] - this->m_origin[i])/this->m_direction[i];
          t_max[i] = (b_max[i] - this->m_origin[i])/this->m_direction[i];
          if (t_min[i] > t_max[i]) {std::swap(t_min[i], t_max[i]);}
        } else if (this->m_origin[i] < b_min[i] || this->m_origin[i] > b_max[i]) {
          return false;
        }
      }
      p1 = this->m_origin + t_min.maxCoeff()*this->m_direction;
      p2 = this->m_origin + t_max.minCoeff()*this->m_direction;
      return (p1 - p2).norm();
    }

  }; // class Ray

} // namespace AABBtree

#endif // AABBTREE_Ray_HXX
