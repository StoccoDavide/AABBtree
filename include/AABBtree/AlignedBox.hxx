/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The AABBtree project is distributed under the BSD 2-Clause License.                           *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// DISCLAIMER: The code in this file is a modified version of the Eigen library.

#pragma once

#ifndef AABBTREE_ALIGNEDBOX_HXX
#define AABBTREE_ALIGNEDBOX_HXX

namespace AABBtree
{

  /**
  * \brief Class container for an axis-aligned box.
  *
  * This class represents an axis-aligned box as a pair of opposite corners, i.e., a minimal corner
  * and a maximal corner.
  * \tparam Real Type of the scalar coefficients
  * \tparam N Dimension of the ambient space.
  */
  template <typename Real, Integer N>
  class AlignedBox
  {
    static_assert(N > 0, "AlignedBox dimension must be positive.");
    static_assert(std::is_floating_point<Real>::value, "AlignedBox real type must be a floating-point type.");

    constexpr static Real EPS{std::numeric_limits<Real>::epsilon()}; /**> Machine epsilon for the scalar type. */
    constexpr static Real MAX{static_cast<Real>(1.0)/EPS}; /**> Maximum value for the scalar type. */
    constexpr static Real MIN{-MAX}; /**> Minimum value for the scalar type. */
    constexpr static Real DUMMY_TOL{EPS*static_cast<Real>(100.0)}; /**> Dummy tolerance for the scalar type. */

  public:
    using Vector = Eigen::Vector<Real, N>; /**> Eigen column vector of real numbers. */
    using Point = Eigen::Vector<Real, N>; /**> Point in the ambient space (Eigen column vector of real numbers). */

  private:
    Point m_min; /**< Minimal corner of the box. */
    Point m_max; /**< Maximal corner of the box. */

  public:
    /**
    * Class destructor for the axis-aligned box.
    */
    ~AlignedBox() {}

    /**
    * Cast the current object to a new scalar type.
    * \tparam NewReal the new scalar type
    * \note If the new real typeis equal to the current scalar type currently used, then this
    * function returns a const reference to the current object.
    */
    template<typename NewReal>
    AlignedBox<NewReal, N> cast() const
    {
      if constexpr (std::is_same<Real, NewReal>::value) {
        return *this;
      } else {
        return AlignedBox<NewReal, N>(this->m_min.template cast<NewReal>(), this->m_max.template cast<NewReal>());
      }
    }

    /**
    * Class constructor for an empty axis-aligned box.
    */
    AlignedBox() {this->set_empty();}

    /**
    * Class constructor for a axis-aligned box given the minimal and maximal corners.
    * \param[in] t_min Minimal corner of the box.
    * \param[in] t_max Maximal corner of the box.
    * \warning If the reordering is not performed, and the minimal corner is greater than the maximal
    * corner in any dimension, undefined behavior may occur.
    */
    AlignedBox(const Point & t_min, const Point & t_max) : m_min(t_min), m_max(t_max) {}

    /**
    * Class constructor for a axis-aligned box containing a single point.
    * \param[in] p Point to be contained in the box.
    */
    AlignedBox(const Point & p) : m_min(p), m_max(p) {}

    /**
    * Class constructor for the 1D axis-aligned bounding box.
    * \param[in] x_min Minimal \f$ x \f$-axis corner of the box.
    * \param[in] x_max Maximal \f$ x \f$-axis corner of the box.
    * \tparam T Type of the scalar coefficients.
    * \note This constructor is only available for 1D boxes.
    */
    template <typename T = Real>
    AlignedBox<T, 1>(Real x_min, Real x_max) : m_min(x_min), m_max(x_max) {}

    /**
    * Class constructor for the 2D axis-aligned bounding box.
    * \param[in] x_min Minimal \f$ x \f$-axis corner of the box.
    * \param[in] y_min Minimal \f$ y \f$-axis corner of the box.
    * \param[in] x_max Maximal \f$ x \f$-axis corner of the box.
    * \param[in] y_max Maximal \f$ y \f$-axis corner of the box.
    * \tparam T Type of the scalar coefficients.
    * \note This constructor is only available for 2D boxes.
    */
    template <typename T = Real>
    AlignedBox<T, 2>(Real x_min, Real y_min, Real x_max, Real y_max)
      : m_min(x_min, y_min), m_max(x_max, y_max) {}

    /**
    * Class constructor for the 3D axis-aligned bounding box.
    * \param[in] x_min Minimal \f$ x \f$-axis corner of the box.
    * \param[in] y_min Minimal \f$ y \f$-axis corner of the box.
    * \param[in] z_min Minimal \f$ z \f$-axis corner of the box.
    * \param[in] x_max Maximal \f$ x \f$-axis corner of the box.
    * \param[in] y_max Maximal \f$ y \f$-axis corner of the box.
    * \param[in] z_max Maximal \f$ z \f$-axis corner of the box.
    * \tparam T Type of the scalar coefficients.
    * \note This constructor is only available for 3D boxes.
    */
    template <typename T = Real>
    AlignedBox<T, 3>(Real x_min, Real y_min, Real z_min, Real x_max, Real y_max, Real z_max)
      : m_min(x_min, y_min, z_min), m_max(x_max, y_max, z_max) {}

    /**
    * Class constructor for the 4D axis-aligned bounding box.
    * \param[in] x_min Minimal \f$ x \f$-axis corner of the box.
    * \param[in] y_min Minimal \f$ y \f$-axis corner of the box.
    * \param[in] z_min Minimal \f$ z \f$-axis corner of the box.
    * \param[in] w_min Minimal \f$ w \f$-axis corner of the box.
    * \param[in] x_max Maximal \f$ x \f$-axis corner of the box.
    * \param[in] y_max Maximal \f$ y \f$-axis corner of the box.
    * \param[in] z_max Maximal \f$ z \f$-axis corner of the box.
    * \param[in] w_max Maximal \f$ w \f$-axis corner of the box.
    * \tparam T Type of the scalar coefficients.
    * \note This constructor is only available for 4D boxes.
    */
    template <typename T = Real>
    AlignedBox<T, 4>(Real x_min, Real y_min, Real z_min, Real w_min, Real x_max, Real y_max, Real
      z_max, Real w_max) : m_min(x_min, y_min, z_min, w_min), m_max(x_max, y_max, z_max, w_max) {}

    /**
    * Copy constructor for a axis-aligned box given another box with a different scalar type.
    * \param[in] b Box to copy.
    * \tparam OtherReal Type of the scalar coefficients of the box to copy.
    */
    template<typename OtherReal>
    explicit AlignedBox(const AlignedBox<OtherReal, N>& b)
    {
      this->m_min = b.m_min.template cast<Real>();
      this->m_max = b.m_max.template cast<Real>();
    }

    /**
    * Reorder the corners of the box such that the minimal corner is less than the maximal corner.
    */
    void reorder()
    {
      for (Integer i{0}; i < N; ++i)
      {
        if (this->m_min[i] > this->m_max[i]) {std::swap(this->m_min[i], this->m_max[i]);}
      }
    }

    /**
    * Check if the current box is approximately equal to another box.
    * \param[in] b Box to compare with.
    * \param[in] tol Tolerance to use for the comparison.
    */
    bool is_approx(const AlignedBox & b, const Real tol = DUMMY_TOL) const
    {
      return this->m_min.isApprox(b.m_min, tol) && this->m_max.isApprox(b.m_max, tol);
    }

    /**
    * Check if the current box is degenrate, i.e., it has zero volume.
    * \param[in] tol Tolerance to use for the comparison.
    * \return True if the box is degenerate, false otherwise.
    */
    bool is_degenerate(const Real tol = DUMMY_TOL) const
    {
      return this->m_min.isApprox(this->m_max, tol);
    }

    /**
    * Set the box to be degenerate.
    * \param[in] p Point to set the box to.
    */
    void set_degenerate(const Point & p)
    {
      this->m_min = p;
      this->m_max = p;
    }

    /**
    * Check if the box is empty, i.e., the minimal corner is greater than the maximal corner in any
    * dimension.
    * \return True if the box is empty, false otherwise.
    */
    bool is_empty() const
    {
      return (this->m_min.array() > this->m_max.array()).any();
    }

    /**
    * Set the box to be empty.
    */
    void set_empty()
    {
      this->m_min.setConstant(MAX);
      this->m_max.setConstant(MIN);
    }

    /**
    * Get the minimal corner.
    * \return The minimal corner of the box.
    */
    const Point & min() const {return this->m_min;}

    /**
    * Get a non const reference to the minimal corner.
    * \return The minimal corner of the box.
    */
    Point & min() {return this->m_min;}

    /**
    * Get the /f$ i /f$-th component of the minimal corner.
    * \param[in] i Index of the component.
    * \return The /f$ i /f$-th component of the minimal corner.
    */
    Real min(Integer i) const {return this->m_min[i];}

    /**
    * Get the maximal corner.
    * \return The maximal corner of the box.
    */
    const Point & max() const {return this->m_max;}

    /**
    * Get a non const reference to the maximal corner.
    * \return The maximal corner of the box.
    */
    Point & max() {return this->m_max;}

    /**
    * Get the /f$ i /f$-th component of the maximal corner.
    * \param[in] i Index of the component.
    * \return The /f$ i /f$-th component of the maximal corner.
    */
    Real max(Integer i) const {return this->m_max[i];}

    /**
    * Get the longest axis of the box.
    * \return The longest axis of the box.
    */
    Integer longest_axis() const
    {
      Vector sizes{this->sizes()};
      return std::distance(sizes.data(), std::max_element(sizes.data(), sizes.data() + N));
    }

    /**
    * Compute the center of the box.
    * \return The center of the box.
    */
    const Point center() const {return (this->m_min + this->m_max)/static_cast<Real>(2.0);}

    /**
    * Compute the lengths of the bounding box's sides.
    * \return The lengths of the bounding box's sides.
    */
    const Vector sizes() const {return this->m_max - this->m_min;}

    /**
    * Compute the volume of the bounding box.
    * \return The volume of the bounding box.
    */
    const Real volume() const {return this->sizes().prod();}

    /**
    * Compute the external surface area of the bounding box.
    * \return The external surface area of the bounding box.
    */
    const Real surface() const
    {
      Vector sizes{this->sizes()};
      Real area{static_cast<Real>(0.0)}, prod;
      for (Integer i{0}, j; i < N; ++i)
      {
        prod = static_cast<Real>(1.0);
        for (j = 0; j < N; ++j) {if (j != i) {prod *= sizes[j];}}
        area += static_cast<Real>(2.0)*prod;
      }
      return area;
    }

    /**
    * Compute the bounding box diagonal vector.
    * \note If the length of the diagonal is needed, \c diagonal().norm() will provide it.
    * \return The diagonal vector of the bounding box.
    */
    const Vector diagonal() const {return this->sizes();}

    /**
    * Check if the point is inside the box.
    * \param[in] p Point to check.
    * \return True if the point is inside the box, false otherwise.
    */
    bool contains(const Point & p) const
    {
      return (this->m_min.array() <= p.array()).all() &&
             (p.array() <= this->m_max.array()).all();
    }

    /**
    * Check if the current box contains a given box.
    * \param[in] b Box to check.
    * \return True if the current box contains the given box, false otherwise.
    */
    bool contains(const AlignedBox & b) const
    {
      return (this->m_min.array() <= b.m_min.array()).all() &&
             (b.m_max.array() <= this->m_max.array()).all();
    }

    /**
    * Check if the current box intersects a given box.
    * \param[in] b Box to check.
    * \return True if the current box intersects the given box, false otherwise.
    */
    bool intersects(const AlignedBox & b) const
    {
      return (this->m_min.array() <= b.m_max.array()).all() &&
             (b.m_min.array() <= this->m_max.array()).all();
    }

    /**
    * Extend the current box such that it contains a given point.
    * \param[in] p Point to extend the box to.
    * \return A reference to the current box.
    */
    template<typename Derived>
    AlignedBox & extend(const Point & p)
    {
      this->m_min = this->m_min.cwiseMin(p);
      this->m_max = this->m_max.cwiseMax(p);
      return *this;
    }

    /**
    * Extend the current box such that it contains the given points.
    * \param[in] p Points to extend the box to.
    * \return A reference to the current box.
    */
    AlignedBox & extend(const std::vector<Point> & p)
    {
      for (const Point & point : p) {this->extend(point);}
      return *this;
    }

    /**
    * Extend the current box such that it contains a given box.
    * \param[in] b Box to extend the box to.
    * \return A reference to the current box.
    * \note Merging with an empty box may result in a box bigger than \c *this.
    */
    AlignedBox & extend(const AlignedBox & b)
    {
      this->m_min = this->m_min.cwiseMin(b.m_min);
      this->m_max = this->m_max.cwiseMax(b.m_max);
      return *this;
    }

    /**
    * Extend the current box such that it contains the given boxes.
    * \param[in] b Boxes to extend the box to.
    * \return A reference to the current box.
    * \note Merging with an empty box may result in a box bigger than \c *this.
    */
    AlignedBox & extend(const std::vector<AlignedBox> & b)
    {
      for (const AlignedBox & box : b) {this->extend(box);}
      return *this;
    }

    /**
    * Clamp the current box by a given box .
    * \param[in] b Box to extend the box to.
    * \return A reference to the current box.
    * \note If the boxes do not intersect, the resulting box is empty.
    */
    AlignedBox & clamp(const AlignedBox & b)
    {
      this->m_min = this->m_min.cwiseMax(b.m_min);
      this->m_max = this->m_max.cwiseMin(b.m_max);
      return *this;
    }

    /**
    * Return the intersection of the current box and another box.
    * \param[in] b Box to intersect with.
    * \return The intersection of the current box and the given box.
    * \note If the boxes do not intersect, the resulting box is empty.
    */
    AlignedBox intersection(const AlignedBox & b) const
    {
      return AlignedBox(this->m_min.cwiseMax(b.m_min), this->m_max.cwiseMin(b.m_max));
    }

    /**
    * Return the union of the current box and another box.
    * \param[in] b Box to merge with.
    * \return The union of the current box and the given box.
    */
    AlignedBox merged(const AlignedBox & b) const
    {
      return AlignedBox(this->m_min.cwiseMin(b.m_min), this->m_max.cwiseMax(b.m_max));
    }

    /**
    * Translate the current box by a given vector.
    * \param[in] t Vector to translate the box by.
    * \return A reference to the current box.
    */
    AlignedBox & translate(const Vector & t)
    {
      this->m_min += t;
      this->m_max += t;
      return *this;
    }

    /**
    * Translate the current box by a given vector.
    * \param[in] t Vector to translate the box by.
    * \return A copy of the current translated box.
    */
    AlignedBox translated(const Vector & t) const
    {
      AlignedBox result(this->m_min, this->m_max);
      result.translate(t);
      return result;
    }

    /**
    * Compute the squared distance between the current box a given point.
    * \param[in] p Point to compute the squared distance to.
    * \return The squared distance between the box and the point.
    * \note The returned value is positive if the point is outside the box, zero otherwise.
    */
    Real squared_distance(const Point & p) const
    {
      if (this->contains(p)) {return static_cast<Real>(0.0);}
      Real dist2{0.0}, aux;
      for (Integer i{0}; i < N; ++i)
      {
        if (this->m_min[i] > p[i]) {
          aux = this->m_min[i] - p[i];
          dist2 += aux*aux;
        } else if (p[i] > this->m_max[i]) {
          aux = p[i] - this->m_max[i];
          dist2 += aux*aux;
        }
      }
      return dist2;
    }

    /**
    * Compute the squared distance between the current box a given point, returning two points at
    * the given distance.
    * \param[in] p Point to compute the squared distance to.
    * \param[out] p_max First point at the distance.
    * \param[out] p_min Second point at the distance.
    * \return The squared distance between the box and the point.
    * \note The returned value is positive if the point is outside the box, zero otherwise.
    */
    Real squared_distance(const Point & p, Point & p_min, Point & p_max) const
    {
      if (this->contains(p)) {return static_cast<Real>(0.0);}
      p_min = p.cwiseMin(this->m_min);
      p_max = p.cwiseMax(this->m_max);
      return (p_max - p_min).squaredNorm();
    }

    /**
    * Compute the distance between the current box a given point.
    * \param[in] p Point to compute the distance to.
    * \return The distance between the box and the point.
    * \note The distance is positive if the point is outside the box, zero otherwise.
    */
    Real distance(const Point & p) const
    {
      return std::sqrt(static_cast<Real>(this->squared_distance(p)));
    }

    /**
    * Compute the distance between the current box a given point, returning two points at the given
    * distance.
    * \param[in] p Point to compute the distance to.
    * \param[out] p_max First point at the distance.
    * \param[out] p_min Second point at the distance.
    * \return The distance between the box and the point.
    * \note The returned value is positive if the point is outside the box, zero otherwise.
    */
    Real distance(const Point & p, Point & p_min, Point & p_max) const
    {
      return std::sqrt(static_cast<Real>(this->squared_distance(p, p_min, p_max)));
    }

    /**
    * Compute the squared \em interior (or \em minimum) distance between the current box a given box.
    * \param[in] b Box to compute the squared distance to.
    * \return The squared distance between the boxes.
    * \note The squared distance is positive if the boxes do not intersect, zero otherwise.
    */
    Real squared_interior_distance(const AlignedBox & b) const
    {
      if (this->intersects(b)) {return static_cast<Real>(0.0);}
      Real dist2{0.0}, aux;
      for (Integer i{0}; i < N; ++i)
      {
        if (this->m_min[i] > b.m_max[i]) {
          aux = this->m_min[i] - b.m_max[i];
          dist2 += aux*aux;
        } else if (b.m_min[i] > this->m_max[i]) {
          aux = b.m_min[i] - this->m_max[i];
          dist2 += aux*aux;
        }
      }
      return dist2;
    }

    /**
    * Compute the squared \em interior (or \em minimum) distance between the current box a given box,
    * returning two points at the minimum distance.
    * \param[in] b Box to compute the squared distance to.
    * \param[out] p_max First point at the minimum distance.
    * \param[out] p_min Second point at the minimum distance.
    * \return The squared distance between the boxes.
    * \note The squared distance is positive if the boxes do not intersect, zero otherwise.
    */
    Real squared_interior_distance(const AlignedBox & b, Point & p_min, Point & p_max) const
    {
      if (this->intersects(b)) {return static_cast<Real>(0.0);}
      for (Integer i{0}; i < N; ++i)
      {
        if (this->m_min[i] > b.m_max[i]) {
          p_min[i] = b.m_max[i];
          p_max[i] = this->m_min[i];
        } else if (b.m_min[i] > this->m_max[i]) {
          p_min[i] = this->m_max[i];
          p_max[i] = b.m_min[i];
        } else {
          p_min[i] = p_max[i] =
            (std::min(this->m_max[i], b.m_max[i]) + std::max(this->m_min[i], b.m_min[i]))/
            static_cast<Real>(2.0);
        }
      }
      return (p_max - p_min).squaredNorm();
    }

    /**
    * Compute the \em interior (or \em minimum) distance between the current box a given box.
    * \param[in] b Box to compute the distance to.
    * \return The distance between the boxes.
    * \note The distance is positive if the boxes do not intersect, zero otherwise.
    */
    Real interior_distance(const AlignedBox & b) const
    {
      return std::sqrt(static_cast<Real>(this->squared_interior_distance(b)));
    }

    /**
    * Compute the \em interior (or \em minimum) distance between the current box a given box,
    * returning two points at the minimum distance.
    * \param[in] b Box to compute the distance to.
    * \param[out] p_max First point at the minimum distance.
    * \param[out] p_min Second point at the minimum distance.
    * \return The distance between the boxes.
    * \note The distance is positive if the boxes do not intersect, zero otherwise.
    */
    Real interior_distance(const AlignedBox & b, Point & p_min, Point & p_max) const
    {
      return std::sqrt(static_cast<Real>(this->squared_interior_distance(b, p_min, p_max)));
    }

    /**
    * Compute the squared \em exterior (or \em maximum) distance between the current box a given box.
    * \param[in] b Box to compute the squared distance to.
    * \return The squared distance between the boxes.
    * \note The squared distance is positive if one box is not contained in the other, zero otherwise.
    */
    Real squared_exterior_distance(const AlignedBox & b) const
    {
      Real dist2{0.0}, aux;
      for (Integer i{0}; i < N; ++i)
      {
        aux = std::max(this->m_max[i], b.m_max[i]) - std::min(this->m_min[i], b.m_min[i]);
        dist2 += aux*aux;
      }
      return dist2;
    }

    /**
    * Compute the squared \em exterior (or \em maximum) distance between the current box a given box,
    * returning two points at the maximum distance.
    * \param[in] b Box to compute the squared distance to.
    * \param[out] p_max First point at the maximum distance.
    * \param[out] p_min Second point at the maximum distance.
    * \return The squared distance between the boxes.
    */
    Real squared_exterior_distance(const AlignedBox & b, Point & p_min, Point & p_max) const
    {
      p_min = b.m_min.cwiseMin(this->m_min);
      p_max = b.m_max.cwiseMax(this->m_max);
      return (p_max - p_min).squaredNorm();
    }

    /**
    * Compute the \em exterior (or \em maximum) distance between the current box a given box.
    * \param[in] b Box to compute the distance to.
    * \return The distance between the boxes.
    * \note The distance is positive if the boxes do not intersect, zero otherwise.
    */
    Real exterior_distance(const AlignedBox & b) const
    {
      return std::sqrt(static_cast<Real>(this->squared_exterior_distance(b)));
    }

    /**
    * Compute the \em exterior (or \em maximum) distance between the current box a given box,
    * and return two points at the maximum distance.
    * \param[in] b Box to compute the distance to.
    * \param[out] p_max First point at the maximum distance.
    * \param[out] p_min Second point at the maximum distance.
    * \return The distance between the boxes.
    * \note The distance is positive if the boxes do not intersect, zero otherwise.
    */
    Real exterior_distance(const AlignedBox & b, Point & p_min, Point & p_max) const
    {
      return std::sqrt(static_cast<Real>(this->squared_exterior_distance(b, p_min, p_max)));
    }

  }; // class AlignedBox

  template class AlignedBox<float, 1>;
  template class AlignedBox<float, 2>;
  template class AlignedBox<float, 3>;
  template class AlignedBox<float, 4>;
  template class AlignedBox<float, 5>;
  template class AlignedBox<double, 1>;
  template class AlignedBox<double, 2>;
  template class AlignedBox<double, 3>;
  template class AlignedBox<double, 4>;
  template class AlignedBox<double, 5>;

} // namespace AABBtree

#endif // AABBTREE_ALIGNEDBOX_HXX
