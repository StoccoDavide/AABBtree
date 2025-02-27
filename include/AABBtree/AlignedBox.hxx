/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
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
  * \tparam Real the type of the scalar coefficients
  * \tparam N the dimension of the ambient space, can be a compile time value or Dynamic.
  */
  template <typename Real, Integer N>
  class AlignedBox
  {
    static_assert(N > 0, "AlignedBox dimension must be positive.");
    static_assert(std::is_floating_point<Real>::value, "AlignedBox real type must be a floating-point type.");

    constexpr static Real EPS{std::numeric_limits<Real>::epsilon()};
    constexpr static Real MAX{static_cast<Real>(1.0)/EPS};
    constexpr static Real MIN{-MAX};
    constexpr static Real DUMMY_TOL{EPS*static_cast<Real>(100.0)};

    using Vector = Eigen::Vector<Real, N>; /**< Eigen column vector of real numbers. */

    Vector m_min; /**< Minimal corner of the box. */
    Vector m_max; /**< Maximal corner of the box. */

  public:

    /**
    * Cast the current object to a new scalar type.
    * \tparam NewReal the new scalar type
    * \note If \a NewReal is equal to the current scalar type currently used, then this function
    * returns a const reference to the current object.
    */
    template<typename NewReal>
    AlignedBox<NewReal, N> cast() const
    {
      if (std::is_same<Real, NewReal>::value) {
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
    * \warning If either component of \a t_min is larger than the same component of \a t_max, the constructed box is empty.
    */
    AlignedBox(Vector const & t_min, Vector const & t_max) : m_min(t_min), m_max(t_max) {}

    /**
    * Class constructor for a axis-aligned box containing a single point.
    * \param[in] p Point to be contained in the box.
    */
    AlignedBox(const Vector & p) : m_min(p), m_max(p) {}

    /**
    * Copy constructor for a axis-aligned box given another box with a different scalar type.
    * \param[in] b Box to copy.
    */
    template<typename OtherReal>
    explicit AlignedBox(const AlignedBox<OtherReal, N>& b)
    {
      this->m_min = b.m_min.template cast<Real>();
      this->m_max = b.m_max.template cast<Real>();
    }

    /**
    * Class destructor for the axis-aligned box.
    */
    ~AlignedBox() {}

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
    * Check if the current box is degenrate.
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
    void set_degenerate(const Vector & p)
    {
      this->m_min = p;
      this->m_max = p;
    }

    /**
    * Check if the box is empty.
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
    const Vector & min() const {return this->m_min;}

    /**
    * Get a non const reference to the minimal corner.
    * \return The minimal corner of the box.
    */
    Vector & min() {return this->m_min;}

    /**
    * Get the maximal corner.
    * \return The maximal corner of the box.
    */
    const Vector & max() const {return this->m_max;}

    /**
    * Get a non const reference to the maximal corner.
    * \return The maximal corner of the box.
    */
    Vector & max() {return this->m_max;}

    /**
    * Compute the center of the box.
    * \return The center of the box.
    */
    const Real center() const {return (this->m_min + this->m_max)/static_cast<Real>(2.0);}

    /**
    * Compute the lengths of the bounding box's sides.
    * \return The lengths of the bounding box's sides.
    */
    const Vector sizes() const {return this->m_max - this->m_min;}

    /** \return the volume of the bounding box */
    const Real volume() const { return this->sizes().prod(); }

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
    bool contains(const Vector & p) const
    {
      return (this->m_min.array() <= p.array()).all() &&
             (p.array() <= this->m_max.array()).all();
    }

    /**
    * Check if the current box contains a given box.
    * \param[in] b Box to check.
    * \return True if the current box contains the box \a b, false otherwise.
    */
    bool contains(const AlignedBox & b) const
    {
      return (this->m_min.array() <= b.m_min.array()).all() &&
             (b.m_max.array() <= this->m_max.array()).all();
    }

    /**
    * Check if the current box intersects a given box.
    * \param[in] b Box to check.
    * \return True if the current box intersects the box \a b, false otherwise.
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
    AlignedBox & extend(const Vector & p)
    {
      this->m_min = this->m_min.cwiseMin(p);
      this->m_max = this->m_max.cwiseMax(p);
      return *this;
    }

    /**
    * Extends the current box such that it contains a given box.
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
    template<typename Derived>
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
    * \note The squared distance is positive if the point is outside the box, zero otherwise.
    */
    Real squared_exterior_distance(const Vector & p) const
    {
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
    * Compute the squared distance between the current box a given box.
    * \param[in] b Box to compute the squared distance to.
    * \return The squared distance between the boxes.
    * \note The squared distance is positive if the boxes do not intersect, zero otherwise.
    */
    Real squared_exterior_distance(const AlignedBox & b) const
    {
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
    * Compute the distance between the current box a given point.
    * \param[in] p Point to compute the distance to.
    * \return The distance between the box and the point.
    * \note The distance is positive if the point is outside the box, zero otherwise.
    */
    Real exterior_distance(const Vector & p) const
    {
      return std::sqrt(static_cast<Real>(this->squared_exterior_distance(p)));
    }

    /**
    * Compute the distance between the current box a given box.
    * \param[in] b Box to compute the distance to.
    * \return The distance between the boxes.
    * \note The distance is positive if the boxes do not intersect, zero otherwise.
    */
    Real exterior_distance(const AlignedBox & b) const
    {
      return std::sqrt(static_cast<Real>(this->squared_exterior_distance(b)));
    }

  }; // class AlignedBox

} // namespace AABBtree

#endif // AABBTREE_ALIGNEDBOX_HXX
