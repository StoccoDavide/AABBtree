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

#ifndef AABBTREE_BOX_HXX
#define AABBTREE_BOX_HXX

#include "AABBtree/Ray.hxx"

namespace AABBtree {

  /**
   * \brief Forward declaration of the Ray class template
   * \tparam Real Type of the scalar coefficients
   * \tparam N Dimension of the ambient space
   */
  template <typename Real, Integer N> class Ray;

  /**
   * \class Box
   * \brief Axis-aligned bounding box (AABB) implementation
   *
   * This class represents an axis-aligned bounding box (AABB) defined by two opposite corners:
   * - Minimum corner (lower bounds in all dimensions)
   * - Maximum corner (upper bounds in all dimensions)
   *
   * The Box class provides various geometric operations including:
   * - Intersection and containment tests
   * - Distance calculations
   * - Transformation operations
   * - Volume and surface area computations
   *
   * \tparam Real Floating-point scalar type (must be float or double)
   * \tparam N Dimension of the space (1, 2, 3 or more)
   *
   * \note The box is considered empty if any component of the minimum corner is greater than
   *       the corresponding component of the maximum corner.
   * \note All methods are thread-safe unless explicitly noted otherwise.
   */
  template <typename Real, Integer N>
  class Box {
    static_assert(std::is_floating_point<Real>::value, "Box Real type must be a floating-point type.");
    static_assert(N > 0, "Box dimension must be positive.");

    // Constants for numerical computations
    constexpr static Real EPS{std::numeric_limits<Real>::epsilon()}; ///< Machine epsilon for the scalar type.
    constexpr static Real MAX{std::numeric_limits<Real>::max()};     ///< Maximum representable value for the scalar type.
    constexpr static Real MIN{-MAX};                                 ///< Minimum representable value for the scalar type.
    constexpr static Real DUMMY_TOL{EPS*static_cast<Real>(100.0)};   ///< Default tolerance for the scalar type

    using Point            = Point<Real, N>;
    using Vector           = Vector<Real, N>;
    using BoxUniquePtr     = BoxUniquePtr<Real, N>;
    using BoxUniquePtrList = BoxUniquePtrList<Real, N>;

    Point m_min; ///< Minimal corner of the box (lower bounds in all dimensions).
    Point m_max; ///< Maximal corner of the box (upper bounds in all dimensions).

  public:
    /**
     * Class destructor for the axis-aligned box.
     */
    ~Box() = default;

    /**
     * \brief Default constructor creating an empty box
     *
     * Creates an empty box where:
     * - Minimum corner is set to MAX in all dimensions
     * - Maximum corner is set to MIN in all dimensions
     */
    Box()
    {
      this->m_min.setConstant(MAX);
      this->m_max.setConstant(MIN);
    }

    /**
     * Copy constructor for a axis-aligned box given another box.
     * \param[in] b Box to copy.
     */
    Box( Box const & b ) : m_min(b.m_min), m_max(b.m_max) {}

    /**
     * Class constructor for a axis-aligned box given the minimal and maximal corners.
     * \param[in] t_min Minimal corner of the box.
     * \param[in] t_max Maximal corner of the box.
     * \warning If the reordering is not performed, and the minimal corner is greater
     *          than the  maximal corner in any dimension, undefined behavior may occur.
     */
    Box( Point const & t_min, Point const & t_max ) : m_min(t_min), m_max(t_max) {}

    /**
     * Class constructor for a axis-aligned box containing a single point.
     * \param[in] p Point to be contained in the box.
     */
    Box(Point const & p) : m_min(p), m_max(p) {}

    /**
     * Class constructor for the 1D axis-aligned bounding box.
     * \param[in] x_min Minimal \f$ x \f$-axis corner of the box.
     * \param[in] x_max Maximal \f$ x \f$-axis corner of the box.
     * \tparam T Type of the scalar coefficients.
     * \note This constructor is only available for 1D boxes.
     */
    template <typename = std::enable_if<N == 1>>
    Box( Real const x_min, Real const x_max ) : m_min(x_min), m_max(x_max) {}

    /**
     * Class constructor for the 2D axis-aligned bounding box.
     * \param[in] x_min Minimal \f$ x \f$-axis corner of the box.
     * \param[in] y_min Minimal \f$ y \f$-axis corner of the box.
     * \param[in] x_max Maximal \f$ x \f$-axis corner of the box.
     * \param[in] y_max Maximal \f$ y \f$-axis corner of the box.
     * \tparam T Type of the scalar coefficients.
     * \note This constructor is only available for 2D boxes.
     */
    template <typename = std::enable_if<N == 2>>
    Box( Real const x_min, Real const y_min, Real const x_max, Real const y_max )
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
    template <typename = std::enable_if<N == 3>>
    Box(
      Real const x_min, Real const y_min, Real const z_min,
      Real const x_max, Real const y_max, Real const z_max
    ) : m_min(x_min, y_min, z_min), m_max(x_max, y_max, z_max) {}

    /**
     * Copy constructor for a axis-aligned box given another box with a different scalar type.
     * \param[in] b Box to copy.
     * \tparam OtherReal Type of the scalar coefficients of the box to copy.
     */
    template<typename OtherReal>
    explicit Box( Box<OtherReal, N> const & b )
    : m_min(b.m_min.template cast<Real>()), m_max(b.m_max.template cast<Real>()) {}

    /**
     * Cast the current object to a new scalar type.
     * \tparam NewReal the new scalar type
     * \note If the new real type is equal to the current scalar type currently used,
     *       then this function returns a const reference to the current object.
     */
    template<typename NewReal>
    Box<NewReal, N> cast() const
    {
      if constexpr (std::is_same<Real, NewReal>::value) return *this;
      return Box<NewReal,N>(this->m_min.template cast<NewReal>(),
                            this->m_max.template cast<NewReal>());
    }

    /**
     * Get the minimal corner.
     * \return The minimal corner of the box.
     */
    Point const & min() const { return this->m_min; }

    /**
     * Get a non const reference to the minimal corner.
     * \return The minimal corner of the box.
     */
    Point & min() { return this->m_min; }

    /**
     * Get the maximal corner.
     * \return The maximal corner of the box.
     */
    Point const & max() const { return this->m_max; }

    /**
     * Get a non const reference to the maximal corner.
     * \return The maximal corner of the box.
     */
    Point & max() { return this->m_max; }

    /**
     * Reorder the corners of the box such that the minimal corner
     * is less than the maximal  corner.
     */
    void reorder() {
      for (Integer i{0}; i < N; ++i) {
        if ( this->m_min[i] > this->m_max[i] )
          std::swap(this->m_min[i], this->m_max[i]);
      }
    }

    /**
     * Check if the box is empty, i.e., the minimal corner is greater than
     * the maximal corner in any dimension.
     * \return True if the box is empty, false otherwise.
     */
    bool is_empty() const { return (this->m_max.array() < this->m_min.array()).any(); }

    /**
     * \brief Check if this box is approximately equal to another box
     * \param[in] b Box to compare with
     * \param[in] tol Tolerance for the comparison
     * \return True if boxes are approximately equal within tolerance
     */
    bool is_approx(Box const & b, Real const tol /*= DUMMY_TOL*/) const
    { return this->m_min.isApprox(b.m_min, tol) && this->m_max.isApprox(b.m_max, tol); }

    /**
     * Check if the current box is degenrate, i.e., it has zero volume.
     * \param[in] tol Tolerance to use for the comparison.
     * \return True if the box is degenerate, false otherwise.
     */
    bool is_degenerate(Real const tol /*= DUMMY_TOL*/) const
    { return this->m_min.isApprox(this->m_max, tol); }

    /**
     * \brief Find the longest axis of the box
     * \return Index of the dimension with the largest extent
     */
    Integer longest_axis() const
    {
      Vector sizes{this->sizes()};
      return std::distance(sizes.data(), std::max_element(sizes.data(), sizes.data() + N));
    }

    /**
     * \brief Find the longest axis and compute its length and midpoint
     * \param[out] max_length Length of the longest axis
     * \param[out] mid_point Midpoint of the longest axis
     * \return Index of the dimension with the largest extent
     */
    Integer longest_axis(Real & max_length, Real & mid_point) const
    {
      Integer axis{-1};
      max_length = -1;
      for (Integer i{0}; i < N; ++i) {
        if (Real const tmp_length{this->m_max(i) - this->m_min(i)}; max_length < tmp_length)
          {max_length = tmp_length; mid_point = 0.5*(this->m_max(i) + this->m_min(i)); axis = i;}
      }
      return axis;
    }

    /**
     * \brief Sort axes by their lengths in descending order
     * \param[out] sizes Vector to store the sizes of each dimension
     * \param[out] ipos Vector to store the sorted dimension indices
     */
    void sort_axes_length(Vector & sizes, Eigen::Vector<Integer, N> & ipos) const
    {
      sizes = this->m_max - this->m_min;
      ipos = Eigen::Vector<Integer, N>::LinSpaced(N, 0, N - 1);
      std::sort(ipos.data(), ipos.data() + N, [&sizes](Integer i, Integer j) {return sizes[i] > sizes[j];});
    }

    /**
     * \brief Compute the barycenter of the box
     * \return Point at the center of the box
     */
    Point baricenter() const { return 0.5*(this->m_min + this->m_max); }

    /**
     * \brief Compute the barycenter coordinate along a specific dimension
     * \param[in] i Dimension index
     * \return Barycenter coordinate along dimension i
     */
    Real baricenter(Integer i) const {return 0.5*(this->m_min(i) + this->m_max(i));}

    /**
     * \brief Compute the sizes of the box in all dimensions
     * \return Vector containing the size along each dimension
     */
    Vector sizes() const {return this->m_max - this->m_min;}

    /**
     * \brief Compute the volume of the box
     * \return Product of sizes in all dimensions
     */
    Real volume() const {return this->sizes().prod();}

    /**
     * \brief Compute the surface area of the box
     * \return Sum of areas of all faces
     */
    Real surface() const
    {
      Vector sizes{this->sizes()};
      Real area{0.0};
      for (Integer i{0}; i < N; ++i) {
        Real prod{1.0};
        for (Integer j{0}; j < N; ++j)
          if (j != i)
            prod *= sizes[j];
        area += 2.0*prod;
      }
      return area;
    }

    /**
     * \brief Compute the diagonal vector of the box
     * \return Vector from min to max corner
     * \note For the length of the diagonal, use diagonal().norm()
     */
    Vector diagonal() const {return this->sizes();}

    /**
     * \brief Set the box to be degenerate (containing a single point)
     * \param[in] p Point to set the box to
     */
    void set_degenerate(Point const & p) {
      this->m_min = p;
      this->m_max = p;
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
     * \brief Check if this box intersects with another box
     * \param[in] b Box to check for intersection
     * \return True if boxes intersect, false otherwise
     */
    bool intersects(Box const & b) const {
      return (this->m_min.array() <= b.m_max.array()).all() &&
             (b.m_min.array() <= this->m_max.array()).all();
    }

    /**
     * \brief Compute the intersection of this box with another box
     * \param[in] b_in Box to intersect with
     * \param[out] b_out Resulting intersection box
     * \return True if intersection exists, false otherwise
     */
    bool intersect(Box const & b_in, Box & b_out) const
    {
      b_out.m_min = this->m_min.cwiseMax(b_in.m_min);
      b_out.m_max = this->m_max.cwiseMin(b_in.m_max);
      return (b_out.m_min.array() <= b_out.m_max.array()).all();
    }

    /**
     * \brief Compute the union of this box with another box
     * \param[in] b Box to merge with
     * \return The smallest box containing both boxes
     */
    Box merged(Box const & b) const
    { return Box(this->m_min.cwiseMin(b.m_min), this->m_max.cwiseMax(b.m_max)); }

    /**
     * \brief Translate the box by a given vector
     * \param[in] t Translation vector
     * \return Reference to this box after translation
     */
    Box & translate(Vector const & t) { this->m_min += t; this->m_max += t; return *this; }

    /**
     * \brief Create a translated copy of this box
     * \param[in] t Translation vector
     * \return New translated box
     */
    Box translated(Vector const & t) const {return Box(this->m_min + t, this->m_max + t);}

    /**
     * \brief Create a transformed copy of this box
     * \param[in] t Transformation to apply
     * \return New transformed box
     * \tparam Transform Type of the transformation
     * \note The resulting box may be larger than the original due to axis-alignment
     */
    template <typename Transform>
    Box transformed(Transform const & t) const
    {
      Box result( this->m_min.transform(t), this->m_max.transform(t) );
      result.reorder();
      return result;
    }

    /**
     * \brief Transform this box in-place
     * \param[in] t Transformation to apply
     * \return Reference to this box after transformation
     * \tparam Transform Type of the transformation
     * \note The resulting box may be larger than the original due to axis-alignment
     */
    template <typename Transform>
    Box & transform(Transform const & t) {
      this->m_min.transform(t);
      this->m_max.transform(t);
      this->reorder();
      return *this;
    }

    /**
     * \enum Side
     * \brief Enumeration for relative position with respect to the box
     */
    enum class Side : Integer {
      LEFT = -1,
      INSIDE = 0,
      RIGHT = +1
    };

    /**
     * \brief Classify a coordinate relative to the box in a specific dimension
     * \param[in] x Coordinate to classify
     * \param[in] tol Tolerance for classification
     * \param[in] dim Dimension to consider
     * \return Side::LEFT if x < min[dim]+tol,
     *         Side::RIGHT if x > max[dim]-tol,
     *         Side::INSIDE otherwise
     */
    Side which_side( Real const x, Real const tol, Integer const dim ) const {
      if (x < this->m_min(dim) + tol) {return Side::LEFT;}
      if (x > this->m_max(dim) - tol) {return Side::RIGHT;}
      return Side::INSIDE;
    }

    /**
     * \brief Check if a point is inside the box
     * \param[in] p Point to check
     * \return True if point is inside or on boundary, false otherwise
     */
    bool contains(Point const & p) const {
      return (this->m_min.array() <= p.array()).all() &&
             (p.array() <= this->m_max.array()).all();
    }

    /**
     * \brief Check if a point intersects the box
     * \param[in] p Point to check
     * \return True if point is inside or on boundary, false otherwise
     * \note Equivalent to contains()
     */
    bool intersects(Point const & p) const {
      return (this->m_min.array() <= p.array()).all() &&
             (p.array() <= this->m_max.array()).all();
    }

    /**
     * \brief Check if this box completely contains another box
     * \param[in] b Box to check
     * \return True if b is entirely within this box, false otherwise
     */
    bool contains(Box const & b) const {
      return (this->m_min.array() <= b.m_min.array()).all() &&
             (b.m_max.array() <= this->m_max.array()).all();
    }

    /**
     * \brief Extend this box to contain a point
     * \param[in] p Point to include
     * \return Reference to this box after extension
     */
    template<typename Derived>
    Box & extend( Point const & p ) {
      this->m_min = this->m_min.cwiseMin(p);
      this->m_max = this->m_max.cwiseMax(p);
      return *this;
    }

    /**
     * \brief Extend this box to contain another box
     * \param[in] b Box to include
     * \return Reference to this box after extension
     * \note Merging with an empty box may result in a box bigger than *this
     */
    Box & extend( Box const & b ) {
      this->m_min = this->m_min.cwiseMin(b.m_min);
      this->m_max = this->m_max.cwiseMax(b.m_max);
      return *this;
    }

    /**
     * Compute the squared \em interior (or \em minimum) distance between the
     * current box a given point.
     * \param[in] p Point to compute the squared distance to.
     * \return The squared distance between the box and the point.
     * \note The returned value is positive if the point is outside the box, zero otherwise.
     */
    Real squared_interior_distance(Point const & p) const
    {
      Real dist2{0.0};
      for (Integer i{0}; i < N; ++i) {
        if (this->m_min[i] > p[i]) {Real aux{this->m_min[i] - p[i]}; dist2 += aux*aux;}
        else if (p[i] > this->m_max[i]) {Real aux{p[i] - this->m_max[i]}; dist2 += aux*aux;}
      }
      return dist2;
    }

    /**
     * Compute the squared \em interior (or \em minimum) distance between
     * the current box a given point,
     * returning a point at the given distance.
     * \param[in] p Point to compute the squared distance to.
     * \param[out] c Closest point to the box.
     * \return The squared distance between the box and the point.
     * \note The returned value is positive if the point is outside the box, zero otherwise.
     */
    Real squared_interior_distance(Point const & p, Point & c) const {
      if (this->contains(p)) return 0;
      c = p.cwiseMax(this->m_min).cwiseMin(this->m_max);
      return (p - c).squaredNorm();
    }

    /**
     * Compute the \em interior (or \em minimum) distance between the current box a given point.
     * \param[in] p Point to compute the distance to.
     * \return The distance between the box and the point.
     * \note The distance is positive if the point is outside the box, zero otherwise.
     */
    Real interior_distance( Point const & p ) const
    { return std::sqrt(this->squared_interior_distance(p)); }

    /**
     * Compute the \em interior (or \em minimum) distance between the current box a given point,
     * returning a point at the given distance.
     * \param[in] p Point to compute the distance to.
     * \param[out] c Closest point to the box.
     * \return The distance between the box and the point.
     * \note The returned value is positive if the point is outside the box, zero otherwise.
     */
    Real interior_distance( Point const & p, Point & c ) const
    { return std::sqrt(this->squared_interior_distance(p, c)); }

    /**
     * Compute the squared \em exterior (or \em maximum) distance between the current box a  given point.
     * \param[in] p Point to compute the squared distance to.
     * \return The squared distance between the box and the point.
     * \note The returned value is positive if the point is outside the box, zero otherwise.
     */
    Real squared_exterior_distance( Point const & p ) const {
      Real dist2{Real(0)};
      for (Integer i{0}; i < N; ++i) {
        Real const aux1{ std::abs(this->m_min[i] - p[i]) };
        Real const aux2{ std::abs(this->m_max[i] - p[i]) };
        Real const aux3{ std::max(aux1,aux2) };
        dist2 += aux3*aux3;
      }
      return dist2;
    }

    /**
     * Compute the squared \em exterior (or \em maximum) distance between the current box a given point,
     * returning a point at the given distance.
     * \param[in] p Point to compute the squared distance to.
     * \param[out] f Farthest point to the box.
     * \return The squared distance between the box and the point.
     * \note The returned value is positive if the point is outside the box, zero otherwise.
     */
    Real squared_exterior_distance( Point const & p, Point & f ) const {
      for (Integer i{0}; i < N; ++i) {
        Real const aux1{ std::abs(this->m_min[i] - p[i]) };
        Real const aux2{ std::abs(this->m_max[i] - p[i]) };
        f[i] = aux1 > aux2 ? this->m_min[i] : this->m_max[i];
      }
      return (f - p).squaredNorm();
    }

    /**
     * Compute the \em exterior (or \em maximum) distance between the current box a given point.
     * \param[in] p Point to compute the distance to.
     * \return The distance between the box and the point.
     * \note The distance is positive if the point is outside the box, zero otherwise.
     */
    Real exterior_distance( Point const & p ) const
    { return std::sqrt(this->squared_exterior_distance(p)); }

    /**
     * Compute the \em exterior (or \em maximum) distance between the current box a given point,
     * returning a point at the given distance.
     * \param[in] p Point to compute the distance to.
     * \param[out] f Farthest point to the box.
     * \return The distance between the box and the point.
     * \note The returned value is positive if the point is outside the box, zero otherwise.
     */
    Real exterior_distance( Point const & p, Point & f ) const
    { return std::sqrt(this->squared_exterior_distance(p, f)); }

    /**
     * Compute the squared \em interior (or \em minimum) distance between the current box a given box.
     * \param[in] b Box to compute the squared distance to.
     * \return The squared distance between the boxes.
     * \note The squared distance is positive if the boxes do not intersect, zero otherwise.
     */
    Real squared_interior_distance( Box const & b ) const {
      if ( this->intersects(b) ) return 0;
      Real dist2{Real(0)};
      for ( Integer i{0}; i < N; ++i ) {
        Real aux{ this->m_min[i] - b.m_max[i] };
        if ( aux < 0 ) {
          aux = b.m_min[i] - this->m_max[i];
          if ( aux < 0 ) continue;
        }
        dist2 += aux*aux;
      }
      return dist2;
    }

    /**
     * Compute the squared \em interior (or \em minimum) distance between the current box a given box,
     * returning two points at the minimum distance.
     * \param[in] b Box to compute the squared distance to.
     * \param[out] p1 First point at the minimum distance.
     * \param[out] p2 Second point at the minimum distance.
     * \return The squared distance between the boxes.
     * \note The squared distance is positive if the boxes do not intersect, zero otherwise.
     */
    Real squared_interior_distance( Box const & b, Point & p1, Point & p2 ) const {
      if ( this->intersects(b) ) return 0;
      for ( Integer i{0}; i < N; ++ i) {
        if      ( this->m_min[i] > b.m_max[i] ) { p1[i] = b.m_max[i];     p2[i] = this->m_min[i]; }
        else if ( b.m_min[i] > this->m_max[i] ) { p1[i] = this->m_max[i]; p2[i] = b.m_min[i];     }
        else {
          p1[i] = p2[i] = 0.5*(std::min(this->m_max[i], b.m_max[i]) + std::max(this->m_min[i], b.m_min[i]));
        }
      }
      return (p2 - p1).squaredNorm();
    }

    /**
     * Compute the \em interior (or \em minimum) distance between the current box a given box.
     * \param[in] b Box to compute the distance to.
     * \return The distance between the boxes.
     * \note The distance is positive if the boxes do not intersect, zero otherwise.
     */
    Real interior_distance( Box const & b ) const { return std::sqrt(this->squared_interior_distance(b)); }

    /**
     * Compute the \em interior (or \em minimum) distance between the current box a given box,
     * returning two points at the minimum distance.
     * \param[in] b Box to compute the distance to.
     * \param[out] p1 First point at the minimum distance.
     * \param[out] p2 Second point at the minimum distance.
     * \return The distance between the boxes.
     * \note The distance is positive if the boxes do not intersect, zero otherwise.
     */
    Real interior_distance( Box const & b, Point & p1, Point & p2 ) const
    { return std::sqrt( this->squared_interior_distance(b, p1, p2) ); }

    /**
     * Compute the squared \em exterior (or \em maximum) distance between the current box a given box.
     * \param[in] b Box to compute the squared distance to.
     * \return The squared distance between the boxes.
     * \note The squared distance is positive if one box is not contained in the other, zero otherwise.
     */
    Real squared_exterior_distance( Box const & b ) const {
      Real dist2{ Real(0) };
      for ( Integer i{0}; i < N; ++i ) {
        Real const aux{ std::max(this->m_max[i], b.m_max[i]) - std::min(this->m_min[i], b.m_min[i]) };
        dist2 += aux*aux;
      }
      return dist2;
    }

    /**
     * Compute the squared \em exterior (or \em maximum) distance between the current box a given box,
     * returning two points at the maximum distance.
     * \param[in] b Box to compute the squared distance to.
     * \param[out] p1 First point at the maximum distance.
     * \param[out] p2 Second point at the maximum distance.
     * \return The squared distance between the boxes.
     */
    Real squared_exterior_distance( Box const & b, Point & p1, Point & p2 ) const {
      p1 = b.m_min.cwiseMin(this->m_min);
      p2 = b.m_max.cwiseMax(this->m_max);
      return (p2 - p1).squaredNorm();
    }

    /**
     * Compute the \em exterior (or \em maximum) distance between the current box a given box.
     * \param[in] b Box to compute the distance to.
     * \return The distance between the boxes.
     * \note The distance is positive if the boxes do not intersect, zero otherwise.
     */
    Real exterior_distance( Box const & b ) const
    { return std::sqrt(this->squared_exterior_distance(b)); }

    /**
     * Compute the \em exterior (or \em maximum) distance between the current box a given box,
     * and return two points at the maximum distance.
     * \param[in] b Box to compute the distance to.
     * \param[out] p1 First point at the maximum distance.
     * \param[out] p2 Second point at the maximum distance.
     * \return The distance between the boxes.
     * \note The distance is positive if the boxes do not intersect, zero otherwise.
     */
    Real exterior_distance( Box const & b, Point & p1, Point & p2 ) const
    { return std::sqrt(this->squared_exterior_distance(b, p1, p2)); }

    /**
     * Check if the current box intersects a given ray.
     * \param[in] r Ray to check.
     * \param[in] tol Tolerance to use for the intersection.
     * \return True if the current box intersects the given ray, false otherwise.
     */
    bool intersects( Ray<Real, N> const & r, Real tol = DUMMY_TOL ) const
    { return r.intersects(*this, tol); }

    /**
     * Check if the current box intersects a given ray.
     * \param[in] r Ray to check.
     * \param[out] c Closest intersection point.
     * \param[out] f Farthest intersection point.
     * \param[in] tol Tolerance to use for the intersection.
     * \return True if the current box intersects the given ray, false otherwise.
     */
    bool intersect( Ray<Real, N> const & r, Point & c, Point & f, Real tol = DUMMY_TOL ) const
    { return r.intersect(*this, c, f, tol); }

    /**
     * Compute the squared \em interior (or \em minimum) distance between the current box a given ray.
     * \param[in] r Ray to compute the squared distance to.
     * \param[in] tol Tolerance to use for the distance computation.
     * \return The squared distance between the ray and the box.
     * \note The squared distance is positive if the ray and the box do not intersect, zero otherwise.
     */
    Real squared_interior_distance(Ray<Real, N> const & r, Real tol = DUMMY_TOL) const
    { return r.squared_interior_distance(*this, tol); }

    /**
     * Compute the squared \em interior (or \em minimum) distance between the current box a given ray,
     * returning two points at the minimum distance.
     * \param[in] r Ray to compute the squared distance to.
     * \param[out] p1 First point at the minimum distance (on the current box).
     * \param[out] p2 Second point at the minimum distance (on the ray).
     * \param[in] tol Tolerance to use for the distance computation.
     * \return The squared distance between the ray and the box.
     * \note The squared distance is positive if the ray and the box do not intersect, zero otherwise.
     */
    Real squared_interior_distance(Ray<Real, N> const & r, Point & p1, Point & p2, Real tol = DUMMY_TOL) const
    { return r.squared_interior_distance(*this, p2, p1, tol); }

    /**
     * Compute the \em interior (or \em minimum) distance between the current box a given ray.
     * \param[in] r Ray to compute the distance to.
     * \param[in] tol Tolerance to use for the distance computation.
     * \return The distance between the ray and the box.
     * \note The distance is positive if the ray and the box do not intersect, zero otherwise.
     */
    Real interior_distance(Ray<Real, N> const & r, Real tol = DUMMY_TOL) const
    { return r.interior_distance(*this, tol); }

    /**
     * Compute the \em interior (or \em minimum) distance between the current box a given ray,
     * returning two points at the minimum distance.
     * \param[in] r Ray to compute the distance to.
     * \param[out] p1 First point at the minimum distance (on the current box).
     * \param[out] p2 Second point at the minimum distance (on the ray).
     * \param[in] tol Tolerance to use for the distance computation.
     * \return The distance between the ray and the box.
     * \note The distance is positive if the ray and the box do not intersect, zero otherwise.
     */
    Real interior_distance( Ray<Real, N> const & r, Point & p1, Point & p2, Real tol = DUMMY_TOL ) const
    { return r.interior_distance(*this, p2, p1, tol); }

    /**
     * Compute the squared \em exterior (or \em maximum) distance between the current box a given ray.
     * \param[in] r Ray to compute the squared distance to.
     * \param[in] tol Tolerance to use for the distance computation.
     * \return The squared distance between the ray and the box.
     * \note The squared distance is positive if one box is not contained in the other, zero otherwise.
     */
    Real squared_exterior_distance( Ray<Real, N> const & r, Real tol = DUMMY_TOL ) const
    { return r.squared_exterior_distance(*this, tol); }

    /**
     * Compute the squared \em exterior (or \em maximum) distance between the current box a given ray,
     * returning two points at the maximum distance.
     * \param[in] r Ray to compute the squared distance to.
     * \param[out] p1 First point at the maximum distance (on the current box).
     * \param[out] p2 Second point at the maximum distance (on the ray).
     * \param[in] tol Tolerance to use for the distance computation.
     * \return The squared distance between the ray and the box.
     */
    Real squared_exterior_distance( Ray<Real, N> const & r, Point & p1, Point & p2, Real tol = DUMMY_TOL ) const
    { return r.squared_exterior_distance(*this, p2, p1, tol); }

    /**
     * Compute the \em exterior (or \em maximum) distance between the current box a given ray.
     * \param[in] r Ray to compute the distance to.
     * \param[in] tol Tolerance to use for the distance computation.
     * \return The distance between the ray and the box.
     * \note The distance is positive if the ray and the box do not intersect, zero otherwise.
     */
    Real exterior_distance( Ray<Real, N> const & r, Real tol = DUMMY_TOL ) const
    { return r.exterior_distance(*this, tol); }

    /**
     * Compute the \em exterior (or \em maximum) distance between the current box a given ray,
     * and return two points at the maximum distance.
     * \param[in] r Ray to compute the distance to.
     * \param[out] p1 First point at the maximum distance (on the current box).
     * \param[out] p2 Second point at the maximum distance (on the ray).
     * \param[in] tol Tolerance to use for the distance computation.
     * \return The distance between the ray and the box.
     * \note The distance is positive if the ray and the box do not intersect, zero otherwise.
     */
    Real exterior_distance( Ray<Real, N> const & r, Point & p1, Point & p2, Real tol = DUMMY_TOL ) const
    { return r.exterior_distance(*this, p2, p1, tol); }

    /**
     * Print the box info to an output stream.
     * \param[in] os Output stream to print the box info to.
     */
    void print(std::ostream & os) const {
      os <<
        "────────────────────────────────────────────────────────────────────────────────\n" <<
        "BOX INFO\n" <<
        "\tmin = " << this->m_min.transpose() << '\n' <<
        "\tmax = " << this->m_max.transpose() << '\n' <<
        "────────────────────────────────────────────────────────────────────────────────\n";
    }

  }; // class Box

  /**
  * Print the box info to an output stream.
  * \param[in] os Output stream to print the box info to.
  * \param[in] b Box to print.
  * \tparam Real Type of the scalar coefficients.
  * \tparam N Dimension of the ambient space.
  */
  template <typename Real, Integer N>
  std::ostream & operator<<(std::ostream & os, Box<Real, N> const & b) { b.print(os); return os; }


} // namespace AABBtree

#endif // AABBTREE_BOX_HXX
