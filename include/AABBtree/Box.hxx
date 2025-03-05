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

#ifndef AABBTREE_BOX_HXX
#define AABBTREE_BOX_HXX

namespace AABBtree {

  /**
  * \brief Class container for an axis-aligned box.
  *
  * This class represents an axis-aligned box as a pair of opposite corners, i.e., a minimal corner
  * and a maximal corner.
  * \tparam Real Type of the scalar coefficients
  * \tparam N Dimension of the ambient space.
  */
  template <typename Real, Integer N>
  class Box {
    static_assert( is_floating_point<Real>::value, "Box Real type must be a floating-point type." );
    static_assert( is_integral<Integer>::value, "Integer must be an integer type." );
    static_assert( N > 0, "Box dimension must be positive." );

    constexpr static Real EPS{numeric_limits<Real>::epsilon()};      /**> Machine epsilon for the scalar type. */
    constexpr static Real MAX{static_cast<Real>(1.0)/EPS};           /**> Maximum value for the scalar type.   */
    constexpr static Real MIN{-MAX};                                 /**> Minimum value for the scalar type.   */
    constexpr static Real DUMMY_TOL{EPS*static_cast<Real>(100.0)};   /**> Dummy tolerance for the scalar type. */

    using Point      = Point<Real,N>;
    using Vector     = Vector<Real,N>;
    using BoxUPtr    = BoxUPtr<Real,N>;
    using BoxUPtrVec = BoxUPtrVec<Real,N>;

    Point m_min; /**< Minimal corner of the box. */
    Point m_max; /**< Maximal corner of the box. */

  public:

    /*
    //   ___        _               _
    //  / __|___ __| |_ _ _ _  _ __| |_ ___ _ _ ___
    // | (__/ _ (_-<  _| '_| || / _|  _/ _ \ '_(_-<
    //  \___\___/__/\__|_|  \_,_\__|\__\___/_| /__/
    */

    /**
     * Class constructor for an empty axis-aligned box.
     */
    Box() {
      m_min.fill(MAX);
      m_max.fill(MIN); 
      this->set_empty();
    }

    /**
     * Class destructor for the axis-aligned box.
     */
    ~Box() = default;

    /**
     * Class constructor for a axis-aligned box given the minimal and maximal corners.
     * \param[in] t_min Minimal corner of the box.
     * \param[in] t_max Maximal corner of the box.
     * \warning If the reordering is not performed, and the minimal corner is greater than the maximal
     * corner in any dimension, undefined behavior may occur.
     */
    Box( Point const & t_min, const Point & t_max ) : m_min(t_min), m_max(t_max) {}

    /**
     * Class constructor for a axis-aligned box containing a single point.
     * \param[in] p Point to be contained in the box.
     */
    Box( Point const & p ) : m_min(p), m_max(p) {}

    /**
     * Class constructor for the 1D axis-aligned bounding box.
     * \param[in] x_min Minimal \f$ x \f$-axis corner of the box.
     * \param[in] x_max Maximal \f$ x \f$-axis corner of the box.
     * \tparam T Type of the scalar coefficients.
     * \note This constructor is only available for 1D boxes.
     */
    template <typename T = Real>
    Box<T, 1>( Real const x_min, Real const x_max ) : m_min(x_min), m_max(x_max) {}

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
    Box<T, 2>(
      Real const x_min,
      Real const y_min,
      Real const x_max,
      Real const y_max
    ) : m_min(x_min, y_min), m_max(x_max, y_max)
    {}

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
    Box<T, 3>(
      Real const x_min,
      Real const y_min,
      Real const z_min,
      Real const x_max,
      Real const y_max,
      Real const z_max
    ) : m_min(x_min, y_min, z_min), m_max(x_max, y_max, z_max)
    {}

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
    Box<T, 4>(
      Real const x_min,
      Real const y_min,
      Real const z_min,
      Real const w_min,
      Real const x_max,
      Real const y_max,
      Real const z_max,
      Real const w_max
    ) : m_min(x_min, y_min, z_min, w_min), m_max(x_max, y_max, z_max, w_max)
    {}

    /**
     * Reorder the corners of the box such that the minimal corner is less than the maximal corner.
     */
    void
    reorder() {
      for ( Integer i{0}; i < N; ++i )
        if ( m_min[i] > m_max[i] )
          swap( m_min[i], m_max[i] );
    }

    /*
    //   ___
    //  / __|___ _ __ _  _
    // | (__/ _ \ '_ \ || |
    //  \___\___/ .__/\_, |
    //          |_|   |__/
    */

    /**
     * Copy constructor for a axis-aligned box given another box with a different scalar type.
     * \param[in] b Box to copy.
     * \tparam OtherReal Type of the scalar coefficients of the box to copy.
     */
    template<typename OtherReal, typename OtherInteger>
    explicit
    Box( Box<OtherReal,N> const & b )
    : m_min(b.m_min.template cast<Real>())
    , m_max(b.m_max.template cast<Real>())
    {}

    /*
    //   ___         _
    //  / __|__ _ __| |_
    // | (__/ _` (_-<  _|
    //  \___\__,_/__/\__|
    */

    /**
     * Cast the current object to a new scalar type.
     * \tparam NewReal the new scalar type
     * \note If the new real typeis equal to the current scalar type currently used, then this
     * function returns a const reference to the current object.
     */
    template<typename NewReal>
    Box<NewReal,N>
    cast() const {
      if constexpr (is_same<Real, NewReal>::value ) return *this;
      return Box<NewReal,N>( m_min.template cast<NewReal>(), m_max.template cast<NewReal>());
    }

    /*
    //  _        __
    // (_)_ __  / _| ___
    // | | '_ \| |_ / _ \
    // | | | | |  _| (_) |
    // |_|_| |_|_|  \___/
    */

    /**
     * Check if the box is empty, i.e., the minimal corner is greater than the maximal corner in any
     * dimension.
     * \return True if the box is empty, false otherwise.
     */
    bool is_empty() const { return !( m_max.array() > m_min.array() ).all(); }

    /**
     * Check if the current box is approximately equal to another box.
     * \param[in] b Box to compare with.
     * \param[in] tol Tolerance to use for the comparison.
     */
    bool
    is_approx( Box const & b, Real const tol = DUMMY_TOL ) const {
      return m_min.isApprox(b.m_min, tol) && 
             m_max.isApprox(b.m_max, tol);
    }

    /**
     * Check if the current box is degenrate, i.e., it has zero volume.
     * \param[in] tol Tolerance to use for the comparison.
     * \return True if the box is degenerate, false otherwise.
     */
    bool is_degenerate( Real const tol = DUMMY_TOL ) const { return this->m_min.isApprox( m_max, tol ); }

    /**
     * Get the minimal corner.
     * \return The minimal corner of the box.
     */
    Point const & min() const { return m_min; }

    /**
     * Get a non const reference to the minimal corner.
     * \return The minimal corner of the box.
     */
    Point & min() { return m_min; }

    /**
     * Get the /f$ i /f$-th component of the minimal corner.
     * \param[in] i Index of the component.
     * \return The /f$ i /f$-th component of the minimal corner.
     */
    Real min( Integer i ) const { return m_min[i]; }

    /**
     * Get the maximal corner.
     * \return The maximal corner of the box.
     */
    Point const & max() const { return m_max; }

    /**
     * Get a non const reference to the maximal corner.
     * \return The maximal corner of the box.
     */
    Point & max() { return m_max; }

    /**
     * Get the /f$ i /f$-th component of the maximal corner.
     * \param[in] i Index of the component.
     * \return The /f$ i /f$-th component of the maximal corner.
     */
    Real max( Integer i ) const { return m_max[i]; }

    /**
     */
    Real length( Integer i ) const { return m_max(i)-m_min(i); }

    /**
     * Get the longest axis of the box.
     * \return The longest axis of the box.
     */
    Integer
    longest_axis( Real & mx, Real & mp ) const {
      Integer ipos{-1};
      mx = -1;
      for ( Integer i{0}; i < N; ++i )
        if ( Real const mx1{ m_max(i) - m_min(i) }; mx < mx1 )
          { mx = mx1; mp = (m_max(i) + m_min(i) ) / 2; ipos = i; }
      return ipos;
    }

    /**
     * Get the longest axis of the box.
     * \return The longest axis of the box.
     */
    void
    sort_axis( Vector & S, Integer ipos[N] ) const {
      std::iota( ipos, ipos + N, 0 );
      S = m_max - m_min;
      std::sort( ipos, ipos+N, [&S]( Integer a, Integer b ) -> bool { return S[a] > S[b]; } );
    }

    /**
     * Compute the center of the box.
     * \return The center of the box.
     */
    Point center() const { return (m_min + m_max)/2; }

    /**
     * Compute the center of the box.
     * \return The center of the box.
     */
    Real center( Integer i ) const { return (m_min(i) + m_max(i))/2; }

    /**
     * Compute the lengths of the bounding box's sides.
     * \return The lengths of the bounding box's sides.
     */
    Vector sizes() const { return m_max - m_min; }

    /**
     * Compute the volume of the bounding box.
     * \return The volume of the bounding box.
     */

    Real volume() const { return this->sizes().prod(); }

    /**
     * Compute the external surface area of the bounding box.
     * \return The external surface area of the bounding box.
     */
    Real
    surface() const {
      Vector sizes { this->sizes() };
      Real   area  { static_cast<Real>(0) };
      for ( Integer i{0}; i < N; ++i ) {
        Real prod{ static_cast<Real>(1) };
        for ( Integer j{0}; j < N; ++j ) if (j != i) prod *= sizes[j];
        area += static_cast<Real>(2.0)*prod;
      }
      return area;
    }

    /**
     * Compute the bounding box diagonal vector.
     * \note If the length of the diagonal is needed, \c diagonal().norm() will provide it.
     * \return The diagonal vector of the bounding box.
     */
    Vector diagonal() const { return this->sizes(); }

    /*
    //                      _ _  __
    //  _ __ ___   ___   __| (_)/ _|_   _
    // | '_ ` _ \ / _ \ / _` | | |_| | | |
    // | | | | | | (_) | (_| | |  _| |_| |
    // |_| |_| |_|\___/ \__,_|_|_|  \__, |
    //                              |___/
    */

    /**
     * Set the box to be degenerate.
     * \param[in] p Point to set the box to.
     */
    void
    set_degenerate( Point const & p ) {
      m_min = p;
      m_max = p;
    }

    /**
     * Set the box to be empty.
     */
    void
    set_empty() {
      m_min.setConstant(MAX);
      m_max.setConstant(MIN);
    }

    /**
     * Clamp the current box by a given box .
     * \param[in] b Box to extend the box to.
     * \return A reference to the current box.
     * \note If the boxes do not intersect, the resulting box is empty.
     */
    Box &
    clamp( Box const & b ) {
      m_min = m_min.cwiseMax(b.m_min);
      m_max = m_max.cwiseMin(b.m_max);
      return *this;
    }

    /**
     * Return the intersection of the current box and another box.
     * \param[in] b Box to intersect with.
     * \return The intersection of the current box and the given box.
     * \note If the boxes do not intersect, the resulting box is empty.
     */
    Box
    intersection( Box const & b ) const {
      return Box( m_min.cwiseMax(b.m_min), m_max.cwiseMin(b.m_max) );
    }

    /**
     * Return the union of the current box and another box.
     * \param[in] b Box to merge with.
     * \return The union of the current box and the given box.
     */
    Box
    get_merged( Box const & b ) const {
      return Box( m_min.cwiseMin(b.m_min), m_max.cwiseMax(b.m_max) );
    }

    /**
     * Translate the current box by a given vector.
     * \param[in] t Vector to translate the box by.
     * \return A reference to the current box.
     */
    Box &
    translate( Vector const & t ) {
      m_min += t;
      m_max += t;
      return *this;
    }

    /**
     * Translate the current box by a given vector.
     * \param[in] t Vector to translate the box by.
     * \return A copy of the current translated box.
     */
    Box
    get_translated( Vector const & t ) const {
      Box result( m_min, m_max );
      result.translate(t);
      return result;
    }

    /*
    //       _               _
    //   ___| |__   ___  ___| | __
    //  / __| '_ \ / _ \/ __| |/ /
    // | (__| | | |  __/ (__|   <
    //  \___|_| |_|\___|\___|_|\_\
    */

    /**
     * Check if the point is inside the box.
     * \param[in] p Point to check.
     * \return True if the point is inside the box, false otherwise.
     */
    bool
    contains( Point const  & p ) const {
      return ( m_min.array() <= p.array()     ).all() &&
             ( p.array()     <= m_max.array() ).all();
    }

    /**
     * Check if the point is inside the interval, on the left or right.
     * \param[in] p Point to check.
     * \return True if the point is inside the box, false otherwise.
     */
    Integer
    classify( Real const x, Real const x_tol, Integer const idim ) const {
      if ( x < m_min(idim)+x_tol ) return -1;
      if ( x > m_max(idim)-x_tol ) return +1;
      return 0;
    }

    /**
     * Check if the current box contains a given box.
     * \param[in] b Box to check.
     * \return True if the current box contains the given box, false otherwise.
     */
    bool
    contains( Box const & b ) const {
      return ( m_min.array()   <= b.m_min.array() ).all() &&
             ( b.m_max.array() <= m_max.array()   ).all();
    }

    /**
     * Check if the current box intersects a given box.
     * \param[in] b Box to check.
     * \return True if the current box intersects the given box, false otherwise.
     */
    bool
    intersects( Box const & b ) const {
      return ( m_min.array()   <= b.m_max.array() ).all() &&
             ( b.m_min.array() <= m_max.array()   ).all();
    }

    /*
    //            _                 _
    //   _____  _| |_ ___ _ __   __| |
    //  / _ \ \/ / __/ _ \ '_ \ / _` |
    // |  __/>  <| ||  __/ | | | (_| |
    //  \___/_/\_\\__\___|_| |_|\__,_|
    */

    /**
     * Extend the current box such that it contains a given point.
     * \param[in] p Point to extend the box to.
     * \return A reference to the current box.
     */
    template<typename Derived>
    Box &
    extend( Point const & p ) {
      m_min = m_min.cwiseMin(p);
      m_max = m_max.cwiseMax(p);
      return *this;
    }

    /**
     * Extend the current box such that it contains a given box.
     * \param[in] b Box to extend the box to.
     * \return A reference to the current box.
     * \note Merging with an empty box may result in a box bigger than \c *this.
     */
    Box &
    extend( Box const & b ) {
      m_min = m_min.cwiseMin(b.m_min);
      m_max = m_max.cwiseMax(b.m_max);
      return *this;
    }

    /**
     * Extend the current box such that it contains a given box.
     * \param[in] b Box to extend the box to.
     * \return A reference to the current box.
     * \note Merging with an empty box may result in a box bigger than \c *this.
     */
    Box &
    extend( BoxUPtr const & b ) {
      m_min = m_min.cwiseMin(b->m_min);
      m_max = m_max.cwiseMax(b->m_max);
      return *this;
    }

    /**
     * Extend the current box such that it contains the given boxes.
     * \param[in] b Boxes to extend the box to.
     * \return A reference to the current box.
     * \note Merging with an empty box may result in a box bigger than \c *this.
     */
    Box &
    extend( BoxUPtrVec const & bvec ) {
      for ( auto const & box : bvec ) this->extend( box );
      return *this;
    }

    /*
    //      _ _     _
    //   __| (_)___| |_ __ _ _ __   ___ ___
    //  / _` | / __| __/ _` | '_ \ / __/ _ \
    // | (_| | \__ \ || (_| | | | | (_|  __/
    //  \__,_|_|___/\__\__,_|_| |_|\___\___|
    */

    /**
     * Compute the squared distance between the current box a given point.
     * \param[in] p Point to compute the squared distance to.
     * \return The squared distance between the box and the point.
     * \note The returned value is positive if the point is outside the box, zero otherwise.
     */
    Real
    squared_distance( Point const & p ) const {
      if ( this->contains(p) ) return 0;
      Real dist2{0};
      for ( Integer i{0}; i < N; ++i ) {
        if      ( m_min[i] > p[i]     ) { Real aux{ m_min[i] - p[i] }; dist2 += aux*aux; }
        else if ( p[i]     > m_max[i] ) { Real aux{ p[i] - m_max[i] }; dist2 += aux*aux; }
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
    Real
    squared_distance( Point const & p, Point & p_min, Point & p_max ) const {
      if ( this->contains(p) ) return static_cast<Real>(0.0);
      p_min = p.cwiseMin(m_min);
      p_max = p.cwiseMax(m_max);
      return (p_max - p_min).squaredNorm();
    }

    /**
     * Compute the distance between the current box a given point.
     * \param[in] p Point to compute the distance to.
     * \return The distance between the box and the point.
     * \note The distance is positive if the point is outside the box, zero otherwise.
     */
    Real distance( Point const & p ) const { return sqrt( squared_distance(p) ); }

    /**
     * Compute the distance between the current box a given point, returning two points at the given
     * distance.
     * \param[in] p Point to compute the distance to.
     * \param[out] p_max First point at the distance.
     * \param[out] p_min Second point at the distance.
     * \return The distance between the box and the point.
     * \note The returned value is positive if the point is outside the box, zero otherwise.
     */
    Real
    distance( Point const & p, Point & p_min, Point & p_max) const
    { return sqrt( squared_distance(p, p_min, p_max) ); }

    /**
     * Compute the squared \em interior (or \em minimum) distance between the current box a given box.
     * \param[in] b Box to compute the squared distance to.
     * \return The squared distance between the boxes.
     * \note The squared distance is positive if the boxes do not intersect, zero otherwise.
     */
    Real
    squared_interior_distance( Box const & b ) const {
      if ( this->intersects(b) ) return 0;
      Real dist2{0.0};
      for ( Integer i{0}; i < N; ++i ) {
        if      ( m_min[i]   > b.m_max[i] ) { Real aux{ m_min[i] - b.m_max[i] }; dist2 += aux*aux; }
        else if ( b.m_min[i] > m_max[i]   ) { Real aux{ b.m_min[i] - m_max[i] }; dist2 += aux*aux; }
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
    Real
    squared_interior_distance( Box const & b, Point & p_min, Point & p_max) const {
      if ( this->intersects(b) ) return 0;
      for ( Integer i{0}; i < N; ++i ) {
        if      ( m_min[i]   > b.m_max[i] ) { p_min[i] = b.m_max[i]; p_max[i] = m_min[i];   }
        else if ( b.m_min[i] > m_max[i]   ) { p_min[i] = m_max[i];   p_max[i] = b.m_min[i]; }
        else    { p_min[i] = p_max[i] = ( std::min(m_max[i], b.m_max[i]) + std::max(m_min[i], b.m_min[i]) ) * static_cast<Real>(0.5); }
      }
      return (p_max - p_min).squaredNorm();
    }

    /**
     * Compute the \em interior (or \em minimum) distance between the current box a given box.
     * \param[in] b Box to compute the distance to.
     * \return The distance between the boxes.
     * \note The distance is positive if the boxes do not intersect, zero otherwise.
     */
    Real interior_distance( Box const & b) const { return sqrt( squared_interior_distance(b) ); }

    /**
     * Compute the \em interior (or \em minimum) distance between the current box a given box,
     * returning two points at the minimum distance.
     * \param[in] b Box to compute the distance to.
     * \param[out] p_max First point at the minimum distance.
     * \param[out] p_min Second point at the minimum distance.
     * \return The distance between the boxes.
     * \note The distance is positive if the boxes do not intersect, zero otherwise.
     */
    Real
    interior_distance( Box const & b, Point & p_min, Point & p_max ) const {
      return sqrt(static_cast<Real>(this->squared_interior_distance(b, p_min, p_max)));
    }

    /**
     * Compute the squared \em exterior (or \em maximum) distance between the current box a given box.
     * \param[in] b Box to compute the squared distance to.
     * \return The squared distance between the boxes.
     * \note The squared distance is positive if one box is not contained in the other, zero otherwise.
     */
    Real
    squared_exterior_distance( Box const & b ) const {
      Real dist2{0};
      for ( Integer i{0}; i < N; ++i ) {
        Real aux{ std::max(m_max[i], b.m_max[i]) - std::min(m_min[i], b.m_min[i]) };
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
    Real
    squared_exterior_distance( Box const & b, Point & p_min, Point & p_max ) const {
      p_min = b.m_min.cwiseMin(m_min);
      p_max = b.m_max.cwiseMax(m_max);
      return (p_max - p_min).squaredNorm();
    }

    /**
     * Compute the \em exterior (or \em maximum) distance between the current box a given box.
     * \param[in] b Box to compute the distance to.
     * \return The distance between the boxes.
     * \note The distance is positive if the boxes do not intersect, zero otherwise.
     */
    Real exterior_distance( Box const & b ) const { return sqrt( squared_exterior_distance(b) ); }

    /**
     * Compute the \em exterior (or \em maximum) distance between the current box a given box,
     * and return two points at the maximum distance.
     * \param[in] b Box to compute the distance to.
     * \param[out] p_max First point at the maximum distance.
     * \param[out] p_min Second point at the maximum distance.
     * \return The distance between the boxes.
     * \note The distance is positive if the boxes do not intersect, zero otherwise.
     */
    Real
    exterior_distance( Box const & b, Point & p_min, Point & p_max ) const
    { return sqrt( squared_exterior_distance(b, p_min, p_max) ); }

  }; // class Box

  template class Box<float,1>;
  template class Box<float,2>;
  template class Box<float,3>;
  template class Box<float,4>;
  template class Box<float,5>;
  template class Box<float,6>;
  template class Box<float,7>;
  template class Box<float,8>;
  template class Box<float,9>;
  template class Box<float,10>;

  template class Box<double,1>;
  template class Box<double,2>;
  template class Box<double,3>;
  template class Box<double,4>;
  template class Box<double,5>;
  template class Box<double,6>;
  template class Box<double,7>;
  template class Box<double,8>;
  template class Box<double,9>;
  template class Box<double,10>;

} // namespace AABBtree

#endif // AABBTREE_Box_HXX
