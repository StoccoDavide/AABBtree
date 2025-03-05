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

#ifndef INCLUDE_AABBTREE_TREE_HXX
#define INCLUDE_AABBTREE_TREE_HXX

namespace AABBtree {

  /**
  * \brief A class representing a \em non-recursive axis-aligned bounding box tree (AABBtree).
  *
  * The Tree class provides an efficient way to store and query a set of axis-aligned bounding
  * boxes. It supports various operations such as building the tree, adding axis-aligned bounding
  * boxes, and performing intersection queries. The Tree class is templated on the type of the
  * real numbers and the dimensions of the bounding boxes and is based. The Tree class is the base
  * class for the Recursive and NonRecursive classes, which are implemented through the curiously
  * recurring template pattern (CRTP).
  * \tparam Real Type of the scalar coefficients
  * \tparam N Dimension of the ambient space.
  */
  template <Integer N, typename Real, typename DerivedTree>
  class Tree {
    static_assert( is_floating_point<Real>::value, "Tree Real type must be a floating-point type." );
    static_assert( is_integral<Integer>::value,    "Tree dimension type must be an integer type."  );
    static_assert( N > 0, "Tree dimension must be positive." );

    DerivedTree       * THIS()       { return static_cast<DerivedTree *>(this); }
    DerivedTree const * THIS() const { return static_cast<DerivedTree const *>(this); }

  protected:
    // Tree parameters
    Integer m_max_nodal_objects{16};    /**< Maximum number of objects per node. */
    Real    m_long_edge_ratio{0.8};     /**< Long edge ratio for bounding boxes. */
    Real    m_collision_tolerance{0.1}; /**< Overlap tolerance for bounding boxes. */
    Real    m_min_size_tolerance{0.0};  /**< Minimum size tolerance for bounding boxes. */

  public:

    using Box        = Box<N,Real>;         /**< Axis-aligned bounding box in N-dimensional space. */
    using BoxUPtr    = BoxUPtr<N,Real>;     /**< Unique pointer to an axis-aligned bounding box. */
    using BoxUPtrVec = BoxUPtrVec<N,Real>;  /**< Vector of unique pointers to an axis-aligned bounding box. */

    using Point      = Point<Real,N>;       /**< Point in the ambient space (Eigen column vector of real numbers). */

    Tree ( Tree const & ) = delete;          /**< Copy constructor. */
    Tree & operator=(Tree const &) = delete; /**< Copy assignment operator. */

    /*!
     * Class destructor for the tree.
     */
    ~Tree() = default;

    /*
    //   ___        _               _
    //  / __|___ __| |_ _ _ _  _ __| |_ ___ _ _ ___
    // | (__/ _ (_-<  _| '_| || / _|  _/ _ \ '_(_-<
    //  \___\___/__/\__|_|  \_,_\__|\__\___/_| /__/
    */

    /*!
     * Class constructor for the tree.
     */
    Tree() = default;

    /*
    //           _
    //  ___  ___| |_
    // / __|/ _ \ __|
    // \__ \  __/ |_
    // |___/\___|\__|
    */

    /**
     * Set the maximum number of objects per node.
     * \param[in] n Maximum number of objects per node.
     */
    void
    set_max_nodes_objects( Integer const n ) {
      AABBTREE_ASSERT( n > 0, "NonRecursive::max_nodes_objects(...): input must be a positive integer." );
      m_max_nodal_objects = n;
    }

    /**
     * Set the long edge ratio for bounding boxes.
     * \param[in] ratio Long edge ratio for bounding boxes.
     */
    void
    set_long_edge_ratio( Real const ratio ) {
      AABBTREE_ASSERT(
        ratio > 0 && ratio < static_cast<Real>(1),
        "NonRecursive::long_edge_ratio(...): input must be in the range [0, 1]."
      );
      m_long_edge_ratio = ratio;
    }

    /**
     * Set the overlap tolerance for bounding boxes.
     * \param[in] tolerance Overlap tolerance for bounding boxes.
     */
    void
    set_collision_tolerance( Real const tolerance ) {
      AABBTREE_ASSERT(
        tolerance > 0 && tolerance < static_cast<Real>(1),
        "NonRecursive::collision_tolerance(...): input must be in the range [0, 1]."
      );
      m_collision_tolerance = tolerance;
    }

    /**
     * Set the minimum size tolerance for bounding boxes.
     * \param[in] tolerance Minimum size tolerance for bounding boxes.
     */
    void
    set_min_size_tolerance( Real const tolerance ) {
      AABBTREE_ASSERT(
        tolerance >= 0,
        "NonRecursive::min_size_tolerance(...): input must be a non-negative real number."
      );
      m_min_size_tolerance = tolerance;
    }

    /*
    //  _        __
    // (_)_ __  / _| ___
    // | | '_ \| |_ / _ \
    // | | | | |  _| (_) |
    // |_|_| |_|_|  \___/
    */

    /**
     * Get the maximum number of objects per node.
     * \return The maximum number of objects per node.
     */
    Integer max_nodes_objects() const { return this->m_max_nodal_objects; }

    /**
     * Get the long edge ratio for bounding boxes.
     * \return The long edge ratio for bounding boxes.
     */
    Real long_edge_ratio() const { return m_long_edge_ratio; }

    /**
     * Get the overlap tolerance for bounding boxes.
     * \return The overlap tolerance for bounding boxes.
     */
    Real collision_tolerance() const { return m_collision_tolerance; }

    /**
     * Get the minimum size tolerance for bounding boxes.
     * \return The minimum size tolerance for bounding boxes.
     */
    Real min_size_tolerance() const { return m_min_size_tolerance; }

    /*
    //  ____  _   _       _        _      _
    // | __ )| \ | |     | |_ _ __(_) ___| | __
    // |  _ \|  \| |_____| __| '__| |/ __| |/ /
    // | |_) | |\  |_____| |_| |  | | (__|   <
    // |____/|_| \_|      \__|_|  |_|\___|_|\_\
    */

    /**
     * Build the tree given the bounding boxes.
     * \param[in] boxes Bounding boxes to build the tree from.
     * \return True if the tree was built successfully, false otherwise.
     */
    void build( BoxUPtrVec const & boxes ) { return THIS()->build_impl(boxes); }

    /**
     * Clear the tree.
     */
    void clear() { THIS()->clear_impl(); }

    /**
     * Check if tree is empty.
     * \return True if the tree is empty, false otherwise.
     */
    bool is_empty() const { return THIS()->is_empty_impl(); }

    /**
     * Print the tree internal structure in an output stream.
     * \param[in] os Output stream to print the tree to.
     */
    void print( ostream_type & os ) const { THIS()->print_impl(os); }

    /*
    //  _       _                          _
    // (_)_ __ | |_ ___ _ __ ___  ___  ___| |_
    // | | '_ \| __/ _ \ '__/ __|/ _ \/ __| __|
    // | | | | | ||  __/ |  \__ \  __/ (__| |_
    // |_|_| |_|\__\___|_|  |___/\___|\___|\__|
    */
    
    /**
     * Intersect the tree with a point.
     * \param[in] point Point to intersect with.
     * \param[out] candidates Intersection result (bounding box indexes).
     * \return True if the point intersects the tree, false otherwise.
     */
    bool intersect( Point const & point, Set & candidates ) const { return THIS()->intersect_impl(point, candidates); }

    /**
     * Intersect the tree with an axis-aligned box.
     * \param[in] box Axis-aligned box to intersect with.
     * \param[out] candidates Intersection result (bounding box indexes).
     * \return True if the point intersects the tree, false otherwise.
     */
    bool intersect( Box const & box, Set & candidates ) const { return THIS()->intersect_impl(box, candidates); }

    /**
     * Intersect the tree with another tree.
     * \param[in] tree Tree to intersect with.
     * \param[out] candidates Intersection result (bounding box indexes).
     * \return True if the point intersects the tree, false otherwise.
     */
    bool intersect( DerivedTree const & tree, Set & candidates ) const { return THIS()->intersect_impl(tree, candidates); }

    /*
    //      _ _     _
    //   __| (_)___| |_ __ _ _ __   ___ ___
    //  / _` | / __| __/ _` | '_ \ / __/ _ \
    // | (_| | \__ \ || (_| | | | | (_|  __/
    //  \__,_|_|___/\__\__,_|_| |_|\___\___|
    */

    /**
     * Minimum distance between a point and the tree.
     * \param[in] point Point to compute the minimum distance to.
     * \param[out] candidates Minimum distance candidates.
     * \return The minimum distance between the point and the tree.
     */
    Real min_distance( Point const & point, Set & candidates ) const { return THIS()->min_distance_impl(point, candidates); }

    /**
    * Minimum distance between an axis-aligned box and the tree.
    * \param[in] box Axis-aligned box to compute the minimum distance to.
    * \param[out] candidates Minimum distance candidates.
    * \return The minimum distance between the box and the tree.
    */
    Real min_distance( Box const & box, Set & candidates ) const { return THIS()->min_distance_impl(box, candidates); }

    /**
     * Minimum distance between an current tree and another tree.
     * \param[in] tree Tree to compute the minimum distance to.
     * \param[out] candidates Minimum distance candidates.
     * \return The minimum distance between the trees.
     */
    Real min_distance( DerivedTree const & tree, Set & candidates ) const { return THIS()->min_distance_impl(tree, candidates); }

    /**
     * Maximum distance between a point and the tree.
     * \param[in] point Point to compute the maximum distance to.
     * \param[out] candidates Maximum distance candidates.
     * \return The maximum distance between the point and the tree.
     */
    Real max_distance( Point const & point, Set & candidates ) const { return THIS()->max_distance_impl(point, candidates); }

    /**
     * Maximum distance between an axis-aligned box and the tree.
     * \param[in] box Axis-aligned box to compute the maximum distance to.
     * \param[out] candidates Maximum distance candidates.
     * \return The maximum distance between the box and the tree.
     */
    Real max_distance( Box const & box, Set & candidates ) const { return THIS()->max_distance_impl(box, candidates); }

    /**
     * Maximum distance between an current tree and another tree.
     * \param[in] tree Tree to compute the maximum distance to.
     * \param[out] candidates Maximum distance candidates.
     * \return The maximum distance between the trees.
     */
    Real max_distance( DerivedTree const & tree, Set & candidates ) const { return THIS()->max_distance_impl(tree, candidates); }

  }; // Tree

} // namespace AABBtree

#endif // INCLUDE_AABBTREE_TREE_HXX
