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
  template <typename Real, Integer N, typename DerivedTree>
  class Tree
  {
    static_assert(N > 0, "Tree dimension must be positive.");
    static_assert(std::is_floating_point<Real>::value, "Tree real type must be a floating-point type.");

    using Box = AlignedBox<Real, N>; /**< Axis-aligned bounding box in N-dimensional space. */
    using BoxUPtr = std::unique_ptr<Box>; /**< Unique pointer to an axis-aligned bounding box. */
    using BoxUPtrVec = std::vector<BoxUPtr>; /**< Vector of unique pointers to an axis-aligned bounding box. */
    using Point = typename Box::Point; /**< Point in the ambient space (Eigen column vector of real numbers). */
    using Set = std::set<Integer>; /**< Set of indexes. */

    // Tree parameters
    Integer m_max_nodal_objects{16}; /**< Maximum number of objects per node. */
    Real m_long_edge_ratio{0.8}; /**< Long edge ratio for bounding boxes. */
    Real m_collision_tolerance{0.1}; /**< Overlap tolerance for bounding boxes. */
    Real m_min_size_tolerance{0.0}; /**< Minimum size tolerance for bounding boxes. */

    Tree(Tree const &) = delete; /**< Copy constructor. */
    Tree & operator=(Tree const &) = delete; /**< Copy assignment operator. */

  public:
    /**
    * Class destructor for the tree.
    */
    ~Tree() {};

    /**
    * Class constructor for the tree.
    */
    Tree() {};

    /**
    * Set the maximum number of objects per node.
    * \param[in] n Maximum number of objects per node.
    */
    void max_nodes_objects(Integer n)
    {
      AABBTREE_ASSERT(n > static_cast<Integer>(0),
        "NonRecursive::max_nodes_objects(...): input must be a positive integer.");
      this->m_max_nodal_objects = n;
    }

    /**
    * Get the maximum number of objects per node.
    * \return The maximum number of objects per node.
    */
    Integer max_nodes_objects() const {return this->m_max_nodal_objects;}

    /**
    * Set the long edge ratio for bounding boxes.
    * \param[in] ratio Long edge ratio for bounding boxes.
    */
    void long_edge_ratio(Real ratio)
    {
      AABBTREE_ASSERT(ratio > static_cast<Real>(0.0) && ratio < static_cast<Real>(1),
        "NonRecursive::long_edge_ratio(...): input must be in the range [0, 1].");
      this->m_long_edge_ratio = ratio;
    }

    /**
    * Get the long edge ratio for bounding boxes.
    * \return The long edge ratio for bounding boxes.
    */
    Real long_edge_ratio() const {return this->m_long_edge_ratio;}

    /**
    * Set the overlap tolerance for bounding boxes.
    * \param[in] tolerance Overlap tolerance for bounding boxes.
    */
    void collision_tolerance(Real tolerance)
    {
      AABBTREE_ASSERT(tolerance > static_cast<Real>(0.0) && tolerance < static_cast<Real>(1),
        "NonRecursive::collision_tolerance(...): input must be in the range [0, 1].");
      this->m_collision_tolerance = tolerance;
    }

    /**
    * Get the overlap tolerance for bounding boxes.
    * \return The overlap tolerance for bounding boxes.
    */
    Real collision_tolerance() const {return this->m_collision_tolerance;}

    /**
    * Set the minimum size tolerance for bounding boxes.
    * \param[in] tolerance Minimum size tolerance for bounding boxes.
    */
    void min_size_tolerance(Real tolerance)
    {
      AABBTREE_ASSERT(tolerance >= static_cast<Real>(0.0),
        "NonRecursive::min_size_tolerance(...): input must be a non-negative real number.");
      this->m_min_size_tolerance = tolerance;
    }

    /**
    * Get the minimum size tolerance for bounding boxes.
    * \return The minimum size tolerance for bounding boxes.
    */
    Real min_size_tolerance() const {return this->m_min_size_tolerance;}

    /**
    * Build the tree given the bounding boxes.
    * \param[in] boxes Bounding boxes to build the tree from.
    * \return True if the tree was built successfully, false otherwise.
    */
    void build(BoxPtrVec const & boxes)
    {
      return static_cast<DerivedTree *>(this)->build_impl(boxes);
    }

    /**
    * Clear the tree.
    */
    void clear()
    {
      static_cast<DerivedTree *>(this)->clear_impl();
    }

    /**
    * Check if tree is empty.
    * \return True if the tree is empty, false otherwise.
    */
    bool is_empty() const
    {
      return static_cast<DerivedTree const *>(this)->is_empty_impl();
    }

    /**
    * Print the tree internal structure in an output stream.
    * \param[in] os Output stream to print the tree to.
    */
    void print(std::ostream & os) const
    {
      static_cast<DerivedTree const *>(this)->print_impl(os);
    }

    /**
    * Build the tree given the bounding boxes.
    * \param[in] boxes Bounding boxes to build the tree from.
    * \return True if the tree was built successfully, false otherwise.
    */
    void build(std::vector<Box> const & boxes)
    {
      return static_cast<DerivedTree *>(this)->build_impl(boxes);
    }

    /**
    * Intersect the tree with a point.
    * \param[in] point Point to intersect with.
    * \param[out] candidates Intersection result (bounding box indexes).
    * \return True if the point intersects the tree, false otherwise.
    */
    bool intersect(Point const & point, Set & candidates) const
    {
      return static_cast<DerivedTree const *>(this)->intersect_impl(point, candidates);
    }

    /**
    * Intersect the tree with an axis-aligned box.
    * \param[in] box Axis-aligned box to intersect with.
    * \param[out] candidates Intersection result (bounding box indexes).
    * \return True if the point intersects the tree, false otherwise.
    */
    bool intersect(Box const & box, Set & candidates) const
    {
      return static_cast<DerivedTree const *>(this)->intersect_impl(box, candidates);
    }

    /**
    * Intersect the tree with another tree.
    * \param[in] tree Tree to intersect with.
    * \param[out] candidates Intersection result (bounding box indexes).
    * \return True if the point intersects the tree, false otherwise.
    */
    bool intersect(DerivedTree const & tree, Set & candidates) const
    {
      return static_cast<DerivedTree const *>(this)->intersect_impl(tree, candidates);
    }

    /**
    * Minimum distance between a point and the tree.
    * \param[in] point Point to compute the minimum distance to.
    * \param[out] candidates Minimum distance candidates.
    * \return The minimum distance between the point and the tree.
    */
    Real min_distance(Point const & point, Set & candidates) const
    {
      return static_cast<DerivedTree const *>(this)->min_distance_impl(point, candidates);
    }

    /**
    * Minimum distance between an axis-aligned box and the tree.
    * \param[in] box Axis-aligned box to compute the minimum distance to.
    * \param[out] candidates Minimum distance candidates.
    * \return The minimum distance between the box and the tree.
    */
    Real min_distance(Box const & box, Set & candidates) const
    {
      return static_cast<DerivedTree const *>(this)->min_distance_impl(box, candidates);
    }

    /**
    * Minimum distance between an current tree and another tree.
    * \param[in] tree Tree to compute the minimum distance to.
    * \param[out] candidates Minimum distance candidates.
    * \return The minimum distance between the trees.
    */
    Real min_distance(DerivedTree const & tree, Set & candidates) const
    {
      return static_cast<DerivedTree const *>(this)->min_distance_impl(tree, candidates);
    }

    /**
    * Maximum distance between a point and the tree.
    * \param[in] point Point to compute the maximum distance to.
    * \param[out] candidates Maximum distance candidates.
    * \return The maximum distance between the point and the tree.
    */
    Real max_distance(Point const & point, Set & candidates) const
    {
      return static_cast<DerivedTree const *>(this)->max_distance_impl(point, candidates);
    }

    /**
    * Maximum distance between an axis-aligned box and the tree.
    * \param[in] box Axis-aligned box to compute the maximum distance to.
    * \param[out] candidates Maximum distance candidates.
    * \return The maximum distance between the box and the tree.
    */
    Real max_distance(Box const & box, Set & candidates) const
    {
      return static_cast<DerivedTree const *>(this)->max_distance_impl(box, candidates);
    }

    /**
    * Maximum distance between an current tree and another tree.
    * \param[in] tree Tree to compute the maximum distance to.
    * \param[out] candidates Maximum distance candidates.
    * \return The maximum distance between the trees.
    */
    Real max_distance(DerivedTree const & tree, Set & candidates) const
    {
      return static_cast<DerivedTree const *>(this)->max_distance_impl(tree, candidates);
    }

  }; // Tree

} // namespace AABBtree

#endif // INCLUDE_AABBTREE_TREE_HXX
