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

#include "AABBtree/Box.hxx"

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
  * \tparam DerivedTree Type of the derived tree class.
  */
  template <typename Real, Integer N, typename DerivedTree>
  class Tree
  {
    static_assert(std::is_floating_point<Real>::value, "Tree Real type must be a floating-point type.");
    static_assert(std::is_integral<Integer>::value, "Tree dimension type must be an integer type." );
    static_assert(N > 0, "Tree dimension must be positive.");

    DerivedTree * THIS() {return static_cast<DerivedTree *>(this);} /**< Cast the current object to the derived tree class. */
    DerivedTree const * THIS() const {return static_cast<DerivedTree const *>(this);} /**< Cast the current object to the derived tree class. */

  protected:
    using Box = Box<Real, N>; /**< Axis-aligned bounding box in N space. */
    using BoxUniquePtr = BoxUniquePtr<Real, N>; /**< Unique pointer to an axis-aligned bounding box. */
    using BoxUniquePtrList = BoxUniquePtrList<Real, N>; /**< Vector of unique pointers to an axis-aligned bounding box. */
    using Point = Point<Real, N>; /**< Point in the ambient space (Eigen column vector of real numbers). */

    // Bounding boxes
    std::unique_ptr<BoxUniquePtrList> m_boxes{nullptr};

    // Tree parameters
    Integer m_max_nodal_objects{1}; /**< Maximum number of objects per node. */
    Real m_separation_ratio_tolerance{0.1}; /**< Tolerance for bounding boxes separation. */
    Real m_balance_ratio_tolerance{0.25}; /**< Tolerance for bounding boxes balance. */
    Real m_min_box_size{0.0}; /**< Minimum size tolerance for bounding boxes. */

    // Statistics
    mutable Integer m_check_counter{0}; /**< Number of collision check (for statistic). */
    mutable Integer m_dump_counter{0}; /**< Number of dumpings (for statistic). */

  public:
    Tree (Tree const &) = delete; /**< Copy constructor. */
    Tree & operator=(Tree const &) = delete; /**< Copy assignment operator. */

    /**
    * Class destructor for the tree.
    */
    ~Tree() = default;

    /**
    * Class constructor for the tree.
    */
    Tree() = default;

    /**
    * Get a look at the vector of unique pointers to the bounding boxes.
    * \return A const reference to the vector of unique pointers to the bounding boxes.
    */
    BoxUniquePtrList const & boxes() const {return *this->m_boxes;}

    /**
    * Get a look at the i-th unique pointer to the bounding box.
    * \param[in] i Index of the unique pointer to the bounding box.
    * \return A const reference to the i-th unique pointer to the bounding box.
    */
    BoxUniquePtr const & box(Integer const i) const {return (*this->m_boxes)[i];}

    /**
    * Set the maximum number of objects per node.
    * \param[in] n Maximum number of objects per node.
    */
    void max_nodal_objects(Integer const n)
    {
      #define CMD "AABBtree::NonRecursive::max_nodal_objects(...): "
      AABBTREE_ASSERT(n > 0, CMD "input must be a positive integer.");
      this->m_max_nodal_objects = n;
      #undef CMD
    }

    /**
    * Get the maximum number of objects per node.
    * \return The maximum number of objects per node.
    */
    Integer max_nodal_objects() const {return this->m_max_nodal_objects;}

    /**
    * Set the balance ratio tolerance for bounding boxes.
    * \param[in] ratio Balance ratio tolerance for bounding boxes.
    */
    void separation_ratio_tolerance(Real const ratio)
    {
      #define CMD "AABBtree::NonRecursive::separation_ratio_tolerance(...): "
      AABBTREE_ASSERT(ratio > 0.0 && ratio < 1.0, CMD "input must be in the range [0, 1].");
      this->m_separation_ratio_tolerance = ratio;
      #undef CMD
    }

    /**
    * Get the balance ratio tolerance for bounding boxes.
    * \return The balance ratio tolerance for bounding boxes.
    */
    Real separation_ratio_tolerance() const {return this->m_separation_ratio_tolerance;}

    /**
    * Set the minimum size for bounding boxes.
    * \param[in] size Minimum size for bounding boxes.
    */
    void min_box_size(Real const size)
    {
      #define CMD "AABBtree::NonRecursive::min_box_size(...): "
      AABBTREE_ASSERT(size >= 0.0, CMD "input must be a non-negative real number.");
      this->m_min_box_size = size;
      #undef CMD
    }

    /**
    * Get the minimum size for bounding boxes.
    * \return The minimum size for bounding boxes.
    */
    Real min_box_size() const {return this->m_min_box_size;}

    /**
    * Check if tree is empty.
    * \return True if the tree is empty, false otherwise.
    */
    bool is_empty() const {return THIS()->is_empty_impl();}

    /**
    * Clear the tree.
    */
    void clear() {THIS()->clear_impl();}

    /**
    * Print the tree internal structure in an output stream.
    * \param[in] os Output stream to print the tree to.
    */
    void print(OutStream & os) const {THIS()->print_impl(os);}

    /**
    * Build the tree given the bounding boxes.
    * \param[in] boxes Bounding boxes to build the tree from.
    */
    void build(std::unique_ptr<BoxUniquePtrList> boxes)
    {this->m_boxes = std::move(boxes); return THIS()->build_impl();}

    /**
    * Intersect the tree with a point.
    * \param[in] point Point to intersect with.
    * \param[out] candidates Intersection result (bounding box indexes).
    * \return True if the point intersects the tree, false otherwise.
    */
    bool intersect(Point const & point, IndexSet & candidates) const
    {return THIS()->intersect_impl(point, candidates);}

    /**
    * Intersect the tree with an axis-aligned box.
    * \param[in] box Axis-aligned box to intersect with.
    * \param[out] candidates Intersection result (bounding box indexes).
    * \return True if the point intersects the tree, false otherwise.
    */
    bool intersect(Box const & box, IndexSet & candidates) const
    {return THIS()->intersect_impl(box, candidates);}

    /**
    * Intersect the tree with another tree.
    * \param[in] tree Tree to intersect with.
    * \param[out] candidates Intersection result (bounding box indexes).
    * \return True if the point intersects the tree, false otherwise.
    */
    bool intersect(DerivedTree const & tree, IndexMap & candidates) const
    {return THIS()->intersect_impl(tree, candidates);}

    /**
    * Minimum distance between a point and the tree.
    * \param[in] point Point to compute the minimum distance to.
    * \param[out] candidates Minimum distance candidates.
    * \return The minimum distance between the point and the tree.
    */
    Real min_distance(Point const & point, IndexSet & candidates) const
    {return THIS()->min_distance_impl(point, candidates);}

    /**
    * Minimum distance between an axis-aligned box and the tree.
    * \param[in] box Axis-aligned box to compute the minimum distance to.
    * \param[out] candidates Minimum distance candidates.
    * \return The minimum distance between the box and the tree.
    */
    Real min_distance(Box const & box, IndexSet & candidates) const
    {return THIS()->min_distance_impl(box, candidates);}

    /**
    * Minimum distance between an current tree and another tree.
    * \param[in] tree Tree to compute the minimum distance to.
    * \param[out] candidates Minimum distance candidates.
    * \return The minimum distance between the trees.
    */
    Real min_distance(DerivedTree const & tree, IndexSet & candidates) const
    {return THIS()->min_distance_impl(tree, candidates);}

    /**
    * Maximum distance between a point and the tree.
    * \param[in] point Point to compute the maximum distance to.
    * \param[out] candidates Maximum distance candidates.
    * \return The maximum distance between the point and the tree.
    */
    Real max_distance(Point const & point, IndexSet & candidates) const
    {return THIS()->max_distance_impl(point, candidates);}

    /**
    * Maximum distance between an axis-aligned box and the tree.
    * \param[in] box Axis-aligned box to compute the maximum distance to.
    * \param[out] candidates Maximum distance candidates.
    * \return The maximum distance between the box and the tree.
    */
    Real max_distance(Box const & box, IndexSet & candidates) const
    {return THIS()->max_distance_impl(box, candidates);}

    /**
    * Maximum distance between an current tree and another tree.
    * \param[in] tree Tree to compute the maximum distance to.
    * \param[out] candidates Maximum distance candidates.
    * \return The maximum distance between the trees.
    */
    Real max_distance(DerivedTree const & tree, IndexSet & candidates) const
    {return THIS()->max_distance_impl(tree, candidates);}

  }; // Tree

} // namespace AABBtree

#endif // INCLUDE_AABBTREE_TREE_HXX
