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

#ifndef INCLUDE_AABBTREE_RECURSIVE_HXX
#define INCLUDE_AABBTREE_RECURSIVE_HXX

#include "AABBtree/Box.hxx"
#include "AABBtree/Tree.hxx"

namespace AABBtree {

  /**
  * \brief A class representing a \em recursive axis-aligned bounding box tree (AABBtree).
  *
  * The Recursive class provides an efficient way to store and query a set of axis-aligned bounding
  * boxes. It supports various operations such as building the tree, adding axis-aligned bounding
  * boxes, and performing intersection queries. The Recursive class is templated on the type of the
  * real numbers and the dimensions of the bounding boxes and is based, as the name suggests, on a
  * novel recursive algorithm.
  * \tparam Real Type of the scalar coefficients
  * \tparam N Dimension of the ambient space.
  */

  template <typename Real, Integer N>
  class Recursive : public Tree<Real, N, Recursive<Real, N>>
  {
    friend class Tree<Real, N, Recursive<Real, N>>; // Allow access to the base class members

    // Basic types definitions
    using Point = Point<Real, N>;
    using Vector = Vector<Real, N>;
    using Box = Box<Real, N>;
    using BoxUniquePtr = BoxUniquePtr<Real, N>;
    using BoxUniquePtrList = BoxUniquePtrList<Real, N>;
    using Children = std::pair<BoxUniquePtr, BoxUniquePtr>;

    Children m_children; /** Pointers to the children (left and right) boxes. */
    BoxUniquePtr m_box; /** Pointer to the current box. */
    BoxUniquePtrList m_objects; /** Object axis-aligned bounding boxes. */

  public:
    Recursive(Recursive const &) = delete; /**< Copy constructor. */
    Recursive & operator=(Recursive const &) = delete; /**< Copy assignment operator. */

    /**
     * Class destructor for the \em recursive axis-aligned bounding box tree.
     */
    ~Recursive() = default;

    /**
     * Class constructor for the \em recursive axis-aligned bounding box tree.
     */
    Recursive() : m_children(nullptr, nullptr) {
      this->m_objects.reserve(this->m_max_nodal_objects);
    };

    /**
     * Clear the tree.
     */
    void
    clear_impl() {
      this->m_objects.clear();
      this->m_children.first.release();
      this->m_children.second.release();
      this->m_children.first = nullptr;
      this->m_children.second = nullptr;
    }

    /**
     * Check if tree is empty.
     * \return True if the tree is empty, false otherwise.
     */
    bool
    is_empty_impl() const {
      return !this->m_children.first && !this->m_children.second && this->m_box.is_empty();
    }

    /**
     * Build the tree given the bounding boxes.
     * \param[in] boxes Bounding boxes to build the tree from.
     * \return True if the tree was built successfully, false otherwise.
     */
    void
    build_impl(BoxUniquePtrList const & boxes) {
      // Check if the input vector is empty
      if (boxes.empty() ) return;

      // Check if the input vector contains only one box
      if (boxes.size() == static_cast<Integer>(1)) {
        this->m_box.emplace_back(boxes.front());
        return;
      }

      // Compute the bounding box of the input vector
      this->m_box = std::make_unique<Box>();
      this->m_box->extend(boxes);

      // Cut on longest edge of the bounding box
      Integer n_cut{ this->m_box->max_dimension()};
      Real    x_min{ this->m_box->min(n_cut)};
      Real    x_max{ this->m_box->max(n_cut)};
      Real    x_cut{ (x_max + x_min)/2.0};

      // Collect the boxes that have a size similar to the current box
      // CHECKME
     this->m_objects.clear();
     for (Box const & box : boxes)
     {
        Real x_mid{(box->min(n_cut) + box->max(n_cut))/2.0};
        if (x_mid > x_cut) {this->m_objects.emplace_back(box);}
      }

      // Split the boxes into left and right
      BoxUniquePtrList left_boxes;
      BoxUniquePtrList right_boxes;
      for (auto const & b : boxes ) {
        Real x_mid{(b.min(n_cut) + b.max(n_cut))*static_cast<Real>(0.5)};
        if (x_mid > x_cut) {left_boxes.emplace_back(b);}
        else {right_boxes.emplace_back(b);}
      }

      // Build the left and right trees
      TreePtr left_tree{std::make_unique<Recursive>()};
      TreePtr right_tree{std::make_unique<Recursive>()};
      left_tree->build(left_boxes);
      right_tree->build(right_boxes);
      if (!left_tree->is_empty()) {m_children.first(left_tree);}
      if (!right_tree->is_empty()) {m_children.second(right_tree);}
    }

    /**
     * Build the tree given the bounding boxes.
     * \param[in] boxes Bounding boxes to build the tree from.
     * \return True if the tree was built successfully, false otherwise.
     */
    void
    build_impl(vector<Box> const & boxes) {
      this->build_impl(BoxUniquePtrList(boxes.begin(), boxes.end()));
    }

    /**
     * Intersect the tree with a point.
     * \param[in] point Point to intersect with.
     * \param[out] candidates Intersection result (bounding box indexes).
     * \return True if the point intersects the tree, false otherwise.
     */
    bool
    intersect_impl(Point const & point, Set & candidates) const {
      return true;
    }

    /**
     * Intersect the tree with an axis-aligned box.
     * \param[in] box Axis-aligned box to intersect with.
     * \param[out] candidates Intersection result (bounding box indexes).
     * \return True if the point intersects the tree, false otherwise.
     */
    bool
    intersect_impl(Box const & box, Set & candidates ) const {
      return true;
    }

    /**
     * Intersect the tree with another tree.
     * \param[in] tree Tree to intersect with.
     * \param[out] candidates Intersection result (bounding box indexes).
     * \return True if the point intersects the tree, false otherwise.
     */
    bool
    intersect_impl(Recursive const & tree, Set & candidates) const {
      return true;
    }

    /**
    * Minimum distance between a point and the tree.
    * \param[in] point Point to compute the minimum distance to.
    * \param[out] candidates Minimum distance candidates.
    * \return The minimum distance between the point and the tree.
    */
    Real min_distance_impl(Point const & point, Set & candidates) const
    {
      return static_cast<Real>(0.0);
    }

    /**
    * Minimum distance between an axis-aligned box and the tree.
    * \param[in] box Axis-aligned box to compute the minimum distance to.
    * \param[out] candidates Minimum distance candidates.
    * \return The minimum distance between the box and the tree.
    */
    Real min_distance_impl(Box const & box, Set & candidates) const
    {
      return static_cast<Real>(0.0);
    }

    /**
    * Minimum distance between an current tree and another tree.
    * \param[in] tree Tree to compute the minimum distance to.
    * \param[out] candidates Minimum distance candidates.
    * \return The minimum distance between the trees.
    */
    Real min_distance_impl(Recursive const & tree, Set & candidates) const
    {
      return static_cast<Real>(0.0);
    }

    /**
    * Maximum distance between a point and the tree.
    * \param[in] point Point to compute the maximum distance to.
    * \param[out] candidates Maximum distance candidates.
    * \return The maximum distance between the point and the tree.
    */
    Real max_distance_impl(Point const & point, Set & candidates) const
    {
      return static_cast<Real>(0.0);
    }

    /**
    * Maximum distance between an axis-aligned box and the tree.
    * \param[in] box Axis-aligned box to compute the maximum distance to.
    * \param[out] candidates Maximum distance candidates.
    * \return The maximum distance between the box and the tree.
    */
    Real max_distance_impl(Box const & box, Set & candidates) const
    {
      return static_cast<Real>(0.0);
    }

    /**
    * Maximum distance between an current tree and another tree.
    * \param[in] tree Tree to compute the maximum distance to.
    * \param[out] candidates Maximum distance candidates.
    * \return The maximum distance between the trees.
    */
    Real max_distance_impl(Recursive const & tree, Set & candidates) const
    {
      return static_cast<Real>(0.0);
    }

  }; // Recursive

} // namespace AABBtree

#endif // INCLUDE_AABBTREE_RECURSIVE_HXX
