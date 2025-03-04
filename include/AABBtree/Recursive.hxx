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
    using Box = typename Tree<Real, N, Recursive<Real, N>>::Box;
    using BoxPtr = typename Tree<Real, N, Recursive<Real, N>>::BoxPtr;
    using BoxUPtrVec = typename Tree<Real, N, Recursive<Real, N>>::BoxUPtrVec;
    using Point = typename Tree<Real, N, Recursive<Real, N>>::Point;
    using Set = typename Tree<Real, N, Recursive<Real, N>>::Set;

    using TreePtr = std::unique_ptr<Recursive<Real, N>>;
    using Children = std::pair<BoxPtr, BoxPtr>;

    Children m_children; /** Pointers to the children (left and right) boxes. */
    BoxPtr m_box; /** Pointer to the current box. */
    BoxUPtrVec m_objects; /** Object axis-aligned bounding boxes. */

    Recursive(Recursive<Real, N> const &) = delete; /**< Copy constructor. */
    Recursive<Real, N> & operator=(Recursive<Real, N> const &) = delete; /**< Copy assignment operator. */

  public:
    /**
    * Class destructor for the \em recursive axis-aligned bounding box tree.
    */
    ~Recursive() {};

    /**
    * Class constructor for the \em recursive axis-aligned bounding box tree.
    */
    Recursive() : m_children(nullptr, nullptr)
    {
      this->m_objects.reserve(this->m_max_nodal_objects);
    };

    /**
    * Clear the tree.
    */
    void clear_impl()
    {
      this->m_objects.clear();
      this->m_children.first->clear();
      this->m_children.second->clear();
      this->m_children.first = nullptr;
      this->m_children.second = nullptr;
    }

    /**
    * Check if tree is empty.
    * \return True if the tree is empty, false otherwise.
    */
    bool is_empty_impl() const
    {
      return !this->m_children.first && !this->m_children.second && this->m_box.is_empty();
    }

    /**
    * Build the tree given the bounding boxes.
    * \param[in] boxes Bounding boxes to build the tree from.
    * \return True if the tree was built successfully, false otherwise.
    */
    void build_impl(const BoxUPtrVec & boxes)
    {
      // Check if the input vector is empty
      if (boxes.empty()) {return;}

      // Check if the input vector contains only one box
      if (boxes.size() == static_cast<Integer>(1))
      {
        this->m_box.emplace_back(boxes.front());
        return;
      }

      // Compute the bounding box of the input vector
      this->m_box = std::make_unique<Box>();
      this->m_box->extend(boxes);

      // Cut on longest edge of the bounding box
      Integer n_cut{this->m_box->max_dimension()};
      Real x_min{this->m_box->min(n_cut)};
      Real x_max{this->m_box->max(n_cut)};
      Real x_cut{(x_max + x_min)/static_cast<Real>(2.0)};

      // Collect the boxes that have a size similar to the current box
      // CHECKME
     this->m_objects.clear();
     for (const Box & box : boxes)
     {
        Real x_mid{(box->min(n_cut) + box->max(n_cut))/static_cast<Real>(2.0)};
        if (x_mid > x_cut) {this->m_objects.emplace_back(box);}
      }

      // Split the boxes into left and right
      std::vector<Box> left_boxes;
      std::vector<Box> right_boxes;
      std::vector<Box>::const_iterator it;
      for (it = boxes.begin(); it != boxes.end(); ++it)
      {
        Real x_mid{((*it)->min(n_cut) + (*it)->max(n_cut))/static_cast<Real>(2.0)};
        if (x_mid > x_cut) {left_boxes.emplace_back(*it);}
        else {right_boxes.emplace_back(*it);}
      }

      // Build the left and right trees
      TreePtr left_tree{std::make_unique<Recursive>()};
      TreePtr right_tree{std::make_unique<Recursive>()};
      left_tree->build(left_boxes);
      right_tree->build(right_boxes);
      if (!left_tree->is_empty()) {this->m_children.first(left_tree);}
      if (!right_tree->is_empty()) {this->m_children.second(right_tree);}
    }

    /**
    * Build the tree given the bounding boxes.
    * \param[in] boxes Bounding boxes to build the tree from.
    * \return True if the tree was built successfully, false otherwise.
    */
    void build_impl(const std::vector<Box> & boxes)
    {
      this->build_impl(BoxUPtrVec(boxes.begin(), boxes.end()));
    }

    /**
    * Intersect the tree with a point.
    * \param[in] point Point to intersect with.
    * \param[out] candidates Intersection result (bounding box indexes).
    * \return True if the point intersects the tree, false otherwise.
    */
    bool intersect_impl(const Point & point, Set & candidates) const
    {
      return true;
    }

    /**
    * Intersect the tree with an axis-aligned box.
    * \param[in] box Axis-aligned box to intersect with.
    * \param[out] candidates Intersection result (bounding box indexes).
    * \return True if the point intersects the tree, false otherwise.
    */
    bool intersect_impl(const Box & box, Set & candidates) const
    {
      return true;
    }

    /**
    * Intersect the tree with another tree.
    * \param[in] tree Tree to intersect with.
    * \param[out] candidates Intersection result (bounding box indexes).
    * \return True if the point intersects the tree, false otherwise.
    */
    bool intersect_impl(const Recursive & tree, Set & candidates) const
    {
      return true;
    }

    /**
    * Minimum distance between a point and the tree.
    * \param[in] point Point to compute the minimum distance to.
    * \param[out] candidates Minimum distance candidates.
    * \return The minimum distance between the point and the tree.
    */
    Real min_distance_impl(const Point & point, Set & candidates) const
    {
      return static_cast<Real>(0.0);
    }

    /**
    * Minimum distance between an axis-aligned box and the tree.
    * \param[in] box Axis-aligned box to compute the minimum distance to.
    * \param[out] candidates Minimum distance candidates.
    * \return The minimum distance between the box and the tree.
    */
    Real min_distance_impl(const Box & box, Set & candidates) const
    {
      return static_cast<Real>(0.0);
    }

    /**
    * Minimum distance between an current tree and another tree.
    * \param[in] tree Tree to compute the minimum distance to.
    * \param[out] candidates Minimum distance candidates.
    * \return The minimum distance between the trees.
    */
    Real min_distance_impl(const Recursive & tree, Set & candidates) const
    {
      return static_cast<Real>(0.0);
    }

    /**
    * Maximum distance between a point and the tree.
    * \param[in] point Point to compute the maximum distance to.
    * \param[out] candidates Maximum distance candidates.
    * \return The maximum distance between the point and the tree.
    */
    Real max_distance_impl(const Point & point, Set & candidates) const
    {
      return static_cast<Real>(0.0);
    }

    /**
    * Maximum distance between an axis-aligned box and the tree.
    * \param[in] box Axis-aligned box to compute the maximum distance to.
    * \param[out] candidates Maximum distance candidates.
    * \return The maximum distance between the box and the tree.
    */
    Real max_distance_impl(const Box & box, Set & candidates) const
    {
      return static_cast<Real>(0.0);
    }

    /**
    * Maximum distance between an current tree and another tree.
    * \param[in] tree Tree to compute the maximum distance to.
    * \param[out] candidates Maximum distance candidates.
    * \return The maximum distance between the trees.
    */
    Real max_distance_impl(const Recursive & tree, Set & candidates) const
    {
      return static_cast<Real>(0.0);
    }

  }; // Recursive

  template class Recursive<float, 1>;
  template class Recursive<float, 2>;
  template class Recursive<float, 3>;
  template class Recursive<float, 4>;
  template class Recursive<double, 1>;
  template class Recursive<double, 2>;
  template class Recursive<double, 3>;
  template class Recursive<double, 4>;

} // namespace AABBtree

#endif // INCLUDE_AABBTREE_RECURSIVE_HXX
