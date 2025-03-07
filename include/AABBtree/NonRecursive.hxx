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

#ifndef INCLUDE_AABBTREE_NONRECURSIVE_HXX
#define INCLUDE_AABBTREE_NONRECURSIVE_HXX

#include "AABBtree/Box.hxx"
#include "AABBtree/Tree.hxx"

namespace AABBtree {

  /**
  * \brief A class representing a \em non-recursive axis-aligned bounding box tree (AABBtree).
  *
  * The NonRecursive class provides an efficient way to store and query a set of axis-aligned bounding
  * boxes. It supports various operations such as building the tree, adding axis-aligned bounding
  * boxes, and performing intersection queries. The NonRecursive class is templated on the type of the
  * real numbers and the dimensions of the bounding boxes and is based, as the name suggests, on a
  * novel non-recursive algorithm.
  * \tparam Real Type of the scalar coefficients
  * \tparam N Dimension of the ambient space.
  */
  template <typename Real, Integer N>
  class NonRecursive : public Tree<Real, N, NonRecursive<Real, N>> {

    friend class Tree<Real, N, NonRecursive<Real, N>>; // Allow access to the base class members

    // Basic types definitions
    using Box = Box<Real, N>;
    using BoxUniquePtr = BoxUniquePtr<Real, N>;
    using BoxUniquePtrList = BoxUniquePtrList<Real, N>;
    using Vector = Vector<Real, N>;
    using Point = Point<Real, N>;

    typedef struct AABBsubtree {
      Box box; /**< Bounding box of the subtree. */
      Box box_long; /**< Bounding box of long boxes. */
      Integer box_ptr; /**< Pointer to the first box in the reordering map of boxes. */
      Integer box_num; /**< Number of boxes in the subtree. */
      Integer parent; /**< Root node of the subtree. */
      Integer child_l; /**< Left child of the subtree. */
      Integer child_r; /**< Right child of the subtree. */
    } AABBnode; /**< Structure representing a node of the AABBtree. */

    // Tree hierarchy
    std::vector<AABBsubtree> m_tree_structure; /**< Tree structure. */
    IndexList m_tree_boxes_map; /**< Reordering between the vector of boxes and the tree internal structure. */
    bool m_dumping_mode{true}; /**< Enable dumping while building the tree. */

    // Stack for non-recursive tree building and query
    mutable IndexList m_stack; /**< Tree stack. */

  public:
    /**
    * Class destructor for the \em non-recursive axis-aligned bounding box tree.
    */
    ~NonRecursive() = default;

    /**
    * Class constructor for the \em non-recursive axis-aligned bounding box tree.
    */
    NonRecursive() = default;

    /**
    * Get the tree structure.
    * \return The tree structure.
    */
    std::vector<AABBsubtree> const & structure() const {return this->m_tree_structure;}

    /**
    * Get the i-th tree node.
    * \param[in] i Index of the node.
    * \return The i-th tree node.
    */
    AABBsubtree const & node(Integer const i) const {return this->m_tree_structure[i];}

    /**
    * Get the number of nodes in the tree.
    * \return The number of nodes in the tree.
    */
    Integer size() const {return this->m_tree_structure.size();}

    /**
    * Enable dumping mode while building the tree.
    */
    void enable_dumping_mode() {this->m_dumping_mode = true;}

    /**
    * Disable dumping mode while building the tree.
    */
    void disable_dumping_mode() {this->m_dumping_mode = false;}

    /**
    * Set dumping mode while building the tree.
    * \param[in] mode Dumping mode.
    */
    void dumping_mode(bool const mode) {this->m_dumping_mode = mode;}

    /**
    * Check if tree is empty.
    * \return True if the tree is empty, false otherwise.
    */
    bool is_empty_impl() const {return this->m_tree_structure.empty();}

    /**
    * Clear the tree.
    */
    void clear_impl()
    {
      this->m_tree_structure.clear();
      this->m_tree_boxes_map.clear();
      this->m_stack.clear();
    }

    /**
    * Print the tree internal structure in an output stream.
    * \param[in] os Output stream to print the tree to.
    */
    void print_impl(OutStream & os) const
    {
      Integer leafs{0};
      for (Integer i{0}; i < this->m_tree_structure.size(); ++i ) {
        if (this->m_tree_structure[i].child_l == -1) {++leafs;}
        if (this->m_tree_structure[i].child_r == -1) {++leafs;}
      }
      os <<
        "AABB tree internal structure" << std::endl <<
        "\tAmbient dimension : " + std::to_string(N) << std::endl <<
        "\tBasic type        : " + std::string(typeid(Real).name()) << std::endl <<
        "\tNumber of nodes   : " + std::to_string(this->m_tree_structure.size()) << std::endl <<
        "\tNumber of leafs   : " + std::to_string(leafs) << std::endl <<
        "\tNumber of objects : " + std::to_string(this->m_boxes.size()) << std::endl;
    }

    /*
                                              m_tree_structure
                     ┌──────┐                  ┌──────┐
                     │      │                  │ CD1  │
                     │      │                  │ CD2  │
                     │      │                  │ CD1  │
    m_nodes_id   --> │ BBOX │    m_children -> │ CD2  │
                     │      │                  │      │
                     │      │                  │      │
                     └──────┘

    */

    /**
    * Build the tree given a set of bounding boxes internally stored.
    */
    void build_impl() {

      // Collect the original object boxes
      BoxUniquePtrList const & boxes {*this->m_boxes};

      // Clear tree structure
      Integer num_boxes{boxes.size()};
      this->m_tree_structure.clear();
      this->m_tree_structure.reserve(static_cast<size_t>(2*num_boxes));

      // Setup the root node
      AABBnode root;
      root.parent  = -1;
      root.child_l = -1;
      root.child_r = -1;
      root.box_ptr = 0;
      root.box_num = 0;
      this->m_tree_structure.emplace_back(root);

      // Setup the boxes map and compute the root box
      Integer depth{std::ceil(std::log2(num_boxes))};
      root.box.set_empty();
      this->m_tree_boxes_map.reserve(2*depth);
      for (auto const & b : boxes) {
        root.box.extend(b);
        this->m_tree_boxes_map.emplace_back(root.box_num++);
      }

      // Setup the stack
      this->m_stack.clear();
      this->m_stack.reserve(2*num_boxes + 1);
      this->m_stack.emplace_back(0);

      // Setup the dump counter (axis to try)
      Integer dump{0};

      // Main loop that divide the nodes iteratively until all constraints satisfied
      while (!this->m_stack.empty())
      {
        // Pop the node from stack
        Integer const id{this->m_stack.back()}; this->m_stack.pop_back();

        // Get the node
        AABBnode & node {this->m_tree_structure[id]};

        // If the node has less than the maximum number of objects, skip it
        if (node.box_num < this->m_max_nodal_objects) {continue;}

        // Compute the separation line and tolerance
        Vector sizes;
        Eigen::Vector<Integer, N> sorting;
        node.box.sort_axes_length(sizes, sorting);
        Integer axis{sorting[dump]};
        Real separation_line{node.box.barycenter(axis)};
        Real separation_tolerance{sizes[axis] * this->m_separation_ratio_tolerance};

        // Separate short and long boxes and compute short boxes barycenter
        Integer n_long{0};
        Integer n_left{0};
        Integer n_right{0};
        Integer id_ini{node.box_ptr};
        Integer id_end{node.box_ptr + node.box_num};
        Real barycenter {0.0};
        while (id_ini < id_end) {
          Box const & box_id{*boxes[this->m_tree_boxes_map[id_ini]]};
          Integer side{static_cast<Integer>(box_id.which_side(separation_line, separation_tolerance, axis))};
          switch (side) {
            case static_cast<Integer>(Box::Side::LEFT): ++n_left; // Left boxes are moved to the end
              std::swap(this->m_tree_boxes_map[id_ini], this->m_tree_boxes_map[--id_end]); break;
            case static_cast<Integer>(Box::Side::RIGHT): ++n_right; // Right boxes are moved to the end
              std::swap(this->m_tree_boxes_map[id_ini], this->m_tree_boxes_map[--id_end]); break;
            default: ++n_long; ++id_ini;
          }
          barycenter += box_id.max(axis) + box_id.min(axis);
        }
        barycenter /= 2.0*node.box_num;

        // Check if the leaves are balanced, if not try to separate boxes on a new separation line
        if (n_long > this->m_max_nodal_objects ||
            std::min(n_left, n_right) < this->m_balance_ratio_tolerance * std::max(n_left, n_right) ) {

          // Perform coputation only if barycenter is well separated from the previous separation line
          if (std::abs(barycenter - separation_line) > separation_tolerance) {
            n_long = n_left = n_right = 0;
            id_ini = node.box_ptr;
            id_end = node.box_ptr + node.box_num;
            while (id_ini < id_end) {
              Box const & box_id{*boxes[this->m_tree_boxes_map[id_ini]]};
              Integer side{static_cast<Integer>(box_id.which_side(barycenter, separation_tolerance, axis))};
              switch (side) {
                case static_cast<Integer>(Box::Side::LEFT): ++n_left;
                  std::swap(this->m_tree_boxes_map[id_ini], this->m_tree_boxes_map[--id_end]); break;
                case static_cast<Integer>(Box::Side::RIGHT): ++n_right;
                  std::swap(this->m_tree_boxes_map[id_ini], this->m_tree_boxes_map[--id_end]); break;
                default: ++n_long; ++id_ini;
              }
            }
          }
        }

        // If the left and right children are yet not balanced, dump the splitting axis
        if (this->m_dumping_mode && n_long > this->m_max_nodal_objects) {
          if (std::min(n_left, n_right) < this->m_balance_ratio_tolerance * std::max(n_left, n_right)) {
            if (++dump < N) {this->m_stack.push_back(id); {continue;}}
        }}

        // Separate the left and right boxes
        n_left = n_right = 0;
        id_ini = node.box_ptr + n_long;
        id_end = node.box_ptr + node.box_num;
        while (id_ini < id_end) {
          Integer ipos{this->m_tree_boxes_map[id_ini]};
          if (boxes[ipos]->barycenter(axis) < separation_line) {
            ++id_ini; ++n_left; // In right position do nothing
          } else {
            --id_end; ++n_right; // In right position swap the current box with the last one
            std::swap(this->m_tree_boxes_map[id_ini], this->m_tree_boxes_map[id_end]);
          }
        }

        // If the left and right children have few boxes, they are leaves
        if (n_left <= this->m_max_nodal_objects && n_right <= this->m_max_nodal_objects) {continue;}

        // Set the left and right children indexes
        node.child_l = static_cast<Integer>(this->m_tree_structure.size() + 0);
        node.child_r = static_cast<Integer>(this->m_tree_structure.size() + 1);

        // Finalize the root node setup (left and right children)
        AABBnode node_l;
        node_l.parent  = id; // Current node
        node_l.child_l = -1; // Default (leaf)
        node_l.child_r = -1; // Default (leaf)
        node_l.box_num = n_left;

        AABBnode node_r;
        node_r.parent  = id; // Current node
        node_r.child_l = -1; // Default (leaf)
        node_r.child_r = -1; // Default (leaf)
        node_r.box_num = n_right;

        // Compute the bounding box of the long boxes, and left and right children
        Integer j{node.box_ptr};
        node.box_num = n_long;
        node.box_long.set_empty();
        for (Integer i{0}; i < n_long; ++i) node.box_long.extend(*boxes[m_tree_boxes_map[j++]]);
        node_l.box.set_empty(); node_l.box_ptr = j;
        for (Integer i{0}; i < n_left; ++i) node_l.box.extend(*boxes[m_tree_boxes_map[j++]]);
        node_r.box.set_empty(); node_r.box_ptr = j;
        for (Integer i{0}; i < n_right; ++i) node_r.box.extend(*boxes[m_tree_boxes_map[j++]]);

        // Push nodes on tree structure
        this->m_tree_structure.emplace_back(node_l);
        this->m_tree_structure.emplace_back(node_r);

        // Push children on stack
        this->m_stack.emplace_back(node.child_l);
        this->m_stack.emplace_back(node.child_r);

        // Reset the dump counter
        dump = 0;
      }

    }

    /**
    * Intersect the tree with a point.
    * \param[in] point Point to intersect with.
    * \param[out] candidates Intersection result (bounding box indexes).
    * \return True if the point intersects the tree, false otherwise.
    */
    bool intersect_impl(Point const & point, IndexSet & candidates) const
    {
      // Reset statistics
      this->m_num_check = 0;

      // Return if the tree is empty
      if (this->is_empty()) {return false;}

      // Collect the original object boxes
      const BoxUniquePtrList & boxes {*this->m_boxes};
      Integer num_boxes{boxes.size()};

      // Setup the stack
      this->m_stack.clear();
      this->m_stack.reserve(2*num_boxes + 1);
      this->m_stack.emplace_back(0);

      // Main loop that checks the intersection iteratively
      while (!this->m_stack.empty())
      {
        // Pop the node from stack
        Integer const id{this->m_stack.back()}; this->m_stack.pop_back();

        // Get the node
        AABBnode & node {this->m_tree_structure[id]};

        // If the point is not in the bounding box, skip the node
        ++this->m_num_check;
        if (!node.box.contains(point)) {continue;}

        // Intersect the point with the long boxes on the leaves
        Integer const id_ini{node.box_ptr};
        Integer const id_end{node.box_ptr + node.box_num};
        for (Integer i{id_ini}; i < id_end; ++i) {
          Integer const pos{this->m_tree_boxes_map[i]};
          Box const & box_i {*boxes[pos]};
          ++this->m_num_check;
          if (box_i.contains(point)) {candidates.insert(pos);}
        }

        // Push children on the stack if thay are not leaves
        if (node.child_l > 0) {this->m_stack.emplace_back(node.child_l);}
        if (node.child_r > 0) {this->m_stack.emplace_back(node.child_r);}
      }
      return !candidates.empty();
    }

    /**
    * Intersect the tree with an axis-aligned box.
    * \param[in] box Axis-aligned box to intersect with.
    * \param[out] candidates Intersection result (bounding box indexes).
    * \return True if the point intersects the tree, false otherwise.
    */
    bool intersect_impl(Box const & box, IndexSet & candidates) const
    {
      // Reset statistics
      this->m_num_check = 0;

      // Collect the original object boxes
      const BoxUniquePtrList & boxes {*this->m_boxes};
      Integer num_boxes{boxes.size()};

      // Return if the tree is empty
      if (this->is_empty()) {return false;}

      // Setup the stack
      this->m_stack.clear();
      this->m_stack.reserve(2*num_boxes + 1);
      this->m_stack.emplace_back(0);

      // Main loop that checks the intersection iteratively
      while (!this->m_stack.empty())
      {
        // Pop the node from stack
        Integer const id{this->m_stack.back()}; this->m_stack.pop_back();

        // Get the node
        AABBnode & node {this->m_tree_structure[id]};

        // If the point is not in the bounding box, skip the node
        ++this->m_num_check;
        if (!node.box.intersects(box)) {continue;}

        // Intersect the box with the long boxes on the leaves
        Integer const id_ini{node.box_ptr};
        Integer const id_end{node.box_ptr + node.box_num};
        for (Integer i{id_ini}; i < id_end; ++i) {
          Integer const pos{this->m_tree_boxes_map[i]};
          Box const & box_i {*boxes[pos]};
          ++this->m_num_check;
          if (box_i.intersects(box)) {candidates.insert(pos);}
        }

        // Push children on the stack if thay are not leaves
        if (node.child_l > 0) {this->m_stack.emplace_back(node.child_l);}
        if (node.child_r > 0) {this->m_stack.emplace_back(node.child_r);}
      }
      return !candidates.empty();
    }

  }; // NonRecursive

} // namespace AABBtree

#endif // INCLUDE_AABBTREE_NONRECURSIVE_HXX
