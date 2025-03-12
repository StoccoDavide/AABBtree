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
      // Count the number of leafs
      Integer leafs{0}, long_boxes{0};
      for (Integer i{0}; i < this->size(); ++i ) {
        if (this->m_tree_structure[i].child_l == -1) {++leafs;}
        if (this->m_tree_structure[i].child_r == -1) {++leafs;}
        if (this->m_tree_structure[i].child_l > 0 && this->m_tree_structure[i].child_r > 0)
          {long_boxes += this->m_tree_structure[i].box_num;}
      }

      // Count the number of left nodes
      Integer left_nodes{0};
      Integer left_leafs{0};
      this->m_stack.clear();
      this->m_stack.reserve(2*this->size() + 1);
      this->m_stack.emplace_back(this->m_tree_structure[0].child_l);
      while (!this->m_stack.empty())
      {
        Integer const id{this->m_stack.back()}; this->m_stack.pop_back();
        AABBnode const & node {this->m_tree_structure[id]};
        if (node.child_l == -1) {++left_leafs;}
        else {this->m_stack.emplace_back(node.child_l);}
        if (node.child_r == -1) {++left_leafs;}
        else {this->m_stack.emplace_back(node.child_r);}
        ++left_nodes;
      }

      // Count the number of right nodes
      Integer right_nodes{0};
      Integer right_leafs{0};
      this->m_stack.clear();
      this->m_stack.reserve(2*this->size() + 1);
      this->m_stack.emplace_back(this->m_tree_structure[0].child_r);
      while (!this->m_stack.empty())
      {
        Integer const id{this->m_stack.back()}; this->m_stack.pop_back();
        AABBnode const & node {this->m_tree_structure[id]};
        if (node.child_l == -1) {++right_leafs;}
        else {this->m_stack.emplace_back(node.child_l);}
        if (node.child_r == -1) {++right_leafs;}
        else {this->m_stack.emplace_back(node.child_r);}
        ++right_nodes;
      }

      os <<
        "AABB tree internal structure ---------" << std::endl <<
        "  Ambient dimension     : " << N << std::endl <<
        "  Basic type            : " << typeid(Real).name() << std::endl <<
        "  Number of nodes       : " << this->size() << std::endl <<
        "  Number of leafs       : " << leafs << std::endl <<
        "  Number of long boxes  : " << long_boxes << std::endl <<
        "  Number of left nodes  : " << left_nodes << std::endl <<
        "  Number of left leafs  : " << left_leafs << std::endl <<
        "  Number of right nodes : " << right_nodes << std::endl <<
        "  Number of right leafs : " << right_leafs << std::endl <<
        "  Balance ratio         : " << static_cast<Real>(left_leafs)/static_cast<Real>(left_nodes) << std::endl <<
        "  Number of objects     : " << this->m_boxes->size() << std::endl <<
        "--------------------------------------" << std::endl;
    }

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

      // Setup the boxes map and compute the root box
      Integer depth{std::ceil(std::log2(num_boxes))};
      root.box.set_empty();
      this->m_tree_boxes_map.reserve(2*depth);
      for (BoxUniquePtr const & box : boxes) {
        root.box.extend(*box);
        this->m_tree_boxes_map.emplace_back(root.box_num++);
      }
      this->m_tree_structure.emplace_back(root);

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
        Real separation_line{node.box.baricenter(axis)};
        Real separation_tolerance{sizes[axis] * this->m_separation_ratio_tolerance};

        // Separate short and long boxes and compute short boxes baricenter
        Integer n_long{0};
        Integer n_left{0};
        Integer n_right{0};
        Integer id_ini{node.box_ptr};
        Integer id_end{node.box_ptr + node.box_num};
        Real baricenter {0.0};
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
          baricenter += box_id.max()[axis] + box_id.min()[axis];
        }
        baricenter /= 2.0*node.box_num;

        // Check if the leafs are balanced, if not try to separate boxes on a new separation line
        if (n_long > this->m_max_nodal_objects ||
            std::min(n_left, n_right) < this->m_balance_ratio_tolerance * std::max(n_left, n_right) ) {

          // Perform coputation only if baricenter is well separated from the previous separation line
          if (std::abs(baricenter - separation_line) > separation_tolerance) {
            n_long = n_left = n_right = 0;
            id_ini = node.box_ptr;
            id_end = node.box_ptr + node.box_num;
            while (id_ini < id_end) {
              Box const & box_id{*boxes[this->m_tree_boxes_map[id_ini]]};
              Integer side{static_cast<Integer>(box_id.which_side(baricenter, separation_tolerance, axis))};
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
            if (dump < N-1) {++dump; this->m_stack.push_back(id); {continue;}}
        }}

        // Reset the dump counter
        dump = 0;

        // Separate the left and right boxes
        n_left = n_right = 0;
        id_ini = node.box_ptr + n_long;
        id_end = node.box_ptr + node.box_num;
        while (id_ini < id_end) {
          Integer ipos{this->m_tree_boxes_map[id_ini]};
          if (boxes[ipos]->baricenter(axis) < separation_line) {
            ++id_ini; ++n_left; // In right position do nothing
          } else {
            --id_end; ++n_right; // In right position swap the current box with the last one
            std::swap(this->m_tree_boxes_map[id_ini], this->m_tree_boxes_map[id_end]);
          }
        }

        // If the left and right children have few boxes, they are leafs
        if (n_left <= this->m_max_nodal_objects && n_right <= this->m_max_nodal_objects) {continue;}

        // Set the left and right children indexes
        node.child_l = static_cast<Integer>(this->size() + 0);
        node.child_r = static_cast<Integer>(this->size() + 1);

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
      }
    }

    /**
    * Intersect the tree with an object.
    * \param[in] obj Object to intersect with.
    * \param[out] candidates Intersection result (bounding box indexes).
    * \return True if the object intersects the tree, false otherwise.
    * \tparam Object Type of the object to intersect with.
    * \note Object must have a method \c intersects that computes the intersection with a box.
    */
    template <typename Object>
    bool intersect_impl(Object const & obj, IndexSet & candidates) const
    {
      // Reset statistics
      this->m_check_counter = 0;

      // Return if the tree is empty
      if (this->is_empty()) {return false;}

      // Collect the original object boxes
      BoxUniquePtrList const & boxes {*this->m_boxes};

      // Setup the stack
      this->m_stack.clear();
      this->m_stack.reserve(2*this->size() + 1);
      this->m_stack.emplace_back(0);

      // Main loop that checks the intersection iteratively
      candidates.clear();
      while (!this->m_stack.empty())
      {
        // Pop the node from stack
        Integer const id{this->m_stack.back()}; this->m_stack.pop_back();

        // Get the node
        AABBnode const & node {this->m_tree_structure[id]};

        // If the object is not in the bounding box, skip the node
        ++this->m_check_counter;
        if (!node.box.intersects(obj)) {continue;}

        // Intersect the object with the long boxes on the node
        Integer const id_ini{node.box_ptr};
        Integer const id_end{node.box_ptr + node.box_num};
        for (Integer i{id_ini}; i < id_end; ++i) {
          Integer const pos{this->m_tree_boxes_map[i]};
          ++this->m_check_counter;
          if (boxes[pos]->intersects(obj)) {candidates.insert(pos);}
        }

        // Push children on the stack if thay are not leafs
        if (node.child_l > 0) {this->m_stack.emplace_back(node.child_l);}
        if (node.child_r > 0) {this->m_stack.emplace_back(node.child_r);}
      }

      // Return true if the object intersects the tree
      return !candidates.empty();
    }

    /**
    * Intersect the tree with another tree.
    * \param[in] tree Tree to intersect with.
    * \param[out] candidates Intersection result (bounding box indexes).
    * \return True if the point intersects the tree, false otherwise.
    */
    bool intersect_impl(NonRecursive const & tree, IndexMap & candidates) const
    {
      // Reset statistics
      this->m_check_counter = 0;

      // Collect the original object boxes
      BoxUniquePtrList const & boxes {*this->m_boxes};

      // Return if the tree is empty
      if (this->is_empty() || tree.is_empty()) {return false;}

      // Setup the stack
      this->m_stack.clear();
      this->m_stack.reserve(this->size() + tree.size() + 2);
      this->m_stack.emplace_back(0);
      this->m_stack.emplace_back(0);

      // Main loop that checks the intersection iteratively
      candidates.clear();
      while (!this->m_stack.empty())
      {
        // TODO
      }

      // Return true if the trees intersect
      return !candidates.empty();
    }

    /**
    * Minimum distance between an object and the tree.
    * \param[in] obj Object to compute the distance to.
    * \param[out] candidates Minimum distance candidates.
    * \return The minimum distance between the Object and the tree.
    * \tparam Object Type of the object to compute the distance to.
    * \note Object must have a method \c interior_distance that computes the distance to a box.
    */
    template <typename Object>
    Real distance_impl(Object const & obj, IndexSet & candidates) const
    {
      // Reset statistics
      this->m_check_counter = 0;

      // Collect the original object boxes
      BoxUniquePtrList const & boxes {*this->m_boxes};

      // Return a negative value if the tree is empty
      if (this->is_empty() || this->is_empty()) {return -1.0;}

      // Setup the stack
      this->m_stack.clear();
      this->m_stack.reserve(2*this->size() + 1);
      this->m_stack.emplace_back(0);
      this->m_stack.emplace_back(0);

      // Main loop that checks the intersection iteratively
      Real distance{std::numeric_limits<Real>::max()};
      candidates.clear();
      while (!this->m_stack.empty())
      {
        // Pop the node from stack
        Integer const id{this->m_stack.back()}; this->m_stack.pop_back();

        // Get the node
        AABBnode const & node {this->m_tree_structure[id]};

        // Compute the distance between the object and the bounding box
        ++this->m_check_counter;
        Real tmp_distance{node.box.interior_distance(obj)};

        // If the distance is greater than the temporary minimum distance, skip the node
        if (tmp_distance > distance) {continue;}

        // Compute the distance between the object and the long boxes on the node
        Integer const id_ini{node.box_ptr};
        Integer const id_end{node.box_ptr + node.box_num};
        for (Integer i{id_ini}; i < id_end; ++i) {
          Integer const pos{this->m_tree_boxes_map[i]};
          ++this->m_check_counter;
          tmp_distance = boxes[pos]->interior_distance(obj);
          if (tmp_distance == distance) {candidates.insert(pos);}
          else if (tmp_distance < distance) {
            candidates.clear(); candidates.insert(pos); distance = tmp_distance;
          }
        }

        // Push children on the stack if thay are not leafs
        if (node.child_l > 0) {this->m_stack.emplace_back(node.child_l);}
        if (node.child_r > 0) {this->m_stack.emplace_back(node.child_r);}
      }

      // Return the distance between the point and the tree
      return distance;
    }

    /**
    * Minimum distance between an current tree and another tree.
    * \param[in] tree Tree to compute the distance to.
    * \param[out] candidates Minimum distance candidates.
    * \return The minimum distance between the trees.
    */
    Real distance_impl(NonRecursive const & tree, IndexSet & candidates) const
    {
      // Reset statistics
      this->m_check_counter = 0;

      // Collect the original object boxes
      BoxUniquePtrList const & boxes {*this->m_boxes};

      // Return if the tree is empty
      if (this->is_empty() || tree.is_empty()) {return false;}

      // Setup the stack
      this->m_stack.clear();
      this->m_stack.reserve(this->size() + tree.size() + 2);
      this->m_stack.emplace_back(0);
      this->m_stack.emplace_back(0);

      // Main loop that checks the intersection iteratively
      Real distance{std::numeric_limits<Real>::max()};
      candidates.clear();
      while (!this->m_stack.empty())
      {
        // TODO
      }

      // Return the distance between the trees
      return distance;
    }

  }; // NonRecursive

} // namespace AABBtree

#endif // INCLUDE_AABBTREE_NONRECURSIVE_HXX
