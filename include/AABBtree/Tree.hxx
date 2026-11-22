/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The AABBtree project is distributed under the BSD 2-Clause License.       *
 *                                                                           *
 * Davide Stocco Enrico Bertolazzi                                           *
 * University of Trento University of Trento                                 *
 * davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef INCLUDE_AABBTREE_TREE_HXX
#define INCLUDE_AABBTREE_TREE_HXX

#include "AABBtree/Box.hxx"

namespace AABBtree {

namespace detail {

/**
 * Access an object stored either by value or through a (smart) pointer.
 */
template <typename T> T const &as_ref(T const &x) { return x; }
template <typename T> T const &as_ref(T const *x) { return *x; }
template <typename T> T const &as_ref(T *x) { return *x; }
template <typename T> T const &as_ref(std::unique_ptr<T> const &x) {
  return *x;
}
template <typename T> T const &as_ref(std::shared_ptr<T> const &x) {
  return *x;
}

/**
 * Integer counter that can be safely written by concurrent readers of a const
 * tree. Relaxed atomics: the value is the one written by the last query.
 */
class RelaxedCounter {
  std::atomic<Integer> m_value{0};

public:
  RelaxedCounter() = default;
  RelaxedCounter(RelaxedCounter const &other) : m_value(other.load()) {}
  RelaxedCounter &operator=(RelaxedCounter const &other) {
    this->store(other.load());
    return *this;
  }
  Integer load() const { return m_value.load(std::memory_order_relaxed); }
  void store(Integer const v) { m_value.store(v, std::memory_order_relaxed); }
};

} // namespace detail

/**
 * \brief A class representing a \em non-recursive axis-aligned bounding box
 * tree (AABB tree).
 *
 * The Tree class provides an efficient way to store and query a set of
 * axis-aligned bounding boxes
 * \em non-recursively. It supports various operations such as building the
 * tree, adding axis-aligned boxes, and performing intersection queries. The
 * Tree class is templated on the type of the real numbers and the dimensions of
 * the bounding boxes and is based, as the name suggests, on a novel
 * non-recursive algorithm.
 *
 * \note Thread safety: once the tree is built, all the \c const methods can be
 * called concurrently from multiple threads. The query statistics returned by
 * \c stats() refer to the last completed query (or batch of queries).
 * \note With OpenMP support (\c _OPENMP defined) the \c *_batch methods and
 * the tree-tree queries (\c intersect(Tree), \c self_intersect and
 * \c distance(Tree)) run in parallel. The number of threads is set with
 * \c num_threads (1 disables parallelism). Tree-tree queries run serially
 * when called from inside a parallel region or when the trees have fewer than
 * \c AABBTREE_PARALLEL_MIN_NODES nodes in total. Results do not depend on the
 * number of threads.
 * \tparam Real Type of the scalar coefficients
 * \tparam N Dimension of the ambient space.
 */
template <typename Real, Integer N> class Tree {
public:
  static_assert(std::is_floating_point<Real>::value,
                "Tree Real type must be a floating-point type.");
  static_assert(std::is_integral<Integer>::value &&
                    std::is_signed<Integer>::value,
                "Tree integer type must be a signed integer type.");
  static_assert(N > 0, "Tree dimension must be positive.");

  /**
   * Structure representing the statistics of the AABB tree.
   */
  struct Statistics {
    // Build statistics
    Integer objects{0};    /**< Number of objects in the tree. */
    Integer nodes{0};      /**< Total number of nodes in the tree. */
    Integer leafs{0};      /**< Number of leaf nodes. */
    Integer long_boxes{0}; /**< Number of long boxes (in non-leaf nodes). */
    Integer depth{0};      /**< Maximum depth of the tree. */
    Integer left_nodes{0}; /**< Number of nodes in left subtree. */
    Integer left_leafs{0}; /**< Number of leafs in left subtree. */
    Integer left_long_boxes{0};  /**< Number of long boxes in left subtree. */
    Integer left_depth{0};       /**< Depth of left subtree. */
    Integer right_nodes{0};      /**< Number of nodes in right subtree. */
    Integer right_leafs{0};      /**< Number of leafs in right subtree. */
    Integer right_long_boxes{0}; /**< Number of long boxes in right subtree. */
    Integer right_depth{0};      /**< Depth of right subtree. */
    Integer dump_counter{0};     /**< Number of dumpings in the last build. */
    Real balance_ratio{0.0};     /**< Size ratio of left/right subtrees. */
    Real depth_ratio{0.0};       /**< Depth ratio of left/right subtree. */

    // Query statistics
    Integer check_counter{0};      /**< Collision checks performed. */
    Integer node_check_counter{0}; /**< Node checks performed. */
    Integer long_check_counter{0}; /**< Long box checks performed. */
    Integer objs_check_counter{0}; /**< Object checks performed. */

    /**
     * Reset all statistics to their initial zero state
     *
     * This method efficiently resets all statistics counters to zero and ratios
     * to 0, returning the structure to its default-constructed state.
     *
     * \note This operation is noexcept and has constant time complexity O(1)
     * \post All counters are zero and ratios are 0
     */
    void reset() noexcept { *this = Statistics{}; }

    /**
     * Print the tree info to an output stream.
     * \param[in] os Output stream to print the tree info to.
     */
    void print(OutStream &os) const {
      // Print the tree info
      os << "AABB TREE INFO" << std::endl
         << "\tAmbient dimension : " << N << std::endl
         << "\tBasic type        : " << typeid(Real).name() << std::endl
         << "\tObjects           : " << objects << std::endl
         << "\tNodes             : " << nodes << std::endl
         << "\tLeafs             : " << leafs << std::endl
         << "\tLong boxes        : " << long_boxes << std::endl
         << "\tDepth             : " << depth << std::endl
         << "\tLeft nodes        : " << left_nodes << std::endl
         << "\tLeft leafs        : " << left_leafs << std::endl
         << "\tLeft long boxes   : " << left_long_boxes << std::endl
         << "\tLeft depth        : " << left_depth << std::endl
         << "\tRight nodes       : " << right_nodes << std::endl
         << "\tRight leafs       : " << right_leafs << std::endl
         << "\tRight long boxes  : " << right_long_boxes << std::endl
         << "\tRight depth       : " << right_depth << std::endl
         << "\tDump counter      : " << dump_counter << std::endl
         << "\tBalance ratio     : " << balance_ratio << std::endl
         << "\tDepth ratio       : " << depth_ratio << std::endl
         << "\tCheck counter     : " << check_counter << std::endl;
    }
  };

private:
  // Basic types definitions
  using Box = AABBtree::Box<Real, N>;
  using BoxUniquePtr = AABBtree::BoxUniquePtr<Real, N>;
  using BoxUniquePtrList = AABBtree::BoxUniquePtrList<Real, N>;
  using Vector = AABBtree::Vector<Real, N>;
  using Point = AABBtree::Point<Real, N>;
  using HeapItem = std::pair<Real, Integer>;

  /**
   * Structure representing a node of the AABB tree.
   */
  struct Node {
    Box box;             /**< Bounding box of the subtree. */
    Box box_long;        /**< Bounding box of long boxes. */
    Integer box_ptr;     /**< Pointer to the first box in the reordering map. */
    Integer box_num;     /**< Number of boxes in the subtree. */
    Integer box_tot_num; /**< Total number of boxes in the subtree. */
    Integer parent;      /**< Root node of the subtree. */
    Integer child_l;     /**< Left child of the subtree. */
    Integer child_r;     /**< Right child of the subtree. */
  };

  /**
   * Counters collected during a single query (or a batch of queries).
   */
  struct QueryCounters {
    Integer node{0}; /**< Node checks. */
    Integer lng{0};  /**< Long boxes checks. */
    Integer objs{0}; /**< Object checks. */
    QueryCounters &operator+=(QueryCounters const &o) {
      node += o.node;
      lng += o.lng;
      objs += o.objs;
      return *this;
    }
  };

  /**
   * Per-thread scratch memory used by the queries.
   */
  struct Workspace {
    IndexList stack;              /**< Traversal stack. */
    std::vector<HeapItem> bstack; /**< Traversal stack (lower bound, node). */
    std::vector<HeapItem> heap;   /**< Heap for the closest-n queries. */
    QueryCounters counters;       /**< Query counters. */
    Workspace() {
      stack.reserve(128);
      bstack.reserve(128);
    }
  };

  /**
   * Default distance function: interior distance between an object and a box.
   */
  struct InteriorDistance {
    template <typename Object>
    Real operator()(Object const &o, Box const &b) const {
      return b.interior_distance(o);
    }
  };

  // Tree hierarchy
  std::unique_ptr<BoxUniquePtrList> m_boxes_ptr{
      nullptr}; /**< Vector of unique pointers to the bounding boxes. */
  std::vector<Node> m_tree_structure; /**< Tree structure. */
  IndexList m_tree_boxes_map; /**< Reordering map between the vector of boxes
                                 and the tree internal structure. */

  // Tree parameters
  Integer m_max_dumpings{N};       /**< Maximum number of dimension dumpings. */
  Integer m_max_nodal_objects{50}; /**< Maximum number of objects per node. */
  Real m_separation_ratio_tolerance{0.1}; /**< Bounding boxes separation. */
  Integer m_num_threads{0}; /**< Threads for batch queries (0 = default). */

  // Statistics
  mutable detail::RelaxedCounter m_node_check_counter; /**< Nodes checks. */
  mutable detail::RelaxedCounter m_long_check_counter; /**< Long box checks. */
  mutable detail::RelaxedCounter m_objs_check_counter; /**< Object checks. */
  Integer m_dump_counter{0}; /**< Number of dumpings in the last build. */

  /**
   * Encode a node id to mark that only its long boxes must be considered.
   * The encoding is an involution: 0 --> -1, 1 --> -2, -1 --> 0, ...
   */
  static constexpr Integer negate(Integer const id) { return -1 - id; }

  /**
   * Push the children of a node on a (lower bound, node) stack so that the
   * nearest child is visited first. Children farther than the bound are not
   * pushed at all.
   */
  template <typename Dist>
  void push_children_nearest_first(std::vector<HeapItem> &stack,
                                   Node const &node, Real const bound,
                                   QueryCounters &c, Dist &&dist) const {
    c.node += 2;
    Real const d_l{dist(m_tree_structure[node.child_l].box)};
    Real const d_r{dist(m_tree_structure[node.child_r].box)};
    // The last pushed child is popped (visited) first
    if (d_l <= d_r) {
      if (d_r <= bound) stack.emplace_back(d_r, node.child_r);
      if (d_l <= bound) stack.emplace_back(d_l, node.child_l);
    } else {
      if (d_l <= bound) stack.emplace_back(d_l, node.child_l);
      if (d_r <= bound) stack.emplace_back(d_r, node.child_r);
    }
  }

  /**
   * Store the counters of the last query.
   */
  void store_counters(QueryCounters const &c) const {
    m_node_check_counter.store(c.node);
    m_long_check_counter.store(c.lng);
    m_objs_check_counter.store(c.objs);
  }

public:
  /**
   * Class destructor for the \em non-recursive axis-aligned bounding box tree.
   */
  ~Tree() = default;

  /**
   * Class constructor for the \em non-recursive axis-aligned bounding box tree.
   */
  Tree() = default;

  /**
   * Get a look at the vector of unique pointers to the bounding boxes.
   * \return A const reference to the vector of unique pointers to the bounding
   * boxes.
   */
  BoxUniquePtrList const &boxes() const { return *m_boxes_ptr; }

  /**
   * Get a look at the i-th unique pointer to the bounding box.
   * \param[in] i Index of the unique pointer to the bounding box.
   * \return A const reference to the i-th unique pointer to the bounding box.
   */
  BoxUniquePtr const &box(Integer const i) const { return (*m_boxes_ptr)[i]; }

  /**
   * Set the maximum number of dimension dumpings.
   * \param[in] max_dumpings Maximum number of dimension dumpings.
   */
  void max_dumpings(Integer const max_dumpings) {
    constexpr char CMD[]{"AABBtree::Tree::max_dumpings(...): "};
    AABBTREE_ASSERT(max_dumpings > 0 && max_dumpings <= N,
                    CMD << "input must be a positive integer in the range [1, "
                        << N << "].");
    m_max_dumpings = max_dumpings;
  }

  /**
   * Get the maximum number of dimension dumpings.
   * \return The maximum number of dimension dumpings.
   */
  Integer max_dumpings() const { return m_max_dumpings; }

  /**
   * Set the maximum number of objects per node.
   * \param[in] max_nodal_objects Maximum number of objects per node.
   */
  void max_nodal_objects(Integer const max_nodal_objects) {
    constexpr char CMD[]{"AABBtree::Tree::max_nodal_objects(...): "};
    AABBTREE_ASSERT(max_nodal_objects > 0,
                    CMD << "input must be a positive integer.");
    m_max_nodal_objects = max_nodal_objects;
  }

  /**
   * Get the maximum number of objects per node.
   * \return The maximum number of objects per node.
   */
  Integer max_nodal_objects() const { return m_max_nodal_objects; }

  /**
   * Set the balance ratio tolerance for bounding boxes.
   * \param[in] ratio Balance ratio tolerance for bounding boxes.
   */
  void separation_ratio_tolerance(Real const ratio) {
    constexpr char CMD[]{"AABBtree::Tree::separation_ratio_tolerance(...): "};
    AABBTREE_ASSERT(ratio > 0.0 && ratio < 1.0,
                    CMD << "input must be in the open range (0, 1).");
    m_separation_ratio_tolerance = ratio;
  }

  /**
   * Get the balance ratio tolerance for bounding boxes.
   * \return The balance ratio tolerance for bounding boxes.
   */
  Real separation_ratio_tolerance() const {
    return m_separation_ratio_tolerance;
  }

  /**
   * Set the number of threads used by the \c *_batch queries.
   * \param[in] num_threads Number of threads (0 means the OpenMP default).
   * \note Without OpenMP support the batch queries always run serially.
   */
  void num_threads(Integer const num_threads) {
    constexpr char CMD[]{"AABBtree::Tree::num_threads(...): "};
    AABBTREE_ASSERT(num_threads >= 0,
                    CMD << "input must be a non-negative integer.");
    m_num_threads = num_threads;
  }

  /**
   * Get the number of threads used by the \c *_batch queries.
   * \return The number of threads (0 means the OpenMP default).
   */
  Integer num_threads() const { return m_num_threads; }

  /**
   * Get the tree structure.
   * \return The tree structure.
   */
  std::vector<Node> const &structure() const { return m_tree_structure; }

  /**
   * Get the i-th tree node.
   * \param[in] i Index of the node.
   * \return The i-th tree node.
   */
  Node const &node(Integer const i) const { return m_tree_structure[i]; }

  /**
   * Get the number of nodes in the tree.
   * \return The number of nodes in the tree.
   */
  Integer size() const { return static_cast<Integer>(m_tree_structure.size()); }

  /**
   * Check if tree is empty.
   * \return True if the tree is empty, false otherwise.
   */
  bool is_empty() const { return m_tree_structure.empty(); }

  /**
   * Clear the tree.
   */
  void clear() {
    m_tree_structure.clear();
    m_tree_boxes_map.clear();
  }

  /**
   * Build the tree given the bounding boxes.
   * \param[in] boxes_ptr Bounding boxes to build the tree from.
   */
  void build(std::unique_ptr<BoxUniquePtrList> boxes_ptr) {
    this->clear();
    m_boxes_ptr = std::move(boxes_ptr);
    this->build();
  }

  /**
   * Build the tree from the bounding boxes already owned by the tree.
   * \warning This method requires the bounding boxes to be already set (see
   * \c build(std::unique_ptr<BoxUniquePtrList>)).
   */
  void build() {

    // Check if the boxes pointer is valid
    AABBTREE_ASSERT(m_boxes_ptr != nullptr,
                    "AABBtree::Tree::build(...): The boxes pointer is "
                    "not valid. Please provide a valid pointer to the "
                    "bounding boxes before building the tree.");

    BoxUniquePtrList &boxes{*m_boxes_ptr};

    // Clear tree structure
    Integer const num_boxes{static_cast<Integer>(boxes.size())};
    m_tree_structure.clear();
    m_tree_boxes_map.clear();
    m_dump_counter = 0;

    // An empty set of boxes produces an empty tree
    if (num_boxes == 0) {
      return;
    }
    m_tree_structure.reserve(2 * num_boxes + 1);

    // Setup the root node
    Node root;
    root.parent = -1;
    root.child_l = -1;
    root.child_r = -1;
    root.box_ptr = 0;
    root.box_num = 0;
    root.box_tot_num = num_boxes;

    // Setup the boxes map and compute the root box
    root.box.set_empty();
    m_tree_boxes_map.reserve(num_boxes);
    for (BoxUniquePtr const &box : boxes) {
      root.box.extend(*box);
      m_tree_boxes_map.emplace_back(root.box_num++);
    }
    // The root is a leaf until it is split: all its boxes are "long" boxes
    root.box_long = root.box;
    m_tree_structure.emplace_back(root);

    // Setup the stack
    IndexList stack;
    stack.reserve(128);
    stack.emplace_back(0);

    // Main loop that divide the nodes iteratively until all constraints
    // satisfied
    IndexList best_map(m_tree_boxes_map.size());
    while (!stack.empty()) {
      // Pop the node from stack
      Integer const id{stack.back()};
      stack.pop_back();

      // Get the node
      Node &node{m_tree_structure[id]};

      // If the node has less than one object, skip it
      if (node.box_num <= 1) {
        continue;
      }

      // Compute the separation line and tolerance
      Vector sizes{node.box.max() - node.box.min()};
      Integer sorting[N];
      std::iota(sorting, sorting + N, 0);
      auto compare = [&sizes](Integer i, Integer j) -> bool {
        return sizes[i] > sizes[j];
      };
      std::sort(sorting, sorting + N, compare);

      Integer axis, n_long, n_left, n_right, id_ini, id_end;
      Real separation_line, separation_tolerance;

      Integer n_long_saved{node.box_num + 1};
      Integer n_diff_saved{node.box_num + 1};
      Integer axis_saved{sorting[0]};
      bool valid_axis_found{false};
      for (Integer dump{0}; dump < m_max_dumpings; ++dump) {
        axis = sorting[dump];
        separation_line = node.box.baricenter(axis);
        separation_tolerance = sizes[axis] * m_separation_ratio_tolerance;

        // Separate short and long boxes
        n_long = n_left = n_right = 0;
        id_ini = node.box_ptr;
        id_end = node.box_ptr + node.box_num;
        while (id_ini < id_end) {
          Box const &box_id{*boxes[m_tree_boxes_map[id_ini]]};
          typename Box::Side const side{
              box_id.which_side(separation_line, separation_tolerance, axis)};
          switch (side) {
          case Box::Side::LEFT: // Left boxes are moved to the end
            ++n_left;
            --id_end;
            std::swap(m_tree_boxes_map[id_ini], m_tree_boxes_map[id_end]);
            break;
          case Box::Side::RIGHT: // Right boxes are moved to the end
            ++n_right;
            --id_end;
            std::swap(m_tree_boxes_map[id_ini], m_tree_boxes_map[id_end]);
            break;
          default:
            ++n_long;
            ++id_ini;
          }
        }

        // Save the first and the best solution found so far
        Integer n_diff{std::abs(n_left - n_right)};
        if ((n_long <= n_long_saved && n_long <= m_max_nodal_objects &&
             n_diff < n_diff_saved) ||
            !valid_axis_found) {
          n_long_saved = n_long;
          n_diff_saved = n_diff;
          axis_saved = axis;
          std::copy_n(m_tree_boxes_map.data() + node.box_ptr, node.box_num,
                      best_map.data() + node.box_ptr);
          valid_axis_found = true;
        }

        ++m_dump_counter;
      }

      // Check if a valid solution has been found
      AABBTREE_ASSERT_WARNING(
          valid_axis_found,
          "AABBtree::Tree::build(...): Failed to find a valid "
          "splitting axis.");

      // Terminate splitting if the number of long boxes if cannot be split
      if (n_long_saved == node.box_num) {
        continue;
      }

      // Use the best solution found
      std::copy_n(best_map.data() + node.box_ptr, node.box_num,
                  m_tree_boxes_map.data() + node.box_ptr);
      n_long = n_long_saved;
      axis = axis_saved;
      separation_line = node.box.baricenter(axis);
      separation_tolerance = sizes[axis] * m_separation_ratio_tolerance;

      // Separate the left and right boxes
      n_left = n_right = 0;
      id_ini = node.box_ptr + n_long;
      id_end = node.box_ptr + node.box_num;
      while (id_ini < id_end) {
        Box const &box_id{*boxes[m_tree_boxes_map[id_ini]]};
        typename Box::Side const side{
            box_id.which_side(separation_line, separation_tolerance, axis)};
        switch (side) {
        case Box::Side::LEFT: // Left boxes stay in place
          ++id_ini;
          ++n_left;
          break;
        case Box::Side::RIGHT: // Right boxes are swapped with the last one
          --id_end;
          ++n_right;
          std::swap(m_tree_boxes_map[id_ini], m_tree_boxes_map[id_end]);
          break;
        default:
          break;
        }
      }

      // Set the left and right children indexes
      node.child_l = static_cast<Integer>(this->size() + 0);
      node.child_r = static_cast<Integer>(this->size() + 1);

      // Finalize the root node setup (left and right children)
      Node node_l;
      node_l.parent = id;  // Current node
      node_l.child_l = -1; // Default (leaf)
      node_l.child_r = -1; // Default (leaf)
      node_l.box_num = n_left;
      node_l.box_tot_num = n_left;

      Node node_r;
      node_r.parent = id;  // Current node
      node_r.child_l = -1; // Default (leaf)
      node_r.child_r = -1; // Default (leaf)
      node_r.box_num = n_right;
      node_r.box_tot_num = n_right;

      // Compute the bounding box of the long boxes, and left and right
      // children
      Integer j{node.box_ptr};

      node.box_num = n_long;
      node.box_long.set_empty();
      for (Integer i{0}; i < n_long; ++i) {
        node.box_long.extend(*boxes[m_tree_boxes_map[j++]]);
      }

      node_l.box.set_empty();
      node_l.box_ptr = j;
      for (Integer i{0}; i < n_left; ++i) {
        node_l.box.extend(*boxes[m_tree_boxes_map[j++]]);
      }
      node_l.box_long = node_l.box;

      node_r.box.set_empty();
      node_r.box_ptr = j;
      for (Integer i{0}; i < n_right; ++i) {
        node_r.box.extend(*boxes[m_tree_boxes_map[j++]]);
      }
      node_r.box_long = node_r.box;

      // Push children on stack (before emplacing, as emplacing may invalidate
      // the reference to the current node)
      Integer const child_l{node.child_l};
      Integer const child_r{node.child_r};

      // Push nodes on tree structure
      m_tree_structure.emplace_back(node_l);
      m_tree_structure.emplace_back(node_r);

      // Push children on stack
      stack.emplace_back(child_l);
      stack.emplace_back(child_r);
    }
  }

  // ---------------------------------------------------------------------------
  // Single queries
  // ---------------------------------------------------------------------------

  /**
   * Intersect the tree with an object.
   * \param[in] obj Object to intersect with.
   * \param[out] candidates Intersection result (boxes indexes).
   * \return True if the object intersects the tree, false otherwise.
   * \tparam Object Type of the object to intersect with.
   * \note Object must have a method \c intersects that computes the
   * intersection with a box.
   */
  template <typename Object>
  bool intersect(Object const &obj, IndexSet &candidates) const {
    Workspace ws;
    this->intersect_impl(obj, candidates, ws);
    this->store_counters(ws.counters);
    return !candidates.empty();
  }

  /**
   * Intersect the tree with another tree.
   * \param[in] tree Tree to intersect with.
   * \param[out] candidates Intersection result (boxes indexes).
   * \return True if the trees intersect, false otherwise.
   */
  bool intersect(Tree const &tree, IndexMap &candidates) const {
    Workspace ws;
    this->intersect_tree_impl(tree, candidates, ws, false);
    this->store_counters(ws.counters);
    return !candidates.empty();
  }

  /**
   * Intersect the tree with itself.
   * \param[out] candidates Intersection result: for each box, the set of the
   * other boxes it intersects (a box is never reported with itself).
   * \return True if at least two boxes intersect, false otherwise.
   */
  bool self_intersect(IndexMap &candidates) const {
    Workspace ws;
    this->intersect_tree_impl(*this, candidates, ws, true);
    this->store_counters(ws.counters);
    return !candidates.empty();
  }

  /**
   * Self-intersect the tree (i.e., intersect the tree with itself to find all
   * the intersecting boxes).
   * \param[out] candidates Intersection result: all the boxes that intersect at
   * least another box.
   * \return True if at least two boxes intersect, false otherwise.
   */
  bool self_intersect(IndexSet &candidates) const {
    IndexMap candidates_map;
    bool const intersects{this->self_intersect(candidates_map)};
    candidates.clear();
    for (auto const &[key, values] : candidates_map) {
      candidates.emplace(key);
      candidates.insert(values.begin(), values.end());
    }
    return intersects;
  }

  /**
   * Compute the minimum distance between an object and the tree.
   * \param[in] obj Object to compute the distance to.
   * \param[out] candidates Minimum distance candidates.
   * \return The minimum distance between the object and the tree, negative if
   * the tree is empty.
   * \tparam Object Type of the object to compute the distance to.
   * \note Object must have a method \c interior_distance that computes the
   * distance to a box.
   */
  template <typename Object>
  Real distance(Object const &obj, IndexSet &candidates) const {
    Workspace ws;
    Real const d{this->distance_impl(obj, candidates, ws)};
    this->store_counters(ws.counters);
    return d;
  }

  /**
   * Compute the minimum distance between the current tree and another tree.
   * \param[in] tree Tree to compute the distance to.
   * \param[out] candidates Minimum distance candidates.
   * \return The minimum distance between the trees, negative if one of the
   * trees is empty.
   */
  Real distance(Tree const &tree, IndexMap &candidates) const {
    Workspace ws;
    Real const d{this->distance_tree_impl(tree, candidates, ws)};
    this->store_counters(ws.counters);
    return d;
  }

  /**
   * Find the first \f$ n \f$ candidates from an object, using a custom
   * distance function.
   * \param[in] obj Object to compute the distance to.
   * \param[in] n Number of candidates to find.
   * \param[out] candidates First \f$ n \f$ candidates.
   * \param[in] distance_function Custom distance function (default is the
   * interior distance).
   * \return The maximum distance between the candidates and the object. The
   * return value is negative if no candidates are found.
   * \tparam Object Type of the object to compute the distance to.
   * \tparam Function Type of the custom distance function.
   * \note The distance function must be a lower bound with respect to box
   * inclusion, i.e. if box \f$ A \f$ contains box \f$ B \f$ then
   * \f$ d(o, A) \leq d(o, B) \f$. This is used to prune the tree.
   */
  template <typename Object, typename Function = InteriorDistance>
  Real closest(Object const &obj, Integer const n, IndexSet &candidates,
               Function distance_function = Function{}) const {
    Workspace ws;
    Real const d{
        this->closest_impl(obj, n, candidates, distance_function, ws)};
    this->store_counters(ws.counters);
    return d;
  }

  /**
   * Find the candidates that are within a given distance from a object, using
   * a custom distance function.
   * \param[in] obj Object to compute the distance to.
   * \param[in] max_distance Maximum distance to consider.
   * \param[out] candidates Minimum distance candidates.
   * \param[in] distance_function Custom distance function (default is the
   * interior distance).
   * \return True if at least one object is within the given distance, false
   * otherwise.
   * \tparam Object Type of the object to compute the distance to.
   * \tparam Function Type of the custom distance function.
   * \note The distance function must be a lower bound with respect to box
   * inclusion (see \c closest).
   */
  template <typename Object, typename Function = InteriorDistance>
  bool within_distance(Object const &obj, Real const max_distance,
                       IndexSet &candidates,
                       Function distance_function = Function{}) const {
    Workspace ws;
    this->within_distance_impl(obj, max_distance, candidates,
                               distance_function, ws);
    this->store_counters(ws.counters);
    return !candidates.empty();
  }

  // ---------------------------------------------------------------------------
  // Batch queries (parallel with OpenMP)
  // ---------------------------------------------------------------------------

  /**
   * Intersect the tree with a batch of objects.
   * \param[in] objs Random-access container of objects (stored by value or by
   * pointer/smart pointer).
   * \param[out] candidates For each object, the intersection result.
   * \return True if at least one object intersects the tree.
   * \tparam Container Type of the container of objects.
   */
  template <typename Container>
  bool intersect_batch(Container const &objs,
                       std::vector<IndexSet> &candidates) const {
    Integer const n{static_cast<Integer>(objs.size())};
    candidates.resize(n);
    this->store_counters(this->parallel_for(n, [&](Integer i, Workspace &ws) {
      this->intersect_impl(detail::as_ref(objs[i]), candidates[i], ws);
    }));
    return std::any_of(candidates.begin(), candidates.end(),
                       [](IndexSet const &c) { return !c.empty(); });
  }

  /**
   * Compute the minimum distance between the tree and a batch of objects.
   * \param[in] objs Random-access container of objects.
   * \param[out] distances For each object, the minimum distance.
   * \param[out] candidates For each object, the minimum distance candidates.
   * \tparam Container Type of the container of objects.
   */
  template <typename Container>
  void distance_batch(Container const &objs, std::vector<Real> &distances,
                      std::vector<IndexSet> &candidates) const {
    Integer const n{static_cast<Integer>(objs.size())};
    distances.resize(n);
    candidates.resize(n);
    this->store_counters(this->parallel_for(n, [&](Integer i, Workspace &ws) {
      distances[i] =
          this->distance_impl(detail::as_ref(objs[i]), candidates[i], ws);
    }));
  }

  /**
   * Find the first \f$ n \f$ candidates for a batch of objects.
   * \param[in] objs Random-access container of objects.
   * \param[in] n Number of candidates to find for each object.
   * \param[out] distances For each object, the maximum candidate distance.
   * \param[out] candidates For each object, the first \f$ n \f$ candidates.
   * \param[in] distance_function Custom distance function (see \c closest).
   * \tparam Container Type of the container of objects.
   * \tparam Function Type of the custom distance function.
   */
  template <typename Container, typename Function = InteriorDistance>
  void closest_batch(Container const &objs, Integer const n,
                     std::vector<Real> &distances,
                     std::vector<IndexSet> &candidates,
                     Function distance_function = Function{}) const {
    Integer const m{static_cast<Integer>(objs.size())};
    distances.resize(m);
    candidates.resize(m);
    this->store_counters(this->parallel_for(m, [&](Integer i, Workspace &ws) {
      distances[i] = this->closest_impl(detail::as_ref(objs[i]), n,
                                        candidates[i], distance_function, ws);
    }));
  }

  /**
   * Find the candidates within a given distance for a batch of objects.
   * \param[in] objs Random-access container of objects.
   * \param[in] max_distance Maximum distance to consider.
   * \param[out] candidates For each object, the candidates within distance.
   * \param[in] distance_function Custom distance function (see \c closest).
   * \return True if at least one candidate is found for at least one object.
   * \tparam Container Type of the container of objects.
   * \tparam Function Type of the custom distance function.
   */
  template <typename Container, typename Function = InteriorDistance>
  bool within_distance_batch(Container const &objs, Real const max_distance,
                             std::vector<IndexSet> &candidates,
                             Function distance_function = Function{}) const {
    Integer const n{static_cast<Integer>(objs.size())};
    candidates.resize(n);
    this->store_counters(this->parallel_for(n, [&](Integer i, Workspace &ws) {
      this->within_distance_impl(detail::as_ref(objs[i]), max_distance,
                                 candidates[i], distance_function, ws);
    }));
    return std::any_of(candidates.begin(), candidates.end(),
                       [](IndexSet const &c) { return !c.empty(); });
  }

  // ---------------------------------------------------------------------------
  // Tree information
  // ---------------------------------------------------------------------------

  /**
   * Compute the depth of the tree from a given node.
   * \param[in] i Index of the node to compute the depth from.
   * \param[out] d Depth of the tree from the given node.
   */
  void depth(Integer const i, Integer &d) const {
    d = 0;
    if (i < 0 || i >= this->size()) {
      return;
    }
    std::vector<std::pair<Integer, Integer>> stack; // (node, depth)
    stack.reserve(128);
    stack.emplace_back(i, 0);
    while (!stack.empty()) {
      auto const [id, level] = stack.back();
      stack.pop_back();
      Node const &node{m_tree_structure[id]};
      if (node.child_l < 0) {
        d = std::max(d, level);
      } else {
        stack.emplace_back(node.child_l, level + 1);
        stack.emplace_back(node.child_r, level + 1);
      }
    }
  }

  /**
   * Compute the number of leafs, nodes, and long boxes of the tree from a
   * given node.
   * \param[in] i Index of the node to start the computation from.
   * \param[out] l Number of leafs of the tree from the given node.
   * \param[out] n Number of nodes of the tree from the given node.
   * \param[out] b Number of long boxes (boxes stored in non-leaf nodes) of the
   * tree from the given node.
   */
  void nodes(Integer const i, Integer &l, Integer &n, Integer &b) const {
    l = n = b = 0;
    if (i < 0 || i >= this->size()) {
      return;
    }
    IndexList stack;
    stack.reserve(128);
    stack.emplace_back(i);
    while (!stack.empty()) {
      Integer const id{stack.back()};
      stack.pop_back();
      Node const &node{m_tree_structure[id]};
      ++n;
      if (node.child_l < 0) {
        ++l;
      } else {
        b += node.box_num;
        stack.emplace_back(node.child_l);
        stack.emplace_back(node.child_r);
      }
    }
  }

  /**
   * Compute some statistics about the current tree.
   * \param[out] stats Statistics about the tree.
   */
  void stats(Statistics &stats) const {
    // Reset statistics
    stats.reset();

    // Compute/copy the build statistics
    stats.objects =
        m_boxes_ptr ? static_cast<Integer>(m_boxes_ptr->size()) : 0;
    stats.dump_counter = m_dump_counter;
    if (!this->is_empty()) {
      Node const &root{m_tree_structure[0]};
      this->nodes(0, stats.leafs, stats.nodes, stats.long_boxes);
      this->depth(0, stats.depth);
      this->nodes(root.child_l, stats.left_leafs, stats.left_nodes,
                  stats.left_long_boxes);
      this->depth(root.child_l, stats.left_depth);
      this->nodes(root.child_r, stats.right_leafs, stats.right_nodes,
                  stats.right_long_boxes);
      this->depth(root.child_r, stats.right_depth);
      auto ratio = [](Integer a, Integer b) -> Real {
        return b != 0 ? static_cast<Real>(a) / static_cast<Real>(b)
                      : static_cast<Real>(0.0);
      };
      stats.balance_ratio = ratio(stats.left_leafs, stats.right_leafs);
      stats.depth_ratio = ratio(stats.left_depth, stats.right_depth);
    }

    // Copy the check counter
    stats.node_check_counter = m_node_check_counter.load();
    stats.long_check_counter = m_long_check_counter.load();
    stats.objs_check_counter = m_objs_check_counter.load();
    stats.check_counter = stats.node_check_counter + stats.long_check_counter +
                          stats.objs_check_counter;
  }

  /**
   * Print the tree info to an output stream.
   * \param[in] os Output stream to print the tree info to.
   */
  void print(OutStream &os) const {
    // Retrieve the statistics
    Statistics stats;
    this->stats(stats);
    // Print the tree info
    stats.print(os);
  }

private:
  // ---------------------------------------------------------------------------
  // Parallel driver
  // ---------------------------------------------------------------------------

  /**
   * Run \c body(i, workspace) for \c i in \f$ [0, n) \f$, in parallel when
   * OpenMP is available. Each thread owns its workspace. The query counters of
   * all the threads are summed and stored as the statistics of the batch.
   * Exceptions thrown by \c body are rethrown on the calling thread.
   */
  template <typename Body>
  QueryCounters parallel_for(Integer const n, Body &&body) const {
    QueryCounters total;
#ifdef _OPENMP
    std::exception_ptr error{nullptr};
    std::atomic<bool> failed{false};
    int const num_threads{this->openmp_threads()};
#pragma omp parallel num_threads(num_threads) if (n > 1 && num_threads > 1)
    {
      Workspace ws;
#pragma omp for schedule(dynamic, 8) nowait
      for (Integer i = 0; i < n; ++i) {
        if (failed.load(std::memory_order_relaxed)) {
          continue;
        }
        try {
          body(i, ws);
        } catch (...) {
#pragma omp critical(aabbtree_parallel_for_error)
          {
            if (!error) {
              error = std::current_exception();
            }
          }
          failed.store(true, std::memory_order_relaxed);
        }
      }
#pragma omp critical(aabbtree_parallel_for_counters)
      total += ws.counters;
    }
    if (error) {
      std::rethrow_exception(error);
    }
#else
    Workspace ws;
    for (Integer i{0}; i < n; ++i) {
      body(i, ws);
    }
    total = ws.counters;
#endif
    return total;
  }

#ifdef _OPENMP
  /**
   * Number of OpenMP threads to use.
   */
  int openmp_threads() const {
    return m_num_threads > 0 ? static_cast<int>(m_num_threads)
                             : omp_get_max_threads();
  }
#endif

  /**
   * Merge the index map \c src into \c dst (\c src is left in a valid but
   * unspecified state). Nodes are moved, not copied.
   */
  static void merge_index_maps(IndexMap &dst, IndexMap &src) {
    dst.merge(src);
    for (auto &[key, values] : src) { // keys already present in dst
      dst[key].merge(values);
    }
  }

  /**
   * Traverse pairs of nodes of this tree and of another tree, starting from the
   * pair of roots. \c step(id_s1, id_s2, result, workspace, stack) processes a
   * pair and pushes the next pairs on the stack.
   *
   * With OpenMP the traversal is first expanded breadth-first until enough
   * independent pairs are available, then the pairs are processed in parallel,
   * each one with its own partial result. Partial results are merged in a
   * fixed order with \c merge(dst, src), so the result is deterministic.
   */
  template <typename Result, typename Step, typename Merge>
  void dual_traversal(Tree const &tree, Result &result, Workspace &ws,
                      Step &&step, Merge &&merge) const {
    IndexList &stack{ws.stack};
    stack.clear();
    stack.emplace_back(0);
    stack.emplace_back(0);
#ifdef _OPENMP
    int const num_threads{this->openmp_threads()};
    if (num_threads > 1 && !omp_in_parallel() &&
        this->size() + tree.size() >= AABBTREE_PARALLEL_MIN_NODES) {
      // Breadth-first expansion until there are enough pairs to balance load
      std::size_t const target{static_cast<std::size_t>(64 * num_threads)};
      IndexList next;
      for (Integer round{0};
           round < 64 && !stack.empty() && stack.size() < target; ++round) {
        next.clear();
        for (std::size_t k{0}; k < stack.size(); k += 2) {
          step(stack[k], stack[k + 1], result, ws, next);
        }
        stack.swap(next);
      }
      if (stack.empty()) {
        return;
      }

      // Depth-first traversal of each pair in parallel
      IndexList const frontier{stack};
      stack.clear();
      Integer const n_pairs{static_cast<Integer>(frontier.size() / 2)};
      std::vector<Result> partial(n_pairs);
      ws.counters += this->parallel_for(n_pairs, [&](Integer i, Workspace &w) {
        IndexList &s{w.stack};
        s.clear();
        s.emplace_back(frontier[2 * i]);
        s.emplace_back(frontier[2 * i + 1]);
        while (!s.empty()) {
          Integer const id_s2{s.back()};
          s.pop_back();
          Integer const id_s1{s.back()};
          s.pop_back();
          step(id_s1, id_s2, partial[i], w, s);
        }
      });
      for (Result &p : partial) {
        merge(result, p);
      }
      return;
    }
#else
    (void)tree;
    (void)merge;
#endif
    while (!stack.empty()) {
      Integer const id_s2{stack.back()};
      stack.pop_back();
      Integer const id_s1{stack.back()};
      stack.pop_back();
      step(id_s1, id_s2, result, ws, stack);
    }
  }

  // ---------------------------------------------------------------------------
  // Query implementations (re-entrant: all the state lives in the workspace)
  // ---------------------------------------------------------------------------

  /**
   * Intersect the tree with an object (see \c intersect).
   */
  template <typename Object>
  void intersect_impl(Object const &obj, IndexSet &candidates,
                      Workspace &ws) const {
    candidates.clear();
    if (this->is_empty()) {
      return;
    }
    BoxUniquePtrList const &boxes{*m_boxes_ptr};
    QueryCounters &c{ws.counters};
    IndexList &stack{ws.stack};
    stack.clear();
    stack.emplace_back(0);
    while (!stack.empty()) {
      Integer const id{stack.back()};
      stack.pop_back();
      Node const &node{m_tree_structure[id]};

      // If the object do not intersects the box, skip the node
      if (++c.node; !node.box.intersects(obj)) {
        continue;
      }

      // Intersect the object with the long boxes on the node. If it is a leaf
      // the long boxes are the boxes of the leaf (and box_long == box).
      bool const leaf{node.child_l < 0};
      if (node.box_num > 0 && (leaf || (++c.lng, node.box_long.intersects(obj)))) {
        Integer const id_end{node.box_ptr + node.box_num};
        for (Integer i{node.box_ptr}; i < id_end; ++i) {
          Integer const pos{m_tree_boxes_map[i]};
          if (++c.objs; boxes[pos]->intersects(obj)) {
            candidates.insert(pos);
          }
        }
      }

      // Push children on the stack
      if (!leaf) {
        stack.emplace_back(node.child_l);
        stack.emplace_back(node.child_r);
      }
    }
  }

  /**
   * Push on the stack the pairs of nodes to be visited after the pair
   * (id_s1, id_s2) in a tree-tree traversal. A negative id means that only the
   * long boxes of the node must be considered.
   */
  static void push_children_pairs(IndexList &stack, Integer const id_s1,
                                  Integer const id_s2, Node const &node_1,
                                  Node const &node_2) {
    Integer const id_1{id_s1 < 0 ? negate(id_s1) : id_s1};
    Integer const id_2{id_s2 < 0 ? negate(id_s2) : id_s2};
    auto push = [&stack](Integer const a, Integer const b) {
      stack.emplace_back(a);
      stack.emplace_back(b);
    };

    // Check if the nodes are leafs
    bool const leaf_1{node_1.child_l < 0};
    bool const leaf_2{node_2.child_l < 0};
    if (leaf_1 && leaf_2) {
      return; // We are done on this branch
    }

    if (leaf_1) {
      push(id_1, node_2.child_l);
      push(id_1, node_2.child_r);
    } else if (leaf_2) {
      push(node_1.child_l, id_2);
      push(node_1.child_r, id_2);
    } else if (node_1.box_tot_num > node_2.box_tot_num) {
      // Split first the larger tree
      if (id_s1 >= 0) {
        push(node_1.child_l, id_s2);
        push(node_1.child_r, id_s2);
        if (node_1.box_num > 0) {
          push(negate(id_1), id_2);
        }
      } else if (id_s2 >= 0) { // and id_s1 < 0
        push(id_s1, node_2.child_l);
        push(id_s1, node_2.child_r);
      }
    } else {
      if (id_s2 >= 0) {
        push(id_s1, node_2.child_l);
        push(id_s1, node_2.child_r);
        if (node_2.box_num > 0) {
          push(id_1, negate(id_2));
        }
      } else if (id_s1 >= 0) { // and id_s2 < 0
        push(node_1.child_l, id_s2);
        push(node_1.child_r, id_s2);
      }
    }
  }

  /**
   * Process one pair of nodes of a tree-tree intersection: intersect their
   * long boxes and push the pairs of children to visit next on the stack.
   * If \c self is true the pairs of a box with itself are discarded.
   */
  void intersect_pair_step(Tree const &tree, Integer const id_s1,
                           Integer const id_s2, IndexMap &candidates,
                           QueryCounters &c, bool const self,
                           IndexList &stack) const {
    BoxUniquePtrList const &boxes_1{*m_boxes_ptr};
    BoxUniquePtrList const &boxes_2{*tree.m_boxes_ptr};
    Integer const id_1{id_s1 < 0 ? negate(id_s1) : id_s1};
    Integer const id_2{id_s2 < 0 ? negate(id_s2) : id_s2};
    Node const &node_1{m_tree_structure[id_1]};
    Node const &node_2{tree.m_tree_structure[id_2]};

    // If the boxes are not intersecting, skip the nodes
    if (++c.node; !node_1.box.intersects(node_2.box)) {
      return;
    }

    // Intersect the long boxes on the nodes
    if (node_1.box_num > 0 && node_2.box_num > 0) {
      if (++c.lng; node_1.box_long.intersects(node_2.box_long)) {
        Integer const id_1_ini{node_1.box_ptr};
        Integer const id_1_end{node_1.box_ptr + node_1.box_num};
        Integer const id_2_ini{node_2.box_ptr};
        Integer const id_2_end{node_2.box_ptr + node_2.box_num};
        if (node_1.box_num < node_2.box_num) {
          for (Integer i{id_1_ini}; i < id_1_end; ++i) {
            Integer const pos_1{m_tree_boxes_map[i]};
            Box const &box_1{*boxes_1[pos_1]};
            if (++c.objs; !box_1.intersects(node_2.box_long)) {
              continue;
            }
            for (Integer j{id_2_ini}; j < id_2_end; ++j) {
              Integer const pos_2{tree.m_tree_boxes_map[j]};
              if (self && pos_1 == pos_2) {
                continue;
              }
              if (++c.objs; box_1.intersects(*boxes_2[pos_2])) {
                candidates[pos_1].insert(pos_2);
              }
            }
          }
        } else {
          for (Integer j{id_2_ini}; j < id_2_end; ++j) {
            Integer const pos_2{tree.m_tree_boxes_map[j]};
            Box const &box_2{*boxes_2[pos_2]};
            if (++c.objs; !box_2.intersects(node_1.box_long)) {
              continue;
            }
            for (Integer i{id_1_ini}; i < id_1_end; ++i) {
              Integer const pos_1{m_tree_boxes_map[i]};
              if (self && pos_1 == pos_2) {
                continue;
              }
              if (++c.objs; boxes_1[pos_1]->intersects(box_2)) {
                candidates[pos_1].insert(pos_2);
              }
            }
          }
        }
      }
    }

    push_children_pairs(stack, id_s1, id_s2, node_1, node_2);
  }

  /**
   * Intersect the tree with another tree (see \c intersect and
   * \c self_intersect).
   */
  void intersect_tree_impl(Tree const &tree, IndexMap &candidates,
                           Workspace &ws, bool const self) const {
    candidates.clear();
    if (this->is_empty() || tree.is_empty()) {
      return;
    }
    this->dual_traversal(
        tree, candidates, ws,
        [this, &tree, self](Integer const id_s1, Integer const id_s2,
                            IndexMap &result, Workspace &w, IndexList &stack) {
          this->intersect_pair_step(tree, id_s1, id_s2, result, w.counters,
                                    self, stack);
        },
        [](IndexMap &dst, IndexMap &src) { merge_index_maps(dst, src); });
  }

  /**
   * Compute the minimum distance between an object and the tree (see
   * \c distance).
   */
  template <typename Object>
  Real distance_impl(Object const &obj, IndexSet &candidates,
                     Workspace &ws) const {
    candidates.clear();
    if (this->is_empty()) {
      return static_cast<Real>(-1.0);
    }
    BoxUniquePtrList const &boxes{*m_boxes_ptr};
    QueryCounters &c{ws.counters};
    auto dist = [&obj](Box const &b) { return b.interior_distance(obj); };

    // Depth-first traversal visiting the nearest child first, so that the
    // pruning bound tightens quickly
    std::vector<HeapItem> &stack{ws.bstack};
    stack.clear();
    ++c.node;
    stack.emplace_back(dist(m_tree_structure[0].box), 0);
    Real distance{std::numeric_limits<Real>::infinity()};
    while (!stack.empty()) {
      auto const [lower_bound, id] = stack.back();
      stack.pop_back();

      // If the node is farther than the temporary minimum distance, skip it
      if (lower_bound > distance) {
        continue;
      }
      Node const &node{m_tree_structure[id]};

      // Compute the distance between the object and the long boxes on the
      // node
      bool const leaf{node.child_l < 0};
      if (node.box_num > 0 &&
          (leaf || (++c.lng, node.box_long.interior_distance(obj) <= distance))) {
        Integer const id_end{node.box_ptr + node.box_num};
        for (Integer i{node.box_ptr}; i < id_end; ++i) {
          Integer const pos{m_tree_boxes_map[i]};
          ++c.objs;
          Real const tmp_distance{boxes[pos]->interior_distance(obj)};
          if (tmp_distance < distance) {
            candidates.clear();
            candidates.insert(pos);
            distance = tmp_distance;
          } else if (tmp_distance == distance) {
            candidates.insert(pos);
          }
        }
      }

      // Push children on the stack (nearest first)
      if (!leaf) {
        this->push_children_nearest_first(stack, node, distance, c, dist);
      }
    }
    return distance;
  }

  /**
   * Partial result of a tree-tree distance query.
   */
  struct DistanceResult {
    Real distance{std::numeric_limits<Real>::infinity()};
    IndexMap candidates;
  };

  /**
   * Atomically set \c a to \f$ \min(a, v) \f$.
   */
  static void atomic_min(std::atomic<Real> &a, Real const v) {
    Real current{a.load(std::memory_order_relaxed)};
    while (v < current &&
           !a.compare_exchange_weak(current, v, std::memory_order_relaxed)) {
    }
  }

  /**
   * Process one pair of nodes of a tree-tree distance query. The pruning bound
   * is the minimum between the local best distance and the global one shared
   * by all the threads.
   */
  void distance_pair_step(Tree const &tree, Integer const id_s1,
                          Integer const id_s2, DistanceResult &r,
                          std::atomic<Real> &global, QueryCounters &c,
                          IndexList &stack) const {
    BoxUniquePtrList const &boxes_1{*m_boxes_ptr};
    BoxUniquePtrList const &boxes_2{*tree.m_boxes_ptr};
    Integer const id_1{id_s1 < 0 ? negate(id_s1) : id_s1};
    Integer const id_2{id_s2 < 0 ? negate(id_s2) : id_s2};
    Node const &node_1{m_tree_structure[id_1]};
    Node const &node_2{tree.m_tree_structure[id_2]};

    Real bound{std::min(r.distance, global.load(std::memory_order_relaxed))};
    auto update = [&r, &global, &bound](Real const d, Integer const pos_1,
                                        Integer const pos_2) {
      if (d > bound) {
        return;
      }
      if (d < r.distance) {
        r.candidates.clear();
        r.candidates[pos_1].insert(pos_2);
        r.distance = bound = d;
        atomic_min(global, d);
      } else { // d == r.distance
        r.candidates[pos_1].insert(pos_2);
      }
    };

    // If the nodes are farther than the minimum distance, skip them
    if (++c.node; node_1.box.interior_distance(node_2.box) > bound) {
      return;
    }

    // Compute the distance between the long boxes on the nodes
    if (node_1.box_num > 0 && node_2.box_num > 0) {
      if (++c.lng; node_1.box_long.interior_distance(node_2.box_long) <= bound) {
        Integer const id_1_ini{node_1.box_ptr};
        Integer const id_1_end{node_1.box_ptr + node_1.box_num};
        Integer const id_2_ini{node_2.box_ptr};
        Integer const id_2_end{node_2.box_ptr + node_2.box_num};
        if (node_1.box_num < node_2.box_num) {
          for (Integer i{id_1_ini}; i < id_1_end; ++i) {
            Integer const pos_1{m_tree_boxes_map[i]};
            Box const &box_1{*boxes_1[pos_1]};
            if (++c.objs; box_1.interior_distance(node_2.box_long) > bound) {
              continue;
            }
            for (Integer j{id_2_ini}; j < id_2_end; ++j) {
              Integer const pos_2{tree.m_tree_boxes_map[j]};
              ++c.objs;
              update(box_1.interior_distance(*boxes_2[pos_2]), pos_1, pos_2);
            }
          }
        } else {
          for (Integer j{id_2_ini}; j < id_2_end; ++j) {
            Integer const pos_2{tree.m_tree_boxes_map[j]};
            Box const &box_2{*boxes_2[pos_2]};
            if (++c.objs; box_2.interior_distance(node_1.box_long) > bound) {
              continue;
            }
            for (Integer i{id_1_ini}; i < id_1_end; ++i) {
              Integer const pos_1{m_tree_boxes_map[i]};
              ++c.objs;
              update(boxes_1[pos_1]->interior_distance(box_2), pos_1, pos_2);
            }
          }
        }
      }
    }

    // Push the next pairs, sorted so that the nearest pair is visited first
    // (and drop the pairs that are already farther than the bound)
    std::size_t const first{stack.size()};
    push_children_pairs(stack, id_s1, id_s2, node_1, node_2);
    std::size_t const n_new{(stack.size() - first) / 2};
    if (n_new > 0) {
      std::array<std::tuple<Real, Integer, Integer>, 3> pairs;
      std::size_t n_kept{0};
      for (std::size_t k{0}; k < n_new; ++k) {
        Integer const a{stack[first + 2 * k]};
        Integer const b{stack[first + 2 * k + 1]};
        Box const &box_a{a < 0 ? m_tree_structure[negate(a)].box_long
                               : m_tree_structure[a].box};
        Box const &box_b{b < 0 ? tree.m_tree_structure[negate(b)].box_long
                               : tree.m_tree_structure[b].box};
        Real const d{box_a.interior_distance(box_b)};
        if (d <= bound) {
          pairs[n_kept++] = std::make_tuple(d, a, b);
        }
      }
      std::sort(pairs.begin(), pairs.begin() + n_kept,
                [](auto const &x, auto const &y) {
                  return std::get<0>(x) > std::get<0>(y);
                });
      stack.resize(first);
      for (std::size_t k{0}; k < n_kept; ++k) {
        stack.emplace_back(std::get<1>(pairs[k]));
        stack.emplace_back(std::get<2>(pairs[k]));
      }
    }
  }

  /**
   * Compute the minimum distance between the tree and another tree (see
   * \c distance).
   */
  Real distance_tree_impl(Tree const &tree, IndexMap &candidates,
                          Workspace &ws) const {
    candidates.clear();
    if (this->is_empty() || tree.is_empty()) {
      return static_cast<Real>(-1.0);
    }
    std::atomic<Real> global{std::numeric_limits<Real>::infinity()};
    DistanceResult result;
    this->dual_traversal(
        tree, result, ws,
        [this, &tree, &global](Integer const id_s1, Integer const id_s2,
                               DistanceResult &r, Workspace &w,
                               IndexList &stack) {
          this->distance_pair_step(tree, id_s1, id_s2, r, global, w.counters,
                                   stack);
        },
        [](DistanceResult &dst, DistanceResult &src) {
          if (src.distance < dst.distance) {
            dst = std::move(src);
          } else if (src.distance == dst.distance) {
            merge_index_maps(dst.candidates, src.candidates);
          }
        });
    candidates = std::move(result.candidates);
    return result.distance;
  }

  /**
   * Find the first \f$ n \f$ candidates from an object (see \c closest).
   */
  template <typename Object, typename Function>
  Real closest_impl(Object const &obj, Integer const n, IndexSet &candidates,
                    Function &distance_function, Workspace &ws) const {
    candidates.clear();
    if (this->is_empty() || n <= 0) {
      return static_cast<Real>(-1.0);
    }
    BoxUniquePtrList const &boxes{*m_boxes_ptr};
    QueryCounters &c{ws.counters};
    auto dist = [&obj, &distance_function](Box const &b) {
      return distance_function(obj, b);
    };

    // Depth-first traversal visiting the nearest child first
    std::vector<HeapItem> &stack{ws.bstack};
    stack.clear();
    ++c.node;
    stack.emplace_back(dist(m_tree_structure[0].box), 0);

    // Max-heap of the best candidates found so far (distance, index)
    std::vector<HeapItem> &heap{ws.heap};
    heap.clear();
    auto cmp = [](HeapItem const &a, HeapItem const &b) {
      return a.first < b.first;
    };

    // Pruning bound: infinite until n candidates have been found, then the
    // maximum distance among the candidates
    Real bound{std::numeric_limits<Real>::infinity()};
    std::size_t const n_size{static_cast<std::size_t>(n)};

    while (!stack.empty()) {
      auto const [lower_bound, id] = stack.back();
      stack.pop_back();

      // If the node is farther than the bound, skip it
      if (lower_bound > bound) {
        continue;
      }
      Node const &node{m_tree_structure[id]};

      // Compute the distance between the object and the long boxes on the
      // node
      bool const leaf{node.child_l < 0};
      if (node.box_num > 0 &&
          (leaf || (++c.lng, distance_function(obj, node.box_long) <= bound))) {
        Integer const id_end{node.box_ptr + node.box_num};
        for (Integer i{node.box_ptr}; i < id_end; ++i) {
          Integer const pos{m_tree_boxes_map[i]};
          ++c.objs;
          Real const tmp_distance{distance_function(obj, *boxes[pos])};
          if (tmp_distance < bound) {
            if (heap.size() == n_size) {
              std::pop_heap(heap.begin(), heap.end(), cmp);
              heap.pop_back();
            }
            heap.emplace_back(tmp_distance, pos);
            std::push_heap(heap.begin(), heap.end(), cmp);
            if (heap.size() == n_size) {
              bound = heap.front().first;
            }
          }
        }
      }

      // Push children on the stack (nearest first)
      if (!leaf) {
        this->push_children_nearest_first(stack, node, bound, c, dist);
      }
    }

    // Extract indices into candidates
    if (heap.empty()) {
      return static_cast<Real>(-1.0);
    }
    for (HeapItem const &item : heap) {
      candidates.insert(item.second);
    }
    return heap.front().first;
  }

  /**
   * Find the candidates within a given distance (see \c within_distance).
   */
  template <typename Object, typename Function>
  void within_distance_impl(Object const &obj, Real const max_distance,
                            IndexSet &candidates, Function &distance_function,
                            Workspace &ws) const {
    candidates.clear();
    if (this->is_empty()) {
      return;
    }
    BoxUniquePtrList const &boxes{*m_boxes_ptr};
    QueryCounters &c{ws.counters};
    IndexList &stack{ws.stack};
    stack.clear();
    stack.emplace_back(0);
    while (!stack.empty()) {
      Integer const id{stack.back()};
      stack.pop_back();
      Node const &node{m_tree_structure[id]};

      // If the node is farther than the maximum distance, skip it
      if (++c.node; distance_function(obj, node.box) > max_distance) {
        continue;
      }

      // Compute the distance between the object and the long boxes on the
      // node
      bool const leaf{node.child_l < 0};
      if (node.box_num > 0 &&
          (leaf ||
           (++c.lng, distance_function(obj, node.box_long) <= max_distance))) {
        Integer const id_end{node.box_ptr + node.box_num};
        for (Integer i{node.box_ptr}; i < id_end; ++i) {
          Integer const pos{m_tree_boxes_map[i]};
          if (++c.objs; distance_function(obj, *boxes[pos]) <= max_distance) {
            candidates.insert(pos);
          }
        }
      }

      // Push children on the stack
      if (!leaf) {
        stack.emplace_back(node.child_l);
        stack.emplace_back(node.child_r);
      }
    }
  }

}; // class Tree

} // namespace AABBtree

#endif // INCLUDE_AABBTREE_TREE_HXX
