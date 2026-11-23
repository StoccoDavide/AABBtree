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

#ifndef INCLUDE_AABBTREE_TREE_HH
#define INCLUDE_AABBTREE_TREE_HH

#include "AABBtree/Box.hh"

namespace AABBtree {

namespace detail {

/**
 * Access an object stored either by value or through a pointer (raw or smart).
 */
template <typename T> decltype(auto) as_ref(T const &x) {
  if constexpr (requires { *x; }) {
    return *x;
  } else {
    return x;
  }
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
 * \brief An object that can be intersected with a box, i.e. the box has a
 * method \c intersects(object).
 */
template <typename Object, typename Real, Integer N>
concept IntersectableWithBox =
    requires(Box<Real, N> const &box, Object const &object) {
      { box.intersects(object) } -> std::convertible_to<bool>;
    };

/**
 * \brief An object whose distance from a box can be computed, i.e. the box has
 * a method \c interior_distance(object).
 */
template <typename Object, typename Real, Integer N>
concept DistanceFromBox =
    requires(Box<Real, N> const &box, Object const &object) {
      { box.interior_distance(object) } -> std::convertible_to<Real>;
    };

/**
 * \brief A function computing a distance between an object and a box, called
 * as \c function(object, box).
 */
template <typename Function, typename Object, typename Real, Integer N>
concept BoxDistanceFunction =
    requires(Function &function, Object const &object, Box<Real, N> const &box) {
      { function(object, box) } -> std::convertible_to<Real>;
    };

/**
 * \brief A random-access range of query objects, stored by value or through
 * pointers (raw or smart).
 */
template <typename Container>
concept QueryRange = std::ranges::random_access_range<Container const> &&
                     std::ranges::sized_range<Container const>;

/**
 * \brief Type of the objects stored in a query range (pointers dereferenced).
 */
template <QueryRange Container>
using QueryObject = std::remove_cvref_t<decltype(detail::as_ref(
    std::declval<std::ranges::range_reference_t<Container const>>()))>;

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
 * <b>How the tree is organized</b>
 * - The boxes are owned by the tree (\c m_boxes_ptr) and never moved: queries
 *   report their original indices. \c m_tree_boxes_map is a permutation of
 *   these indices such that every node owns a contiguous range of it, starting
 *   at \c box_ptr and \c box_num long.
 * - A node is split along one axis: the boxes entirely on the left (right) of
 *   the separation line go to the left (right) child, while the boxes crossing
 *   it, called \em long boxes, stay in the node. So every node stores its long
 *   boxes (bounded by \c box_long), and a leaf stores all its boxes as long
 *   boxes (and \c box_long is equal to \c box).
 * - All the traversals use explicit stacks, never recursion.
 *
 * \note Thread safety: once the tree is built, all the \c const methods can be
 * called concurrently from multiple threads. The query statistics returned by
 * \c stats() refer to the last completed query (or batch of queries).
 * \note With OpenMP support (\c _OPENMP defined) \c build, the \c *_batch
 * methods and the tree-tree queries (\c intersect(Tree), \c self_intersect and
 * \c distance(Tree)) run in parallel. The tree built does not depend on the
 * number of threads. The number of threads is set with
 * \c num_threads (1 disables parallelism). Tree-tree queries run serially
 * when called from inside a parallel region or when the trees have fewer than
 * \c AABBTREE_PARALLEL_MIN_NODES nodes in total. Results do not depend on the
 * number of threads.
 * \tparam Real Type of the scalar coefficients
 * \tparam N Dimension of the ambient space.
 */
template <std::floating_point Real, Integer N> class Tree {
public:
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
  using DistanceIndex = std::pair<Real, Integer>; /**< (distance, index). */

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
    Integer node_checks{0};     /**< Checks against node bounding boxes. */
    Integer long_box_checks{0}; /**< Checks against long boxes bounds. */
    Integer object_checks{0};   /**< Checks against single boxes. */
    QueryCounters &operator+=(QueryCounters const &other) {
      node_checks += other.node_checks;
      long_box_checks += other.long_box_checks;
      object_checks += other.object_checks;
      return *this;
    }
  };

  /**
   * Per-thread scratch memory used by the queries.
   */
  struct Workspace {
    IndexList stack; /**< Traversal stack of nodes (or of pairs of nodes). */
    std::vector<DistanceIndex> nearest_stack; /**< Stack of (distance, node). */
    std::vector<DistanceIndex> heap; /**< Best candidates of \c closest. */
    QueryCounters counters;          /**< Query counters. */
    Workspace() {
      stack.reserve(128);
      nearest_stack.reserve(128);
    }
  };

  /**
   * Default distance function: interior distance between an object and a box.
   */
  struct InteriorDistance {
    template <DistanceFromBox<Real, N> Object>
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

  // In the tree-tree traversals a node is referenced either as a whole subtree
  // (by its index, which is >= 0) or as "only its long boxes" (by the negative
  // value -1 - index). The functions below hide this encoding.

  /** Reference to only the long boxes of the node \c id. */
  static constexpr Integer only_long_boxes(Integer const id) { return -1 - id; }

  /** Check if a node reference is to only the long boxes of the node. */
  static constexpr bool is_only_long_boxes(Integer const ref) { return ref < 0; }

  /** Index of the node a reference points to. */
  static constexpr Integer node_index(Integer const ref) {
    return ref < 0 ? -1 - ref : ref;
  }

  /** Check if a node is a leaf. */
  static bool is_leaf(Node const &node) { return node.child_l < 0; }

  /**
   * Push the children of a node on a (lower bound, node) stack so that the
   * nearest child is visited first. Children farther than the bound are not
   * pushed at all.
   */
  template <typename Dist>
  void push_children_nearest_first(std::vector<DistanceIndex> &stack,
                                   Node const &node, Real const bound,
                                   QueryCounters &c, Dist &&dist) const {
    c.node_checks += 2;
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
    m_node_check_counter.store(c.node_checks);
    m_long_check_counter.store(c.long_box_checks);
    m_objs_check_counter.store(c.object_checks);
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

    BoxUniquePtrList const &boxes{*m_boxes_ptr};

    // Clear tree structure
    Integer const num_boxes{static_cast<Integer>(boxes.size())};
    m_tree_structure.clear();
    m_tree_boxes_map.clear();
    m_dump_counter = 0;

    // An empty set of boxes produces an empty tree
    if (num_boxes == 0) {
      return;
    }

    // Setup the root node
    Node root{.box = Box{},
              .box_long = Box{},
              .box_ptr = 0,
              .box_num = num_boxes,
              .box_tot_num = num_boxes,
              .parent = -1,
              .child_l = -1,
              .child_r = -1};

    // Setup the boxes map and compute the root box
    root.box.set_empty();
    m_tree_boxes_map.resize(num_boxes);
    for (Integer i{0}; i < num_boxes; ++i) {
      root.box.extend(*boxes[i]);
      m_tree_boxes_map[i] = i;
    }
    // The root is a leaf until it is split: all its boxes are "long" boxes
    root.box_long = root.box;

    // Build the subtrees. Every node only reorders its own range of the boxes
    // map, so different subtrees can be built concurrently. Large subtrees
    // are split into OpenMP tasks, small ones are built serially.
    IndexList best_map(num_boxes);
    BuildPart top;
#ifdef _OPENMP
    int const num_threads{this->openmp_threads()};
    bool const parallel{num_threads > 1 && !omp_in_parallel() &&
                        num_boxes > AABBTREE_PARALLEL_BUILD_GRAIN};
#pragma omp parallel num_threads(num_threads) if (parallel)
#pragma omp single
#endif
    this->build_part(&top, root, &best_map, 0);

    // Assemble the parts in a fixed order (the result does not depend on the
    // number of threads)
    m_tree_structure.reserve(2 * num_boxes + 1);
    this->assemble_part(top, -1);
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
  template <IntersectableWithBox<Real, N> Object>
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
   * Intersect the tree with another tree.
   * \param[in] tree Tree to intersect with.
   * \param[out] pairs Intersection result: the pairs (i, j) of intersecting
   * boxes, with i a box of this tree and j a box of the other tree, sorted and
   * without duplicates.
   * \return True if the trees intersect, false otherwise.
   * \note This is faster than the \c IndexMap overload for large results.
   */
  bool intersect(Tree const &tree, IndexPairList &pairs) const {
    Workspace ws;
    this->intersect_tree_impl(tree, pairs, ws, false);
    this->store_counters(ws.counters);
    return !pairs.empty();
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
   * Intersect the tree with itself.
   * \param[out] pairs Intersection result: the pairs (i, j) of intersecting
   * boxes, with i < j (each pair is reported once), sorted.
   * \return True if at least two boxes intersect, false otherwise.
   * \note This is faster than the \c IndexMap overload for large results.
   */
  bool self_intersect(IndexPairList &pairs) const {
    Workspace ws;
    this->intersect_tree_impl(*this, pairs, ws, true);
    this->store_counters(ws.counters);
    return !pairs.empty();
  }

  /**
   * Self-intersect the tree (i.e., intersect the tree with itself to find all
   * the intersecting boxes).
   * \param[out] candidates Intersection result: all the boxes that intersect at
   * least another box.
   * \return True if at least two boxes intersect, false otherwise.
   */
  bool self_intersect(IndexSet &candidates) const {
    IndexPairList pairs;
    bool const intersects{this->self_intersect(pairs)};
    IndexList ids;
    ids.reserve(2 * pairs.size());
    for (IndexPair const &p : pairs) {
      ids.emplace_back(p.first);
      ids.emplace_back(p.second);
    }
    std::ranges::sort(ids);
    candidates.clear();
    for (Integer const id : ids) {
      candidates.emplace_hint(candidates.end(), id);
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
  template <DistanceFromBox<Real, N> Object>
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
    requires BoxDistanceFunction<Function, Object, Real, N>
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
    requires BoxDistanceFunction<Function, Object, Real, N>
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
  template <QueryRange Container>
    requires IntersectableWithBox<QueryObject<Container>, Real, N>
  bool intersect_batch(Container const &objs,
                       std::vector<IndexSet> &candidates) const {
    Integer const n{static_cast<Integer>(objs.size())};
    candidates.resize(n);
    this->store_counters(this->parallel_for(n, [&](Integer i, Workspace &ws) {
      this->intersect_impl(detail::as_ref(objs[i]), candidates[i], ws);
    }));
    return std::ranges::any_of(candidates,
                               [](IndexSet const &c) { return !c.empty(); });
  }

  /**
   * Compute the minimum distance between the tree and a batch of objects.
   * \param[in] objs Random-access container of objects.
   * \param[out] distances For each object, the minimum distance.
   * \param[out] candidates For each object, the minimum distance candidates.
   * \tparam Container Type of the container of objects.
   */
  template <QueryRange Container>
    requires DistanceFromBox<QueryObject<Container>, Real, N>
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
  template <QueryRange Container, typename Function = InteriorDistance>
    requires BoxDistanceFunction<Function, QueryObject<Container>, Real, N>
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
  template <QueryRange Container, typename Function = InteriorDistance>
    requires BoxDistanceFunction<Function, QueryObject<Container>, Real, N>
  bool within_distance_batch(Container const &objs, Real const max_distance,
                             std::vector<IndexSet> &candidates,
                             Function distance_function = Function{}) const {
    Integer const n{static_cast<Integer>(objs.size())};
    candidates.resize(n);
    this->store_counters(this->parallel_for(n, [&](Integer i, Workspace &ws) {
      this->within_distance_impl(detail::as_ref(objs[i]), max_distance,
                                 candidates[i], distance_function, ws);
    }));
    return std::ranges::any_of(candidates,
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
  // Tree building
  // ---------------------------------------------------------------------------

  /**
   * A part of the tree built by a single task: either a subtree built
   * serially (\c nodes, with local indices) or a node split into two parts.
   */
  struct BuildPart {
    Node node;                       /**< Split node (if left/right set). */
    std::vector<Node> nodes;         /**< Serially built subtree. */
    std::unique_ptr<BuildPart> left; /**< Left part. */
    std::unique_ptr<BuildPart> right; /**< Right part. */
    Integer dumps{0};                /**< Dumpings performed in this part. */
  };

  /**
   * Split a node into its long boxes and two children.
   * \param[in,out] node Node to split (its long boxes are updated).
   * \param[in,out] best_map Scratch copy of the boxes map (only the node
   * range is used).
   * \param[in,out] dumps Dumpings counter.
   * \param[out] node_l Left child (parent and child indices not set).
   * \param[out] node_r Right child (parent and child indices not set).
   * \return True if the node has been split, false if it is a leaf.
   */
  bool split_node(Node &node, IndexList &best_map, Integer &dumps,
                  Node &node_l, Node &node_r) {
    BoxUniquePtrList const &boxes{*m_boxes_ptr};

    // If the node has less than one object, do not split it
    if (node.box_num <= 1) {
      return false;
    }

    // Compute the separation line and tolerance
    Vector const sizes{node.box.max() - node.box.min()};
    Integer sorting[N];
    std::iota(sorting, sorting + N, 0);
    std::ranges::sort(sorting, std::ranges::greater{},
                      [&sizes](Integer i) { return sizes[i]; });

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
      Integer const n_diff{std::abs(n_left - n_right)};
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

      ++dumps;
    }

    // Terminate splitting if the long boxes cannot be split
    if (n_long_saved == node.box_num) {
      return false;
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

    // Compute the bounding box of the long boxes, and left and right children
    Integer j{node.box_ptr};

    node.box_num = n_long;
    node.box_long.set_empty();
    for (Integer i{0}; i < n_long; ++i) {
      node.box_long.extend(*boxes[m_tree_boxes_map[j++]]);
    }

    // A new leaf owning the next `count` boxes of the map (the parent is set by
    // the caller)
    auto make_leaf = [&](Integer const count) {
      Integer const first{j};
      Box bound; // Empty box
      for (Integer i{0}; i < count; ++i) {
        bound.extend(*boxes[m_tree_boxes_map[j++]]);
      }
      return Node{.box = bound,
                  .box_long = bound, // A leaf: all its boxes are long boxes
                  .box_ptr = first,
                  .box_num = count,
                  .box_tot_num = count,
                  .parent = -1,
                  .child_l = -1,
                  .child_r = -1};
    };
    node_l = make_leaf(n_left);
    node_r = make_leaf(n_right);
    return true;
  }

  /**
   * Build serially (and non-recursively) the subtree rooted at a node.
   * \param[in] root Root of the subtree.
   * \param[out] nodes Nodes of the subtree, with local indices (the root is
   * the first node, its parent is -1).
   * \param[in,out] best_map Scratch copy of the boxes map.
   * \param[in,out] dumps Dumpings counter.
   */
  void build_serial(Node const &root, std::vector<Node> &nodes,
                    IndexList &best_map, Integer &dumps) {
    nodes.clear();
    nodes.reserve(2 * root.box_num + 1);
    nodes.emplace_back(root);
    nodes.back().parent = -1;
    nodes.back().child_l = nodes.back().child_r = -1;
    IndexList stack;
    stack.reserve(128);
    stack.emplace_back(0);
    while (!stack.empty()) {
      Integer const id{stack.back()};
      stack.pop_back();
      Node node_l, node_r;
      if (!this->split_node(nodes[id], best_map, dumps, node_l, node_r)) {
        continue;
      }
      Integer const child_l{static_cast<Integer>(nodes.size())};
      nodes[id].child_l = child_l;
      nodes[id].child_r = child_l + 1;
      node_l.parent = node_r.parent = id;
      nodes.emplace_back(node_l); // May invalidate references to nodes
      nodes.emplace_back(node_r);
      stack.emplace_back(child_l);
      stack.emplace_back(child_l + 1);
    }
  }

  /**
   * Build the part of the tree rooted at a node. Nodes with more than
   * \c AABBTREE_PARALLEL_BUILD_GRAIN boxes are split and their children are
   * built as independent (OpenMP) tasks; smaller subtrees, or subtrees deeper
   * than 64 levels, are built serially.
   */
  void build_part(BuildPart *part, Node node, IndexList *best_map,
                  Integer const level) {
    if (node.box_num <= AABBTREE_PARALLEL_BUILD_GRAIN || level >= 64) {
      this->build_serial(node, part->nodes, *best_map, part->dumps);
      return;
    }
    Node node_l, node_r;
    if (!this->split_node(node, *best_map, part->dumps, node_l, node_r)) {
      part->nodes.emplace_back(node);
      return;
    }
    part->node = node;
    part->left = std::make_unique<BuildPart>();
    part->right = std::make_unique<BuildPart>();
    BuildPart *left{part->left.get()};
    BuildPart *right{part->right.get()};
#ifdef _OPENMP
#pragma omp task firstprivate(left, node_l, best_map, level)
#endif
    this->build_part(left, node_l, best_map, level + 1);
#ifdef _OPENMP
#pragma omp task firstprivate(right, node_r, best_map, level)
#endif
    this->build_part(right, node_r, best_map, level + 1);
#ifdef _OPENMP
#pragma omp taskwait
#endif
  }

  /**
   * Append a built part to the tree structure.
   * \param[in] part Part to append.
   * \param[in] parent Index of the parent node (-1 for the root).
   * \return The index of the root node of the part.
   */
  Integer assemble_part(BuildPart const &part, Integer const parent) {
    m_dump_counter += part.dumps;
    Integer const offset{this->size()};
    if (!part.left) {
      for (Node node : part.nodes) {
        node.parent = node.parent < 0 ? parent : node.parent + offset;
        if (node.child_l >= 0) {
          node.child_l += offset;
          node.child_r += offset;
        }
        m_tree_structure.emplace_back(node);
      }
      return offset;
    }
    m_tree_structure.emplace_back(part.node);
    m_tree_structure.back().parent = parent;
    Integer const child_l{this->assemble_part(*part.left, offset)};
    Integer const child_r{this->assemble_part(*part.right, offset)};
    m_tree_structure[offset].child_l = child_l;
    m_tree_structure[offset].child_r = child_r;
    return offset;
  }

  // ---------------------------------------------------------------------------
  // Parallel driver
  // ---------------------------------------------------------------------------

  /**
   * Run \c body(i, workspace) for \c i in \f$ [0, n) \f$, in parallel when
   * OpenMP is available. Each thread owns its workspace.
   * \return The sum of the query counters of all the threads.
   * \note Exceptions thrown by \c body are rethrown on the calling thread.
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
          continue; // An exception has been thrown: skip the remaining work
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
    dst.merge(src); // Moves the keys that are not in dst yet
    for (auto &[key, values] : src) { // The keys already in dst are left
      dst[key].merge(values);
    }
  }

  // ---------------------------------------------------------------------------
  // Tree-tree traversal
  // ---------------------------------------------------------------------------

  /**
   * Pop a pair of node references from a tree-tree traversal stack.
   */
  static void pop_pair(IndexList &stack, Integer &ref_1, Integer &ref_2) {
    ref_2 = stack.back();
    stack.pop_back();
    ref_1 = stack.back();
    stack.pop_back();
  }

  /**
   * Push a pair of node references on a tree-tree traversal stack.
   */
  static void push_pair(IndexList &stack, Integer const ref_1,
                        Integer const ref_2) {
    stack.emplace_back(ref_1);
    stack.emplace_back(ref_2);
  }

  /**
   * Traverse pairs of nodes of this tree and of another tree, starting from the
   * pair of roots. \c step(ref_1, ref_2, result, workspace, stack) processes a
   * pair of node references and pushes the pairs to visit next on the stack.
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
    push_pair(stack, 0, 0);
#ifdef _OPENMP
    int const num_threads{this->openmp_threads()};
    bool const parallel{num_threads > 1 && !omp_in_parallel() &&
                        this->size() + tree.size() >= AABBTREE_PARALLEL_MIN_NODES};
    if (parallel) {
      // Breadth-first expansion until there are enough pairs to balance the
      // load (a few rounds are enough, the number of pairs grows quickly)
      std::size_t const target_size{static_cast<std::size_t>(2 * 64 * num_threads)};
      IndexList next;
      for (Integer round{0}; round < 64 && !stack.empty() && stack.size() < target_size;
           ++round) {
        next.clear();
        for (std::size_t k{0}; k < stack.size(); k += 2) {
          step(stack[k], stack[k + 1], result, ws, next);
        }
        stack.swap(next);
      }
      if (stack.empty()) {
        return;
      }

      // Depth-first traversal from each pair, in parallel
      IndexList const frontier{stack};
      stack.clear();
      Integer const n_pairs{static_cast<Integer>(frontier.size() / 2)};
      std::vector<Result> partial(n_pairs);
      ws.counters += this->parallel_for(n_pairs, [&](Integer i, Workspace &w) {
        IndexList &local_stack{w.stack};
        local_stack.clear();
        push_pair(local_stack, frontier[2 * i], frontier[2 * i + 1]);
        Integer ref_1, ref_2;
        while (!local_stack.empty()) {
          pop_pair(local_stack, ref_1, ref_2);
          step(ref_1, ref_2, partial[i], w, local_stack);
        }
      });
      for (Result &p : partial) {
        merge(result, p);
      }
      return;
    }
#else
    static_cast<void>(tree);
    static_cast<void>(merge);
#endif
    Integer ref_1, ref_2;
    while (!stack.empty()) {
      pop_pair(stack, ref_1, ref_2);
      step(ref_1, ref_2, result, ws, stack);
    }
  }

  /**
   * Push the pairs of node references to visit after (ref_1, ref_2) in a
   * tree-tree traversal. The node with more boxes is split first; when a node
   * is split, its long boxes are paired with the other node separately.
   */
  static void push_children_pairs(IndexList &stack, Integer const ref_1,
                                  Integer const ref_2, Node const &node_1,
                                  Node const &node_2) {
    Integer const id_1{node_index(ref_1)};
    Integer const id_2{node_index(ref_2)};
    bool const leaf_1{is_leaf(node_1)};
    bool const leaf_2{is_leaf(node_2)};

    if (leaf_1 && leaf_2) {
      return; // Nothing below: this branch is done
    }
    if (leaf_1) {
      push_pair(stack, id_1, node_2.child_l);
      push_pair(stack, id_1, node_2.child_r);
      return;
    }
    if (leaf_2) {
      push_pair(stack, node_1.child_l, id_2);
      push_pair(stack, node_1.child_r, id_2);
      return;
    }

    // Both nodes have children: split the larger one (if it is referenced as
    // a whole subtree), otherwise split the other one
    bool const split_1{node_1.box_tot_num > node_2.box_tot_num
                           ? !is_only_long_boxes(ref_1)
                           : is_only_long_boxes(ref_2)};
    if (split_1 && !is_only_long_boxes(ref_1)) {
      push_pair(stack, node_1.child_l, ref_2);
      push_pair(stack, node_1.child_r, ref_2);
      if (node_1.box_num > 0 && !is_only_long_boxes(ref_2)) {
        push_pair(stack, only_long_boxes(id_1), id_2);
      }
    } else if (!is_only_long_boxes(ref_2)) {
      push_pair(stack, ref_1, node_2.child_l);
      push_pair(stack, ref_1, node_2.child_r);
      if (node_2.box_num > 0 && !is_only_long_boxes(ref_1)) {
        push_pair(stack, id_1, only_long_boxes(id_2));
      }
    }
  }

  /**
   * Push the pairs to visit after (ref_1, ref_2) in a self-intersection. A
   * node paired with itself, (x, x), is split without generating mirrored
   * pairs: (xl, xl), (xr, xr), (xl, xr), (long(x), xl), (long(x), xr). The long
   * boxes of x against themselves are already handled at (x, x). Pairs of
   * distinct nodes follow the general rule.
   */
  static void push_self_children_pairs(IndexList &stack, Integer const ref_1,
                                       Integer const ref_2, Node const &node_1,
                                       Node const &node_2) {
    bool const same_node{ref_1 == ref_2 && !is_only_long_boxes(ref_1)};
    if (!same_node) {
      push_children_pairs(stack, ref_1, ref_2, node_1, node_2);
      return;
    }
    if (is_leaf(node_1)) {
      return; // A leaf with itself: nothing below
    }
    Integer const l{node_1.child_l};
    Integer const r{node_1.child_r};
    push_pair(stack, l, l);
    push_pair(stack, r, r);
    push_pair(stack, l, r);
    if (node_1.box_num > 0) {
      push_pair(stack, only_long_boxes(ref_1), l);
      push_pair(stack, only_long_boxes(ref_1), r);
    }
  }

  /**
   * Visit the pairs made of a long box of \c node_1 (this tree) and a long box
   * of \c node_2 (the other tree).
   *
   * The boxes of the node with fewer long boxes are visited in the outer loop
   * and, as a cheap pre-filter, each of them is first tested against the bound
   * of the long boxes of the other node: if \c near(box, other_box_long) is
   * false, the inner loop is skipped.
   * \c visit(pos_1, box_1, pos_2, box_2) is then called for each pair.
   */
  template <typename Near, typename Visit>
  void for_each_long_box_pair(Tree const &tree, Node const &node_1,
                              Node const &node_2, QueryCounters &c, Near &&near,
                              Visit &&visit) const {
    BoxUniquePtrList const &boxes_1{*m_boxes_ptr};
    BoxUniquePtrList const &boxes_2{*tree.m_boxes_ptr};
    Integer const end_1{node_1.box_ptr + node_1.box_num};
    Integer const end_2{node_2.box_ptr + node_2.box_num};
    if (node_1.box_num < node_2.box_num) {
      for (Integer i{node_1.box_ptr}; i < end_1; ++i) {
        Integer const pos_1{m_tree_boxes_map[i]};
        Box const &box_1{*boxes_1[pos_1]};
        ++c.object_checks;
        if (!near(box_1, node_2.box_long)) {
          continue;
        }
        for (Integer j{node_2.box_ptr}; j < end_2; ++j) {
          Integer const pos_2{tree.m_tree_boxes_map[j]};
          visit(pos_1, box_1, pos_2, *boxes_2[pos_2]);
        }
      }
    } else {
      for (Integer j{node_2.box_ptr}; j < end_2; ++j) {
        Integer const pos_2{tree.m_tree_boxes_map[j]};
        Box const &box_2{*boxes_2[pos_2]};
        ++c.object_checks;
        if (!near(box_2, node_1.box_long)) {
          continue;
        }
        for (Integer i{node_1.box_ptr}; i < end_1; ++i) {
          Integer const pos_1{m_tree_boxes_map[i]};
          visit(pos_1, *boxes_1[pos_1], pos_2, box_2);
        }
      }
    }
  }

  // ---------------------------------------------------------------------------
  // Tree-tree intersection
  // ---------------------------------------------------------------------------

  /**
   * Store an intersecting pair. In a self-intersection each pair is found
   * once: pair lists store it as (min, max), maps in both directions.
   */
  static void store_pair(IndexPairList &pairs, Integer const a, Integer const b,
                         bool const self) {
    if (self && b < a) {
      pairs.emplace_back(b, a);
    } else {
      pairs.emplace_back(a, b);
    }
  }
  static void store_pair(IndexMap &map, Integer const a, Integer const b,
                         bool const self) {
    map[a].insert(b);
    if (self) {
      map[b].insert(a);
    }
  }

  /**
   * Merge two partial results of a tree-tree intersection.
   */
  static void merge_results(IndexPairList &dst, IndexPairList &src) {
    dst.insert(dst.end(), src.begin(), src.end());
    IndexPairList().swap(src); // Release the memory
  }
  static void merge_results(IndexMap &dst, IndexMap &src) {
    merge_index_maps(dst, src);
  }

  /**
   * Finalize the result of a tree-tree intersection: sort the pairs and remove
   * the duplicates (the same pair may be found more than once, e.g. the long
   * boxes of a node against themselves in a self-intersection).
   */
  static void finalize_result(IndexPairList &pairs) {
    std::ranges::sort(pairs);
    auto const duplicates{std::ranges::unique(pairs)};
    pairs.erase(duplicates.begin(), duplicates.end());
  }
  static void finalize_result(IndexMap &) {}

  /**
   * Process one pair of node references of a tree-tree intersection: intersect
   * their long boxes and push the pairs to visit next on the stack. If \c self
   * is true the pairs of a box with itself are discarded.
   */
  template <typename Result>
  void intersect_pair_step(Tree const &tree, Integer const ref_1,
                           Integer const ref_2, Result &result,
                           QueryCounters &c, bool const self,
                           IndexList &stack) const {
    Node const &node_1{m_tree_structure[node_index(ref_1)]};
    Node const &node_2{tree.m_tree_structure[node_index(ref_2)]};

    // If the nodes do not intersect, nothing below them can
    ++c.node_checks;
    if (!node_1.box.intersects(node_2.box)) {
      return;
    }

    // Intersect the long boxes of the two nodes
    if (node_1.box_num > 0 && node_2.box_num > 0) {
      ++c.long_box_checks;
      if (node_1.box_long.intersects(node_2.box_long)) {
        auto near = [](Box const &box, Box const &other_box_long) {
          return box.intersects(other_box_long);
        };
        auto visit = [&](Integer const pos_1, Box const &box_1,
                         Integer const pos_2, Box const &box_2) {
          if (self && pos_1 == pos_2) {
            return; // A box always intersects itself
          }
          ++c.object_checks;
          if (box_1.intersects(box_2)) {
            store_pair(result, pos_1, pos_2, self);
          }
        };
        this->for_each_long_box_pair(tree, node_1, node_2, c, near, visit);
      }
    }

    if (self) {
      push_self_children_pairs(stack, ref_1, ref_2, node_1, node_2);
    } else {
      push_children_pairs(stack, ref_1, ref_2, node_1, node_2);
    }
  }

  /**
   * Intersect the tree with another tree (see \c intersect and
   * \c self_intersect).
   */
  template <typename Result>
  void intersect_tree_impl(Tree const &tree, Result &result, Workspace &ws,
                           bool const self) const {
    result.clear();
    if (this->is_empty() || tree.is_empty()) {
      return;
    }
    auto step = [this, &tree, self](Integer const ref_1, Integer const ref_2,
                                    Result &r, Workspace &w, IndexList &stack) {
      this->intersect_pair_step(tree, ref_1, ref_2, r, w.counters, self, stack);
    };
    auto merge = [](Result &dst, Result &src) { merge_results(dst, src); };
    this->dual_traversal(tree, result, ws, step, merge);
    finalize_result(result);
  }

  // ---------------------------------------------------------------------------
  // Tree-tree distance
  // ---------------------------------------------------------------------------

  /**
   * Partial result of a tree-tree distance query.
   */
  struct DistanceResult {
    Real distance{std::numeric_limits<Real>::infinity()}; /**< Best distance. */
    IndexMap candidates; /**< Pairs of boxes at the best distance. */
  };

  /**
   * Atomically set \c a to \f$ \min(a, v) \f$.
   */
  static void atomic_min(std::atomic<Real> &a, Real const v) {
    Real current{a.load(std::memory_order_relaxed)};
    while (v < current &&
           !a.compare_exchange_weak(current, v, std::memory_order_relaxed)) {
      // On failure current is reloaded: try again while v is still smaller
    }
  }

  /**
   * Bounding box of a node reference (the bound of only its long boxes, or of
   * the whole subtree).
   */
  Box const &referenced_box(Integer const ref) const {
    Node const &node{m_tree_structure[node_index(ref)]};
    return is_only_long_boxes(ref) ? node.box_long : node.box;
  }

  /**
   * Reorder the pairs pushed on the stack from position \c first so that the
   * nearest pair is popped (visited) first, and drop the pairs that are
   * already farther than \c bound.
   */
  void sort_pairs_nearest_first(IndexList &stack, std::size_t const first,
                                Tree const &tree, Real const bound) const {
    // At most 3 pairs are pushed per step: an insertion sort is the simplest
    // (std::sort on such small arrays triggers false -Warray-bounds with GCC)
    std::array<std::tuple<Real, Integer, Integer>, 3> pairs;
    std::size_t n_kept{0};
    for (std::size_t k{first}; k < stack.size(); k += 2) {
      Integer const ref_1{stack[k]};
      Integer const ref_2{stack[k + 1]};
      Real const d{this->referenced_box(ref_1).interior_distance(
          tree.referenced_box(ref_2))};
      if (d <= bound) {
        pairs[n_kept++] = std::make_tuple(d, ref_1, ref_2);
      }
    }
    for (std::size_t k{1}; k < n_kept; ++k) { // Decreasing distance
      for (std::size_t h{k};
           h > 0 && std::get<0>(pairs[h - 1]) < std::get<0>(pairs[h]); --h) {
        std::swap(pairs[h - 1], pairs[h]);
      }
    }
    stack.resize(first);
    for (std::size_t k{0}; k < n_kept; ++k) {
      push_pair(stack, std::get<1>(pairs[k]), std::get<2>(pairs[k]));
    }
  }

  /**
   * Process one pair of node references of a tree-tree distance query. The
   * pruning bound is the minimum between the best distance found by this task
   * and the best distance found by all the threads (\c global).
   */
  void distance_pair_step(Tree const &tree, Integer const ref_1,
                          Integer const ref_2, DistanceResult &result,
                          std::atomic<Real> &global, QueryCounters &c,
                          IndexList &stack) const {
    Node const &node_1{m_tree_structure[node_index(ref_1)]};
    Node const &node_2{tree.m_tree_structure[node_index(ref_2)]};
    Real bound{std::min(result.distance, global.load(std::memory_order_relaxed))};

    // If the nodes are farther than the bound, nothing below them is closer
    ++c.node_checks;
    if (node_1.box.interior_distance(node_2.box) > bound) {
      return;
    }

    // Compare the long boxes of the two nodes
    if (node_1.box_num > 0 && node_2.box_num > 0) {
      ++c.long_box_checks;
      if (node_1.box_long.interior_distance(node_2.box_long) <= bound) {
        auto near = [&bound](Box const &box, Box const &other_box_long) {
          return box.interior_distance(other_box_long) <= bound;
        };
        auto visit = [&](Integer const pos_1, Box const &box_1,
                         Integer const pos_2, Box const &box_2) {
          ++c.object_checks;
          Real const d{box_1.interior_distance(box_2)};
          if (d < result.distance) { // New best pair
            result.candidates.clear();
            result.candidates[pos_1].insert(pos_2);
            result.distance = d;
            bound = std::min(bound, d);
            atomic_min(global, d);
          } else if (d == result.distance) { // Tie with the best pair
            result.candidates[pos_1].insert(pos_2);
          }
        };
        this->for_each_long_box_pair(tree, node_1, node_2, c, near, visit);
      }
    }

    // Push the next pairs, the nearest one on top
    std::size_t const first{stack.size()};
    push_children_pairs(stack, ref_1, ref_2, node_1, node_2);
    this->sort_pairs_nearest_first(stack, first, tree, bound);
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
    auto step = [this, &tree, &global](Integer const ref_1, Integer const ref_2,
                                       DistanceResult &r, Workspace &w,
                                       IndexList &stack) {
      this->distance_pair_step(tree, ref_1, ref_2, r, global, w.counters, stack);
    };
    auto merge = [](DistanceResult &dst, DistanceResult &src) {
      if (src.distance < dst.distance) {
        dst = std::move(src);
      } else if (src.distance == dst.distance) {
        merge_index_maps(dst.candidates, src.candidates);
      }
    };
    this->dual_traversal(tree, result, ws, step, merge);
    candidates = std::move(result.candidates);
    return result.distance;
  }

  // ---------------------------------------------------------------------------
  // Object queries
  // ---------------------------------------------------------------------------

  /**
   * Check whether the long boxes of a node can satisfy a query. For a leaf the
   * long boxes are all the boxes of the node, and \c box_long is equal to
   * \c box, already checked by the caller: the check is skipped.
   * \param[in] node Node to check.
   * \param[in,out] c Query counters.
   * \param[in] check Check of the bound of the long boxes, \c check(box_long).
   */
  template <typename Check>
  static bool long_boxes_may_match(Node const &node, QueryCounters &c,
                                   Check &&check) {
    if (node.box_num == 0) {
      return false;
    }
    if (is_leaf(node)) {
      return true;
    }
    ++c.long_box_checks;
    return check(node.box_long);
  }

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
    auto intersects = [&obj](Box const &box) { return box.intersects(obj); };

    IndexList &stack{ws.stack};
    stack.clear();
    stack.emplace_back(0);
    while (!stack.empty()) {
      Node const &node{m_tree_structure[stack.back()]};
      stack.pop_back();

      // If the object does not intersect the node, it cannot intersect
      // anything below it
      ++c.node_checks;
      if (!intersects(node.box)) {
        continue;
      }

      // Intersect the object with the long boxes of the node
      if (long_boxes_may_match(node, c, intersects)) {
        Integer const end{node.box_ptr + node.box_num};
        for (Integer i{node.box_ptr}; i < end; ++i) {
          Integer const pos{m_tree_boxes_map[i]};
          ++c.object_checks;
          if (intersects(*boxes[pos])) {
            candidates.insert(pos);
          }
        }
      }

      if (!is_leaf(node)) {
        stack.emplace_back(node.child_l);
        stack.emplace_back(node.child_r);
      }
    }
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
    auto distance_to = [&obj](Box const &box) {
      return box.interior_distance(obj);
    };
    Real distance{std::numeric_limits<Real>::infinity()}; // Best so far
    auto may_be_closer = [&](Box const &box) {
      return distance_to(box) <= distance;
    };

    // Depth-first traversal visiting the nearest child first, so that the
    // best distance (used to prune) decreases quickly. Each node is stored
    // with its distance, a lower bound of the distance of its boxes.
    std::vector<DistanceIndex> &stack{ws.nearest_stack};
    stack.clear();
    ++c.node_checks;
    stack.emplace_back(distance_to(m_tree_structure[0].box), 0);
    while (!stack.empty()) {
      auto const [node_distance, id] = stack.back();
      stack.pop_back();
      if (node_distance > distance) {
        continue; // Nothing below this node can be closer
      }
      Node const &node{m_tree_structure[id]};

      if (long_boxes_may_match(node, c, may_be_closer)) {
        Integer const end{node.box_ptr + node.box_num};
        for (Integer i{node.box_ptr}; i < end; ++i) {
          Integer const pos{m_tree_boxes_map[i]};
          ++c.object_checks;
          Real const d{distance_to(*boxes[pos])};
          if (d < distance) { // New closest box
            candidates.clear();
            candidates.insert(pos);
            distance = d;
          } else if (d == distance) { // Tie with the closest box
            candidates.insert(pos);
          }
        }
      }

      if (!is_leaf(node)) {
        this->push_children_nearest_first(stack, node, distance, c,
                                          distance_to);
      }
    }
    return distance;
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
    auto distance_to = [&obj, &distance_function](Box const &box) {
      return distance_function(obj, box);
    };

    // The best candidates found so far, in a max-heap: the farthest on top
    std::vector<DistanceIndex> &best{ws.heap};
    best.clear();
    auto farther = [](DistanceIndex const &a, DistanceIndex const &b) {
      return a.first < b.first;
    };
    std::size_t const n_wanted{static_cast<std::size_t>(n)};

    // Pruning bound: infinite until n candidates have been found, then the
    // distance of the farthest candidate
    Real bound{std::numeric_limits<Real>::infinity()};
    auto may_be_closer = [&](Box const &box) { return distance_to(box) <= bound; };

    // Depth-first traversal visiting the nearest child first
    std::vector<DistanceIndex> &stack{ws.nearest_stack};
    stack.clear();
    ++c.node_checks;
    stack.emplace_back(distance_to(m_tree_structure[0].box), 0);
    while (!stack.empty()) {
      auto const [node_distance, id] = stack.back();
      stack.pop_back();
      if (node_distance > bound) {
        continue; // Nothing below this node can be closer
      }
      Node const &node{m_tree_structure[id]};

      if (long_boxes_may_match(node, c, may_be_closer)) {
        Integer const end{node.box_ptr + node.box_num};
        for (Integer i{node.box_ptr}; i < end; ++i) {
          Integer const pos{m_tree_boxes_map[i]};
          ++c.object_checks;
          Real const d{distance_to(*boxes[pos])};
          if (d >= bound) {
            continue;
          }
          if (best.size() == n_wanted) { // Replace the farthest candidate
            std::pop_heap(best.begin(), best.end(), farther);
            best.pop_back();
          }
          best.emplace_back(d, pos);
          std::push_heap(best.begin(), best.end(), farther);
          if (best.size() == n_wanted) {
            bound = best.front().first;
          }
        }
      }

      if (!is_leaf(node)) {
        this->push_children_nearest_first(stack, node, bound, c, distance_to);
      }
    }

    if (best.empty()) {
      return static_cast<Real>(-1.0);
    }
    for (DistanceIndex const &candidate : best) {
      candidates.insert(candidate.second);
    }
    return best.front().first; // Distance of the farthest candidate
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
    auto within = [&obj, &distance_function, max_distance](Box const &box) {
      return distance_function(obj, box) <= max_distance;
    };

    IndexList &stack{ws.stack};
    stack.clear();
    stack.emplace_back(0);
    while (!stack.empty()) {
      Node const &node{m_tree_structure[stack.back()]};
      stack.pop_back();

      // If the node is too far, everything below it is too far
      ++c.node_checks;
      if (!within(node.box)) {
        continue;
      }

      if (long_boxes_may_match(node, c, within)) {
        Integer const end{node.box_ptr + node.box_num};
        for (Integer i{node.box_ptr}; i < end; ++i) {
          Integer const pos{m_tree_boxes_map[i]};
          ++c.object_checks;
          if (within(*boxes[pos])) {
            candidates.insert(pos);
          }
        }
      }

      if (!is_leaf(node)) {
        stack.emplace_back(node.child_l);
        stack.emplace_back(node.child_r);
      }
    }
  }

}; // class Tree

} // namespace AABBtree

#endif // INCLUDE_AABBTREE_TREE_HH
