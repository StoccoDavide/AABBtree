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
  template <unsigned N, typename Real=double>
  class NonRecursive {
  public:
    // Basic types for managing the axis-aligned bounding boxes
    using Box = Box<Real, N>; /**< Axis-aligned bounding box in N-dimensional space. */
    using Boxes = vector<unique_ptr<Box>>; /**< Vector of unique pointers to an axis-aligned bounding box. */

    // Basic types for managing the tree
    using Indexes = std::vector<Integer>; /**< Vector of indexes. */

    // Basic types for retrieving indexes
    using Set = set<Integer>; /**< Set of indexes. */
    using Map = map<Integer, Set>; /**< Map of indexes. */

  private:
    // Tree structure
    Boxes m_objs_boxes; /**< Object axis-aligned bounding boxes. */
    Boxes m_tree_boxes; /**< Tree axis-aligned bounding boxes. */

    Integer m_num_objects{0};    /**< Number of objects. */
    Integer m_num_tree_nodes{0}; /**< Number of nodes in the tree. */

    // Tree hierarchy
    Indexes m_parent; /**< Father nodes. */
    Indexes m_children; /**< Child nodes. */
    Indexes m_nodes_id; /**< Nodes ID. */

    // Tree parameters
    Integer m_max_nodes_objects{16}; /**< Maximum number of objects per node. */
    Real m_long_edge_ratio{0.8}; /**< Long edge ratio for bounding boxes. */
    Real m_collision_tolerance{0.1}; /**< Overlap tolerance for bounding boxes. */
    Real m_min_size_tolerance{0.0}; /**< Minimum size tolerance for bounding boxes. */

    // Cache and statistics
    mutable std::vector<Integer> m_stack; /**< Tree stack. */
    mutable Integer m_check_counter{0}; /**< Number of overlap check (for statistic). */

  public:
    /**
    * Class destructor for the \em non-recursive axis-aligned bounding box tree.
    */
    ~NonRecursive() {};

    /**
    * Class constructor for the \em non-recursive axis-aligned bounding box tree.
    */
    NonRecursive() {};

    /**
    * Class constructor for the \em non-recursive axis-aligned bounding box tree.
    * \param[in] tree Axis-aligned bounding box tree to copy.
    */
    NonRecursive(NonRecursive<Real, N> const & tree)
    {
      this->m_max_nodes_objects = tree.m_max_nodes_objects;
      this->m_long_edge_ratio = tree.m_long_edge_ratio;
      this->m_collision_tolerance = tree.m_collision_tolerance;
      this->m_min_size_tolerance = tree.m_min_size_tolerance;
    }

    /**
    * Allocate memory for the \em non-recursive AABB tree given the number of bounding boxes.
    * \param[in] n Number of axis-aligned boxes to allocate.
    */
    void allocate(Integer n)
    {
      this->m_objs_boxes.clear();
      this->m_tree_boxes.clear();
      this->m_parent.clear();
      this->m_children.clear();
      this->m_nodes_id.clear();

      this->m_objs_boxes.reserve(n); // CHECK
      this->m_tree_boxes.reserve(2*n); // CHECK
      this->m_parent.reserve(n); // CHECK
      this->m_children.reserve(n); // CHECK
      this->m_nodes_id.reserve(n); // CHECK
    }

    /**
    * Set the maximum number of objects per node.
    * \param[in] n Maximum number of objects per node.
    */
    void max_nodes_objects(Integer n)
    {
      AABBTREE_ASSERT(n > static_cast<Integer>(0),
        "NonRecursive::max_nodes_objects(...): input must be a positive integer.");
      this->m_max_nodes_objects = n;
    }

    /**
    * Get the maximum number of objects per node.
    * \return The maximum number of objects per node.
    */
    Integer max_nodes_objects() const {return this->m_max_nodes_objects;}

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
    * Add an axis-aligned bounding box to the AABB tree.
    * \param[in] box Pointer to the box to add.
    */
    void add_box( unique_ptr<Box> & box)
    {
        AABBTREE_ASSERT(!box->is_empty(), "NonRecursive::add_box(...): empty box detected.");
        this->m_objs_boxes.push_back(std::move(box));
    }

//    /**
//    * Add an axis-aligned bounding box to the AABB tree.
//    * \param[in] box Box to add.
//    */
//    void add_box(Box const & box)
//    {
//        AABBTREE_ASSERT(!box.is_empty(), "NonRecursive::add_box(...): empty box detected.");
//        this->m_objs_boxes.emplace_back(box);
//    }//

//    /**
//    * Add axis-aligned bounding boxes to the AABB tree.
//    * \param[in] boxes Bounding boxes to add.
//    */
//    void add_boxes(Boxes & boxes)
//    {
//      for (auto & box : boxes) {this->add_box(box);}
//    }//

//    /**
//    * Add axis-aligned bounding boxes to the AABB tree.
//    * \param[in] boxes Bounding boxes to add.
//    */
//    void add_boxes(std::vector<Box> const & boxes)
//    {
//      for (auto const & box : boxes) {this->add_box(box);}
//    }//

//  /**
//  * Build the AABB tree with the already internally stored axis-aligned bounding boxes.
//  */
//  void build()
//  {
//    #define CMD "NonRecursive::build(...): "
//
//    Integer num_objects{this->m_objs_boxes.size()};
//
//    // Build the root box that contains all tree boxes
//    this->m_bbox_tree.emplace_back(*this->m_tree_boxes);
//    for (Integer i{1}; i < num_objects; ++i) {this->m_bbox_tree[0].extend(this->m_tree_boxes[i]);}
//
//    // Clear the stack and reserve memory
//    this->m_stack.clear();
//    this->m_stack.reserve(2*m_num_objects+1);
//    this->m_stack.emplace_back(0);
//    this->m_num_tree_nodes = 1;
//
//    // Divide nodes until all constraints satisfied
//    while (!this->m_stack.empty())
//    {
//
//      // Pop node from stack
//      Integer id_father{this->m_stack.back()};
//      this->m_stack.pop_back();
//
//      // Set no children for the moment
//      this->m_child[id_father] = -1;
//
//      // Get rectangles id in parent
//      Integer num{this->m_num_nodes[id_father]};
//
//      // If there are few boxes stop splitting
//      if (num < this->m_max_num_objects_per_node) {continue;}
//
//      Integer iptr{this->m_ptr_nodes[id_father]};
//      Integer * ptr{this->m_id_nodes + iptr};
//
//      // Split plane on longest axis (euristic)
//      Real const * father_min{this->m_bbox_tree + 2*N*id_father};
//      Real const * father_max{father_min + N};
//
//      Integer idim = 0;
//      Real    mx   = father_max[0] - father_min[0];
//      for (Integer i{1}; i < N; ++i) {
//        Real mx1 = father_max[i] - father_min[i];
//        if (mx < mx1) {mx = mx1; idim = i;}
//      }
//
//      // If too small bbox stop splitting
//      if (mx < this->m_bbox_min_size_tolerance) {continue;}
//
//      Real tol_len{this->m_bbox_long_edge_ratio * mx};
//      Real sp{0.0};
//
//      // Separate short/long and accumulate short baricenter
//      Integer n_long{0};
//      Integer n_short{0};
//      while (n_long + n_short < num)
//      {
//        Integer id = ptr[n_long];
//        Real const * id_min = this->m_bbox_objs + 2*N*id;
//        Real const * id_max = id_min + N;
//        Real id_len = id_max[idim] - id_min[idim];
//        if (id_len > tol_len) {
//          // Found long bbox, increment n_long and update position
//          ++n_long;
//        } else {
//          // Found short bbox, increment n_short and exchange with bottom
//          ++n_short;
//          swap(ptr[n_long], ptr[num-n_short]);
//          sp += id_max[idim] + id_min[idim];
//        }
//      }
//
//      // If split rectangles do not improve search, stop split at this level
//      if (n_short < 2) {continue;}
//
//      // Select the split position: take the mean of the set of (non-"long")
//      // rectangle centers along axis idim
//      sp /= 2*n_short;
//
//      // Partition based on centers
//      Integer n_left{0};
//      Integer n_right{0};
//
//      while (n_long + n_left + n_right < num)
//      {
//        Integer id = ptr[n_long+n_left];
//        Real const * id_min = this->m_bbox_objs + id * this->m_2dim;
//        Real const * id_max = id_min + N;
//        Real id_mid = (id_max[idim] + id_min[idim])/2;
//        if (id_mid < sp) {
//          ++n_left; // In right position do nothing
//        } else {
//          ++n_right;
//          swap(ptr[n_long+n_left], ptr[num-n_right]);
//        }
//      }
//
//      // If cannot improve bbox, stop splitting at this level
//      if (n_left == 0 || n_right == 0) {continue;}
//
//      // Child indexing
//      Integer id_left{this->m_num_tree_nodes + 0};
//      Integer id_right{this->m_num_tree_nodes + 1};
//
//      // Compute bbox of left and right child
//      Box bb_left;
//      Real * bb_left_min = this->m_bbox_tree + id_left * this->m_2dim;
//      Real * bb_left_max = bb_left_min + N;
//      for (Integer i{0}; i < n_left; ++i) {
//        Integer id = ptr[n_long+i];
//        Real const * bb_id_min = this->m_bbox_objs + id * this->m_2dim;
//        Real const * bb_id_max = bb_id_min + N;
//        if (i == 0) {
//          copy_n(bb_id_min, this->m_2dim, bb_left_min);
//          //copy_n(bb_id_max, N, left_max);
//        } else {
//          for (Integer j{0}; j < N; ++j) {
//            if (bb_left_min[j] > bb_id_min[j]) {bb_left_min[j] = bb_id_min[j];}
//            if (bb_left_max[j] < bb_id_max[j]) {bb_left_max[j] = bb_id_max[j];}
//          }
//        }
//      }
//
//      Box bb_right;
//      Real * bb_right_min = this->m_bbox_tree + id_right * this->m_2dim;
//      Real * bb_right_max = bb_right_min + N;
//      for (Integer i{0}; i < n_right; ++i) {
//        Integer id = ptr[n_long+n_left+i];
//        Real const * bb_id_min = this->m_bbox_objs + id * this->m_2dim;
//        Real const * bb_id_max = bb_id_min + N;
//        if (i == 0) {
//          copy_n(bb_id_min, this->m_2dim, bb_right_min);
//        } else {
//          for (Integer j{0}; j < N; ++j) {
//            if (bb_right_min[j] > bb_id_min[j]) {bb_right_min[j] = bb_id_min[j];}
//            if (bb_right_max[j] < bb_id_max[j]) {bb_right_max[j] = bb_id_max[j];}
//          }
//        }
//      }
//
//      // Check again if split improve the tree, otherwise stop the exploration
//      if (n_left < this->m_max_num_objects_per_node || n_right < this->m_max_num_objects_per_node) {
//        // Few nodes, check if improve volume
//        Real vo{bb_left.intersection(bb_right).volume()};
//        Real vol_left{bb_left.volume()};
//        Real vol_right{bb_right.volume()};
//        // If do not improve volume, stop splitting
//        Real otol{std::pow(this->m_bbox_overlap_tolerance, N)};
//        if (vo > (vol_left+vol_right-vo)*otol) {continue;}
//      }
//
//      // Push child nodes onto the stack
//      this->m_father[id_left]  = id_father;
//      this->m_father[id_right] = id_father;
//      this->m_child[id_father] = id_left;
//
//      this->m_num_nodes[id_father] = n_long;
//
//      this->m_ptr_nodes[id_left]  = iptr + n_long;
//      this->m_num_nodes[id_left]  = n_left;
//
//      this->m_ptr_nodes[id_right] = iptr + n_long + n_left;
//      this->m_num_nodes[id_right] = n_right;
//
//      this->m_stack.emplace_back(id_left);
//      this->m_stack.emplace_back(id_right);
//      this->m_num_tree_nodes += 2;
//    }
//
//    #undef CMD
//  }
//
//  /**
//  * Build the AABB tree given the bounding boxes.
//  * \param[in] boxes Bounding boxes to build the tree from.
//  */
//  void build(Boxes const & boxes)
//  {
//    this->allocate(boxes.size());
//    this->add_bboxes(boxes);
//    this->build();
//  }
//
//  /**
//  * Build the AABB tree given the bounding boxes.
//  * \param[in] boxes Bounding boxes to build the tree from.
//  */
//  void build(std::vector<Box> const & boxes)
//  {
//    this->allocate(boxes.size());
//    this->add_bboxes(boxes);
//    this->build();
//  }
//
//  /**
//  * Intersect the AABB tree with a point.
//  * \param[in] point Point to intersect with.
//  * \param[out] candidates Intersection result (bounding box indexes).
//  */
//  void intersect_with_one_point(Point const & point, Set & candidates) const
//  {
//    this->m_num_check = 0;
//
//    // Quick return on empty inputs
//    if (this->m_num_tree_nodes == 0) return;
//
//    // Descend tree from root
//    this->m_stack.clear();
//    this->m_stack.reserve(2*m_num_tree_nodes+1);
//    this->m_stack.emplace_back(0);
//    while (!this->m_stack.empty()) {
//      // Pop node from the stack
//      Integer id_father = this->m_stack.back();
//      this->m_stack.pop_back();
//
//      // Get bbox
//      Real const * bb_father = this->m_bbox_tree + id_father * this->m_2dim;
//
//      ++this->m_num_check;
//      bool overlap = this->m_check_overlap_with_point(point, bb_father, N);
//
//      // If do not overlap skip
//      if (!overlap) {continue;}
//
//      // Get rectangles id in parent
//      this->get_bbox_indexes_of_a_node(id_father, candidates);
//
//      Integer nn = this->m_child[id_father];
//      if (nn > 0) { // root == 0, children > 0
//        // Push children on stack
//        this->m_stack.emplace_back(nn);
//        this->m_stack.emplace_back(nn+1);
//      }
//    }
//  }

  }; // NonRecursive

  template class NonRecursive<float, 1>;
  template class NonRecursive<float, 2>;
  template class NonRecursive<float, 3>;
  template class NonRecursive<float, 4>;
  template class NonRecursive<double, 1>;
  template class NonRecursive<double, 2>;
  template class NonRecursive<double, 3>;
  template class NonRecursive<double, 4>;

} // namespace AABBtree

#endif // INCLUDE_AABBTREE_NONRECURSIVE_HXX
