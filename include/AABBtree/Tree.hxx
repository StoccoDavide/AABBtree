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

#include <unordered_set>

#include "AABBtree/Box.hxx"
#include "AABBtree/Tree.hxx"

namespace AABBtree {

  /**
   * \brief A class representing a \em non-recursive axis-aligned bounding box tree (AABB tree).
   *
   * The Tree class provides an efficient way to store and query a set of axis-aligned bounding boxes
   * \em non-recursively. It supports various operations such as building the tree, adding axis-aligned
   * boxes, and performing intersection queries. The Tree class is templated on the type of the
   * real numbers and the dimensions of the bounding boxes and is based, as the name suggests, on a
   * novel non-recursive algorithm.
   * \tparam Real Type of the scalar coefficients
   * \tparam N Dimension of the ambient space.
   */
  template <typename Real, Integer N>
  class Tree {
  public:
    static_assert(std::is_floating_point<Real>::value, "Tree Real type must be a floating-point type.");
    static_assert(std::is_integral<Integer>::value, "Tree dimension type must be an integer type." );
    static_assert(N > 0, "Tree dimension must be positive.");

    /**
     * \brief Structure representing the statistics of the AABB tree.
     */
    struct Statistics {
      // Build statistics
      Integer objects{0};          //< Number of objects in the tree.
      Integer nodes{0};            //< Total number of nodes in the tree.
      Integer leafs{0};            //< Number of leaf nodes (nodes with objects).
      Integer long_boxes{0};       //< Number of boxes with aspect ratio > threshold.
      Integer depth{0};            //< Maximum depth of the tree.
      Integer left_nodes{0};       //< Number of nodes in left subtree.
      Integer left_leafs{0};       //< Number of leafs in left subtree.
      Integer left_long_boxes{0};  //< Number of long boxes in left subtree.
      Integer left_depth{0};       //< Depth of left subtree.
      Integer right_nodes{0};      //< Number of nodes in right subtree.
      Integer right_leafs{0};      //< Number of leafs in right subtree.
      Integer right_long_boxes{0}; //< Number of long boxes in right subtree.
      Integer right_depth{0};      //< Depth of right subtree.
      Integer dump_counter{0};     //< Number of tree rebuilds (dumps).
      Real    balance_ratio{0.0};  //< Ratio of left/right subtree sizes (0=perfect balance).
      Real    depth_ratio{0.0};    //< Ratio of left/right subtree depths (0=perfect balance).
  
      // Query statistics
      Integer check_counter{0}; //< Total number of collision checks performed.
  
      /**
       * \brief Reset all statistics to their initial zero state
       * 
       * This method efficiently resets all statistics counters to zero and ratios to 0,
       * returning the structure to its default-constructed state.
       * 
       * \note This operation is noexcept and has constant time complexity O(1)
       * \post All counters are zero and ratios are 0
       */
      void reset() noexcept { *this = Statistics{}; }

      /**
       * \brief Print the tree info to an output stream.
       * \param[in] os Output stream to print the tree info to.
       */
      void print( OutStream & os ) const {
        // Print the tree info
        os << "────────────────────────────────────────────────────────────────────────────────\n"
           << "AABB TREE INFO"
           << "\n\tAmbient dimension : " << N
           << "\n\tBasic type        : " << typeid(Real).name()
           << "\n\tObjects           : " << objects
           << "\n\tNodes             : " << nodes
           << "\n\tLeafs             : " << leafs
           << "\n\tLong boxes        : " << long_boxes
           << "\n\tDepth             : " << depth
           << "\n\tLeft nodes        : " << left_nodes
           << "\n\tLeft leafs        : " << left_leafs
           << "\n\tLeft long boxes   : " << left_long_boxes
           << "\n\tLeft depth        : " << left_depth
           << "\n\tRight nodes       : " << right_nodes
           << "\n\tRight leafs       : " << right_leafs
           << "\n\tRight long boxes  : " << right_long_boxes
           << "\n\tRight depth       : " << right_depth
           << "\n\tDump counter      : " << dump_counter
           << "\n\tBalance ratio     : " << balance_ratio
           << "\n\tDepth ratio       : " << depth_ratio
           << "\n\tCheck counter     : " << check_counter
           << "\n────────────────────────────────────────────────────────────────────────────────\n";
      }

    };

  private:
    // Basic types definitions
    using Box              = Box<Real, N>;
    using BoxUniquePtr     = BoxUniquePtr<Real, N>;
    using BoxUniquePtrList = BoxUniquePtrList<Real, N>;
    using Vector           = Vector<Real, N>;
    using Point            = Point<Real, N>;

    /**
     * \brief Structure representing a node of the AABB tree.
     */
    struct Node {
      Box     box;         //< Bounding box of the subtree.
      Box     box_long;    //< Bounding box of long boxes.
      Integer box_ptr;     //< Pointer to the first box in the reordering map of boxes.
      Integer box_num;     //< Number of boxes in the subtree.
      Integer box_tot_num; //< Total number of boxes in the subtree.
      Integer parent;      //< Root node of the subtree.
      Integer child_l;     //< Left child of the subtree.
      Integer child_r;     //< Right child of the subtree.
    };

    // Tree hierarchy
    std::unique_ptr<BoxUniquePtrList> m_boxes_ptr{nullptr};
    std::vector<Node>                 m_tree_structure;        //< Tree structure.
    IndexList                         m_tree_boxes_map;        //< Reordering between the vector of boxes and the tree internal structure.
    bool                              m_dumping_mode{true};    //< Enable dumping while building the tree.

    // Tree parameters
    Integer m_max_nodal_objects{10};            //< Maximum number of objects per node.
    Real    m_separation_ratio_tolerance{0.25}; //< Tolerance for bounding boxes separation.
    Real    m_balance_ratio_tolerance{0.1};     //< Tolerance for bounding boxes balance.
    Real    m_min_box_size{0.0};                //< Minimum size tolerance for bounding boxes.

    // Statistics
    mutable Integer m_check_counter{0}; //< Number of collision check.
    mutable Integer m_dump_counter{0};  //< Number of dumpings.

    // Stack for non-recursive tree building and query
    mutable IndexList m_stack; //< Tree stack.

  public:
    /**
     * \brief Class destructor for the \em non-recursive axis-aligned bounding box tree.
     */
    ~Tree() = default;

    /**
     * \brief Class constructor for the \em non-recursive axis-aligned bounding box tree.
     */
    Tree() = default;

    /**
     * \brief Get a look at the vector of unique pointers to the bounding boxes.
     * \return A const reference to the vector of unique pointers to the bounding boxes.
     */
    BoxUniquePtrList const & boxes() const { return *m_boxes_ptr; }

    /**
     * \brief Get a look at the i-th unique pointer to the bounding box.
     * \param[in] i Index of the unique pointer to the bounding box.
     * \return A const reference to the i-th unique pointer to the bounding box.
     */
    BoxUniquePtr const & box(Integer const i) const { return (*m_boxes_ptr)[i]; }

    /**
     * \brief Set the maximum number of objects per node.
     * \param[in] n Maximum number of objects per node.
     */
    void max_nodal_objects( Integer const n ) {
      constexpr char CMD[]{"AABBtree::Tree::max_nodal_objects(...): "};
      AABBTREE_ASSERT(n > 0, CMD << "input must be a positive integer.");
      m_max_nodal_objects = n;
    }

    /**
     * \brief Get the maximum number of objects per node.
     * \return The maximum number of objects per node.
     */
    Integer max_nodal_objects() const { return m_max_nodal_objects; }

    /**
     * \brief Set the balance ratio tolerance for bounding boxes.
     * \param[in] ratio Balance ratio tolerance for bounding boxes.
     */
    void separation_ratio_tolerance( Real const ratio ) {
      constexpr char CMD[]{ "AABBtree::Tree::separation_ratio_tolerance(...): " };
      AABBTREE_ASSERT(ratio > 0.0 && ratio < 1.0, CMD << "input must be in the range [0, 1].");
      m_separation_ratio_tolerance = ratio;
    }

    /**
     * \brief Get the balance ratio tolerance for bounding boxes.
     * \return The balance ratio tolerance for bounding boxes.
     */
    Real separation_ratio_tolerance() const { return m_separation_ratio_tolerance; }

    /**
     * \brief Set the minimum size for bounding boxes.
     * \param[in] size Minimum size for bounding boxes.
     */
    void min_box_size( Real const size ) {
      constexpr char CMD[]{ "AABBtree::Tree::min_box_size(...): " };
      AABBTREE_ASSERT( size >= 0.0, CMD << "input must be a non-negative real number." );
      m_min_box_size = size;
    }

    /**
     * \brief Get the minimum size for bounding boxes.
     * \return The minimum size for bounding boxes.
     */
    Real min_box_size() const { return m_min_box_size; }

    /**
     * \brief Get the tree structure.
     * \return The tree structure.
     */
    std::vector<Node> const & structure() const { return m_tree_structure; }

    /**
     * \brief Get the i-th tree node.
     * \param[in] i Index of the node.
     * \return The i-th tree node.
     */
    Node const & node( Integer const i ) const { return m_tree_structure[i]; }

    /**
     * \brief Get the number of nodes in the tree.
     * \return The number of nodes in the tree.
     */
    Integer size() const { return static_cast<Integer>( m_tree_structure.size() ); }

    /**
     * \brief Enable dumping mode while building the tree.
     */
    void enable_dumping_mode() { m_dumping_mode = true; }

    /**
     * \brief Disable dumping mode while building the tree.
     */
    void disable_dumping_mode() { m_dumping_mode = false; }

    /**
     * \brief Set dumping mode while building the tree.
     * \param[in] mode Dumping mode.
     */
    void dumping_mode( bool const mode ) { m_dumping_mode = mode; }

    /**
     * \brief Check if tree is empty.
     * \return True if the tree is empty, false otherwise.
     */
    bool is_empty() const { return m_tree_structure.empty(); }

    /**
     * \brief Clear the tree.
     */
    void clear() {
      m_tree_structure.clear();
      m_tree_boxes_map.clear();
      m_stack.clear();
    }

    /**
     * \brief Build the tree given the bounding boxes.
     * \param[in] boxes Bounding boxes to build the tree from.
     */
    void
    build( std::unique_ptr<BoxUniquePtrList> boxes_ptr ) {

      // Collect the original object boxes
      m_boxes_ptr = std::move(boxes_ptr);

      BoxUniquePtrList & boxes{ *m_boxes_ptr };

      // Clear tree structure
      Integer num_boxes{ static_cast<Integer>( boxes.size() ) };
      m_tree_structure.clear();
      m_tree_structure.reserve(20*num_boxes + 1);

      // Setup the root node
      Node root;
      root.parent      = -1;
      root.child_l     = -1;
      root.child_r     = -1;
      root.box_ptr     = 0;
      root.box_num     = 0;
      root.box_tot_num = num_boxes;

      // Setup the boxes map and compute the root box
      Integer depth{ std::ceil( std::log2(num_boxes) ) };
      root.box.set_empty();
      m_tree_boxes_map.reserve(2*depth);
      for ( BoxUniquePtr const & box : boxes ) {
        root.box.extend(*box);
        m_tree_boxes_map.emplace_back(root.box_num++);
      }
      m_tree_structure.emplace_back(root);

      IndexList m_map(m_tree_boxes_map.size());

      // Setup the stack
      m_stack.clear();
      m_stack.reserve(2*num_boxes + 1);
      m_stack.emplace_back(0);

      // Main loop that divide the nodes iteratively until all constraints satisfied
      while ( !m_stack.empty() ) {
        // Pop the node from stack
        Integer const id{ m_stack.back() }; m_stack.pop_back();

        // Get the node
        Node & node{ m_tree_structure[id] };

        // If the node has less than the maximum number of objects, skip it
        if ( node.box_num < m_max_nodal_objects ) continue;

        // Compute the separation line and tolerance
        Vector  sizes{ node.box.max() - node.box.min() };
        Integer sorting[N]; std::iota(sorting, sorting+N, 0);
        auto compare = [&sizes](Integer i, Integer j) -> bool { return sizes[i] > sizes[j]; };
        std::sort( sorting, sorting+N, compare );
        //std::make_heap( sorting, sorting+N, compare );

        Integer axis, n_long, n_left, n_right, id_ini, id_end;
        Real    separation_line;

        Integer n_long_saved{node.box_num+1};
        Integer n_diff_saved{node.box_num+1};
        Integer axis_saved{sorting[0]};

        // salva la porzione di indici nella zona di lavoro
        //std::copy_n( m_tree_boxes_map.data()+node.box_ptr, node.box_num+1, m_map.data()+node.box_ptr );

        for ( Integer dump{0}; dump < N; ++dump ) {
          axis            = sorting[dump];
          separation_line = node.box.baricenter(axis);

          // Separate short and long boxes and compute short boxes baricenter
          n_long = n_left = n_right = 0;
          id_ini = node.box_ptr;
          id_end = node.box_ptr + node.box_num;

          Real baricenter           { 0 };
          Real separation_tolerance { sizes[axis] * m_separation_ratio_tolerance };
          while ( id_ini < id_end ) {
            Box const & box_id{ *boxes[m_tree_boxes_map[id_ini]] };
            typename Box::Side const side{ box_id.which_side(separation_line, separation_tolerance, axis) };
            switch (side) {
              case Box::Side::LEFT: // Left boxes are moved to the end
                ++n_left;  --id_end; std::swap( m_tree_boxes_map[id_ini], m_tree_boxes_map[id_end] );
                break;
              case Box::Side::RIGHT: // Right boxes are moved to the end
                ++n_right; --id_end; std::swap( m_tree_boxes_map[id_ini], m_tree_boxes_map[id_end] );
                break;
              default:
                ++n_long; ++id_ini;
            }
            baricenter += box_id.max()[axis] + box_id.min()[axis];
          }
          baricenter /= 2.0*node.box_num;

          // salva soluzione migliore se migliora n_long o n_diff
          Integer n_diff{ std::abs( n_left - n_right ) };
          if ( n_long_saved > n_long || n_diff_saved > n_diff ) {
            n_long_saved = n_long;
            n_diff_saved = n_diff;
            axis_saved   = axis;
            std::copy_n( m_tree_boxes_map.data()+node.box_ptr, node.box_num, m_map.data()+node.box_ptr );
          }

          //if ( m_dumping_mode ) break;
          if ( n_long > m_max_nodal_objects ) continue;

          // If the left and right children are yet not balanced, dump the splitting axis
          if ( std::min(n_left, n_right) >= m_balance_ratio_tolerance * std::max(n_left, n_right) ) break;
        }

        // usa migliore salvato
        std::copy_n( m_map.data()+node.box_ptr, node.box_num, m_tree_boxes_map.data()+node.box_ptr );
        n_long = n_long_saved;
        axis   = axis_saved;
        separation_line = node.box.baricenter(axis);

        // Separate the left and right boxes
        n_left = n_right = 0;
        id_ini = node.box_ptr + n_long;
        id_end = node.box_ptr + node.box_num;
        while ( id_ini < id_end ) {
          Integer ipos{ m_tree_boxes_map[id_ini] };
          if ( boxes[ipos]->baricenter(axis) < separation_line ) {
            ++id_ini; ++n_left; // In right position do nothing
          } else {
            --id_end; ++n_right; // In right position swap the current box with the last one
            std::swap( m_tree_boxes_map[id_ini], m_tree_boxes_map[id_end] );
          }
        }

        // If the left and right children have few boxes, they are leafs
        if ( n_left <= m_max_nodal_objects && n_right <= m_max_nodal_objects ) continue;

        // Set the left and right children indexes
        node.child_l = static_cast<Integer>(this->size() + 0);
        node.child_r = static_cast<Integer>(this->size() + 1);

        // Finalize the root node setup (left and right children)
        Node node_l;
        node_l.parent      = id; // Current node
        node_l.child_l     = -1; // Default (leaf)
        node_l.child_r     = -1; // Default (leaf)
        node_l.box_num     = n_left;
        node_l.box_tot_num = n_left;

        Node node_r;
        node_r.parent      = id; // Current node
        node_r.child_l     = -1; // Default (leaf)
        node_r.child_r     = -1; // Default (leaf)
        node_r.box_num     = n_right;
        node_r.box_tot_num = n_right;

        // Compute the bounding box of the long boxes, and left and right children
        Integer j{node.box_ptr};

        node.box_num = n_long;
        node.box_long.set_empty();
        for ( Integer i{0}; i < n_long; ++i ) 
          node.box_long.extend( *boxes[ m_tree_boxes_map[j++] ] );
        
        node_l.box.set_empty();
        node_l.box_ptr = j;
        for (Integer i{0}; i < n_left; ++i)
          node_l.box.extend( *boxes[m_tree_boxes_map[j++] ] );
        node_l.box_long = node_l.box;

        node_r.box.set_empty();
        node_r.box_ptr = j;
        for (Integer i{0}; i < n_right; ++i)
          node_r.box.extend( *boxes[m_tree_boxes_map[j++] ] );
        node_r.box_long = node_r.box;

        // Push nodes on tree structure
        m_tree_structure.emplace_back(node_l);
        m_tree_structure.emplace_back(node_r);

        // Push children on stack
        m_stack.emplace_back(node.child_l);
        m_stack.emplace_back(node.child_r);
      }
    }

    /**
     * \brief Intersect the tree with an object.
     * \param[in] obj Object to intersect with.
     * \param[out] candidates Intersection result (boxes indexes).
     * \return True if the object intersects the tree, false otherwise.
     * \tparam Object Type of the object to intersect with.
     * \note Object must have a method \c intersects that computes the intersection with a box.
     */
    template <typename Object>
    bool
    intersect( Object const & obj, IndexSet & candidates ) const {
      // Reset statistics
      m_check_counter = 0;

      // Return if the tree is empty
      if ( this->is_empty() ) return false;

      // Collect the original object boxes
      BoxUniquePtrList const & boxes{ *m_boxes_ptr };

      // Setup the stack
      m_stack.clear();
      m_stack.reserve(2*this->size() + 1);
      m_stack.emplace_back(0);

      // Main loop that checks the intersection iteratively
      candidates.clear();
      while (!m_stack.empty()) {

        // Pop the node from stack
        Integer const id{ m_stack.back() }; m_stack.pop_back();

        // Get the node
        Node const & node{ m_tree_structure[id] };

        // If the object do not intersects the box, skip the node
        if ( ++m_check_counter; !node.box.intersects(obj) ) continue;

        // Intersect the object with the long boxes on the node.
        // If it is a leaf long boxes are also the nodes of the leaf.
        if ( node.box_num > 0 ) {
          if ( ++m_check_counter; node.box_long.intersects(obj) ) {
            Integer const id_ini{node.box_ptr};
            Integer const id_end{node.box_ptr + node.box_num};
            for ( Integer i{id_ini}; i < id_end; ++i ) {
              Integer const pos{ m_tree_boxes_map[i] };
              if ( ++m_check_counter; boxes[pos]->intersects(obj) ) candidates.insert(pos);
            }
          }
        }

        // Push children on the stack if they are not leafs
        if ( node.child_l > 0 ) m_stack.emplace_back( node.child_l );
        if ( node.child_r > 0 ) m_stack.emplace_back( node.child_r );
      }

      // Return true if the object intersects the tree
      return !candidates.empty();
    }

    /**
     * \brief Intersect the tree with another tree.
     * \param[in] tree Tree to intersect with.
     * \param[out] candidates Intersection result (boxes indexes).
     * \return True if the point intersects the tree, false otherwise.
     */
    bool
    intersect( Tree const & tree, IndexMap & candidates ) const {
      // Reset statistics
      m_check_counter = 0;

      // Return if the tree is empty
      if ( this->is_empty() || tree.is_empty() ) return false;

      // Collect the original object boxes
      BoxUniquePtrList const & boxes_1{ *m_boxes_ptr      };
      BoxUniquePtrList const & boxes_2{ *tree.m_boxes_ptr };

      // Setup the stack
      m_stack.clear();
      m_stack.reserve( this->size() + tree.size() + 2 );
      m_stack.emplace_back(0);
      m_stack.emplace_back(0);

      // Negate the id of the node to distinguish the tree
      // if 0 --> -1, 1 --> -2, 2 --> -3
      auto negate = [] ( Integer const id ) { return -1-id; };

      // Main loop that checks the intersection iteratively
      candidates.clear();
      while ( !m_stack.empty() ) {

        // Pop the node from stack (reversed order)
        Integer const id_s2 { m_stack.back() }; m_stack.pop_back();
        Integer const id_2  { id_s2 < 0 ? negate(id_s2) : id_s2 };
        Integer const id_s1 { m_stack.back() }; m_stack.pop_back();
        Integer const id_1  { id_s1 < 0 ? negate(id_s1) : id_s1 };

        // Get the node
        Node const & node_1{ m_tree_structure[id_1] };
        Node const & node_2{ tree.m_tree_structure[id_2] };

        // If the boxes are not intersecting, skip the nodes
        if ( ++m_check_counter; !node_1.box.intersects(node_2.box) ) continue;

        // Intersect the long boxes on the nodes
        // If both are leaf then intersect the corresponding boxes
        if ( node_1.box_num > 0 && node_2.box_num > 0 ) {
          if ( ++m_check_counter; node_1.box_long.intersects(node_2.box_long) ) {
            Integer const id_1_ini{ node_1.box_ptr                  };
            Integer const id_1_end{ node_1.box_ptr + node_1.box_num };
            Integer const id_2_ini{ node_2.box_ptr                  };
            Integer const id_2_end{ node_2.box_ptr + node_2.box_num };
            if ( node_1.box_num < node_2.box_num ) {
              for ( Integer i{id_1_ini}; i < id_1_end; ++i ) {
                Integer const pos_1{ m_tree_boxes_map[i] };
                if ( ++m_check_counter; !boxes_1[pos_1]->intersects(node_2.box_long) ) continue;
                for ( Integer j{id_2_ini}; j < id_2_end; ++j ) {
                  Integer const pos_2{ tree.m_tree_boxes_map[j] };
                  if ( ++m_check_counter; boxes_1[pos_1]->intersects(*boxes_2[pos_2]) ) candidates[pos_1].insert(pos_2);
                }
              }
            } else {
              for ( Integer j{id_2_ini}; j < id_2_end; ++j ) {
                Integer const pos_2{ tree.m_tree_boxes_map[j] };
                if ( ++m_check_counter; !boxes_2[pos_2]->intersects(node_1.box_long) ) continue;
                for ( Integer i{id_1_ini}; i < id_1_end; ++i ) {
                  Integer const pos_1{ m_tree_boxes_map[i] };
                    if ( ++m_check_counter; boxes_1[pos_1]->intersects(*boxes_2[pos_2]) ) candidates[pos_1].insert(pos_2);
                }
              }
            }
          }
        }
        
        // controlla se entrambe foglie
        bool const leaf_1{ node_1.child_l < 0 };
        bool const leaf_2{ node_2.child_l < 0 };
        if ( leaf_1 && leaf_2 ) continue; // finita esplorazione questo ramo

        // Push children of both trees on the stack if they are not leafs
        auto stack_emplace_back = [this]( Integer const node_1, Integer const node_2 ) {
          m_stack.emplace_back(node_1);
          m_stack.emplace_back(node_2);
        };
        
        if ( leaf_1 ) {
          stack_emplace_back( id_1, node_2.child_l );
          stack_emplace_back( id_1, node_2.child_r );
        } else if ( leaf_2 ) {
          stack_emplace_back( node_1.child_l, id_2 );
          stack_emplace_back( node_1.child_r, id_2 );
        } else {
          if ( node_1.box_tot_num > node_2.box_tot_num ) { // split first larger tree
            if ( id_s1 >= 0 ) {
              stack_emplace_back( node_1.child_l, id_s2 );
              stack_emplace_back( node_1.child_r, id_s2 );
              if ( node_1.box_num > 0 ) stack_emplace_back( negate(id_1), id_2 );
            } else if ( id_s2 >= 0 ) { // and id_s1 < 0
              stack_emplace_back( id_s1, node_2.child_l );
              stack_emplace_back( id_s1, node_2.child_r );
            }
          } else {
            if ( id_s2 >= 0 ) {
              stack_emplace_back( id_s1, node_2.child_l );
              stack_emplace_back( id_s1, node_2.child_r );
              if ( node_2.box_num > 0 ) stack_emplace_back( id_1, negate(id_2) );
            } else if ( id_s1 >= 0 ) { // and id_s2 < 2
              stack_emplace_back( node_1.child_l, id_s2 );
              stack_emplace_back( node_1.child_r, id_s2 );
            }
          }
        }
      }

      // Return true if the trees intersect
      return !candidates.empty();
    }

    /**
     * \brief Self-intersect the tree (i.e., intersect the tree with itself to find all the intersecting boxes).
     * \param[out] candidates Intersection result (boxes indexes).
     * \return True if the tree intersects itself, false otherwise.
     */
    bool
    self_intersect( IndexSet & candidates ) const {
      IndexMap candidates_map;
      bool intersects{this->intersect(*this, candidates_map)};
      candidates.clear();
      for (const auto & [key, values] : candidates_map) {
        for (int value : values) {candidates.emplace(key); candidates.emplace(value);}
      }
      return intersects;
    }

    /**
     * \brief Compute the minimum distance between an object and the tree.
     * \param[in] obj Object to compute the distance to.
     * \param[out] candidates Minimum distance candidates.
     * \return The minimum distance between the object and the tree.
     * \tparam Object Type of the object to compute the distance to.
     * \note Object must have a method \c interior_distance that computes the distance to a box.
     */
    template <typename Object>
    Real distance( Object const & obj, IndexSet & candidates ) const {
      // Reset statistics
      m_check_counter = 0;

      // Return a negative value if the tree is empty
      if ( this->is_empty() || this->is_empty() ) return -1;

      // Collect the original object boxes
      BoxUniquePtrList const & boxes{ *m_boxes_ptr };

      // Setup the stack
      m_stack.clear();
      m_stack.reserve(2*this->size() + 1);
      m_stack.emplace_back(0);
      m_stack.emplace_back(0);

      // Main loop that checks the intersection iteratively
      Real distance{std::numeric_limits<Real>::max()};
      candidates.clear();
      while (!m_stack.empty())
      {
        // Pop the node from stack
        Integer const id{ m_stack.back() }; m_stack.pop_back();

        // Get the node
        Node const & node { m_tree_structure[id] };

        // Compute the distance between the object and the bounding box
        ++m_check_counter;
        Real tmp_distance{node.box.interior_distance(obj)};

        // If the distance is greater than the temporary minimum distance, skip the node
        if (tmp_distance > distance) {continue;}

        // Compute the distance between the object and the long boxes on the node
        if (node.box_num > 0) {
          
          if ( ++m_check_counter; node.box_long.interior_distance(obj) <= distance) {
            Integer const id_ini{ node.box_ptr                };
            Integer const id_end{ node.box_ptr + node.box_num };
            for (Integer i{id_ini}; i < id_end; ++i) {
              Integer const pos{ m_tree_boxes_map[i] };
              ++m_check_counter;
              tmp_distance = boxes[pos]->interior_distance(obj);
              if (tmp_distance < distance) {
                candidates.clear(); candidates.insert(pos); distance = tmp_distance;
              } else if (tmp_distance == distance) {
                candidates.insert(pos);
              }
            }
          }
        }

        // Push children on the stack if thay are not leafs
        if ( node.child_l > 0 ) m_stack.emplace_back(node.child_l);
        if ( node.child_r > 0 ) m_stack.emplace_back(node.child_r);
      }

      // Return the distance between the point and the tree
      return distance;
    }

    /**
     * \brief Compute the minimum distance between the current tree and another tree.
     * \param[in] tree Tree to compute the distance to.
     * \param[out] candidates Minimum distance candidates.
     * \return The minimum distance between the trees.
     */
    Real distance( Tree const & tree, IndexMap & candidates ) const {
      // Reset statistics
      m_check_counter = 0;

      // Return if the tree is empty
      if ( this->is_empty() || tree.is_empty() ) return -1;

      // Collect the original object boxes
      BoxUniquePtrList const & boxes_1{ *m_boxes_ptr      };
      BoxUniquePtrList const & boxes_2{ *tree.m_boxes_ptr };

      // Setup the stack
      m_stack.clear();
      m_stack.reserve( this->size() + tree.size() + 2 );
      m_stack.emplace_back(0);
      m_stack.emplace_back(0);

      // Negate the id of the node to distinguish the tree
      auto negate = [] (Integer const id) { return -1-id; };

      // Main loop that checks the intersection iteratively
      Real distance{std::numeric_limits<Real>::max()};
      candidates.clear();
      while (!m_stack.empty())
      {
        // Pop the node from stack (reversed order)
        Integer const id_s2 { m_stack.back() }; m_stack.pop_back();
        Integer const id_2  { id_s2 >= 0 ? id_s2 : negate(id_s2)};
        Integer const id_s1 { m_stack.back() }; m_stack.pop_back();
        Integer const id_1  { id_s1 >= 0 ? id_s1 : negate(id_s1)};

        // Get the node
        Node const & node_1{ m_tree_structure[id_1] };
        Node const & node_2{ tree.m_tree_structure[id_2] };

        // Compute the distance between the bounding boxes
        ++m_check_counter;
        Real tmp_distance{node_1.box.interior_distance(node_2.box)};

        // If the distance is greater than the temporary minimum distance, skip the nodes
        if (tmp_distance > distance) {continue;}

        // Compute the distance between the long boxes on the nodes
        if (node_1.box_num > 0 && node_2.box_num > 0) {
          if (++m_check_counter; node_1.box_long.interior_distance(node_2.box_long) <= distance) {
            Integer const id_1_ini{ node_1.box_ptr                  };
            Integer const id_1_end{ node_1.box_ptr + node_1.box_num };
            Integer const id_2_ini{ node_2.box_ptr                  };
            Integer const id_2_end{ node_2.box_ptr + node_2.box_num };
            for (Integer i{id_1_ini}; i < id_1_end; ++i) {
              Integer const pos_1{ m_tree_boxes_map[i] };
              for (Integer j{id_2_ini}; j < id_2_end; ++j) {
                Integer const pos_2{ tree.m_tree_boxes_map[j] };
                ++m_check_counter;
                tmp_distance = boxes_1[pos_1]->interior_distance(*boxes_2[pos_2]);
                if (tmp_distance < distance) {
                  candidates.clear();
                  candidates[pos_1].insert(pos_2);
                  distance = tmp_distance;
                } else if (tmp_distance == distance) {
                  candidates[pos_1].insert(pos_2);
                }
              }
            }
          }
        }

        // Push children of both trees on the stack if thay are not leafs
        auto stack_emplace_back = [this](Integer const node_1, Integer const node_2) {
          m_stack.emplace_back(node_1); m_stack.emplace_back(node_2);
        };
        if (id_s1 >= 0) {
          stack_emplace_back(node_1.child_l, id_s2);
          stack_emplace_back(node_1.child_r, id_s2);
          if (node_1.box_num > 0) {stack_emplace_back(negate(id_1), id_2);}
        } else if (id_s2 >= 0) {
          stack_emplace_back(id_s1, node_2.child_l);
          stack_emplace_back(id_s1, node_2.child_r);
        }
      }

      // Return the distance between the trees
      return distance;
    }

    /**
    * Find the first \f$ n \f$ candidates from an object, using a custom distance function.
    * \param[in] obj Object to compute the distance to.
    * \param[in] n Number of candidates to find.
    * \param[out] candidates First \f$ n \f$ candidates.
    * \param[in] distance_function Custom distance function (default is the interior distance).
    * \return The maximum distance between the candidates and the point. The return value is negative
    * if no candidates are found.
    * \tparam Object Type of the object to compute the distance to.
    * \tparam Function Type of the custom distance function.
    */
    template <typename Object, typename Function = std::function<Real(Object const &, Box const &)>>
    Real
    closest(
      Object  const & obj,
      Integer const   n,
      IndexSet      & candidates,
      Function        distance_function = [] (Object const & o, Box const & b) { return b.interior_distance(o); }
    ) const {
      // Reset statistics
      m_check_counter = 0;

      // Return if the tree is empty
      if ( this->is_empty() ) return -1;

      // Collect the original object boxes
      BoxUniquePtrList const & boxes{ *m_boxes_ptr };

      // Setup the stack
      m_stack.clear();
      m_stack.reserve(2*this->size() + 1);
      m_stack.emplace_back(0);

      // Candidate vector distance and index
      using Pair = std::pair<Real, Integer>;
      auto cmp = [](const Pair & a, const Pair & b) {return a.first < b.first;};
      std::priority_queue<Pair, std::vector<Pair>, decltype(cmp)> queue(cmp);

      // Main loop that checks the intersection iteratively
      Real max_distance{std::numeric_limits<Real>::max()}; // Maximum distance in the queue
      candidates.clear();
      while (!m_stack.empty())
      {
        // Pop the node from stack
        Integer const id{ m_stack.back() }; m_stack.pop_back();

        // Get the node
        Node const & node{ m_tree_structure[id] };

        // Compute the distance between the object and the bounding box
        ++m_check_counter;
        Real tmp_distance{distance_function(obj, node.box)};

        // If the distance is greater than the maximum distance, skip the node
        if (tmp_distance > max_distance) {continue;}

        // Compute the distance between the object and the long boxes on the node
        if (node.box_num > 0) {
          Integer const id_ini{ node.box_ptr                };
          Integer const id_end{ node.box_ptr + node.box_num };
          for (Integer i{id_ini}; i < id_end; ++i) {
            Integer const pos{ m_tree_boxes_map[i] };
            ++m_check_counter;
            tmp_distance = distance_function(obj, *boxes[pos]);
            if (tmp_distance < max_distance) {
              if (static_cast<Integer>(queue.size()) < n) {
                queue.emplace(tmp_distance, pos);
              } else {
                queue.pop(); queue.emplace(tmp_distance, pos);
              }
              max_distance = queue.top().first;
            }
          }
        }

        // Push children on the stack if thay are not leafs
        if ( node.child_l > 0 ) m_stack.emplace_back(node.child_l);
        if ( node.child_r > 0 ) m_stack.emplace_back(node.child_r);
      }

      // Extract indices into candidates
      Real min_distance{queue.empty() ? max_distance : queue.top().first};
      while (!queue.empty()) {candidates.insert(queue.top().second); queue.pop();}
      return min_distance;
    }

    /**
     * \brief Find the candidates that are within a given distance from a object, using a custom distance function.
     * \param[in] obj Object to compute the distance to.
     * \param[in] max_distance Maximum distance to consider.
     * \param[out] candidates Minimum distance candidates.
     * \param[in] distance_function Custom distance function (default is the interior distance).
     * \return True if at least one object is within the given distance, false otherwise.
     * \tparam Object Type of the object to compute the distance to.
     * \tparam Function Type of the custom distance function.
     */
    template <typename Object, typename Function = std::function<Real(Object const &, Box const &)>>
    bool
    within_distance(
      Object const & obj,
      Real   const   max_distance,
      IndexSet     & candidates,
      Function       distance_function = [] ( Object const & o, Box const & b ) { return b.interior_distance(o); }
    ) const {
      // Reset statistics
      m_check_counter = 0;

      // Return if the tree is empty
      if ( this->is_empty() ) return false;

      // Collect the original object boxes
      BoxUniquePtrList const & boxes{ *m_boxes_ptr };

      // Setup the stack
      m_stack.clear();
      m_stack.reserve(2*this->size() + 1);
      m_stack.emplace_back(0);

      // Main loop that checks the intersection iteratively
      candidates.clear();
      while (!m_stack.empty())
      {
        // Pop the node from stack
        Integer const id{m_stack.back()}; m_stack.pop_back();

        // Get the node
        Node const & node{m_tree_structure[id]};

        // Compute the distance between the object and the bounding box
        ++m_check_counter;
        Real distance{distance_function(obj, node.box)};

        // If the distance is greater than the maximum distance, skip the node
        if (distance > max_distance) {continue;}

        // Compute the distance between the object and the long boxes on the node
        if (node.box_num > 0) {
          Integer const id_ini{ node.box_ptr                };
          Integer const id_end{ node.box_ptr + node.box_num };
          for (Integer i{id_ini}; i < id_end; ++i) {
            Integer const pos{ m_tree_boxes_map[i] };
            ++m_check_counter;
            Real tmp_distance{distance_function(obj, *boxes[pos])};
            if (tmp_distance <= max_distance) {candidates.insert(pos);}
          }
        }

        // Push children on the stack if thay are not leafs
        if ( node.child_l > 0 ) m_stack.emplace_back(node.child_l);
        if ( node.child_r > 0 ) m_stack.emplace_back(node.child_r);
      }

      // Return true if at least one candidate is within the given distance
      return !candidates.empty();
    }

    /**
     * \brief Compute the depth of the tree from a given node.
     * \param[in] i Index of the node to compute the depth from.
     * \param[out] d Depth of the tree from the given node.
     */
    void depth( Integer const i, Integer & d ) const {
      d = 0;
      if (i < 0) return;
      m_stack.clear();
      m_stack.reserve(2*this->size() + 1);
      m_stack.emplace_back(i);
      std::vector<Integer> depth_stack;
      depth_stack.reserve(2*this->size() + 1);
      depth_stack.emplace_back(0);
      Integer depth{0};
      while (!m_stack.empty())
      {
        Integer const id{ m_stack.back() }; m_stack.pop_back();
        Node const & node { m_tree_structure[id] };
        depth = static_cast<Integer>(depth_stack.back()); depth_stack.pop_back();
        if (node.child_l == -1) { d = std::max(d, depth); }
        else                    { m_stack.emplace_back(node.child_l); depth_stack.emplace_back(depth+1); }
        if (node.child_r == -1) { d = std::max(d, depth); }
        else                    { m_stack.emplace_back(node.child_r); depth_stack.emplace_back(depth+1); }
      }
    }

    /**
     * \brief Compute the number of leafs, nodes, and long boxes of the tree from a given node.
     * \param[in] i Index of the node to start the computation from.
     * \param[out] l Number of leafs of the tree from the given node.
     * \param[out] n Number of nodes of the tree from the given node.
     * \param[out] b Number of long boxes of the tree from the given node.
     */
    void nodes( Integer const i, Integer & l, Integer & n, Integer & b ) const {
      l = n = b = 0;
      if (i < 0) {return;}
      m_stack.clear();
      m_stack.reserve(2*this->size() + 1);
      m_stack.emplace_back(i);
      while (!m_stack.empty())
      {
        Integer const id{ m_stack.back() }; m_stack.pop_back();
        Node const & node { m_tree_structure[id] };
        if ( node.child_l == -1 ) ++l; else m_stack.emplace_back(node.child_l);
        if ( node.child_r == -1 ) ++l; else m_stack.emplace_back(node.child_r);
        ++n; b += node.box_num;
      }
    }

    /**
     * \brief Compute some statistics about the current tree.
     * \param[out] stats Statistics about the tree.
     */
    void
    stats( Statistics & stats ) const {
      // Reset statistics
      stats.reset();

      // Compute/copy the build statistics
      stats.objects = static_cast<Integer>(m_boxes_ptr->size());
      this->nodes( 0, stats.leafs, stats.nodes, stats.long_boxes );
      this->depth( 0, stats.depth );
      this->nodes( m_tree_structure[0].child_l, stats.left_leafs, stats.left_nodes, stats.left_long_boxes );
      this->depth( m_tree_structure[0].child_l, stats.left_depth );
      this->nodes( m_tree_structure[0].child_r, stats.right_leafs, stats.right_nodes, stats.right_long_boxes );
      this->depth( m_tree_structure[0].child_r, stats.right_depth );
      stats.dump_counter  = m_dump_counter;
      stats.balance_ratio = static_cast<Real>(stats.left_leafs)/static_cast<Real>(stats.right_leafs);
      stats.depth_ratio   = static_cast<Real>(stats.left_depth)/static_cast<Real>(stats.right_depth);

      // Copy the check counter
      stats.check_counter = m_check_counter;
    }

    /**
     * \brief Print the tree info to an output stream.
     * \param[in] os Output stream to print the tree info to.
     */
    void print( OutStream & os ) const  {
      // Retrieve the statistics
      Statistics stats; this->stats(stats);
      // Print the tree info
      stats.print( os );
    }

  }; // Tree

} // namespace AABBtree

#endif // INCLUDE_AABBTREE_TREE_HXX
