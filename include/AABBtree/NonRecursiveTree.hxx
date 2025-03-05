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
  template <typename Real, Integer N>
  class NonRecursive : public Tree<Real,N,Recursive<Real,N>> {
  public:
    // Basic types for managing the axis-aligned bounding boxes

    using Box        = Box<Real,N>;         /**< Axis-aligned bounding box in N-dimensional space. */
    using BoxUPtr    = BoxUPtr<Real,N>;     /**< Unique pointer to an axis-aligned bounding box. */
    using BoxUPtrVec = BoxUPtrVec<Real,N>;  /**< Vector of unique pointers to an axis-aligned bounding box. */

    using Vector     = Vector<Real,N>;       /**< Point in the ambient space (Eigen column vector of real numbers). */
    using Point      = Point<Real,N>;       /**< Point in the ambient space (Eigen column vector of real numbers). */

    typedef struct AABB_sub_tree {
      Box     bbox;     //< bounding box of the subtree
      Box     bbox3;    //< bounding box of "middle" boxes
      Integer parent;   //< root node (bbox) of the subtree
      Integer child_L;  //< children of the root of the subtree
      Integer child_R;  //< children of the root of the subtree
      Integer bb_ptr;   //< position of the first bbox in the m_map_bbox reordering of bboxes
      Integer bb_num;   //
    } AABBnode;

  private:

  // Tree hierarchy
    vector<AABB_sub_tree> m_AABB_structure;
    Indexes               m_map_bbox;

    // Cache and statistics
    mutable vector<Integer> m_stack; /**< Tree stack. */

    using Tree<Real,N,Recursive<Real,N>>::m_boxes;
    using Tree<Real,N,Recursive<Real,N>>::m_check_counter;
    using Tree<Real,N,Recursive<Real,N>>::m_separation_ratio_tolerance;
    using Tree<Real,N,Recursive<Real,N>>::m_balance_ratio_LR;

    using Tree<Real,N,Recursive<Real,N>>::m_max_nodal_objects;
    using Tree<Real,N,Recursive<Real,N>>::m_collision_tolerance;
    using Tree<Real,N,Recursive<Real,N>>::m_min_size_tolerance;

  public:

    /**
    * Class destructor for the \em non-recursive axis-aligned bounding box tree.
    */
    ~NonRecursive() = default;

    /**
    * Class constructor for the \em non-recursive axis-aligned bounding box tree.
    */
    NonRecursive() = default;

    /*
                                              m_AABB_structure
                     ┌──────┐                  ┌──────┐
                     │      │                  │ CD1  │
                     │      │                  │ CD2  │
                     │      │                  │ CD1  │
    m_nodes_id   --> │ BBOX │    m_children -> │ CD2  │
                     │      │                  │      │
                     │      │                  │      │
                     └──────┘

    */

    void
    build_impl() {
      BoxUPtrVec const & boxes { *m_boxes };
      Integer nbox  { static_cast<Integer>(boxes.size()) };
      Integer depth { static_cast<Integer>(ceil( log2( static_cast<Real>(nbox) ) )) };
      Integer dump  { 0 };

      m_AABB_structure.clear();
      m_AABB_structure.reserve( static_cast<size_t>( 2*nbox) );

      // setup root node
      AABBnode n;
      n.parent  = -1;
      n.child_L = -1;
      n.child_R = -1;
      n.bb_ptr  = 0;
      n.bb_num  = 0;
      m_AABB_structure.emplace_back( n );

      // setup 
      n.bbox.set_empty();
      m_map_bbox.reserve( static_cast<size_t>( 2*depth ) );
      for ( auto const & b : boxes ) {
        n.bbox.extend( b );
        m_map_bbox.emplace_back( n.bb_num++ );
      }

      // main loop: divide nodes until all constraints satisfied
      m_stack.clear();
      m_stack.reserve(2*boxes.size()+1);
      m_stack.emplace_back(0);

      while ( !m_stack.empty() ) {

        // pop node from stack
        Integer const id{ m_stack.back() }; m_stack.pop_back();

        AABBnode & n { m_AABB_structure[id] };
  
        // if few bbox stop splitting
        if ( n.bb_num < m_max_nodal_objects ) continue;
  
        //
        // The boxes are at the position m_ptr_nodes[ id_parent .. id_parent + num - 1 ]
        //
        Vector  mx_vec;
        Integer ipos_vec[N];
        n.bbox.sort_axis( mx_vec, ipos_vec );
        Integer idim                 { ipos_vec[dump] };
        Real    mx                   { mx_vec[idim] };
        Real    separation_line      { n.bbox.center(idim) };
        Real    separation_tolerance { mx * m_separation_ratio_tolerance };
  
        /*
        //                              _
        //  ___ ___ _ __  __ _ _ _ __ _| |_ ___
        // (_-</ -_) '_ \/ _` | '_/ _` |  _/ -_)
        // /__/\___| .__/\__,_|_| \__,_|\__\___|
        //         |_|
        */
        // separate short/long and accumulate short baricenter
        //
        // classified left and right are pushed to the end of the list
        //
        Integer n_long     { 0 };
        Integer n_left     { 0 };
        Integer n_right    { 0 };
        Integer ib         { n.bb_ptr};
        Integer ie         { n.bb_ptr+n.bb_num };
        Real    baricenter { 0 };
        while ( ib < ie ) {
          Box const & bb { *boxes[ m_map_bbox[ib] ] };
          Integer     ic { bb.classify( separation_line, separation_tolerance, idim  )};
          switch ( ic ) {
            case -1: ++n_left;  swap( m_map_bbox[ib], m_map_bbox[--ie] ); break;
            case +1: ++n_right; swap( m_map_bbox[ib], m_map_bbox[--ie] ); break;
            default: ++n_long; ++ib;
          }
          baricenter += bb.max(idim)+bb.min(idim);

        }
        baricenter /= 2*n.bb_num;

        // check if balanced, if not try to separate boxes on a new separation line
        if ( n_long > m_max_nodal_objects ||
             std::min(n_left,n_right) < m_balance_ratio_LR * std::max(n_left,n_right) ) {

          // try computation only if baricenter is well separated from separation_line
          if ( abs(baricenter-separation_line) > separation_tolerance ) {
            n_long = n_left = n_right = 0;
            ib     = n.bb_ptr;
            ie     = n.bb_ptr+n.bb_num;
            while ( ib < ie ) {
              Box const & bb { *boxes[ m_map_bbox[ib] ] };
              Integer     ic { bb.classify( baricenter, separation_tolerance, idim  )};
              switch ( ic ) {
                case -1: ++n_left;  swap( m_map_bbox[ib], m_map_bbox[--ie] ); break;
                case +1: ++n_right; swap( m_map_bbox[ib], m_map_bbox[--ie] ); break;
                default: ++n_long; ++ib;
              }
            }
          }
        }

        if ( n_long > m_max_nodal_objects &&
             std::min(n_left,n_right) < m_balance_ratio_LR * std::max(n_left,n_right) ) {
          if ( ++dump < N ) m_stack.push_back(id); // try in another direction
          continue;
        }

        // separe left e right
        n_left = n_right = 0;
        ib = n.bb_ptr + n_long;
        ie = n.bb_ptr + n.bb_num;
        while ( ib < ie ) {
          Integer ipos{ m_map_bbox[ib] };
          if ( boxes[ipos]->center(idim) < separation_line ) {
            ++ib; ++n_left; // in right position do nothing
          } else {
            --ie; ++n_right; swap( m_map_bbox[ib], m_map_bbox[ie] );
          }
        }

        // if cannot improve bbox, stop split at this level!
        if ( n_left <= m_max_nodal_objects && n_right <= m_max_nodal_objects ) continue;

        // child indexing
        n.child_L = static_cast<Integer>( m_AABB_structure.size()+0 ); 
        n.child_R = static_cast<Integer>( m_AABB_structure.size()+1 ); 

        // setup root node
        AABBnode L, R;
        L.parent  = id;
        L.child_L = -1;
        L.child_R = -1;
        L.bb_num  = n_left;

        R.parent  = id;
        R.child_L = -1;
        R.child_R = -1;
        R.bb_num  = n_right;

        // compute bbox of left and right child
        Integer ii{n.bb_ptr};
        n.bbox3.set_empty();
        n.bb_num = n_long;
        for ( Integer i{0}; i < n_long; ++i ) n.bbox3.extend( *boxes[m_map_bbox[ii++]] );
        L.bbox.set_empty(); L.bb_ptr = ii;
        for ( Integer i{0}; i < n_left; ++i ) L.bbox.extend( *boxes[m_map_bbox[ii++]] );
        R.bbox.set_empty(); R.bb_ptr = ii;
        for ( Integer i{0}; i < n_right; ++i ) R.bbox.extend( *boxes[m_map_bbox[ii++]] );

        m_AABB_structure.emplace_back( L );
        m_AABB_structure.emplace_back( R );

        m_stack.emplace_back( n.child_L ); 
        m_stack.emplace_back( n.child_R );
      }

    }

  }; // NonRecursive

  template class NonRecursive<float,1>;
  template class NonRecursive<float,2>;
  template class NonRecursive<float,3>;
  template class NonRecursive<float,4>;
  template class NonRecursive<float,5>;
  template class NonRecursive<float,6>;
  template class NonRecursive<float,7>;
  template class NonRecursive<float,8>;
  template class NonRecursive<float,9>;
  template class NonRecursive<float,10>;

  template class NonRecursive<double,1>;
  template class NonRecursive<double,2>;
  template class NonRecursive<double,3>;
  template class NonRecursive<double,4>;
  template class NonRecursive<double,5>;
  template class NonRecursive<double,6>;
  template class NonRecursive<double,7>;
  template class NonRecursive<double,8>;
  template class NonRecursive<double,9>;
  template class NonRecursive<double,10>;

} // namespace AABBtree

#endif // INCLUDE_AABBTREE_NONRECURSIVE_HXX
