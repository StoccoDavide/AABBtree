/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef INCLUDE_AABBTREE_HXX
#define INCLUDE_AABBTREE_HXX

namespace AABBtree {

  /**
  * \brief A class representing an axis-aligned bounding box tree (AABBtree).
  *
  * The AABBtree class provides an efficient way to store and query a set of axis-aligned bounding
  * boxes. It supports various operations such as building the tree, adding axis-aligned bounding
  * boxes, and performing intersection queries. The AABBtree class is templated on the type of the
  * real numbers and the dimensions of the bounding boxes.
  * \tparam Real Type of the real numbers.
  * \tparam N Dimensions of the space.
  */
  template <typename Real, Integer N>
  class AABBtree {
  public:
    // Basic types for the AABBtree
    using Point = Eigen::Vector<Real, N>; /**< Point in N-dimensional space (Eigen column vector). */
    using PointPtr = std::unique_ptr<Point>; /**< Unique pointer to a point. */
    using AlignedBox = AlignedBox<Real, 2>; /**< Axis-aligned bounding box in N-dimensional space. */
    using AlignedBoxPtr = std::unique_ptr<AlignedBox>; /**< Unique pointer to an axis-aligned bounding box. */

    // Set and map types for retrieving indexes
    using Set = std::set<Integer>; /**< Set of indexes. */
    using Map = std::map<Integer, Set>; /**< Map of indexes. */

  private:
    Real *
    std::vector<Real> m_rmem; /**< Memory allocator for reals. */
    std::vector<Integer> m_imem; /**< Memory allocator for integers. */

    // Structure
    RealPtr m_tree_bbox; /**< Tree bounding boxes (size: N*m_2dim). */
    RealPtr m_objs_bbox; /**< Object bounding boxes (size: m_num_objects*m_2dim). */

    SizePtr m_father; /**< Father nodes (size: N). */
    SizePtr m_child; /**< Child nodes (size: N). */
    SizePtr m_ptr_nodes; /**< Pointer to nodes (size: N). */
    SizePtr m_num_nodes; /**< Number of nodes (size: N). */
    SizePtr m_id_nodes; /**< Nodes id (size: m_num_objects). */

    mutable std::vector<Integer> m_stack; /**< Tree stack. */

    // Parameters
    Integer m_max_nodal_objects{16}; /**< Maximum number of objects per node. */
    Real m_long_edge_ratio{0.8}; /**< Long edge ratio for bounding boxes. */
    Real m_collision_tolerance{0.1}; /**< Overlap tolerance for bounding boxes. */
    Real m_min_size_tolerance{0}; /**< Minimum size tolerance for bounding boxes. */

    // Statistics
    mutable Integer m_check_counter{0}; /**< Number of overlap check (for statistic). */

  public:
    /**
    * Class constructor for the axis-aligned bounding box tree.
    */
    AABBtree(void) {};

    /**
    * Class constructor for the axis-aligned bounding box tree.
    * \param[in] tree Axis-aligned bounding box tree to copy.
    */
    AABBtree(AABBtree<Real, N> const & tree)
    {

      this->m_max_nodal_objects   = tree.m_max_nodal_objects;
      this->m_long_edge_ratio     = tree.m_long_edge_ratio;
      this->m_collision_tolerance = tree.m_collision_tolerance;
      this->m_min_size_tolerance  = tree.m_min_size_tolerance;
    }

    /**
    * Collision detection function between axis-aligned bounding boxes.
    * \param[in] bbox_1 First axis-aligned bounding box.
    * \param[in] bbox_2 Second axis-aligned bounding box.
    * \return True if the axis-aligned bounding boxes collide, false otherwise.
    */
    bool AABBtree::do_collide(BBox const bbox_1, BBox const bbox_2);

    /**
    * Overlap detection function between an axis-aligned bounding box and a point.
    * \param[in] bbox Axis-aligned bounding box.
    * \param[in] pnt Point.
    * \return True if the axis-aligned bounding box and the point overlap, false otherwise.
    */
    bool AABBtree::do_collide(BBox const bbox, Point const pnt);


  }; // AABBtree


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename Real>
  bool AABBtree<Real, 1>::do_collide(BBox const bbox_1, BBox const bbox_2)
  {
    return bbox_1[0] <= bbox_2[1] && bbox_1[1] >= bbox_2[0];
  }

  template <typename Real>
  bool AABBtree<Real, 2>::do_collide(BBox const bbox_1, BBox const bbox_2)
  {
    return bbox_1[0] <= bbox_2[2] && bbox_1[2] >= bbox_2[0] &&
           bbox_1[1] <= bbox_2[3] && bbox_1[3] >= bbox_2[1];
  }

  template <typename Real>
  bool AABBtree<Real, 3>::do_collide(BBox const bbox_1, BBox const bbox_2)
  {
    return bbox_1[0] <= bbox_2[3] && bbox_1[3] >= bbox_2[0] &&
           bbox_1[1] <= bbox_2[4] && bbox_1[4] >= bbox_2[1] &&
           bbox_1[2] <= bbox_2[5] && bbox_1[5] >= bbox_2[2];
  }

  template <typename Real>
  bool AABBtree<Real, 4>::do_collide(BBox const bbox_1, BBox const bbox_2)
  {
    return bbox_1[0] <= bbox_2[4] && bbox_1[4] >= bbox_2[0] &&
           bbox_1[1] <= bbox_2[5] && bbox_1[5] >= bbox_2[1] &&
           bbox_1[2] <= bbox_2[6] && bbox_1[6] >= bbox_2[2] &&
           bbox_1[3] <= bbox_2[7] && bbox_1[7] >= bbox_2[3];
  }

  template <typename Real>
  bool AABBtree<Real, 5>::do_collide(BBox const bbox_1, BBox const bbox_2)
  {
    return bbox_1[0] <= bbox_2[5] && bbox_1[5] >= bbox_2[0] &&
           bbox_1[1] <= bbox_2[6] && bbox_1[6] >= bbox_2[1] &&
           bbox_1[2] <= bbox_2[7] && bbox_1[7] >= bbox_2[2] &&
           bbox_1[3] <= bbox_2[8] && bbox_1[8] >= bbox_2[3] &&
           bbox_1[4] <= bbox_2[9] && bbox_1[9] >= bbox_2[4];
  }

  template <typename Real>
  bool AABBtree<Real, 6>::do_collide(BBox const bbox_1, BBox const bbox_2)
  {
    return bbox_1[0] <= bbox_2[6]  && bbox_1[6]  >= bbox_2[0] &&
           bbox_1[1] <= bbox_2[7]  && bbox_1[7]  >= bbox_2[1] &&
           bbox_1[2] <= bbox_2[8]  && bbox_1[8]  >= bbox_2[2] &&
           bbox_1[3] <= bbox_2[9]  && bbox_1[9]  >= bbox_2[3] &&
           bbox_1[4] <= bbox_2[10] && bbox_1[10] >= bbox_2[4] &&
           bbox_1[5] <= bbox_2[11] && bbox_1[11] >= bbox_2[5];
  }

  template <typename Real>
  bool AABBtree<Real, 7>::do_collide(BBox const bbox_1, BBox const bbox_2)
  {
    return bbox_1[0] <= bbox_2[7]  && bbox_1[7]  >= bbox_2[0] &&
           bbox_1[1] <= bbox_2[8]  && bbox_1[8]  >= bbox_2[1] &&
           bbox_1[2] <= bbox_2[9]  && bbox_1[9]  >= bbox_2[2] &&
           bbox_1[3] <= bbox_2[10] && bbox_1[10] >= bbox_2[3] &&
           bbox_1[4] <= bbox_2[11] && bbox_1[11] >= bbox_2[4] &&
           bbox_1[5] <= bbox_2[12] && bbox_1[12] >= bbox_2[5] &&
           bbox_1[6] <= bbox_2[13] && bbox_1[14] >= bbox_2[6];
  }

  template <typename Real>
  bool AABBtree<Real, 8>::do_collide(BBox const bbox_1, BBox const bbox_2)
  {
    return bbox_1[0] <= bbox_2[8]  && bbox_1[8]  >= bbox_2[0] &&
           bbox_1[1] <= bbox_2[9]  && bbox_1[9]  >= bbox_2[1] &&
           bbox_1[2] <= bbox_2[10] && bbox_1[10] >= bbox_2[2] &&
           bbox_1[3] <= bbox_2[11] && bbox_1[11] >= bbox_2[3] &&
           bbox_1[4] <= bbox_2[12] && bbox_1[12] >= bbox_2[4] &&
           bbox_1[5] <= bbox_2[13] && bbox_1[13] >= bbox_2[5] &&
           bbox_1[6] <= bbox_2[14] && bbox_1[14] >= bbox_2[6] &&
           bbox_1[7] <= bbox_2[15] && bbox_1[15] >= bbox_2[7];
  }

  template <typename Real>
  bool AABBtree<Real, 8>::do_collide(BBox const bbox_1, BBox const bbox_2)
  {
    return bbox_1[0] <= bbox_2[8]  && bbox_1[8]  >= bbox_2[0] &&
           bbox_1[1] <= bbox_2[9]  && bbox_1[9]  >= bbox_2[1] &&
           bbox_1[2] <= bbox_2[10] && bbox_1[10] >= bbox_2[2] &&
           bbox_1[3] <= bbox_2[11] && bbox_1[11] >= bbox_2[3] &&
           bbox_1[4] <= bbox_2[12] && bbox_1[12] >= bbox_2[4] &&
           bbox_1[5] <= bbox_2[13] && bbox_1[13] >= bbox_2[5] &&
           bbox_1[6] <= bbox_2[14] && bbox_1[14] >= bbox_2[6] &&
           bbox_1[7] <= bbox_2[15] && bbox_1[15] >= bbox_2[7];
  }

  template <typename Real>
  bool AABBtree<Real, 9>::do_collide(BBox const bbox_1, BBox const bbox_2)
  {
    return bbox_1[0] <= bbox_2[9]  && bbox_1[9]  >= bbox_2[0] &&
           bbox_1[1] <= bbox_2[10] && bbox_1[10] >= bbox_2[1] &&
           bbox_1[2] <= bbox_2[11] && bbox_1[11] >= bbox_2[2] &&
           bbox_1[3] <= bbox_2[12] && bbox_1[12] >= bbox_2[3] &&
           bbox_1[4] <= bbox_2[13] && bbox_1[13] >= bbox_2[4] &&
           bbox_1[5] <= bbox_2[14] && bbox_1[14] >= bbox_2[5] &&
           bbox_1[6] <= bbox_2[15] && bbox_1[15] >= bbox_2[6] &&
           bbox_1[7] <= bbox_2[16] && bbox_1[16] >= bbox_2[7] &&
           bbox_1[8] <= bbox_2[17] && bbox_1[17] >= bbox_2[8];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename Real>
  bool AABBtree<Real, 1>::do_collide(BBox const bbox, Point const pnt)
  {
    return pnt[0] <= bbox[1] && pnt[0] >= bbox[0];
  }

  template <typename Real>
  bool AABBtree<Real, 2>::do_collide(BBox const bbox, Point const pnt)
  {
    return pnt[0] <= bbox[2] && pnt[0] >= bbox[0] &&
           pnt[1] <= bbox[3] && pnt[1] >= bbox[1];
  }

  template <typename Real>
  bool AABBtree<Real, 3>::do_collide(BBox const bbox, Point const pnt)
  {
    return pnt[0] <= bbox[3] && pnt[0] >= bbox[0] &&
           pnt[1] <= bbox[4] && pnt[1] >= bbox[1] &&
           pnt[2] <= bbox[5] && pnt[2] >= bbox[2];
  }

  template <typename Real>
  bool AABBtree<Real, 4>::do_collide(BBox const bbox, Point const pnt)
  {
    return pnt[0] <= bbox[4] && pnt[0] >= bbox[0] &&
           pnt[1] <= bbox[5] && pnt[1] >= bbox[1] &&
           pnt[2] <= bbox[6] && pnt[2] >= bbox[2] &&
           pnt[3] <= bbox[7] && pnt[3] >= bbox[3];
  }

  template <typename Real>
  bool AABBtree<Real, 5>::do_collide(BBox const bbox, Point const pnt)
  {
    return pnt[0] <= bbox[5] && pnt[0] >= bbox[0] &&
           pnt[1] <= bbox[6] && pnt[1] >= bbox[1] &&
           pnt[2] <= bbox[7] && pnt[2] >= bbox[2] &&
           pnt[3] <= bbox[8] && pnt[3] >= bbox[3] &&
           pnt[4] <= bbox[9] && pnt[4] >= bbox[4];
  }

  template <typename Real>
  bool AABBtree<Real, 6>::do_collide(BBox const bbox, Point const pnt)
  {
    return pnt[0] <= bbox[6]  && pnt[0] >= bbox[0] &&
           pnt[1] <= bbox[7]  && pnt[1] >= bbox[1] &&
           pnt[2] <= bbox[8]  && pnt[2] >= bbox[2] &&
           pnt[3] <= bbox[9]  && pnt[3] >= bbox[3] &&
           pnt[4] <= bbox[10] && pnt[4] >= bbox[4] &&
           pnt[5] <= bbox[11] && pnt[5] >= bbox[5];
  }

  template <typename Real>
  bool AABBtree<Real, 7>::do_collide(BBox const bbox, Point const pnt)
  {
    return pnt[0] <= bbox[7]  && pnt[0] >= bbox[0] &&
           pnt[1] <= bbox[8]  && pnt[1] >= bbox[1] &&
           pnt[2] <= bbox[9]  && pnt[2] >= bbox[2] &&
           pnt[3] <= bbox[10] && pnt[3] >= bbox[3] &&
           pnt[4] <= bbox[11] && pnt[4] >= bbox[4] &&
           pnt[5] <= bbox[12] && pnt[5] >= bbox[5] &&
           pnt[6] <= bbox[13] && pnt[6] >= bbox[6];
  }

  template <typename Real>
  bool AABBtree<Real, 8>::do_collide(BBox const bbox, Point const pnt)
  {
    return pnt[0] <= bbox[8]  && pnt[0] >= bbox[0] &&
           pnt[1] <= bbox[9]  && pnt[1] >= bbox[1] &&
           pnt[2] <= bbox[10] && pnt[2] >= bbox[2] &&
           pnt[3] <= bbox[11] && pnt[3] >= bbox[3] &&
           pnt[4] <= bbox[12] && pnt[4] >= bbox[4] &&
           pnt[5] <= bbox[13] && pnt[5] >= bbox[5] &&
           pnt[6] <= bbox[14] && pnt[6] >= bbox[6] &&
           pnt[7] <= bbox[15] && pnt[7] >= bbox[7];
  }

  template <typename Real>
  bool AABBtree<Real, 9>::do_collide(BBox const bbox, Point const pnt)
  {
    return pnt[0] <= bbox[9]  && pnt[0] >= bbox[0] &&
           pnt[1] <= bbox[10] && pnt[1] >= bbox[1] &&
           pnt[2] <= bbox[11] && pnt[2] >= bbox[2] &&
           pnt[3] <= bbox[12] && pnt[3] >= bbox[3] &&
           pnt[4] <= bbox[13] && pnt[4] >= bbox[4] &&
           pnt[5] <= bbox[14] && pnt[5] >= bbox[5] &&
           pnt[6] <= bbox[15] && pnt[6] >= bbox[6] &&
           pnt[7] <= bbox[16] && pnt[7] >= bbox[7] &&
           pnt[8] <= bbox[17] && pnt[8] >= bbox[8];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

} // namespace AABBtree

#endif // INCLUDE_AABBTREE_HXX
