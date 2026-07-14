/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The AABBtree project is distributed under the BSD 2-Clause License.       *
 *                                                                           *
 * Davide Stocco Enrico Bertolazzi                                           *
 * University of Trento University of Trento                                 *
 * davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// C++17 standard libraries
#include <cmath>
#include <memory>
#include <vector>

// AABBtree library
#include "AABBtree.hh"
using namespace AABBtree;

// Catch2 library
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
using namespace Catch::Matchers;

// Test utilities
#include "TestUtilities.hh"
using namespace TestUtilities;

namespace {

// Ambient space dimension and problem size used by this test. Kept small so
// that all the (max_dumpings, max_nodal_objects) combinations run fast.
constexpr Integer N_DIM{3};
constexpr Integer N_OBJECTS_1{6000};
constexpr Integer N_OBJECTS_2{6000};
constexpr unsigned SEED_1{123456u};
constexpr unsigned SEED_2{654321u};

// Build a reproducible set of random N_DIM-dimensional boxes. The domain and
// box sizes are chosen so that boxes are close to each other (and often
// overlapping), which is a meaningful setup for distance queries.
template <typename Real>
std::unique_ptr<BoxUniquePtrList<Real, N_DIM>> BuildBoxes(Integer n_objects,
                                                          unsigned seed) {
  using Vector = AABBtree::Vector<Real, N_DIM>;
  using Box = AABBtree::Box<Real, N_DIM>;

  std::vector<RandomBox<Real, N_DIM>> random_boxes{
      GenerateRandomBoxes<Real, N_DIM, Integer>(
          n_objects, seed, static_cast<Real>(0.0), static_cast<Real>(10.0),
          static_cast<Real>(0.5), static_cast<Real>(2.0))};

  auto boxes{std::make_unique<BoxUniquePtrList<Real, N_DIM>>()};
  boxes->reserve(n_objects);
  for (Integer i{0}; i < n_objects; ++i) {
    Vector box_min, box_max;
    for (Integer d{0}; d < N_DIM; ++d) {
      box_min(d) = random_boxes[i].min[d];
      box_max(d) = random_boxes[i].max[d];
    }
    boxes->push_back(std::make_unique<Box>(box_min, box_max));
    boxes->back()->reorder();
  }
  return boxes;
}

// Tolerance used to compare distances computed by trees built with different
// (max_dumpings, max_nodal_objects) parameters. It should not depend on the
// tree parameters, only on the floating point type used.
template <typename Real> Real distance_tolerance() {
  return static_cast<Real>(100.0) *
         std::sqrt(std::numeric_limits<Real>::epsilon());
}

} // namespace

TEMPLATE_TEST_CASE("Tree distance results are invariant to build parameters",
                   "[template]", float, double) {
  using Tree = AABBtree::Tree<TestType, N_DIM>;
  using Box = AABBtree::Box<TestType, N_DIM>;

  // Reference results computed with default build parameters
  Tree ref_tree_1, ref_tree_2;
  ref_tree_1.build(BuildBoxes<TestType>(N_OBJECTS_1, SEED_1));
  ref_tree_2.build(BuildBoxes<TestType>(N_OBJECTS_2, SEED_2));

  IndexMap ref_candidates_tt;
  TestType const ref_distance_tt{
      ref_tree_1.distance(ref_tree_2, ref_candidates_tt)};

  Box const ref_query_box{*ref_tree_2.box(0)};
  IndexSet ref_candidates_tb;
  TestType const ref_distance_tb{
      ref_tree_1.distance(ref_query_box, ref_candidates_tb)};

  TestType const tolerance{distance_tolerance<TestType>()};

  // Vary the tree build parameters and check that the distance results
  // (both the minimum distance and the set of candidates achieving it) do
  // not change.
  Integer const max_dumpings{GENERATE(1, 2, 3)};
  Integer const max_nodal_objects{GENERATE(1, 10, 100, 1000, 10000)};

  CAPTURE(max_dumpings, max_nodal_objects);

  Tree tree_1, tree_2;
  tree_1.max_dumpings(max_dumpings);
  tree_2.max_dumpings(max_dumpings);
  tree_1.max_nodal_objects(max_nodal_objects);
  tree_2.max_nodal_objects(max_nodal_objects);
  tree_1.build(BuildBoxes<TestType>(N_OBJECTS_1, SEED_1));
  tree_2.build(BuildBoxes<TestType>(N_OBJECTS_2, SEED_2));

  // Tree-tree distance
  IndexMap candidates_tt;
  TestType const distance_tt{tree_1.distance(tree_2, candidates_tt)};
  CHECK_THAT(distance_tt, WithinAbs(ref_distance_tt, tolerance));
  CHECK(candidates_tt == ref_candidates_tt);

  // Tree-box distance
  IndexSet candidates_tb;
  TestType const distance_tb{tree_1.distance(ref_query_box, candidates_tb)};
  CHECK_THAT(distance_tb, WithinAbs(ref_distance_tb, tolerance));
  CHECK(candidates_tb == ref_candidates_tb);
}
