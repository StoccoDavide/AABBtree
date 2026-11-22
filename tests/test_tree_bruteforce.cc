/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The AABBtree project is distributed under the BSD 2-Clause License.       *
 *                                                                           *
 * Davide Stocco Enrico Bertolazzi                                           *
 * University of Trento University of Trento                                 *
 * davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Compare every tree query against a brute-force O(n^2) reference.

// C++17 standard libraries
#include <algorithm>
#include <cmath>
#include <memory>
#include <thread>
#include <vector>

// AABBtree library (force the parallel tree-tree path also for small trees)
#define AABBTREE_PARALLEL_MIN_NODES 8
#include "AABBtree.hh"
using namespace AABBtree;

// Catch2 library
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

// Test utilities
#include "TestUtilities.hh"
using namespace TestUtilities;

namespace {

constexpr Integer N_DIM{3};
constexpr Integer N_OBJECTS_1{1500};
constexpr Integer N_OBJECTS_2{1200};
constexpr Integer N_QUERIES{40};

template <typename Real>
std::unique_ptr<BoxUniquePtrList<Real, N_DIM>> BuildBoxes(Integer n_objects,
                                                          unsigned seed) {
  using Vector = AABBtree::Vector<Real, N_DIM>;
  using Box = AABBtree::Box<Real, N_DIM>;
  std::vector<RandomBox<Real, N_DIM>> random_boxes{
      GenerateRandomBoxes<Real, N_DIM, Integer>(
          n_objects, seed, static_cast<Real>(0.0), static_cast<Real>(30.0),
          static_cast<Real>(0.1), static_cast<Real>(2.0))};
  auto boxes{std::make_unique<BoxUniquePtrList<Real, N_DIM>>()};
  boxes->reserve(n_objects);
  for (Integer i{0}; i < n_objects; ++i) {
    Vector box_min, box_max;
    for (Integer d{0}; d < N_DIM; ++d) {
      box_min(d) = random_boxes[i].min[d];
      box_max(d) = random_boxes[i].max[d];
    }
    boxes->push_back(std::make_unique<Box>(box_min, box_max));
  }
  return boxes;
}

template <typename Real>
using BoxList = std::vector<AABBtree::Box<Real, N_DIM>>;

template <typename Real>
BoxList<Real> CopyBoxes(BoxUniquePtrList<Real, N_DIM> const &boxes) {
  BoxList<Real> out;
  out.reserve(boxes.size());
  for (auto const &b : boxes) {
    out.push_back(*b);
  }
  return out;
}

} // namespace

TEMPLATE_TEST_CASE("Tree queries match brute force", "[template]", float,
                   double) {
  using Real = TestType;
  using Tree = AABBtree::Tree<Real, N_DIM>;
  using Box = AABBtree::Box<Real, N_DIM>;
  using Point = AABBtree::Point<Real, N_DIM>;

  Integer const max_dumpings{GENERATE(1, 3)};
  Integer const max_nodal_objects{GENERATE(1, 7, 50, 5000)};
  Integer const num_threads{GENERATE(1, 0)};
  CAPTURE(max_dumpings, max_nodal_objects, num_threads);

  auto boxes_1_ptr{BuildBoxes<Real>(N_OBJECTS_1, 1234u)};
  auto boxes_2_ptr{BuildBoxes<Real>(N_OBJECTS_2, 4321u)};
  BoxList<Real> const b1{CopyBoxes<Real>(*boxes_1_ptr)};
  BoxList<Real> const b2{CopyBoxes<Real>(*boxes_2_ptr)};
  Integer const n1{static_cast<Integer>(b1.size())};
  Integer const n2{static_cast<Integer>(b2.size())};

  Tree tree_1, tree_2;
  tree_1.max_dumpings(max_dumpings);
  tree_1.max_nodal_objects(max_nodal_objects);
  tree_2.max_dumpings(max_dumpings);
  tree_2.max_nodal_objects(max_nodal_objects);
  tree_1.num_threads(num_threads);
  tree_1.build(std::move(boxes_1_ptr));
  tree_2.build(std::move(boxes_2_ptr));

  // Query objects: some boxes of the second set and some points
  std::vector<Box> query_boxes;
  std::vector<Point> query_points;
  for (Integer q{0}; q < N_QUERIES; ++q) {
    query_boxes.push_back(b2[(q * 29) % n2]);
    query_points.push_back(b2[(q * 31) % n2].baricenter());
  }

  SECTION("Object intersection") {
    for (Box const &q : query_boxes) {
      IndexSet expected;
      for (Integer i{0}; i < n1; ++i) {
        if (b1[i].intersects(q)) {
          expected.insert(i);
        }
      }
      IndexSet candidates;
      bool const hit{tree_1.intersect(q, candidates)};
      CHECK(hit == !expected.empty());
      CHECK(candidates == expected);
    }
    for (Point const &q : query_points) {
      IndexSet expected;
      for (Integer i{0}; i < n1; ++i) {
        if (b1[i].intersects(q)) {
          expected.insert(i);
        }
      }
      IndexSet candidates;
      tree_1.intersect(q, candidates);
      CHECK(candidates == expected);
    }
  }

  SECTION("Tree-tree intersection") {
    IndexMap expected;
    for (Integer i{0}; i < n1; ++i) {
      for (Integer j{0}; j < n2; ++j) {
        if (b1[i].intersects(b2[j])) {
          expected[i].insert(j);
        }
      }
    }
    IndexMap candidates;
    bool const hit{tree_1.intersect(tree_2, candidates)};
    CHECK(hit == !expected.empty());
    CHECK(candidates == expected);
  }

  SECTION("Self intersection") {
    IndexMap expected_map;
    IndexSet expected_set;
    for (Integer i{0}; i < n1; ++i) {
      for (Integer j{0}; j < n1; ++j) {
        if (i != j && b1[i].intersects(b1[j])) {
          expected_map[i].insert(j);
          expected_set.insert(i);
        }
      }
    }
    IndexMap candidates_map;
    CHECK(tree_1.self_intersect(candidates_map) == !expected_map.empty());
    CHECK(candidates_map == expected_map);

    IndexSet candidates_set;
    CHECK(tree_1.self_intersect(candidates_set) == !expected_set.empty());
    CHECK(candidates_set == expected_set);
  }

  SECTION("Object distance") {
    for (Point const &q : query_points) {
      Real expected{std::numeric_limits<Real>::infinity()};
      for (Integer i{0}; i < n1; ++i) {
        expected = std::min(expected, b1[i].interior_distance(q));
      }
      IndexSet expected_set;
      for (Integer i{0}; i < n1; ++i) {
        if (b1[i].interior_distance(q) == expected) {
          expected_set.insert(i);
        }
      }
      IndexSet candidates;
      CHECK(tree_1.distance(q, candidates) == expected);
      CHECK(candidates == expected_set);
    }
  }

  SECTION("Tree-tree distance") {
    // Shift the second set so that the trees are disjoint
    Real expected{std::numeric_limits<Real>::infinity()};
    for (Integer i{0}; i < n1; ++i) {
      for (Integer j{0}; j < n2; ++j) {
        expected = std::min(expected, b1[i].interior_distance(b2[j]));
      }
    }
    IndexMap expected_map;
    for (Integer i{0}; i < n1; ++i) {
      for (Integer j{0}; j < n2; ++j) {
        if (b1[i].interior_distance(b2[j]) == expected) {
          expected_map[i].insert(j);
        }
      }
    }
    IndexMap candidates;
    CHECK(tree_1.distance(tree_2, candidates) == expected);
    CHECK(candidates == expected_map);
  }

  SECTION("Closest n objects") {
    for (Integer const k : {1, 5, 37}) {
      for (Point const &q : query_points) {
        std::vector<Real> dists(n1);
        for (Integer i{0}; i < n1; ++i) {
          dists[i] = b1[i].interior_distance(q);
        }
        std::vector<Real> sorted{dists};
        std::sort(sorted.begin(), sorted.end());
        Real const kth{sorted[k - 1]};
        // -ffast-math may round the same distance differently at different
        // call sites, so allow a few ulps
        Real const tol{static_cast<Real>(16.0) *
                       std::numeric_limits<Real>::epsilon() *
                       std::max(static_cast<Real>(1.0), kth)};

        IndexSet candidates;
        Real const max_dist{tree_1.closest(q, k, candidates)};
        CAPTURE(k);
        REQUIRE(static_cast<Integer>(candidates.size()) == k);
        CHECK(std::abs(max_dist - kth) <= tol);
        for (Integer c : candidates) {
          CHECK(dists[c] <= kth + tol);
        }
      }
    }
  }

  SECTION("Within distance") {
    Real const max_distance{static_cast<Real>(1.5)};
    for (Point const &q : query_points) {
      IndexSet expected;
      for (Integer i{0}; i < n1; ++i) {
        if (b1[i].interior_distance(q) <= max_distance) {
          expected.insert(i);
        }
      }
      IndexSet candidates;
      CHECK(tree_1.within_distance(q, max_distance, candidates) ==
            !expected.empty());
      CHECK(candidates == expected);
    }
  }

  SECTION("Statistics") {
    typename Tree::Statistics stats;
    tree_1.stats(stats);
    Integer leafs{0}, long_boxes{0}, objects{0};
    for (Integer i{0}; i < tree_1.size(); ++i) {
      auto const &node{tree_1.node(i)};
      objects += node.box_num;
      if (node.child_l < 0) {
        ++leafs;
      } else {
        long_boxes += node.box_num;
      }
    }
    CHECK(objects == n1);
    CHECK(stats.objects == n1);
    CHECK(stats.nodes == tree_1.size());
    CHECK(stats.leafs == leafs);
    CHECK(stats.long_boxes == long_boxes);
  }
}

TEMPLATE_TEST_CASE("Tree handles degenerate inputs", "[template]", float,
                   double) {
  using Real = TestType;
  using Tree = AABBtree::Tree<Real, N_DIM>;
  using Box = AABBtree::Box<Real, N_DIM>;
  using Point = AABBtree::Point<Real, N_DIM>;

  SECTION("Empty box list") {
    Tree tree;
    tree.build(std::make_unique<BoxUniquePtrList<Real, N_DIM>>());
    IndexSet candidates;
    Point const p{Point::Zero()};
    CHECK_FALSE(tree.intersect(p, candidates));
    CHECK(candidates.empty());
    CHECK(tree.distance(p, candidates) < 0);
    CHECK(tree.closest(p, 3, candidates) < 0);
    typename Tree::Statistics stats;
    CHECK_NOTHROW(tree.stats(stats));
    CHECK(stats.objects == 0);
  }

  SECTION("Single box") {
    auto boxes{std::make_unique<BoxUniquePtrList<Real, N_DIM>>()};
    boxes->push_back(
        std::make_unique<Box>(Point::Zero().eval(), Point::Ones().eval()));
    Tree tree;
    tree.build(std::move(boxes));
    IndexSet candidates;
    CHECK(tree.intersect(Point::Constant(0.5).eval(), candidates));
    CHECK(candidates == IndexSet{0});
    IndexMap self;
    CHECK_FALSE(tree.self_intersect(self));
    IndexSet closest;
    CHECK(tree.closest(Point::Constant(3.0).eval(), 4, closest) > 0);
    CHECK(closest == IndexSet{0});
  }
}

TEMPLATE_TEST_CASE("Batch queries match single queries", "[template]", float,
                   double) {
  using Real = TestType;
  using Tree = AABBtree::Tree<Real, N_DIM>;
  using Box = AABBtree::Box<Real, N_DIM>;
  using Point = AABBtree::Point<Real, N_DIM>;

  Integer const num_threads{GENERATE(1, 0, 3)};
  CAPTURE(num_threads);

  Tree tree;
  tree.num_threads(num_threads);
  tree.build(BuildBoxes<Real>(N_OBJECTS_1, 1234u));
  auto const query_ptrs{BuildBoxes<Real>(500, 777u)};
  std::vector<Box> query_boxes;
  std::vector<Point> query_points;
  for (auto const &b : *query_ptrs) {
    query_boxes.push_back(*b);
    query_points.push_back(b->baricenter());
  }

  SECTION("intersect_batch") {
    std::vector<IndexSet> batch;
    bool const any{tree.intersect_batch(query_boxes, batch)};
    REQUIRE(batch.size() == query_boxes.size());
    bool any_single{false};
    for (std::size_t i{0}; i < query_boxes.size(); ++i) {
      IndexSet single;
      any_single |= tree.intersect(query_boxes[i], single);
      CHECK(batch[i] == single);
    }
    CHECK(any == any_single);
    // Containers of unique pointers are accepted as well
    std::vector<IndexSet> batch_ptr;
    tree.intersect_batch(*query_ptrs, batch_ptr);
    CHECK(batch_ptr == batch);
  }

  SECTION("distance_batch") {
    std::vector<Real> distances;
    std::vector<IndexSet> batch;
    tree.distance_batch(query_points, distances, batch);
    for (std::size_t i{0}; i < query_points.size(); ++i) {
      IndexSet single;
      CHECK(distances[i] == tree.distance(query_points[i], single));
      CHECK(batch[i] == single);
    }
  }

  SECTION("closest_batch") {
    auto fun = [](Point const &p, Box const &b) {
      return b.squared_interior_distance(p);
    };
    std::vector<Real> distances;
    std::vector<IndexSet> batch;
    tree.closest_batch(query_points, 9, distances, batch, fun);
    for (std::size_t i{0}; i < query_points.size(); ++i) {
      IndexSet single;
      CHECK(distances[i] == tree.closest(query_points[i], 9, single, fun));
      CHECK(batch[i] == single);
    }
  }

  SECTION("within_distance_batch") {
    std::vector<IndexSet> batch;
    tree.within_distance_batch(query_points, static_cast<Real>(1.0), batch);
    for (std::size_t i{0}; i < query_points.size(); ++i) {
      IndexSet single;
      tree.within_distance(query_points[i], static_cast<Real>(1.0), single);
      CHECK(batch[i] == single);
    }
  }

  SECTION("Exceptions are propagated") {
    auto fun = [](Point const &, Box const &) -> Real {
      throw std::runtime_error("boom");
    };
    std::vector<IndexSet> batch;
    CHECK_THROWS_AS(tree.within_distance_batch(query_points,
                                               static_cast<Real>(1.0), batch,
                                               fun),
                    std::runtime_error);
  }
}

TEST_CASE("Concurrent const queries are safe", "[threads]") {
  using Real = double;
  using Tree = AABBtree::Tree<Real, N_DIM>;
  using Point = AABBtree::Point<Real, N_DIM>;

  Tree tree;
  tree.build(BuildBoxes<Real>(N_OBJECTS_1, 1234u));
  Tree other;
  other.build(BuildBoxes<Real>(300, 99u));

  auto const query_ptrs{BuildBoxes<Real>(200, 777u)};
  std::vector<Point> points;
  for (auto const &b : *query_ptrs) {
    points.push_back(b->baricenter());
  }

  // Serial reference
  std::vector<IndexSet> ref_int(points.size()), ref_cls(points.size());
  for (std::size_t i{0}; i < points.size(); ++i) {
    tree.intersect(*(*query_ptrs)[i], ref_int[i]);
    tree.closest(points[i], 5, ref_cls[i]);
  }
  IndexMap ref_tt;
  tree.intersect(other, ref_tt);

  constexpr int n_threads{4};
  std::vector<int> errors(n_threads, 0);
  std::vector<std::thread> threads;
  for (int t{0}; t < n_threads; ++t) {
    threads.emplace_back([&, t]() {
      for (int rep{0}; rep < 3; ++rep) {
        for (std::size_t i{0}; i < points.size(); ++i) {
          IndexSet c;
          tree.intersect(*(*query_ptrs)[i], c);
          errors[t] += (c != ref_int[i]);
          tree.closest(points[i], 5, c);
          errors[t] += (c != ref_cls[i]);
        }
        IndexMap tt;
        tree.intersect(other, tt);
        errors[t] += (tt != ref_tt);
        typename Tree::Statistics stats;
        tree.stats(stats);
      }
    });
  }
  for (auto &th : threads) {
    th.join();
  }
  for (int t{0}; t < n_threads; ++t) {
    CHECK(errors[t] == 0);
  }
}
