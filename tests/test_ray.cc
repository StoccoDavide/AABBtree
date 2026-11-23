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
#include <iostream>
#include <random>
#include <vector>

// AABBtree library
#include "AABBtree.hh"
using namespace AABBtree;

// Catch2 library
#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
using namespace Catch::Matchers;

// Matplot++ library
#ifdef AABBTREE_ENABLE_PLOTTING
#include "Plot2D.hh"
#endif

// Test utilities
#include "TestUtilities.hh"
using namespace TestUtilities;

TEMPLATE_TEST_CASE("Ray", "[template]", float, double) {

  using Point = Eigen::Matrix<TestType, 2, 1>;
  using Box = AABBtree::Box<TestType, 2>;
  using Ray = AABBtree::Ray<TestType, 2>;

  std::vector<std::string> colors = {"r", "g", "b", "c", "m", "y", "k", "w"};
  TestType const tol{1.0e-8};

  SECTION("Intersection") {
    Ray ray(-1.0, 2.0, 0.0, -1.0);
    Box box(-2.0, -2.0, 1.0, 1.0);
    Point c, f;
    ray.intersect(box, c, f, tol);
    TestType d_int{ray.interior_distance(box, tol)};
#ifdef AABBTREE_ENABLE_PLOTTING
    Plot2D P;
    P.xlim(-3.0, 3.0);
    P.xlabel("x");
    P.ylim(-3.0, 3.0);
    P.ylabel("y");
    P.grid(true);
    P.title("Intersection");
    P.plot_ray<TestType, 2>(ray, colors[0], 1.0);
    P.hold(true);
    P.plot_box<TestType, 2>(box, colors[1], 1.0);
    P.plot_point<TestType, 2>(c, colors[2], 4.0);
    P.plot_point<TestType, 2>(f, colors[4], 2.0);
    P.hold(false);
    P.show();
#endif
    REQUIRE(ray.intersects(box, tol));
    REQUIRE_THAT(1.0, WithinAbs((ray.origin() - c).norm(), tol));
    REQUIRE_THAT(4.0, WithinAbs((ray.origin() - f).norm(), tol));
    REQUIRE_THAT(d_int, WithinAbs(ray.interior_distance(box, tol), tol));
  }

  SECTION("Point distance") {
    Ray ray(-1.0, 2.0, 0.0, -1.0);
    Point pnt(2.0, -2.0);
    Point c;
    ray.distance(pnt, c);
    TestType d{ray.distance(pnt)};
#ifdef AABBTREE_ENABLE_PLOTTING
    Plot2D P;
    P.xlim(-3.0, 3.0);
    P.xlabel("x");
    P.ylim(-3.0, 3.0);
    P.ylabel("y");
    P.grid(true);
    P.title("Point distance");
    P.plot_ray<TestType, 2>(ray, colors[0], 1.0);
    P.hold(true);
    P.plot_point<TestType, 2>(pnt, colors[1], 2.0);
    P.plot_point<TestType, 2>(c, colors[2], 2.0);
    P.plot_segment<TestType, 2>(pnt, c, colors[4], 2.0);
    P.hold(false);
    P.show();
#endif
    REQUIRE_THAT(4.0, WithinAbs((ray.origin() - c).norm(), tol));
    REQUIRE_THAT(d, WithinAbs(ray.distance(pnt), tol));
  }

  SECTION("Interior distance") {
    Ray ray(2.0, 2.0, 0.0, -1.0);
    Box box(-2.0, -2.0, 1.0, 1.0);
    Point c, f;
    ray.interior_distance(box, c, f, tol);
    TestType d_int{ray.interior_distance(box, tol)};
#ifdef AABBTREE_ENABLE_PLOTTING
    Plot2D P;
    P.xlim(-3.0, 3.0);
    P.xlabel("x");
    P.ylim(-3.0, 3.0);
    P.ylabel("y");
    P.grid(true);
    P.title("Interior distance");
    P.plot_ray<TestType, 2>(ray, colors[0], 1.0);
    P.hold(true);
    P.plot_box<TestType, 2>(box, colors[1], 1.0);
    P.plot_point<TestType, 2>(c, colors[2], 4.0);
    P.plot_point<TestType, 2>(f, colors[4], 2.0);
    P.hold(false);
    P.show();
#endif
    REQUIRE(!ray.intersects(box, tol));
    REQUIRE_THAT(d_int, WithinAbs(ray.interior_distance(box, tol), tol));
  }

  SECTION("Exterior distance") {
    Ray ray(2.0, 2.0, 0.5, 1.0);
    Box box(-2.0, -2.0, 1.0, 1.0);
    Point c, f;
    ray.exterior_distance(box, c, f, tol);
    TestType d_int{ray.exterior_distance(box, tol)};
#ifdef AABBTREE_ENABLE_PLOTTING
    Plot2D P;
    P.xlim(-3.0, 3.0);
    P.xlabel("x");
    P.ylim(-3.0, 3.0);
    P.ylabel("y");
    P.grid(true);
    P.title("Exterior distance");
    P.plot_ray<TestType, 2>(ray, colors[0], 1.0);
    P.hold(true);
    P.plot_box<TestType, 2>(box, colors[1], 1.0);
    P.plot_point<TestType, 2>(c, colors[2], 4.0);
    P.plot_point<TestType, 2>(f, colors[4], 2.0);
    P.hold(false);
    P.show();
#endif
    REQUIRE(!ray.intersects(box, tol));
    REQUIRE_THAT(d_int, WithinAbs(ray.exterior_distance(box, tol), tol));
  }
}

namespace {

// Brute-force minimum over sampled ray points of a point-box distance
template <typename Ray, typename Box, typename Fun>
double sampled_minimum(Ray const &ray, Box const &box, Fun const &fun) {
  double best{std::numeric_limits<double>::infinity()};
  constexpr int n_samples{100000};
  constexpr double t_max{80.0};
  for (int k{0}; k <= n_samples; ++k) {
    double const t{t_max * k / n_samples};
    auto const p{ray.point(static_cast<decltype(ray.direction()[0] + 0)>(t))};
    best = std::min(best, static_cast<double>(fun(box, p)));
  }
  return best;
}

} // namespace

TEMPLATE_TEST_CASE("Ray-box queries match brute force", "[template]", float,
                   double) {
  using Real = TestType;
  constexpr AABBtree::Integer D{3};
  using Point = AABBtree::Point<Real, D>;
  using Box = AABBtree::Box<Real, D>;
  using Ray = AABBtree::Ray<Real, D>;

  std::mt19937_64 gen(2026);
  std::uniform_real_distribution<Real> coord(-10, 10), size(0.1, 4);
  // Sampling resolution (80 / 1e5 along a unit direction) plus float rounding
  double const tol{std::is_same<Real, float>::value ? 5.0e-3 : 2.0e-3};

  for (int test{0}; test < 60; ++test) {
    CAPTURE(test);
    Point const o(coord(gen), coord(gen), coord(gen));
    Point d(coord(gen), coord(gen), coord(gen));
    if (test % 5 == 0) {
      d[test % 3] = 0; // Parallel to an axis
    }
    Ray const ray(o, d.normalized());
    Point const c(coord(gen), coord(gen), coord(gen));
    Point const h(size(gen), size(gen), size(gen));
    Box const box(Point(c - h), Point(c + h));

    // Interior distance: exact minimum, never above any sample
    Point p1, p2;
    double const d_int{ray.interior_distance(box, p1, p2)};
    double const d_int_brute{sampled_minimum(
        ray, box, [](Box const &b, Point const &p) {
          return b.interior_distance(p);
        })};
    CHECK(d_int <= d_int_brute + 1.0e-4);
    CHECK(d_int >= d_int_brute - tol);
    CHECK(std::abs(d_int - (p1 - p2).norm()) < 1.0e-4);
    CHECK(ray.distance(p1) < 1.0e-4); // p1 is on the ray
    CHECK(box.interior_distance(p2) < 1.0e-4); // p2 is in the box

    // Exterior distance
    double const d_ext{ray.exterior_distance(box)};
    double const d_ext_brute{sampled_minimum(
        ray, box, [](Box const &b, Point const &p) {
          return b.exterior_distance(p);
        })};
    CHECK(d_ext <= d_ext_brute + 1.0e-4);
    CHECK(d_ext >= d_ext_brute - tol);
    CHECK(d_ext >= d_int);

    // Intersection is consistent with the interior distance
    Point entry, exit;
    bool const hit{ray.intersect(box, entry, exit)};
    CHECK(hit == ray.intersects(box));
    if (d_int > 1.0e-3) {
      CHECK_FALSE(hit);
    }
    if (d_int == 0) {
      CHECK(hit);
    }
    if (hit) {
      CHECK(box.interior_distance(entry) < 1.0e-4);
      CHECK(box.interior_distance(exit) < 1.0e-4);
      CHECK((entry - o).dot(d) >= -1.0e-4); // Not behind the origin
    }
  }
}

TEMPLATE_TEST_CASE("Ray utilities", "[template]", float, double) {
  using Real = TestType;
  using Point2 = AABBtree::Point<Real, 2>;
  using Ray2 = AABBtree::Ray<Real, 2>;
  using Box2 = AABBtree::Box<Real, 2>;

  Ray2 const ray(1.0, 1.0, 2.0, 0.0);

  SECTION("contains") {
    CHECK(ray.contains(Point2(1.0, 1.0)));      // Origin
    CHECK(ray.contains(Point2(5.0, 1.0)));      // Ahead
    CHECK_FALSE(ray.contains(Point2(0.0, 1.0))); // Behind
    CHECK_FALSE(ray.contains(Point2(5.0, 2.0))); // Off the ray
  }

  SECTION("point distance") {
    CHECK(ray.distance(Point2(4.0, 3.0)) == Catch::Approx(2.0));
    CHECK(ray.distance(Point2(-2.0, 5.0)) == Catch::Approx(5.0)); // Behind
  }

  SECTION("intersect from inside") {
    Box2 const box(0.0, 0.0, 3.0, 2.0);
    Point2 entry, exit;
    REQUIRE(ray.intersect(box, entry, exit));
    CHECK(entry.isApprox(Point2(1.0, 1.0)));
    CHECK(exit.isApprox(Point2(3.0, 1.0)));
  }

  SECTION("transform") {
    Eigen::Transform<Real, 2, Eigen::Affine> t{Eigen::Transform<Real, 2, Eigen::Affine>::Identity()};
    t.rotate(Real(EIGEN_PI / 2)).pretranslate(Point2(1.0, 0.0));
    Ray2 const r{ray.transformed(t)};
    CHECK(r.origin().isApprox(Point2(0.0, 1.0)));
    CHECK(r.direction().isApprox(Point2(0.0, 2.0)));

    // A box rotated by 90 degrees swaps its sizes around the new center
    Box2 const box(0.0, 0.0, 4.0, 2.0);
    Box2 const b{box.transformed(t)};
    CHECK(b.min()[0] == Catch::Approx(-1.0));
    CHECK(b.max()[0] == Catch::Approx(1.0));
    CHECK(b.min()[1] == Catch::Approx(0.0));
    CHECK(b.max()[1] == Catch::Approx(4.0));
  }
}

TEMPLATE_TEST_CASE("Tree ray queries match brute force", "[template]", float,
                   double) {
  using Real = TestType;
  constexpr AABBtree::Integer D{3};
  using Point = AABBtree::Point<Real, D>;
  using Box = AABBtree::Box<Real, D>;
  using Ray = AABBtree::Ray<Real, D>;
  using Tree = AABBtree::Tree<Real, D>;

  std::mt19937_64 gen(7);
  std::uniform_real_distribution<Real> coord(0, 50), size(0.1, 1.5), dir(-1, 1);
  auto boxes{std::make_unique<AABBtree::BoxUniquePtrList<Real, D>>()};
  std::vector<Box> copy;
  for (int i{0}; i < 3000; ++i) {
    Point const c(coord(gen), coord(gen), coord(gen));
    Point const h(size(gen), size(gen), size(gen));
    copy.emplace_back(Point(c - h), Point(c + h));
    boxes->push_back(std::make_unique<Box>(copy.back()));
  }
  Tree tree;
  tree.max_nodal_objects(GENERATE(1, 16));
  tree.build(std::move(boxes));

  std::vector<Ray> rays;
  for (int k{0}; k < 30; ++k) {
    Point const o(coord(gen), coord(gen), coord(gen));
    rays.emplace_back(o, Point(dir(gen), dir(gen), dir(gen)).normalized());
  }
  for (Ray const &ray : rays) {
    AABBtree::IndexSet expected_hits, hits;
    Real min_dist{std::numeric_limits<Real>::infinity()};
    std::vector<Real> dist(copy.size());
    for (std::size_t i{0}; i < copy.size(); ++i) {
      if (ray.intersects(copy[i])) {
        expected_hits.insert(static_cast<AABBtree::Integer>(i));
      }
      dist[i] = copy[i].interior_distance(ray);
      min_dist = std::min(min_dist, dist[i]);
    }
    tree.intersect(ray, hits);
    CHECK(hits == expected_hits);

    AABBtree::IndexSet nearest;
    CHECK(tree.distance(ray, nearest) == Catch::Approx(min_dist).margin(1e-5));
    for (AABBtree::Integer i : nearest) {
      CHECK(dist[i] == Catch::Approx(min_dist).margin(1e-5));
    }

    std::vector<Real> sorted{dist};
    std::sort(sorted.begin(), sorted.end());
    AABBtree::IndexSet closest;
    Real const kth{tree.closest(ray, 7, closest)};
    CHECK(closest.size() == 7);
    CHECK(kth == Catch::Approx(sorted[6]).margin(1e-5));
  }
}
