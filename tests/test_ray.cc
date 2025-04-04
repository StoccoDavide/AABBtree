/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The AABBtree project is distributed under the BSD 2-Clause License.                           *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// C++17 standard libraries
#include <vector>
#include <iostream>

// AABBtree library
#include "AABBtree.hh"
using namespace AABBtree;

// Catch2 library
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/catch_template_test_macros.hpp>
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
  using Box   = AABBtree::Box<TestType, 2>;
  using Ray   = AABBtree::Ray<TestType, 2>;

  std::vector<std::string> colors = { "r", "g", "b", "c", "m", "y", "k", "w" };
  TestType const tol{1.0e-8};

  SECTION("Intersection") {
    Ray ray(-1.0, 2.0, 0.0, -1.0);
    Box box(-2.0, -2.0, 1.0, 1.0);
    Point c, f; ray.intersect(box, c, f, tol);
    TestType d_int{ray.interior_distance(box, tol)};
    #ifdef AABBTREE_ENABLE_PLOTTING
    Plot2D P;
    P.xlim(-3.0, 3.0); P.xlabel("x");
    P.ylim(-3.0, 3.0); P.ylabel("y");
    P.grid( true);
    P.title("Intersection");
    P.plot_ray<TestType, 2>( ray, colors[0], 1.0);
    P.hold(true);
    P.plot_box<TestType, 2>( box, colors[1], 1.0);
    P.plot_point<TestType, 2>( c, colors[2], 4.0);
    P.plot_point<TestType, 2>( f, colors[4], 2.0);
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
    Point c; ray.distance(pnt, c);
    TestType d{ray.distance(pnt)};
    #ifdef AABBTREE_ENABLE_PLOTTING
    Plot2D P;
    P.xlim(-3.0, 3.0); P.xlabel("x");
    P.ylim(-3.0, 3.0); P.ylabel("y");
    P.grid( true);
    P.title("Point distance");
    P.plot_ray<TestType, 2>( ray, colors[0], 1.0);
    P.hold(true);
    P.plot_point<TestType, 2>( pnt, colors[1], 2.0);
    P.plot_point<TestType, 2>( c, colors[2], 2.0);
    P.plot_segment<TestType, 2>( pnt, c, colors[4], 2.0);
    P.hold(false);
    P.show();
    #endif
    REQUIRE_THAT(4.0, WithinAbs((ray.origin() - c).norm(), tol));
    REQUIRE_THAT(d, WithinAbs(ray.distance(pnt), tol));
  }

  SECTION("Interior distance") {
    Ray ray(2.0, 2.0, 0.0, -1.0);
    Box box(-2.0, -2.0, 1.0, 1.0);
    Point c, f; ray.interior_distance(box, c, f, tol);
    TestType d_int{ray.interior_distance(box, tol)};
    #ifdef AABBTREE_ENABLE_PLOTTING
    Plot2D P;
    P.xlim(-3.0, 3.0); P.xlabel("x");
    P.ylim(-3.0, 3.0); P.ylabel("y");
    P.grid( true);
    P.title("Interior distance");
    P.plot_ray<TestType, 2>( ray, colors[0], 1.0);
    P.hold(true);
    P.plot_box<TestType, 2>( box, colors[1], 1.0);
    P.plot_point<TestType, 2>( c, colors[2], 4.0);
    P.plot_point<TestType, 2>( f, colors[4], 2.0);
    P.hold(false);
    P.show();
    #endif
    REQUIRE(!ray.intersects(box, tol));
    REQUIRE_THAT(d_int, WithinAbs(ray.interior_distance(box, tol), tol));
  }

  SECTION("Exterior distance") {
    Ray ray(2.0, 2.0, 0.5, 1.0);
    Box box(-2.0, -2.0, 1.0, 1.0);
    Point c, f; ray.exterior_distance(box, c, f, tol);
    TestType d_int{ray.exterior_distance(box, tol)};
    #ifdef AABBTREE_ENABLE_PLOTTING
    Plot2D P;
    P.xlim(-3.0, 3.0); P.xlabel("x");
    P.ylim(-3.0, 3.0); P.ylabel("y");
    P.grid( true);
    P.title("Exterior distance");
    P.plot_ray<TestType, 2>( ray, colors[0], 1.0);
    P.hold(true);
    P.plot_box<TestType, 2>( box, colors[1], 1.0);
    P.plot_point<TestType, 2>( c, colors[2], 4.0);
    P.plot_point<TestType, 2>( f, colors[4], 2.0);
    P.hold(false);
    P.show();
    #endif
    REQUIRE(!ray.intersects(box, tol));
    REQUIRE_THAT(d_int, WithinAbs(ray.exterior_distance(box, tol), tol));
  }
}
