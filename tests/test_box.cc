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

// Matplot++ library
#ifdef AABBTREE_ENABLE_PLOTTING
#include "Plot2D.hh"
#endif

using namespace AABBtree;

// Catch2 library
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
using namespace Catch::Matchers;


// Test utilities
#include "TestUtilities.hh"
using namespace TestUtilities;

TEMPLATE_TEST_CASE("Box", "[template]", float, double) {

  using Vector = Eigen::Matrix<TestType, 2, 1>;
  using Box = AABBtree::Box<TestType, 2>;

  std::vector<std::string> colors = {"r", "g", "b", "c", "m", "y", "k", "w"};
  TestType const tol{1.0e-8};

  SECTION("Intersection") {
    Box box_1(-1.0, -1.0, 2.0, 2.0);
    Box box_2(-2.0, -2.0, 1.0, 1.0);
    Box box_3; box_1.intersect(box_2, box_3);
    #ifdef AABBTREE_ENABLE_PLOTTING
    Plot2D P;
    P.xlim({-3.0, 3.0}); P.xlabel("x");
    P.ylim({-3.0, 3.0}); P.ylabel("y");
    P.grid(true);
    P.title( "Intersection" );
    P.plot_box<TestType, 2>( box_1, colors[0], 1.0);
    P.hold(true);
    P.plot_box<TestType, 2>( box_2, colors[1], 1.0);
    P.plot_box<TestType, 2>( box_3, colors[2], 2.0);
    P.hold(false);
    P.show();
    #endif
    REQUIRE(box_3.min().isApprox(Vector(-1.0, -1.0)));
    REQUIRE(box_3.max().isApprox(Vector(1.0, 1.0)));
    REQUIRE(box_1.intersects(box_2));
    REQUIRE(box_2.intersects(box_1));
    REQUIRE_FALSE(box_1.contains(box_2));
    REQUIRE_FALSE(box_2.contains(box_1));
  }

  SECTION("Point distance") {
    Vector pnt(0.0, 0.0);
    Box box(-2.5, 0.5, 2.0, 2.0);
    Vector c, f;
    TestType d_int{box.interior_distance(pnt, c)};
    TestType d_ext{box.exterior_distance(pnt, f)};
    #ifdef AABBTREE_ENABLE_PLOTTING
    Plot2D P;
    P.xlim({-3.0, 3.0}); P.xlabel("x");
    P.ylim({-3.0, 3.0}); P.ylabel("y");
    P.grid(true);
    P.title( "Point distance" );
    P.plot_box<TestType, 2>( box, colors[0], 1.0); P.hold(true);
    P.plot_segment<TestType, 2>( pnt, c, colors[1], 2.0);
    P.plot_segment<TestType, 2>( pnt, f, colors[2], 2.0);
    P.hold(false);
    P.show();
    #endif
    REQUIRE_THAT(d_int, WithinAbs(box.interior_distance(pnt), tol));
    REQUIRE_THAT(d_ext, WithinAbs(box.exterior_distance(pnt), tol));
  }

  SECTION("Interior distance") {
    Box box_1(-2.0, -2.0, -0.5, -0.5);
    Box box_2(-2.5, 0.5, 2.0, 2.0);
    Vector p_1, p_2;
    TestType d{box_1.interior_distance(box_2, p_1, p_2)};
    #ifdef AABBTREE_ENABLE_PLOTTING
    Plot2D P;
    P.xlim({-3.0, 3.0}); P.xlabel("x");
    P.ylim({-3.0, 3.0}); P.ylabel("y");
    P.grid(true);
    P.title( "Interior distance" );
    P.plot_box<TestType, 2>( box_1, colors[0], 1.0);
    P.hold(true);
    P.plot_box<TestType, 2>( box_2, colors[1], 1.0);
    P.plot_segment<TestType, 2>( p_1, p_2, colors[2], 2.0);
    P.hold(false);
    P.show();
    #endif
    REQUIRE_THAT(d, WithinAbs(box_1.interior_distance(box_2), tol));
  }

  SECTION("Exterior distance") {
    Box box_1(-2.0, -2.0, 0.0, 0.0);
    Box box_2(0.0, 0.0, 2.0, 2.0);
    Vector p_1, p_2;
    TestType d{box_1.exterior_distance(box_2, p_1, p_2)};
    #ifdef AABBTREE_ENABLE_PLOTTING
    Plot2D P;
    P.xlim({-3.0, 3.0}); P.xlabel("x");
    P.ylim({-3.0, 3.0}); P.ylabel("y");
    P.grid(true);
    P.title( "Exterior distance" );
    P.plot_box<TestType, 2>( box_1, colors[0], 1.0);
    P.hold(true);
    P.plot_box<TestType, 2>( box_2, colors[1], 1.0);
    P.plot_segment<TestType, 2>( p_1, p_2, colors[2], 2.0);
    P.hold(false);
    P.show();
    #endif
    REQUIRE_THAT(d, WithinAbs(box_1.exterior_distance(box_2), tol));
    REQUIRE(p_1.isApprox(Vector(-2.0, -2.0)));
    REQUIRE(p_2.isApprox(Vector(2.0, 2.0)));
    REQUIRE(box_1.intersects(box_2));
  }
}
