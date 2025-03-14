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
#include <matplot/matplot.h>
using namespace matplot;
static auto fig{figure(true)};
static axes_handle ax{fig->current_axes()};
#endif

// Test utilities
#include "TestUtilities.hh"
using namespace TestUtilities;

#ifdef AABBTREE_ENABLE_PLOTTING
#ifndef SET_PLOT
#define SET_PLOT \
xlim(ax, {-3.0, 3.0}); xlabel(ax, "x"); \
ylim(ax, {-3.0, 3.0}); ylabel(ax, "y"); \
grid(ax, true);
#endif
#endif

TEMPLATE_TEST_CASE("Ray", "[template]", float, double) {

  using Point = Eigen::Matrix<TestType, 2, 1>;
  using Box = AABBtree::Box<TestType, 2>;
  using Ray = AABBtree::Ray<TestType, 2>;

  std::vector<std::string> colors = {"r", "g", "b", "c", "m", "y", "k", "w"};
  TestType const tol{1.0e-8};

  SECTION("Intersection") {
    Ray ray(-1.0, 2.0, 0.0, -1.0);
    Box box(-2.0, -2.0, 1.0, 1.0);
    Point c, f; ray.intersect(box, c, f, tol);
    TestType d_int{ray.interior_distance(box, tol)};
    #ifdef AABBTREE_ENABLE_PLOTTING
    SET_PLOT
    title(ax, "Intersection");
    plot_ray<TestType, 2>(ray, colors[0], 1.0); ax->hold(true);
    plot_box<TestType, 2>(box, colors[1], 1.0);
    plot_point<TestType, 2>(c, colors[2], 4.0);
    plot_point<TestType, 2>(f, colors[4], 2.0); ax->hold(false); show(fig);
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
    SET_PLOT
    title(ax, "Point distance");
    plot_ray<TestType, 2>(ray, colors[0], 1.0); ax->hold(true);
    plot_point<TestType, 2>(pnt, colors[1], 2.0);
    plot_point<TestType, 2>(c, colors[2], 2.0);
    plot_segment<TestType, 2>(pnt, c, colors[4], 2.0); ax->hold(false); show(fig);
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
    SET_PLOT
    title(ax, "Interior distance");
    plot_ray<TestType, 2>(ray, colors[0], 1.0); ax->hold(true);
    plot_box<TestType, 2>(box, colors[1], 1.0);
    plot_point<TestType, 2>(c, colors[2], 4.0);
    plot_point<TestType, 2>(f, colors[4], 2.0); ax->hold(false); show(fig);
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
    SET_PLOT
    title(ax, "Exterior distance");
    plot_ray<TestType, 2>(ray, colors[0], 1.0); ax->hold(true);
    plot_box<TestType, 2>(box, colors[1], 1.0);
    plot_point<TestType, 2>(c, colors[2], 4.0);
    plot_point<TestType, 2>(f, colors[4], 2.0); ax->hold(false); show(fig);
    #endif
    REQUIRE(!ray.intersects(box, tol));
    REQUIRE_THAT(d_int, WithinAbs(ray.exterior_distance(box, tol), tol));
  }
}
