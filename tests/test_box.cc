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
#include "test_utilities.hh"

#ifdef AABBTREE_ENABLE_PLOTTING
#ifndef SET_PLOT
#define SET_PLOT \
xlim(ax, {-3.0, 3.0}); xlabel(ax, "x"); \
ylim(ax, {-3.0, 3.0}); ylabel(ax, "y"); \
grid(ax, true);
#endif
#endif

TEMPLATE_TEST_CASE("Box", "[template]", float, double) {

  std::vector<std::string> colors = {"r", "g", "b", "c", "m", "y", "k", "w"};

  SECTION("Intersection") {
    using Vector = Eigen::Matrix<TestType, 2, 1>;
    using Box = AABBtree::Box<TestType, 2>;
    Box box_1(-1.0, -1.0, 2.0, 2.0);
    Box box_2(-2.0, -2.0, 1.0, 1.0);
    Box box_3 = box_1.intersection(box_2);
    SET_PLOT
    title(ax, "Intersection");
    plot_box<TestType, 2>(box_1, colors[0], 1.0); ax->hold(true);
    plot_box<TestType, 2>(box_2, colors[1], 1.0);
    plot_box<TestType, 2>(box_3, colors[2], 2.0); ax->hold(false); show(fig);
    REQUIRE(box_3.min().isApprox(Vector(-1.0, -1.0)));
    REQUIRE(box_3.max().isApprox(Vector(1.0, 1.0)));
    REQUIRE(box_1.intersects(box_2));
    REQUIRE(box_2.intersects(box_1));
    REQUIRE_FALSE(box_1.contains(box_2));
    REQUIRE_FALSE(box_2.contains(box_1));
  }

  SECTION("Interior distance") {
    using Vector = Eigen::Matrix<TestType, 2, 1>;
    using Box = AABBtree::Box<TestType, 2>;
    Box box_1(-2.0, -2.0, -0.5, -0.5);
    Box box_2(-2.5, 0.5, 2.0, 2.0);
    Vector p_1, p_2;
    TestType d{box_1.interior_distance(box_2, p_1, p_2)};
    SET_PLOT
    title(ax, "Interior distance");
    plot_box<TestType, 2>(box_1, colors[0], 1.0); ax->hold(true);
    plot_box<TestType, 2>(box_2, colors[1], 1.0);
    plot_segment<TestType, 2>(p_1, p_2, colors[2], 2.0); ax->hold(false); show(fig);
    REQUIRE_THAT(d, WithinAbs(box_1.interior_distance(box_2), 1.0e-8));
  }

  SECTION("Exterior distance") {
    using Vector = Eigen::Matrix<TestType, 2, 1>;
    using Box = AABBtree::Box<TestType, 2>;
    Box box_1(-2.0, -2.0, 0.0, 0.0);
    Box box_2(0.0, 0.0, 2.0, 2.0);
    Vector p_1, p_2;
    TestType d{box_1.exterior_distance(box_2, p_1, p_2)};
    SET_PLOT
    title(ax, "Exterior distance");
    plot_box<TestType, 2>(box_1, colors[0], 1.0); ax->hold(true);
    plot_box<TestType, 2>(box_2, colors[1], 1.0);
    plot_segment<TestType, 2>(p_1, p_2, colors[2], 2.0); ax->hold(false); show(fig);
    REQUIRE_THAT(d, WithinAbs(box_1.exterior_distance(box_2), 1.0e-8));
    REQUIRE(p_1.isApprox(Vector(-2.0, -2.0)));
    REQUIRE(p_2.isApprox(Vector(2.0, 2.0)));
    REQUIRE(box_1.intersects(box_2));
  }
}
