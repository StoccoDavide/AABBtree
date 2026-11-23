/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The AABBtree project is distributed under the BSD 2-Clause License.       *
 *                                                                           *
 * Davide Stocco Enrico Bertolazzi                                           *
 * University of Trento University of Trento                                 *
 * davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef INCLUDE_AABBTREE_HH
#define INCLUDE_AABBTREE_HH

// C++20 standard libraries
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <concepts>
#include <cstdlib>
#include <exception>
#include <functional>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <queue>
#include <ranges>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <typeinfo>
#include <utility>
#include <vector>

// Eigen library
#include <Eigen/Dense>
#include <Eigen/Geometry>

// OpenMP (optional, used by the batch queries)
#ifdef _OPENMP
#include <omp.h>
#endif

// Print AABBtree errors
#ifndef AABBTREE_ERROR
#define AABBTREE_ERROR(MSG)                                                    \
  {                                                                            \
    std::ostringstream os;                                                     \
    os << MSG;                                                                 \
    throw std::runtime_error(os.str());                                        \
  }
#endif

// Assert for AABBtree
#ifndef AABBTREE_ASSERT
#define AABBTREE_ASSERT(COND, MSG)                                             \
  if (!(COND)) {                                                               \
    AABBTREE_ERROR(MSG);                                                       \
  }
#endif

// Warning for AABBtree
#ifndef AABBTREE_WARNING
#define AABBTREE_WARNING(MSG)                                                  \
  {                                                                            \
    std::cout << MSG << std::endl;                                             \
  }
#endif

// Warning assert for AABBtree
#ifndef AABBTREE_ASSERT_WARNING
#define AABBTREE_ASSERT_WARNING(COND, MSG)                                     \
  if (!(COND)) {                                                               \
    AABBTREE_WARNING(MSG);                                                     \
  }
#endif

// Minimum number of nodes (sum over the two trees) for a tree-tree query to
// run in parallel with OpenMP
#ifndef AABBTREE_PARALLEL_MIN_NODES
#define AABBTREE_PARALLEL_MIN_NODES 1024
#endif

// Subtrees with more boxes than this are built as separate OpenMP tasks
#ifndef AABBTREE_PARALLEL_BUILD_GRAIN
#define AABBTREE_PARALLEL_BUILD_GRAIN 4096
#endif

// Default integer AABBtree type
#ifndef AABBTREE_DEFAULT_INTEGER_TYPE
#define AABBTREE_DEFAULT_INTEGER_TYPE int
#endif

/**
 * \brief Namespace for the AABBtree library.
 *
 * The AABBtree namespace contains all the classes and functions of the AABBtree
 * library.
 */
namespace AABBtree {

/**
 * \brief The Integer type used in the AABBtree class.
 *
 * The Integer type, is defined by the preprocessor symbol \c
 * AABBTREE_DEFAULT_INTEGER_TYPE. The default value is \c int.
 */
using Integer = AABBTREE_DEFAULT_INTEGER_TYPE;
static_assert(std::signed_integral<Integer>,
              "AABBTREE_DEFAULT_INTEGER_TYPE must be a signed integral type.");

using IndexSet = std::set<Integer>;
using IndexMap = std::map<Integer, IndexSet>;
using IndexList = std::vector<Integer>;
using IndexPair = std::pair<Integer, Integer>;
using IndexPairList = std::vector<IndexPair>;
using OutStream = std::basic_ostream<char>;

// Type aliases
template <std::floating_point Real, Integer N> class Box;
template <typename Real, Integer N>
using BoxUniquePtr = std::unique_ptr<Box<Real, N>>;
template <typename Real, Integer N>
using BoxUniquePtrList = std::vector<BoxUniquePtr<Real, N>>;
template <typename Real, Integer N> using Vector = Eigen::Vector<Real, N>;
template <typename Real, Integer N> using Point = Eigen::Vector<Real, N>;

} // namespace AABBtree

#include "AABBtree/Box.hh"
#include "AABBtree/Ray.hh"
#include "AABBtree/Tree.hh"

#endif // INCLUDE_AABBTREE_HH
