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

// C++17 standard libraries
#include <algorithm>
#include <array>
#include <functional>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <queue>
#include <set>
#include <type_traits>
#include <utility>
#include <vector>

// Eigen library
#include <Eigen/Dense>
#include <Eigen/Geometry>

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
static_assert(std::is_integral<Integer>::value,
              "AABBTREE_DEFAULT_INTEGER_TYPE must be an integral type.");

using IndexSet = std::set<Integer>;
using IndexMap = std::map<Integer, IndexSet>;
using IndexList = std::vector<Integer>;
using OutStream = std::basic_ostream<char>;

// Type aliases
template <typename Real, Integer N> class Box;
template <typename Real, Integer N>
using BoxUniquePtr = std::unique_ptr<Box<Real, N>>;
template <typename Real, Integer N>
using BoxUniquePtrList = std::vector<BoxUniquePtr<Real, N>>;
template <typename Real, Integer N> using Vector = Eigen::Vector<Real, N>;
template <typename Real, Integer N> using Point = Eigen::Vector<Real, N>;

} // namespace AABBtree

#include "AABBtree/Box.hxx"
#include "AABBtree/Ray.hxx"
#include "AABBtree/Tree.hxx"

#endif // INCLUDE_AABBTREE_HH
