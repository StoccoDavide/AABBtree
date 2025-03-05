/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The AABBtree project is distributed under the BSD 2-Clause License.                           *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef INCLUDE_AABBTREE_HH
#define INCLUDE_AABBTREE_HH

// C++17 standard libraries
#include <limits>
#include <type_traits>
#include <array>
#include <vector>
#include <set>
#include <map>
#include <iostream>
#include <utility>
#include <memory>
#include <algorithm>

// Eigen library
#include <Eigen/Dense>

// Print AABBtree errors
#ifndef AABBTREE_ERROR
#define AABBTREE_ERROR(MSG) {         \
  std::ostringstream os; os << MSG;   \
  throw std::runtime_error(os.str()); \
}
#endif

// Assert for AABBtree
#ifndef AABBTREE_ASSERT
#define AABBTREE_ASSERT(COND, MSG) if (!(COND)) {AABBTREE_ERROR(MSG);}
#endif

// Warning for AABBtree
#ifndef AABBTREE_WARNING
#define AABBTREE_WARNING(MSG) {std::cout << MSG << std::endl;}
#endif

// Warning assert for AABBtree
#ifndef AABBTREE_ASSERT_WARNING
#define AABBTREE_ASSERT_WARNING(COND, MSG) if (!(COND)) {AABBTREE_WARNING(MSG);}
#endif

// Default integer AABBtree type
#ifndef AABBTREE_DEFAULT_INTEGER_TYPE
#define AABBTREE_DEFAULT_INTEGER_TYPE int
#endif

/**
* \brief Namespace for the AABBtree library.
*
* The AABBtree namespace contains all the classes and functions of the AABBtree library.
*/
namespace AABBtree {

  using std::sqrt;
  using std::swap;
  using std::numeric_limits;
  using std::is_floating_point;
  using std::is_integral;
  using std::is_same;
  using std::distance;
  using std::max_element;
  using std::min_element;

  using std::vector;
  using std::set;
  using std::map;

  using std::unique_ptr;
  using std::make_unique;

  /**
  * \brief The Integer type used in the AABBtree class.
  *
  * The Integer type, is defined by the preprocessor symbol \c AABBTREE_DEFAULT_INTEGER_TYPE. The
  * default value is \c int.
  */
  using Integer = AABBTREE_DEFAULT_INTEGER_TYPE;
  using Set     = set<Integer>;        /**< Set of indexes. */

  using istream_type = std::basic_istream<char>;
  using ostream_type = std::basic_ostream<char>;

  static_assert( is_integral<Integer>::value, "Integer must be an integer type." );

  template <Integer N, typename Real>                       class Box;
  template <Integer N, typename Real, typename DerivedTree> class Tree;

  template <Integer N, typename Real> using BoxUPtr    = unique_ptr<Box<N,Real>>;
  template <Integer N, typename Real> using BoxUPtrVec = vector<BoxUPtr<N,Real>>;
  template <typename Real, Integer N> using Vector     = Eigen::Vector<Real,N>; /**> Eigen column vector of real numbers. */
  template <typename Real, Integer N> using Point      = Eigen::Vector<Real,N>; /**> Point in the ambient space (Eigen column vector of real numbers). */

} // namespace AABBtree

#include "AABBtree/Box.hxx"
#include "AABBtree/Tree.hxx"
#include "AABBtree/RecursiveTree.hxx"
//#include "AABBtree/NonRecursiveTree.hxx"

#endif // INCLUDE_AABBTREE_HH
