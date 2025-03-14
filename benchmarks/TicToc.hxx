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

#ifndef INCLUDE_AABBTREE_UTILITIES_TICTOC_HXX
#define INCLUDE_AABBTREE_UTILITIES_TICTOC_HXX

#include "AABBtree/Box.hxx"

namespace AABBtree {


    /**
    * \brief Class container for a timer.
    *
    * Class container for a timer that uses the C++ standard library chrono to measure the elapsed time
    * between two points in the code.
    * \tparam Real Scalar number type.
    */
    template <typename Real>
    class TicToc {

      // Basic types definitions
      using Clock = std::chrono::high_resolution_clock; /**< Clock type. */
      using Time  = std::chrono::microseconds; /**< Time type. */

      // Timer variables
      Clock::time_point m_start_time; /**< Start time. */
      Clock::time_point m_stop_time; /**< Stop time. */
      Time m_elapsed_time; /**< Elapsed time. */

    public:
      /**
      * Copy constructor for the timer.
      */
      TicToc(const TicToc &) = delete;

      /**
      * Assignment operator for the timer.
      */
      TicToc & operator=(TicToc const &) = delete;

      /**
      * Class constructor for the timer.
      */
      TicToc() : m_elapsed_time(0) {}

      /**
      * Start the timer.
      */
      void tic() {this->m_start_time = Clock::now();}

      /**
      * Stop the timer.
      */
      void toc()
      {
        this->m_stop_time    = Clock::now();
        this->m_elapsed_time = std::chrono::duration_cast<Time>(this->m_stop_time - this->m_start_time);
      }

      /**
      * Get the elapsed time in seconds.
      * \return The elapsed time in seconds.
      */
      Real elapsed_s() const {return static_cast<Real>(1.0e-6)*this->m_elapsed_time.count();}

      /**
      * Get the elapsed time in milliseconds.
      * \return The elapsed time in milliseconds.
      */
      Real elapsed_ms() const {return static_cast<Real>(1.0e-3)*this->m_elapsed_time.count();}

      /**
      * Get the elapsed time in microseconds.
      * \return The elapsed time in microseconds.
      */
      Real elapsed_us() const {return this->m_elapsed_time.count();}

    }; // class TicToc

  } // namespace Utilities

} // namespace AABBtree

#endif // INCLUDE_AABBTREE_UTILITIES_TICTOC_HXX