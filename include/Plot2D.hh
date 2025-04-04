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

#ifndef INCLUDE_PLOT2D_HH
#define INCLUDE_PLOT2D_HH

#include "AABBtree.hh"
#include <matplot/matplot.h>

namespace AABBtree {

  typedef matplot::axes_handle                  plot_obj;
  typedef std::shared_ptr<matplot::figure_type> fig_obj;

  class Plot2D {

    plot_obj m_ax;
    fig_obj  m_fig;

  public:

    Plot2D() {
      m_fig = matplot::figure(true);
      m_ax  = m_fig->current_axes();
    }

    void hold   ( bool const YN ) { m_ax->hold( YN ); }
    void xlabel ( std::string const & s ) { m_ax->xlabel( s ); }
    void ylabel ( std::string const & s ) { m_ax->ylabel( s ); }
    void xlim   ( double const mi, double const ma ) { return m_ax->xlim( {mi,ma} ); }
    void ylim   ( double const mi, double const ma ) { return m_ax->ylim( {mi,ma} ); }
    void title  ( std::string const & s ) { m_ax->title( s ); }  
    void show()  { m_fig->show(); }  
    void clear() { m_ax->clear(); }

    template <typename... Args> auto plot  ( Args... args ) { return m_ax->plot( args... ); }
    template <typename... Args> auto grid  ( Args... args ) { return m_ax->grid( args... ); }
    template <typename... Args> auto loglog( Args... args ) { return m_ax->loglog( args... ); }
    template <typename... Args> auto legend( Args... args ) { return m_ax->legend( args... ); }

    // Plot a box
    template <typename Real, Integer N>
    void
    plot_box(
      Box<Real, N> const & box,
      std::string  const & color,
      Real         const   line_width = 1
    ) {
      static_assert(N == 2, "Plotting is only supported for 2D objects.");
      std::vector<Real> min(box.min().data(), box.min().data() + N);
      std::vector<Real> max(box.max().data(), box.max().data() + N);
      if ( m_ax ) {
        m_ax->plot(
          {min[0], max[0], max[0], min[0], min[0]},
          {min[1], min[1], max[1], max[1], min[1]}
        )->line_width(line_width).color(color);
      }
    }

    // Plot a tree
    template <typename Real, Integer N>
    void
    plot_tree(
      Tree<Real, N> const & tree,
      std::string   const & color,
      Real          const   line_width = 1
    ) {
      static_assert(N == 2, "Plotting is only supported for 2D objects.");
      for (Integer i{0}; i < static_cast<Integer>(tree.size()); ++i) {
        this->plot_box<Real, N>( tree.node(i).box,      color, line_width );
        this->plot_box<Real, N>( tree.node(i).box_long, color, 2.0*line_width );
      }
    }

    // Plot a point
    template <typename Real, Integer N>
    void
    plot_point(
      AABBtree::Vector<Real, N> const & point,
      std::string               const & color,
      Real                      const   marker_size = 0.5
    ) {
      static_assert(N == 2, "Plotting is only supported for 2D objects.");
      std::vector<Real> p(point.data(), point.data() + N);
      if (m_ax) {
        m_ax->plot({p[0]}, {p[1]}, "o")
            ->marker_size(marker_size)
            .color(color)
            .marker_face(true);
      }
    }

    // Plot a segment
    template <typename Real, Integer N>
    void
    plot_segment(
      AABBtree::Vector<Real, N> const & p_1,
      AABBtree::Vector<Real, N> const & p_2,
      std::string               const & color,
      Real                      const   line_width = 1
    ) {
      static_assert(N == 2, "Plotting is only supported for 2D objects.");
      std::vector<Real> v_1(p_1.data(), p_1.data() + N);
      std::vector<Real> v_2(p_2.data(), p_2.data() + N);
      if (m_ax) {
        m_ax->plot({v_1[0], v_2[0]}, {v_1[1], v_2[1]}, "-o")
            ->color(color)
            .line_width(line_width)
            .marker_size(2.0*line_width)
            .marker_face(true);
      }
    }

    // Plot a ray
    template <typename Real, Integer N>
    void
    plot_ray(
      AABBtree::Ray<Real, N> const & ray,
      std::string            const & color,
      Real                   const   line_width = 1
    ) {
      static_assert(N == 2, "Plotting is only supported for 2D objects.");
      Real t{1000.0};
      this->plot_segment<Real, N>( ray.origin(), ray.origin() + t*ray.direction(), color, line_width );
    }

    // Plot a circle
    template <typename Real, Integer N>
    void
    plot_circle(
      Vector<Real, 2> const & center,
      Real            const   radius,
      std::string     const & color,
      Real            const   line_width = 1,
      Integer         const n_points = 100
    ) {
      std::vector<Real> x, y;
      Real const d_theta{2.0*M_PI/n_points};
      for (Integer i{0}; i < n_points; ++i) {
        Real const theta{d_theta*i};
        x.push_back(center[0] + radius*std::cos(theta));
        y.push_back(center[1] + radius*std::sin(theta));
      }
      if (m_ax) {
        m_ax->plot(x, y)->line_width(line_width).color(color);
      }
    }
  };
}

#endif
