/*********************************************************************************
* Copyright (c) 2022 Rob McDonald <ramcdona@calpoly.edu>
*
* All rights reserved. This program and the accompanying materials
* are made available under the terms of the Eclipse Public License v1.0
* which accompanies this distribution, and is available at
* http://www.eclipse.org/legal/epl-v10.html
*
* Contributors:
*    Rob McDonald - initial code and implementation
********************************************************************************/

#ifndef eli_geom_intersect_one_d_curve_solver_hpp
#define eli_geom_intersect_one_d_curve_solver_hpp

#include <cmath>
#include <vector>
#include <list>
#include <algorithm>

#include "eli/code_eli.hpp"

#include "eli/mutil/nls/newton_raphson_method.hpp"

#include "eli/geom/curve/piecewise.hpp"

namespace eli
{
  namespace geom
  {
    namespace curve
    {
      template<template<typename, unsigned short, typename> class curve__, typename data__, unsigned short dim__, typename tol__ >
      class piecewise;
    }

    namespace intersect
    {
      namespace internal
      {
        template <typename onedbezcurve__>
        struct onedbezcurve_g_functor
        {
          const onedbezcurve__ *pc;

          typename onedbezcurve__::data_type operator()(const typename onedbezcurve__::data_type &t) const
          {
            typename onedbezcurve__::data_type tt(t);

            if ( !(tt>=pc->get_t0()) )
            {
              std::cout << "One D Bezier curve g_functor, tt less than minimum.  tt: " << tt << " t0: " << pc->get_t0() << std::endl;
              tt=pc->get_t0();
            }
            if ( !(tt<=pc->get_tmax()) )
            {
              std::cout << "One D Bezier curve g_functor, tt greater than maximum.  tt: " << tt << " tmax: " << pc->get_tmax() << std::endl;
              tt=pc->get_tmax();
            }

            assert((tt>=pc->get_t0()) && (tt<=pc->get_tmax()));

            return pc->f(tt)(0);
          }
        };

        template <typename onedbezcurve__>
        struct onedbezcurve_gp_functor
        {
          const onedbezcurve__ *pc;

          typename onedbezcurve__::data_type operator()(const typename onedbezcurve__::data_type &t) const
          {
            typename onedbezcurve__::data_type tt(t);

            if ( !(tt>=pc->get_t0()) )
            {
              std::cout << "One D Bezier curve gp_functor, tt less than minimum.  tt: " << tt << " t0: " << pc->get_t0() << std::endl;
              tt=pc->get_t0();
            }
            if ( !(tt<=pc->get_tmax()) )
            {
              std::cout << "One D Bezier curve gp_functor, tt greater than maximum.  tt: " << tt << " tmax: " << pc->get_tmax() << std::endl;
              tt=pc->get_tmax();
            }

            assert((tt>=pc->get_t0()) && (tt<=pc->get_tmax()));

            return pc->fp(tt)(0);
          }
        };
      }

      template<typename onedbezcurve__>
      typename onedbezcurve__::data_type find_zero( typename onedbezcurve__::data_type &t, const onedbezcurve__ &c, const typename onedbezcurve__::data_type &t0, const typename onedbezcurve__::data_type &tmin, const typename onedbezcurve__::data_type &tmax )
      {
        eli::mutil::nls::newton_raphson_method<typename onedbezcurve__::data_type> nrm;
        internal::onedbezcurve_g_functor<onedbezcurve__> g;
        internal::onedbezcurve_gp_functor<onedbezcurve__> gp;
        typename onedbezcurve__::data_type val0, val;
        typename onedbezcurve__::tolerance_type tol;

        // setup the functors
        g.pc=&c;
        gp.pc=&c;

        // setup the solver
        nrm.set_absolute_f_tolerance(tol.get_absolute_tolerance());
        nrm.set_max_iteration(10);

        nrm.set_lower_condition( tmin, eli::mutil::nls::iterative_root_base_constrained<typename onedbezcurve__::data_type>::IRC_EXCLUSIVE);
        nrm.set_upper_condition( tmax, eli::mutil::nls::iterative_root_base_constrained<typename onedbezcurve__::data_type>::IRC_EXCLUSIVE);

        // set the initial guess
        nrm.set_initial_guess(t0);
        val0 = c.f(t0)(0);

        // find the root
        nrm.find_root(t, g, gp, 0);

        // if root is within bounds and is closer than initial guess
        {
          assert((t>=c.get_t0()) && (t<=c.get_tmax()));

          val = c.f(t)(0);
          if ( std::abs(val) <= std::abs(val0) )
          {
            return val;
          }
        }

        // couldn't find better answer so return initial guess
        t=t0;
        return val0;
      }

      template<typename onedbezcurve__>
      typename onedbezcurve__::data_type find_zero( typename onedbezcurve__::data_type &t, const onedbezcurve__ &c, const typename onedbezcurve__::data_type &t0 )
      {
        return find_zero( t, c, t0, c.get_t0(), c.get_tmax() );
      }

      template<typename onedbezcurve__>
      void findzeros( std::vector< typename onedbezcurve__::data_type > & optpts,
              const typename onedbezcurve__::data_type &tstart,
              const typename onedbezcurve__::data_type &tend,
              const onedbezcurve__ &objcurve, const typename onedbezcurve__::index_type &nsplit )
      {
        typedef typename onedbezcurve__::data_type data_type;

        data_type tmid = ( tstart + tend ) * 0.5;

        onedbezcurve__ low, high;

        objcurve.split( low, high, 0.5 );

        int nlow = low.numzerocrossings();

        if ( nlow == 0 ) // Curve does not cross zero
        {
        }
        else if ( nlow == -1 ) // Entire curve is zero
        {
          optpts.push_back( ( tstart + tmid ) * 0.5 );
        }
        else  // Curve may have zero crossing
        {
          if ( nsplit <= 0 ) // Reached end of targeted bisection
          { // Check end and mid points, choose closest as candidate start point.
            data_type t = tstart;
            data_type obj = std::abs( low.f( 0 )(0) );
            data_type obj2 = std::abs( low.f( 0.5 )(0) );
            if ( obj2 < obj )
            {
              obj = obj2;
              t = ( tstart + tmid ) * 0.5;
            }
            obj2 = std::abs( low.f( 1 )(0) );
            if ( obj2 < obj )
            {
              t = tmid;
            }
            optpts.push_back( t );
          }
          else // Recurse to next level
          {
            findzeros( optpts, tstart, tmid, low, nsplit - 1 );
          }
        }

        int nhigh = high.numzerocrossings();

        if ( nhigh == 0 ) // Curve does not cross zero
        {
        }
        else if ( nhigh == -1 ) // Entire curve is zero
        {
          optpts.push_back( ( tmid + tend ) * 0.5 );
        }
        else  // Curve may have zero crossing
        {
          if ( nsplit <= 0 ) // Reached end of targeted bisection
          { // Check end and mid points, choose closest as candidate start point.
            data_type t = tmid;
            data_type obj = std::abs( high.f( 0 )(0) );
            data_type obj2 = std::abs( high.f( 0.5 )(0) );
            if ( obj2 < obj )
            {
              obj = obj2;
              t = ( tmid + tend ) * 0.5;
            }
            obj2 = std::abs( high.f( 1 )(0) );
            if ( obj2 < obj )
            {
              t = tend;
            }
            optpts.push_back( t );
          }
          else // Recurse to next level
          {
            findzeros( optpts, tmid, tend, high, nsplit - 1 );
          }
        }
      }

    }
  }
}
#endif
