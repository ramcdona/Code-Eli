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
        struct onedcurve_g_gp_functor
        {
          const onedbezcurve__ *pc;

          void operator()(typename onedbezcurve__::data_type &g, typename onedbezcurve__::data_type &gp, const typename onedbezcurve__::data_type &t) const
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

            g = pc->f(tt)(0);
            gp = pc->fp(tt)(0);
          }
        };
      }

      template<typename onedcurve__>
      typename onedcurve__::data_type find_zero(typename onedcurve__::data_type &t, int &ret, const onedcurve__ &c, const typename onedcurve__::data_type &t0, const typename onedcurve__::data_type &tmin, const typename onedcurve__::data_type &tmax )
      {
        eli::mutil::nls::newton_raphson_method<typename onedcurve__::data_type> nrm;
        internal::onedcurve_g_gp_functor<onedcurve__> ggp;
        typedef typename onedcurve__::data_type data_type;
        data_type val0, val;
        typename onedcurve__::tolerance_type tol;

        // setup the functors
        ggp.pc=&c;

        // setup the solver
        nrm.set_absolute_f_tolerance(tol.get_absolute_tolerance());
        nrm.set_max_iteration(10);

        nrm.set_lower_condition( tmin, eli::mutil::nls::iterative_root_base_constrained<typename onedcurve__::data_type>::IRC_EXCLUSIVE);
        nrm.set_upper_condition( tmax, eli::mutil::nls::iterative_root_base_constrained<typename onedcurve__::data_type>::IRC_EXCLUSIVE);

        // set the initial guess
        nrm.set_initial_guess(t0);
        val0 = c.f(t0)(0);

        // Set up initial guess as worst case scenario.
        t = t0;
        val = val0;

        // find the root
        data_type tnrm = t0;
        ret = nrm.find_root(tnrm, ggp, 0);

        if ( (ret == nrm.converged) && (t>=tmin) && (t<=tmax))
        {
          val0 = c.f(tnrm)(0);
          if ( std::abs(val0) <= std::abs(val) )
          {
            t = tnrm;
            val = val0;
          }
        }

        // Try tmin.
        val0 = c.f(tmin)(0);
        if ( std::abs(val0) <= std::abs(val) )
        {
          t = tmin;
          val = val0;
        }

        // Try tmax.
        val0 = c.f(tmax)(0);
        if ( std::abs(val0) <= std::abs(val) )
        {
          t = tmax;
          val = val0;
        }

        return val;
      }

      template<typename onedcurve__>
      typename onedcurve__::data_type find_zero(typename onedcurve__::data_type &t, const onedcurve__ &c, const typename onedcurve__::data_type &t0, const typename onedcurve__::data_type &tmin, const typename onedcurve__::data_type &tmax )
      {
        int ret;
        return find_zero( t, ret, c, t0, tmin, tmax );
      }

      template<typename onedcurve__>
      typename onedcurve__::data_type find_zero(typename onedcurve__::data_type &t, const onedcurve__ &c, const typename onedcurve__::data_type &t0 )
      {
        return find_zero( t, c, t0, c.get_t0(), c.get_tmax() );
      }

      template<typename onedcurve__>
      typename onedcurve__::data_type find_zero(typename onedcurve__::data_type &t, int &ret, const onedcurve__ &c, const typename onedcurve__::data_type &t0 )
      {
        return find_zero( t, ret, c, t0, c.get_t0(), c.get_tmax(), ret );
      }

      template<typename onedcurve__>
      void findzeros(std::vector< typename onedcurve__::data_type > & optpts,
                     const typename onedcurve__::data_type &tstart,
                     const typename onedcurve__::data_type &tend,
                     const onedcurve__ &objcurve, const typename onedcurve__::index_type &nsplit )
      {
        typedef typename onedcurve__::data_type data_type;

        data_type tmid = ( tstart + tend ) * 0.5;

        onedcurve__ low, high;

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
