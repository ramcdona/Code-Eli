/*********************************************************************************
* Copyright (c) 2017 Rob McDonald <ramcdona@calpoly.edu>
*
* All rights reserved. This program and the accompanying materials
* are made available under the terms of the Eclipse Public License v1.0
* which accompanies this distribution, and is available at
* http://www.eclipse.org/legal/epl-v10.html
*
* Contributors:
*    Rob McDonald - initial code and implementation
********************************************************************************/

#ifndef eli_geom_intersect_intersect_plane_curve_hpp
#define eli_geom_intersect_intersect_plane_curve_hpp

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
              std::cout << "Intersect plane curve g_functor, tt less than minimum.  tt: " << tt << " t0: " << pc->get_t0() << std::endl;
              tt=pc->get_t0();
            }
            if ( !(tt<=pc->get_tmax()) )
            {
              std::cout << "Intersect plane curve g_functor, tt greater than maximum.  tt: " << tt << " tmax: " << pc->get_tmax() << std::endl;
              tt=pc->get_tmax();
            }

            assert((tt>=pc->get_t0()) && (tt<=pc->get_tmax()));

            return pc->f(tt);
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
              std::cout << "Intersect plane curve gp_functor, tt less than minimum.  tt: " << tt << " t0: " << pc->get_t0() << std::endl;
              tt=pc->get_t0();
            }
            if ( !(tt<=pc->get_tmax()) )
            {
              std::cout << "Intersect plane curve gp_functor, tt greater than maximum.  tt: " << tt << " tmax: " << pc->get_tmax() << std::endl;
              tt=pc->get_tmax();
            }

            assert((tt>=pc->get_t0()) && (tt<=pc->get_tmax()));

            return pc->fp(tt);
          }
        };
      }

      template<typename onedbezcurve__>
      typename onedbezcurve__::data_type find_zero( typename onedbezcurve__::data_type &t, const onedbezcurve__ &c, const typename onedbezcurve__::data_type &t0 )
      {
        eli::mutil::nls::newton_raphson_method<typename onedbezcurve__::data_type> nrm;
        internal::curve_g_functor<onedbezcurve__> g;
        internal::curve_gp_functor<onedbezcurve__> gp;
        typename onedbezcurve__::data_type val0, val;
        typename onedbezcurve__::tolerance_type tol;

        // setup the functors
        g.pc=&c;
        gp.pc=&c;

        // setup the solver
        nrm.set_absolute_f_tolerance(tol.get_absolute_tolerance());
        nrm.set_max_iteration(10);
        if (c.open())
        {
          nrm.set_lower_condition(c.get_t0(), eli::mutil::nls::iterative_root_base_constrained<typename onedbezcurve__::data_type>::IRC_EXCLUSIVE);
          nrm.set_upper_condition(c.get_tmax(), eli::mutil::nls::iterative_root_base_constrained<typename onedbezcurve__::data_type>::IRC_EXCLUSIVE);
        }
        else
        {
          nrm.set_periodic_condition(c.get_t0(), c.get_tmax());
        }

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
      void findzeros( std::vector< typename onedbezcurve__::data_type > & optpts,
              const typename onedbezcurve__::data_type &tstart,
              const typename onedbezcurve__::data_type &tend,
              const onedbezcurve__ &objcurve, const typename onedbezcurve__::index_type &nsplit )
      {
        typedef typename onedbezcurve__::data_type data_type;

        data_type tmid = ( tstart + tend ) * 0.5;

        data_type smallpos = 10000 * std::numeric_limits< data_type >::epsilon();

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

      template<typename curve__>
      typename curve__::data_type intersect_plane(typename curve__::data_type &t, const curve__ &c, const typename curve__::point_type &pt, const typename curve__::point_type &nvec, typename curve__::data_type &t0 )
      {
        typedef typename curve__::onedbezcurve objcurve;
        typedef typename curve__::data_type data_type;

        data_type val;
        objcurve obj = c.signeddistcurve( pt, nvec );

        val = find_zero( t, obj, t0 );

        return val;
      }

      template<typename curve__>
      typename curve__::data_type intersect_plane(typename curve__::data_type &t, const curve__ &c, const typename curve__::point_type &pt, const typename curve__::point_type &nvec )
      {
        typedef typename curve__::onedbezcurve objcurve;
        typedef typename curve__::data_type data_type;
        typename std::vector< data_type >::size_type i;
        data_type tt, vv;

        data_type val = std::numeric_limits<data_type>::max();

        data_type start = 0;
        data_type end = 1;

        objcurve obj = c.signeddistcurve( pt, nvec );

        int nzero = obj.numzerocrossings();

        if ( nzero <= 0 ) // Entire curve is zero or no zeros.
        {
          t = 0.5;
          val = obj.f( t )(0);
          return val;
        }

        std::vector< data_type > optpts;
        findzeros( optpts, start, end, obj, 6 );

        // No candidate zeros after subdivision.
        if ( optpts.empty() )
        {
          t = 0.5;
          val = obj.f( t )(0);
          return val;
        }

        for ( i = 0; i < optpts.size(); i++ )
        {
          data_type t0 = optpts[i];

          vv = find_zero( tt, obj, t0 );

          if ( std::abs( vv ) < std::abs( val ) )
          {
            val = vv;
            t = tt;
          }
        }

        return val;
      }

      template<template<typename, unsigned short, typename> class curve__, typename data__, unsigned short dim__, typename tol__>
      typename curve::piecewise<curve__, data__, dim__, tol__>::data_type intersect_plane(typename curve::piecewise<curve__, data__, dim__, tol__>::data_type &t,
                                                                                    const curve::piecewise<curve__, data__, dim__, tol__> &pc,
                                                                                    const typename curve::piecewise<curve__, data__, dim__, tol__>::point_type &pt,
                                                                                    const typename curve::piecewise<curve__, data__, dim__, tol__>::point_type &nvec)
      {
        typedef curve::piecewise<curve__, data__, dim__, tol__> piecewise_type;
        typedef typename piecewise_type::data_type data_type;
        typedef typename piecewise_type::bounding_box_type bounding_box_type;

        data_type val = std::numeric_limits<data_type>::max();
        t = pc.get_t0();

        typedef typename piecewise_type::segment_collection_type::const_iterator segit;

        // Check piecewise segment bounding boxes for candidate intersections.
        for (segit seg=pc.segments.begin(); seg!=pc.segments.end(); ++seg)
        {
          bounding_box_type bb_local;
          seg->second.get_bounding_box(bb_local);

          if ( bb_local.intersect_plane( pt, nvec ) ) // Segment is candidate
          {
            data_type tlocal, vv;
            vv = intersect_plane( tlocal, seg->second, pt, nvec );
            if ( std::abs( vv ) < std::abs( val ) )
            {
              data_type tstart( seg->first );
              data_type dt( pc.get_delta_t( seg ) );

              val = vv;
              t = tstart + tlocal * dt;
            }
          }
        }

        return val;
      }
    }
  }
}
#endif
