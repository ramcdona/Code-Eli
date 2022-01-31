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
#include "eli/geom/intersect/one_d_curve_solver.hpp"

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
      template<typename curve__>
      typename curve__::data_type intersect_plane(typename curve__::data_type &t, const curve__ &c, const typename curve__::point_type &pt, const typename curve__::point_type &nvec, typename curve__::data_type &t0 )
      {
        typedef typename curve__::onedbezcurve objcurve;
        typedef typename curve__::data_type data_type;

        data_type val;
        objcurve obj = c.signedcurveplanedistcurve(pt, nvec);

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

        objcurve obj = c.signedcurveplanedistcurve(pt, nvec);

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
