/*********************************************************************************
* Copyright (c) 2013 David D. Marshall <ddmarsha@calpoly.edu>
*
* All rights reserved. This program and the accompanying materials
* are made available under the terms of the Eclipse Public License v1.0
* which accompanies this distribution, and is available at
* http://www.eclipse.org/legal/epl-v10.html
*
* Contributors:
*    Rob McDonald - initial code and implementation
********************************************************************************/

#ifndef eli_geom_intersect_minimum_dimension_curve_hpp
#define eli_geom_intersect_minimum_dimension_curve_hpp

#include <cmath>
#include <vector>
#include <list>
#include <algorithm>

#include "eli/code_eli.hpp"

#include "eli/mutil/nls/newton_raphson_method.hpp"

#include "eli/geom/point/distance.hpp"
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
        template <typename curve__>
        struct curve_g_gp_dim_functor
        {
          const curve__ *pc;
          typename curve__::data_type idim;

          void operator()(typename curve__::data_type &g, typename curve__::data_type &gp, const typename curve__::data_type &t) const
          {
            typename curve__::data_type tt(t);

            if ( !(tt>=pc->get_t0()) )
            {
              std::cout << "Minimum dimension curve g_functor, tt less than minimum.  tt: " << tt << " t0: " << pc->get_t0() << std::endl;
              tt=pc->get_t0();
            }
            if ( !(tt<=pc->get_tmax()) )
            {
              std::cout << "Minimum dimension curve g_functor, tt greater than maximum.  tt: " << tt << " tmax: " << pc->get_tmax() << std::endl;
              tt=pc->get_tmax();
            }

            assert((tt>=pc->get_t0()) && (tt<=pc->get_tmax()));

            g = pc->fp(tt)[idim];

            gp = pc->fpp(tt)[idim];
          }
        };
      }

      template<typename curve__>
      typename curve__::data_type minimum_dimension(typename curve__::data_type &t, const curve__ &c, const typename curve__::data_type &idim, const typename curve__::data_type &t0)
      {
        eli::mutil::nls::newton_raphson_method<typename curve__::data_type> nrm;
        internal::curve_g_gp_dim_functor<curve__> ggp;

        typename curve__::data_type r0, r;
        typename curve__::tolerance_type tol;

        // setup the functors
        ggp.pc=&c;
        ggp.idim=idim;

        // setup the solver
        nrm.set_absolute_f_tolerance(tol.get_absolute_tolerance());
        nrm.set_max_iteration(10);
        if (c.open())
        {
          nrm.set_lower_condition(c.get_t0(), eli::mutil::nls::iterative_root_base_constrained<typename curve__::data_type>::IRC_EXCLUSIVE);
          nrm.set_upper_condition(c.get_tmax(), eli::mutil::nls::iterative_root_base_constrained<typename curve__::data_type>::IRC_EXCLUSIVE);
        }
        else
        {
          nrm.set_periodic_condition(c.get_t0(), c.get_tmax());
        }

        // set the initial guess
        nrm.set_initial_guess(t0);
        r0=c.f(t0)[idim];

        // find the root
        nrm.find_root(t, ggp, 0);

        // if root is within bounds and is closer than initial guess
        {
          assert((t>=c.get_t0()) && (t<=c.get_tmax()));

          r = c.f(t)[idim];
          if  (r<=r0)
          {
            return r;
          }
        }

        // couldn't find better answer so return initial guess
        t=t0;
        return r0;
      }

      template<typename curve__>
      typename curve__::data_type minimum_dimension(typename curve__::data_type &t, const curve__ &c, const typename curve__::index_type &idim)
      {
        typedef typename curve__::onedbezcurve objcurve;
        typedef typename curve__::data_type data_type;
        typename std::vector< data_type >::size_type i;
        data_type tt, dd;

        data_type dist = std::numeric_limits<data_type>::max();

        data_type start = 0;
        data_type end = 1;

        // Make dimension derivative squared curve.
        objcurve obj = c.singledimensioncurve( idim );
        objcurve * objderiv = obj.getderiv();
        objcurve osq;
        osq.square( *objderiv );

        // Identify possible zeros from control points.
        std::vector< data_type > optpts;
        findnonpos( optpts, start, end, osq, 6 );

        if ( optpts.empty() )
        {
          optpts.push_back( 0.5 );
        }

        for ( i = 0; i < optpts.size(); i++ )
        {
          data_type t0 = optpts[i];

          dd = minimum_dimension( tt, c, idim, t0 );

          if ( dd < dist )
          {
            dist = dd;
            t = tt;
          }
        }

        std::vector< data_type > forcepts;
        forcepts.push_back( 0 );
        forcepts.push_back( 1 );

        for ( i = 0; i < forcepts.size(); i++ )
        {
          data_type t0 = forcepts[i];

          dd = c.f(t0)( idim );

          if ( dd < dist )
          {
            dist = dd;
            t = t0;
          }
        }

        return dist;
      }

      template<template<typename, unsigned short, typename> class curve__, typename data__, unsigned short dim__, typename tol__>
      typename curve::piecewise<curve__, data__, dim__, tol__>::data_type minimum_dimension(typename curve::piecewise<curve__, data__, dim__, tol__>::data_type &t,
                                                                                           const curve::piecewise<curve__, data__, dim__, tol__> &pc,
                                                                                           const typename curve::piecewise<curve__, data__, dim__, tol__>::index_type &idim)
      {
        typedef curve::piecewise<curve__, data__, dim__, tol__> piecewise_type;
        typedef typename piecewise_type::curve_type curve_type;
        typedef typename piecewise_type::data_type data_type;
        typedef typename piecewise_type::bounding_box_type bounding_box_type;

        typedef typename piecewise_type::segment_collection_type::const_iterator segit;

        typedef std::vector< std::pair<data_type,segit> > dvec;
        dvec minbbdist;

        // Find closest corner of bounding boxes, add them to vector
        // Simple linear search, would be more efficient with some sort of tree.
        for (segit seg=pc.segments.begin(); seg!=pc.segments.end(); ++seg)
        {
          bounding_box_type bb_local;
          seg->second.get_bounding_box(bb_local);

          minbbdist.push_back(std::make_pair(bb_local.get_min()(idim), seg));
        }

        // Sort by nearest distance.
        std::sort( minbbdist.begin(), minbbdist.end(), pairfirstcompare<data_type, segit> );

        // Iterate over segments, starting with nearest bounding box
        data_type dmin(std::numeric_limits<data_type>::max());
        typename dvec::const_iterator it;
        for (it=minbbdist.begin(); it!=minbbdist.end(); ++it)
        {
          // If nearest bb distance is farther than current best, we're done.
          if(it->first < dmin )
          {
            segit seg = it->second;

            curve_type c(seg->second);

            data_type tlocal, d;
            d=minimum_dimension(tlocal,c,idim);

            if(d < dmin)
            {
              data_type tstart(seg->first);
              data_type dt(pc.get_delta_t(seg));

              dmin = d;
              t=tstart+tlocal*dt;
            }
          }
          else
          {
            break;
          }

        }

        return dmin;
      }

    }
  }
}
#endif
