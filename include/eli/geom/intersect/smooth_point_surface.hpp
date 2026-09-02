/*********************************************************************************
* Copyright (c) 2013 David D. Marshall <ddmarsha@calpoly.edu>
*
* All rights reserved. This program and the accompanying materials
* are made available under the terms of the Eclipse Public License v1.0
* which accompanies this distribution, and is available at
* http://www.eclipse.org/legal/epl-v10.html
*
* Contributors:
*    Rob McDonald
********************************************************************************/

#ifndef eli_geom_intersect_smooth_point_surface_hpp
#define eli_geom_intersect_smooth_point_surface_hpp

#include <cmath>
#include <limits>
#include <vector>

#include "eli/code_eli.hpp"
#include "eli/geom/point/distance.hpp"

namespace eli
{
  namespace geom
  {
    namespace intersect
    {
      // Where to move a mesh node so the triangles around it are well shaped.
      //
      // The obvious objective is no use.  Least sum of squared distances to the neighbours has
      // gradient 2n(S-c).Su and 2n(S-c).Sv, with c the centroid, so its minimum is exactly the
      // point of the surface nearest the centroid -- Laplacian smoothing and the projection of
      // the average are the same thing, and being the same thing they share the same fault:
      // the centroid of a one sided ring of neighbours lies off to one side, so a node against
      // a boundary is dragged inward, away from it.
      //
      // Aim at the edge lengths instead.  What the mesher wants of a node is that the edges
      // leaving it are the length it asked for, all of them:
      //
      //     E(u,v) = sum over neighbours of ( |S(u,v) - q_i| - L )^2
      //
      // A ring that is one sided pulls no harder than a ring that surrounds, since each
      // neighbour only asks to be its own distance away.  Nothing in it refers to a centroid.
      //
      // Gauss-Newton, from where the node already is.  Two or three steps, because this is a
      // smoothing pass and not a solve -- the mesh is going to move again anyway.
      template<typename surface__>
      void smooth_point(typename surface__::data_type &u, typename surface__::data_type &v,
                        const surface__ &s,
                        const std::vector<typename surface__::point_type> &q,
                        const typename surface__::data_type &len,
                        const typename surface__::data_type &u0,
                        const typename surface__::data_type &v0,
                        int nstep = 3)
      {
        typedef typename surface__::data_type data_type;
        typedef typename surface__::point_type point_type;

        data_type umin, umax, vmin, vmax;
        s.get_parameter_min(umin,vmin);
        s.get_parameter_max(umax,vmax);

        u = std::min(std::max(u0, umin), umax);
        v = std::min(std::max(v0, vmin), vmax);

        if ( q.empty() )
        {
          return;
        }

        for ( int it = 0; it < nstep; ++it )
        {
          point_type S( s.f(u, v) );
          point_type Su( s.f_u(u, v) );
          point_type Sv( s.f_v(u, v) );

          // Normal equations of the linearized residuals r_i = |S - q_i| - len.
          data_type Auu(0), Auv(0), Avv(0), bu(0), bv(0);

          for ( std::size_t i = 0; i < q.size(); ++i )
          {
            point_type d( S - q[i] );
            data_type dist( d.norm() );

            if ( !(dist > 0) )
            {
              continue;
            }

            data_type r( dist - len );
            data_type ju( d.dot(Su) / dist );      // d(dist)/du
            data_type jv( d.dot(Sv) / dist );

            Auu += ju * ju;
            Auv += ju * jv;
            Avv += jv * jv;
            bu  += ju * r;
            bv  += jv * r;
          }

          data_type det( Auu * Avv - Auv * Auv );

          // A ring that is nearly collinear in parameter space, or a node with too few
          // neighbours, leaves the normal equations singular or nearly so.  Solving them
          // anyway throws the node somewhere absurd, and the surface is then asked for a
          // parameter it does not have.
          data_type scale( Auu * Avv );

          if ( !( det > std::numeric_limits<data_type>::epsilon() * scale ) )
          {
            return;
          }

          data_type du( -( Avv * bu - Auv * bv ) / det );
          data_type dv( -( Auu * bv - Auv * bu ) / det );

          if ( !( std::abs( du ) < std::numeric_limits<data_type>::max() ) ||
               !( std::abs( dv ) < std::numeric_limits<data_type>::max() ) )
          {
            return;                                  // catches NaN as well as overflow
          }

          // A smoothing step has no business moving a node a long way.  Held to a quarter of
          // the patch, which is far more than a well shaped ring ever asks for and stops a
          // badly conditioned one from launching the node across the surface.
          data_type ulim( ( umax - umin ) / static_cast<data_type>(4) );
          data_type vlim( ( vmax - vmin ) / static_cast<data_type>(4) );

          du = std::min( std::max( du, -ulim ), ulim );
          dv = std::min( std::max( dv, -vlim ), vlim );

          u = std::min(std::max(u + du, umin), umax);
          v = std::min(std::max(v + dv, vmin), vmax);
        }
      }
    }
  }
}
#endif
