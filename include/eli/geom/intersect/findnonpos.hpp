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

#ifndef eli_geom_intersect_findnonpos_hpp
#define eli_geom_intersect_findnonpos_hpp

#include <cmath>
#include <vector>
#include <list>
#include <algorithm>
#include <limits>

#include "eli/code_eli.hpp"


namespace eli
{
  namespace geom
  {
    namespace intersect
    {

      template<typename onedsurf__>
      void findnonpos( std::vector< std::pair< typename onedsurf__::data_type, typename onedsurf__::data_type > > & optpts,
              const std::pair< typename onedsurf__::data_type, typename onedsurf__::data_type > &uvstart,
              const std::pair< typename onedsurf__::data_type, typename onedsurf__::data_type > &uvend,
              const onedsurf__ &objsurf, const typename onedsurf__::index_type &nsplit )
      {
        typedef typename onedsurf__::data_type data_type;
        typedef std::pair< data_type, data_type > uvpair;

        uvpair uvmid = std::make_pair( (uvstart.first + uvend.first) * 0.5, (uvstart.second + uvend.second) * 0.5 );

        data_type smallpos = 10000 * std::numeric_limits< data_type >::epsilon();

        onedsurf__ ulow, uhigh;
        onedsurf__ ulowvlow, ulowvhigh, uhighvlow, uhighvhigh;

        objsurf.split_u( ulow, uhigh, 0.5 );
        ulow.split_v( ulowvlow, ulowvhigh, 0.5 );
        uhigh.split_v( uhighvlow, uhighvhigh, 0.5 );

        if ( !ulowvlow.allpos( smallpos ) )
        {
          if ( nsplit <= 0 )
          {
            uvpair uv = std::make_pair( (uvstart.first + uvmid.first) * 0.5, (uvstart.second + uvmid.second) * 0.5 );
            optpts.push_back( uv );
          }
          else
          {
            findnonpos( optpts, uvstart, uvmid, ulowvlow, nsplit - 1 );
          }
        }

        if ( !ulowvhigh.allpos( smallpos ) )
        {
          if ( nsplit <= 0 )
          {
            uvpair uv = std::make_pair( (uvstart.first + uvmid.first) * 0.5, (uvmid.second + uvend.second) * 0.5 );
            optpts.push_back( uv );
          }
          else
          {
            uvpair uvone = std::make_pair( uvstart.first, uvmid.second );
            uvpair uvtwo = std::make_pair( uvmid.first, uvend.second );

            findnonpos( optpts, uvone, uvtwo, ulowvhigh, nsplit - 1 );
          }
        }

        if ( !uhighvlow.allpos( smallpos ) )
        {
          if ( nsplit <= 0 )
          {
            uvpair uv = std::make_pair( (uvmid.first + uvend.first) * 0.5, (uvstart.second + uvmid.second) * 0.5 );
            optpts.push_back( uv );
          }
          else
          {
            uvpair uvone = std::make_pair( uvmid.first, uvstart.second );
            uvpair uvtwo = std::make_pair( uvend.first, uvmid.second );

            findnonpos( optpts, uvone, uvtwo, uhighvlow, nsplit - 1 );
          }
        }

        if ( !uhighvhigh.allpos( smallpos ) )
        {
          if ( nsplit <= 0 )
          {
            uvpair uv = std::make_pair( (uvmid.first + uvend.first) * 0.5, (uvmid.second + uvend.second) * 0.5 );
            optpts.push_back( uv );
          }
          else
          {
            findnonpos( optpts, uvmid, uvend, uhighvhigh, nsplit - 1 );
          }
        }

      }



      // Divides out of f a factor u, 1 - u, v or 1 - v wherever f's control points along that edge are all within
      // tol of zero, as they are along an edge collapsed to a point.  The quotient has the zeros f has inside the
      // patch but not those along the edge, which would otherwise keep every cell beside it.  f zero throughout is
      // left as it is.
      template<typename onedsurf__>
      void deflate_zero_edges( onedsurf__ &f, const typename onedsurf__::data_type &tol )
      {
        typedef typename onedsurf__::data_type data_type;
        typedef typename onedsurf__::index_type index_type;
        typedef typename onedsurf__::point_type point_type;

        if ( f.allnear( tol ) )
        {
          return;
        }

        bool changed = true;
        while ( changed )
        {
          changed = false;
          index_type n = f.degree_u();
          index_type m = f.degree_v();

          for ( int edge = 0; edge < 4 && !changed; edge++ )
          {
            bool alongu = ( edge < 2 );
            index_type deg = m;
            index_type len = n;
            if ( alongu )
            {
              deg = n;
              len = m;
            }
            if ( deg < 1 )
            {
              continue;
            }
            index_type row = 0;
            if ( edge == 1 || edge == 3 )
            {
              row = deg;
            }

            bool zero = true;
            for ( index_type k = 0; k <= len && zero; k++ )
            {
              point_type p;
              if ( alongu )
              {
                p = f.get_control_point( row, k );
              }
              else
              {
                p = f.get_control_point( k, row );
              }
              zero = std::abs( p(0) ) <= tol;
            }
            if ( !zero )
            {
              continue;
            }

            onedsurf__ g;
            if ( alongu )
            {
              g.resize( n - 1, m );
            }
            else
            {
              g.resize( n, m - 1 );
            }
            for ( index_type i = 0; i < deg; i++ )
            {
              data_type c;
              index_type from;
              if ( row == 0 )
              {
                c = static_cast< data_type >( deg ) / static_cast< data_type >( i + 1 );
                from = i + 1;
              }
              else
              {
                c = static_cast< data_type >( deg ) / static_cast< data_type >( deg - i );
                from = i;
              }
              for ( index_type k = 0; k <= len; k++ )
              {
                if ( alongu )
                {
                  g.set_control_point( c * f.get_control_point( from, k ), i, k );
                }
                else
                {
                  g.set_control_point( c * f.get_control_point( k, from ), k, i );
                }
              }
            }
            f = g;
            changed = true;
          }
        }
      }

      // The cells, as their corners, of [uvstart, uvend] where Su . d and Sv . d, given as du and dv over the cell,
      // may both vanish, the distance to the point being stationary there, found by splitting the cell nsplit times.  A cell is passed
      // over where either's control points are all of one sign.  Where one is zero to within tol over the whole
      // cell, the distance does not change in its direction, so every point across it is as good as any other: the
      // cell is split only in the other direction, and the middle stands for them all.  A cell where both are zero
      // is one point, its middle.
      template<typename onedsurf__>
      void findstationary( std::vector< std::pair< std::pair< typename onedsurf__::data_type, typename onedsurf__::data_type >,
                                                   std::pair< typename onedsurf__::data_type, typename onedsurf__::data_type > > > & cells,
                           const std::pair< typename onedsurf__::data_type, typename onedsurf__::data_type > &uvstart,
                           const std::pair< typename onedsurf__::data_type, typename onedsurf__::data_type > &uvend,
                           const onedsurf__ &du, const onedsurf__ &dv, const typename onedsurf__::data_type &tol,
                           const typename onedsurf__::index_type &nsplit )
      {
        typedef typename onedsurf__::data_type data_type;
        typedef std::pair< data_type, data_type > uvpair;

        bool zu = du.allnear( tol );
        bool zv = dv.allnear( tol );
        if ( ( !zu && !du.straddles( tol ) ) || ( !zv && !dv.straddles( tol ) ) )
        {
          return;
        }

        uvpair uvmid = std::make_pair( (uvstart.first + uvend.first) * 0.5, (uvstart.second + uvend.second) * 0.5 );
        if ( ( zu && zv ) || nsplit <= 0 )
        {
          cells.push_back( std::make_pair( uvstart, uvend ) );
          return;
        }

        onedsurf__ dulo, duhi, dvlo, dvhi;
        if ( zv )
        {
          du.split_u( dulo, duhi, 0.5 );
          dv.split_u( dvlo, dvhi, 0.5 );
          findstationary( cells, uvstart, std::make_pair( uvmid.first, uvend.second ), dulo, dvlo, tol, nsplit - 1 );
          findstationary( cells, std::make_pair( uvmid.first, uvstart.second ), uvend, duhi, dvhi, tol, nsplit - 1 );
          return;
        }
        if ( zu )
        {
          du.split_v( dulo, duhi, 0.5 );
          dv.split_v( dvlo, dvhi, 0.5 );
          findstationary( cells, uvstart, std::make_pair( uvend.first, uvmid.second ), dulo, dvlo, tol, nsplit - 1 );
          findstationary( cells, std::make_pair( uvstart.first, uvmid.second ), uvend, duhi, dvhi, tol, nsplit - 1 );
          return;
        }

        onedsurf__ dull, dulh, duhl, duhh, dvll, dvlh, dvhl, dvhh;
        du.split_u( dulo, duhi, 0.5 );
        dulo.split_v( dull, dulh, 0.5 );
        duhi.split_v( duhl, duhh, 0.5 );
        dv.split_u( dvlo, dvhi, 0.5 );
        dvlo.split_v( dvll, dvlh, 0.5 );
        dvhi.split_v( dvhl, dvhh, 0.5 );

        findstationary( cells, uvstart, uvmid, dull, dvll, tol, nsplit - 1 );
        findstationary( cells, std::make_pair( uvstart.first, uvmid.second ), std::make_pair( uvmid.first, uvend.second ),
                        dulh, dvlh, tol, nsplit - 1 );
        findstationary( cells, std::make_pair( uvmid.first, uvstart.second ), std::make_pair( uvend.first, uvmid.second ),
                        duhl, dvhl, tol, nsplit - 1 );
        findstationary( cells, uvmid, uvend, duhh, dvhh, tol, nsplit - 1 );
      }
    }
  }
}
#endif
