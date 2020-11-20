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


    }
  }
}
#endif
