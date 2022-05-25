/*********************************************************************************
* Copyright (c) 2022 Rob McDonald <rob.a.mcdonald@gmail.com>
*
* All rights reserved. This program and the accompanying materials
* are made available under the terms of the Eclipse Public License v1.0
* which accompanies this distribution, and is available at
* http://www.eclipse.org/legal/epl-v10.html
*
* Contributors:
*    Rob McDonald - initial code and implementation
********************************************************************************/

#ifndef ELI_INTERSECT_SEGMENT_TRIANGLE_HPP
#define ELI_INTERSECT_SEGMENT_TRIANGLE_HPP

#include <cmath>
#include <limits>

#include "eli/code_eli.hpp"

namespace eli
{
  namespace geom
  {
    namespace intersect
    {
      // ==== Triangle - Line Segment Intersection ====//
      // ==== A - Base Point on Triangle
      // ==== B - Vector for one   Side of Tri
      // ==== C - Vector for other Side of Tri
      // ==== D - Base Point for Line Seg
      // ==== E - Vector for Line Seg
      // ==============================================//
      template < typename data__, int dim__, typename tol__=eli::util::tolerance < data__ > >
      bool seg_tri_intersect( const Eigen::Matrix < data__, 1, dim__ > &A, const Eigen::Matrix < data__, 1, dim__ > &B,
                              const Eigen::Matrix < data__, 1, dim__ > &C, const Eigen::Matrix < data__, 1, dim__ > &D,
                              const Eigen::Matrix < data__, 1, dim__ > &E,
                              data__ &u, data__ &w, data__ &t )
      {
          typedef data__ data_type;
          typedef Eigen::Matrix < data__, 1, dim__ > point_type;

          data_type zero = -sqrt( std::numeric_limits < data_type >::epsilon() );
          data_type one = 1.0 - zero;

          point_type cs = B.cross( C );
          data_type denom = cs.dot( E );

          if ( std::abs( denom ) <= std::numeric_limits < data_type >::epsilon() )
          {
              return false;
          }

          t = ( cs.dot( A ) - cs.dot( D )) / denom;

          if (( t < zero ) || ( t > one ))
          {
              return false;
          }

          cs = C.cross( E );
          denom = cs.dot( B );

          if ( std::abs( denom ) <= std::numeric_limits < data_type >::epsilon() )
          {
              return false;
          }

          u = ( cs.dot( D ) - cs.dot( A )) / denom;

          if (( u < zero ) || ( u > one ))
          {
              return false;
          }

          cs = B.cross( E );
          denom = cs.dot( C );

          if ( std::abs( denom ) <= std::numeric_limits < data_type >::epsilon() )
          {
              return false;
          }

          w = ( cs.dot( D ) - cs.dot( A )) / denom;

          if (( w < zero ) || ( w > one ))
          {
              return false;
          }

          if (( w + u ) > one )
          {
              return false;
          }

          return true;
      }

    }
  }
}



#endif //ELI_INTERSECT_SEGMENT_TRIANGLE_HPP
