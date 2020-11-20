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

#ifndef eli_geom_intersect_intersect_axis_bounding_box_hpp
#define eli_geom_intersect_intersect_axis_bounding_box_hpp

#include "eli/code_eli.hpp"

#include "eli/geom/point/distance.hpp"
#include "eli/geom/general/bounding_box.hpp"


namespace eli
{
  namespace geom
  {
    namespace intersect
    {
      template<typename data__, unsigned short dim__, typename tol__>
      data__ minimum_distance(const eli::geom::general::bounding_box<data__, dim__, tol__> &bb,
                              const typename eli::geom::general::bounding_box<data__, dim__, tol__>::point_type &pt,
                              const typename eli::geom::general::bounding_box<data__, dim__, tol__>::index_type iax)
      {
        data__ dist2(0), len;

        bool below_min, above_max;

        // for each dimension that is outside corresponding min/max add that to distance
        below_min = pt(0, iax) < bb.get_min()(0, iax);
        above_max = pt(0, iax) > bb.get_max()(0, iax);

        if (below_min)
        {
          len = bb.get_min()(0,iax) - pt(0,iax);
        }
        else if (above_max)
        {
          len = pt(0,iax) - bb.get_max()(0,iax);
        }
        else
        {
          len=0;
        }

        return len;
      }

      template<typename data__, unsigned short dim__, typename tol__>
      data__ maximum_distance(const eli::geom::general::bounding_box<data__, dim__, tol__> &bb,
                              const typename eli::geom::general::bounding_box<data__, dim__, tol__>::point_type &pt,
                              const typename eli::geom::general::bounding_box<data__, dim__, tol__>::index_type iax)
      {
        data__ dist2(0), len;
        bool below_min, above_max;

        // for each dimension that is outside corresponding min/max add that to distance
        below_min = pt(0, iax) < bb.get_min()(0, iax);
        above_max = pt(0, iax) > bb.get_max()(0, iax);

        if (below_min)
        {
          len = bb.get_max()(0,iax) - pt(0,iax);
        }
        else if (above_max)
        {
          len = pt(0,iax) - bb.get_min()(0,iax);
        }
        else
        {
          data__ l1, l2;
          l1 = pt(0,iax) - bb.get_min()(0,iax);
          l2 = bb.get_max()(0,iax) - pt(0,iax);

          if ( l1 >= l2 )
            len = l1;
          else
            len = l2;
        }

        return len;
      }

      template<typename data__, unsigned short dim__, typename tol__>
      void minmax_distance(const eli::geom::general::bounding_box<data__, dim__, tol__> &bb,
                           const typename eli::geom::general::bounding_box<data__, dim__, tol__>::point_type &pt,
                           const typename eli::geom::general::bounding_box<data__, dim__, tol__>::index_type iax,
                           data__ &dmin,
                           data__ &dmax)
      {
        bool below_min, above_max;

        // for each dimension that is outside corresponding min/max add that to distance
        below_min = pt(0, iax) < bb.get_min()(0, iax);
        above_max = pt(0, iax) > bb.get_max()(0, iax);

        if (below_min)
        {
          dmin = bb.get_min()(0,iax) - pt(0,iax);
          dmax = bb.get_max()(0,iax) - pt(0,iax);
        }
        else if (above_max)
        {
          dmin = pt(0,iax) - bb.get_max()(0,iax);
          dmax = pt(0,iax) - bb.get_min()(0,iax);
        }
        else
        {
          dmin=0;

          data__ l1, l2;
          l1 = pt(0,iax) - bb.get_min()(0,iax);
          l2 = bb.get_max()(0,iax) - pt(0,iax);

          if( l1 >= l2 )
            dmax = l1;
          else
            dmax = l2;
        }

      }

    }
  }
}
#endif
