/*********************************************************************************
* Copyright (c) 2013 David D. Marshall <ddmarsha@calpoly.edu>
*
* All rights reserved. This program and the accompanying materials
* are made available under the terms of the Eclipse Public License v1.0
* which accompanies this distribution, and is available at
* http://www.eclipse.org/legal/epl-v10.html
*
* Contributors:
*    Rob McDonald - implementation of binary cubic creator
********************************************************************************/

#ifndef eli_geom_curve_piecewise_binary_cubic_creator_abstract_hpp
#define eli_geom_curve_piecewise_binary_cubic_creator_abstract_hpp

#include <iterator>
#include <vector>
#include <cmath>

#include "eli/code_eli.hpp"

#include "eli/geom/curve/piecewise.hpp"
#include "eli/geom/curve/bezier.hpp"
#include "eli/geom/point/distance.hpp"

namespace eli
{
  namespace geom
  {
    namespace curve
    {
      template<typename fun__, typename data__, unsigned short dim__, typename tol__>
      class piecewise_binary_cubic_creator_abstract
      {
        public:
          typedef fun__ function_type;
          typedef data__  data_type;
          typedef Eigen::Matrix<data_type, 1, dim__> point_type;
          typedef typename point_type::Index index_type;
          typedef tol__ tolerance_type;

          typedef piecewise<bezier, data_type, dim__, tolerance_type> piecewise_curve_type;
          typedef typename piecewise_curve_type::curve_type curve_type;

          piecewise_binary_cubic_creator_abstract() {}

          void setup(const function_type &f, const data_type &t0, const data_type &t1, const data_type &t, const index_type &dmin, const index_type &dmax)
          {
            fun = f;

            tmin = t0;
            tmax = t1;
            ttol = t;

            min_depth = dmin;
            max_depth = dmax;
          }

          virtual index_type create(piecewise<bezier, data_type, dim__, tolerance_type> &pc) const
          {
            index_type d = -1;
            data_type t0, t1;
            t0 = tmin;
            t1 = tmax;

            point_type p0, m0, p1, m1;
            fun(t0, p0, m0);
            fun(t1, p1, m1);

            pc.clear();

            // set the start parameter
            pc.set_t0( t0 );

            // Build approximate curve.
            curve_type c;
            c = make_curve_point_slope(p0, m0, p1, m1, t1-t0);

            pc.push_back(c, t1-t0);

            d = adapt_pc( pc, t0, p0, m0, t1, p1, m1);

            return d;
          }

        protected:

          index_type adapt_pc( piecewise<bezier, data_type, dim__, tolerance_type> &pc, const data_type &t0, const point_type &p0, const point_type &m0,
                         const data_type &t1, const point_type &p1, const point_type &m1, index_type depth = 0 ) const
          {
            data_type tmid = ( t0 + t1 ) / 2.0;
            point_type pmid;
            point_type mmid;

            fun( tmid, pmid, mmid );

            bool adapt = false;

            if ( depth < min_depth )
            {
              adapt = true;
            }
            else if ( eli::geom::point::distance( pmid, pc.f( tmid ) ) > ttol )
            {
              adapt = true;
            }

            if ( adapt && depth < max_depth )
            {
              piecewise<bezier, data_type, dim__, tolerance_type> insert;

              insert.set_t0( t0 );

              curve_type c1;
              c1 = make_curve_point_slope(p0, m0, pmid, mmid, tmid-t0);
              insert.push_back(c1, tmid-t0);

              curve_type c2;
              c2 = make_curve_point_slope(pmid, mmid, p1, m1, t1-tmid);
              insert.push_back(c2, t1-tmid);

              pc.replace_t( insert, t0 );

              index_type d1 = adapt_pc( pc, t0, p0, m0, tmid, pmid, mmid, depth + 1 );
              index_type d2 = adapt_pc( pc, tmid, pmid, mmid, t1, p1, m1, depth + 1 );
              return std::max( d1, d2 );
            }
            else
            {
              return depth;
            }
          }

          static curve_type make_curve_point_slope(const point_type &p0, const point_type &m0,
                                                   const point_type &p1, const point_type &m1, const data_type &dt)
          {
            curve_type c(3);
            point_type cp[4];

            cp[0]=p0;
            cp[1]=p0+(dt*m0/3.0);
            cp[2]=p1-(dt*m1/3.0);
            cp[3]=p1;

            for (index_type i=0; i<4; ++i)
            {
              c.set_control_point(cp[i], i);
            }

            return c;
          }

        private:
          function_type fun;

          data_type ttol;

          data_type tmin;
          data_type tmax;

          index_type min_depth;
          index_type max_depth;

      };
    }
  }
}
#endif
