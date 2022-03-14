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

#ifndef eli_geom_curve_piecewise_binary_cubic_cylinder_projector_hpp
#define eli_geom_curve_piecewise_binary_cubic_cylinder_projector_hpp

#include <iterator>
#include <vector>

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
      template<typename data__, unsigned short dim__, typename tol__>
      class piecewise_binary_cubic_cylinder_projector
      {
        public:
          typedef data__  data_type;
          typedef Eigen::Matrix<data_type, 1, dim__> point_type;
          typedef typename point_type::Index index_type;
          typedef tol__ tolerance_type;

          typedef piecewise<bezier, data_type, dim__, tolerance_type> piecewise_curve_type;
          typedef typename piecewise_curve_type::curve_type curve_type;

          piecewise_binary_cubic_cylinder_projector() {}

          void setup(const piecewise_curve_type &pc, const data_type &rad, const data_type &t, const index_type &dmin, const index_type &dmax)
          {
            parent_curve = pc;
            parent_curve.get_pmap( parent_pmap );

            ttol = t;
            atol = -1.0;
            r = rad;

            min_depth = dmin;
            max_depth = dmax;
          }

          void setup(const piecewise_curve_type &pc, const data_type &rad, const data_type &t, const data_type &a, const index_type &dmin, const index_type &dmax)
          {
            setup( pc, rad, t, dmin, dmax );
            atol = a;
          }

          virtual index_type corner_create(piecewise<bezier, data_type, dim__, tolerance_type> &pc) const
          {
            std::vector<data_type> tdisc;

            // Sometimes includes first/last point.  Sometimes doesn't.
            parent_curve.find_discontinuities( atol, tdisc );
            // Always append last point, repeated points don't hurt corner_create.
            // also doesn't matter whether corner_create has first point or not.
            tdisc.push_back( parent_curve.get_tmax() );

            return corner_create( pc, tdisc );
          }

          virtual index_type corner_create(piecewise<bezier, data_type, dim__, tolerance_type> &pc, const std::vector<data_type> &tdisc ) const
          {
            index_type d = -1;
            point_type p0, m01, m02, p1, m11, m12;
            data_type theta0, theta1;
            point_type p0t, m01t, m02t, p1t, m11t, m12t;

            data_type t0, t1;
            t0 = parent_curve.get_t0();

            pc.clear();

            // set the start parameter
            pc.set_t0( t0 );

            p0 = parent_curve.f(t0);
            parent_curve.fps(t0, m01, m02);
            theta0 = calc_theta( p0 );
            transform_point_slopes( p0, theta0, m01, m02, p0t, m01t, m02t );

            for ( typename std::vector< data_type>::size_type i = 0; i < tdisc.size(); i++ )
            {
              t1 = tdisc[i];
              if ( t1 > t0 )
              {
                p1 = parent_curve.f(t1);
                parent_curve.fps(t1, m11, m12);
                theta1 = calc_theta( p1 );
                transform_point_slopes( p1, theta1, m11, m12, p1t, m11t, m12t );

                // Build approximate curve.
                curve_type c;
                c = make_curve_point_slope(p0t, m02t, p1t, m11t, t1-t0);

                pc.push_back(c, t1-t0);

                d = adapt_pc( pc, t0, p0t, m02t, theta0, t1, p1t, m11t, theta1 );

                t0 = t1;
                p0t = p1t;
                theta0 = theta1;
                m01t = m11t;
                m02t = m12t;
              }
            }
            return d;
          }

          virtual index_type create(piecewise<bezier, data_type, dim__, tolerance_type> &pc) const
          {
            index_type d = -1;
            data_type t0, t1;
            t0 = parent_curve.get_t0();
            t1 = parent_curve.get_tmax();

            point_type p0, m0, p1, m1;
            data_type theta0, theta1;
            point_type p0t, m0t, p1t, m1t;
            p0 = parent_curve.f(t0);
            m0 = parent_curve.fp(t0);
            theta0 = calc_theta( p0 );
            transform_point_slope( p0, theta0, m0, p0t, m0t );
            p1 = parent_curve.f(t1);
            m1 = parent_curve.fp(t1);
            theta1 = calc_theta( p1 );
            transform_point_slope( p1, theta1, m1, p1t, m1t );

            pc.clear();

            // set the start parameter
            pc.set_t0( t0 );

            // Build approximate curve.
            curve_type c;
            c = make_curve_point_slope(p0t, m0t, p1t, m1t, t1-t0);

            pc.push_back(c, t1-t0);

            d = adapt_pc( pc, t0, p0t, m0t, theta0, t1, p1t, m1t, theta1 );

            return d;
          }

        protected:

          index_type adapt_pc( piecewise<bezier, data_type, dim__, tolerance_type> &pc, const data_type &t0, const point_type &p0, const point_type &m0, const data_type &theta0,
                         const data_type &t1, const point_type &p1, const point_type &m1, const data_type &theta1, index_type depth = 0 ) const
          {
            data_type tmid = ( t0 + t1 ) / 2.0;
            data_type t25 = t0 + 0.25 * ( t1 - t0 ) ;
            data_type t75 = t0 + 0.75 * ( t1 - t0 ) ;
            point_type pmid0 = parent_curve.f( tmid );
            data_type thetamid = calc_theta( pmid0 );
            point_type pmid = transform_point( pmid0, thetamid );

            bool adapt = false;

            if ( depth < min_depth )
            {
              adapt = true;
            }
            else if ( eli::geom::point::distance( pmid, pc.f( tmid ) ) > ttol )
            {
              adapt = true;
            }
            else if ( eli::geom::point::distance( transform_point( parent_curve.f( t25 ) ), pc.f( t25 ) ) > ttol )
            {
                adapt = true;
            }
            else if ( eli::geom::point::distance( transform_point( parent_curve.f( t75 ) ), pc.f( t75 ) ) > ttol )
            {
                adapt = true;
            }
            else if ( !test_match( pc, t0, t1 ) )
            {
              adapt = true;
            }

            if ( adapt && depth < max_depth )
            {
              point_type mmid = transform_slope( parent_curve.fp( tmid ), thetamid );

              piecewise<bezier, data_type, dim__, tolerance_type> insert;

              insert.set_t0( t0 );

              curve_type c1;
              c1 = make_curve_point_slope(p0, m0, pmid, mmid, tmid-t0);
              insert.push_back(c1, tmid-t0);

              curve_type c2;
              c2 = make_curve_point_slope(pmid, mmid, p1, m1, t1-tmid);
              insert.push_back(c2, t1-tmid);

              pc.replace_t( insert, t0 );

              index_type d1 = adapt_pc( pc, t0, p0, m0, theta0, tmid, pmid, mmid, thetamid, depth + 1 );
              index_type d2 = adapt_pc( pc, tmid, pmid, mmid, thetamid, t1, p1, m1, theta1, depth + 1 );
              return std::max( d1, d2 );
            }
            else
            {
              return depth;
            }
          }

          bool test_match( const piecewise<bezier, data_type, dim__, tolerance_type> &pc, const data_type &tstart, const data_type &tend ) const
          {
            data_type tcheck;

            // Check start/end and midpoints of piecewise parent curve that overlap range.
            for ( typename std::vector< data_type>::size_type i = 0; i < parent_pmap.size() - 1; i++ )
            {
              if ( parent_pmap[ i ] <= tend && parent_pmap[ i + 1 ] >= tstart )
              {
                tcheck = parent_pmap[ i ];
                if ( tcheck >= tstart && tcheck <= tend )
                {
                  if ( eli::geom::point::distance( transform_point( parent_curve.f( tcheck ) ), pc.f( tcheck ) ) > ttol )
                  {
                    return false;
                  }
                }

                tcheck = ( parent_pmap[ i ] + parent_pmap[ i + 1 ] ) / 2.0;
                if ( tcheck >= tstart && tcheck <= tend )
                {
                  if ( eli::geom::point::distance( transform_point( parent_curve.f( tcheck ) ), pc.f( tcheck ) ) > ttol )
                  {
                    return false;
                  }
                }
              }
            }

            return true;
          }

          static curve_type make_curve_point_slope( const point_type &p0, const point_type &m0,
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

          data_type calc_theta( const point_type &p ) const
          {
              return p.z() / r;
          }

          point_type transform_point( const point_type &p, const data_type &theta ) const
          {
            return point_type( p.x(),
                               r * std::cos(theta),
                               r * std::sin(theta) );
          }

          point_type transform_point( const point_type &p ) const
          {
            return transform_point( p, calc_theta( p ) );
          }

          point_type transform_slope( const point_type &m, const data_type &theta ) const
          {
            data_type cth = std::cos(theta);
            data_type sth = std::sin(theta);

            return point_type( m.x(), - sth * m.z(), cth * m.z() );
          }

          void transform_point_slope( const point_type &p, const data_type &theta, const point_type &m, point_type &pt, point_type &mt ) const
          {
            data_type cth = std::cos(theta);
            data_type sth = std::sin(theta);
            pt << p.x(), r * cth, r * sth;
            mt << m.x(), - sth * m.z(), cth * m.z();
          }

          void transform_point_slopes( const point_type &p0, const data_type &theta, const point_type &m0, const point_type &m1, point_type &pt, point_type &m0t, point_type &m1t ) const
          {
            data_type cth = std::cos(theta);
            data_type sth = std::sin(theta);
            pt << p0.x(), r * cth, r * sth;
            m0t << m0.x(), - sth * m0.z(), cth * m0.z();
            m1t << m1.x(), - sth * m1.z(), cth * m1.z();
          }

        private:
          piecewise_curve_type parent_curve;

          std::vector < data_type > parent_pmap;

          data_type ttol;
          data_type atol;

          data_type r;

          index_type min_depth;
          index_type max_depth;

      };
    }
  }
}
#endif
