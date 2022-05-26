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

#ifndef intersect_segment_surface_test_suite_hpp
#define intersect_segment_surface_test_suite_hpp

#include <cmath>
#include <typeinfo>
#include "eli/util/tolerance.hpp"

#include "eli/geom/surface/piecewise.hpp"
#include "eli/geom/surface/piecewise_body_of_revolution_creator.hpp"
#include "eli/geom/intersect/intersect_segment_triangle.hpp"
#include "eli/geom/intersect/intersect_segment_surface.hpp"

template<typename data__>
class intersect_segment_surface_test_suite : public Test::Suite
{
  private:
    typedef eli::geom::surface::piecewise<eli::geom::surface::bezier, data__, 3> piecewise_surface_type;
    typedef typename piecewise_surface_type::surface_type surface_type;
    typedef typename piecewise_surface_type::point_type point_type;
    typedef typename piecewise_surface_type::data_type data_type;
    typedef typename piecewise_surface_type::index_type index_type;
    typedef typename piecewise_surface_type::tolerance_type tolerance_type;
    typedef typename piecewise_surface_type::piecewise_curve_type piecewise_curve_type;
    typedef typename piecewise_curve_type::curve_type curve_type;
    typedef typename curve_type::control_point_type curve_control_point_type;

    tolerance_type tol;

  protected:
    void AddTests(const float &)
    {
      // add the tests
      TEST_ADD( intersect_segment_surface_test_suite<float>::segment_surface_test);
      TEST_ADD( intersect_segment_surface_test_suite<float>::surface_inside_test);
    }
    void AddTests(const double &)
    {
      // add the tests
      TEST_ADD( intersect_segment_surface_test_suite<double>::segment_surface_test);
      TEST_ADD( intersect_segment_surface_test_suite<double>::surface_inside_test);
    }
    void AddTests(const long double &)
    {
      // add the tests
      TEST_ADD( intersect_segment_surface_test_suite<long double>::segment_surface_test);
      TEST_ADD( intersect_segment_surface_test_suite<long double>::surface_inside_test);
    }

  public:
    intersect_segment_surface_test_suite()
    {
      AddTests(data__());
    }
    ~intersect_segment_surface_test_suite()
    {
    }

  private:
    void segment_surface_test()
    {

      piecewise_curve_type pc;
      curve_type c(3);
      curve_control_point_type cp[4];
      data_type k=eli::constants::math<data_type>::cubic_bezier_circle_const()*(eli::constants::math<data_type>::sqrt_two()-1);
      index_type i;

      // create curve
      cp[0] << 1, 0, 0;
      cp[1] << 1, k, 0;
      cp[2] << k, 1, 0;
      cp[3] << 0, 1, 0;
      for (i=0; i<4; ++i)
      {
        c.set_control_point(cp[i], i);
      }
      TEST_ASSERT(pc.push_back(c, 0.25)==piecewise_curve_type::NO_ERRORS);

      // set 2nd quadrant curve
      cp[0] <<  0, 1, 0;
      cp[1] << -k, 1, 0;
      cp[2] << -1, k, 0;
      cp[3] << -1, 0, 0;
      for (i=0; i<4; ++i)
      {
        c.set_control_point(cp[i], i);
      }
      TEST_ASSERT(pc.push_back(c, 0.25)==piecewise_curve_type::NO_ERRORS);

      piecewise_surface_type ps;

      TEST_ASSERT(eli::geom::surface::create_body_of_revolution(ps, pc, 0, true));

      std::vector < data_type > tvec;

      point_type p0, vec;
      data_type t_ref;
      p0 << 0, 0, 0;
      vec << 0, 10, 0;

      eli::geom::intersect::intersect_segment( tvec, ps, p0, vec );

      // std::cout << std::scientific;
      // std::cout << "First test size: " << tvec.size() << std::endl;
      // for ( i = 0; i < tvec.size(); i++ )
      // {
      //     std::cout << i << " " << tvec[i] << std::endl;
      // }

      TEST_ASSERT( tvec.size() == 1 );
      if ( tvec.size() >= 1 )
      {
        t_ref = 0.1;
        TEST_ASSERT( tol.approximately_equal( tvec[0], t_ref ) );
      }

      p0 << 0, -2, 0;

      eli::geom::intersect::intersect_segment( tvec, ps, p0, vec );

      // std::cout << "Second test size: " << tvec.size() << std::endl;
      // for ( i = 0; i < tvec.size(); i++ )
      // {
      //   std::cout << i << " " << tvec[i] << std::endl;
      // }

      TEST_ASSERT( tvec.size() == 2 );
      if ( tvec.size() >= 1 )
      {
        t_ref = 0.1;
        TEST_ASSERT( tol.approximately_equal( tvec[ 0 ], t_ref ));
      }

      if ( tvec.size() >= 2 )
      {
        t_ref = 0.3;
        TEST_ASSERT( tol.approximately_equal( tvec[ 1 ], t_ref ));
      }
    }

    void surface_inside_test()
    {
      piecewise_curve_type pc;
      curve_type c(3);
      curve_control_point_type cp[4];
      data_type k=eli::constants::math<data_type>::cubic_bezier_circle_const()*(eli::constants::math<data_type>::sqrt_two()-1);
      index_type i;

      // create curve
      cp[0] << 1, 0, 0;
      cp[1] << 1, k, 0;
      cp[2] << k, 1, 0;
      cp[3] << 0, 1, 0;
      for (i=0; i<4; ++i)
      {
        c.set_control_point(cp[i], i);
      }
      TEST_ASSERT(pc.push_back(c, 0.25)==piecewise_curve_type::NO_ERRORS);

      // set 2nd quadrant curve
      cp[0] <<  0, 1, 0;
      cp[1] << -k, 1, 0;
      cp[2] << -1, k, 0;
      cp[3] << -1, 0, 0;
      for (i=0; i<4; ++i)
      {
        c.set_control_point(cp[i], i);
      }
      TEST_ASSERT(pc.push_back(c, 0.25)==piecewise_curve_type::NO_ERRORS);

      piecewise_surface_type ps;

      TEST_ASSERT(eli::geom::surface::create_body_of_revolution(ps, pc, 0, true));

      point_type p0;
      p0 << 0, 0, 0;

      TEST_ASSERT( eli::geom::intersect::inside( ps, p0 ) );

      data_type r0, s0, t0;

      index_type nr = 11;
      index_type ns = 9;
      index_type nt = 7;

      data_type rmin( -1.5 ), rmax( 1.5 );
      data_type smin( -1.5 ), smax( 1.5 );
      data_type tmin( -1.5 ), tmax( 1.5 );

      data_type dr = ( rmax - rmin ) / ( nr - 1 );
      data_type ds = ( smax - smin) /  ( ns - 1 );
      data_type dt = ( tmax - tmin) /  ( nt - 1 );

      index_type in_cnt = 0;
      index_type out_cnt = 0;

      // Find array of points, starting from center every time.
      r0 = rmin;
      for ( index_type ir = 0; ir < nr; ir++ )
      {
        s0 = smin;
        for ( index_type is = 0; is < ns; is++ )
        {
          t0 = tmin;
          for ( index_type it = 0; it < nt; it++ )
          {
            data_type rad = sqrt( r0*r0 + s0*s0 + t0*t0 );

            // Avoid points _on_ the sphere as our piecewise surface is not exactly spherical.
            if ( std::abs(rad-1.0) > 0.01 )
            {
                bool in_ref = rad < 1.0;

                if ( in_ref )
                    in_cnt++;
                else
                    out_cnt++;

                p0 << r0, s0, t0;
                TEST_ASSERT( in_ref == eli::geom::intersect::inside( ps, p0 ));
            }
            t0 = t0 + dt;
            if ( t0 > tmax ) t0 = tmax;
          }
          s0 = s0 + ds;
          if ( s0 > smax ) s0 = smax;
        }
        r0 = r0 + dr;
        if ( r0 > rmax ) r0 = rmax;
      }

      // std::cout << "in: " << in_cnt << std::endl;
      // std::cout << "out: " << out_cnt << std::endl;

      // Arbitrary point on surface
      data_type u = 0.1234;
      data_type v = 0.5678;

      p0 = ps.f( u, v );
      point_type n = ps.normal( u, v );

      // Plane approximation tolerance is 1e-5
      for ( data_type i = 1.0; i < 4; i++ )
      {
        point_type p;
        p = p0 + std::pow( 10.0, -1.0 * i ) * n;
        // Outside
        TEST_ASSERT( eli::geom::intersect::inside( ps, p ) );

        p = p0 - std::pow( 10.0, -1.0 * i ) * n;
        // Inside
        TEST_ASSERT( !eli::geom::intersect::inside( ps, p ) );
      }


    }
};

#endif

