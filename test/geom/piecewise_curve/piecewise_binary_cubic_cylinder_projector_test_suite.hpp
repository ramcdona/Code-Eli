/*********************************************************************************
* Copyright (c) 2021 Rob McDonald
*
* All rights reserved. This program and the accompanying materials
* are made available under the terms of the Eclipse Public License v1.0
* which accompanies this distribution, and is available at
* http://www.eclipse.org/legal/epl-v10.html
*
********************************************************************************/

#ifndef piecewise_binary_cubic_circle_projector_test_suite_hpp
#define piecewise_binary_cubic_circle_projector_test_suite_hpp

#include <cmath>    // std::pow, std::exp

#include <typeinfo> // typeid
#include <string>   // std::string
#include <sstream>  // std::stringstream
#include <iomanip>  // std::setw
#include <limits>   // std::numeric_limits

#include "eli/constants/math.hpp"
#include "eli/geom/curve/piecewise.hpp"
#include "eli/geom/curve/piecewise_binary_cubic_cylinder_projector.hpp"
#include "eli/geom/intersect/minimum_distance_curve.hpp"

#include "eli/geom/curve/piecewise_four_digit_creator.hpp"
#include "eli/geom/curve/pseudo/four_digit.hpp"

#include "octave_helpers.hpp"

template<typename data__>
class piecewise_binary_cubic_cylinder_projector_test_suite : public Test::Suite
{
  private:
    typedef eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 3> piecewise_curve_type;
    typedef typename piecewise_curve_type::curve_type curve_type;
    typedef typename piecewise_curve_type::point_type point_type;
    typedef typename piecewise_curve_type::data_type data_type;
    typedef typename piecewise_curve_type::index_type index_type;
    typedef typename piecewise_curve_type::tolerance_type tolerance_type;

    typedef eli::geom::curve::piecewise_binary_cubic_cylinder_projector<data__, 3, tolerance_type> binary_projector_type;
    typedef eli::geom::curve::piecewise_binary_cubic_creator<data__, 3, tolerance_type> binary_creator_type;
    typedef eli::geom::curve::piecewise_four_digit_creator<data__, 3, tolerance_type> four_digit_type;
    typedef eli::geom::curve::piecewise_circle_creator<data__, 3, tolerance_type> circle_creator_type;

    typedef typename piecewise_curve_type::onedpiecewisecurve oned_type;
    typedef typename piecewise_curve_type::onedbezcurve onedbezcurve;


    tolerance_type tol;

  protected:
    void AddTests(const float &)
    {
      // add the tests
      TEST_ADD(piecewise_binary_cubic_cylinder_projector_test_suite<float>::circle_test);
      TEST_ADD(piecewise_binary_cubic_cylinder_projector_test_suite<float>::spiral_test);
      TEST_ADD(piecewise_binary_cubic_cylinder_projector_test_suite<float>::foil_test);
      TEST_ADD(piecewise_binary_cubic_cylinder_projector_test_suite<float>::circle_test2);
    }
    void AddTests(const double &)
    {
      // add the tests
      TEST_ADD(piecewise_binary_cubic_cylinder_projector_test_suite<double>::circle_test);
      TEST_ADD(piecewise_binary_cubic_cylinder_projector_test_suite<double>::spiral_test);
      TEST_ADD(piecewise_binary_cubic_cylinder_projector_test_suite<double>::foil_test);
      TEST_ADD(piecewise_binary_cubic_cylinder_projector_test_suite<double>::circle_test2);
    }
    void AddTests(const long double &)
    {
      // add the tests
      TEST_ADD(piecewise_binary_cubic_cylinder_projector_test_suite<long double>::circle_test);
      TEST_ADD(piecewise_binary_cubic_cylinder_projector_test_suite<long double>::spiral_test);
      TEST_ADD(piecewise_binary_cubic_cylinder_projector_test_suite<long double>::foil_test);
      TEST_ADD(piecewise_binary_cubic_cylinder_projector_test_suite<long double>::circle_test2);
    }

  public:
    piecewise_binary_cubic_cylinder_projector_test_suite() : tol()
    {
      AddTests(data__());
    }
    ~piecewise_binary_cubic_cylinder_projector_test_suite()
    {
    }

  private:

    void circle_test()
    {
      binary_projector_type bp;
      piecewise_curve_type pc, pcout;
      curve_type c;
      point_type cp[2];

      data_type r(1.0);

      data_type tol(1e-6);

      cp[0] << 0., 1., 0.;
      cp[1] << 0., 1., 2.*eli::constants::math<data__>::pi();

      c.resize(1);
      for (index_type i=0; i<2; ++i)
      {
          c.set_control_point( cp[i], i);
      }

      pc.push_back( c, 1.0 );

      bp.setup( pc, r, tol, 2, 15 );
      bp.create( pcout );

      // pcout.octave_print( 1 );

      // Check evenly around circle.
      index_type n(101);
      for ( index_type i = 0; i <= n; i++ )
      {
        data_type t = (data_type) i / (data_type) n;
        point_type p0 = pc.f(t);
        data_type theta = p0.z() / r;
        point_type p0t;
        p0t << p0.x(), r * std::cos(theta), r * std::sin(theta);
        point_type p0c = pcout.f(t);

        TEST_ASSERT( ( p0t - p0c ).norm() < tol );
        if ( ( p0t - p0c ).norm() >= tol  )
        {
          std::cout << i << " " << t << " " << ( p0t - p0c ).norm() << std::endl <<
          " x " << p0t.x() << " " << p0c.x() << std::endl <<
          " y " << p0t.y() << " " << p0c.y() << std::endl <<
          " z " << p0t.z() << " " << p0c.z() << std::endl << std::endl;
        }
      }

      data_type t0 = pcout.get_t0();
      data_type dt;
      pcout.get( c, dt, 0 );

      // Check just one segment.  Use same tolerance as was used to adapt curve.
      // Hopefully it presented a worst case scenario.
      for ( index_type i = 0; i <= n; i++ )
      {
        data_type t = t0 + dt * (data_type) i / (data_type) n;
        point_type p0 = pc.f(t);
        data_type theta = p0.z() / r;
        point_type p0t;
        p0t << p0.x(), r * std::cos(theta), r * std::sin(theta);
        point_type p0c = pcout.f(t);

        TEST_ASSERT( ( p0t - p0c ).norm() < tol );
        if ( ( p0t - p0c ).norm() >= tol  )
        {
          std::cout << i << " " << t << " " << ( p0t - p0c ).norm() << std::endl <<
          " x " << p0t.x() << " " << p0c.x() << std::endl <<
          " y " << p0t.y() << " " << p0c.y() << std::endl <<
          " z " << p0t.z() << " " << p0c.z() << std::endl << std::endl;
        }
      }
    }

    void spiral_test()
    {
      binary_projector_type bp;
      piecewise_curve_type pc, pcout;
      curve_type c;
      point_type cp[2];

      data_type r(1.0);

      data_type tol(1e-3);

      cp[0] << 0., 1., 0.;
      cp[1] << 3., 1., 3.*2.*eli::constants::math<data__>::pi();

      c.resize(1);
      for (index_type i=0; i<2; ++i)
      {
          c.set_control_point( cp[i], i);
      }

      pc.push_back( c, 1.0 );

      bp.setup( pc, r, tol, 1, 15 );
      bp.create( pcout );

      // Check evenly around curve.  Use very fine tolerance because curve is exact at critical points.
      index_type n(101);
      for ( index_type i = 0; i <= n; i++ )
      {
        data_type t = (data_type) i / (data_type) n;
        point_type p0 = pc.f(t);
        data_type theta = p0.z() / r;
        point_type p0t;
        p0t << p0.x(), r * std::cos(theta), r * std::sin(theta);
        point_type p0c = pcout.f(t);

        TEST_ASSERT( ( p0t - p0c ).norm() < tol );
        if ( ( p0t - p0c ).norm() >= tol  )
        {
          std::cout << i << " " << t << " " << ( p0t - p0c ).norm() << std::endl <<
          " x " << p0t.x() << " " << p0c.x() << std::endl <<
          " y " << p0t.y() << " " << p0c.y() << std::endl <<
          " z " << p0t.z() << " " << p0c.z() << std::endl << std::endl;
        }
      }

      data_type t0 = pcout.get_t0();
      data_type dt;
      pcout.get( c, dt, 0 );

      // Check just one segment.  Use same tolerance as was used to adapt curve.
      // Hopefully it presented a worst case scenario.
      for ( index_type i = 0; i <= n; i++ )
      {
        data_type t = t0 + dt * (data_type) i / (data_type) n;
        point_type p0 = pc.f(t);
        data_type theta = p0.z() / r;
        point_type p0t;
        p0t << p0.x(), r * std::cos(theta), r * std::sin(theta);
        point_type p0c = pcout.f(t);

        TEST_ASSERT( ( p0t - p0c ).norm() < tol );
        if ( ( p0t - p0c ).norm() >= tol  )
        {
          std::cout << i << " " << t << " " << ( p0t - p0c ).norm() << std::endl <<
          " x " << p0t.x() << " " << p0c.x() << std::endl <<
          " y " << p0t.y() << " " << p0c.y() << std::endl <<
          " z " << p0t.z() << " " << p0c.z() << std::endl << std::endl;
        }
      }
    }

    void foil_test()
    {
      binary_projector_type bp;
      binary_creator_type bc;
      piecewise_curve_type pc, pcc, pcout;
      curve_type c;

      four_digit_type af;

      data_type tc, cam, cam_loc;
      bool rtn;

      data_type tol(1e-6);

      // set airfoil thickness
      tc = 0.24;
      rtn = af.set_thickness( tc );
      TEST_ASSERT( rtn );

      // set airfoil camber
      cam = .04;
      cam_loc = .2;
      rtn = af.set_camber( cam, cam_loc );
      TEST_ASSERT( rtn );

      af.set_sharp_trailing_edge( true );

      bool fit_success;
      fit_success = af.create( pc );
      TEST_ASSERT( fit_success );

      typename piecewise_curve_type::rotation_matrix_type rmat;
      data_type th = 0.5 * eli::constants::math<data__>::pi();
      // Rotate 90 about Y.
      rmat << std::cos(th), 0, -std::sin(th),
              0, 1., 0,
              std::sin(th), 0, std::cos(th);
      pc.rotate( rmat );
      // Rotate 90 about Z.
      rmat << std::cos(th), -std::sin(th), 0,
              std::sin(th), std::cos(th), 0,
              0, 0, 1.;
      pc.rotate( rmat );

      point_type dx;
      dx << 1.0, 0.0, 1.0;
      pc.translate( dx );

      // Convert to binary cubic curve.
      bc.setup( pc, tol, 1, 15 );
      index_type dpcc = bc.create( pcc );
      TEST_ASSERT( dpcc <= 15 );

      data_type r(1.0);

      // pcc.octave_print( 1 );

      bp.setup( pcc, r, tol, 1, 15 );
      index_type d = bp.create( pcout );

      TEST_ASSERT( d <= 15 );

      // std::cout << "Recursion depth: " << dpcc << " and " << d << std::endl;
      // std::cout << "Planar nseg: " << pcc.number_segments() << " cylindrical nseg: " << pcout.number_segments() << std::endl;

      // pcout.octave_print( 2 );

      // Check evenly around curve.
      index_type n(101);
      for ( index_type i = 0; i <= n; i++ )
      {
        data_type t = (data_type) i / (data_type) n;
        point_type p0 = pcc.f(t);
        data_type theta = p0.z() / r;
        point_type p0t;
        p0t << p0.x(), r * std::cos(theta), r * std::sin(theta);
        point_type p0c = pcout.f(t);

        TEST_ASSERT( ( p0t - p0c ).norm() < tol );
        if ( ( p0t - p0c ).norm() >= tol  )
        {
          std::cout << i << " " << t << " " << ( p0t - p0c ).norm() << std::endl <<
          " x " << p0t.x() << " " << p0c.x() << std::endl <<
          " y " << p0t.y() << " " << p0c.y() << std::endl <<
          " z " << p0t.z() << " " << p0c.z() << std::endl << std::endl;
        }
      }

      data_type t0 = pcout.get_t0();
      data_type dt;
      pcout.get( c, dt, 0 );

      // Check just one segment.  Use same tolerance as was used to adapt curve.
      // Hopefully it presented a worst case scenario.
      for ( index_type i = 0; i <= n; i++ )
      {
        data_type t = t0 + dt * (data_type) i / (data_type) n;
        point_type p0 = pcc.f(t);
        data_type theta = p0.z() / r;
        point_type p0t;
        p0t << p0.x(), r * std::cos(theta), r * std::sin(theta);
        point_type p0c = pcout.f(t);

        TEST_ASSERT( ( p0t - p0c ).norm() < tol );
        if ( ( p0t - p0c ).norm() >= tol  )
        {
          std::cout << i << " " << t << " " << ( p0t - p0c ).norm() << std::endl <<
          " x " << p0t.x() << " " << p0c.x() << std::endl <<
          " y " << p0t.y() << " " << p0c.y() << std::endl <<
          " z " << p0t.z() << " " << p0c.z() << std::endl << std::endl;
        }
      }

    }

    void circle_test2()
    {
      binary_projector_type projector;
      binary_creator_type approximator;
      circle_creator_type circle_creator;
      piecewise_curve_type c, pc, pc2;
      point_type origin, x, y;
      data_type radius;

      // set the parameters for circle
      origin << 1, 1, 1;
      x << 0, 0, 1;
      y << 1, 0, 0;
      radius=1;

      circle_creator.set( origin, x, y, radius );

      // create the circle
      circle_creator.create(c);

      typename piecewise_curve_type::rotation_matrix_type rmat;

      data_type th = 0.321 * eli::constants::math<data__>::pi();
      // Rotate 90 about Y.
      rmat << std::cos(th), 0, -std::sin(th),
              0, 1., 0,
              std::sin(th), 0, std::cos(th);
      c.rotate( rmat );

      data_type tol(1e-6);

      approximator.setup( c, tol, 1, 15 );
      index_type dpcc = approximator.create( pc );
      TEST_ASSERT( dpcc <= 15 );

      // pc.octave_print( 1 );

      data_type r(1.0);

      projector.setup( pc, r, tol, 1, 15 );
      index_type d = projector.create( pc2 );
      TEST_ASSERT( d <= 15 );

      // pc2.octave_print( 2 );

      // std::cout << "Recursion depth: " << d << std::endl;
      // std::cout << "Planar nseg: " << pc.number_segments() << " cylindrical nseg: " << pc2.number_segments() << std::endl;
    }

};

#endif

