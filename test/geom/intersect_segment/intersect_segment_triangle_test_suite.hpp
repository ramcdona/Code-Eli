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

#ifndef intersect_segment_triangle_test_suite_hpp
#define intersect_segment_triangle_test_suite_hpp

#include <cmath>    // cos(), sin()

#include <typeinfo> // typeid

#include "eli/util/tolerance.hpp"

#include "eli/geom/intersect/intersect_segment_triangle.hpp"

template<typename data__>
class intersect_segment_triangle_test_suite : public Test::Suite
{
  private:
    typedef data__ data_type;
    typedef Eigen::Matrix<data_type, 1, 2> point_type2;
    typedef Eigen::Matrix<data_type, 1, 3> point_type3;
    typedef eli::util::tolerance<data_type> tolerance_type;

    tolerance_type tol;

  protected:
    void AddTests(const float &)
    {
      // add the tests
      TEST_ADD( intersect_segment_triangle_test_suite<float>::segment_triangle_test);
    }
    void AddTests(const double &)
    {
      // add the tests
      TEST_ADD( intersect_segment_triangle_test_suite<double>::segment_triangle_test);
    }
    void AddTests(const long double &)
    {
      // add the tests
      TEST_ADD( intersect_segment_triangle_test_suite<long double>::segment_triangle_test);
    }

  public:
    intersect_segment_triangle_test_suite()
    {
      AddTests(data__());
    }
    ~intersect_segment_triangle_test_suite()
    {
    }

  private:
    void segment_triangle_test()
    {
      point_type3 a0, a1, a2, p0, p1;
      data_type t, u, w;
      data_type t_ref, u_ref, w_ref;
      bool isect;

      a0 << 0, 0, 0;
      a1 << 0, 1, 0;
      a2 << 1, 0, 0;

      p0 << 0.25, 0.25, -1;
      p1 << 0, 0, 2;

      isect = eli::geom::intersect::seg_tri_intersect( a0, a1, a2, p0, p1, u, w, t );

      TEST_ASSERT( isect );
      t_ref = 0.5;
      u_ref = 0.25;
      w_ref = 0.25;
      TEST_ASSERT( tol.approximately_equal( t, t_ref ) );
      TEST_ASSERT( tol.approximately_equal( u, u_ref ) );
      TEST_ASSERT( tol.approximately_equal( w, w_ref ) );

      p0 << 1, 1, -1;

      isect = eli::geom::intersect::seg_tri_intersect( a0, a1, a2, p0, p1, u, w, t );

      TEST_ASSERT( !isect );
      t_ref = 0.5;
      u_ref = 1;
      w_ref = 1;
      TEST_ASSERT( tol.approximately_equal( t, t_ref ) );
      TEST_ASSERT( tol.approximately_equal( u, u_ref ) );
      TEST_ASSERT( tol.approximately_equal( w, w_ref ) );


      p0 << 0.5, 0.5, -1;

      isect = eli::geom::intersect::seg_tri_intersect( a0, a1, a2, p0, p1, u, w, t );

      TEST_ASSERT( isect );
      t_ref = 0.5;
      u_ref = 0.5;
      w_ref = 0.5;
      TEST_ASSERT( tol.approximately_equal( t, t_ref ) );
      TEST_ASSERT( tol.approximately_equal( u, u_ref ) );
      TEST_ASSERT( tol.approximately_equal( w, w_ref ) );

      // std::cout << std::scientific;
      // std::cout << isect << std::endl;
      // std::cout << t-t_ref << std::endl;
      // std::cout << u-u_ref << std::endl;
      // std::cout << w-w_ref << std::endl;
    }
};

#endif

