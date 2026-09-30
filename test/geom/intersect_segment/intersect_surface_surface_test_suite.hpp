/*********************************************************************************
* Copyright (c) 2026 Rob McDonald <rob.a.mcdonald@gmail.com>
*
* All rights reserved. This program and the accompanying materials
* are made available under the terms of the Eclipse Public License v1.0
* which accompanies this distribution, and is available at
* http://www.eclipse.org/legal/epl-v10.html
*
* Contributors:
*    Rob McDonald - initial code and implementation
********************************************************************************/

#ifndef intersect_surface_surface_test_suite_hpp
#define intersect_surface_surface_test_suite_hpp

#include <cmath>
#include "eli/util/tolerance.hpp"

#include "eli/geom/surface/piecewise.hpp"
#include "eli/geom/intersect/intersect_surface_surface.hpp"

template<typename data__>
class intersect_surface_surface_test_suite : public Test::Suite
{
  private:
    typedef eli::geom::surface::piecewise<eli::geom::surface::bezier, data__, 3> piecewise_surface_type;
    typedef typename piecewise_surface_type::surface_type surface_type;
    typedef typename piecewise_surface_type::point_type point_type;
    typedef typename piecewise_surface_type::data_type data_type;
    typedef typename piecewise_surface_type::index_type index_type;

    // Model size: large enough that the product of the two surfaces' normals is far from one
    data_type L;

  protected:
    void AddTests(const double &)
    {
      TEST_ADD( intersect_surface_surface_test_suite<double>::crossing_test);
      TEST_ADD( intersect_surface_surface_test_suite<double>::edge_crossing_test);
      TEST_ADD( intersect_surface_surface_test_suite<double>::held_test);
    }

  public:
    intersect_surface_surface_test_suite()
    {
      L = 1000;
      AddTests(data__());
    }
    ~intersect_surface_surface_test_suite()
    {
    }

  private:
    // The curved floor z = 0.2 x^2 / L over 0 <= x, y <= L, with u along x and v along y
    void make_floor(piecewise_surface_type &pws)
    {
      surface_type s(3, 3);
      point_type cp;
      data_type zc[4] = { 0, 0, static_cast<data_type>(1)/3, 1 };
      for (index_type i=0; i<=3; ++i)
      {
        for (index_type j=0; j<=3; ++j)
        {
          cp << L*i/3, L*j/3, static_cast<data_type>(0.2)*L*zc[i];
          s.set_control_point(cp, i, j);
        }
      }
      pws.init_uv(1, 1);
      pws.set(s, 0, 0);
    }

    // The plane x = x0 + 0.3 z + sy y over 0 <= y <= L and -L/2 <= z <= L/2, with u along y and v
    // along z
    void make_wall(piecewise_surface_type &pws, const data_type &x0, const data_type &sy)
    {
      surface_type s(3, 3);
      point_type cp;
      for (index_type i=0; i<=3; ++i)
      {
        for (index_type j=0; j<=3; ++j)
        {
          data_type z = -L/2 + L*j/3;
          cp << x0 + static_cast<data_type>(0.3)*z + sy*L*i/3, L*i/3, z;
          s.set_control_point(cp, i, j);
        }
      }
      pws.init_uv(1, 1);
      pws.set(s, 0, 0);
    }

    // The wall crosses the floor where x = L/2 + 0.3 z(x)
    void crossing_test()
    {
      piecewise_surface_type s1, s2;
      make_floor(s1);
      make_wall(s2, L/2, 0);

      data_type u1, v1, u2, v2, d;
      point_type pt;
      pt << L*static_cast<data_type>(0.51), L*static_cast<data_type>(0.4), L*static_cast<data_type>(0.01);

      index_type ret = eli::geom::intersect::intersect(u1, v1, u2, v2, d, s1, s2, pt,
                                                       static_cast<data_type>(0.51), static_cast<data_type>(0.4),
                                                       static_cast<data_type>(0.4), static_cast<data_type>(0.51));
      TEST_ASSERT(ret==0);
      TEST_ASSERT(d<1e-9*L);
      data_type x = L*u1;
      TEST_ASSERT(std::abs(x-L/2-static_cast<data_type>(0.06)*x*x/L)<1e-9*L);
      // Held to the plane through pt across the curve, which runs along y
      TEST_ASSERT(std::abs(v1-static_cast<data_type>(0.4))<1e-12);
    }

    // The intersection leaves the floor through its u = 1 edge at y = L/2.  From a guess past
    // the edge at y = 0.7 L the solve stops against the edge, short of where the curve crosses
    // it.
    void edge_crossing_test()
    {
      piecewise_surface_type s1, s2;
      make_floor(s1);
      make_wall(s2, static_cast<data_type>(0.89)*L, static_cast<data_type>(0.1));

      data_type u1, v1, u2, v2, d;
      point_type pt;
      pt << L, L*static_cast<data_type>(0.7), L*static_cast<data_type>(0.2);

      index_type ret = eli::geom::intersect::intersect(u1, v1, u2, v2, d, s1, s2, pt,
                                                       static_cast<data_type>(1), static_cast<data_type>(0.7),
                                                       static_cast<data_type>(0.7), static_cast<data_type>(0.7));
      TEST_ASSERT(ret==0);
      TEST_ASSERT(d<1e-9*L);
      TEST_ASSERT(u1==1);
      TEST_ASSERT(std::abs(v1-static_cast<data_type>(0.5))<1e-12);
      TEST_ASSERT(std::abs(u2-static_cast<data_type>(0.5))<1e-12);
      TEST_ASSERT(std::abs(v2-static_cast<data_type>(0.7))<1e-12);
    }

    // The crossing of the intersection with the floor's line v = 0.3, which a free solve from the
    // same guess, held to the plane through pt, does not reach
    void held_test()
    {
      piecewise_surface_type s1, s2;
      make_floor(s1);
      make_wall(s2, L/2, 0);

      data_type u1, v1, u2, v2, d;
      point_type pt;
      pt << L*static_cast<data_type>(0.52), L*static_cast<data_type>(0.3), L*static_cast<data_type>(0.01);

      index_type ret = eli::geom::intersect::intersect(u1, v1, u2, v2, d, s1, s2, pt,
                                                       static_cast<data_type>(0.52), static_cast<data_type>(0.3),
                                                       static_cast<data_type>(0.33), static_cast<data_type>(0.51), 1);
      TEST_ASSERT(ret==0);
      TEST_ASSERT(d<1e-9*L);
      TEST_ASSERT(v1==static_cast<data_type>(0.3));
      TEST_ASSERT(std::abs(u2-static_cast<data_type>(0.3))<1e-12);
      data_type x = L*u1;
      TEST_ASSERT(std::abs(x-L/2-static_cast<data_type>(0.06)*x*x/L)<1e-9*L);
    }
};

#endif
