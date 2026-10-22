/*********************************************************************************
* Copyright (c) 2026 Rob McDonald <rob.a.mcdonald@gmail.com>
*
* All rights reserved. This program and the accompanying materials
* are made available under the terms of the Eclipse Public License v1.0
* which accompanies this distribution, and is available at
* http://www.eclipse.org/legal/epl-v10.html
********************************************************************************/

#ifndef piecewise_cubic_spline_skinning_surface_creator_test_suite_hpp
#define piecewise_cubic_spline_skinning_surface_creator_test_suite_hpp

#include <cmath>
#include <limits>
#include <vector>

#include <cpptest.h>

#include "eli/geom/surface/piecewise_cubic_spline_skinning_surface_creator.hpp"

template<typename data__>
class piecewise_cubic_spline_skinning_surface_creator_test_suite : public Test::Suite
{
  private:
    typedef eli::geom::surface::piecewise<eli::geom::surface::bezier, data__, 3> piecewise_surface_type;
    typedef eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 3> piecewise_curve_type;
    typedef typename piecewise_curve_type::curve_type curve_type;
    typedef typename piecewise_surface_type::point_type point_type;
    typedef typename piecewise_surface_type::data_type data_type;
    typedef typename piecewise_surface_type::index_type index_type;
    typedef typename piecewise_surface_type::tolerance_type tolerance_type;
    typedef typename eli::geom::surface::piecewise_cubic_spline_skinning_surface_creator<data__, 3, tolerance_type> creator_type;
    typedef typename creator_type::rib_data_type rib_data_type;

  protected:
    void AddTests(const float &)
    {
      TEST_ADD(piecewise_cubic_spline_skinning_surface_creator_test_suite<float>::accumulated_parameter_test);
    }
    void AddTests(const double &)
    {
      TEST_ADD(piecewise_cubic_spline_skinning_surface_creator_test_suite<double>::accumulated_parameter_test);
    }
    void AddTests(const long double &)
    {
      TEST_ADD(piecewise_cubic_spline_skinning_surface_creator_test_suite<long double>::accumulated_parameter_test);
    }

  public:
    piecewise_cubic_spline_skinning_surface_creator_test_suite()
    {
      AddTests(data__());
    }
    ~piecewise_cubic_spline_skinning_surface_creator_test_suite()
    {
    }

  private:
    // A straight rib along x at height y.
    rib_data_type make_rib(const data_type &y)
    {
      piecewise_curve_type pc;
      curve_type c(1);
      typename curve_type::control_point_type cp;

      cp << 0, y, 0;
      c.set_control_point(cp, 0);
      cp << 1, y, 0;
      c.set_control_point(cp, 1);

      pc.set_t0(0);
      pc.push_back(c, 1);

      rib_data_type rib;
      rib.set_f(pc);
      return rib;
    }

    // Ribs placed the way a propeller blade places them: three per interval between the
    // stations, every station a discontinuity, all normalized by the station span.  For
    // stations 0.1, 0.25 and 1 in double, the rib parameter differences summed again from
    // the start come to 0.9999999999999999, short of the last discontinuity at 1.
    void accumulated_parameter_test()
    {
      std::vector<data_type> station(3);
      station[0] = static_cast<data_type>(0.1);
      station[1] = static_cast<data_type>(0.25);
      station[2] = static_cast<data_type>(1);

      data_type rfirst(station[0]), rlast(station[2]);

      std::vector<data_type> param, tdisc;
      for (size_t k=0; k<station.size()-1; ++k)
      {
        data_type dt((station[k+1]-station[k])/3);
        param.push_back((station[k]-rfirst)/(rlast-rfirst));
        param.push_back((station[k]+dt-rfirst)/(rlast-rfirst));
        param.push_back((station[k]+2*dt-rfirst)/(rlast-rfirst));
      }
      param.push_back((station.back()-rfirst)/(rlast-rfirst));

      for (size_t k=0; k<station.size(); ++k)
      {
        tdisc.push_back((station[k]-rfirst)/(rlast-rfirst));
      }

      index_type nribs(static_cast<index_type>(param.size())), i;

      std::vector<rib_data_type> ribs(nribs);
      for (i=0; i<nribs; ++i)
      {
        ribs[i] = make_rib(10*param[i]*param[i]);
      }

      std::vector<index_type> maxd(nribs-1, 0);
      creator_type sc;
      TEST_ASSERT(sc.set_conditions(ribs, maxd, false));

      sc.set_u0(param[0]);
      for (i=0; i<nribs-1; ++i)
      {
        sc.set_segment_du(param[i+1]-param[i], i);
      }
      sc.set_tdisc(tdisc);

      piecewise_surface_type ps;
      TEST_ASSERT(sc.create(ps));
      TEST_ASSERT(ps.number_u_patches()==nribs-1);

      data_type umin, vmin, umax, vmax;
      ps.get_parameter_min(umin, vmin);
      ps.get_parameter_max(umax, vmax);

      // The surface passes through every rib.
      data_type tol(std::sqrt(std::numeric_limits<data_type>::epsilon()));
      data_type u(umin);
      for (i=0; i<nribs; ++i)
      {
        point_type p(ps.f(u, static_cast<data_type>(0.5)));
        TEST_ASSERT(std::abs(p.y()-10*param[i]*param[i])<tol);
        if (i<nribs-1)
        {
          u += param[i+1]-param[i];
          if (u>umax)
          {
            u = umax;
          }
        }
      }
    }
};

#endif
