/*********************************************************************************
* Copyright (c) 2026 Rob McDonald <rob.a.mcdonald@gmail.com>
*
* All rights reserved. This program and the accompanying materials
* are made available under the terms of the Eclipse Public License v1.0
* which accompanies this distribution, and is available at
* http://www.eclipse.org/legal/epl-v10.html
********************************************************************************/

#ifndef piecewise_uniform_skinning_surface_creator_test_suite_hpp
#define piecewise_uniform_skinning_surface_creator_test_suite_hpp

#include <cmath>
#include <limits>
#include <vector>

#include <cpptest.h>

#include "eli/geom/surface/piecewise_general_skinning_surface_creator.hpp"
#include "eli/geom/surface/piecewise_uniform_skinning_surface_creator.hpp"

// Parity tests: the uniform structure skinning creator must reproduce the general
// creator's surfaces to numerical precision for every rib configuration expressible
// through the current rib API.
template<typename data__>
class piecewise_uniform_skinning_surface_creator_test_suite : public Test::Suite
{
  private:
    typedef eli::geom::surface::piecewise<eli::geom::surface::bezier, data__, 3> piecewise_surface_type;
    typedef eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 3> piecewise_curve_type;
    typedef typename piecewise_curve_type::curve_type curve_type;
    typedef typename piecewise_surface_type::surface_type surface_type;
    typedef typename piecewise_surface_type::point_type point_type;
    typedef typename piecewise_surface_type::data_type data_type;
    typedef typename piecewise_surface_type::index_type index_type;
    typedef typename piecewise_surface_type::tolerance_type tolerance_type;
    typedef typename eli::geom::surface::connection_data<data__, 3, tolerance_type> rib_data_type;
    typedef typename eli::geom::surface::piecewise_general_skinning_surface_creator<data__, 3, tolerance_type> general_creator_type;
    typedef typename eli::geom::surface::piecewise_uniform_skinning_surface_creator<data__, 3, tolerance_type> uniform_creator_type;

    int rseed;

  protected:
    void AddTests(const float &)
    {
      TEST_ADD(piecewise_uniform_skinning_surface_creator_test_suite<float>::parity_test);
    }
    void AddTests(const double &)
    {
      TEST_ADD(piecewise_uniform_skinning_surface_creator_test_suite<double>::parity_test);
    }
    void AddTests(const long double &)
    {
      TEST_ADD(piecewise_uniform_skinning_surface_creator_test_suite<long double>::parity_test);
    }

  public:
    piecewise_uniform_skinning_surface_creator_test_suite() : rseed(7)
    {
      AddTests(data__());
    }
    ~piecewise_uniform_skinning_surface_creator_test_suite()
    {
    }

  private:
    data_type rnd()
    {
      rseed = ( 1103515245 * rseed + 12345 ) & 0x7fffffff;
      return static_cast<data_type>( rseed % 2000 ) / 1000 - 1;
    }

    piecewise_curve_type make_section(index_type nseg, const data_type &xoff, const data_type &scale)
    {
      piecewise_curve_type pc;
      pc.set_t0(0);
      for (index_type s=0; s<nseg; ++s)
      {
        curve_type c;
        c.resize(3);
        for (index_type i=0; i<=3; ++i)
        {
          typename curve_type::control_point_type cp;
          cp << xoff+rnd()/10, scale*(s+static_cast<data_type>(i)/3)+rnd()/10, scale*rnd();
          c.set_control_point(cp, i);
        }
        pc.push_back(c, static_cast<data_type>(1)/nseg);
      }
      return pc;
    }

    piecewise_curve_type make_asym_spine(index_type nseg)
    {
      piecewise_curve_type pc;
      pc.set_t0(0);
      for (index_type s=0; s<nseg; ++s)
      {
        curve_type c;
        c.resize(3);
        data_type side((s<nseg/2) ? 10 : static_cast<data_type>(-0.05));
        for (index_type i=0; i<=3; ++i)
        {
          typename curve_type::control_point_type cp;
          cp << side*(1+3*rnd()/10), side*rnd()/2, rnd()/5;
          c.set_control_point(cp, i);
        }
        pc.push_back(c, static_cast<data_type>(1)/nseg);
      }
      return pc;
    }

    bool build_and_compare(std::vector<rib_data_type> &ribs, std::vector<index_type> maxdeg = std::vector<index_type>())
    {
      if (maxdeg.empty())
      {
        maxdeg.assign(ribs.size()-1, 0);
      }

      general_creator_type gc;
      uniform_creator_type uc;
      piecewise_surface_type psg, psu;
      index_type nseg(static_cast<index_type>(ribs.size())-1), i;

      if (!gc.set_conditions(ribs, maxdeg, false))
        return false;
      if (!uc.set_conditions(ribs, maxdeg, false))
        return false;

      gc.set_u0(0);
      uc.set_u0(0);
      for (i=0; i<nseg; ++i)
      {
        gc.set_segment_du(1+static_cast<data_type>(i)/2, i);
        uc.set_segment_du(1+static_cast<data_type>(i)/2, i);
      }

      if (!gc.create(psg))
        return false;
      if (!uc.create(psu))
        return false;

      if ( (psg.number_u_patches()!=psu.number_u_patches()) ||
           (psg.number_v_patches()!=psu.number_v_patches()) )
        return false;

      data_type err(0), scale(0);
      for (index_type iu=0; iu<psg.number_u_patches(); ++iu)
      {
        for (index_type iv=0; iv<psg.number_v_patches(); ++iv)
        {
          const surface_type *sg(psg.get_patch(iu, iv)), *su(psu.get_patch(iu, iv));
          if ( (sg->degree_u()!=su->degree_u()) || (sg->degree_v()!=su->degree_v()) )
            return false;
          for (index_type ii=0; ii<=sg->degree_u(); ++ii)
          {
            for (index_type jj=0; jj<=sg->degree_v(); ++jj)
            {
              err=std::max(err, static_cast<data_type>((sg->get_control_point(ii, jj)-su->get_control_point(ii, jj)).norm()));
              scale=std::max(scale, static_cast<data_type>(sg->get_control_point(ii, jj).norm()));
            }
          }
        }
      }
      if (scale==0)
      {
        scale=1;
      }

      return err<=500*std::numeric_limits<data_type>::epsilon()*scale;
    }

    void parity_test()
    {
      // plain C0 ribs
      {
        std::vector<rib_data_type> ribs(3);
        for (index_type i=0; i<3; ++i)
        {
          ribs[i].set_f(make_section(4, static_cast<data_type>(i), 1));
        }
        TEST_ASSERT(build_and_compare(ribs));
      }

      // mixed continuity with side-varying interior derivative spines
      {
        std::vector<rib_data_type> ribs(4);
        for (index_type i=0; i<4; ++i)
        {
          ribs[i].set_f(make_section(4, static_cast<data_type>(i), 1+static_cast<data_type>(i)/5));
        }
        ribs[1].set_continuity(rib_data_type::C1);
        ribs[2].set_left_fp(make_asym_spine(4));
        ribs[2].set_right_fp(make_asym_spine(4));
        TEST_ASSERT(build_and_compare(ribs));
      }

      // end tangent and curvature spines
      {
        std::vector<rib_data_type> ribs(3);
        for (index_type i=0; i<3; ++i)
        {
          ribs[i].set_f(make_section(4, static_cast<data_type>(i), 1));
        }
        ribs[0].set_right_fp(make_asym_spine(4));
        ribs[0].set_right_fpp(make_asym_spine(4));
        ribs[2].set_left_fp(make_asym_spine(4));
        TEST_ASSERT(build_and_compare(ribs));
      }

      // mismatched rib segmentation (split/promote path) with C2 and interior spines
      {
        std::vector<rib_data_type> ribs(4);
        ribs[0].set_f(make_section(2, 0, 1));
        ribs[1].set_f(make_section(4, 1, static_cast<data_type>(1.2)));
        ribs[2].set_f(make_section(3, 2, static_cast<data_type>(0.8)));
        ribs[3].set_f(make_section(6, 3, 1));
        ribs[1].set_continuity(rib_data_type::C2);
        ribs[2].set_continuity(rib_data_type::C1);
        ribs[2].set_left_fp(make_asym_spine(3));
        ribs[2].set_right_fp(make_asym_spine(3));
        TEST_ASSERT(build_and_compare(ribs));
      }

      // one sided C1/C2 derivative mirroring and a bare C2 joint
      {
        std::vector<rib_data_type> ribs(5);
        for (index_type i=0; i<5; ++i)
        {
          ribs[i].set_f(make_section(3, static_cast<data_type>(i), 1));
        }
        ribs[1].set_continuity(rib_data_type::C1);
        ribs[1].set_left_fp(make_asym_spine(3));
        ribs[2].set_continuity(rib_data_type::C2);
        ribs[2].set_right_fp(make_asym_spine(3));
        ribs[2].set_right_fpp(make_asym_spine(3));
        ribs[3].set_continuity(rib_data_type::C2);
        TEST_ASSERT(build_and_compare(ribs));
      }

      // maximum degree caps interacting with continuity degree bumps
      {
        std::vector<rib_data_type> ribs(4);
        for (index_type i=0; i<4; ++i)
        {
          ribs[i].set_f(make_section(3, static_cast<data_type>(i), 1));
        }
        ribs[1].set_continuity(rib_data_type::C2);
        ribs[2].set_continuity(rib_data_type::C1);
        std::vector<index_type> maxdeg(3);
        maxdeg[0]=4;
        maxdeg[1]=5;
        maxdeg[2]=4;
        TEST_ASSERT(build_and_compare(ribs, maxdeg));
      }
    }
};

#endif
