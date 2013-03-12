/*********************************************************************************
* Copyright (c) 2013 David D. Marshall <ddmarsha@calpoly.edu>
*
* All rights reserved. This program and the accompanying materials
* are made available under the terms of the Eclipse Public License v1.0
* which accompanies this distribution, and is available at
* http://www.eclipse.org/legal/epl-v10.html
*
* Contributors:
*    David D. Marshall - initial code and implementation
********************************************************************************/

#ifndef piecewise_surface_test_suite_hpp
#define piecewise_surface_test_suite_hpp

#include "eli/code_eli.hpp"

#include "eli/util/tolerance.hpp"

#include "eli/constants/math.hpp"
#include "eli/geom/point/distance.hpp"
#include "eli/geom/surface/bezier.hpp"
// #include "eli/geom/surface/area.hpp"
#include "eli/geom/surface/curvature.hpp"
#include "eli/geom/surface/piecewise.hpp"

#include <cmath>    // std::pow, std::exp
#include <cassert>  // assert()

#include <typeinfo> // typeid
#include <string>   // std::string
#include <sstream>  // std::stringstream
#include <iomanip>  // std::setw
#include <limits>   // std::numeric_limits

template<typename data__>
class piecewise_surface_test_suite : public Test::Suite
{
  private:
    typedef eli::geom::surface::piecewise<eli::geom::surface::bezier, data__, 3> piecewise_surface_type;
    typedef typename piecewise_surface_type::surface_type surface_type;
    typedef typename piecewise_surface_type::point_type point_type;
    typedef typename piecewise_surface_type::data_type data_type;
    typedef typename piecewise_surface_type::index_type index_type;
    typedef typename piecewise_surface_type::tolerance_type tolerance_type;

    tolerance_type tol;

  protected:
    void AddTests(const float &)
    {
      // add the tests
      TEST_ADD(piecewise_surface_test_suite<float>::creation_test);
      TEST_ADD(piecewise_surface_test_suite<float>::bounding_box_test);
      TEST_ADD(piecewise_surface_test_suite<float>::reverse_test);
      TEST_ADD(piecewise_surface_test_suite<float>::replace_test);
      TEST_ADD(piecewise_surface_test_suite<float>::evaluation_test);
      TEST_ADD(piecewise_surface_test_suite<float>::split_test);
      TEST_ADD(piecewise_surface_test_suite<float>::area_test);
    }
    void AddTests(const double &)
    {
      // add the tests
      TEST_ADD(piecewise_surface_test_suite<double>::creation_test);
      TEST_ADD(piecewise_surface_test_suite<double>::bounding_box_test);
      TEST_ADD(piecewise_surface_test_suite<double>::reverse_test);
      TEST_ADD(piecewise_surface_test_suite<double>::replace_test);
      TEST_ADD(piecewise_surface_test_suite<double>::evaluation_test);
      TEST_ADD(piecewise_surface_test_suite<double>::split_test);
      TEST_ADD(piecewise_surface_test_suite<double>::area_test);
    }
    void AddTests(const long double &)
    {
      // add the tests
      TEST_ADD(piecewise_surface_test_suite<long double>::creation_test);
      TEST_ADD(piecewise_surface_test_suite<long double>::bounding_box_test);
      TEST_ADD(piecewise_surface_test_suite<long double>::reverse_test);
      TEST_ADD(piecewise_surface_test_suite<long double>::replace_test);
      TEST_ADD(piecewise_surface_test_suite<long double>::evaluation_test);
      TEST_ADD(piecewise_surface_test_suite<long double>::split_test);
      TEST_ADD(piecewise_surface_test_suite<long double>::area_test);
    }

#ifdef ELI_QD_FOUND
    void AddTests(const dd_real &)
    {
      // add the tests
      TEST_ADD(piecewise_surface_test_suite<dd_real>::creation_test);
      TEST_ADD(piecewise_surface_test_suite<dd_real>::bounding_box_test);
      TEST_ADD(piecewise_surface_test_suite<dd_real>::reverse_test);
      TEST_ADD(piecewise_surface_test_suite<dd_real>::replace_test);
      TEST_ADD(piecewise_surface_test_suite<dd_real>::evaluation_test);
      TEST_ADD(piecewise_surface_test_suite<dd_real>::split_test);
      TEST_ADD(piecewise_surface_test_suite<dd_real>::area_test);
    }

    void AddTests(const qd_real &)
    {
      // add the tests
      TEST_ADD(piecewise_surface_test_suite<qd_real>::creation_test);
      TEST_ADD(piecewise_surface_test_suite<qd_real>::bounding_box_test);
      TEST_ADD(piecewise_surface_test_suite<qd_real>::reverse_test);
      TEST_ADD(piecewise_surface_test_suite<qd_real>::replace_test);
      TEST_ADD(piecewise_surface_test_suite<qd_real>::evaluation_test);
      TEST_ADD(piecewise_surface_test_suite<qd_real>::split_test);
      TEST_ADD(piecewise_surface_test_suite<qd_real>::area_test);
    }
#endif

  public:
    piecewise_surface_test_suite()
    {
      AddTests(data__());
    }
    ~piecewise_surface_test_suite()
    {
    }

  private:
    void octave_print(int figno, /*const point_type pts[][], */const piecewise_surface_type &ps) const
    {
      index_type i, j, pp, qq, nup, nvp;
      data_type umin, vmin, umax, vmax;

      nup=ps.number_u_patches();
      nvp=ps.number_v_patches();
      ps.get_parameter_min(umin, vmin);
      ps.get_parameter_max(umax, vmax);

      std::cout << "figure(" << figno << ");" << std::endl;
      std::cout << "cp_x=[" << std::endl;
      for (pp=0; pp<nup; ++pp)
      {
        for (qq=0; qq<nvp; ++qq)
        {
          surface_type bez;
          ps.get(bez, pp, qq);
          for (i=0; i<=bez.degree_u(); ++i)
          {
            std::cout << bez.get_control_point(i, 0).x();
            for (j=1; j<bez.degree_v(); ++j)
            {
              std::cout << ", " << bez.get_control_point(i, j).x();
            }
            j=bez.degree_v();
            std::cout << ", " << bez.get_control_point(i, j).x();
            if (i<bez.degree_u())
              std::cout << "; ";
            else if (pp<nup-1)
              std::cout << "; ";
            else if (qq<nvp-1)
              std::cout << "; ";
          }
          std::cout << std::endl;
        }
      }
      std::cout << "];" << std::endl;

      std::cout << "cp_y=[";
      for (pp=0; pp<nup; ++pp)
      {
        for (qq=0; qq<nvp; ++qq)
        {
          surface_type bez;
          ps.get(bez, pp, qq);
          for (i=0; i<=bez.degree_u(); ++i)
          {
            std::cout << bez.get_control_point(i, 0).y();
            for (j=1; j<bez.degree_v(); ++j)
            {
              std::cout << ", " << bez.get_control_point(i, j).y();
            }
            j=bez.degree_v();
            std::cout << ", " << bez.get_control_point(i, j).y();
            if (i<bez.degree_u())
              std::cout << "; ";
            else if (pp<nup-1)
              std::cout << "; ";
            else if (qq<nvp-1)
              std::cout << "; ";
          }
          std::cout << std::endl;
        }
      }
      std::cout << "];" << std::endl;

      std::cout << "cp_z=[";
      for (pp=0; pp<nup; ++pp)
      {
        for (qq=0; qq<nvp; ++qq)
        {
          surface_type bez;
          ps.get(bez, pp, qq);
          for (i=0; i<=bez.degree_u(); ++i)
          {
            std::cout << bez.get_control_point(i, 0).z();
            for (j=1; j<bez.degree_v(); ++j)
            {
              std::cout << ", " << bez.get_control_point(i, j).z();
            }
            j=bez.degree_v();
            std::cout << ", " << bez.get_control_point(i, j).z();
            if (i<bez.degree_u())
              std::cout << "; ";
            else if (pp<nup-1)
              std::cout << "; ";
            else if (qq<nvp-1)
              std::cout << "; ";
          }
          std::cout << std::endl;
        }
      }
      std::cout << "];" << std::endl;

      // initialize the u & v parameters
      std::vector<data__> u(11), v(11);
      for (i=0; i<static_cast<index_type>(u.size()); ++i)
      {
        u[i]=umin+(umax-umin)*static_cast<data__>(i)/(u.size()-1);
      }
      for (j=0; j<static_cast<index_type>(v.size()); ++j)
      {
        v[j]=vmin+(vmax-vmin)*static_cast<data__>(j)/(v.size()-1);
      }

      // set the surface points
      std::cout << "surf_x=[";
      for (i=0; i<static_cast<index_type>(u.size()); ++i)
      {
        std::cout << ps.f(u[i], v[0]).x();
        for (j=1; j<static_cast<index_type>(v.size()-1); ++j)
        {
          std::cout << ", " << ps.f(u[i], v[j]).x();
        }
        j=static_cast<index_type>(v.size()-1);
        std::cout << ", " << ps.f(u[i], v[j]).x();
        if (i<static_cast<index_type>(u.size()-1))
          std::cout << "; " << std::endl;
      }
      std::cout << "];" << std::endl;

      std::cout << "surf_y=[";
      for (i=0; i<static_cast<index_type>(u.size()); ++i)
      {
        std::cout << ps.f(u[i], v[0]).y();
        for (j=1; j<static_cast<index_type>(v.size()-1); ++j)
        {
          std::cout << ", " << ps.f(u[i], v[j]).y();
        }
        j=static_cast<index_type>(v.size()-1);
        std::cout << ", " << ps.f(u[i], v[j]).y();
        if (i<static_cast<index_type>(u.size()-1))
          std::cout << "; " << std::endl;
      }
      std::cout << "];" << std::endl;

      std::cout << "surf_z=[";
      for (i=0; i<static_cast<index_type>(u.size()); ++i)
      {
        std::cout << ps.f(u[i], v[0]).z();
        for (j=1; j<static_cast<index_type>(v.size()-1); ++j)
        {
          std::cout << ", " << ps.f(u[i], v[j]).z();
        }
        j=static_cast<index_type>(v.size()-1);
        std::cout << ", " << ps.f(u[i], v[j]).z();
        if (i<static_cast<index_type>(u.size()-1))
          std::cout << "; " << std::endl;
      }
      std::cout << "];" << std::endl;

      std::cout << "setenv('GNUTERM', 'x11');" << std::endl;
      std::cout << "mesh(surf_x, surf_y, surf_z, zeros(size(surf_z)), 'EdgeColor', [0 0 0]);" << std::endl;
      std::cout << "hold on;" << std::endl;
      std::cout << "plot3(cp_x, cp_y, cp_z, 'ok', 'MarkerFaceColor', [0 0 0]);" << std::endl;
      std::cout << "hold off;" << std::endl;
    }

    void creation_test()
    {
      surface_type s, s1, s2, s3, s4, s5, s6, s_patches[6];
      piecewise_surface_type ps1, ps2;
      typename piecewise_surface_type::error_code err;

      // create 3x2 patches with unit spacing
      ps1.resize(3, 2);
      TEST_ASSERT(ps1.number_u_patches()==3);
      TEST_ASSERT(ps1.number_v_patches()==2);

      // create and set each surface
      index_type i, j, n(3), m(3);
      point_type pt[3+1][3+1], pt_out;

      // create surface with specified control points
      pt[0][0] << -15, 0,  15;
      pt[1][0] <<  -5, 5,  15;
      pt[2][0] <<   5, 5,  15;
      pt[3][0] <<  15, 0,  15;
      pt[0][1] << -15, 5,   5;
      pt[1][1] <<  -5, 5,   5;
      pt[2][1] <<   5, 5,   5;
      pt[3][1] <<  15, 5,   5;
      pt[0][2] << -15, 5,  -5;
      pt[1][2] <<  -5, 5,  -5;
      pt[2][2] <<   5, 5,  -5;
      pt[3][2] <<  15, 5,  -5;
      pt[0][3] << -15, 0, -15;
      pt[1][3] <<  -5, 5, -15;
      pt[2][3] <<   5, 5, -15;
      pt[3][3] <<  15, 0, -15;
      s.resize(n, m);
      for (i=0; i<=n; ++i)
      {
        for (j=0; j<=m; ++j)
        {
          s.set_control_point(pt[i][j], i, j);
        }
      }

      s.split_v(s1, s2, 0.5);  // this splits surface into lower and upper
      s1.split_u(s3, s4, 0.5); // this splits lower into first segment and last two
      err=ps1.replace(s3, 0, 0); s_patches[0]=s3;
      TEST_ASSERT(err==piecewise_surface_type::NO_ERROR);
      s2.split_u(s5, s6, 0.5); // this splits upper into first segment and last two
      err=ps1.replace(s5, 0, 1); s_patches[3]=s5;
      TEST_ASSERT(err==piecewise_surface_type::NO_ERROR);
      s4.split_u(s1, s2, 0.5); // this splits lower end into final two pieces
      err=ps1.replace(s1, 1, 0); s_patches[1]=s1;
      TEST_ASSERT(err==piecewise_surface_type::NO_ERROR);
      err=ps1.replace(s2, 2, 0); s_patches[2]=s2;
      TEST_ASSERT(err==piecewise_surface_type::NO_ERROR);
      s6.split_u(s1, s2, 0.5); // this splits the upper end into final two pieces
      err=ps1.replace(s1, 1, 1); s_patches[4]=s1;
      TEST_ASSERT(err==piecewise_surface_type::NO_ERROR);
      err=ps1.replace(s2, 2, 1); s_patches[5]=s2;
      TEST_ASSERT(err==piecewise_surface_type::NO_ERROR);

      // print out
//       if (typeid(data_type)==typeid(double))
//         octave_print(2, ps1);

      // test getting parameter max
      data_type umax, vmax;
      ps1.get_parameter_max(umax, vmax);
      TEST_ASSERT(tol.approximately_equal(umax, 3));
      TEST_ASSERT(tol.approximately_equal(vmax, 2));

      // test copy ctr
      piecewise_surface_type ps1c(ps1);
      TEST_ASSERT(ps1==ps1c);

      // test changing u0 and v0
      ps1c.set_u0(-1.5);
      TEST_ASSERT(ps1c.get_u0()==-1.5);
      ps1c.set_v0(-1.5);
      TEST_ASSERT(ps1c.get_v0()==-1.5);

      // test assignment operator
      ps2=ps1;
      TEST_ASSERT(ps2==ps1);
    }

    void bounding_box_test()
    {
      surface_type s, s1, s2, s3, s4, s5, s6, s_patches[6];
      piecewise_surface_type ps1, ps2;
      typename piecewise_surface_type::error_code err;

      // create 3x2 patches with unit spacing
      ps1.resize(3, 2);
      TEST_ASSERT(ps1.number_u_patches()==3);
      TEST_ASSERT(ps1.number_v_patches()==2);

      // create and set each surface
      index_type i, j, n(3), m(3);
      point_type pt[3+1][3+1], pt_out;

      // create surface with specified control points
      pt[0][0] << -15, 0,  15;
      pt[1][0] <<  -5, 5,  15;
      pt[2][0] <<   5, 5,  15;
      pt[3][0] <<  15, 0,  15;
      pt[0][1] << -15, 5,   5;
      pt[1][1] <<  -5, 5,   5;
      pt[2][1] <<   5, 5,   5;
      pt[3][1] <<  15, 5,   5;
      pt[0][2] << -15, 5,  -5;
      pt[1][2] <<  -5, 5,  -5;
      pt[2][2] <<   5, 5,  -5;
      pt[3][2] <<  15, 5,  -5;
      pt[0][3] << -15, 0, -15;
      pt[1][3] <<  -5, 5, -15;
      pt[2][3] <<   5, 5, -15;
      pt[3][3] <<  15, 0, -15;
      s.resize(n, m);
      for (i=0; i<=n; ++i)
      {
        for (j=0; j<=m; ++j)
        {
          s.set_control_point(pt[i][j], i, j);
        }
      }

      s.split_v(s1, s2, 0.5);  // this splits surface into lower and upper
      s1.split_u(s3, s4, 0.5); // this splits lower into first segment and last two
      err=ps1.replace(s3, 0, 0); s_patches[0]=s3;
      TEST_ASSERT(err==piecewise_surface_type::NO_ERROR);
      s2.split_u(s5, s6, 0.5); // this splits upper into first segment and last two
      err=ps1.replace(s5, 0, 1); s_patches[3]=s5;
      TEST_ASSERT(err==piecewise_surface_type::NO_ERROR);
      s4.split_u(s1, s2, 0.5); // this splits lower end into final two pieces
      err=ps1.replace(s1, 1, 0); s_patches[1]=s1;
      TEST_ASSERT(err==piecewise_surface_type::NO_ERROR);
      err=ps1.replace(s2, 2, 0); s_patches[2]=s2;
      TEST_ASSERT(err==piecewise_surface_type::NO_ERROR);
      s6.split_u(s1, s2, 0.5); // this splits the upper end into final two pieces
      err=ps1.replace(s1, 1, 1); s_patches[4]=s1;
      TEST_ASSERT(err==piecewise_surface_type::NO_ERROR);
      err=ps1.replace(s2, 2, 1); s_patches[5]=s2;
      TEST_ASSERT(err==piecewise_surface_type::NO_ERROR);

      // test the bounding box
      point_type pmin, pmax, pmin_ref, pmax_ref;

      ps1.get_bounding_box(pmin, pmax);
      pmin_ref << -15, 0, -15;
      pmax_ref << 15, 4.6875, 15;
      TEST_ASSERT(pmin==pmin_ref);
      TEST_ASSERT(pmax==pmax_ref);
    }

    void reverse_test()
    {
#if 0
      piecewise_surface_type c1, c2;
      curve_type bc[3], bc1_out, bc2_out;
      data_type dt[3], dt1_out, dt2_out;
      index_type i;
      typename curve_type::control_point_type cntrl1_in[4], cntrl2_in[5], cntrl3_in[3];
      typename piecewise_surface_type::error_code err;

      // create bezier curves
      cntrl1_in[0] << 2.0, 2.0, 0.0;
      cntrl1_in[1] << 1.0, 1.5, 0.0;
      cntrl1_in[2] << 3.5, 0.0, 0.0;
      cntrl1_in[3] << 4.0, 1.0, 0.0;
      dt[0]=0.5;
      bc[0].resize(3);
      for (i=0; i<4; ++i)
      {
        bc[0].set_control_point(cntrl1_in[i], i);
      }
      cntrl2_in[0] << 4.0, 1.0, 0.0;
      cntrl2_in[1] << 5.0, 2.5, 0.0;
      cntrl2_in[2] << 5.5, 1.0, 0.0;
      cntrl2_in[3] << 6.0, 0.0, 0.0;
      cntrl2_in[4] << 6.5,-0.5, 0.0;
      dt[1]=2.0;
      bc[1].resize(4);
      for (i=0; i<5; ++i)
      {
        bc[1].set_control_point(cntrl2_in[i], i);
      }
      cntrl3_in[0] << 6.5,-0.5, 0.0;
      cntrl3_in[1] << 6.0,-1.0, 0.0;
      cntrl3_in[2] << 5.5,-2.0, 0.0;
      dt[2]=1.5;
      bc[2].resize(2);
      for (i=0; i<3; ++i)
      {
        bc[2].set_control_point(cntrl3_in[i], i);
      }

      // initialize by passing iterators to curve collection
      err=c1.set(bc, bc+3, dt);
      TEST_ASSERT(err==piecewise_surface_type::NO_ERROR);

      c1=c2;
      c2.reverse();
      for (i=0; i<c1.number_segments(); ++i)
      {
        c1.get(bc1_out, dt1_out, i);
        c2.get(bc2_out, dt2_out, c1.number_segments()-i-1);
        for (index_type ii=0; ii<=bc1_out.degree(); ++ii)
        {
          TEST_ASSERT(bc1_out.get_control_point(ii)==bc2_out.get_control_point(bc1_out.degree()-ii));
        }
      }
#endif
    }

    void replace_test()
    {
#if 0
      piecewise_surface_type c1, c1c, c2, c3;
      curve_type bc[3], bc2, bc_out;
      data_type dt[3], dt2, dt_out;
      index_type i;
      typename curve_type::control_point_type cntrl1_in[4], cntrl2_in[5], cntrl3_in[3], cntrl2a_in[2];
      typename piecewise_surface_type::error_code err;

      // create bezier curves
      cntrl1_in[0] << 2.0, 2.0, 0.0;
      cntrl1_in[1] << 1.0, 1.5, 0.0;
      cntrl1_in[2] << 3.5, 0.0, 0.0;
      cntrl1_in[3] << 4.0, 1.0, 0.0;
      dt[0]=0.5;
      bc[0].resize(3);
      for (i=0; i<4; ++i)
      {
        bc[0].set_control_point(cntrl1_in[i], i);
      }
      cntrl2_in[0] << 4.0, 1.0, 0.0;
      cntrl2_in[1] << 5.0, 2.5, 0.0;
      cntrl2_in[2] << 5.5, 1.0, 0.0;
      cntrl2_in[3] << 6.0, 0.0, 0.0;
      cntrl2_in[4] << 6.5,-0.5, 0.0;
      dt[1]=2.0;
      bc[1].resize(4);
      for (i=0; i<5; ++i)
      {
        bc[1].set_control_point(cntrl2_in[i], i);
      }
      cntrl3_in[0] << 6.5,-0.5, 0.0;
      cntrl3_in[1] << 6.0,-1.0, 0.0;
      cntrl3_in[2] << 5.5,-2.0, 0.0;
      dt[2]=1.5;
      bc[2].resize(2);
      for (i=0; i<3; ++i)
      {
        bc[2].set_control_point(cntrl3_in[i], i);
      }

      // initialize by passing iterators to curve collection
      err=c1.set(bc, bc+3, dt);
      TEST_ASSERT(err==piecewise_surface_type::NO_ERROR);

      // replace segment correctly
      cntrl2a_in[0] << 4.0, 1.0, 0.0;
      cntrl2a_in[1] << 6.5,-0.5, 0.0;
      bc2.resize(1);
      for (i=0; i<2; ++i)
      {
        bc2.set_control_point(cntrl2a_in[i], i);
      }
      err=c1.replace(bc2, 1);
      TEST_ASSERT(err==piecewise_surface_type::NO_ERROR);
      err=c1.get(bc_out, dt_out, 1);
      TEST_ASSERT(err==piecewise_surface_type::NO_ERROR);
      TEST_ASSERT(bc2==bc_out);
      TEST_ASSERT(dt[1]==dt_out);

      // replace segment with segment that doesn't connect to neighbor(s)
      cntrl2a_in[0] << 4.0, 1.0, 0.0;
      cntrl2a_in[1] << 5.5,-0.5, 0.0;
      bc2.resize(1);
      for (i=0; i<2; ++i)
      {
        bc2.set_control_point(cntrl2a_in[i], i);
      }
      err=c1.replace(bc2, 1);
      TEST_ASSERT(err==piecewise_surface_type::SEGMENT_NOT_CONNECTED);
      cntrl2a_in[0] << 4.0, 2.0, 0.0;
      cntrl2a_in[1] << 6.5,-0.5, 0.0;
      bc2.resize(1);
      for (i=0; i<2; ++i)
      {
        bc2.set_control_point(cntrl2a_in[i], i);
      }
      err=c1.replace(bc2, 1);
      TEST_ASSERT(err==piecewise_surface_type::SEGMENT_NOT_CONNECTED);

      // replace one segment with several segments
      cntrl1_in[0] << 4.0, 1.0, 0.0;
      cntrl1_in[1] << 3.0, 1.5, 0.0;
      cntrl1_in[2] << 3.5, 2.0, 0.0;
      cntrl1_in[3] << 5.5, 2.5, 0.0;
      dt[0]=0.25;
      bc[0].resize(3);
      for (i=0; i<4; ++i)
      {
        bc[0].set_control_point(cntrl1_in[i], i);
      }
      cntrl2_in[0] << 5.5, 2.5, 0.0;
      cntrl2_in[1] << 5.0, 2.5, 0.0;
      cntrl2_in[2] << 5.5, 2.0, 0.0;
      cntrl2_in[3] << 5.0, 2.0, 0.0;
      cntrl2_in[4] << 5.0, 0.5, 0.0;
      dt[1]=1.0;
      bc[1].resize(4);
      for (i=0; i<5; ++i)
      {
        bc[1].set_control_point(cntrl2_in[i], i);
      }
      cntrl3_in[0] << 5.0, 0.5, 0.0;
      cntrl3_in[1] << 5.5, 0.0, 0.0;
      cntrl3_in[2] << 6.5,-0.5, 0.0;
      dt[2]=3.0;
      bc[2].resize(2);
      for (i=0; i<3; ++i)
      {
        bc[2].set_control_point(cntrl3_in[i], i);
      }
      err=c2.set(bc, bc+3, dt);
      TEST_ASSERT(err==piecewise_surface_type::NO_ERROR);
      c1c=c1;
      err=c1.replace(c2, 1);
      TEST_ASSERT(err==piecewise_surface_type::NO_ERROR);
      TEST_ASSERT(c1.number_segments()==5);
      err=c1.get(bc_out, dt_out, 1);
      TEST_ASSERT(err==piecewise_surface_type::NO_ERROR);
      TEST_ASSERT(bc[0]==bc_out);
      TEST_ASSERT(dt[0]==dt_out);
      err=c1.get(bc_out, dt_out, 2);
      TEST_ASSERT(err==piecewise_surface_type::NO_ERROR);
      TEST_ASSERT(bc[1]==bc_out);
      TEST_ASSERT(dt[1]==dt_out);
      err=c1.get(bc_out, dt_out, 3);
      TEST_ASSERT(err==piecewise_surface_type::NO_ERROR);
      TEST_ASSERT(bc[2]==bc_out);
      TEST_ASSERT(dt[2]==dt_out);
      c3=c1;

      // replace segment with piecewise that doesn't connect to neighbor(s)
      c1=c1c;
      cntrl1_in[0] << 5.0, 1.0, 0.0;
      cntrl1_in[1] << 3.0, 1.5, 0.0;
      cntrl1_in[2] << 3.5, 2.0, 0.0;
      cntrl1_in[3] << 5.5, 2.5, 0.0;
      dt[0]=0.25;
      bc[0].resize(3);
      for (i=0; i<4; ++i)
      {
        bc[0].set_control_point(cntrl1_in[i], i);
      }
      cntrl2_in[0] << 5.5, 2.5, 0.0;
      cntrl2_in[1] << 5.0, 2.5, 0.0;
      cntrl2_in[2] << 5.5, 2.0, 0.0;
      cntrl2_in[3] << 5.0, 2.0, 0.0;
      cntrl2_in[4] << 5.0, 0.5, 0.0;
      dt[1]=1.0;
      bc[1].resize(4);
      for (i=0; i<5; ++i)
      {
        bc[1].set_control_point(cntrl2_in[i], i);
      }
      cntrl3_in[0] << 5.0, 0.5, 0.0;
      cntrl3_in[1] << 5.5, 0.0, 0.0;
      cntrl3_in[2] << 6.5,-0.5, 0.0;
      dt[2]=3.0;
      bc[2].resize(2);
      for (i=0; i<3; ++i)
      {
        bc[2].set_control_point(cntrl3_in[i], i);
      }
      err=c2.set(bc, bc+3, dt);
      TEST_ASSERT(err==piecewise_surface_type::NO_ERROR);
      err=c1.replace(c2, 1);
      TEST_ASSERT(err==piecewise_surface_type::SEGMENT_NOT_CONNECTED);
      cntrl1_in[0] << 4.0, 1.0, 0.0;
      cntrl1_in[1] << 3.0, 1.5, 0.0;
      cntrl1_in[2] << 3.5, 2.0, 0.0;
      cntrl1_in[3] << 5.5, 2.5, 0.0;
      dt[0]=0.25;
      bc[0].resize(3);
      for (i=0; i<4; ++i)
      {
        bc[0].set_control_point(cntrl1_in[i], i);
      }
      cntrl2_in[0] << 5.5, 2.5, 0.0;
      cntrl2_in[1] << 5.0, 2.5, 0.0;
      cntrl2_in[2] << 5.5, 2.0, 0.0;
      cntrl2_in[3] << 5.0, 2.0, 0.0;
      cntrl2_in[4] << 5.0, 0.5, 0.0;
      dt[1]=1.0;
      bc[1].resize(4);
      for (i=0; i<5; ++i)
      {
        bc[1].set_control_point(cntrl2_in[i], i);
      }
      cntrl3_in[0] << 5.0, 0.5, 0.0;
      cntrl3_in[1] << 5.5, 0.0, 0.0;
      cntrl3_in[2] << 7.5,-0.5, 0.0;
      dt[2]=3.0;
      bc[2].resize(2);
      for (i=0; i<3; ++i)
      {
        bc[2].set_control_point(cntrl3_in[i], i);
      }
      err=c2.set(bc, bc+3, dt);
      TEST_ASSERT(err==piecewise_surface_type::NO_ERROR);
      err=c1.replace(c2, 1);
      TEST_ASSERT(err==piecewise_surface_type::SEGMENT_NOT_CONNECTED);

      // replace several segments with one segment
      c1=c3;
      c1c=c1;
      cntrl2a_in[0] << 4.0, 1.0, 0.0;
      cntrl2a_in[1] << 5.0, 0.5, 0.0;
      dt2=1.25;
      bc2.resize(1);
      for (i=0; i<2; ++i)
      {
        bc2.set_control_point(cntrl2a_in[i], i);
      }
      err=c1.replace(bc2, dt2, 1, 3);
      TEST_ASSERT(err==piecewise_surface_type::NO_ERROR);
      TEST_ASSERT(c1.number_segments()==4);
      err=c1.get(bc_out, dt_out, 1);
      TEST_ASSERT(err==piecewise_surface_type::NO_ERROR);
      TEST_ASSERT(bc2==bc_out);
      TEST_ASSERT(dt2==dt_out);
      c3=c1;

      // replace segments with segment that doesn't connect to neighbor(s)
      c1=c1c;
      cntrl2a_in[0] << 4.5, 1.0, 0.0;
      cntrl2a_in[1] << 5.0, 0.5, 0.0;
      dt2=1.25;
      bc2.resize(1);
      for (i=0; i<2; ++i)
      {
        bc2.set_control_point(cntrl2a_in[i], i);
      }
      err=c1.replace(bc2, dt2, 1, 3);
      TEST_ASSERT(err==piecewise_surface_type::SEGMENT_NOT_CONNECTED);
      cntrl2a_in[0] << 4.0, 1.0, 0.0;
      cntrl2a_in[1] << 5.5, 0.5, 0.0;
      dt2=1.25;
      bc2.resize(1);
      for (i=0; i<2; ++i)
      {
        bc2.set_control_point(cntrl2a_in[i], i);
      }
      err=c1.replace(bc2, dt2, 1, 3);
      TEST_ASSERT(err==piecewise_surface_type::SEGMENT_NOT_CONNECTED);

      // replace several segments with several segments
      c1=c3;
      c1c=c1;
      cntrl1_in[0] << 4.0, 1.0, 0.0;
      cntrl1_in[1] << 3.0, 1.5, 0.0;
      cntrl1_in[2] << 3.5, 2.0, 0.0;
      cntrl1_in[3] << 5.5, 2.5, 0.0;
      dt[0]=0.25;
      bc[0].resize(3);
      for (i=0; i<4; ++i)
      {
        bc[0].set_control_point(cntrl1_in[i], i);
      }
      cntrl2_in[0] << 5.5, 2.5, 0.0;
      cntrl2_in[1] << 5.0, 2.5, 0.0;
      cntrl2_in[2] << 4.5, 2.0, 0.0;
      cntrl2_in[3] << 4.0, 2.0, 0.0;
      cntrl2_in[4] << 3.5, 0.5, 0.0;
      dt[1]=1.0;
      bc[1].resize(4);
      for (i=0; i<5; ++i)
      {
        bc[1].set_control_point(cntrl2_in[i], i);
      }
      cntrl3_in[0] << 3.5, 0.5, 0.0;
      cntrl3_in[1] << 4.5, 0.0, 0.0;
      cntrl3_in[2] << 6.5,-0.5, 0.0;
      dt[2]=3.0;
      bc[2].resize(2);
      for (i=0; i<3; ++i)
      {
        bc[2].set_control_point(cntrl3_in[i], i);
      }
      err=c2.set(bc, bc+3, dt);
      TEST_ASSERT(err==piecewise_surface_type::NO_ERROR);
      err=c1.replace(c2, 1, 3);
      TEST_ASSERT(err==piecewise_surface_type::NO_ERROR);
      TEST_ASSERT(c1.number_segments()==5);
      err=c1.get(bc_out, dt_out, 1);
      TEST_ASSERT(err==piecewise_surface_type::NO_ERROR);
      TEST_ASSERT(bc[0]==bc_out);
      TEST_ASSERT(dt[0]==dt_out);
      err=c1.get(bc_out, dt_out, 2);
      TEST_ASSERT(err==piecewise_surface_type::NO_ERROR);
      TEST_ASSERT(bc[1]==bc_out);
      TEST_ASSERT(dt[1]==dt_out);
      err=c1.get(bc_out, dt_out, 3);
      TEST_ASSERT(err==piecewise_surface_type::NO_ERROR);
      TEST_ASSERT(bc[2]==bc_out);
      TEST_ASSERT(dt[2]==dt_out);

      // replace segments with piecewise that doesn't connect to neighbor(s)
      c1=c1c;
      cntrl1_in[0] << 3.0, 1.0, 0.0;
      cntrl1_in[1] << 3.0, 1.5, 0.0;
      cntrl1_in[2] << 3.5, 2.0, 0.0;
      cntrl1_in[3] << 5.5, 2.5, 0.0;
      dt[0]=1.25;
      bc[0].resize(3);
      for (i=0; i<4; ++i)
      {
        bc[0].set_control_point(cntrl1_in[i], i);
      }
      cntrl2_in[0] << 5.5, 2.5, 0.0;
      cntrl2_in[1] << 5.0, 2.5, 0.0;
      cntrl2_in[2] << 4.5, 2.0, 0.0;
      cntrl2_in[3] << 4.0, 2.0, 0.0;
      cntrl2_in[4] << 3.5, 0.5, 0.0;
      dt[1]=1.5;
      bc[1].resize(4);
      for (i=0; i<5; ++i)
      {
        bc[1].set_control_point(cntrl2_in[i], i);
      }
      cntrl3_in[0] << 3.5, 0.5, 0.0;
      cntrl3_in[1] << 4.5, 0.0, 0.0;
      cntrl3_in[2] << 6.5,-0.5, 0.0;
      dt[2]=3.5;
      bc[2].resize(2);
      for (i=0; i<3; ++i)
      {
        bc[2].set_control_point(cntrl3_in[i], i);
      }
      err=c2.set(bc, bc+3, dt);
      TEST_ASSERT(err==piecewise_surface_type::NO_ERROR);
      err=c1.replace(c2, 1, 3);
      TEST_ASSERT(err==piecewise_surface_type::SEGMENT_NOT_CONNECTED);
      cntrl1_in[0] << 4.0, 1.0, 0.0;
      cntrl1_in[1] << 3.0, 1.5, 0.0;
      cntrl1_in[2] << 3.5, 2.0, 0.0;
      cntrl1_in[3] << 5.5, 2.5, 0.0;
      dt[0]=1.25;
      bc[0].resize(3);
      for (i=0; i<4; ++i)
      {
        bc[0].set_control_point(cntrl1_in[i], i);
      }
      cntrl2_in[0] << 5.5, 2.5, 0.0;
      cntrl2_in[1] << 5.0, 2.5, 0.0;
      cntrl2_in[2] << 4.5, 2.0, 0.0;
      cntrl2_in[3] << 4.0, 2.0, 0.0;
      cntrl2_in[4] << 3.5, 0.5, 0.0;
      dt[1]=1.5;
      bc[1].resize(4);
      for (i=0; i<5; ++i)
      {
        bc[1].set_control_point(cntrl2_in[i], i);
      }
      cntrl3_in[0] << 3.5, 0.5, 0.0;
      cntrl3_in[1] << 4.5, 0.0, 0.0;
      cntrl3_in[2] << 5.5,-0.5, 0.0;
      dt[2]=3.5;
      bc[2].resize(2);
      for (i=0; i<3; ++i)
      {
        bc[2].set_control_point(cntrl3_in[i], i);
      }
      err=c2.set(bc, bc+3, dt);
      TEST_ASSERT(err==piecewise_surface_type::NO_ERROR);
      err=c1.replace(c2, 1, 3);
      TEST_ASSERT(err==piecewise_surface_type::SEGMENT_NOT_CONNECTED);
#endif
    }

    void evaluation_test()
    {
#if 0
      data_type eps(std::numeric_limits<data__>::epsilon());
#ifdef ELI_QD_FOUND
      if ( (typeid(data_type)==typeid(dd_real)) || (typeid(data_type)==typeid(qd_real)) )
        eps=std::numeric_limits<double>::epsilon();
#endif

      // test two curves with delta t=1
      {
        piecewise_surface_type pwc;
        typename curve_type::control_point_type cntrl_in[4];
        curve_type bc1, bc1l, bc1r;
        point_type eval_out, eval_ref;
        data_type t, tl, tr, ts;
        tl = static_cast<data__>(0.3);
        tr = static_cast<data__>(0.87);
        ts = static_cast<data__>(0.586);

        // set control points
        cntrl_in[0] << 0, 0, 0;
        cntrl_in[1] << 0, 2, 0;
        cntrl_in[2] << 8, 2, 0;
        cntrl_in[3] << 4, 0, 0;
        bc1.resize(3);
        for (index_type i=0; i<4; ++i)
        {
          bc1.set_control_point(cntrl_in[i], i);
        }

        // split curve and create piecewise
        bc1.split(bc1l, bc1r, ts);
        pwc.push_back(bc1l);
        pwc.push_back(bc1r);
        TEST_ASSERT(pwc.number_segments()==2);

        // check the left curve
        t=tl*ts;
        eval_out=pwc.f(tl);
        eval_ref=bc1.f(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<1e2*eps);
        eval_out=pwc.fp(tl);
        eval_ref=bc1.fp(t)*ts;
        TEST_ASSERT((eval_out-eval_ref).norm()<1e2*eps);
        eval_out=pwc.fpp(tl);
        eval_ref=bc1.fpp(t)*ts*ts;
        TEST_ASSERT((eval_out-eval_ref).norm()<1e2*eps);
        eval_out=pwc.fppp(tl);
        eval_ref=bc1.fppp(t)*ts*ts*ts;
        TEST_ASSERT((eval_out-eval_ref).norm()<1e2*eps);

        // check the right curve
        t=ts+tr*(1-ts);
        eval_out=pwc.f(1+tr);
        eval_ref=bc1.f(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<1e2*eps);
        eval_out=pwc.fp(1+tr);
        eval_ref=bc1.fp(t)*(1-ts);
        TEST_ASSERT((eval_out-eval_ref).norm()<1e2*eps);
        eval_out=pwc.fpp(1+tr);
        eval_ref=bc1.fpp(t)*(1-ts)*(1-ts);
        TEST_ASSERT((eval_out-eval_ref).norm()<1e2*eps);
        eval_out=pwc.fppp(1+tr);
        eval_ref=bc1.fppp(t)*(1-ts)*(1-ts)*(1-ts);
        TEST_ASSERT((eval_out-eval_ref).norm()<1.21e2*eps);
      }

      // test two curves with delta t!=1
      {
        piecewise_surface_type pwc;
        typename curve_type::control_point_type cntrl_in[4];
        curve_type bc1, bc1l, bc1r;
        point_type eval_out, eval_ref;
        data_type tl, tr, ts;
        tl = static_cast<data__>(0.3);
        tr = static_cast<data__>(0.87);
        ts = static_cast<data__>(0.586);

        // set control points
        cntrl_in[0] << 0, 0, 0;
        cntrl_in[1] << 0, 2, 0;
        cntrl_in[2] << 8, 2, 0;
        cntrl_in[3] << 4, 0, 0;
        bc1.resize(3);
        for (index_type i=0; i<4; ++i)
        {
          bc1.set_control_point(cntrl_in[i], i);
        }

        // split curve and create piecewise
        bc1.split(bc1l, bc1r, ts);
        pwc.push_back(bc1l, ts);
        pwc.push_back(bc1r, 1-ts);
        TEST_ASSERT(pwc.number_segments()==2);

        // check the left curve
        eval_out=pwc.f(tl);
        eval_ref=bc1.f(tl);
        TEST_ASSERT((eval_out-eval_ref).norm()<1e2*eps);
        eval_out=pwc.fp(tl);
        eval_ref=bc1.fp(tl);
        TEST_ASSERT((eval_out-eval_ref).norm()<1e2*eps);
        eval_out=pwc.fpp(tl);
        eval_ref=bc1.fpp(tl);
        TEST_ASSERT((eval_out-eval_ref).norm()<1e2*eps);
        eval_out=pwc.fppp(tl);
        eval_ref=bc1.fppp(tl);
        TEST_ASSERT((eval_out-eval_ref).norm()<1e2*eps);

        // check the right curve
        eval_out=pwc.f(tr);
        eval_ref=bc1.f(tr);
        TEST_ASSERT((eval_out-eval_ref).norm()<1e2*eps);
        eval_out=pwc.fp(tr);
        eval_ref=bc1.fp(tr);
        TEST_ASSERT((eval_out-eval_ref).norm()<1e2*eps);
        eval_out=pwc.fpp(tr);
        eval_ref=bc1.fpp(tr);
        TEST_ASSERT((eval_out-eval_ref).norm()<2.6e2*eps);
        eval_out=pwc.fppp(tr);
        eval_ref=bc1.fppp(tr);
        TEST_ASSERT((eval_out-eval_ref).norm()<1.8e3*eps);
      }
#endif
    }

    void split_test()
    {
#if 0
      data_type eps(std::numeric_limits<data__>::epsilon());
#ifdef ELI_QD_FOUND
      if ( (typeid(data_type)==typeid(dd_real)) || (typeid(data_type)==typeid(qd_real)) )
        eps=std::numeric_limits<double>::epsilon();
#endif
        piecewise_surface_type pwc0, pwc1;
        typename curve_type::control_point_type cntrl_in[4];
        typename piecewise_surface_type::error_code err;
        curve_type bc;
        point_type eval_out, eval_ref;
        data_type ts, t;
        ts=static_cast<data__>(1.56);

        // build piecewise curve
        cntrl_in[0] << 0, 0, 0;
        cntrl_in[1] << 0, 2, 0;
        cntrl_in[2] << 8, 2, 0;
        cntrl_in[3] << 4, 0, 0;
        bc.resize(3);
        for (index_type i=0; i<4; ++i)
        {
          bc.set_control_point(cntrl_in[i], i);
        }
        err=pwc0.push_back(bc);
        TEST_ASSERT(err==piecewise_surface_type::NO_ERROR);
        cntrl_in[0] << 4,  0,   0;
        cntrl_in[1] << 3, -0.5, 0;
        cntrl_in[2] << 2, -1,   0;
        cntrl_in[3] << 1, -1,   0;
        bc.resize(3);
        for (index_type i=0; i<4; ++i)
        {
          bc.set_control_point(cntrl_in[i], i);
        }
        err=pwc0.push_back(bc);
        TEST_ASSERT(err==piecewise_surface_type::NO_ERROR);
        TEST_ASSERT(pwc0.number_segments()==2);

        // split curve and create piecewise
        pwc1=pwc0;
        err=pwc0.split(ts);
        TEST_ASSERT(err==piecewise_surface_type::NO_ERROR);
        TEST_ASSERT(pwc0.number_segments()==pwc1.number_segments()+1);

        t=0;
        eval_out=pwc0.f(t);
        eval_ref=pwc1.f(t);
        TEST_ASSERT(eval_out==eval_ref);
        eval_out=pwc0.fp(t);
        eval_ref=pwc1.fp(t);
        TEST_ASSERT(eval_out==eval_ref);
        eval_out=pwc0.fpp(t);
        eval_ref=pwc1.fpp(t);
        TEST_ASSERT(eval_out==eval_ref);
        eval_out=pwc0.fppp(t);
        eval_ref=pwc1.fppp(t);
        TEST_ASSERT(eval_out==eval_ref);

        t=0.5;
        eval_out=pwc0.f(t);
        eval_ref=pwc1.f(t);
        TEST_ASSERT(eval_out==eval_ref);
        eval_out=pwc0.fp(t);
        eval_ref=pwc1.fp(t);
        TEST_ASSERT(eval_out==eval_ref);
        eval_out=pwc0.fpp(t);
        eval_ref=pwc1.fpp(t);
        TEST_ASSERT(eval_out==eval_ref);
        eval_out=pwc0.fppp(t);
        eval_ref=pwc1.fppp(t);
        TEST_ASSERT(eval_out==eval_ref);

        t=1;
        eval_out=pwc0.f(t);
        eval_ref=pwc1.f(t);
        TEST_ASSERT(eval_out==eval_ref);
        eval_out=pwc0.fp(t);
        eval_ref=pwc1.fp(t);
        TEST_ASSERT(eval_out==eval_ref);
        eval_out=pwc0.fpp(t);
        eval_ref=pwc1.fpp(t);
        TEST_ASSERT(eval_out==eval_ref);
        eval_out=pwc0.fppp(t);
        eval_ref=pwc1.fppp(t);
        TEST_ASSERT(eval_out==eval_ref);

        t=1.25;
        eval_out=pwc0.f(t);
        eval_ref=pwc1.f(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<3*eps);
        eval_out=pwc0.fp(t);
        eval_ref=pwc1.fp(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<3*eps);
        eval_out=pwc0.fpp(t);
        eval_ref=pwc1.fpp(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<18*eps);
        eval_out=pwc0.fppp(t);
        eval_ref=pwc1.fppp(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<138*eps);

        t=1.5;
        eval_out=pwc0.f(t);
        eval_ref=pwc1.f(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<3*eps);
        eval_out=pwc0.fp(t);
        eval_ref=pwc1.fp(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<7*eps);
        eval_out=pwc0.fpp(t);
        eval_ref=pwc1.fpp(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<35*eps);
        eval_out=pwc0.fppp(t);
        eval_ref=pwc1.fppp(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<138*eps);

        t=1.75;
        eval_out=pwc0.f(t);
        eval_ref=pwc1.f(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<2*eps);
        eval_out=pwc0.fp(t);
        eval_ref=pwc1.fp(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<3*eps);
        eval_out=pwc0.fpp(t);
        eval_ref=pwc1.fpp(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<28*eps);
        eval_out=pwc0.fppp(t);
        eval_ref=pwc1.fppp(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<142*eps);

        t=2;
        eval_out=pwc0.f(t);
        eval_ref=pwc1.f(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<1*eps);
        eval_out=pwc0.fp(t);
        eval_ref=pwc1.fp(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<9*eps);
        eval_out=pwc0.fpp(t);
        eval_ref=pwc1.fpp(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<63*eps);
        eval_out=pwc0.fppp(t);
        eval_ref=pwc1.fppp(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<142*eps);
#endif
    }

    void area_test()
    {
#if 0
      data_type eps(std::numeric_limits<data__>::epsilon());
#ifdef ELI_QD_FOUND
      if ( (typeid(data_type)==typeid(dd_real)) || (typeid(data_type)==typeid(qd_real)) )
        eps=std::numeric_limits<double>::epsilon();
#endif
      piecewise_surface_type c1;
      curve_type bc[3];
      data_type dt[3], len, bc_len[3], ref_len, t0, t1, temp0, temp1;
      typename curve_type::control_point_type cntrl1_in[4], cntrl2_in[5], cntrl3_in[3];
      typename piecewise_surface_type::error_code err;
      data_type tol(std::sqrt(eps));

      // create piecewise curve
      cntrl1_in[0] << 2.0, 2.0, 0.0;
      cntrl1_in[1] << 1.0, 1.5, 0.0;
      cntrl1_in[2] << 3.5, 0.0, 0.0;
      cntrl1_in[3] << 4.0, 1.0, 0.0;
      dt[0]=0.5;
      bc[0].resize(3);
      for (index_type i=0; i<4; ++i)
      {
        bc[0].set_control_point(cntrl1_in[i], i);
      }
      eli::geom::curve::length(bc_len[0], bc[0], tol);
      cntrl2_in[0] << 4.0, 1.0, 0.0;
      cntrl2_in[1] << 5.0, 2.5, 0.0;
      cntrl2_in[2] << 5.5, 1.0, 0.0;
      cntrl2_in[3] << 6.0, 0.0, 0.0;
      cntrl2_in[4] << 6.5,-0.5, 0.0;
      dt[1]=2.0;
      bc[1].resize(4);
      for (index_type i=0; i<5; ++i)
      {
        bc[1].set_control_point(cntrl2_in[i], i);
      }
      eli::geom::curve::length(bc_len[1], bc[1], tol);
      cntrl3_in[0] << 6.5,-0.5, 0.0;
      cntrl3_in[1] << 6.0,-1.0, 0.0;
      cntrl3_in[2] << 5.5,-2.0, 0.0;
      dt[2]=1.5;
      bc[2].resize(2);
      for (index_type i=0; i<3; ++i)
      {
        bc[2].set_control_point(cntrl3_in[i], i);
      }
      eli::geom::curve::length(bc_len[2], bc[2], tol);
      err=c1.set(bc, bc+3, dt);
      TEST_ASSERT(err==piecewise_surface_type::NO_ERROR);
      TEST_ASSERT(c1.number_segments()==3);

      // create two segment curve calc length of each segment to compare
      eli::geom::curve::length(len, c1, tol);
      ref_len=bc_len[0]+bc_len[1]+bc_len[2];
#ifdef ELI_QD_FOUND
      if ( (typeid(data_type)==typeid(dd_real)) || (typeid(data_type)==typeid(qd_real)) )
      {
        TEST_ASSERT(std::abs(len-ref_len)<3*eps);
      }
      else
#endif
      {
        TEST_ASSERT(len==ref_len);
      }

      // choose part of first segment to calc length and compare
      t0=0.125;
      t1=0.375;
      eli::geom::curve::length(len, c1, t0, t1, tol);
      eli::geom::curve::length(ref_len, bc[0], 0.25, 0.75, tol);
      TEST_ASSERT(len==ref_len);

      // choose part of second segment to calc length and compare
      t0=0.25;
      t1=1.5;
      eli::geom::curve::length(len, c1, t0, t1, tol);
      eli::geom::curve::length(temp0, bc[0], 0.5, 1, tol);
      eli::geom::curve::length(temp1, bc[1], 0, 0.5, tol);
      ref_len=temp0+temp1;
      TEST_ASSERT(std::abs(len-ref_len)<2*tol);

      // choose part of third segment to calc length and compare
      t0=0.25;
      t1=3.25;
      eli::geom::curve::length(len, c1, t0, t1, tol);
      eli::geom::curve::length(temp0, bc[0], 0.5, 1, tol);
      eli::geom::curve::length(temp1, bc[2], 0, 0.5, tol);
      ref_len=temp0+bc_len[1]+temp1;
      TEST_ASSERT(std::abs(len-ref_len)<166*tol);
#endif
    }
};

#endif
