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

#ifndef bezier_curve_test_suite_hpp
#define bezier_curve_test_suite_hpp

#include <cmath>    // std::pow, std::exp

#include <typeinfo> // typeid
#include <string>   // std::string
#include <sstream>  // std::stringstream
#include <iomanip>  // std::setw
#include <limits>   // std::numeric_limits

#include "eli/constants/math.hpp"
#include "eli/geom/point/distance.hpp"
#include "eli/geom/curve/bezier.hpp"
#include "eli/geom/curve/length.hpp"
#include "eli/geom/curve/curvature.hpp"

template<typename data__>
class bezier_curve_test_suite : public Test::Suite
{
  private:
    typedef eli::geom::curve::bezier<data__, 3> bezier_type;
    typedef typename bezier_type::fit_container_type fit_container_type;
    typedef typename fit_container_type::constraint_point_type constraint_point_type;
    typedef typename bezier_type::index_type index_type;
    typedef typename bezier_type::point_type point_type;
    typedef typename bezier_type::data_type data_type;
    typedef typename bezier_type::onedbezcurve oned_type;

  protected:
    void AddTests(const float &)
    {
      // add the tests
      TEST_ADD(bezier_curve_test_suite<float>::assignment_test);
      TEST_ADD(bezier_curve_test_suite<float>::bounding_box_test);
      TEST_ADD(bezier_curve_test_suite<float>::transformation_test);
      TEST_ADD(bezier_curve_test_suite<float>::evaluation_test);
      TEST_ADD(bezier_curve_test_suite<float>::derivative_1_test);
      TEST_ADD(bezier_curve_test_suite<float>::derivative_2_test);
      TEST_ADD(bezier_curve_test_suite<float>::derivative_3_test);
      TEST_ADD(bezier_curve_test_suite<float>::frenet_serret_test);
      TEST_ADD(bezier_curve_test_suite<float>::reverse_test);
      TEST_ADD(bezier_curve_test_suite<float>::reflection_test);
      TEST_ADD(bezier_curve_test_suite<float>::promotion_test);
      TEST_ADD(bezier_curve_test_suite<float>::promotion_to_test);
      TEST_ADD(bezier_curve_test_suite<float>::demotion_test);
      TEST_ADD(bezier_curve_test_suite<float>::degree_to_cubic_test);
      TEST_ADD(bezier_curve_test_suite<float>::distance_bound_test);
      TEST_ADD(bezier_curve_test_suite<float>::split_test);
      TEST_ADD(bezier_curve_test_suite<float>::length_test);
      TEST_ADD(bezier_curve_test_suite<float>::math_test);
      TEST_ADD(bezier_curve_test_suite<float>::integral_test);
      TEST_ADD(bezier_curve_test_suite<float>::demo_pt_project_test);
    }
    void AddTests(const double &)
    {
      // add the tests
      TEST_ADD(bezier_curve_test_suite<double>::assignment_test);
      TEST_ADD(bezier_curve_test_suite<double>::bounding_box_test);
      TEST_ADD(bezier_curve_test_suite<double>::transformation_test);
      TEST_ADD(bezier_curve_test_suite<double>::evaluation_test);
      TEST_ADD(bezier_curve_test_suite<double>::derivative_1_test);
      TEST_ADD(bezier_curve_test_suite<double>::derivative_2_test);
      TEST_ADD(bezier_curve_test_suite<double>::derivative_3_test);
      TEST_ADD(bezier_curve_test_suite<double>::frenet_serret_test);
      TEST_ADD(bezier_curve_test_suite<double>::reverse_test);
      TEST_ADD(bezier_curve_test_suite<double>::reflection_test);
      TEST_ADD(bezier_curve_test_suite<double>::promotion_test);
      TEST_ADD(bezier_curve_test_suite<double>::promotion_to_test);
      TEST_ADD(bezier_curve_test_suite<double>::demotion_test);
      TEST_ADD(bezier_curve_test_suite<double>::degree_to_cubic_test);
      TEST_ADD(bezier_curve_test_suite<double>::distance_bound_test);
      TEST_ADD(bezier_curve_test_suite<double>::split_test);
      TEST_ADD(bezier_curve_test_suite<double>::length_test);
      TEST_ADD(bezier_curve_test_suite<double>::math_test);
      TEST_ADD(bezier_curve_test_suite<double>::integral_test);
      TEST_ADD(bezier_curve_test_suite<double>::demo_pt_project_test);
    }
    void AddTests(const long double &)
    {
      // add the tests
      TEST_ADD(bezier_curve_test_suite<long double>::assignment_test);
      TEST_ADD(bezier_curve_test_suite<long double>::bounding_box_test);
      TEST_ADD(bezier_curve_test_suite<long double>::transformation_test);
      TEST_ADD(bezier_curve_test_suite<long double>::evaluation_test);
      TEST_ADD(bezier_curve_test_suite<long double>::derivative_1_test);
      TEST_ADD(bezier_curve_test_suite<long double>::derivative_2_test);
      TEST_ADD(bezier_curve_test_suite<long double>::derivative_3_test);
      TEST_ADD(bezier_curve_test_suite<long double>::frenet_serret_test);
      TEST_ADD(bezier_curve_test_suite<long double>::reverse_test);
      TEST_ADD(bezier_curve_test_suite<long double>::reflection_test);
      TEST_ADD(bezier_curve_test_suite<long double>::promotion_test);
      TEST_ADD(bezier_curve_test_suite<long double>::promotion_to_test);
      TEST_ADD(bezier_curve_test_suite<long double>::demotion_test);
      TEST_ADD(bezier_curve_test_suite<long double>::degree_to_cubic_test);
      TEST_ADD(bezier_curve_test_suite<long double>::distance_bound_test);
      TEST_ADD(bezier_curve_test_suite<long double>::split_test);
      TEST_ADD(bezier_curve_test_suite<long double>::length_test);
      TEST_ADD(bezier_curve_test_suite<long double>::math_test);
      TEST_ADD(bezier_curve_test_suite<long double>::integral_test);
      TEST_ADD(bezier_curve_test_suite<long double>::demo_pt_project_test);
    }

  public:
    bezier_curve_test_suite()
    {
      AddTests(data__());
    }
    ~bezier_curve_test_suite()
    {
    }

  private:
    void octave_print(int figno, const std::vector<point_type> &pts, const bezier_type &bez) const
    {
      size_t i;

      std::cout << "figure(" << figno << ");" << std::endl;
      std::cout << "xpts=[" << pts[0].x();
      for (i=1; i<pts.size(); ++i)
        std::cout << ", " << pts[i].x();
      std::cout << "];" << std::endl;
      std::cout << "ypts=[" << pts[0].y();
      for (i=1; i<pts.size(); ++i)
        std::cout << ", " << pts[i].y();
      std::cout << "];" << std::endl;

      std::vector<data__> t(101);
      for (i=0; i<t.size(); ++i)
        t[i]=static_cast<data__>(i)/(t.size()-1);

      std::cout << "xint=[" << bez.f(t[0])(0);
      for (i=1; i<t.size(); ++i)
        std::cout << ", " << bez.f(t[i])(0);
      std::cout << "];" << std::endl;
      std::cout << "yint=[" << bez.f(t[0])(1);
      for (i=1; i<t.size(); ++i)
        std::cout << ", " << bez.f(t[i])(1);
      std::cout << "];" << std::endl;

      std::cout << "plot(xpts, ypts, 'bo', xint, yint, 'k-');" << std::endl;
    }

    void create_circle(std::vector<point_type> &pts)
    {
      // NOTE: This will not create a closed circle, the last point will be
      //       the point just prior to the 2*pi point
      size_t n=pts.size();
      for (size_t i=0; i<n; ++i)
      {
        data__ theta(eli::constants::math<data__>::two_pi()*static_cast<data__>(i)/n);
        pts[i](0)=std::cos(theta);
        pts[i](1)=std::sin(theta);
        pts[i](2)=0;
      }
    }

    void assignment_test()
    {
      bezier_type bc1, bc2;
      point_type cntrl_in[4];
      index_type i;

      // test default constructor then set control points
      cntrl_in[0] << 2, 2, 0;
      cntrl_in[1] << 1, 1, 0;
      cntrl_in[2] << 3, 0, 0;
      cntrl_in[3] << 4, 1, 0;
      bc1.resize(3);
      for (i=0; i<4; ++i)
      {
        bc1.set_control_point(cntrl_in[i], i);
      }
      for (i=0; i<4; ++i)
      {
        TEST_ASSERT(cntrl_in[i]==bc1.get_control_point(i));
      }

      // test copy ctr
      bezier_type bcc1(bc1);
      for (i=0; i<4; ++i)
      {
        TEST_ASSERT(bc1.get_control_point(i)==bcc1.get_control_point(i));
      }

      // test assignment operator
      bc2=bc1;
      for (i=0; i<4; ++i)
      {
        TEST_ASSERT(bc1.get_control_point(i)==bc2.get_control_point(i));
      }

      // test order
      TEST_ASSERT(bc2.degree()==4-1);
    }

    void bounding_box_test()
    {
      bezier_type bc;
      point_type cntrl_in[4];
      index_type i;

      // test default constructor then set control points
      cntrl_in[0] << 2, 2, 0;
      cntrl_in[1] << 1, 1, 0;
      cntrl_in[2] << 3, 0, 0;
      cntrl_in[3] << 4, 1, 0;
      bc.resize(3);
      for (i=0; i<4; ++i)
      {
        bc.set_control_point(cntrl_in[i], i);
      }

      // test bounding box
      typename bezier_type::bounding_box_type bb;
      point_type pmin_ref, pmax_ref;

      bc.get_bounding_box(bb);
      pmin_ref << 1, 0, 0;
      pmax_ref << 4, 2, 0;
      TEST_ASSERT(bb.get_min()==pmin_ref);
      TEST_ASSERT(bb.get_max()==pmax_ref);
    }

    void transformation_test()
    {
      point_type cntrl_in[4];
      data_type eps(std::numeric_limits<data__>::epsilon());

      // set control points
      cntrl_in[0] << 2.0, 2.0, 0.0;
      cntrl_in[1] << 1.0, 1.5, 0.0;
      cntrl_in[2] << 3.5, 0.0, 0.0;
      cntrl_in[3] << 4.0, 1.0, 0.0;

      bezier_type bc1(3), bc2;
      point_type  eval_out, eval_ref;
      data_type  t;

      // set control points
      for (index_type i=0; i<4; ++i)
      {
        bc1.set_control_point(cntrl_in[i], i);
      }

      // test translation
      {
        point_type trans;

        // set up translation vector and apply
        bc2=bc1;
        trans << 2, 1, 3;
        bc2.translate(trans);

        // test evaluations
        t=0;
        eval_out=bc2.f(t);
        eval_ref=trans+cntrl_in[0];
        TEST_ASSERT(eval_out==eval_ref);
        t=1;
        eval_out=bc2.f(t);
        eval_ref=trans+cntrl_in[3];
        TEST_ASSERT(eval_out==eval_ref);

        // test evaluation at interior point
        t=static_cast<data__>(0.45);
        eval_out=bc2.f(t);
        eval_ref << static_cast<data__>(2.2750625), static_cast<data__>(1.0364375), static_cast<data__>(0);
        eval_ref+=trans;
        TEST_ASSERT((eval_out-eval_ref).norm()<5e3*eps);
      }

      // test rotation about origin
      {
        typename bezier_type::rotation_matrix_type rmat;

        // set up rotation and apply
        data_type one(1);
        bc2=bc1;
        rmat << std::cos(one), 0,  -std::sin(one),
                0,             one, 0,
                std::sin(one), 0,   std::cos(one);
        bc2.rotate(rmat);

        // test evaluations
        t=0;
        eval_out=bc2.f(t);
        eval_ref=(cntrl_in[0]*rmat.transpose());
        TEST_ASSERT(eval_out==eval_ref);
        t=1;
        eval_out=bc2.f(t);
        eval_ref=(cntrl_in[3]*rmat.transpose());
        TEST_ASSERT(eval_out==eval_ref);

        // test evaluation at interior point
        t=static_cast<data__>(0.45);
        eval_out=bc2.f(t);
        eval_ref << static_cast<data__>(2.2750625), static_cast<data__>(1.0364375), static_cast<data__>(0);
        eval_ref*=rmat.transpose();
        TEST_ASSERT((eval_out-eval_ref).norm()<5e3*eps);
      }

      // test rotation about point
      {
        point_type rorig;
        typename bezier_type::rotation_matrix_type rmat;

        // set up rotation and apply
        data_type one(1);
        bc2=bc1;
        rorig << 2, 1, 3;
        rmat << std::cos(one), 0,  -std::sin(one),
                0,             one, 0,
                std::sin(one), 0,   std::cos(one);
        bc2.rotate(rmat, rorig);

        // test evaluations
        t=0;
        eval_out=bc2.f(t);
        eval_ref=rorig+(cntrl_in[0]-rorig)*rmat.transpose();
        TEST_ASSERT(eval_out==eval_ref);
        t=1;
        eval_out=bc2.f(t);
        eval_ref=rorig+(cntrl_in[3]-rorig)*rmat.transpose();
        TEST_ASSERT(eval_out==eval_ref);

        // test evaluation at interior point
        t=static_cast<data__>(0.45);
        eval_out=bc2.f(t);
        eval_ref << static_cast<data__>(2.2750625), static_cast<data__>(1.0364375), static_cast<data__>(0);
        eval_ref=rorig+(eval_ref-rorig)*rmat.transpose();
        TEST_ASSERT((eval_out-eval_ref).norm()<5e3*eps);
      }
    }

    void evaluation_test()
    {
      point_type cntrl_in[4];
      data_type eps(std::numeric_limits<data__>::epsilon());

      // set control points
      cntrl_in[0] << 2.0, 2.0, 0.0;
      cntrl_in[1] << 1.0, 1.5, 0.0;
      cntrl_in[2] << 3.5, 0.0, 0.0;
      cntrl_in[3] << 4.0, 1.0, 0.0;

      bezier_type bc1(3);
      point_type  eval_out, eval_ref;
      data_type  t;

      // set control points
      for (index_type i=0; i<4; ++i)
      {
        bc1.set_control_point(cntrl_in[i], i);
      }

      // test evaluation at end points
      t=0;
      eval_out=bc1.f(t);
      TEST_ASSERT(eval_out==cntrl_in[0]);
      t=1;
      eval_out=bc1.f(t);
      TEST_ASSERT(eval_out==cntrl_in[3]);

      // test evaluation at interior point
      t=static_cast<data__>(0.45);
      eval_out=bc1.f(t);
      eval_ref << static_cast<data__>(2.2750625), static_cast<data__>(1.0364375), static_cast<data__>(0);
      TEST_ASSERT((eval_out-eval_ref).norm()<5e3*eps);
    }

    void derivative_1_test()
    {
      point_type cntrl_in[4];
      data_type eps(std::numeric_limits<data__>::epsilon());

      // set control points
      cntrl_in[0] << 2.0, 2.0, 0.0;
      cntrl_in[1] << 1.0, 1.5, 0.0;
      cntrl_in[2] << 3.5, 0.0, 0.0;
      cntrl_in[3] << 4.0, 1.0, 0.0;

      bezier_type bc1(3);
      bezier_type bp;
      point_type eval_out, eval_ref;
      data_type t;

      // set control points
      for (index_type i=0; i<4; ++i)
      {
        bc1.set_control_point(cntrl_in[i], i);
      }

      bc1.fp( bp );

      // test 1st derivative at end points
      t=0;
      eval_out=bc1.fp(t);
      eval_ref << -3, -1.5, 0;
      TEST_ASSERT(eval_out==eval_ref);
      eval_out=bp.f(t);
      TEST_ASSERT(eval_out==eval_ref);
      t=1;
      eval_out=bc1.fp(t);
      eval_ref << 1.5, 3, 0;
      TEST_ASSERT(eval_out==eval_ref);
      eval_out=bp.f(t);
      TEST_ASSERT(eval_out==eval_ref);

      // test 1st derivative at interior point
      t=static_cast<data__>(0.45);
      eval_out=bc1.fp(t);
      eval_ref << static_cast<data__>(3.10875), static_cast<data__>(-2.07375), static_cast<data__>(0);
      TEST_ASSERT((eval_out-eval_ref).norm() < 5e3*eps);
      eval_out=bp.f(t);
      TEST_ASSERT((eval_out-eval_ref).norm() < 5e3*eps);
    }

    void derivative_2_test()
    {
      point_type cntrl_in[4];
      data_type eps(std::numeric_limits<data__>::epsilon());

      // set control points
      cntrl_in[0] << 2.0, 2.0, 0.0;
      cntrl_in[1] << 1.0, 1.5, 0.0;
      cntrl_in[2] << 3.5, 0.0, 0.0;
      cntrl_in[3] << 4.0, 1.0, 0.0;

      bezier_type bc1(3);
      point_type eval_out, eval_ref;
      data_type t;

      // set control points
      for (index_type i=0; i<4; ++i)
      {
        bc1.set_control_point(cntrl_in[i], i);
      }

      // test 2nd derivative at end points
      t=0;
      eval_out=bc1.fpp(t);
      eval_ref << 21, -6, 0;
      TEST_ASSERT(eval_out==eval_ref);
      t=1;
      eval_out=bc1.fpp(t);
      eval_ref << -12, 15, 0;
      TEST_ASSERT(eval_out==eval_ref);

      // test 2nd derivative at interior point
      t=static_cast<data__>(0.45);
      eval_out=bc1.fpp(t);
      eval_ref << static_cast<data__>(6.15), static_cast<data__>(3.45), static_cast<data__>(0);
      TEST_ASSERT_DELTA((eval_out-eval_ref).norm(), 0, 1e4*eps);
    }

    void derivative_3_test()
    {
      point_type cntrl_in[4];
      data_type eps(std::numeric_limits<data__>::epsilon());

      // set control points
      cntrl_in[0] << 2.0, 2.0, 0.0;
      cntrl_in[1] << 1.0, 1.5, 0.0;
      cntrl_in[2] << 3.5, 0.0, 0.0;
      cntrl_in[3] << 4.0, 1.0, 0.0;

      bezier_type bc1(3);
      point_type eval_out, eval_ref;
      data_type t;

      // set control points
      for (index_type i=0; i<4; ++i)
      {
        bc1.set_control_point(cntrl_in[i], i);
      }

      // test 3rd derivative at end points
      t=0;
      eval_out=bc1.fppp(t);
      eval_ref << -33, 21, 0;
      TEST_ASSERT(eval_out==eval_ref);
      t=1;
      eval_out=bc1.fppp(t);
      eval_ref << -33, 21, 0;
      TEST_ASSERT(eval_out==eval_ref);

      // test 3rd derivative at interior point
      t=static_cast<data__>(0.45);
      eval_out=bc1.fppp(t);
      eval_ref << -33, 21, 0;
      TEST_ASSERT(eval_out==eval_ref);

      // test curvature at end points
      data_type curv_out, curv_ref;
      t=0;
      eli::geom::curve::curvature(curv_out, bc1, t);
      curv_ref=static_cast<data__>(1.31182654679988);
      TEST_ASSERT_DELTA(curv_out, curv_ref, std::max(static_cast<data__>(1e-13), static_cast<data__>(1e2)*eps));
      t=1;
      eli::geom::curve::curvature(curv_out, bc1, t);
      curv_ref=static_cast<data__>(1.55034046439985);
      TEST_ASSERT_DELTA(curv_out, curv_ref, std::max(static_cast<data__>(1e-13), static_cast<data__>(1e2)*eps));

      // test curvature at interior point
      t=static_cast<data__>(0.45);
      eli::geom::curve::curvature(curv_out, bc1, t);
      curv_ref=static_cast<data__>(0.449908807121445);
      TEST_ASSERT_DELTA(curv_out, curv_ref, std::max(static_cast<data__>(1e-13), static_cast<data__>(1e2)*eps));
    }

    void frenet_serret_test()
    {
      point_type cntrl_in[4];
      data_type eps(std::numeric_limits<data__>::epsilon());

      // set control points
      cntrl_in[0] << 2.0, 2.0, 0.0;
      cntrl_in[1] << 1.0, 1.5, 0.0;
      cntrl_in[2] << 3.5, 0.0, 0.0;
      cntrl_in[3] << 4.0, 1.0, 0.0;

      bezier_type bc(3);
      point_type  t_dir, n_dir, b_dir, t_dir_ref;
      data_type  t;

      // set control points
      for (index_type i=0; i<4; ++i)
      {
        bc.set_control_point(cntrl_in[i], i);
      }

      t=0.5;
      t_dir=bc.tangent(t);
      t_dir_ref << static_cast<data_type>(0.874157), static_cast<data_type>(-0.485643), 0;
      TEST_ASSERT((t_dir-t_dir_ref).norm()<1e-5);

      bc.frenet_serret_frame(t_dir, n_dir, b_dir, t);
      TEST_ASSERT((t_dir.cross(n_dir)-b_dir).norm()<2*eps);
    }

    void reverse_test()
    {
      point_type cntrl_in[4];

      // set control points
      cntrl_in[0] << 2.0, 2.0, 0.0;
      cntrl_in[1] << 1.0, 1.5, 0.0;
      cntrl_in[2] << 3.5, 0.0, 0.0;
      cntrl_in[3] << 4.0, 1.0, 0.0;

      bezier_type bc1(3);

      // set control points
      for (index_type i=0; i<4; ++i)
      {
        bc1.set_control_point(cntrl_in[i], i);
      }

      // reverse curve parameterization
      bc1.reverse();

      // get control points and check them
      for (index_type i=0; i<4; ++i)
      {
        TEST_ASSERT(bc1.get_control_point(i)==cntrl_in[3-i]);
      }
    }

    void reflection_test()
    {
      point_type cntrl_in[4], cntrl_out[4];
      index_type i;

      // set control points
      cntrl_in[0] << 2.0, 2.0, 1.0;
      cntrl_in[1] << 1.0, 1.5, 2.0;
      cntrl_in[2] << 3.5, 0.0, 3.0;
      cntrl_in[3] << 4.0, 1.0, 4.0;

      bezier_type bc(3), bcc;

      // set control points
      for (index_type i=0; i<4; ++i)
      {
        bc.set_control_point(cntrl_in[i], i);
      }
      bcc=bc;

      // reflect about xy-plane
      {
        bc.reflect_xy();
        for (i=0; i<4; ++i)
        {
          cntrl_out[i].x()= cntrl_in[i].x();
          cntrl_out[i].y()= cntrl_in[i].y();
          cntrl_out[i].z()=-cntrl_in[i].z();
        }

        // check control points
        for (i=0; i<4; ++i)
        {
          TEST_ASSERT(bc.get_control_point(i)==cntrl_out[i]);
        }
      }

      // reflect about xz-plane
      bc=bcc;
      {
        bc.reflect_xz();
        for (i=0; i<4; ++i)
        {
          cntrl_out[i].x()= cntrl_in[i].x();
          cntrl_out[i].y()=-cntrl_in[i].y();
          cntrl_out[i].z()= cntrl_in[i].z();
        }

        // check control points
        for (i=0; i<4; ++i)
        {
          TEST_ASSERT(bc.get_control_point(i)==cntrl_out[i]);
        }
      }

      // reflect about yz-plane
      bc=bcc;
      {
        bc.reflect_yz();
        for (i=0; i<4; ++i)
        {
          cntrl_out[i].x()=-cntrl_in[i].x();
          cntrl_out[i].y()= cntrl_in[i].y();
          cntrl_out[i].z()= cntrl_in[i].z();
        }

        // check control points
        for (i=0; i<4; ++i)
        {
          TEST_ASSERT(bc.get_control_point(i)==cntrl_out[i]);
        }
      }

      // reflect about arbitrary plane
      bc=bcc;
      {
        data_type eps(std::numeric_limits<data__>::epsilon());
        point_type normal, n;

        normal << 1, 2, 3;
        n=normal;
        n.normalize();

        bc.reflect(normal);
        for (i=0; i<4; ++i)
        {
          cntrl_out[i]=cntrl_in[i]-2*cntrl_in[i].dot(n)*n;
        }

        // check control points
        for (i=0; i<4; ++i)
        {
          TEST_ASSERT((bc.get_control_point(i)-cntrl_out[i]).norm()<10*eps);
        }
      }

      // reflect about arbitrary plane not through origin
      bc=bcc;
      {
        data_type eps(std::numeric_limits<data__>::epsilon());
        point_type normal, n;
        data_type d(4);

        normal << 1, 2, 3;
        n=normal;
        n.normalize();

        bc.reflect(normal, d);
        for (i=0; i<4; ++i)
        {
          cntrl_out[i]=cntrl_in[i]-2*(cntrl_in[i].dot(n)-d)*n;
        }

        // check control points
        for (i=0; i<4; ++i)
        {
          TEST_ASSERT((bc.get_control_point(i)-cntrl_out[i]).norm()<10*eps);
        }
      }
    }

    void promotion_test()
    {
      point_type cntrl_in[5];
      data_type eps(std::numeric_limits<data__>::epsilon());

      // set control points
      cntrl_in[0] << 0,   0, 0;
      cntrl_in[1] << 0,   4, 0;
      cntrl_in[2] << 2,   4, 0;
      cntrl_in[3] << 2,   3, 0;
      cntrl_in[4] << 1.5, 3, 0;

      bezier_type bc1(4), bc2;
      point_type eval_out, eval_ref;
      data_type t, curv_out, curv_ref;

      // set control points
      for (index_type i=0; i<5; ++i)
      {
        bc1.set_control_point(cntrl_in[i], i);
      }
      bc2=bc1;

      // promote curve
      bc2.degree_promote();

      // test to see if degree has increased
      TEST_ASSERT(bc1.degree()+1==bc2.degree());

      // test to see if get expected polygon
      if (bc1.degree()+1==bc2.degree())
      {
        point_type cntrl_ref[6];

        cntrl_ref[0] << 0.0,                      0.0,                      0.0;
        cntrl_ref[1] << 0.0,                      static_cast<data__>(3.2), 0.0;
        cntrl_ref[2] << static_cast<data__>(1.2), 4.0,                      0.0;
        cntrl_ref[3] << 2.0,                      static_cast<data__>(3.6), 0.0;
        cntrl_ref[4] << static_cast<data__>(1.9), 3.0,                      0.0;
        cntrl_ref[5] << static_cast<data__>(1.5), 3.0,                      0.0;
        for (index_type i=0; i<6; ++i)
        {
          if (typeid(data__)==typeid(long double))
          {
            TEST_ASSERT((bc2.get_control_point(i)-cntrl_ref[i]).norm()<std::sqrt(eps));
          }
          else
          {
            TEST_ASSERT(bc2.get_control_point(i)==cntrl_ref[i]);
          }
        }
      }

      // test evaluation at end points
      t=0;
      eval_out=bc2.f(t);
      eval_ref=bc1.f(t);
      TEST_ASSERT(eval_out==eval_ref);
      t=1;
      eval_out=bc2.f(t);
      eval_ref=bc1.f(t);
      TEST_ASSERT(eval_out==eval_ref);

      // test evaluation at interior point
      t=static_cast<data__>(0.45);
      eval_out=bc2.f(t);
      eval_ref=bc1.f(t);
      TEST_ASSERT((eval_out-eval_ref).norm()<5*eps);

      // test 1st derivative at end points
      t=0;
      eval_out=bc2.fp(t);
      eval_ref=bc1.fp(t);
      TEST_ASSERT(eval_out==eval_ref);
      t=1;
      eval_out=bc2.fp(t);
      eval_ref=bc1.fp(t);
      TEST_ASSERT((eval_out-eval_ref).norm()<5*eps);

      // test 1st derivative at interior point
      t=static_cast<data__>(0.45);
      eval_out=bc2.fp(t);
      eval_ref=bc1.fp(t);
      TEST_ASSERT((eval_out-eval_ref).norm()<5*eps);

      // test 2nd derivative at end points
      t=0;
      eval_out=bc2.fpp(t);
      eval_ref=bc1.fpp(t);
      TEST_ASSERT((eval_out-eval_ref).norm()<33*eps);
      t=1;
      eval_out=bc2.fpp(t);
      eval_ref=bc1.fpp(t);
      TEST_ASSERT((eval_out-eval_ref).norm()<18*eps);

      // test 2nd derivative at interior point
      t=static_cast<data__>(0.45);
      eval_out=bc2.fpp(t);
      eval_ref=bc1.fpp(t);
      TEST_ASSERT((eval_out-eval_ref).norm()<33*eps);

      // test 3rd derivative at end points
      t=0;
      eval_out=bc2.fppp(t);
      eval_ref=bc1.fppp(t);
      TEST_ASSERT((eval_out-eval_ref).norm()<390*eps);
      t=1;
      eval_out=bc2.fppp(t);
      eval_ref=bc1.fppp(t);
      TEST_ASSERT((eval_out-eval_ref).norm()<390*eps);

      // test 3rd derivative at interior point
      t=static_cast<data__>(0.45);
      eval_out=bc2.fppp(t);
      eval_ref=bc1.fppp(t);
      TEST_ASSERT((eval_out-eval_ref).norm()<203*eps);

      // test curvature at end points
      t=0;
      eli::geom::curve::curvature(curv_out, bc2, t);
      eli::geom::curve::curvature(curv_ref, bc1, t);
      TEST_ASSERT(curv_out==curv_ref);
      t=1;
      eli::geom::curve::curvature(curv_out, bc2, t);
      eli::geom::curve::curvature(curv_ref, bc1, t);
      TEST_ASSERT(std::abs(curv_out-curv_ref)<203*eps);

      // test curvature at interior point
      t=static_cast<data__>(0.45);
      eli::geom::curve::curvature(curv_out, bc2, t);
      eli::geom::curve::curvature(curv_ref, bc1, t);
      TEST_ASSERT(std::abs(curv_out-curv_ref)<203*eps);
    }

    void promotion_to_test()
    {
      point_type cntrl_in[5];
      data_type eps(std::numeric_limits<data__>::epsilon());

      // set control points
      cntrl_in[0] << 0,   0, 0;
      cntrl_in[1] << 0,   4, 0;
      cntrl_in[2] << 2,   4, 0;
      cntrl_in[3] << 2,   3, 0;
      cntrl_in[4] << 1.5, 3, 0;

      bezier_type bc1(4), bc2, bc3, bc4;
      point_type eval_out, eval_ref;
      data_type t, curv_out, curv_ref;

      // set control points
      for (index_type i=0; i<5; ++i)
      {
        bc1.set_control_point(cntrl_in[i], i);
      }
      bc2=bc1;
      bc3=bc1;
      bc4=bc1;

      // promote curve single degree
      bc2.degree_promote();
      // promote_to curve
      bc3.degree_promote_to(bc1.degree()+1);

      // test to see if degree has increased
      TEST_ASSERT(bc3.degree()==bc2.degree());

      // test to see if get exact behavior of single-degree promotion
      if (bc3.degree()==bc2.degree())
      {
        for (index_type i=0; i<6; ++i)
        {
          TEST_ASSERT(bc3.get_control_point(i)==bc2.get_control_point(i));
        }
      }

      // promote_to curve high order
      bc4.degree_promote_to(bc1.degree()+5);

      // test to see if degree has increased
      TEST_ASSERT(bc4.degree()==bc1.degree()+5);

      // test evaluation at end points
      t=0;
      eval_out=bc4.f(t);
      eval_ref=bc1.f(t);
      TEST_ASSERT(eval_out==eval_ref);
      t=1;
      eval_out=bc4.f(t);
      eval_ref=bc1.f(t);
      TEST_ASSERT(eval_out==eval_ref);

      // test evaluation at interior point
      t=static_cast<data__>(0.45);
      eval_out=bc4.f(t);
      eval_ref=bc1.f(t);
      TEST_ASSERT((eval_out-eval_ref).norm()<5*eps);

      // test 1st derivative at end points
      t=0;
      eval_out=bc4.fp(t);
      eval_ref=bc1.fp(t);
      TEST_ASSERT((eval_out-eval_ref).norm()<20*eps);
      t=1;
      eval_out=bc4.fp(t);
      eval_ref=bc1.fp(t);
      TEST_ASSERT((eval_out-eval_ref).norm()<5*4*eps);

      // test 1st derivative at interior point
      t=static_cast<data__>(0.45);
      eval_out=bc4.fp(t);
      eval_ref=bc1.fp(t);
      TEST_ASSERT((eval_out-eval_ref).norm()<5*eps);

      // test 2nd derivative at end points
      t=0;
      eval_out=bc4.fpp(t);
      eval_ref=bc1.fpp(t);
      TEST_ASSERT((eval_out-eval_ref).norm()<200*eps);
      t=1;
      eval_out=bc4.fpp(t);
      eval_ref=bc1.fpp(t);
      TEST_ASSERT((eval_out-eval_ref).norm()<18*16*eps);

      // test 2nd derivative at interior point
      t=static_cast<data__>(0.45);
      eval_out=bc4.fpp(t);
      eval_ref=bc1.fpp(t);
      TEST_ASSERT((eval_out-eval_ref).norm()<35*eps);

      // test 3rd derivative at end points
      t=0;
      eval_out=bc4.fppp(t);
      eval_ref=bc1.fppp(t);
      TEST_ASSERT((eval_out-eval_ref).norm()<390*16*eps);
      t=1;
      eval_out=bc4.fppp(t);
      eval_ref=bc1.fppp(t);
      TEST_ASSERT((eval_out-eval_ref).norm()<390*4*10*eps);

      // test 3rd derivative at interior point
      t=static_cast<data__>(0.45);
      eval_out=bc4.fppp(t);
      eval_ref=bc1.fppp(t);
      TEST_ASSERT((eval_out-eval_ref).norm()<203*4*eps);

      // test curvature at end points
      t=0;
      eli::geom::curve::curvature(curv_out, bc4, t);
      eli::geom::curve::curvature(curv_ref, bc1, t);
      TEST_ASSERT(std::abs(curv_out-curv_ref)<5*eps);
      t=1;
      eli::geom::curve::curvature(curv_out, bc4, t);
      eli::geom::curve::curvature(curv_ref, bc1, t);
      TEST_ASSERT(std::abs(curv_out-curv_ref)<203*eps);

      // test curvature at interior point
      t=static_cast<data__>(0.45);
      eli::geom::curve::curvature(curv_out, bc4, t);
      eli::geom::curve::curvature(curv_ref, bc1, t);
      TEST_ASSERT(std::abs(curv_out-curv_ref)<203*eps);
    }

    void demotion_test()
    {
      data_type eps(std::numeric_limits<data__>::epsilon());

      {
        // hand checked with "Degree Reduction of Bezier Curves" by Dave Morgan
        point_type cntrl_in[7], cntrl_ref[6];

        cntrl_in[0] <<  0, 0, 0;
        cntrl_in[1] <<  2, 6, 0;
        cntrl_in[2] <<  3, 0, 0;
        cntrl_in[3] <<  5, 4, 0;
        cntrl_in[4] <<  7, 1, 0;
        cntrl_in[5] <<  5, 5, 0;
        cntrl_in[6] << 10, 6, 0;

        bezier_type bc1(6), bc2;

        // set control points
        for (index_type i=0; i<7; ++i)
        {
          bc1.set_control_point(cntrl_in[i], i);
        }

        bc2=bc1;
        bc1.degree_demote(eli::geom::general::NOT_CONNECTED);

        // test if the degree is correct
        TEST_ASSERT(bc1.degree()+1==bc2.degree());

        if (bc1.degree()+1==bc2.degree())
        {
          // test of get the correct control points
          cntrl_ref[0] << static_cast<data__>(-0.00878906), static_cast<data__>(0.0610352),  0;
          cntrl_ref[1] << static_cast<data__>( 2.5177734),  static_cast<data__>(6.3821289),  0;
          cntrl_ref[2] << static_cast<data__>( 2.8060547), static_cast<data__>(-0.16982422), 0;
          cntrl_ref[3] << static_cast<data__>( 8.0060547), static_cast<data__>( 2.5301758),  0;
          cntrl_ref[4] << static_cast<data__>( 4.1177734), static_cast<data__>( 3.9821289),  0;
          cntrl_ref[5] << static_cast<data__>( 9.9912109), static_cast<data__>( 6.0610352),  0;
          for (index_type i=0; i<6; ++i)
          {
            TEST_ASSERT_DELTA((bc1.get_control_point(i)-cntrl_ref[i]).norm(), 0, 1e-6);
          }
        }
      }

      // test whether can promote and demote and get the same control points
      {
        point_type cntrl_in[5];

        // set control points
        cntrl_in[0] << 0,   0, 0;
        cntrl_in[1] << 0,   4, 0;
        cntrl_in[2] << 2,   4, 0;
        cntrl_in[3] << 2,   3, 0;
        cntrl_in[4] << 1.5, 3, 0;

        bezier_type bc1(4), bc2;

        // set control points
        for (index_type i=0; i<5; ++i)
        {
          bc1.set_control_point(cntrl_in[i], i);
        }
        bc2=bc1;
        bc2.degree_promote();
        bc2.degree_demote(eli::geom::general::NOT_CONNECTED);

        // test to see if degree has changed
        TEST_ASSERT(bc1.degree()==bc2.degree());

        // test to see if get expected polygon
        for (index_type i=0; i<5; ++i)
        {
          if (typeid(data__)==typeid(long double))
          {
            TEST_ASSERT((bc2.get_control_point(i)-cntrl_in[i]).norm()<std::sqrt(eps));
          }
          else
          {
            TEST_ASSERT((bc2.get_control_point(i)-cntrl_in[i]).norm()<4*eps);
          }
        }
      }

      // test whether can demote and maintain continuity at end points
      {
        point_type cntrl_in[10];

        cntrl_in[0] <<  0, 0, 0,
        cntrl_in[1] <<  2, 6, 0,
        cntrl_in[2] <<  3, 0, 0,
        cntrl_in[3] <<  5, 4, 0,
        cntrl_in[4] <<  7, 1, 0,
        cntrl_in[5] <<  5, 5, 0,
        cntrl_in[6] <<  6, 6, 0,
        cntrl_in[7] <<  7, 5, 0,
        cntrl_in[8] <<  5, 5, 0,
        cntrl_in[9] << 10, 6, 0;

        bezier_type bc1(9);

        // set control points
        for (index_type i=0; i<10; ++i)
        {
          bc1.set_control_point(cntrl_in[i], i);
        }

        // No End Constraint
        {
          bezier_type bc2(bc1);
          point_type eval_out[3], eval_ref[3];
          data_type d[3], tv[3]={0, 1, static_cast<data__>(0.45)};
          bool success;

          success=bc2.degree_demote(eli::geom::general::NOT_CONNECTED);
          TEST_ASSERT(success);

          for (size_t k=0; k<3; ++k)
          {
            data_type t(tv[k]);
            eval_out[0]=bc2.f(t);
            eval_out[1]=bc2.fp(t);
            eval_out[2]=bc2.fpp(t);
            eval_ref[0]=bc1.f(t);
            eval_ref[1]=bc1.fp(t);
            eval_ref[2]=bc1.fpp(t);
            d[0]=eli::geom::point::distance(eval_out[0], eval_ref[0]);
            d[1]=eli::geom::point::distance(eval_out[1], eval_ref[1]);
            d[2]=eli::geom::point::distance(eval_out[2], eval_ref[2]);

            if (k<2)
            {
              TEST_ASSERT_DELTA(d[0], 0.00435372, 5e-4);
              TEST_ASSERT_DELTA(d[1], 0.70530, 5e-3);
              TEST_ASSERT_DELTA(d[2], 37.6161, 5e-2);
            }
            else
            {
              TEST_ASSERT_DELTA(d[0], 0.003414, 1e-5);
              TEST_ASSERT_DELTA(d[1], 0.04887, 5e-4);
              TEST_ASSERT_DELTA(d[2], 1.10759, 5e-4);
            }
          }
        }

        // Point End Constraint
        {
          bezier_type bc2(bc1);
          point_type eval_out[3], eval_ref[3];
          data_type d[3], tv[3]={0, 1, static_cast<data__>(0.45)};
          bool success;

          success=bc2.degree_demote(eli::geom::general::C0);
          TEST_ASSERT(success);

          for (size_t k=0; k<3; ++k)
          {
            data_type t(tv[k]);
            eval_out[0]=bc2.f(t);
            eval_out[1]=bc2.fp(t);
            eval_out[2]=bc2.fpp(t);
            eval_ref[0]=bc1.f(t);
            eval_ref[1]=bc1.fp(t);
            eval_ref[2]=bc1.fpp(t);
            d[0]=eli::geom::point::distance(eval_out[0], eval_ref[0]);
            d[1]=eli::geom::point::distance(eval_out[1], eval_ref[1]);
            d[2]=eli::geom::point::distance(eval_out[2], eval_ref[2]);

            if (k<2)
            {
              TEST_ASSERT_DELTA(d[0], 0.0, 5e-6);
              TEST_ASSERT_DELTA(d[1], 0.64553, 5e-5);
              TEST_ASSERT_DELTA(d[2], 37.4409, 5e-4);
            }
            else
            {
              TEST_ASSERT_DELTA(d[0], 0.00305362, 5e-6);
              TEST_ASSERT_DELTA(d[1], 0.042752, 5e-5);
              TEST_ASSERT_DELTA(d[2], 1.043406, 1e-4);
            }
          }
        }

        // Slope End Constraint
        {
          bezier_type bc2(bc1);
          point_type eval_out[3], eval_ref[3];
          data_type d[3], tv[3]={0, 1, static_cast<data__>(0.45)};
          bool success;

          success=bc2.degree_demote(eli::geom::general::C1);
          TEST_ASSERT(success);

          for (size_t k=0; k<3; ++k)
          {
            data_type t(tv[k]);
            eval_out[0]=bc2.f(t);
            eval_out[1]=bc2.fp(t);
            eval_out[2]=bc2.fpp(t);
            eval_ref[0]=bc1.f(t);
            eval_ref[1]=bc1.fp(t);
            eval_ref[2]=bc1.fpp(t);
            d[0]=eli::geom::point::distance(eval_out[0], eval_ref[0]);
            d[1]=eli::geom::point::distance(eval_out[1], eval_ref[1]);
            d[2]=eli::geom::point::distance(eval_out[2], eval_ref[2]);

            if (k<2)
            {
              TEST_ASSERT(d[0]==0);
              TEST_ASSERT(d[1]==0);
              TEST_ASSERT_DELTA(d[2], 16.7838, 5e-4);
            }
            else
            {
              TEST_ASSERT_DELTA(d[0], 0.0057941, 5e-6);
              TEST_ASSERT_DELTA(d[1], 0.086370, 5e-5);
              TEST_ASSERT_DELTA(d[2], 1.69369, 5e-4);
            }
          }
        }

        // Second Derivative Constraint
        {
          bezier_type bc2(bc1);
          point_type eval_out[3], eval_ref[3];
          data_type d[3], tv[3]={0, 1, static_cast<data__>(0.45)};
          bool success;

          success=bc2.degree_demote(eli::geom::general::C2);
          TEST_ASSERT(success);

          for (size_t k=0; k<3; ++k)
          {
            data_type t(tv[k]);
            eval_out[0]=bc2.f(t);
            eval_out[1]=bc2.fp(t);
            eval_out[2]=bc2.fpp(t);
            eval_ref[0]=bc1.f(t);
            eval_ref[1]=bc1.fp(t);
            eval_ref[2]=bc1.fpp(t);
            d[0]=eli::geom::point::distance(eval_out[0], eval_ref[0]);
            d[1]=eli::geom::point::distance(eval_out[1], eval_ref[1]);
            d[2]=eli::geom::point::distance(eval_out[2], eval_ref[2]);

            if (k<2)
            {
              TEST_ASSERT(d[0]==0);
              TEST_ASSERT(d[1]==0);
              TEST_ASSERT(d[2] < 129*eps);
            }
            else
            {
              TEST_ASSERT_DELTA(d[0], 0.0180029, 5e-5);
              TEST_ASSERT_DELTA(d[1], 0.294979, 5e-5);
              TEST_ASSERT_DELTA(d[2], 3.78229, 5e-4);
            }
          }
        }
      }
    }

    void degree_to_cubic_test()  // degree_to_cubic();
    {
      data_type eps(std::numeric_limits<data__>::epsilon());

      // Test linear promotion case.
      {
        point_type cntrl_in[2];

        cntrl_in[0] <<  0, 0, 0;
        cntrl_in[1] <<  2, 6, 0;

        bezier_type bc1(1), bc2;
        point_type eval_out, eval_ref;
        data_type t, curv_out, curv_ref;

        // set control points
        for (index_type i=0; i<2; ++i)
        {
          bc1.set_control_point(cntrl_in[i], i);
        }

        bc2=bc1;
        bc1.degree_to_cubic();

        // test if the degree is correct
        TEST_ASSERT(bc1.degree()==3);

        // test evaluation at end points
        t=0;
        eval_out=bc2.f(t);
        eval_ref=bc1.f(t);
        TEST_ASSERT(eval_out==eval_ref);
        t=1;
        eval_out=bc2.f(t);
        eval_ref=bc1.f(t);
        TEST_ASSERT(eval_out==eval_ref);

        // test evaluation at interior point
        t=static_cast<data__>(0.45);
        eval_out=bc2.f(t);
        eval_ref=bc1.f(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<5*eps);

        // test 1st derivative at end points
        t=0;
        eval_out=bc2.fp(t);
        eval_ref=bc1.fp(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<5*eps);
        t=1;
        eval_out=bc2.fp(t);
        eval_ref=bc1.fp(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<5*eps);

        // test 1st derivative at interior point
        t=static_cast<data__>(0.45);
        eval_out=bc2.fp(t);
        eval_ref=bc1.fp(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<5*eps);

        // test 2nd derivative at end points
        t=0;
        eval_out=bc2.fpp(t);
        eval_ref=bc1.fpp(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<33*eps);
        t=1;
        eval_out=bc2.fpp(t);
        eval_ref=bc1.fpp(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<18*eps);

        // test 2nd derivative at interior point
        t=static_cast<data__>(0.45);
        eval_out=bc2.fpp(t);
        eval_ref=bc1.fpp(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<33*eps);

        // test 3rd derivative at end points
        t=0;
        eval_out=bc2.fppp(t);
        eval_ref=bc1.fppp(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<390*eps);
        t=1;
        eval_out=bc2.fppp(t);
        eval_ref=bc1.fppp(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<390*eps);

        // test 3rd derivative at interior point
        t=static_cast<data__>(0.45);
        eval_out=bc2.fppp(t);
        eval_ref=bc1.fppp(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<203*eps);

        // test curvature at end points
        t=0;
        eli::geom::curve::curvature(curv_out, bc2, t);
        eli::geom::curve::curvature(curv_ref, bc1, t);
        TEST_ASSERT(curv_out==curv_ref);
        t=1;
        eli::geom::curve::curvature(curv_out, bc2, t);
        eli::geom::curve::curvature(curv_ref, bc1, t);
        TEST_ASSERT(std::abs(curv_out-curv_ref)<203*eps);

        // test curvature at interior point
        t=static_cast<data__>(0.45);
        eli::geom::curve::curvature(curv_out, bc2, t);
        eli::geom::curve::curvature(curv_ref, bc1, t);
        TEST_ASSERT(std::abs(curv_out-curv_ref)<203*eps);
      }

      // Test quadratic promotion case.
      {
        point_type cntrl_in[3];

        cntrl_in[0] <<  0, 0, 0;
        cntrl_in[1] <<  2, 6, 0;
        cntrl_in[2] <<  3, 0, 0;

        bezier_type bc1(2), bc2;
        point_type eval_out, eval_ref;
        data_type t, curv_out, curv_ref;

        // set control points
        for (index_type i=0; i<3; ++i)
        {
          bc1.set_control_point(cntrl_in[i], i);
        }

        bc2=bc1;
        bc1.degree_to_cubic();

        // test if the degree is correct
        TEST_ASSERT(bc1.degree()==3);

        // test evaluation at end points
        t=0;
        eval_out=bc2.f(t);
        eval_ref=bc1.f(t);
        TEST_ASSERT(eval_out==eval_ref);
        t=1;
        eval_out=bc2.f(t);
        eval_ref=bc1.f(t);
        TEST_ASSERT(eval_out==eval_ref);

        // test evaluation at interior point
        t=static_cast<data__>(0.45);
        eval_out=bc2.f(t);
        eval_ref=bc1.f(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<5*eps);

        // test 1st derivative at end points
        t=0;
        eval_out=bc2.fp(t);
        eval_ref=bc1.fp(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<5*eps);
        t=1;
        eval_out=bc2.fp(t);
        eval_ref=bc1.fp(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<5*eps);

        // test 1st derivative at interior point
        t=static_cast<data__>(0.45);
        eval_out=bc2.fp(t);
        eval_ref=bc1.fp(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<5*eps);

        // test 2nd derivative at end points
        t=0;
        eval_out=bc2.fpp(t);
        eval_ref=bc1.fpp(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<33*eps);
        t=1;
        eval_out=bc2.fpp(t);
        eval_ref=bc1.fpp(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<18*eps);

        // test 2nd derivative at interior point
        t=static_cast<data__>(0.45);
        eval_out=bc2.fpp(t);
        eval_ref=bc1.fpp(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<33*eps);

        // test 3rd derivative at end points
        t=0;
        eval_out=bc2.fppp(t);
        eval_ref=bc1.fppp(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<390*eps);
        t=1;
        eval_out=bc2.fppp(t);
        eval_ref=bc1.fppp(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<390*eps);

        // test 3rd derivative at interior point
        t=static_cast<data__>(0.45);
        eval_out=bc2.fppp(t);
        eval_ref=bc1.fppp(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<203*eps);

        // test curvature at end points
        t=0;
        eli::geom::curve::curvature(curv_out, bc2, t);
        eli::geom::curve::curvature(curv_ref, bc1, t);
        TEST_ASSERT((curv_out-curv_ref)<10*eps);
        t=1;
        eli::geom::curve::curvature(curv_out, bc2, t);
        eli::geom::curve::curvature(curv_ref, bc1, t);
        TEST_ASSERT(std::abs(curv_out-curv_ref)<203*eps);

        // test curvature at interior point
        t=static_cast<data__>(0.45);
        eli::geom::curve::curvature(curv_out, bc2, t);
        eli::geom::curve::curvature(curv_ref, bc1, t);
        TEST_ASSERT(std::abs(curv_out-curv_ref)<203*eps);
      }

      // Test do-nothing cubic case.
      {
        point_type cntrl_in[4];

        cntrl_in[0] <<  0, 0, 0;
        cntrl_in[1] <<  2, 6, 0;
        cntrl_in[2] <<  3, 0, 0;
        cntrl_in[3] <<  5, 4, 0;

        bezier_type bc1(3), bc2;

        // set control points
        for (index_type i=0; i<4; ++i)
        {
          bc1.set_control_point(cntrl_in[i], i);
        }

        bc2=bc1;
        bc1.degree_to_cubic();

        // test if the degree is correct
        TEST_ASSERT(bc1.degree()==3);

        // test to see that nothing happened
        if (bc1.degree()==3)
        {
          for (index_type i=0; i<4; ++i)
          {
            TEST_ASSERT(bc1.get_control_point(i)==bc2.get_control_point(i));
          }
        }
      }

      // Test demotion case.
      {
        point_type cntrl_in[7];

        cntrl_in[0] <<  0, 0, 0;
        cntrl_in[1] <<  2, 6, 0;
        cntrl_in[2] <<  3, 0, 0;
        cntrl_in[3] <<  5, 4, 0;
        cntrl_in[4] <<  7, 1, 0;
        cntrl_in[5] <<  5, 5, 0;
        cntrl_in[6] << 10, 6, 0;

        bezier_type bc1(6), bc2;

        // set control points
        for (index_type i=0; i<7; ++i)
        {
          bc1.set_control_point(cntrl_in[i], i);
        }

        bc2=bc1;
        bc1.degree_to_cubic();

        // test if the degree is correct
        TEST_ASSERT(bc1.degree()==3);

        // test that start/end point/derivatives are equal.
        data_type t;
        t = 0;
        TEST_ASSERT(eli::geom::point::distance(bc2.f(t), bc1.f(t))==0);
        TEST_ASSERT(eli::geom::point::distance(bc2.fp(t), bc1.fp(t))==0);

        t = 1;
        TEST_ASSERT(eli::geom::point::distance(bc2.f(t), bc1.f(t))==0);
        TEST_ASSERT(eli::geom::point::distance(bc2.fp(t), bc1.fp(t))==0);
      }
    }

    void distance_bound_test()
    {
      // Test that curve has zero distance itself
      {
        point_type cntrl_in[5];

        // set control points
        cntrl_in[0] << 0,   0, 0;
        cntrl_in[1] << 0,   4, 0;
        cntrl_in[2] << 2,   4, 0;
        cntrl_in[3] << 2,   3, 0;
        cntrl_in[4] << 1.5, 3, 0;

        bezier_type bc1(4);

        // set control points
        for (index_type i=0; i<5; ++i)
        {
          bc1.set_control_point(cntrl_in[i], i);
        }

        data_type d1 = bc1.eqp_distance_bound(bc1);
        TEST_ASSERT(d1==0);
      }

      // Test that curve has zero distance from promoted-self (and vis-versa)
      {
        point_type cntrl_in[5];

        // set control points
        cntrl_in[0] << 0,   0, 0;
        cntrl_in[1] << 0,   4, 0;
        cntrl_in[2] << 2,   4, 0;
        cntrl_in[3] << 2,   3, 0;
        cntrl_in[4] << 1.5, 3, 0;

        bezier_type bc1(4), bc2;

        // set control points
        for (index_type i=0; i<5; ++i)
        {
          bc1.set_control_point(cntrl_in[i], i);
        }
        bc2=bc1;

        // promote_to curve high order
        bc2.degree_promote_to(bc1.degree()+5);

        data_type d1 = bc2.eqp_distance_bound(bc1);
        TEST_ASSERT(d1==0);

        data_type d2 = bc1.eqp_distance_bound(bc2);
        TEST_ASSERT(d2==0);
      }

      // Test that curve has known distance from offset & promoted-self (and vis-versa)
      {
        point_type cntrl_in[5];

        // set control points
        cntrl_in[0] << 0,   0, 0;
        cntrl_in[1] << 0,   4, 0;
        cntrl_in[2] << 2,   4, 0;
        cntrl_in[3] << 2,   3, 0;
        cntrl_in[4] << 1.5, 3, 0;

        data_type dz(3);
        point_type offset;
        offset << 0, 0, dz;

        bezier_type bc1(4), bc2(4);

        // set control points
        for (index_type i=0; i<5; ++i)
        {
          bc1.set_control_point(cntrl_in[i], i);
          bc2.set_control_point(cntrl_in[i]+offset, i);
        }

        // promote_to curve high order
        bc2.degree_promote_to(bc1.degree()+5);

        data_type d1 = bc2.eqp_distance_bound(bc1);
        TEST_ASSERT(d1==dz);

        data_type d2 = bc1.eqp_distance_bound(bc2);
        TEST_ASSERT(d2==dz);
      }

      // Test that curve has known distance from reversed & promoted-self (and vis-versa)
      {
        point_type cntrl_in[5];

        // set control points
        // curve constructed such that greatest distance from reverse will be the endpoints.
        cntrl_in[0] << 0,   0,   0;
        cntrl_in[1] << 1,   1,   0;
        cntrl_in[2] << 2,   0.5, 0;
        cntrl_in[3] << 3,   2.5, 0;
        cntrl_in[4] << 4,   6.0, 0;

        bezier_type bc1(4), bc2(4);

        // set control points
        for (index_type i=0; i<5; ++i)
        {
          bc1.set_control_point(cntrl_in[i], i);
        }
        bc2=bc1;

        bc2.reverse();

        data_type dendpts = (bc1.get_control_point(4)-bc1.get_control_point(0) ).norm();

        // promote_to curve high order
        bc2.degree_promote_to(bc1.degree()+5);

        data_type d1 = bc2.eqp_distance_bound(bc1);
        TEST_ASSERT(d1==dendpts);

        data_type d2 = bc1.eqp_distance_bound(bc2);
        TEST_ASSERT(d2==dendpts);
      }
    }

    void split_test()
    {
      point_type cntrl_in[4], cntrl_ref[4];
      data_type eps(std::numeric_limits<data__>::epsilon());

      // set control points
      cntrl_in[0] << 0, 0, 0;
      cntrl_in[1] << 0, 2, 0;
      cntrl_in[2] << 8, 2, 0;
      cntrl_in[3] << 4, 0, 0;

      bezier_type bc1(3), bc1l, bc1r;
      point_type eval_out, eval_ref;
      data_type t;

      // set control points
      for (index_type i=0; i<4; ++i)
      {
        bc1.set_control_point(cntrl_in[i], i);
      }

      // test split with known control points
      t=0.5;
      bc1.split(bc1l, bc1r, t);
      cntrl_ref[0] << 0.0, 0.0, 0.0;
      cntrl_ref[1] << 0.0, 1.0, 0.0;
      cntrl_ref[2] << 2.0, 1.5, 0.0;
      cntrl_ref[3] << 3.5, 1.5, 0.0;
      for (index_type i=0; i<4; ++i)
      {
        TEST_ASSERT(bc1l.get_control_point(i)==cntrl_ref[i]);
      }
      cntrl_ref[0] << 3.5, 1.5, 0.0;
      cntrl_ref[1] << 5.0, 1.5, 0.0;
      cntrl_ref[2] << 6.0, 1.0, 0.0;
      cntrl_ref[3] << 4.0, 0.0, 0.0;
      for (index_type i=0; i<4; ++i)
      {
        TEST_ASSERT(bc1r.get_control_point(i)==cntrl_ref[i]);
      }

      // split the curve and check the evaluations
      data_type tl, tr, ts;
      tl=static_cast<data__>(0.3);
      tr=static_cast<data__>(0.87);
      ts=static_cast<data__>(0.586);

      bc1.split(bc1l, bc1r, ts);

      // check the left curve
      t=tl*ts;
      eval_out=bc1l.f(tl);
      eval_ref=bc1.f(t);
      TEST_ASSERT((eval_out-eval_ref).norm()<1e2*eps);
      eval_out=bc1l.fp(tl);
      eval_ref=bc1.fp(t)*ts;
      TEST_ASSERT((eval_out-eval_ref).norm()<1e2*eps);
      eval_out=bc1l.fpp(tl);
      eval_ref=bc1.fpp(t)*ts*ts;
      TEST_ASSERT((eval_out-eval_ref).norm()<1e2*eps);
      eval_out=bc1l.fppp(tl);
      eval_ref=bc1.fppp(t)*ts*ts*ts;
      TEST_ASSERT((eval_out-eval_ref).norm()<1e2*eps);

      // check the right curve
      t=ts+tr*(1-ts);
      eval_out=bc1r.f(tr);
      eval_ref=bc1.f(t);
      TEST_ASSERT((eval_out-eval_ref).norm()<1e2*eps);
      eval_out=bc1r.fp(tr);
      eval_ref=bc1.fp(t)*(1-ts);
      TEST_ASSERT((eval_out-eval_ref).norm()<1e2*eps);
      eval_out=bc1r.fpp(tr);
      eval_ref=bc1.fpp(t)*(1-ts)*(1-ts);
      TEST_ASSERT((eval_out-eval_ref).norm()<1e2*eps);
      eval_out=bc1r.fppp(tr);
      eval_ref=bc1.fppp(t)*(1-ts)*(1-ts)*(1-ts);
      TEST_ASSERT((eval_out-eval_ref).norm()<1.21e2*eps);
    }

    void length_test()
    {
      point_type  cntrl_in[4];
      data_type eps(std::numeric_limits<data__>::epsilon());

      // set control points
      cntrl_in[0] << 0, 0, 0;
      cntrl_in[1] << 0, 2, 0;
      cntrl_in[2] << 8, 2, 0;
      cntrl_in[3] << 4, 0, 0;

      bezier_type bc1(3);
      data_type length_cal, length_ref;

      // set control points
      for (index_type i=0; i<4; ++i)
      {
        bc1.set_control_point(cntrl_in[i], i);
      }

      // calculate the reference length using simpson's rule
      size_t i, npts(1001), n(npts-1);
      std::vector<data__> speed(npts);
      for (i=0; i<npts; ++i)
      {
        data_type t(i/static_cast<data__>(n));
        speed[i]=bc1.fp(t).norm();
      }
      length_ref=speed[0];
      for (i=1; i<=n-1; i+=2)
      {
        length_ref+=4*speed[i];
        length_ref+=2*speed[i+1];
      }
      length_ref-=speed[n];
      length_ref*=static_cast<data__>(1-0)/n/3;

      // compute length and compare
      data_type tol(std::sqrt(eps));
      length(length_cal, bc1, tol);
      TEST_ASSERT_DELTA(1, length_cal/length_ref, tol);

      // test computing some segment length
      data_type length01_cal, length12_cal, t0, t1, t2;
      t0=0;
      t1=static_cast<data__>(0.3);
      t2=1;

      length(length01_cal, bc1, t0, t1, tol);
      length(length12_cal, bc1, t1, t2, tol);
      TEST_ASSERT_DELTA(1, (length01_cal+length12_cal)/length_cal, tol);
    }

    void math_test()
    {
      data_type eps(std::numeric_limits<data__>::epsilon());

      point_type cntrl_in[7];

      cntrl_in[0] <<  0, 0, 0;
      cntrl_in[1] <<  2, 6, 0;
      cntrl_in[2] <<  3, 0, 0;
      cntrl_in[3] <<  5, 4, 0;
      cntrl_in[4] <<  7, 1, 0;
      cntrl_in[5] <<  5, 5, 0;
      cntrl_in[6] << 10, 6, 0;

      bezier_type bc1(6);

      // set control points
      for (index_type i=0; i<7; ++i)
      {
        bc1.set_control_point(cntrl_in[i], i);
      }

      point_type  cntrl_in2[4];

      // set control points
      cntrl_in2[0] << 0, 0, 0;
      cntrl_in2[1] << 0, 2, 0;
      cntrl_in2[2] << 8, 2, 0;
      cntrl_in2[3] << 4, 0, 0;

      bezier_type bc2(3);

      // set control points
      for (index_type i=0; i<4; ++i)
      {
        bc2.set_control_point(cntrl_in2[i], i);
      }


      bezier_type bc3;

      bc3.product( bc1, bc2 );

      point_type eval_out, eval_ref, a, b;
      data_type x_ref, t;

      t = 0.5;

      eval_out = bc3.f( t );
      eval_ref = bc1.f(t).cwiseProduct( bc2.f(t) );

      TEST_ASSERT( (eval_out-eval_ref).norm() < 5*eps );

      size_t i, n(20);

      for ( i = 0; i < n; i++ )
      {
        t = i * 1.0 / (n-1);
        eval_out = bc3.f( t );
        eval_ref = bc1.f(t).cwiseProduct( bc2.f(t) );
        TEST_ASSERT( (eval_out-eval_ref).norm() < 100*eps );
      }

      bezier_type delta, deltasq;

      oned_type rsq, r;
      typename oned_type::point_type x;

      bezier_type bc2neg( bc2 );
      bc2neg.scale( -1.0 );
      delta.sum( bc1, bc2neg );
      deltasq.square( delta );
      rsq = deltasq.sumcompcurve();


      t = 0.5;

      eval_out = delta.f( t );
      eval_ref = bc1.f(t) - bc2.f(t);
      TEST_ASSERT( (eval_out-eval_ref).norm() < 10*eps );

      eval_out = deltasq.f( t );
      eval_ref = eval_ref.cwiseProduct( eval_ref );
      TEST_ASSERT( (eval_out-eval_ref).norm() < 10*eps );

      x = rsq.f( t );
      x_ref = eval_ref(0) + eval_ref(1) + eval_ref(2);
      TEST_ASSERT( (x(0)-x_ref) < 10*eps );

      std::cout << std::endl;

      for ( i = 0; i < n; i++ )
      {
         t = i * 1.0 / (n-1);

         eval_out = delta.f( t );
         eval_ref = bc1.f(t) - bc2.f(t);
         TEST_ASSERT( (eval_out-eval_ref).norm() < 100*eps );

         eval_out = deltasq.f( t );
         eval_ref = eval_ref.cwiseProduct( eval_ref );
         TEST_ASSERT( (eval_out-eval_ref).norm() < 500*eps );

         x = rsq.f( t );
         x_ref = eval_ref(0) + eval_ref(1) + eval_ref(2);
         TEST_ASSERT( (x(0)-x_ref) < 100*eps );
      }
      std::cout << std::endl;
      std::cout << std::endl;
      std::cout << std::endl;

    }

    void integral_test()
    {
      point_type cntrl_in[4];
      data_type eps(std::numeric_limits<data__>::epsilon());

      // set control points
      cntrl_in[0] << 2.0, 2.0, 0.0;
      cntrl_in[1] << 1.0, 1.5, 0.0;
      cntrl_in[2] << 3.5, 0.0, 0.0;
      cntrl_in[3] << 4.0, 1.0, 0.0;

      bezier_type bc(3);
      bezier_type bi, bip;
      bezier_type bp, bpi;

      // set control points
      for (index_type i=0; i<4; ++i)
      {
        bc.set_control_point(cntrl_in[i], i);
      }

      bc.fi( bi );  // Integrate bc to form bi.
      bi.fp( bip ); // Differentiate bi to form bip -- which should equal bc.

      bc.fp( bp );  // Differentiate bc to bp
      bp.fi( bpi ); // Integrate bp to form bpi -- which should equal bc - with zero starting point...
      bpi.translate( bc.get_control_point(0) );

      TEST_ASSERT(bc.degree()==bip.degree())

      TEST_ASSERT(bc.degree()==bpi.degree())

      for (index_type i=0; i<4; ++i)
      {
        TEST_ASSERT( bip.get_control_point(i)==bc.get_control_point(i) );
      }

      for (index_type i=0; i<4; ++i)
      {
        TEST_ASSERT( (bpi.get_control_point(i)-bc.get_control_point(i)).norm() < 1300*eps );
      }
    }

    // Test that demonstrates the steps in projecting a point to a curve.  Used for demonstrating
    // algebraic operators on Bezier curves.
    void demo_pt_project_test()
    {
      bool output = false;
      point_type cntrl_in[4];

      point_type x;

      x << -0.6, -0.25, -0.0;

      // set control points
      cntrl_in[0] << 0.0, 0.0, 0.0;
      cntrl_in[1] << 0.0, 1.0, 0.0;
      cntrl_in[2] << 1.4, 0.9, 0.0;
      cntrl_in[3] << 1.0, 0.0, 0.0;

      bezier_type bc(3);

      // set control points
      for (index_type i=0; i<4; ++i)
      {
        bc.set_control_point(cntrl_in[i], i);
      }

      if ( output )
        bc.octave_print( 1 );

      bc.translate( x );

      if ( output )
        bc.octave_print( 2 );

      bezier_type bcprod;

      bcprod.product( bc, bc );

      if ( output )
        bcprod.octave_print( 3 );

      oned_type rsq, dr;
      rsq = bcprod.sumcompcurve();

      if ( output )
        rsq.octave_print( 4 );

      rsq.fp( dr );

      if ( output )
        dr.octave_print( 5 );

      oned_type o1, o2;
      oned_type o11, o12, o21, o22;

      dr.split( o1, o2, 0.5 );
      o1.split( o11, o12, 0.5 );
      o2.split( o21, o22, 0.5 );

      if ( output )
      {
        o11.octave_print( 6 );
        o12.octave_print( 7 );
        o21.octave_print( 8 );
        o22.octave_print( 9 );

        o1.octave_print( 10 );
        o2.octave_print( 11 );
      }

    }

};

#endif

