/*********************************************************************************
* Copyright (c) 2015 Rob McDonald <ramcdona@calpoly.edu>
*
* All rights reserved. This program and the accompanying materials
* are made available under the terms of the Eclipse Public License v1.0
* which accompanies this distribution, and is available at
* http://www.eclipse.org/legal/epl-v10.html
*
* Contributors:
*    Rob McDonald - initial code and implementation
********************************************************************************/

#ifndef eli_geom_intersect_specified_thickness_curve_hpp
#define eli_geom_intersect_specified_thickness_curve_hpp

#include <cmath>
#include <vector>
#include <list>
#include <algorithm>

#include "eli/code_eli.hpp"

#include "eli/mutil/nls/newton_raphson_system_method.hpp"

#include "eli/geom/point/distance.hpp"
#include "eli/geom/curve/piecewise.hpp"
#include "eli/geom/intersect/minimum_distance_bounding_box.hpp"

namespace eli
{
  namespace geom
  {
    namespace curve
    {
      template<template<typename, unsigned short, typename> class curve__, typename data__, unsigned short dim__, typename tol__ >
      class piecewise;
    }

    namespace intersect
    {
      namespace internal
      {
        template <typename curve__>
        struct curve_thick_g_gp_functor
        {
          const curve__ *pc;
          typename curve__::point_type pt;
          typename curve__::data_type thick;

          typedef typename Eigen::Matrix<typename curve__::data_type, 2, 1> vec;
          typedef typename Eigen::Matrix<typename curve__::data_type, 2, 2> mat;

          void operator()(vec &g, mat &gp, const vec &ts) const
          {
            typename curve__::data_type t1(ts[0]), t2(ts[1]);

            if ( !(t1>=pc->get_t0()) )
            {
              std::cout << "Specified thickness curve g_functor, t1 less than minimum.  t1: " << t1 << " t0: " << pc->get_t0() << std::endl;
              t1=pc->get_t0();
            }
            if ( !(t1<=pc->get_tmax()) )
            {
              std::cout << "Specified thickness curve g_functor, t1 greater than maximum.  t1: " << t1 << " tmax: " << pc->get_tmax() << std::endl;
              t1=pc->get_tmax();
            }

            if ( !(t2>=pc->get_t0()) )
            {
              std::cout << "Specified thickness curve g_functor, t2 less than minimum.  t2: " << t2 << " t0: " << pc->get_t0() << std::endl;
              t2=pc->get_t0();
            }
            if ( !(t2<=pc->get_tmax()) )
            {
              std::cout << "Specified thickness curve g_functor, t2 greater than maximum.  t2: " << t2 << " tmax: " << pc->get_tmax() << std::endl;
              t2=pc->get_tmax();
            }

            assert((t1>=pc->get_t0()) && (t1<=pc->get_tmax()));
            assert((t2>=pc->get_t0()) && (t2<=pc->get_tmax()));

            typename curve__::point_type u = pc->f(t1) - pt;
            typename curve__::point_type v = pc->f(t2) - pt;
            typename curve__::point_type w = u - v;

            g(0) = u.dot(u) - v.dot(v);
            g(1) = w.dot(w) - thick * thick;

            typename curve__::point_type du = pc->fp(t1);
            typename curve__::point_type dv = pc->fp(t2);

            gp(0,0) = 2.0 * u.dot(du);
            gp(0,1) = - 2.0 * v.dot(dv);
            gp(1,0) = 2.0 * w.dot(du);
            gp(1,1) = - 2.0 * w.dot(dv);
          }
        };
      }

      template<typename curve__>
      typename curve__::data_type specified_thickness(typename curve__::data_type &t1, typename curve__::data_type &t2, const curve__ &c, const typename curve__::point_type &pt, const typename curve__::data_type &d, const typename curve__::data_type &t10, const typename curve__::data_type &t20)
      {
        typedef eli::mutil::nls::newton_raphson_system_method<typename curve__::data_type, 2, 1> nonlinear_solver_type;
        nonlinear_solver_type nrm;
        internal::curve_thick_g_gp_functor<curve__> ggp;
        typename curve__::data_type r1, r2, dist, dist0;
        typename curve__::tolerance_type tol;

        // setup the functors
        ggp.pc=&c;
        ggp.pt=pt;
        ggp.thick=d;

        typename curve__::data_type tmin(c.get_t0());
        typename curve__::data_type tmax(c.get_tmax());
        typename curve__::data_type tmid((tmin+tmax)/2.0);

        // setup the solver
        nrm.set_absolute_f_tolerance(tol.get_absolute_tolerance());
        nrm.set_max_iteration(10);

        nrm.set_lower_condition(0,tmin, nonlinear_solver_type::IRC_EXCLUSIVE);
        nrm.set_upper_condition(0,tmid, nonlinear_solver_type::IRC_EXCLUSIVE);
        nrm.set_lower_condition(1,tmid, nonlinear_solver_type::IRC_EXCLUSIVE);
        nrm.set_upper_condition(1,tmax, nonlinear_solver_type::IRC_EXCLUSIVE);

        dist0=eli::geom::point::distance(c.f(t10), c.f(t20))-d;

        typename nonlinear_solver_type::solution_matrix uinit, rhs, ans;

        // set the initial guess
        uinit(0)=t10;
        uinit(1)=t20;
        nrm.set_initial_guess(uinit);

        rhs.setZero();

        // find the root
        nrm.find_root(ans, ggp, rhs);
        t1=ans(0);
        t2=ans(1);

        // if root is within bounds and is closer than initial guess
        {
          assert((t1>=tmin) && (t1<=tmax));
          assert((t2>=tmin) && (t2<=tmax));

          dist=eli::geom::point::distance(c.f(t1), c.f(t2))-d;
          if ( std::abs(dist) <= std::abs(dist0) )
          {
            return dist;
          }
        }

        // couldn't find better answer so return initial guess
        t1=t10;
        t2=t20;
        return dist0;
      }

    }
  }
}
#endif
