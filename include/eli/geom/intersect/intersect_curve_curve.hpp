/*********************************************************************************
* Copyright (c) 2013 David D. Marshall <ddmarsha@calpoly.edu>
*
* All rights reserved. This program and the accompanying materials
* are made available under the terms of the Eclipse Public License v1.0
* which accompanies this distribution, and is available at
* http://www.eclipse.org/legal/epl-v10.html
*
* Contributors:
*    Rob McDonald - initial code and implementation
********************************************************************************/

#ifndef eli_geom_intersect_intersect_curve_curve_hpp
#define eli_geom_intersect_intersect_curve_curve_hpp

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
    namespace intersect
    {
      namespace internal
      {
        template <typename curve__>
        struct curve_curve_g_gp_functor
        {
          const curve__ *c1;
          const curve__ *c2;

          typedef typename Eigen::Matrix<typename curve__::data_type, 2, 1> vec;
          typedef typename Eigen::Matrix<typename curve__::data_type, 2, 2> mat;

          void operator()(vec &g, mat &gp, const vec &x) const
          {
            typename curve__::data_type t1(x[0]);
            typename curve__::data_type t2(x[1]);

            typename curve__::data_type tmin1, tmax1, tmin2, tmax2;

            tmin1 = c1->get_t0();
            tmax1 = c1->get_tmax();
            tmin2 = c2->get_t0();
            tmax2 = c2->get_tmax();

            t1 = std::min(std::max(t1, tmin1), tmax1);
            t2 = std::min(std::max(t2, tmin2), tmax2);

            typename curve__::point_type disp;
            typename curve__::point_type Ct1, Ct2;

            disp = c1->f(t1) - c2->f(t2);

            g(0) = disp(0);
            g(1) = disp(1);

            Ct1 = c1->fp(t1);
            Ct2 = c2->fp(t2);

            gp(0,0)=Ct1(0);
            gp(0,1)=-Ct2(0);

            gp(1,0)=Ct1(1);
            gp(1,1)=-Ct2(1);

            // TODO: What to do if matrix becomes singular?

          }
        };
      }



      template<typename curve__>
      typename curve__::data_type intersect(typename curve__::data_type &t1,
                                            typename curve__::data_type &t2,
                                            const curve__ &c1, const curve__ &c2,
                                            const typename curve__::data_type &t01,
                                            const typename curve__::data_type &t02,
                                            int & ret)
      {
        typedef eli::mutil::nls::newton_raphson_system_method<typename curve__::data_type, 2, 1> nonlinear_solver_type;
        nonlinear_solver_type nrm;
        internal::curve_curve_g_gp_functor<curve__> ggp;
        typename curve__::data_type dist0, dist;
        typename curve__::tolerance_type tol;

        typename curve__::data_type tmin1, tmax1, tmin2, tmax2;

        tmin1 = c1.get_t0();
        tmax1 = c1.get_tmax();
        tmin2 = c2.get_t0();
        tmax2 = c2.get_tmax();

        t1 = std::min(std::max(t1, tmin1), tmax1);
        t2 = std::min(std::max(t2, tmin2), tmax2);

        typename curve__::point_type p1, p2;

        p1 = c1.f(t01);
        p2 = c2.f(t02);

        // setup the functors
        ggp.c1 = &c1;
        ggp.c2 = &c2;

        // setup the solver
        nrm.set_absolute_f_tolerance(tol.get_absolute_tolerance());
        nrm.set_max_iteration(20);
        nrm.set_norm_type(nonlinear_solver_type::max_norm);

        if (c1.open())
        {
          nrm.set_lower_condition(0, tmin1, nonlinear_solver_type::IRC_INCLUSIVE);
          nrm.set_upper_condition(0, tmax1, nonlinear_solver_type::IRC_INCLUSIVE);
        }
        else
        {
          nrm.set_periodic_condition(0, tmin1, tmax1);
        }

        if (c2.open())
        {
          nrm.set_lower_condition(1, tmin2, nonlinear_solver_type::IRC_INCLUSIVE);
          nrm.set_upper_condition(1, tmax2, nonlinear_solver_type::IRC_INCLUSIVE);
        }
        else
        {
          nrm.set_periodic_condition(1, tmin2, tmax2);
        }

        // set the initial guess
        typename nonlinear_solver_type::solution_matrix xinit, rhs, ans;

        xinit(0) = t01;
        xinit(1) = t02;

        nrm.set_initial_guess(xinit);
        rhs.setZero();

        dist0=eli::geom::point::distance(p1, p2);

        // find the root
        ret = nrm.find_root(ans, ggp, rhs);
        t1 = ans(0);
        t2 = ans(1);

        // if root is within bounds and is closer than initial guess
        assert((t1 >= tmin1) && (t1 <= tmax1));
        assert((t2 >= tmin2) && (t2 <= tmax2));

        dist = eli::geom::point::distance(c1.f(t1), c2.f(t2));

        if  (dist<=dist0)
        {
          return dist;
        }

        // couldn't find better answer so return initial guess
        t1 = t01;
        t2 = t02;

        return dist0;
      }


    }
  }
}
#endif
