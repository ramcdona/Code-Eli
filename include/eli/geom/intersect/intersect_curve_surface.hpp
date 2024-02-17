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

#ifndef eli_geom_intersect_intersect_curve_surface_hpp
#define eli_geom_intersect_intersect_curve_surface_hpp

#include <cmath>
#include <vector>
#include <list>
#include <algorithm>

#include "eli/code_eli.hpp"

#include "eli/mutil/nls/newton_raphson_method.hpp"

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
        template <typename surface__, typename curve__>
        struct surf_curve_g_gp_functor
        {
          const surface__ *s;
          const curve__ *c;

          typedef typename Eigen::Matrix<typename surface__::data_type, 3, 1> vec;
          typedef typename Eigen::Matrix<typename surface__::data_type, 3, 3> mat;

          void operator()(vec &g, mat &gp, const vec &x) const
          {
            typename surface__::data_type u(x[0]), v(x[1]);
            typename surface__::data_type t(x[2]);

            typename surface__::data_type umin, umax, vmin, vmax;
            typename curve__::data_type tmin, tmax;

            s->get_parameter_min(umin,vmin);
            s->get_parameter_max(umax,vmax);
            tmin = c->get_t0();
            tmax = c->get_tmax();

            u=std::min(std::max(u, umin), umax);
            v=std::min(std::max(v, vmin), vmax);
            t=std::min(std::max(t, tmin), tmax);

            typename surface__::point_type disp;
            typename surface__::point_type Su, Sv;
            typename curve__::point_type Ct;

            disp = c->f(t) - s->f(u,v);

            g(0)=disp(0);
            g(1)=disp(1);
            g(2)=disp(2);

            Su=s->f_u(u, v);
            Sv=s->f_v(u, v);

            Ct=c->fp(t);

            gp(0,0)=-Su(0);
            gp(0,1)=-Sv(0);
            gp(0,2)=Ct(0);

            gp(1,0)=-Su(1);
            gp(1,1)=-Sv(1);
            gp(1,2)=Ct(1);

            gp(2,0)=-Su(2);
            gp(2,1)=-Sv(2);
            gp(2,2)=Ct(2);

            // TODO: What to do if matrix becomes singular?

          }
        };
      }



      template<typename surface__, typename curve__>
      typename surface__::data_type intersect(typename surface__::data_type &u, typename surface__::data_type &v,
                                              typename curve__::data_type &t,
                                              const surface__ &s, const curve__ &c,
                                              const typename surface__::data_type &u0, const typename surface__::data_type &v0,
                                              const typename curve__::data_type &t0 )
      {
        typedef eli::mutil::nls::newton_raphson_system_method<typename surface__::data_type, 3, 1> nonlinear_solver_type;
        nonlinear_solver_type nrm;
        internal::surf_curve_g_gp_functor<surface__, curve__> ggp;
        typename surface__::data_type dist0, dist;
        typename surface__::tolerance_type tol;

        typename surface__::data_type umin, umax, vmin, vmax;
        typename curve__::data_type tmin, tmax;

        s.get_parameter_min(umin,vmin);
        s.get_parameter_max(umax,vmax);
        tmin = c.get_t0();
        tmax = c.get_tmax();

        u=std::min(std::max(u, umin), umax);
        v=std::min(std::max(v, vmin), vmax);
        t=std::min(std::max(t, tmin), tmax);

        typename surface__::point_type p1, p2;

        p1=s.f(u0,v0);
        p2=c.f(t0);

        // setup the functors
        ggp.s=&s;
        ggp.c=&c;

        // setup the solver
        nrm.set_absolute_f_tolerance(tol.get_absolute_tolerance());
        nrm.set_max_iteration(20);
        nrm.set_norm_type(nonlinear_solver_type::max_norm);

        if (s.open_u())
        {
          nrm.set_lower_condition(0, umin, nonlinear_solver_type::IRC_EXCLUSIVE);
          nrm.set_upper_condition(0, umax, nonlinear_solver_type::IRC_EXCLUSIVE);
        }
        else
        {
          nrm.set_periodic_condition(0, umin, umax);
        }

        if (s.open_v())
        {
          nrm.set_lower_condition(1, vmin, nonlinear_solver_type::IRC_EXCLUSIVE);
          nrm.set_upper_condition(1, vmax, nonlinear_solver_type::IRC_EXCLUSIVE);
        }
        else
        {
          nrm.set_periodic_condition(1, vmin, vmax);
        }

        if (c.open())
        {
          nrm.set_lower_condition(2, tmin, nonlinear_solver_type::IRC_EXCLUSIVE);
          nrm.set_upper_condition(2, tmax, nonlinear_solver_type::IRC_EXCLUSIVE);
        }
        else
        {
          nrm.set_periodic_condition(2, tmin, tmax);
        }

        // set the initial guess
        typename nonlinear_solver_type::solution_matrix xinit, rhs, ans;

        xinit(0)=u0;
        xinit(1)=v0;
        xinit(2)=t0;

        nrm.set_initial_guess(xinit);
        rhs.setZero();

        dist0=eli::geom::point::distance(p1, p2);

        // find the root
        nrm.find_root(ans, ggp, rhs);
        u=ans(0);
        v=ans(1);
        t=ans(2);

        // if root is within bounds and is closer than initial guess
        assert((u>=umin) && (u<=umax));
        assert((v>=vmin) && (v<=vmax));
        assert((t>=tmin) && (t<=tmax));

        dist = eli::geom::point::distance(s.f(u, v), c.f(t));

        if  (dist<=dist0)
        {
          return dist;
        }

        // couldn't find better answer so return initial guess
        u=u0;
        v=v0;
        t=t0;

        return dist0;
      }


    }
  }
}
#endif
