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

#ifndef eli_geom_intersect_intersect_surface_surface_hpp
#define eli_geom_intersect_intersect_surface_surface_hpp

#include <cmath>
#include <vector>
#include <list>
#include <algorithm>

#include "eli/code_eli.hpp"

#include "eli/mutil/nls/newton_raphson_method.hpp"
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
        template <typename surface__>
        struct surf_surf_g_gp_functor
        {
          const surface__ *s1;
          const surface__ *s2;
          typename surface__::point_type p0;

          // Unit direction of the intersection curve at the initial guess.  The solution is held
          // to the plane through p0 normal to it, which keeps the fourth equation in units of
          // length, like the other three, and linear in the unknowns.
          typename surface__::point_type t0;

          typedef typename Eigen::Matrix<typename surface__::data_type, 4, 1> vec;
          typedef typename Eigen::Matrix<typename surface__::data_type, 4, 4> mat;

          void operator()(vec &g, mat &gp, const vec &x) const
          {
            typename surface__::data_type u1(x[0]), v1(x[1]);
            typename surface__::data_type u2(x[2]), v2(x[3]);

            typename surface__::data_type u1min, u1max, v1min, v1max;
            typename surface__::data_type u2min, u2max, v2min, v2max;

            s1->get_parameter_min(u1min,v1min);
            s1->get_parameter_max(u1max,v1max);
            s2->get_parameter_min(u2min,v2min);
            s2->get_parameter_max(u2max,v2max);

            u1=std::min(std::max(u1, u1min), u1max);
            v1=std::min(std::max(v1, v1min), v1max);
            u2=std::min(std::max(u2, u2min), u2max);
            v2=std::min(std::max(v2, v2min), v2max);

            typename surface__::point_type p1, p2, pave, disp;

            typename surface__::point_type S1u, S1v;
            typename surface__::point_type S2u, S2v;

            p1=s1->f(u1,v1,p0);
            p2=s2->f(u2,v2,p0);

            pave=(p1+p2)*0.5;
            disp=p2-p1;

            S1u=s1->f_u(u1, v1);
            S1v=s1->f_v(u1, v1);

            S2u=s2->f_u(u2, v2);
            S2v=s2->f_v(u2, v2);

            g(0)=disp(0);
            g(1)=disp(1);
            g(2)=disp(2);
            g(3)=t0.dot(pave);

            gp(0,0)=-S1u(0);
            gp(0,1)=-S1v(0);
            gp(0,2)=S2u(0);
            gp(0,3)=S2v(0);

            gp(1,0)=-S1u(1);
            gp(1,1)=-S1v(1);
            gp(1,2)=S2u(1);
            gp(1,3)=S2v(1);

            gp(2,0)=-S1u(2);
            gp(2,1)=-S1v(2);
            gp(2,2)=S2u(2);
            gp(2,3)=S2v(2);

            gp(3,0)=t0.dot( S1u * 0.5 );
            gp(3,1)=t0.dot( S1v * 0.5 );
            gp(3,2)=t0.dot( S2u * 0.5 );
            gp(3,3)=t0.dot( S2v * 0.5 );

            // TODO: What to do if matrix becomes singular?
          }
        };

        // The same point on both surfaces with one of the four parameters held at a value: where
        // the intersection crosses that parameter's bound.  The held parameter's column is left out.
        template <typename surface__>
        struct surf_surf_held_g_gp_functor
        {
          const surface__ *s1;
          const surface__ *s2;
          typename surface__::point_type p0;
          int held;
          typename surface__::data_type held_val;

          typedef typename Eigen::Matrix<typename surface__::data_type, 3, 1> vec;
          typedef typename Eigen::Matrix<typename surface__::data_type, 3, 3> mat;

          void full(typename surface__::data_type x4[4], const vec &y) const
          {
            int j = 0;
            for (int i=0; i<4; ++i)
            {
              if (i==held)
              {
                x4[i]=held_val;
              }
              else
              {
                x4[i]=y[j];
                ++j;
              }
            }
          }

          void operator()(vec &g, mat &gp, const vec &y) const
          {
            typename surface__::data_type x4[4];
            full(x4, y);

            typename surface__::data_type u1min, u1max, v1min, v1max;
            typename surface__::data_type u2min, u2max, v2min, v2max;

            s1->get_parameter_min(u1min,v1min);
            s1->get_parameter_max(u1max,v1max);
            s2->get_parameter_min(u2min,v2min);
            s2->get_parameter_max(u2max,v2max);

            typename surface__::data_type u1(std::min(std::max(x4[0], u1min), u1max));
            typename surface__::data_type v1(std::min(std::max(x4[1], v1min), v1max));
            typename surface__::data_type u2(std::min(std::max(x4[2], u2min), u2max));
            typename surface__::data_type v2(std::min(std::max(x4[3], v2min), v2max));

            typename surface__::point_type disp, cols[4];

            disp=s2->f(u2,v2,p0)-s1->f(u1,v1,p0);

            cols[0]=-s1->f_u(u1, v1);
            cols[1]=-s1->f_v(u1, v1);
            cols[2]=s2->f_u(u2, v2);
            cols[3]=s2->f_v(u2, v2);

            g(0)=disp(0);
            g(1)=disp(1);
            g(2)=disp(2);

            int j = 0;
            for (int i=0; i<4; ++i)
            {
              if (i!=held)
              {
                gp(0,j)=cols[i](0);
                gp(1,j)=cols[i](1);
                gp(2,j)=cols[i](2);
                ++j;
              }
            }
          }
        };

        // Newton solve for the point on both surfaces with x[held] fixed, from the other three in
        // x.  Returns whether it converged, with the answer in x.
        template<typename surface__>
        bool intersect_held(typename surface__::data_type x[4], const surface__ &s1, const surface__ &s2,
                            const typename surface__::point_type &pt, int held)
        {
          typedef eli::mutil::nls::newton_raphson_system_method<typename surface__::data_type, 3, 1> held_solver_type;
          held_solver_type hnrm;
          surf_surf_held_g_gp_functor<surface__> hggp;
          typename surface__::tolerance_type tol;

          typename surface__::data_type pmin[4], pmax[4];
          s1.get_parameter_min(pmin[0],pmin[1]);
          s1.get_parameter_max(pmax[0],pmax[1]);
          s2.get_parameter_min(pmin[2],pmin[3]);
          s2.get_parameter_max(pmax[2],pmax[3]);

          bool open[4] = { s1.open_u(), s1.open_v(), s2.open_u(), s2.open_v() };

          hggp.s1=&s1;
          hggp.s2=&s2;
          hggp.p0=pt;
          hggp.held=held;
          hggp.held_val=x[held];

          hnrm.set_absolute_f_tolerance(tol.get_absolute_tolerance());
          hnrm.set_max_iteration(20);
          hnrm.set_norm_type(held_solver_type::max_norm);

          typename held_solver_type::solution_matrix yinit, hrhs, yans;
          int j = 0;
          for (int i=0; i<4; ++i)
          {
            if (i!=held)
            {
              if (open[i])
              {
                hnrm.set_lower_condition(j, pmin[i], held_solver_type::IRC_EXCLUSIVE);
                hnrm.set_upper_condition(j, pmax[i], held_solver_type::IRC_EXCLUSIVE);
              }
              else
              {
                hnrm.set_periodic_condition(j, pmin[i], pmax[i]);
              }
              yinit(j)=x[i];
              ++j;
            }
          }
          hnrm.set_initial_guess(yinit);
          hrhs.setZero();

          if ( hnrm.find_root(yans, hggp, hrhs) != hnrm.converged )
          {
            return false;
          }

          hggp.full(x, yans);
          return true;
        }
      }

      // The point on both surfaces with one of the four parameters, u1 v1 u2 v2 in that order,
      // held at the value given for it: where the intersection crosses that parameter line.
      // Returns 0 where it converges, with the answer and the distance between the surfaces
      // there; otherwise the guess is returned.
      template<typename surface__>
      typename surface__::index_type intersect(typename surface__::data_type &u1, typename surface__::data_type &v1,
                                              typename surface__::data_type &u2, typename surface__::data_type &v2,
                                              typename surface__::data_type &dist,
                                              const surface__ &s1, const surface__ &s2, const typename surface__::point_type &pt,
                                              const typename surface__::data_type &u01, const typename surface__::data_type &v01,
                                              const typename surface__::data_type &u02, const typename surface__::data_type &v02,
                                              int held )
      {
        typename surface__::data_type x[4] = { u01, v01, u02, v02 };

        if ( internal::intersect_held( x, s1, s2, pt, held ) )
        {
          u1=x[0];
          v1=x[1];
          u2=x[2];
          v2=x[3];
          dist = eli::geom::point::distance(s1.f(u1, v1, pt), s2.f(u2, v2, pt));
          return 0;
        }

        u1=u01;
        v1=v01;
        u2=u02;
        v2=v02;
        dist = eli::geom::point::distance(s1.f(u1, v1, pt), s2.f(u2, v2, pt));
        return 1;
      }

      template<typename surface__>
      typename surface__::index_type intersect(typename surface__::data_type &u1, typename surface__::data_type &v1,
                                              typename surface__::data_type &u2, typename surface__::data_type &v2,
                                              typename surface__::data_type &dist,
                                              const surface__ &s1, const surface__ &s2, const typename surface__::point_type &pt,
                                              const typename surface__::data_type &u01, const typename surface__::data_type &v01,
                                              const typename surface__::data_type &u02, const typename surface__::data_type &v02 )
      {
        typedef eli::mutil::nls::newton_raphson_system_method<typename surface__::data_type, 4, 1> nonlinear_solver_type;
        nonlinear_solver_type nrm;
        internal::surf_surf_g_gp_functor<surface__> ggp;
        typename surface__::data_type dist0;
        typename surface__::tolerance_type tol;

        typename surface__::data_type u1min, u1max, v1min, v1max;
        typename surface__::data_type u2min, u2max, v2min, v2max;

        s1.get_parameter_min(u1min,v1min);
        s1.get_parameter_max(u1max,v1max);
        s2.get_parameter_min(u2min,v2min);
        s2.get_parameter_max(u2max,v2max);

        typename surface__::point_type p1, p2;

        // Use offset function evaluation here and in functors to shift surfaces to be centered near initial intersection point.
        // This forces all coordinates to be close to zero thereby increasing available precision for the calculations.
        p1=s1.f(u01,v01,pt);
        p2=s2.f(u02,v02,pt);

        // The direction the intersection runs at the guess.  Where the surfaces are tangent
        // there is none, and any plane through the guess will do.
        typename surface__::point_type tvec;
        tvec=(s1.f_u(u01, v01).cross(s1.f_v(u01, v01))).cross(s2.f_u(u02, v02).cross(s2.f_v(u02, v02)));
        if (tvec.norm()>0)
        {
          tvec.normalize();
        }
        else
        {
          tvec=s1.f_u(u01, v01);
          if (tvec.norm()>0)
          {
            tvec.normalize();
          }
        }

        // setup the functors
        ggp.s1=&s1;
        ggp.s2=&s2;
        ggp.p0=pt;
        ggp.t0=tvec;

        // setup the solver
        nrm.set_absolute_f_tolerance(tol.get_absolute_tolerance());
        nrm.set_max_iteration(20);
        nrm.set_norm_type(nonlinear_solver_type::max_norm);

        if (s1.open_u())
        {
          nrm.set_lower_condition(0, u1min, nonlinear_solver_type::IRC_EXCLUSIVE);
          nrm.set_upper_condition(0, u1max, nonlinear_solver_type::IRC_EXCLUSIVE);
        }
        else
        {
          nrm.set_periodic_condition(0, u1min, u1max);
        }

        if (s1.open_v())
        {
          nrm.set_lower_condition(1, v1min, nonlinear_solver_type::IRC_EXCLUSIVE);
          nrm.set_upper_condition(1, v1max, nonlinear_solver_type::IRC_EXCLUSIVE);
        }
        else
        {
          nrm.set_periodic_condition(1, v1min, v1max);
        }

        if (s2.open_u())
        {
          nrm.set_lower_condition(2, u2min, nonlinear_solver_type::IRC_EXCLUSIVE);
          nrm.set_upper_condition(2, u2max, nonlinear_solver_type::IRC_EXCLUSIVE);
        }
        else
        {
          nrm.set_periodic_condition(2, u2min, u2max);
        }

        if (s2.open_v())
        {
          nrm.set_lower_condition(3, v2min, nonlinear_solver_type::IRC_EXCLUSIVE);
          nrm.set_upper_condition(3, v2max, nonlinear_solver_type::IRC_EXCLUSIVE);
        }
        else
        {
          nrm.set_periodic_condition(3, v2min, v2max);
        }

        // set the initial guess
        typename nonlinear_solver_type::solution_matrix xinit, rhs, ans;

        xinit(0)=u01;
        xinit(1)=v01;
        xinit(2)=u02;
        xinit(3)=v02;
        nrm.set_initial_guess(xinit);
        rhs.setZero();

        dist0=eli::geom::point::distance(p1, p2);

        // find the root
        typename surface__::index_type ret = nrm.find_root(ans, ggp, rhs);

        // Stopped against a bound: the intersection crosses that edge of a surface.  The clamp
        // leaves the parameter exactly on the bound, where it is held while the other three put
        // the point on both surfaces.  Where the edge lies along the other surface there is no
        // one crossing, and the held solve finds none.
        if ( ret == nrm.hit_constraint )
        {
          typename surface__::data_type pmin[4] = { u1min, v1min, u2min, v2min };
          typename surface__::data_type pmax[4] = { u1max, v1max, u2max, v2max };

          int held = -1;
          for (int i=0; i<4 && held<0; ++i)
          {
            if (ans(i)==pmin[i] || ans(i)==pmax[i])
            {
              held = i;
            }
          }

          if ( held >= 0 )
          {
            typename surface__::data_type x4[4] = { ans(0), ans(1), ans(2), ans(3) };
            if ( internal::intersect_held( x4, s1, s2, pt, held ) )
            {
              for (int i=0; i<4; ++i)
              {
                ans(i)=x4[i];
              }
              ret = nrm.converged;
            }
          }
        }

        if ( ret == nrm.converged )
        {
          u1=ans(0);
          v1=ans(1);
          u2=ans(2);
          v2=ans(3);

          dist = eli::geom::point::distance(s1.f(u1, v1, pt), s2.f(u2, v2, pt));

//        if( dist > 1e-6 )
//        {
//          printf("d0: %g d: %g\n", dist0, dist );
//          printf(" u01: %f u1: %f\n", u01, u1 );
//          printf(" v01: %f v1: %f\n", v01, v1 );
//          printf(" u02: %f u2: %f\n", u02, u2 );
//          printf(" v02: %f v2: %f\n", v02, v2 );
//        }

          if  (dist<=dist0)
          {
            return ret;
          }
          ret = 3; // Converged, but worse answer than initial guess.
        }

        // couldn't find better answer so return initial guess
        u1=u01;
        v1=v01;
        u2=u02;
        v2=v02;
        dist=dist0;
        return ret;
      }


    }
  }
}
#endif
