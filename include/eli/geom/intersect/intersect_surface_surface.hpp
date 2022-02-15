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
        struct surf_surf_g_functor
        {
          const surface__ *s1;
          const surface__ *s2;
          typename surface__::data_type k;

          typedef typename Eigen::Matrix<typename surface__::data_type, 4, 1> vec;

          vec operator()(const vec &x) const
          {
            typename surface__::data_type u1(x[0]), v1(x[1]);
            typename surface__::data_type u2(x[2]), v2(x[3]);
            vec rtn;

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

            typename surface__::point_type p1, p2, pave, disp, tvec;

            p1=s1->f(u1,v1);
            p2=s2->f(u2,v2);

            pave=(p1+p2)*0.5;
            disp=p2-p1;

            tvec=(s1->f_u(u1,v1).cross(s1->f_v(u1,v1))).cross(s2->f_u(u2,v2).cross(s2->f_v(u2,v2)));

            rtn(0)=disp(0);
            rtn(1)=disp(1);
            rtn(2)=disp(2);
            rtn(3)=k*tvec.dot(pave);
            return rtn;
          }
        };

        template <typename surface__>
        struct surf_surf_gp_functor
        {
          const surface__ *s1;
          const surface__ *s2;
          typename surface__::data_type k;

          typedef typename Eigen::Matrix<typename surface__::data_type, 4, 1> vec;
          typedef typename Eigen::Matrix<typename surface__::data_type, 4, 4> mat;

          mat operator()(const vec &x) const
          {
            typename surface__::data_type u1(x[0]), v1(x[1]);
            typename surface__::data_type u2(x[2]), v2(x[3]);
            mat rtn;

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

            typename surface__::point_type S1u, S1v, S1uu, S1uv, S1vv;
            typename surface__::point_type S2u, S2v, S2uu, S2uv, S2vv;

            typename surface__::point_type p1, p2, pave, tvec, dist;

            p1=s1->f(u1,v1);
            p2=s2->f(u2,v2);

            pave=(p1+p2)*0.5;
            dist=pave;


            S1u=s1->f_u(u1, v1);
            S1v=s1->f_v(u1, v1);
            S1uu=s1->f_uu(u1, v1);
            S1uv=s1->f_uv(u1, v1);
            S1vv=s1->f_vv(u1, v1);

            S2u=s2->f_u(u2, v2);
            S2v=s2->f_v(u2, v2);
            S2uu=s2->f_uu(u2, v2);
            S2uv=s2->f_uv(u2, v2);
            S2vv=s2->f_vv(u2, v2);

            tvec=(S1u.cross(S1v)).cross(S2u.cross(S2v));


            rtn(0,0)=-S1u(0);
            rtn(0,1)=-S1v(0);
            rtn(0,2)=S2u(0);
            rtn(0,3)=S2v(0);

            rtn(1,0)=-S1u(1);
            rtn(1,1)=-S1v(1);
            rtn(1,2)=S2u(1);
            rtn(1,3)=S2v(1);

            rtn(2,0)=-S1u(2);
            rtn(2,1)=-S1v(2);
            rtn(2,2)=S2u(2);
            rtn(2,3)=S2v(2);

            // tvec.dot(dist);

            rtn(3,0)= k*(dist.dot( ( S1uu.cross(S1v)+S1u.cross(S1uv) ).cross(S2u.cross(S2v)) ) + tvec.dot( S1u * 0.5 ));
            rtn(3,1)= k*(dist.dot( ( S1uv.cross(S1v)+S1u.cross(S1vv) ).cross(S2u.cross(S2v)) ) + tvec.dot( S1v * 0.5 ));
            rtn(3,2)= k*(dist.dot( (S1u.cross(S1v)).cross( S2uu.cross(S2v)+S2u.cross(S2uv) ) ) + tvec.dot( S2u * 0.5 ));
            rtn(3,3)= k*(dist.dot( (S1u.cross(S1v)).cross( S2uv.cross(S2v)+S2u.cross(S2vv) ) ) + tvec.dot( S2v * 0.5 ));

            // TODO: What to do if matrix becomes singular?

            return rtn;
          }
        };
      }



      template<typename surface__>
      typename surface__::index_type intersect(typename surface__::data_type &u1, typename surface__::data_type &v1,
                                              typename surface__::data_type &u2, typename surface__::data_type &v2,
                                              typename surface__::data_type &dist,
                                              const surface__ &s1in, const surface__ &s2in, const typename surface__::point_type &pt,
                                              const typename surface__::data_type &u01, const typename surface__::data_type &v01,
                                              const typename surface__::data_type &u02, const typename surface__::data_type &v02 )
      {
        typedef eli::mutil::nls::newton_raphson_system_method<typename surface__::data_type, 4, 1> nonlinear_solver_type;
        nonlinear_solver_type nrm;
        internal::surf_surf_g_functor<surface__> g;
        internal::surf_surf_gp_functor<surface__> gp;
        typename surface__::data_type dist0;
        typename surface__::tolerance_type tol;

        typename surface__::data_type u1min, u1max, v1min, v1max;
        typename surface__::data_type u2min, u2max, v2min, v2max;

        surface__ s1 = s1in;
        surface__ s2 = s2in;

        // Shift surfaces to be centered at initial intersection point.  This forces all coordinates to be close to zero
        // thereby increasing available precision for the calculations.
        s1.translate( -pt );
        s2.translate( -pt );

        s1.get_parameter_min(u1min,v1min);
        s1.get_parameter_max(u1max,v1max);
        s2.get_parameter_min(u2min,v2min);
        s2.get_parameter_max(u2max,v2max);

        typename surface__::point_type p1, p2;

        p1=s1.f(u01,v01);
        p2=s2.f(u02,v02);

        // Relative importance of nearness to base point.
        typename surface__::data_type k = 1.0e-3;

        // setup the functors
        g.s1=&s1;
        g.s2=&s2;
        g.k=k;
        gp.s1=&s1;
        gp.s2=&s2;
        gp.k=k;

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
        typename surface__::index_type ret = nrm.find_root(ans, g, gp, rhs);

        if ( ret == nrm.converged )
        {
          u1=ans(0);
          v1=ans(1);
          u2=ans(2);
          v2=ans(3);

          dist = eli::geom::point::distance(s1.f(u1, v1), s2.f(u2,v2));

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
