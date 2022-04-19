/*********************************************************************************
* Copyright (c) 2022 Rob McDonald <rob.a.mcdonald@gmail.com>
*
* All rights reserved. This program and the accompanying materials
* are made available under the terms of the Eclipse Public License v1.0
* which accompanies this distribution, and is available at
* http://www.eclipse.org/legal/epl-v10.html
*
* Contributors:
*    Rob McDonald - initial code and implementation
********************************************************************************/

#ifndef eli_geom_intersect_find_rst_surface_hpp
#define eli_geom_intersect_find_rst_surface_hpp

#include <cmath>
#include <vector>
#include <list>
#include <algorithm>
#include <limits>

#include "eli/code_eli.hpp"

#include "eli/mutil/nls/iterative_system_root_base_constrained.hpp"
#include "eli/mutil/nls/newton_raphson_system_method.hpp"

#include "eli/geom/point/distance.hpp"

namespace eli
{
  namespace geom
  {
    namespace intersect
    {
      namespace internal
      {

        template <typename surface__>
        struct rst_g_functor
        {
          const surface__ *ps;
          typename surface__::point_type pt;
          typedef typename Eigen::Matrix<typename surface__::data_type, 3, 1> vec;

          vec operator()(const vec &x) const
          {
            typename surface__::data_type r( x(0) ), s( x(1) ), t( x(2) );
            vec rtn;

            typename surface__::data_type rmin( 0.0 ), rmax( 1.0 );
            typename surface__::data_type smin( 0.0 ), smax( 0.5 );
            typename surface__::data_type tmin( 0.0 ), tmax( 1.0 );

            if ( !(r>=rmin) )
            {
              std::cout << "rst_g_functor, r less than minimum.  r: " << r << " rmin: " << rmin << std::endl;
              r=rmin;
            }
            if ( !(r<=rmax) )
            {
              std::cout << "rst_g_functor, r greater than maximum.  r: " << r << " ramx: " << rmax << std::endl;
              r=rmax;
            }

            if ( !(s>=smin) )
            {
              std::cout << "rst_g_functor, s less than minimum.  s: " << s << " smin: " << smin << std::endl;
              s=smin;
            }
            if ( !(s<=smax) )
            {
              std::cout << "rst_g_functor, s greater than maximum.  s: " << s << " smax: " << smax << std::endl;
              s=smax;
            }

            if ( !(t>=tmin) )
            {
              std::cout << "rst_g_functor, t less than minimum.  t: " << t << " tmin: " << tmin << std::endl;
              t=tmin;
            }
            if ( !(t<=tmax) )
            {
              std::cout << "rst_g_functor, t greater than maximum.  t: " << t << " tmax: " << tmax << std::endl;
              t=tmax;
            }
            typename surface__::point_type tmp;

            tmp = ps->fRST( r, s, t ) - pt;
            rtn( 0 ) = tmp( 0 );
            rtn( 1 ) = tmp( 1 );
            rtn( 2 ) = tmp( 2 );

            return rtn;
          }
        };

        template <typename surface__>
        struct rst_gp_functor
        {
          const surface__ *ps;
          typename surface__::point_type pt;
          typedef typename Eigen::Matrix<typename surface__::data_type, 3, 1> vec;
          typedef typename Eigen::Matrix<typename surface__::data_type, 3, 3> mat;

          mat operator()(const vec &x) const
          {
            typename surface__::data_type r( x(0) ), s( x(1) ), t( x(2) );
            mat rtn;
            typename surface__::data_type rmin( 0.0 ), rmax( 1.0 );
            typename surface__::data_type smin( 0.0 ), smax( 0.5 );
            typename surface__::data_type tmin( 0.0 ), tmax( 1.0 );

            if ( !(r>=rmin) )
            {
              std::cout << "rst_g_functor, r less than minimum.  r: " << r << " rmin: " << rmin << std::endl;
              r=rmin;
            }
            if ( !(r<=rmax) )
            {
              std::cout << "rst_g_functor, r greater than maximum.  r: " << r << " ramx: " << rmax << std::endl;
              r=rmax;
            }

            if ( !(s>=smin) )
            {
              std::cout << "rst_g_functor, s less than minimum.  s: " << s << " smin: " << smin << std::endl;
              s=smin;
            }
            if ( !(s<=smax) )
            {
              std::cout << "rst_g_functor, s greater than maximum.  s: " << s << " smax: " << smax << std::endl;
              s=smax;
            }

            if ( !(t>=tmin) )
            {
              std::cout << "rst_g_functor, t less than minimum.  t: " << t << " tmin: " << tmin << std::endl;
              t=tmin;
            }
            if ( !(t<=tmax) )
            {
              std::cout << "rst_g_functor, t greater than maximum.  t: " << t << " tmax: " << tmax << std::endl;
              t=tmax;
            }

            typename surface__::point_type Sr, Ss, St;

            Sr = ps->f_R( r, s, t );
            Ss = ps->f_S( r, s, t );
            St = ps->f_T( r, s, t );

            rtn( 0, 0 ) = Sr( 0 );
            rtn( 0, 1 ) = Ss( 0 );
            rtn( 0, 2 ) = St( 0 );

            rtn( 1, 0 ) = Sr( 1 );
            rtn( 1, 1 ) = Ss( 1 );
            rtn( 1, 2 ) = St( 1 );

            rtn( 2, 0 ) = Sr( 2 );
            rtn( 2, 1 ) = Ss( 2 );
            rtn( 2, 2 ) = St( 2 );

            return rtn;
          }
        };
      }

      template<typename surface__>
      typename surface__::data_type find_rst( typename surface__::data_type &r, typename surface__::data_type &s, typename surface__::data_type &t,
                                              const surface__ &surf, const typename surface__::point_type &pt,
                                              const typename surface__::data_type &r0, const typename surface__::data_type &s0, const typename surface__::data_type &t0,
                                              int & ret )
      {
        typedef eli::mutil::nls::newton_raphson_system_method<typename surface__::data_type, 3, 1> nonlinear_solver_type;
        nonlinear_solver_type nrm;
        internal::rst_g_functor<surface__> g;
        internal::rst_gp_functor<surface__> gp;
        typename surface__::data_type dist0, dist;
        typename surface__::tolerance_type tol;

        typename surface__::data_type rmin( 0.0 ), rmax( 1.0 );
        typename surface__::data_type smin( 0.0 ), smax( 0.5 );
        typename surface__::data_type tmin( 0.0 ), tmax( 1.0 );

          // setup the functors
        g.ps = &surf;
        g.pt = pt;
        gp.ps = &surf;
        gp.pt = pt;

        // setup the solver
        nrm.set_absolute_f_tolerance( tol.get_absolute_tolerance() );
        nrm.set_max_iteration( 20 );
        nrm.set_norm_type( nonlinear_solver_type::max_norm );

        nrm.set_lower_condition( 0, rmin, nonlinear_solver_type::IRC_INCLUSIVE );
        nrm.set_upper_condition( 0, rmax, nonlinear_solver_type::IRC_INCLUSIVE );

        nrm.set_lower_condition( 1, smin, nonlinear_solver_type::IRC_INCLUSIVE );
        nrm.set_upper_condition( 1, smax, nonlinear_solver_type::IRC_INCLUSIVE );

        nrm.set_lower_condition( 2, tmin, nonlinear_solver_type::IRC_INCLUSIVE );
        nrm.set_upper_condition( 2, tmax, nonlinear_solver_type::IRC_INCLUSIVE );

        // set the initial guess
        typename nonlinear_solver_type::solution_matrix xinit, rhs, ans;
        typename nonlinear_solver_type::solution_matrix resid;
        typename nonlinear_solver_type::jacobian_matrix jacob;

        xinit(0) = r0;
        xinit(1) = s0;
        xinit(2) = t0;

//        resid = g( xinit );
//        std::cout << "Initial residual " << std::endl;
//        std::cout << resid << std::endl;
//
//        jacob = gp( xinit );
//        std::cout << "Initial Jacobian " << std::endl;
//        std::cout << jacob << std::endl;

        nrm.set_initial_guess( xinit );
        rhs.setZero();
        dist0 = eli::geom::point::distance( surf.fRST( r0, s0, t0 ), pt );

        // find the root
        ret = nrm.find_root( ans, g, gp, rhs );
        r = ans( 0 );
        s = ans( 1 );
        t = ans( 2 );

//         std::cout << "niter: " << nrm.get_iteration_count() << std::endl;

        // if root is within bounds and is closer than initial guess
        {
          assert( ( r >= rmin ) && (r <= rmax) );
          assert( ( s >= smin ) && (s <= smax) );
          assert( ( t >= tmin ) && (t <= tmax) );

          dist = eli::geom::point::distance( surf.fRST( r, s, t ), pt );

          // std::cout << r << " " << s << " " << t << std::endl;
          // std::cout << "dist0 " << dist0 << " dist " << dist << std::endl;

          if  ( dist <= dist0 )
          {
            return dist;
          }
        }
//         else
//         {
//             std::cout << "% not converged";
//             if (stat==nonlinear_solver_type::hit_constraint)
//               std::cout << " because hit constraint" << std::endl;
//             else if (stat==nonlinear_solver_type::max_iteration)
//               std::cout << " reached max iteration" << std::endl;
//             else
//               std::cout << " for out of range parameters (" << ans(0) << ", " << ans(1) << ")" << std::endl;
//         }

        // couldn't find better answer so return initial guess
        r = r0;
        s = s0;
        t = t0;
        return dist0;
      }

    }
  }
}
#endif
