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

#ifndef eli_geom_intersect_minimum_distance_surface_hpp
#define eli_geom_intersect_minimum_distance_surface_hpp

#include <cmath>
#include <iostream>
#include <vector>
#include <list>
#include <algorithm>
#include <limits>

#include "eli/code_eli.hpp"

#include "eli/mutil/nls/iterative_system_root_base_constrained.hpp"
#include "eli/mutil/nls/newton_raphson_system_method.hpp"

#include "eli/geom/intersect/minimum_distance_curve.hpp"
#include "eli/geom/point/distance.hpp"
#include "eli/geom/intersect/minimum_distance_bounding_box.hpp"
#include "eli/geom/intersect/findnonpos.hpp"
#include "eli/util/clamp.hpp"

namespace eli
{
  namespace geom
  {
    namespace intersect
    {
      namespace internal
      {

        template<typename surface__, size_t N__, size_t NSOL__=1>
        class tangent_plane_method : public mutil::nls::iterative_system_root_base_constrained<typename surface__::data_type, N__, NSOL__>
        {
          public:
            static const int hit_constraint = 101;
            const surface__ *s;
            typename surface__::point_type pt;
            typename surface__::data_type xtol;
            typename surface__::index_type maxit;

          public:
            tangent_plane_method()
            : mutil::nls::iterative_system_root_base_constrained<typename surface__::data_type, N__, NSOL__>()
            {
              x0.setConstant(static_cast<typename surface__::data_type>(0));
            }

            tangent_plane_method(const tangent_plane_method<surface__, N__, NSOL__> &tpm)
            : mutil::nls::iterative_system_root_base_constrained<typename surface__::data_type, N__, NSOL__>(tpm), x0(tpm.x0)
            {
            }

            ~tangent_plane_method()
            {
            }

            void set_initial_guess(const typename mutil::nls::iterative_system_root_base<typename surface__::data_type, N__, NSOL__>::solution_matrix &xg)
            {
              x0=xg;
            }

            const typename mutil::nls::iterative_system_root_base<typename surface__::data_type, N__, NSOL__>::solution_matrix & get_initial_guess() const
            {
              return x0;
            }

            int find_root(typename mutil::nls::iterative_system_root_base<typename surface__::data_type, N__, NSOL__>::solution_matrix &root) const
            {
              typename mutil::nls::iterative_system_root_base<typename surface__::data_type, N__, NSOL__>::solution_matrix dx, x(x0);
              typename surface__::data_type abs_x_norm, prev_dx;
              typename surface__::index_type count;

              typename surface__::point_type q;
              typename surface__::point_type Su, Sv, r, A, B, norm;

              typename surface__::data_type umin, umax, vmin, vmax;
              s->get_parameter_min(umin,vmin);
              s->get_parameter_max(umax,vmax);

              bool divflag = false;

              // Divergence is judged on the distance, which is what the iteration is trying to
              // reduce, not on the length of the step.  The step norm is the wrong quantity: it is
              // a max norm over both parameters, so as one converges the norm becomes whichever
              // component is still moving, and a drift of parts in 1e7 in that component reads as
              // two growing steps and the whole iteration is abandoned.  Traced on this surface,
              // that discarded a converged answer at d=0.607 and returned the seed at d=1.969 --
              // and of 1670 steps traced across failing cases, exactly one iteration ever reached
              // the converged return at all.
              typename surface__::data_type dprev( std::numeric_limits<typename surface__::data_type>::max() );
              typename surface__::data_type dnow( dprev );

              // Best point visited.  On divergence the iteration used to hand back the starting
              // guess, discarding everything it had found: one trace reaches d=0.948 and then
              // returns the seed at d=2.453.  The test for divergence is now right; the response
              // to it was still wrong.
              typename mutil::nls::iterative_system_root_base<typename surface__::data_type, N__, NSOL__>::solution_matrix xbest( x0 );
              typename surface__::data_type dbest( std::numeric_limits<typename surface__::data_type>::max() );

              // Trust region radius, as a fraction of the parameter span.  See the clamp below.
              typename surface__::data_type trust( static_cast<typename surface__::data_type>(0.25) );

              bool all_zero = false;

              abs_x_norm = std::numeric_limits<typename surface__::data_type>::max();
              prev_dx = abs_x_norm;
              count = 0;
              while ( count < maxit && abs_x_norm > xtol && !all_zero)
              {
                s->f_pt_derivs( x(0), x(1), q, Su, Sv );

                r = q - pt;

                dprev = dnow;
                dnow = r.norm();

                // Widen the trust region while the steps are working and shrink it when they are
                // not.  A fixed clamp is what produced the limit cycle this replaces: the raw
                // steps were +1.26 and -1.89, both longer than the limit, so both came out at
                // exactly the limit and v ping-ponged between 1.175 and 2.175 for ever.  The
                // divergence test could not see it either, since it wants two worsening steps in a
                // row and a two-cycle alternates.  Halving on a worse step turns the cycle into a
                // spiral that closes on the answer.
                //
                // It only ever shrinks.  Restoring the radius after a good step looks like the
                // usual trust region rule and is wrong here: a wing gave a four-cycle in which
                // every other step improved, so the radius doubled back up as fast as it came
                // down and the cycle survived.  There is nothing to restore in any case -- a
                // radius that is not binding costs nothing, so the only radius worth having is
                // one large enough to have been binding once.
                //
                // This costs nothing.  The radius is updated from the evaluation the iteration was
                // going to make anyway, and a step already shorter than the radius is untouched --
                // so a run that never gets clamped never notices any of it.
                if ( dnow > dprev )
                {
                  trust /= 2;
                }

                if ( dnow < dbest )
                {
                  dbest = dnow;
                  xbest = x;
                }

                A = Sv.cross( r );
                B = Su.cross( r );
                norm = Su.cross( Sv );

                typename surface__::data_type N = norm.dot( norm );

                if( std::abs( N ) > std::numeric_limits< typename surface__::data_type >::min() )
                {
                    dx(0) = A.dot( norm ) / N;
                    dx(1) = -B.dot( norm ) / N;
                }
                else
                {
                    dx(0) = 0.0;
                    dx(1) = 0.0;
                }

                // Keep the step inside the range the tangent plane can speak for.  The step is the
                // solution of a linearisation about the current point; asking it to move most of
                // the way round a closed surface is asking it about geometry it knows nothing of.
                // 26% of steps traced here wanted more than half the v period, and the periodic
                // wrap then turns such a step into a small one in the opposite direction -- or
                // into nothing at all, which reads as a converged iteration sitting at a point
                // that is not a solution.  A quarter of the span is well inside where a tangent
                // plane means anything, and clamping is not a line search: no extra evaluation,
                // and a step already short enough is untouched.
                typename surface__::data_type ulim( trust*( umax - umin ) );
                typename surface__::data_type vlim( trust*( vmax - vmin ) );

                // Each parameter on its own, not the whole step scaled to preserve its direction.
                // Preserving the direction sounds like the careful choice and it is the wrong one
                // here: near a crease Sv is small, so the v part of the step blows up -- traced
                // wanting one and a half periods -- and scaling the pair to fit drags a perfectly
                // good u step down with it, to 6.6e-3 where the tangent plane had asked for 0.32.
                // The two parameters are not commensurable in the first place, and a box is the
                // usual shape for a trust region on unknowns that are not.
                // ...except where one of the parameter directions is nearly collapsed.  There the
                // step in that direction is a ratio with a vanishing denominator and means nothing
                // -- near a pod's nose Sv is a fortieth of Su and the tangent plane asks to travel
                // most of the way round the ring -- so a box, which grants it the whole allowance,
                // takes the iteration somewhere it has no reason to go.  Scaling the pair together
                // there keeps the step pointing where the reliable component says.
                typename surface__::data_type sumag( Su.norm() ), svmag( Sv.norm() );
                typename surface__::data_type flat( static_cast<typename surface__::data_type>(0.05) );

                bool unreliable = ( sumag < flat*svmag ) || ( svmag < flat*sumag );

                if ( std::abs( dx(0) ) > ulim )
                {
                  if ( unreliable )
                  {
                    dx(1) *= ulim/std::abs( dx(0) );
                  }
                  dx(0) *= ulim/std::abs( dx(0) );
                }
                if ( std::abs( dx(1) ) > vlim )
                {
                  if ( unreliable )
                  {
                    dx(0) *= vlim/std::abs( dx(1) );
                  }
                  dx(1) *= vlim/std::abs( dx(1) );
                }

                dx = this->calculate_delta_factor(x, dx);
                x+=dx;

                prev_dx = abs_x_norm;
                abs_x_norm = this->calculate_norm(dx);

                if ( divflag && ( dnow > dprev ) ) // Getting further away twice in a row.
                {
                  // Solution diverging, return initial guess.
                  root = xbest;
                  return this->no_root_found;
                }
                else if ( dnow > dprev ) // Further away, first time
                {
                  divflag = true;
                }
                else // Not diverging.
                {
                  divflag = false;
                }

                all_zero = true;
                for (size_t i=0; i<N__; ++i)
                {
                  // check if stuck and cannot move x anymore
                  if ( std::abs( dx(i) ) > std::numeric_limits<typename surface__::data_type>::epsilon() )
                  {
                    all_zero = false;
                    break;
                  }
                }

                ++count;
              }

              // The loop evaluates at the top and steps at the bottom, so the point it leaves in
              // x has never been measured.  Measure it once and keep whichever of it and the best
              // point visited is actually closer.  Reading xbest only on the divergence return
              // would throw away the best answer the run found on every other exit -- and
              // max_iteration, not convergence, is the exit that dominates here.  One f()
              // against up to maxit f_pt_derivs() is not a cost worth the wrong answer.
              q = s->f( x(0), x(1) );
              r = q - pt;

              if ( dbest < r.norm() )
              {
                root = xbest;
              }
              else
              {
                root = x;
              }

              if ( all_zero )
              {
                return this->hit_constraint;
              }

              if ( count >= maxit )
              {
                return this->max_iteration;
              }

              return this->converged;
            }

          private:
            typename mutil::nls::iterative_system_root_base<typename surface__::data_type, N__, NSOL__>::solution_matrix x0;
        };

        template <typename surface__>
        struct surface_g_gp_functor
        {
          const surface__ *ps;
          typename surface__::point_type pt;
          typedef typename Eigen::Matrix<typename surface__::data_type, 2, 1> vec;
          typedef typename Eigen::Matrix<typename surface__::data_type, 2, 2> mat;

          void operator()(vec &g, mat &gp, const vec &u) const
          {
            typename surface__::data_type uu(u[0]), vv(u[1]);

            typename surface__::data_type umin, umax, vmin, vmax;
            ps->get_parameter_min(umin,vmin);
            ps->get_parameter_max(umax,vmax);

            if ( !(uu>=umin) )
            {
#ifdef ELI_DEBUG_CLAMP
              std::cout << "Minimum distance surface g_functor, u less than minimum.  uu: " << uu << " umin: " << umin << std::endl;
#endif
              uu=umin;
            }
            if ( !(uu<=umax) )
            {
#ifdef ELI_DEBUG_CLAMP
              std::cout << "Minimum distance surface g_functor, u greater than maximum.  uu: " << uu << " uamx: " << umax << std::endl;
#endif
              uu=umax;
            }

            if ( !(vv>=vmin) )
            {
#ifdef ELI_DEBUG_CLAMP
              std::cout << "Minimum distance surface g_functor, v less than minimum.  vv: " << vv << " vmin: " << vmin << std::endl;
#endif
              vv=vmin;
            }
            if ( !(vv<=vmax) )
            {
#ifdef ELI_DEBUG_CLAMP
              std::cout << "Minimum distance surface g_functor, v greater than maximum.  vv: " << vv << " vmax: " << vmax << std::endl;
#endif
              vv=vmax;
            }

            assert((uu>=umin) && (uu<=umax));
            assert((vv>=vmin) && (vv<=vmax));

            uu=std::min(std::max(uu, static_cast<typename surface__::data_type>(umin)), static_cast<typename surface__::data_type>(umax));
            vv=std::min(std::max(vv, static_cast<typename surface__::data_type>(vmin)), static_cast<typename surface__::data_type>(vmax));

            typename surface__::point_type tmp, Su, Sv, Suu, Suv, Svv;

            tmp=ps->f(uu, vv)-pt;
            Su=ps->f_u(uu, vv);
            Sv=ps->f_v(uu, vv);

            g(0)=tmp.dot(Su);
            g(1)=tmp.dot(Sv);

            Suu=ps->f_uu(uu, vv);
            Suv=ps->f_uv(uu, vv);
            Svv=ps->f_vv(uu, vv);

            gp(0,0)=Su.dot(Su)+tmp.dot(Suu);
            gp(0,1)=Su.dot(Sv)+tmp.dot(Suv);
            gp(1,0)=gp(0,1);
            gp(1,1)=Sv.dot(Sv)+tmp.dot(Svv);

            // TODO: What to do if matrix becomes singular?
          }
        };
      }

      template<typename surface__>
      typename surface__::data_type minimum_distance_tan(typename surface__::data_type &u, typename surface__::data_type &v, const surface__ &s, const typename surface__::point_type &pt,
                                                     const typename surface__::data_type &u0, const typename surface__::data_type &v0, int & ret,
                                                     const typename surface__::data_type uminc = 0, const typename surface__::data_type umaxc = 0,
                                                     const typename surface__::data_type vminc = 0, const typename surface__::data_type vmaxc = 0 )
      {
        typedef internal::tangent_plane_method<surface__, 2, 1> nonlinear_solver_type;

        typename surface__::data_type umin, umax, vmin, vmax;

        bool user_con = false;

        if ( umaxc == 0 && uminc == 0 && vmaxc == 0 && vminc == 0 )
        {
          s.get_parameter_min(umin,vmin);
          s.get_parameter_max(umax,vmax);
        }
        else
        {
          umin = uminc;
          umax = umaxc;
          vmin = vminc;
          vmax = vmaxc;
          user_con = true;
        }

        nonlinear_solver_type tpsolve;

        typename surface__::tolerance_type tol;


        if ( user_con || s.open_u() ) // user_con first: skip the closed_u() boundary-curve sweep when the caller supplied bounds
        {
          tpsolve.set_lower_condition(0, umin, nonlinear_solver_type::IRC_EXCLUSIVE);
          tpsolve.set_upper_condition(0, umax, nonlinear_solver_type::IRC_EXCLUSIVE);
        }
        else
        {
          tpsolve.set_periodic_condition(0, umin, umax);
        }

        if ( user_con || s.open_v() ) // user_con first: skip the closed_v() boundary-curve sweep when the caller supplied bounds
        {
          tpsolve.set_lower_condition(1, vmin, nonlinear_solver_type::IRC_EXCLUSIVE);
          tpsolve.set_upper_condition(1, vmax, nonlinear_solver_type::IRC_EXCLUSIVE);
        }
        else
        {
          tpsolve.set_periodic_condition(1, vmin, vmax);
        }

        // setup the solver
        tpsolve.set_absolute_f_tolerance(tol.get_absolute_tolerance());
        tpsolve.set_max_iteration(20);
        tpsolve.set_norm_type(nonlinear_solver_type::max_norm);

        typename nonlinear_solver_type::solution_matrix x, x0;

        typename surface__::point_type q;
        typename surface__::data_type dist, dist0;



        assert((u0>=umin) && (u0<=umax));
        assert((v0>=vmin) && (v0<=vmax));

        x0(0) = u0;
        x0(1) = v0;
        tpsolve.set_initial_guess( x0 );
        tpsolve.pt = pt;
        tpsolve.s = &s;
        tpsolve.maxit = 20;
        tpsolve.xtol = tol.get_absolute_tolerance();

        q = s.f(x0(0), x0(1));
        dist0 = eli::geom::point::distance(q, pt);

        ret = tpsolve.find_root( x );

        q = s.f(x(0), x(1));
        dist = eli::geom::point::distance(q, pt);

        if ( dist > dist0 )
        {
          x = x0; // Genuinely no progress; the starting guess is the best there is.
          dist = dist0;
        }

        u = x(0);
        v = x(1);
        return dist;
      }

      template<typename surface__>
      typename surface__::data_type minimum_distance_nrm(typename surface__::data_type &u, typename surface__::data_type &v, const surface__ &s, const typename surface__::point_type &pt,
                                                     const typename surface__::data_type &u0, const typename surface__::data_type &v0, int & ret,
                                                     const typename surface__::data_type uminc = 0, const typename surface__::data_type umaxc = 0,
                                                     const typename surface__::data_type vminc = 0, const typename surface__::data_type vmaxc = 0 )
      {
        typedef eli::mutil::nls::newton_raphson_system_method<typename surface__::data_type, 2, 1> nonlinear_solver_type;
        nonlinear_solver_type nrm;
        internal::surface_g_gp_functor<surface__> ggp;
        typename surface__::data_type dist0, dist;
        typename surface__::tolerance_type tol;

        bool user_con = false;

        typename surface__::data_type umin, umax, vmin, vmax;
        if ( umaxc == 0 && uminc == 0 && vmaxc == 0 && vminc == 0 )
        {
          s.get_parameter_min(umin,vmin);
          s.get_parameter_max(umax,vmax);
        }
        else
        {
          umin = uminc;
          umax = umaxc;
          vmin = vminc;
          vmax = vmaxc;
          user_con = true;
        }

        // setup the functors
        ggp.ps=&s;
        ggp.pt=pt;

        // setup the solver
        nrm.set_absolute_f_tolerance(tol.get_absolute_tolerance());
        nrm.set_max_iteration(20);
        nrm.set_norm_type(nonlinear_solver_type::max_norm);

        if ( user_con || s.open_u() ) // user_con first: skip the closed_u() boundary-curve sweep when the caller supplied bounds
        {
          nrm.set_lower_condition(0, umin, nonlinear_solver_type::IRC_EXCLUSIVE);
          nrm.set_upper_condition(0, umax, nonlinear_solver_type::IRC_EXCLUSIVE);
        }
        else
        {
          nrm.set_periodic_condition(0, umin, umax);
        }

        if ( user_con || s.open_v() ) // user_con first: skip the closed_v() boundary-curve sweep when the caller supplied bounds
        {
          nrm.set_lower_condition(1, vmin, nonlinear_solver_type::IRC_EXCLUSIVE);
          nrm.set_upper_condition(1, vmax, nonlinear_solver_type::IRC_EXCLUSIVE);
        }
        else
        {
          nrm.set_periodic_condition(1, vmin, vmax);
        }

        // set the initial guess
        typename nonlinear_solver_type::solution_matrix uinit, rhs, ans;

        uinit(0)=u0;
        uinit(1)=v0;
        nrm.set_initial_guess(uinit);
        rhs.setZero();
        dist0=eli::geom::point::distance(s.f(u0, v0), pt);

        // find the root
        ret = nrm.find_root(ans, ggp, rhs);
        u=ans(0);
        v=ans(1);

        // if root is within bounds and is closer than initial guess
        {
          assert((u>=umin) && (u<=umax));
          assert((v>=vmin) && (v<=vmax));

          u = util::clamp( u, umin, umax );
          v = util::clamp( v, vmin, vmax );

          dist = eli::geom::point::distance(s.f(u, v), pt);
          if  (dist<=dist0)
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
        u=u0;
        v=v0;
        return dist0;
      }

      template<typename surface__>
      typename surface__::data_type minimum_distance_core(typename surface__::data_type &u, typename surface__::data_type &v, const surface__ &s, const typename surface__::point_type &pt,
                                                     const typename surface__::data_type &u0, const typename surface__::data_type &v0,
                                                     const typename surface__::data_type uminc = 0, const typename surface__::data_type umaxc = 0,
                                                     const typename surface__::data_type vminc = 0, const typename surface__::data_type vmaxc = 0 )
      {
        internal::tangent_plane_method<surface__, 2, 1> tan_solver;
        typename surface__::data_type dist_tan, dist_nrmt, dist_nrm0;

        int rett = -1;
        dist_tan = minimum_distance_tan( u, v, s, pt, u0, v0, rett, uminc, umaxc, vminc, vmaxc );

        if ( rett == tan_solver.converged || rett == tan_solver.hit_constraint )
        {
          return dist_tan;
        }

        typename surface__::data_type u0t, v0t;
        u0t = u;
        v0t = v;

        int retn = -1;
        dist_nrmt = minimum_distance_nrm( u, v, s, pt, u0t, v0t, retn, uminc, umaxc, vminc, vmaxc );

        if ( retn == tan_solver.converged )
        {
          if ( dist_nrmt <= dist_tan )
          {
            return dist_nrmt;
          }
        }

        dist_nrm0 = minimum_distance_nrm( u, v, s, pt, u0, v0, retn, uminc, umaxc, vminc, vmaxc );

//        if ( retn != tan_solver.converged )
//        {
//          printf("Nothing converged.\n" );
//          printf("%g %g %g\n", dist_tan, dist_nrmt, dist_nrm0 );
//        }
//        else
//        {
//          printf("Newton x0 converged.\n" );
//          printf("%g %g %g\n", dist_tan, dist_nrmt, dist_nrm0 );
//        }

        if ( dist_nrm0 < dist_tan )
        {
            return dist_nrm0;
        }

        u = u0t;
        v = v0t;
        return dist_tan;
      }


      // Judge the answer the iteration produced, and do something about it if it is not one.
      //
      // The core above is Gauss-Newton, and Gauss-Newton is drawn to any point where the residual
      // is perpendicular to the surface -- the far side ridge of a body as readily as the near
      // side.  It need not even arrive: on a ridge crest the step is small and *grows*, which is
      // what escaping an unstable stationary point looks like, and the iteration budget runs out
      // long before it gets anywhere.  It can also stop at a degenerate pole, where the whole v
      // line collapses to a point so there is no way to move around it, or at a boundary it had no
      // business stopping at.  None of this shows in the return code.
      //
      // Second derivatives settle it.  The objective is 0.5*|r|^2 with r = S - pt, so its Hessian
      // is J^T J + r.S'', which is exact and cheap now that f_pt_derivs2 fetches the whole set in
      // one patch lookup.  A minimum needs it positive definite; a ridge crest fails on r.Svv, a
      // pole fails on the collapsed row, and a point found exactly gives r = 0 and H = J^T J,
      // which is positive definite -- so a correct answer never triggers any of what follows, and
      // pays only for the test.
      //
      // Three things follow from it, in order: a Newton polish with that exact Hessian, for the
      // case where the answer is a minimum but the iteration had not reached it; the same test
      // again on the polished point, together with a gradient check at any active boundary; and a
      // restart from a line of samples when the answer is not a minimum at all.
      template<typename surface__>
      typename surface__::data_type minimum_distance(typename surface__::data_type &u, typename surface__::data_type &v, const surface__ &s, const typename surface__::point_type &pt,
                                                     const typename surface__::data_type &u0, const typename surface__::data_type &v0,
                                                     const typename surface__::data_type uminc = 0, const typename surface__::data_type umaxc = 0,
                                                     const typename surface__::data_type vminc = 0, const typename surface__::data_type vmaxc = 0 )
      {
        typename surface__::data_type dist = minimum_distance_core( u, v, s, pt, u0, v0, uminc, umaxc, vminc, vmaxc );

        typename surface__::data_type umin, umax, vmin, vmax;
        s.get_parameter_min( umin, vmin );
        s.get_parameter_max( umax, vmax );

        typename surface__::data_type edge( std::sqrt( std::numeric_limits< typename surface__::data_type >::epsilon() ) );

        // Is the answer in hand a minimum at all?
        //
        // The tangent plane iteration is Gauss-Newton on the residual, so it is drawn to any point
        // where the residual is perpendicular to the surface -- the far side ridge of a body of
        // revolution as readily as the near side.  Worse, it need not even arrive: traced on a
        // ridge crest the v step is 5e-4 and *grows* by half a percent each pass, which is what
        // escaping an unstable stationary point looks like, and twenty iterations move it by 0.01
        // before the budget runs out.  Both end at a point that is not the answer, and neither
        // shows up in the return code.
        //
        // Second derivatives settle it outright.  The objective is 0.5*|r|^2 with r = S - pt, so
        //
        //     H = [ Su.Su + r.Suu   Su.Sv + r.Suv ]
        //         [ Su.Sv + r.Suv   Sv.Sv + r.Svv ]
        //
        // and a minimum needs H positive definite.  On a ridge crest r points outward while the
        // surface curves away from it, r.Svv is negative and large, and H fails -- exactly the
        // case above.  At a pole one row collapses and it fails too, so this subsumes the
        // collapsed-derivative test it replaces.  Sitting on the surface gives r = 0 and H = J^T J,
        // which is positive definite, so a point found exactly never triggers anything.
        //
        // Three extra surface evaluations, once, after the iteration has finished -- and only ever
        // followed by more work when the answer is demonstrably not a minimum.
        typename surface__::point_type q, su, sv, suu, suv, svv;
        s.f_pt_derivs2( u, v, q, su, sv, suu, suv, svv );

        typename surface__::point_type r( q - pt );

        typename surface__::data_type h00( su.dot( su ) + r.dot( suu ) );
        typename surface__::data_type h01( su.dot( sv ) + r.dot( suv ) );
        typename surface__::data_type h11( sv.dot( sv ) + r.dot( svv ) );

        // Whether the point the iteration handed over was a minimum.  Kept separately from the
        // test after the polish below, because the polish is local: it will happily settle a point
        // that was sitting on a ridge into the nearest minimum on that ridge, which is a better
        // answer and still the wrong one.  A restart is wanted if either test fails.
        bool minimum0 = ( h00 > 0 ) && ( h00*h11 - h01*h01 > 0 );

        {
          // Is this the answer the iteration was heading for?
          //
          // Gauss-Newton drops the r.S'' term from the Hessian, which is what makes the tangent
          // plane cheap and what makes it stall.  Where that term is comparable to J^T J the step
          // is built on the wrong curvature and convergence goes linear at a rate near one:
          // traced here as twenty iterations of a v step that stays at 8e-4 and never shrinks,
          // moving v by 0.019 in total while the distance improves in the sixth decimal.  It ends
          // at max_iteration, short of the answer, and the return code says only that it ran out.
          //
          // The Hessian just computed is the exact one, and it is positive definite, so a full
          // Newton step is a descent step and lands on the answer in a few passes rather than
          // hundreds.  Whether to take one is decided by the gradient, scaled to a cosine so the
          // test means the same thing at any distance: an iteration that genuinely converged has
          // nothing here to trigger it, and pays for the test alone.
          for ( int i = 0; i < 8; ++i )
          {
            typename surface__::data_type g0( r.dot( su ) ), g1( r.dot( sv ) );
            typename surface__::data_type rmag( r.norm() ), sumag( su.norm() ), svmag( sv.norm() );

            // A parameter on its bound with the gradient pushing it outward is held there.  The
            // answer is a constrained one in that parameter, so its gradient need not vanish, and
            // the step is the one dimensional Newton step in the other parameter.  Solving the
            // full system and then clamping would leave the free parameter with a step computed
            // for a move the bound does not allow.
            bool ufix = false;
            bool vfix = false;

            if ( ( u <= umin + edge*( umax - umin ) ) && ( g0 > 0 ) )
            {
              ufix = true;
            }
            else if ( ( u >= umax - edge*( umax - umin ) ) && ( g0 < 0 ) )
            {
              ufix = true;
            }

            if ( s.open_v() )
            {
              if ( ( v <= vmin + edge*( vmax - vmin ) ) && ( g1 > 0 ) )
              {
                vfix = true;
              }
              else if ( ( v >= vmax - edge*( vmax - vmin ) ) && ( g1 < 0 ) )
              {
                vfix = true;
              }
            }

            bool stationary = true;

            if ( rmag > 0 )
            {
              if ( !ufix && ( sumag > 0 ) && ( std::abs( g0 ) > edge*sumag*rmag ) )
              {
                stationary = false;
              }
              if ( !vfix && ( svmag > 0 ) && ( std::abs( g1 ) > edge*svmag*rmag ) )
              {
                stationary = false;
              }
            }

            if ( stationary )
            {
              break;
            }

            typename surface__::data_type det( h00*h11 - h01*h01 );

            typename surface__::data_type du( 0 ), dv( 0 );

            if ( ( svmag <= 0 ) && ( h00 > 0 ) )
            {
              // The surface has no v extent here at all.  A sharp trailing edge is built this
              // way -- the closing strip has zero width, so Sv vanishes along the whole crease --
              // and the nearest point to anything off the back of a wing lies on that crease,
              // where the distance has a corner in v rather than a stationary point.  The tangent
              // plane cannot converge to a corner: traced here it ping-ponged across the seam,
              // 3.958, 3.990, 0.021, 0.005, 3.990, with the trust region halving each time until
              // the step underflowed.
              //
              // v is not a free variable on a crease.  Hold it and solve the one dimensional
              // problem along the edge, which is smooth and has an ordinary minimum.
              du = -g0/h00;
            }
            else if ( ( sumag <= 0 ) && ( h11 > 0 ) )
            {
              dv = -g1/h11;
            }
            else if ( ufix && ( h11 > 0 ) )
            {
              dv = -g1/h11;
            }
            else if ( vfix && ( h00 > 0 ) )
            {
              du = -g0/h00;
            }
            else if ( ( det > 0 ) && ( h00 > 0 ) )
            {
              du = ( -h11*g0 + h01*g1 )/det;
              dv = ( h01*g0 - h00*g1 )/det;
            }
            else
            {
              break;
            }

            // Bounded, for the same reason the tangent plane step is.  A positive but nearly
            // singular determinant gives an arbitrarily long step, and a step of 1e30 in a
            // periodic parameter is not merely useless -- wrapping it back into range one period
            // at a time does not finish in any usable time.
            typename surface__::data_type ulim( ( umax - umin )/4 );
            typename surface__::data_type vlim( ( vmax - vmin )/4 );

            if ( std::abs( du ) > ulim )
            {
              dv *= ulim/std::abs( du );
              du *= ulim/std::abs( du );
            }
            if ( std::abs( dv ) > vlim )
            {
              du *= vlim/std::abs( dv );
              dv *= vlim/std::abs( dv );
            }

            typename surface__::data_type unew( u + du ), vnew( v + dv );

            // Same domain rules the iteration itself works under: u is clamped, v is wrapped where
            // the surface closes and clamped where it does not.
            if ( unew < umin ) { unew = umin; }
            if ( unew > umax ) { unew = umax; }

            if ( s.open_v() )
            {
              if ( vnew < vmin ) { vnew = vmin; }
              if ( vnew > vmax ) { vnew = vmax; }
            }
            else
            {
              if ( vnew < vmin ) { vnew += ( vmax - vmin ); }
              if ( vnew > vmax ) { vnew -= ( vmax - vmin ); }
            }

            typename surface__::point_type qn, sun, svn, suun, suvn, svvn;
            s.f_pt_derivs2( unew, vnew, qn, sun, svn, suun, suvn, svvn );

            typename surface__::point_type rn( qn - pt );
            typename surface__::data_type dn( rn.norm() );

            if ( !( dn < dist ) )
            {
              // Newton overshot.  The point in hand is the better one and there is nothing more
              // this step can offer.
              break;
            }

            u = unew;
            v = vnew;
            dist = dn;
            q = qn;
            su = sun;
            sv = svn;
            r = rn;

            h00 = su.dot( su ) + r.dot( suun );
            h01 = su.dot( sv ) + r.dot( suvn );
            h11 = sv.dot( sv ) + r.dot( svvn );

          }
        }

        // Now judge the answer in hand.  The Hessian is the exact one, and after the polish above
        // it describes the point actually being returned.
        bool minimum = minimum0 && ( h00 > 0 ) && ( h00*h11 - h01*h01 > 0 );

        // A boundary answer is a constrained one, and the Hessian says nothing about whether the
        // constraint is the reason it stopped.  The gradient does: at u = umin the answer is only
        // a minimum if moving u up makes things worse.  Traced on a wing, an iteration ran into
        // the root edge and stopped at d = 22.5 with the gradient still pointing inboard, where
        // the answer was 13.7 -- and the Hessian test passed it, since a constrained minimum need
        // not have a positive definite Hessian and this one was not even constrained.
        //
        // Free quantities; the residual and both derivatives are already in hand.
        if ( u <= umin + edge*( umax - umin ) )
        {
          if ( r.dot( su ) < 0 )
          {
            minimum = false;
          }
        }
        else if ( u >= umax - edge*( umax - umin ) )
        {
          if ( r.dot( su ) > 0 )
          {
            minimum = false;
          }
        }

        if ( s.open_v() )
        {
          if ( v <= vmin + edge*( vmax - vmin ) )
          {
            if ( r.dot( sv ) < 0 )
            {
              minimum = false;
            }
          }
          else if ( v >= vmax - edge*( vmax - vmin ) )
          {
            if ( r.dot( sv ) > 0 )
            {
              minimum = false;
            }
          }
        }


        if ( !minimum )
        {
          // Where to restart, and along which parameter.
          //
          // Half a period round in v is the obvious guess and it is wrong at a pole, which is
          // where most of these end up: every v names the same physical point there, so the v the
          // iteration produced carries no information and adding to it carries none either.
          // Traced on a point just past the nose and off the axis, the half-period restart lands
          // on the far side and walks straight back to the pole.  One line of samples finds the
          // right value outright.
          //
          // Which parameter to sample is the same question again.  A collapsed derivative is
          // precisely what "carries no information" means, so sample the collapsed one: at a pod's
          // nose Sv vanishes and the samples run round the ring, while at a wing's tip cap it is
          // Su that vanishes -- traced at |Su| = 4e-5 against |Sv| = 4e-2, with the distance flat
          // to six figures over the whole cap -- and the samples must run along the span instead.
          // A ring is no use there; it stays inside the cap.  Where neither has collapsed, v is
          // the periodic one and the better bet.
          typename surface__::data_type sumag( su.norm() ), svmag( sv.norm() );
          typename surface__::data_type flat( static_cast<typename surface__::data_type>(0.01) );

          // Not less than, but not greater than.  A wing's tip closing patch is degenerate in
          // both directions at once -- Su and Sv both vanish, and so does the Hessian entirely --
          // and a strict test then picks v, whose samples all stay inside the cap.  Neither
          // parameter carries information there, and it is u that leads out.
          bool scan_u = ( sumag <= flat*svmag );

          typename surface__::data_type unudge( u ), vnudge( v );

          if ( u <= umin + edge*( umax - umin ) )
          {
            unudge = umin + static_cast<typename surface__::data_type>(0.01)*( umax - umin );
          }
          else if ( u >= umax - edge*( umax - umin ) )
          {
            unudge = umax - static_cast<typename surface__::data_type>(0.01)*( umax - umin );
          }

          const int nscan( 12 );

          typename surface__::data_type dscan( std::numeric_limits< typename surface__::data_type >::max() );

          for ( int i = 0; i < nscan; ++i )
          {
            // Cell centres.  Starting at the lower bound would put a sample on the parameter
            // boundary, and a restart from exactly there spends its first steps arguing with the
            // constraint rather than descending.
            typename surface__::data_type t( ( i + static_cast<typename surface__::data_type>(0.5) )/nscan );

            typename surface__::data_type ui( unudge ), vi( v );

            if ( scan_u )
            {
              ui = umin + ( umax - umin )*t;
            }
            else
            {
              vi = vmin + ( vmax - vmin )*t;
            }

            typename surface__::data_type di( ( s.f( ui, vi ) - pt ).norm() );

            if ( di < dscan )
            {
              dscan = di;
              unudge = ui;
              vnudge = vi;
            }
          }

          typename surface__::data_type u2, v2;
          typename surface__::data_type dist2 = minimum_distance_core( u2, v2, s, pt, unudge, vnudge, uminc, umaxc, vminc, vmaxc );

          if ( dist2 < dist )
          {
            u = u2;
            v = v2;
            dist = dist2;
          }
        }

        return dist;
      }

      template<typename surface__>
      typename surface__::data_type minimum_distance_old(typename surface__::data_type &u, typename surface__::data_type &v, const surface__ &s, const typename surface__::point_type &pt)
      {
        typename surface__::tolerance_type tol;

        // possible that end points are closest, so start by checking them
        typename surface__::data_type dist, uu, vv, dd;

        typename surface__::data_type umin, umax, vmin, vmax, uspan, vspan;
        s.get_parameter_min(umin,vmin);
        s.get_parameter_max(umax,vmax);
        uspan=umax-umin;
        vspan=vmax-vmin;

        typename surface__::index_type i, j, nu, nv;
        typename surface__::data_type du, dv;

        nu=2*s.degree_u()+1;
        nv=2*s.degree_v()+1;

        // Evenly spaced in parameter, don't repeat 0/1 if closed curve.
        if (s.open_u())
        {
          du = uspan/(nu-1);
        }
        else
        {
          du = uspan/nu;
        }

        if (s.open_v())
        {
          dv = vspan/(nv-1);
        }
        else
        {
          dv = vspan/nv;
        }

        // Find closest of evenly spaced points.
        uu=umin;
        dist = std::numeric_limits<typename surface__::data_type>::max();
        for (i = 0; i < nu; i++)
        {
          vv=vmin;
          for (j = 0; j < nv; j++)
          {
            dd=eli::geom::point::distance(s.f(uu,vv), pt);

            if( dd < dist )
            {
              u=uu;
              v=vv;
              dist=dd;
            }
            vv+=dv;
            if(vv>=vmax)
            {
              vv=vmax;
            }
          }
          uu+=du;
          if(uu>=umax)
          {
            uu=umax;
          }
        }

        // Polish best point with Newton's method search.
        dd=minimum_distance(uu, vv, s, pt, u, v);

        if ((uu>=umin) && (uu<=umax) && (vv>=vmin) && (vv<=vmax))
        {
          if (dd<dist)
          {
            u=uu;
            v=vv;
            dist=dd;
          }
        }

        // next check edges
        // Since these are always edges, we could implement an edge curve extraction routine
        // that returned the control points directly instead of performing an arbitrary curve
        // extraction calculation.
        typename surface__::curve_type bc;
        if(u<=(umin+std::abs(umin)*2*std::numeric_limits<typename surface__::data_type>::epsilon()))
        {
          s.get_umin_bndy_curve(bc);
          dd=eli::geom::intersect::minimum_distance(vv, bc, pt, v);

          if (dd<dist)
          {
            u=umin;
            v=vv;
            dist=dd;
          }
        }

        if(u>=(umax-std::abs(umax)*2*std::numeric_limits<typename surface__::data_type>::epsilon()))
        {
          s.get_umax_bndy_curve(bc);
          dd=eli::geom::intersect::minimum_distance(vv, bc, pt, v);

          if (dd<dist)
          {
            u=umax;
            v=vv;
            dist=dd;
          }
        }

        if(v<=(vmin+std::abs(vmin)*2*std::numeric_limits<typename surface__::data_type>::epsilon()))
        {
          s.get_vmin_bndy_curve(bc);
          dd=eli::geom::intersect::minimum_distance(uu, bc, pt, u);

          if (dd<dist)
          {
            u=uu;
            v=vmin;
            dist=dd;
          }
        }

        if(v>=(vmax-std::abs(vmax)*2*std::numeric_limits<typename surface__::data_type>::epsilon()))
        {
          s.get_vmax_bndy_curve(bc);
          dd=eli::geom::intersect::minimum_distance(uu, bc, pt, u);

          if (dd<dist)
          {
            u=uu;
            v=vmax;
            dist=dd;
          }
        }

        return dist;

      }

      template<typename surface__>
      typename surface__::data_type minimum_distance_new(typename surface__::data_type &u, typename surface__::data_type &v, const surface__ &s, const typename surface__::point_type &pt)
      {
        typedef typename surface__::onedbezsurf objsurf;
        typedef typename surface__::data_type data_type;
        typedef std::pair< data_type, data_type > uvpair;
        typename std::vector< uvpair >::size_type i;
        data_type uu, vv, dd;

        data_type dist = std::numeric_limits<data_type>::max();

        uvpair start = std::make_pair( 0, 0 );
        uvpair end = std::make_pair( 1, 1 );

        objsurf obj = s.mindistsurf( pt );

        std::vector< uvpair > optpts;
        findnonpos( optpts, start, end, obj, 6 );

        if ( optpts.empty() )
        {
          optpts.push_back( std::make_pair( 0.5, 0.5 ) );
        }

        for ( i = 0; i < optpts.size(); i++ )
        {
          uvpair uv = optpts[i];
          int ret = -1;

          dd = minimum_distance_nrm( uu, vv, s, pt, uv.first, uv.second, ret );

          if ( dd < dist )
          {
            dist = dd;
            u = uu;
            v = vv;
          }
        }

        // next check edges
        typename surface__::data_type umin(0), umax(1), vmin(0), vmax(1);

        typename surface__::curve_type bc;

        s.get_umin_bndy_curve(bc);
        dd=eli::geom::intersect::minimum_distance(vv, bc, pt);

        if ( dd < dist )
        {
          u = umin;
          v = vv;
          dist = dd;
        }

        s.get_umax_bndy_curve(bc);
        dd=eli::geom::intersect::minimum_distance(vv, bc, pt);

        if ( dd < dist )
        {
          u = umax;
          v  =vv;
          dist = dd;
        }

        s.get_vmin_bndy_curve(bc);
        dd=eli::geom::intersect::minimum_distance(uu, bc, pt);

        if ( dd < dist )
        {
          u = uu;
          v = vmin;
          dist = dd;
        }

        s.get_vmax_bndy_curve(bc);
        dd=eli::geom::intersect::minimum_distance(uu, bc, pt);

        if ( dd < dist )
        {
          u = uu;
          v = vmax;
          dist = dd;
        }

        return dist;
      }

      template<typename surface__>
      typename surface__::data_type minimum_distance(typename surface__::data_type &u, typename surface__::data_type &v, const surface__ &s, const typename surface__::point_type &pt)
      {
        return minimum_distance_new( u, v, s, pt );
      }

// Defined for minimum_distance_curve.  Could be moved to util somewhere.
//      template< typename first__, typename second__>
//      bool pairfirstcompare( const std::pair < first__, second__ > &a, const std::pair < first__, second__ > &b )
//      {
//          return ( a.first < b.first );
//      }

      template<template<typename, unsigned short, typename> class surface__, typename data__, unsigned short dim__, typename tol__ >
      typename surface::piecewise<surface__, data__, dim__, tol__>::data_type minimum_distance(
          typename surface::piecewise<surface__, data__, dim__, tol__>::data_type &u,
          typename surface::piecewise<surface__, data__, dim__, tol__>::data_type &v,
          const surface::piecewise<surface__, data__, dim__, tol__> &ps,
          const typename surface::piecewise<surface__, data__, dim__, tol__>::point_type &pt)
      {
        typedef surface::piecewise<surface__, data__, dim__, tol__> piecewise_type;
        typedef typename piecewise_type::index_type index_type;
        typedef typename piecewise_type::data_type data_type;
        typedef typename piecewise_type::bounding_box_type bounding_box_type;

        typedef typename piecewise_type::keymap_type keymap_type;
        typedef typename keymap_type::const_iterator keyit;

        typedef std::pair<keyit, keyit> itpair;
        typedef std::vector< std::pair<data_type, itpair > > dvec;
        dvec minbbdist;

        // Find closest corner of bounding boxes, add them to vector
        // Simple linear search, would be more efficient with some sort of tree.
        for(keyit uit = ps.ukey.key.begin(); uit != ps.ukey.key.end(); ++uit)
        {
          for(keyit vit = ps.vkey.key.begin(); vit != ps.vkey.key.end(); ++vit)
          {
            index_type uk = uit->second;
            index_type vk = vit->second;

            bounding_box_type bb_local;
            ps.patches[uk][vk].get_bounding_box(bb_local);

            data_type dbbmin;
            dbbmin = minimum_distance(bb_local, pt);

            minbbdist.push_back(std::make_pair(dbbmin, std::make_pair(uit, vit)));

          }
        }

        // Sort by nearest distance.
        std::sort( minbbdist.begin(), minbbdist.end(), pairfirstcompare<data_type, itpair > );


        // Iterate over segments, starting with nearest bounding box
        data_type dist(std::numeric_limits<data_type>::max());

        typename dvec::const_iterator it;
        for (it=minbbdist.begin(); it!=minbbdist.end(); ++it)
        {
          // If nearest bb distance is farther than current best, we're done.
          if(it->first < dist )
          {
            itpair itp = it->second;
            keyit uit = itp.first;
            keyit vit = itp.second;

            index_type uk = uit->second;
            index_type vk = vit->second;

            data_type uu, vv, d;
            d=minimum_distance(uu, vv, ps.patches[uk][vk], pt);

            if(d < dist)
            {
              data_type du(ps.ukey.get_delta_parm(uit));
              data_type dv(ps.vkey.get_delta_parm(vit));

              data_type ustart(uit->first);
              data_type vstart(vit->first);

              dist = d;
              u=ustart+uu*du;
              v=vstart+vv*dv;
            }
          }
          else
          {
            break;
          }

        }
        return dist;
      }

    }
  }
}
#endif
