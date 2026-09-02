/*********************************************************************************
* Copyright (c) 2013 David D. Marshall <ddmarsha@calpoly.edu>
*
* All rights reserved. This program and the accompanying materials
* are made available under the terms of the Eclipse Public License v1.0
* which accompanies this distribution, and is available at
* http://www.eclipse.org/legal/epl-v10.html
*
* Contributors:
*    Rob McDonald
********************************************************************************/

#ifndef eli_geom_intersect_equidistant_surface_hpp
#define eli_geom_intersect_equidistant_surface_hpp

#include <cmath>
#include <iostream>
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
        // The point on a surface that is equidistant from two given points, and of all such
        // points the one nearest to them.
        //
        // Halving an edge that runs between two points on a surface cannot be done by taking
        // the midpoint of the two and projecting it: the midpoint need not lie anywhere near
        // the surface, and where the surface curves away -- around a wing tip, say -- the
        // nearest point to it is back beside one of the ends.  Nor does least sum of squared
        // distances help, since d0^2 + d1^2 is 2|S-m|^2 plus a constant and so has its
        // minimum at exactly that projection.
        //
        // What is wanted is the least of d0 subject to d0 = d1.  Written out:
        //
        //   d0^2 - d1^2 is linear in S, and vanishes on the plane that bisects the two
        //   points, so the constraint is simply
        //
        //     F0 = ( S - m ) . n                                with m the midpoint, n = p1 - p0
        //
        //   and at the constrained minimum the gradient of d0^2 is parallel to the gradient
        //   of F0, which with the multiplier eliminated is
        //
        //     F1 = ( (S-p0).Su )( Sv.n ) - ( (S-p0).Sv )( Su.n )
        //
        // Two equations in u and v, solved by Newton with the Jacobian written out below.
        template <typename surface__>
        struct surface_equidistant_functor
        {
          const surface__ *ps;
          typename surface__::point_type p0;
          typename surface__::point_type p1;

          typedef typename Eigen::Matrix<typename surface__::data_type, 2, 1> vec;
          typedef typename Eigen::Matrix<typename surface__::data_type, 2, 2> mat;

          void operator()(vec &f, mat &J, const vec &u) const
          {
            typedef typename surface__::data_type data_type;
            typedef typename surface__::point_type point_type;

            data_type uu(u[0]), vv(u[1]);

            data_type umin, umax, vmin, vmax;
            ps->get_parameter_min(umin,vmin);
            ps->get_parameter_max(umax,vmax);

            uu=std::min(std::max(uu, umin), umax);
            vv=std::min(std::max(vv, vmin), vmax);

            point_type m( (p0+p1)/static_cast<data_type>(2) );
            point_type n( p1-p0 );

            point_type S( ps->f(uu, vv) );
            point_type Su( ps->f_u(uu, vv) );
            point_type Sv( ps->f_v(uu, vv) );
            point_type Suu( ps->f_uu(uu, vv) );
            point_type Suv( ps->f_uv(uu, vv) );
            point_type Svv( ps->f_vv(uu, vv) );

            point_type r( S-p0 );

            data_type A( r.dot(Su) );          // d(d0^2)/du / 2
            data_type B( r.dot(Sv) );          // d(d0^2)/dv / 2
            data_type C( Sv.dot(n) );          // dF0/dv
            data_type D( Su.dot(n) );          // dF0/du

            f(0) = (S-m).dot(n);
            f(1) = A*C - B*D;

            data_type dAdu( Su.dot(Su) + r.dot(Suu) );
            data_type dAdv( Sv.dot(Su) + r.dot(Suv) );
            data_type dBdu( Su.dot(Sv) + r.dot(Suv) );
            data_type dBdv( Sv.dot(Sv) + r.dot(Svv) );

            data_type dCdu( Suv.dot(n) );
            data_type dCdv( Svv.dot(n) );
            data_type dDdu( Suu.dot(n) );
            data_type dDdv( Suv.dot(n) );

            J(0,0) = D;
            J(0,1) = C;
            J(1,0) = dAdu*C + A*dCdu - dBdu*D - B*dDdu;
            J(1,1) = dAdv*C + A*dCdv - dBdv*D - B*dDdv;
          }
        };
      }

      // The equidistant point along the straight line in parameter space between two
      // parameter locations.
      //
      // Restricting the search to that line turns the two by two system into one equation in
      // one unknown, and a well behaved one: at t=0 the point is p0's parameter, so d0-d1 is
      // negative, and at t=1 it is p1's, so it is positive.  A root therefore always exists in
      // between and bisection always finds it.  There is no seed to get wrong and no
      // stationary point to be caught on.
      //
      // The answer is not in general the nearest equidistant point on the surface -- it is
      // only the nearest one on this line -- but it is equidistant, which the parametric
      // midpoint is not, and that makes it a far better place to start a two dimensional
      // solve from.
      //
      // Composing the surface with the line would make xyz(t) a polynomial of degree m+n and
      // the root an analytic one.  Code-Eli cannot yet compose a Bezier with a Bezier, so the
      // root is bisected instead; the interface would not change if it could.
      template<typename surface__>
      void equidistant_uwline(typename surface__::data_type &u, typename surface__::data_type &v,
                              const surface__ &s,
                              const typename surface__::point_type &p0,
                              const typename surface__::point_type &p1,
                              const typename surface__::data_type &u0,
                              const typename surface__::data_type &v0,
                              const typename surface__::data_type &u1,
                              const typename surface__::data_type &v1)
      {
        typedef typename surface__::data_type data_type;
        typedef typename surface__::point_type point_type;

        data_type tlo(0), thi(1);

        // f(t) = d0^2 - d1^2 along the line.  Negative at one end, positive at the other.
        //
        // Few steps on purpose.  This is a seed, not the answer: what follows it either
        // solves properly or searches, and both start from wherever this leaves off.  Each
        // step costs a surface evaluation and is paid on every split, so a tight bracket here
        // is worth less than it costs -- eight halvings put it inside a two hundredth of the
        // edge, which is plenty to start from.
        for ( int i = 0; i < 8; ++i )
        {
          data_type tm( ( tlo + thi ) / static_cast<data_type>(2) );

          point_type q( s.f( u0 + ( u1 - u0 ) * tm, v0 + ( v1 - v0 ) * tm ) );

          data_type f( ( q - p0 ).squaredNorm() - ( q - p1 ).squaredNorm() );

          if ( f < static_cast<data_type>(0) )
          {
            tlo = tm;
          }
          else
          {
            thi = tm;
          }
        }

        data_type t( ( tlo + thi ) / static_cast<data_type>(2) );

        u = u0 + ( u1 - u0 ) * t;
        v = v0 + ( v1 - v0 ) * t;
      }

      // Find (u,v) on s equidistant from p0 and p1, and as near to them as such a point can
      // be.  u0, v0 seed the search; the answer is held inside [ullim,uulim] x [vllim,vulim].
      //
      // Returns the distance from the answer to p0.  ret carries the solver's return code.
      template<typename surface__>
      typename surface__::data_type equidistant(typename surface__::data_type &u, typename surface__::data_type &v,
                                                const surface__ &s,
                                                const typename surface__::point_type &p0,
                                                const typename surface__::point_type &p1,
                                                const typename surface__::data_type &u0,
                                                const typename surface__::data_type &v0,
                                                const typename surface__::data_type &ullim,
                                                const typename surface__::data_type &uulim,
                                                const typename surface__::data_type &vllim,
                                                const typename surface__::data_type &vulim,
                                                int &ret)
      {
        typedef eli::mutil::nls::newton_raphson_system_method<typename surface__::data_type, 2, 1> nonlinear_solver_type;
        typedef typename surface__::data_type data_type;

        nonlinear_solver_type nrm;
        internal::surface_equidistant_functor<surface__> feq;
        typename surface__::tolerance_type tol;

        feq.ps=&s;
        feq.p0=p0;
        feq.p1=p1;

        nrm.set_absolute_f_tolerance(tol.get_absolute_tolerance());
        nrm.set_max_iteration(20);
        nrm.set_norm_type(nonlinear_solver_type::max_norm);

        nrm.set_lower_condition(0, ullim, nonlinear_solver_type::IRC_EXCLUSIVE);
        nrm.set_upper_condition(0, uulim, nonlinear_solver_type::IRC_EXCLUSIVE);
        nrm.set_lower_condition(1, vllim, nonlinear_solver_type::IRC_EXCLUSIVE);
        nrm.set_upper_condition(1, vulim, nonlinear_solver_type::IRC_EXCLUSIVE);

        typename nonlinear_solver_type::solution_matrix uinit, rhs, ans;
        uinit(0)=u0;
        uinit(1)=v0;
        nrm.set_initial_guess(uinit);
        rhs(0)=static_cast<data_type>(0);
        rhs(1)=static_cast<data_type>(0);

        ret = nrm.find_root(ans, feq, rhs);

        u=ans(0);
        v=ans(1);

        // The seed is kept when the solve does not beat it.  What matters to a caller
        // splitting an edge is the larger of the two distances, since that is the child edge
        // it will have to deal with, so that is what the two are judged on.
        typename surface__::point_type pa( s.f(u, v) );
        typename surface__::point_type pb( s.f(u0, v0) );

        data_type fa( std::max( eli::geom::point::distance(pa, p0), eli::geom::point::distance(pa, p1) ) );
        data_type fb( std::max( eli::geom::point::distance(pb, p0), eli::geom::point::distance(pb, p1) ) );

        if ( !( fa <= fb ) )        // catches NaN as well as a worse answer
        {
          u=u0;
          v=v0;
          return eli::geom::point::distance(pb, p0);
        }

        return eli::geom::point::distance(pa, p0);
      }
    }
  }
}
#endif
