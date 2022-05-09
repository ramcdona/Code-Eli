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

#include "eli/geom/intersect.hpp"

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
                                              const typename surface__::data_type &rmin, const typename surface__::data_type &rmax,
                                              const typename surface__::data_type &smin, const typename surface__::data_type &smax,
                                              typename surface__::index_type & ret )
      {
        typedef eli::mutil::nls::newton_raphson_system_method<typename surface__::data_type, 3, 1> nonlinear_solver_type;
        nonlinear_solver_type nrm;
        internal::rst_g_functor<surface__> g;
        internal::rst_gp_functor<surface__> gp;
        typename surface__::data_type dist0, dist;
        typename surface__::tolerance_type tol;

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

      template<typename surface__>
      typename surface__::data_type find_rst( typename surface__::data_type &r, typename surface__::data_type &s, typename surface__::data_type &t,
                                              const surface__ &surf, const typename surface__::point_type &pt,
                                              const typename surface__::data_type &r0, const typename surface__::data_type &s0, const typename surface__::data_type &t0,
                                              typename surface__::index_type & ret )
      {
        typedef eli::mutil::nls::newton_raphson_system_method<typename surface__::data_type, 3, 1> nonlinear_solver_type;
        nonlinear_solver_type nrm;
        internal::rst_g_functor<surface__> g;
        internal::rst_gp_functor<surface__> gp;
        typename surface__::data_type dist0, dist;
        typename surface__::tolerance_type tol;

        typename surface__::data_type rmin( 0.0 ), rmax( 1.0 );
        typename surface__::data_type smin( 0.0 ), smax( 0.5 );

        return find_rst( r, s, t, surf, pt, r0, s0, t0, rmin, rmax, smin, smax, ret );
      }

      template<template<typename, unsigned short, typename> class surface__, typename data__, unsigned short dim__, typename tol__ >
      typename surface::piecewise<surface__, data__, dim__, tol__>::data_type find_rst(
          typename surface::piecewise<surface__, data__, dim__, tol__>::data_type &r,
          typename surface::piecewise<surface__, data__, dim__, tol__>::data_type &s,
          typename surface::piecewise<surface__, data__, dim__, tol__>::data_type &t,
          const surface::piecewise<surface__, data__, dim__, tol__> &ps,
          const typename surface::piecewise<surface__, data__, dim__, tol__>::point_type &pt,
          typename surface::piecewise<surface__, data__, dim__, tol__>::index_type &ret )
      {
        typedef surface::piecewise<surface__, data__, dim__, tol__> piecewise_type;
        typedef typename piecewise_type::index_type index_type;
        typedef typename piecewise_type::data_type data_type;
        typedef typename piecewise_type::bounding_box_type bounding_box_type;
        typedef typename piecewise_type::surface_type patch_type;
        typedef typename surface::piecewise<surface__, data__, dim__, tol__>::point_type point_type;

        typedef std::pair<data_type, data_type> pdata;
        typedef std::pair<index_type, index_type> iuivdata;
        typedef std::vector< std::pair<data_type, iuivdata > > dvec;
        dvec minbbdist;

        bool dbg = false;

        data_type vmin, vmid, vmax;
        data_type umin, umax;
        data_type urng, vrng;
        piecewise_type lower, upper;

        ps.get_parameter_min( umin, vmin );
        ps.get_parameter_max( umax, vmax );
        urng = umax - umin;
        vrng = vmax - vmin;

        vmid = 0.5 * ( vmin + vmax );

        ps.split_v( lower, upper, vmid );
        upper.reverse_v();
        upper.set_v0( vmin );

        piecewise_type::parm_match_u( lower, upper );
        piecewise_type::parm_match_v( lower, upper );

        // Not needed as no in-depth control point manipulations happen.
        // piecewise_type::order_match_u( lower, upper );
        // piecewise_type::order_match_v( lower, upper );

        index_type nu = lower.number_u_patches();
        index_type nv = lower.number_v_patches();

        // std::cout << "Number of patches " << nu * nv << std::endl;

        // Find closest corner of bounding boxes, add them to vector
        // Simple linear search, would be more efficient with some sort of tree.
        for( index_type i = 0; i < nu; i++ )
        {
          for( index_type j = 0; j < nv; j++ )
          {
            data_type ustart, du, vstart, dv;
            patch_type *plow = lower.get_patch( i, j, ustart, du, vstart, dv );
            patch_type *pup = upper.get_patch( i, j );

            bounding_box_type bb, bb_up;
            plow->get_bounding_box( bb );
            pup->get_bounding_box( bb_up );
            bb.add( bb_up );

            data_type dbbmin;
            dbbmin = minimum_distance( bb, pt );

            if ( dbg )
            {
              std::cout << "i " << i << " j " << j << " d " << dbbmin << std::endl;
            }

            minbbdist.push_back( std::make_pair( dbbmin, std::make_pair( i, j ) ) );
          }
        }

        // std::cout << std::endl << std::endl;

        // Sort by nearest distance.
        std::sort( minbbdist.begin(), minbbdist.end(), pairfirstcompare< data_type, iuivdata > );

        if ( dbg )
        {
            std::cout << "Number of trial points " << minbbdist.size() << std::endl;
            std::cout << "Largest d0 " << (minbbdist.back()).first << std::endl;
        }

        data_type dist( std::numeric_limits< data_type >::max() );

        index_type count = 0;
        typename dvec::const_iterator it;
        for ( it = minbbdist.begin(); it != minbbdist.end(); ++it )
        {
          // std::cout << "count " << count << " d0 " << it->first << " dist " << dist << std::endl;

          // If current best is farther than the next nearest bb distance, try again
          // Place sqrt(eps) tolerance on bb distance so we don't keep trying when we've achieved
          // 1e-16, but the bb distance is 0.0.
          // It is likely that the first few bb distances will be 0.0 as multiple bounding boxes
          // may contain our point of interest.
          if( dist > ( it->first + std::sqrt( std::numeric_limits< data_type >::epsilon() ) ) )
          {
            data_type dbb = it->first;
            iuivdata iuiv = it->second;

            index_type iu = iuiv.first;
            index_type iv = iuiv.second;

            data_type ustart, du, vstart, dv;

            patch_type *plow = lower.get_patch( iu, iv, ustart, du, vstart, dv );
            // patch_type *pup = upper.get_patch( iu, iv );

            data_type u = ustart + 0.5 * du;
            data_type v = vstart + 0.5 * dv;

            data_type r0 = ( u - umin ) / urng;
            data_type s0 = ( v - vmin ) / vrng;
            data_type t0 = 0.5;

            data_type rmin = ( ustart - umin ) / urng;
            data_type rmax = ( ustart + du - umin ) / urng;

            data_type smin = ( vstart - vmin ) / vrng;
            data_type smax = ( vstart + dv - vmin ) / vrng;

            data_type rr, ss, tt;
            data_type d;
            index_type rret;

            d = find_rst( rr, ss, tt, ps, pt, r0, s0, t0, rmin, rmax, smin, smax, rret );

            if ( dbg )
            {
                std::cout << "iu " << iu << " iv " << iv << " dbb " << dbb << std::endl;
                std::cout << "Target x " << pt.x() << " y " << pt.y() << " z " << pt.z() << std::endl;
                std::cout << "rmin: " << rmin << " r0: " << r0 << " rmax: " << rmax << std::endl;
                std::cout << "smin: " << smin << " s0: " << s0 << " smax: " << smax << std::endl;
                std::cout << "t0: " << t0 << std::endl;
                std::cout << "d " << d << " rret " << rret << std::endl << std::endl << std::endl;
            }

            if( d < dist ) // 0 means converged.
            {
              dist = d;
              r = rr;
              s = ss;
              t = tt;
              ret = rret;
            }
          }
          else
          {
            break;
          }
          count++;
        }

        if ( ret != 0 ) // Not converged.  Try surface projection.
        {
          data_type u, v, d;
          d = minimum_distance( u, v, ps, pt );

          if ( dbg )
          {
            std::cout << "Resorted to surface projection." << std::endl;
          }

          if( d < dist ) // Check for improvement.
          {
            data_type rr = ( u - umin ) / urng;
            data_type ss = ( v - vmin ) / vrng;
            data_type tt = 0.0; // Assume on bottom surface.

            if ( ss > 0.5 ) // Actually on top surface.
            {
              ss = 1.0 - ss;
              tt = 1.0;
            }
            index_type rret = 0;

            if ( dbg )
            {
              std::cout << "Surface projection success." << std::endl;
              std::cout << "d " << d << std::endl << std::endl << std::endl;
            }

            dist = d;
            r = rr;
            s = ss;
            t = tt;
            ret = rret;
          }
          else
          {
            if ( dbg )
            {
              std::cout << "Surface projection failure." << std::endl;
            }
          }
        }

        return dist;
      }

    }
  }
}
#endif
