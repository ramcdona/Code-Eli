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

#ifndef eli_geom_intersect_intersect_axis_surface_hpp
#define eli_geom_intersect_intersect_axis_surface_hpp

#include <cmath>
#include <vector>
#include <list>
#include <algorithm>

#include "eli/code_eli.hpp"

#include "eli/mutil/nls/iterative_system_root_base_constrained.hpp"
#include "eli/mutil/nls/newton_raphson_system_method.hpp"

#include "eli/geom/intersect/minimum_distance_curve.hpp"
#include "eli/geom/point/distance.hpp"
#include "eli/geom/intersect/intersect_axis_bounding_box.hpp"
#include "eli/geom/intersect/findnonpos.hpp"

namespace eli
{
  namespace geom
  {
    namespace intersect
    {
      namespace internal
      {
        enum DIR
        {
          X_DIR = 0,
          Y_DIR = 1,
          Z_DIR = 2
        };

        template <typename surface__>
        struct surf_axis_g_functor
        {
          const surface__ *s;
          typename surface__::index_type dir1, dir2;
          typename surface__::point_type p0;

          typedef typename Eigen::Matrix<typename surface__::data_type, 2, 1> vec;


          vec operator()(const vec &x) const
          {
            typename surface__::data_type u(x[0]), v(x[1]);
            vec rtn;

            typename surface__::data_type umin, umax, vmin, vmax;

            s->get_parameter_min(umin,vmin);
            s->get_parameter_max(umax,vmax);

            u=std::min(std::max(u, umin), umax);
            v=std::min(std::max(v, vmin), vmax);

            typename surface__::point_type disp;

            disp = p0 - s->f(u,v);

            rtn(0) = disp(dir1);
            rtn(1) = disp(dir2);

            return rtn;
          }
        };

        template <typename surface__>
        struct surf_axis_gp_functor
        {
          const surface__ *s;
          typename surface__::index_type dir1, dir2;
          typename surface__::point_type p0;

          typedef typename Eigen::Matrix<typename surface__::data_type, 2, 1> vec;
          typedef typename Eigen::Matrix<typename surface__::data_type, 2, 2> mat;

          mat operator()(const vec &x) const
          {
            typename surface__::data_type u(x[0]), v(x[1]);
            mat rtn;

            typename surface__::data_type umin, umax, vmin, vmax;

            s->get_parameter_min(umin,vmin);
            s->get_parameter_max(umax,vmax);

            u=std::min(std::max(u, umin), umax);
            v=std::min(std::max(v, vmin), vmax);

            typename surface__::point_type Su, Sv;

            Su=s->f_u(u, v);
            Sv=s->f_v(u, v);

            rtn(0,0)=-Su(dir1);
            rtn(0,1)=-Sv(dir1);

            rtn(1,0)=-Su(dir2);
            rtn(1,1)=-Sv(dir2);

            // TODO: What to do if matrix becomes singular?

            return rtn;
          }
        };
      }


      template<typename surface__>
      typename surface__::data_type intersect(typename surface__::data_type &u, typename surface__::data_type &v, typename surface__::point_type &p,
                                              const surface__ &s, const typename surface__::point_type &p0, const typename surface__::index_type &iproj,
                                              const typename surface__::data_type &u0, const typename surface__::data_type &v0 )
      {
        typedef eli::mutil::nls::newton_raphson_system_method<typename surface__::data_type, 2, 1> nonlinear_solver_type;
        nonlinear_solver_type nrm;
        internal::surf_axis_g_functor<surface__> g;
        internal::surf_axis_gp_functor<surface__> gp;
        typename surface__::data_type idist;
        typename surface__::tolerance_type tol;

        typename surface__::data_type umin, umax, vmin, vmax;

        typename surface__::index_type dir1, dir2;
        if ( iproj == internal::X_DIR )
        {
          dir1 = internal::Y_DIR;
          dir2 = internal::Z_DIR;
        }
        else if ( iproj == internal::Y_DIR )
        {
          dir1 = internal::Z_DIR;
          dir2 = internal::X_DIR;
        }
        else
        {
          dir1 = internal::X_DIR;
          dir2 = internal::Y_DIR;
        }

        s.get_parameter_min(umin,vmin);
        s.get_parameter_max(umax,vmax);

        u=std::min(std::max(u, umin), umax);
        v=std::min(std::max(v, vmin), vmax);

        // setup the functors
        g.s=&s;
        g.p0=p0;
        g.dir1=dir1;
        g.dir2=dir2;

        gp.s=&s;
        gp.p0=p0;
        gp.dir1=dir1;
        gp.dir2=dir2;

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

        // set the initial guess
        typename nonlinear_solver_type::solution_matrix xinit;
        typename nonlinear_solver_type::solution_matrix rhs;
        typename nonlinear_solver_type::solution_matrix ans;

        xinit(0)=u0;
        xinit(1)=v0;

        nrm.set_initial_guess(xinit);
        rhs.setZero();

        // find the root
        typename surface__::index_type ret = nrm.find_root(ans, g, gp, rhs);
        u=ans(0);
        v=ans(1);

        // if root is within bounds and is a solution
        assert((u>=umin) && (u<=umax));
        assert((v>=vmin) && (v<=vmax));

        p = s.f(u, v);
        idist = std::abs( p(iproj) - p0(iproj) );

        if ( ret == nrm.converged )
        {
          // Solution u, v, p passed back as references.
          return idist;
        }

        // couldn't find better answer so return original point and invalid parameters
        u = -1;
        v = -1;
        p = p0;

        return -1;  // idist is always positive, so this signals non-convergence
      }

      template<typename surface__>
      typename surface__::data_type intersect(typename surface__::data_type &u, typename surface__::data_type &v, typename surface__::point_type &p,
                                              const surface__ &s, const typename surface__::point_type &p0, const typename surface__::index_type &iproj)
      {
        typename surface__::data_type umin, umax, vmin, vmax;

        typedef typename surface__::onedbezsurf objsurf;
        typedef typename surface__::data_type data_type;
        typedef std::pair< data_type, data_type > uvpair;
        typename std::vector< uvpair >::size_type i;
        data_type uu, vv, idd;
        typename surface__::point_type pp;

        data_type idist = std::numeric_limits<data_type>::max();

        uvpair start = std::make_pair( 0, 0 );
        uvpair end = std::make_pair( 1, 1 );

        objsurf obj = s.intaxissurf( p0, iproj );

        std::vector< uvpair > optpts;
        findnonpos( optpts, start, end, obj, 6 );

        if ( optpts.empty() )
        {
          optpts.push_back( std::make_pair( 0.5, 0.5 ) );
        }

        bool converged = false;
        for ( i = 0; i < optpts.size(); i++ )
        {
          uvpair uv = optpts[i];

          idd = intersect( uu, vv, pp, s, p0, iproj, uv.first, uv.second );

          if ( idd >= 0 && idd < idist )
          {
            idist = idd;
            u = uu;
            v = vv;
            p = pp;
            converged = true;
          }
        }

        if ( converged )  // At least one point converged to a solution.
        {
          // Returning u, v, p, idist from best solution.
          return idist;
        }

        // Did not converge, return miss.
        p = p0;
        u = -1;
        v = -1;
        idist = -1;
        return idist;
      }

      template<template<typename, unsigned short, typename> class surface__, typename data__, unsigned short dim__, typename tol__ >
      typename surface::piecewise<surface__, data__, dim__, tol__>::data_type intersect(
          typename surface::piecewise<surface__, data__, dim__, tol__>::data_type &u,
          typename surface::piecewise<surface__, data__, dim__, tol__>::data_type &v,
          typename surface::piecewise<surface__, data__, dim__, tol__>::point_type &p,
          const surface::piecewise<surface__, data__, dim__, tol__> &ps,
          const typename surface::piecewise<surface__, data__, dim__, tol__>::point_type &p0,
          const typename surface::piecewise<surface__, data__, dim__, tol__>::index_type &iproj)
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

        // Find closest distance to bbox in iproj direction
        // Simple linear search, would be more efficient with some sort of tree.
        for(keyit uit = ps.ukey.key.begin(); uit != ps.ukey.key.end(); ++uit)
        {
          for(keyit vit = ps.vkey.key.begin(); vit != ps.vkey.key.end(); ++vit)
          {
            index_type uk = uit->second;
            index_type vk = vit->second;

            bounding_box_type bb_local;
            ps.patches[uk][vk].get_bounding_box(bb_local);

            if ( bb_local.intersect_axis( p0, iproj ) )  // bounding box hit, worth considering.
            {
              data_type dbbmin;
              dbbmin = minimum_distance( bb_local, p0, iproj );  // Find lower bound on projected distance

              minbbdist.push_back( std::make_pair( dbbmin, std::make_pair( uit, vit ) ) );
            }
          }
        }

        // Sort by nearest distance.
        std::sort( minbbdist.begin(), minbbdist.end(), pairfirstcompare<data_type, itpair > );


        // Iterate over segments, starting with nearest bounding box
        data_type idist( std::numeric_limits<data_type>::max() );

        bool converged = false;
        typename dvec::const_iterator it;
        for ( it = minbbdist.begin(); it != minbbdist.end(); ++it )
        {
          // If nearest bb distance is farther than current best, we're done.
          if( it->first < idist )
          {
            itpair itp = it->second;
            keyit uit = itp.first;
            keyit vit = itp.second;

            index_type uk = uit->second;
            index_type vk = vit->second;

            data_type uu, vv, d;
            typename piecewise_type::point_type pp;

            d = intersect( uu, vv, pp, ps.patches[uk][vk], p0, iproj );

            if( d >=0 && d < idist )
            {
              data_type du( ps.ukey.get_delta_parm(uit) );
              data_type dv( ps.vkey.get_delta_parm(vit) );

              data_type ustart( uit->first );
              data_type vstart( vit->first );

              idist = d;
              u = ustart + uu * du;
              v = vstart + vv * dv;
              p = pp;
              converged = true;
            }
          }
          else
          {
            break;
          }
        }

        if ( converged )  // At least one point converged to a solution.
        {
          // Returning u, v, p, idist from best solution.
          return idist;
        }

        // Did not converge, return miss.
        p = p0;
        u = -1;
        v = -1;
        idist = -1;
        return idist;
      }

    }
  }
}
#endif
