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

#ifndef eli_geom_surface_piecewise_hpp
#define eli_geom_surface_piecewise_hpp

#include <vector>
#include <iterator>
#include <utility>
#include <algorithm>

#include "eli/code_eli.hpp"

#include "eli/util/tolerance.hpp"
#include "eli/geom/curve/piecewise.hpp"

#include "eli/geom/general/continuity.hpp"
#include "eli/geom/curve/equivalent_curves.hpp"
#include "eli/geom/intersect/specified_distance_curve.hpp"

namespace eli
{
  namespace geom
  {

    namespace surface
    {
      template<template<typename, unsigned short, typename> class surface__, typename data__, unsigned short dim__, typename tol__ >
      class piecewise;
    }

    namespace intersect
    {
      template<template<typename, unsigned short, typename> class surface__, typename data__, unsigned short dim__, typename tol__ >
      typename surface::piecewise<surface__, data__, dim__, tol__>::data_type
        find_rst(
          typename surface::piecewise<surface__, data__, dim__, tol__>::data_type &r,
          typename surface::piecewise<surface__, data__, dim__, tol__>::data_type &s,
          typename surface::piecewise<surface__, data__, dim__, tol__>::data_type &t,
          const surface::piecewise<surface__, data__, dim__, tol__> &ps,
          const typename surface::piecewise<surface__, data__, dim__, tol__>::point_type &pt,
          typename surface::piecewise<surface__, data__, dim__, tol__>::index_type &ret );

      template<template<typename, unsigned short, typename> class surface1__, typename data1__, unsigned short dim1__, typename tol1__ >
      typename surface::piecewise<surface1__, data1__, dim1__, tol1__>::data_type
        minimum_distance(
          typename surface::piecewise<surface1__, data1__, dim1__, tol1__>::data_type &u,
	      typename surface::piecewise<surface1__, data1__, dim1__, tol1__>::data_type &v,
          const surface::piecewise<surface1__, data1__, dim1__, tol1__> &ps,
          const typename surface::piecewise<surface1__, data1__, dim1__, tol1__>::point_type &pt);

      template<template<typename, unsigned short, typename> class surface1__, typename data1__, unsigned short dim1__, typename tol1__ >
      typename surface::piecewise<surface1__, data1__, dim1__, tol1__>::data_type
        intersect(
          typename surface::piecewise<surface1__, data1__, dim1__, tol1__>::data_type &u,
          typename surface::piecewise<surface1__, data1__, dim1__, tol1__>::data_type &v,
          typename surface::piecewise<surface1__, data1__, dim1__, tol1__>::point_type &p,
          const surface::piecewise<surface1__, data1__, dim1__, tol1__> &ps,
          const typename surface::piecewise<surface1__, data1__, dim1__, tol1__>::point_type &p0,
          const typename surface::piecewise<surface1__, data1__, dim1__, tol1__>::index_type &iproj);

      template<template<typename, unsigned short, typename> class surface1__, typename data1__, unsigned short dim1__, typename tol1__ >
      void
        intersect_segment( std::vector<data1__> &tvec,
          const surface::piecewise<surface1__, data1__, dim1__, tol1__> &ps,
          const typename surface::piecewise<surface1__, data1__, dim1__, tol1__>::point_type &pt,
          const typename surface::piecewise<surface1__, data1__, dim1__, tol1__>::point_type &vec,
          const typename surface::piecewise<surface1__, data1__, dim1__, tol1__>::bounding_box_type &bbox);
    }

    namespace surface
    {
      template<template<typename, unsigned short, typename> class surface__, typename data__, unsigned short dim__, typename tol__=eli::util::tolerance<data__> >
      class piecewise
      {
        public:
          typedef surface__<data__, dim__, tol__> surface_type;
          typedef typename surface_type::index_type index_type;
          typedef typename surface_type::point_type point_type;
          typedef typename surface_type::control_point_type control_point_type;
          typedef typename surface_type::rotation_matrix_type rotation_matrix_type;
          typedef typename surface_type::bounding_box_type bounding_box_type;
          typedef typename surface_type::curve_type curve_type;
          typedef eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, dim__, tol__> piecewise_curve_type;

          typedef std::pair< index_type, index_type > patch_boundary_code_type;

          typedef data__ data_type;
          typedef unsigned short dimension_type;
          typedef tol__ tolerance_type;
          enum error_code
          {
            NO_ERRORS=0,
            INVALID_INDEX=1,
            INDEX_NOT_FOUND=2,
            INVALID_PARAM=50,
            INVALID_PARAM_DIFFERENCE=51,
            PATCH_NOT_CONNECTED=100,
            UNKNOWN_ERROR=999
          };

        public:
          piecewise() : nu(0), nv(0), uclosecache(UNKNOWN), vclosecache(UNKNOWN) {}
          piecewise(const piecewise<surface__, data_type, dim__, tol__> &p)
            : patches(p.patches), ukey(p.ukey), vkey(p.vkey), nu(p.nu), nv(p.nv),
            uclosecache(p.uclosecache), vclosecache(p.vclosecache) {}
          ~piecewise() {}

          piecewise & operator=(const piecewise<surface__, data_type, dim__> &p)
          {
            if (this==&p)
              return (*this);

            patches=p.patches;
            ukey=p.ukey;
            vkey=p.vkey;
            nu=p.nu;
            nv=p.nv;
            uclosecache=p.uclosecache;
            vclosecache=p.vclosecache;

            return (*this);
          }

          bool operator==(const piecewise<surface__, data_type, dim__> &p) const
          {
            if (this==&p)
              return true;
            if (nu!=p.nu)
              return false;
            if (nv!=p.nv)
              return false;
            if (ukey!=p.ukey)
              return false;
            if (vkey!=p.vkey)
              return false;
            if (number_u_patches()!=p.number_u_patches())
              return false;
            if (number_v_patches()!=p.number_v_patches())
              return false;
            if (uclosecache!=p.uclosecache)
              return false;
            if (vclosecache!=p.vclosecache)
              return false;
            typename patch_collection_type::const_iterator scit, it;
            for (scit=patches.begin(), it=p.patches.begin(); scit!=patches.end(); ++scit, ++it)
            {
              if ((*it)!=(*scit))
                return false;
            }

            return true;
          }

          bool operator!=(const piecewise<surface__, data_type, dim__> &p) const
          {
            return !operator==(p);
          }

          static dimension_type dimension() {return dim__;}

          data_type get_u0() const {return ukey.get_pmin();}
          void set_u0(const data_type &u0_in) {ukey.set_pmin(u0_in);}

          data_type get_v0() const {return vkey.get_pmin();}
          void set_v0(const data_type &v0_in) {vkey.set_pmin(v0_in);}

          data_type get_umax() const {return ukey.get_pmax();}
          data_type get_vmax() const {return vkey.get_pmax();}

          index_type number_u_patches() const {return nu;}
          index_type number_v_patches() const {return nv;}

          surface_type * get_patch( const index_type &ui, const index_type &vi)
          {
              index_type uk, vk;
              find_patch(uk, vk, ui, vi);
              return &patches[uk][vk];
          }

          const surface_type * get_patch( const index_type &ui, const index_type &vi) const
          {
              index_type uk, vk;
              find_patch(uk, vk, ui, vi);
              return &patches[uk][vk];
          }

          surface_type * get_patch( const index_type &ui, const index_type &vi, data_type &ustart, data_type &du, data_type &vstart, data_type &dv)
          {
              index_type uk, vk;
              typename keymap_type::const_iterator uit, vit;

              find_patch( uk, vk, uit, vit, ui, vi );

              ustart = uit->first;
              du = ukey.get_delta_parm( uit );

              vstart = vit->first;
              dv = vkey.get_delta_parm( vit );

              return &patches[uk][vk];
          }

          const surface_type * get_patch( const index_type &ui, const index_type &vi, data_type &ustart, data_type &du, data_type &vstart, data_type &dv) const
          {
              index_type uk, vk;
              typename keymap_type::const_iterator uit, vit;

              find_patch( uk, vk, uit, vit, ui, vi );

              ustart = uit->first;
              du = ukey.get_delta_parm( uit );

              vstart = vit->first;
              dv = vkey.get_delta_parm( vit );

              return &patches[uk][vk];
          }

          surface_type * get_patch_unordered( const index_type &uk, const index_type &vk)
          {
              return &patches[uk][vk];
          }

          const surface_type * get_patch_unordered( const index_type &uk, const index_type &vk) const
          {
              return &patches[uk][vk];
          }

          void get_parameter_min(data_type &umin, data_type &vmin) const
          {
            umin=ukey.get_pmin();
            vmin=vkey.get_pmin();
          }

          void get_parameter_max(data_type &umax, data_type &vmax) const
          {
            umax=ukey.get_pmax();
            vmax=vkey.get_pmax();
          }

          index_type find_u( const data_type &u_in )
          {
            index_type uk;
            typename keymap_type::iterator uit;
            data_type uu(0);

            ukey.find_segment( uk, uit, uu, u_in );

            return ukey.find_index( uit );
          }

          index_type find_v( const data_type &v_in )
          {
            index_type vk;
            typename keymap_type::iterator vit;
            data_type vv(0);

            vkey.find_segment( vk, vit, vv, v_in );

            return vkey.find_index( vit );
          }

          void parameter_report() const
          {
            printf("U parameter:\n");
            ukey.parameter_report();
            printf("V parameter:\n");
            vkey.parameter_report();
          }

          void octave_print(int figno ) const
          {
            index_type i, j, pp, qq, nup, nvp;
            data_type umin, vmin, umax, vmax;

            nup = number_u_patches();
            nvp = number_v_patches();
            get_parameter_min(umin, vmin);
            get_parameter_max(umax, vmax);

            std::cout << "figure(" << figno << ");" << std::endl;

            // initialize the u & v parameters
            std::vector<data__> u(31), v(31);
            for (i=0; i<static_cast<index_type>(u.size()); ++i)
            {
              u[i]=umin+(umax-umin)*static_cast<data__>(i)/(u.size()-1);
            }
            for (j=0; j<static_cast<index_type>(v.size()); ++j)
            {
              v[j]=vmin+(vmax-vmin)*static_cast<data__>(j)/(v.size()-1);
            }

            // set the surface points
            std::cout << "surf_x=[";
            for (i=0; i<static_cast<index_type>(u.size()); ++i)
            {
              std::cout << this->f(u[i], v[0]).x();
              for (j=1; j<static_cast<index_type>(v.size()-1); ++j)
              {
                std::cout << ", " << this->f(u[i], v[j]).x();
              }
              j=static_cast<index_type>(v.size()-1);
              std::cout << ", " << this->f(u[i], v[j]).x();
              if (i<static_cast<index_type>(u.size()-1))
                std::cout << "; " << std::endl;
            }
            std::cout << "];" << std::endl;

            std::cout << "surf_y=[";
            for (i=0; i<static_cast<index_type>(u.size()); ++i)
            {
              std::cout << f(u[i], v[0]).y();
              for (j=1; j<static_cast<index_type>(v.size()-1); ++j)
              {
                std::cout << ", " << f(u[i], v[j]).y();
              }
              j=static_cast<index_type>(v.size()-1);
              std::cout << ", " << f(u[i], v[j]).y();
              if (i<static_cast<index_type>(u.size()-1))
                std::cout << "; " << std::endl;
            }
            std::cout << "];" << std::endl;

            std::cout << "surf_z=[";
            for (i=0; i<static_cast<index_type>(u.size()); ++i)
            {
              std::cout << f(u[i], v[0]).z();
              for (j=1; j<static_cast<index_type>(v.size()-1); ++j)
              {
                std::cout << ", " << f(u[i], v[j]).z();
              }
              j=static_cast<index_type>(v.size()-1);
              std::cout << ", " << f(u[i], v[j]).z();
              if (i<static_cast<index_type>(u.size()-1))
                std::cout << "; " << std::endl;
            }
            std::cout << "];" << std::endl;

            std::cout << "setenv('GNUTERM', 'x11');" << std::endl;
            std::cout << "mesh(surf_x, surf_y, surf_z, zeros(size(surf_z)), 'EdgeColor', [0 0 0]);" << std::endl;
            std::cout << "axis equal" << std::endl;
            std::cout << "axis off" << std::endl;
          }

          void get_pmap_u(std::vector<data_type> &pmap) const
          {
            ukey.get_pmap(pmap);
          }

          void get_pmap_v(std::vector<data_type> &pmap) const
          {
            vkey.get_pmap(pmap);
          }

          void get_pmap_uv(std::vector<data_type> &upmap, std::vector<data_type> &vpmap) const
          {
            ukey.get_pmap(upmap);
            vkey.get_pmap(vpmap);
          }

          void init_u(const index_type &nsegu, const data_type &du = 1, const data_type &u0 = 0)
          {
            patches.clear();
            resize_store(nsegu, nv);
            ukey.init(nsegu, du, u0);
            uclosecache=UNKNOWN;
            vclosecache=UNKNOWN;
          }

          void init_v(const index_type &nsegv, const data_type &dv = 1, const data_type &v0 = 0)
          {
            patches.clear();
            resize_store(nu, nsegv);
            vkey.init(nsegv, dv, v0);
            uclosecache=UNKNOWN;
            vclosecache=UNKNOWN;
          }

          void init_uv(const index_type &nsegu, const index_type &nsegv, const data_type &du = 1, const data_type &dv = 1, const data_type &u0 = 0, const data_type &v0 = 0)
          {
            patches.clear();
            resize_store(nsegu, nsegv);
            ukey.init(nsegu, du, u0);
            vkey.init(nsegv, dv, v0);
            uclosecache=UNKNOWN;
            vclosecache=UNKNOWN;
          }

          template<typename it__>
          void init_u(const it__ &dus, const it__ &due, const data_type &u0 = 0)
          {
            patches.clear();
            ukey.init(dus, due, u0);
            resize_store(ukey.key.size(), nv);
            uclosecache=UNKNOWN;
            vclosecache=UNKNOWN;
          }

          template<typename it__>
          void init_v(const it__ &dvs, const it__ &dve, const data_type &v0 = 0)
          {
            patches.clear();
            vkey.init(dvs, dve, v0);
            resize_store(nu, vkey.key.size());
            uclosecache=UNKNOWN;
            vclosecache=UNKNOWN;
          }

          template<typename it__>
          void init_uv(const it__ &dus, const it__ &due, const it__ &dvs, const it__ &dve, const data_type &u0 = 0, const data_type &v0 = 0)
          {
            patches.clear();
            ukey.init(dus, due, u0);
            vkey.init(dvs, dve, v0);
            resize_store(ukey.key.size(), vkey.key.size());
            uclosecache=UNKNOWN;
            vclosecache=UNKNOWN;
          }

          template<typename it__>
          void init_uv(const it__ &dus, const it__ &due, const index_type &nsegv, const data_type &dv = 1, const data_type &u0 = 0, const data_type &v0 = 0)
          {
            patches.clear();
            ukey.init(dus, due, u0);
            vkey.init(nsegv, dv, v0);
            resize_store(ukey.key.size(), vkey.key.size());
            uclosecache=UNKNOWN;
            vclosecache=UNKNOWN;
          }

          template<typename it__>
          void init_uv(const index_type &nsegu, const it__ &dvs, const it__ &dve, const data_type &du = 1, const data_type &u0 = 0, const data_type &v0 = 0)
          {
            patches.clear();
            ukey.init(nsegu, du, u0);
            vkey.init(dvs, dve, v0);
            resize_store(ukey.key.size(), vkey.key.size());
            uclosecache=UNKNOWN;
            vclosecache=UNKNOWN;
          }

          void init_uv( const std::vector<data_type> &umap, const std::vector<data_type> &vmap )
          {
            patches.clear();
            ukey.init( umap );
            vkey.init( vmap );
            resize_store(ukey.key.size(), vkey.key.size());
            uclosecache=UNKNOWN;
            vclosecache=UNKNOWN;
          }

          void degree_u(index_type &mind, index_type &maxd)
          {
            typename patch_collection_type::iterator uit;
            typename patch_strip_type::iterator vit;

            uit=patches.begin();
            vit=(*uit).begin();

            index_type d = vit->degree_u();
            mind = d;
            maxd = d;

            for (uit=patches.begin(); uit!=patches.end(); ++uit)
            {
              for (vit=(*uit).begin(); vit!=(*uit).end(); ++vit)
              {
                d = vit->degree_u();

                if(d<mind)
                {
                  mind = d;
                }
                if(d>maxd)
                {
                  maxd=d;
                }
              }
            }
          }

          void degree_u( std::vector<index_type> &deg)
          {
            index_type vk, j;
            typename keymap_type::const_iterator vit;

            deg.resize( nv );

            for ( j = 0, vit = vkey.key.begin(); vit != vkey.key.end(); ++vit, ++j )
            {
              vk = vit->second;

              deg[j] = patches[0][vk].degree_u();
            }
          }

          void degree_v(index_type &mind, index_type &maxd)
          {
            typename patch_collection_type::iterator uit;
            typename patch_strip_type::iterator vit;

            uit=patches.begin();
            vit=(*uit).begin();

            index_type d = vit->degree_v();
            mind = d;
            maxd = d;

            for (uit=patches.begin(); uit!=patches.end(); ++uit)
            {
              for (vit=(*uit).begin(); vit!=(*uit).end(); ++vit)
              {
                d = vit->degree_v();

                if(d<mind)
                {
                  mind = d;
                }
                if(d>maxd)
                {
                  maxd=d;
                }
              }
            }
          }

          void degree_v( std::vector<index_type> &deg)
          {
            index_type uk, i;
            typename keymap_type::const_iterator uit;

            deg.resize( nu );

            for ( i = 0, uit = ukey.key.begin(); uit != ukey.key.end(); ++uit, ++i )
            {
              uk = uit->second;

              deg[i] = patches[uk][0].degree_v();
            }
          }

          bool open_u() const
          {
            return !closed_u();
          }
          bool closed_u() const
          {
            if (uclosecache != UNKNOWN)
            {
              if (uclosecache == CLOSED)
                return true;

              return false;
            }

            index_type ifirst, ilast, j;
            typename surface_type::curve_type bc0, bc1;

            ifirst = ukey.key.begin()->second;
            ilast = ukey.key.rbegin()->second;

            for (j=0; j<nv; ++j)
            {
              patches[ifirst][j].get_umin_bndy_curve(bc0);
              patches[ilast][j].get_umax_bndy_curve(bc1);
              if (!eli::geom::curve::equivalent_curves(bc0, bc1))
              {
                uclosecache = OPEN;
                return false;
              }
            }

            uclosecache = CLOSED;
            return true;
          }

          bool open_v() const
          {
            return !closed_v();
          }
          bool closed_v() const
          {
            if (vclosecache != UNKNOWN)
            {
              if (vclosecache == CLOSED)
                return true;

              return false;
            }
            index_type i, jfirst, jlast;
            typename surface_type::curve_type bc0, bc1;

            jfirst = vkey.key.begin()->second;
            jlast = vkey.key.rbegin()->second;

            for (i=0; i<nu; ++i)
            {
              patches[i][jfirst].get_vmin_bndy_curve(bc0);
              patches[i][jlast].get_vmax_bndy_curve(bc1);
              if (!eli::geom::curve::equivalent_curves(bc0, bc1))
              {
                vclosecache = OPEN;
                return false;
              }
            }

            vclosecache = CLOSED;
            return true;
          }

          void get_bounding_box(bounding_box_type &bb) const
          {
            typename patch_collection_type::const_iterator uit;
            typename patch_strip_type::const_iterator vit;
            bounding_box_type bb_local;

            bb.clear();

            // cycle through all patches to get each bounding box to compare
            for (uit=patches.begin(); uit!=patches.end(); ++uit)
            {
              for (vit=(*uit).begin(); vit!=(*uit).end(); ++vit)
              {
                vit->get_bounding_box(bb_local);
                bb.add(bb_local);
              }
            }
          }

          void get_bounding_box( bounding_box_type &bb, const data_type U0, const data_type Uf, const data_type W0, const data_type Wf ) const
          {
              bounding_box_type bb_local;
              index_type uk, vk;
              typename keymap_type::const_iterator uit, vit;
              tolerance_type tol;

              if ( tol.approximately_less_than( U0, ukey.get_pmin() ) || tol.approximately_less_than( W0, vkey.get_pmin() ) ||
                   tol.approximately_less_than( ukey.get_pmax(), Uf ) || tol.approximately_less_than( vkey.get_pmax(), Wf ) )
                  return;

              // Create a temporary copy to use for splitting
              piecewise<surface__, data_type, dim__, tol__> s( *this );

              bb.clear();

              // Split to create patch borders at U/W limits. No split is performed if patch border at U/W value already exists
              s.split_u( U0 );
              s.split_u( Uf );
              s.split_v( W0 );
              s.split_v( Wf );

              // cycle through ukey and vkey within input limits to get each patch bounding box
              for ( uit = s.ukey.key.begin(); uit != s.ukey.key.end(); ++uit )
              {
                  if ( ( tol.approximately_equal( uit->first, U0 ) || tol.approximately_less_than( U0, uit->first ) ) && 
                       tol.approximately_less_than( uit->first, Uf ) )
                  {
                      uk = uit->second;
                      for ( vit = s.vkey.key.begin(); vit != s.vkey.key.end(); ++vit )
                      {
                          if ( ( tol.approximately_equal( vit->first, W0 ) || tol.approximately_less_than( W0, vit->first ) ) && 
                               tol.approximately_less_than( vit->first, Wf ) )
                          {
                              vk = vit->second;
                              s.patches[uk][vk].get_bounding_box( bb_local );
                              bb.add( bb_local );
                          }
                      }
                  }
              }
          }

          void rotate(const rotation_matrix_type &rmat)
          {
            typename patch_collection_type::iterator uit;
            typename patch_strip_type::iterator vit;

            for (uit=patches.begin(); uit!=patches.end(); ++uit)
            {
              for (vit=(*uit).begin(); vit!=(*uit).end(); ++vit)
              {
                vit->rotate(rmat);
              }
            }
          }

          void rotate(const rotation_matrix_type &rmat, const point_type &rorig)
          {
            typename patch_collection_type::iterator uit;
            typename patch_strip_type::iterator vit;

            for (uit=patches.begin(); uit!=patches.end(); ++uit)
            {
              for (vit=(*uit).begin(); vit!=(*uit).end(); ++vit)
              {
                vit->rotate(rmat, rorig);
              }
            }
          }

          void translate(const point_type &trans)
          {
            typename patch_collection_type::iterator uit;
            typename patch_strip_type::iterator vit;

            for (uit=patches.begin(); uit!=patches.end(); ++uit)
            {
              for (vit=(*uit).begin(); vit!=(*uit).end(); ++vit)
              {
                vit->translate(trans);
              }
            }
          }

          void reverse_u()
          {
            typename patch_collection_type::iterator uit;
            typename patch_strip_type::iterator vit;

            for (uit=patches.begin(); uit!=patches.end(); ++uit)
            {
              for (vit=(*uit).begin(); vit!=(*uit).end(); ++vit)
              {
                vit->reverse_u();
              }
            }
            ukey.reverse_keymap();
          }

          void reverse_v()
          {
            typename patch_collection_type::iterator uit;
            typename patch_strip_type::iterator vit;

            for (uit=patches.begin(); uit!=patches.end(); ++uit)
            {
              for (vit=(*uit).begin(); vit!=(*uit).end(); ++vit)
              {
                vit->reverse_v();
              }
            }
            vkey.reverse_keymap();
          }

          void roll_u( const index_type &index )
          {
            ukey.roll_keymap( index );
          }

          void roll_v( const index_type &index )
          {
            vkey.roll_keymap( index );
          }

          void scale(const data_type &s)
          {
            typename patch_collection_type::iterator uit;
            typename patch_strip_type::iterator vit;

            for (uit=patches.begin(); uit!=patches.end(); ++uit)
            {
              for (vit=(*uit).begin(); vit!=(*uit).end(); ++vit)
              {
                vit->scale(s);
              }
            }
          }

          void scale_x(const data_type &s)
          {
            typename patch_collection_type::iterator uit;
            typename patch_strip_type::iterator vit;

            for (uit=patches.begin(); uit!=patches.end(); ++uit)
            {
              for (vit=(*uit).begin(); vit!=(*uit).end(); ++vit)
              {
                vit->scale_x(s);
              }
            }
          }

          void scale_y(const data_type &s)
          {
            typename patch_collection_type::iterator uit;
            typename patch_strip_type::iterator vit;

            for (uit=patches.begin(); uit!=patches.end(); ++uit)
            {
              for (vit=(*uit).begin(); vit!=(*uit).end(); ++vit)
              {
                vit->scale_y(s);
              }
            }
          }

          void scale_z(const data_type &s)
          {
            typename patch_collection_type::iterator uit;
            typename patch_strip_type::iterator vit;

            for (uit=patches.begin(); uit!=patches.end(); ++uit)
            {
              for (vit=(*uit).begin(); vit!=(*uit).end(); ++vit)
              {
                vit->scale_z(s);
              }
            }
          }

          void swap_uv()
          {
            patch_collection_type old_patches;
            old_patches.swap(patches);

            index_type nu_old(nu), nv_old(nv);

            // Resizes patches and also assigns nu, nv.
            resize_store(nv_old, nu_old);

            for (index_type i=0; i<nu; ++i)
            {
              for (index_type j=0; j<nv; ++j)
              {
                patches[i][j]=old_patches[j][i];
                patches[i][j].swap_uv();
              }
            }

            data_type pmaxtmp;
            pmaxtmp = ukey.pmax;
            ukey.pmax = vkey.pmax;
            vkey.pmax = pmaxtmp;

            swap(ukey.key, vkey.key);

            index_type tmp( uclosecache );
            uclosecache = vclosecache;
            vclosecache = tmp;
          }

          void clear()
          {
            nu=0;
            nv=0;
            patches.clear();
            ukey.clear();
            vkey.clear();
            uclosecache=UNKNOWN;
            vclosecache=UNKNOWN;
          }

          error_code get(surface_type &surf, const index_type &ui, const index_type &vi) const
          {
            data_type du, dv;
            return get(surf, du, dv, ui, vi);
          }

          error_code get(surface_type &surf, data_type &du, data_type &dv, const index_type &ui, const index_type &vi) const
          {
            if ((ui>=number_u_patches()) || (vi>=number_v_patches()))
              return INVALID_INDEX;

            index_type uk, vk;
            typename keymap_type::const_iterator uit, vit;
            find_patch(uk, vk, uit, vit, ui, vi);

            du = ukey.get_delta_parm(uit);
            dv = vkey.get_delta_parm(vit);
            surf = patches[uk][vk];

            return NO_ERRORS;
          }

          data_type get_du( const index_type &ui ) const
          {
            if ( ui >= number_u_patches() )
              return INVALID_INDEX;

            index_type uk, vk;
            typename keymap_type::const_iterator uit, vit;
            find_patch( uk, vk, uit, vit, ui, 0 );

            return ukey.get_delta_parm( uit );
          }

          data_type get_dv( const index_type &vi ) const
          {
            if ( vi >= number_v_patches() )
              return INVALID_INDEX;

            index_type uk, vk;
            typename keymap_type::const_iterator uit, vit;
            find_patch(uk, vk, uit, vit, 0, vi);

            return vkey.get_delta_parm(vit);
          }

          error_code set(const surface_type &surf, const index_type &ui, const index_type &vi)
          {
            if ((ui>=number_u_patches()) || (vi>=number_v_patches()))
              return INVALID_INDEX;

            index_type uk, vk;
            typename keymap_type::const_iterator uit, vit;
            find_patch(uk, vk, uit, vit, ui, vi);

            // set the new surf
            patches[uk][vk]=surf;

            uclosecache=UNKNOWN;
            vclosecache=UNKNOWN;

            return NO_ERRORS;
          }

          error_code replace(const surface_type &surf, const index_type &ui, const index_type &vi)
          {
            if ((ui>=number_u_patches()) || (vi>=number_v_patches()))
              return INVALID_INDEX;

            // advance to desired index
            index_type uk, vk;
            typename keymap_type::const_iterator uit, vit;
            find_patch(uk, vk, uit, vit, ui, vi);

            // set the new surf
            patches[uk][vk]=surf;

            uclosecache=UNKNOWN;
            vclosecache=UNKNOWN;

            return NO_ERRORS;
          }

          error_code split_u(const data_type &u_in)
          {
            index_type uk, vk;
            typename keymap_type::iterator uit, vit;
            data_type uu(0), vv(0);
            data_type vmin = vkey.get_pmin();

            find_patch(uk, vk, uit, vit, uu, vv, u_in, vmin);

            // check for out of range input
            if ((uk == -1) || (vk == -1))
              return INVALID_PARAM;

            // check if no need to split
            tolerance_type tol;
            if (tol.approximately_equal(uu, 0))
              return NO_ERRORS;
            if (tol.approximately_equal(uu, 1))
              return NO_ERRORS;

            return split_u(uk, uit, u_in, uu);
          }

          error_code split_u(piecewise<surface__, data_type, dim__, tol__> &before, piecewise<surface__, data_type, dim__, tol__> &after, const data_type &u_in) const
          {
            before.clear();
            after.clear();

            if (u_in <= ukey.get_pmin())
            {
              after=(*this);
              return NO_ERRORS;
            }

            if (u_in >= ukey.get_pmax())
            {
              before=(*this);
              return NO_ERRORS;
            }

            piecewise<surface__, data_type, dim__, tol__> s(*this);

            error_code spliterr = s.split_u( u_in );

            typename keymap_type::const_iterator uit, vit;
            data_type uu(0), vv(0);
            data_type vmin = vkey.get_pmin();

            s.find_patch(uit, vit, uu, vv, u_in, vmin);

            s.subsurf(before, s.ukey.key.begin(), uit, s.vkey.key.begin(), s.vkey.key.end());
            s.subsurf(after, uit, s.ukey.key.end(), s.vkey.key.begin(), s.vkey.key.end());

            return spliterr;
          }

          error_code split_v(const data_type &v_in)
          {
            index_type uk, vk;
            typename keymap_type::iterator uit, vit;
            data_type uu(0), vv(0);
            data_type umin = ukey.get_pmin();

            find_patch(uk, vk, uit, vit, uu, vv, umin, v_in);

            // check for out of range input
            if ((uk == -1) || (vk == -1))
              return INVALID_PARAM;

            // check if no need to split
            tolerance_type tol;
            if (tol.approximately_equal(vv, 0))
              return NO_ERRORS;
            if (tol.approximately_equal(vv, 1))
              return NO_ERRORS;

            return split_v(vk, vit, v_in, vv);
          }

          error_code split_v(piecewise<surface__, data_type, dim__, tol__> &before, piecewise<surface__, data_type, dim__, tol__> &after, const data_type &v_in) const
          {
            before.clear();
            after.clear();

            if (v_in <= vkey.get_pmin())
            {
              after=(*this);
              return NO_ERRORS;
            }

            if (v_in >= vkey.get_pmax())
            {
              before=(*this);
              return NO_ERRORS;
            }

            piecewise<surface__, data_type, dim__, tol__> s(*this);

            error_code spliterr = s.split_v( v_in );

            typename keymap_type::const_iterator uit, vit;
            data_type uu(0), vv(0);
            data_type umin = ukey.get_pmin();

            s.find_patch(uit, vit, uu, vv, umin, v_in);

            s.subsurf(before, s.ukey.key.begin(), s.ukey.key.end(), s.vkey.key.begin(), vit);
            s.subsurf(after, s.ukey.key.begin(), s.ukey.key.end(), vit, s.vkey.key.end());

            return spliterr;
          }

          void join_u( const piecewise<surface__, data_type, dim__, tol__> &a, const piecewise<surface__, data_type, dim__, tol__> &b )
          {
            tolerance_type tol;

            std::vector<data_type> upmap1, upmap2, vpmap1, vpmap2, upmap, vpmap;
            a.get_pmap_uv( upmap1, vpmap1 );
            b.get_pmap_uv( upmap2, vpmap2 );

            index_type nu1, nv1, nu2, nv2;
            nu1 = a.number_u_patches();
            nv1 = a.number_v_patches();
            nu2 = b.number_u_patches();
            nv2 = b.number_v_patches();

            upmap.resize( nu1 + nu2 + 1 );

            // Append U parameter values
            index_type j = 0;
            for ( index_type i = 0; i < nu1 + 1; i++, j++ )
            {
              upmap[j] = upmap1[i];
            }
            for ( index_type i = 0; i < nu2; i++, j++ )
            {
              data_type du = upmap2[i+1] - upmap2[i];
              upmap[j] = upmap[j-1] + du;
            }

            // V should be identical.
            vpmap = vpmap1;

            init_uv( upmap, vpmap );

            index_type ju = 0;
            for ( index_type iu = 0; iu < nu1; iu++, ju++ )
            {
              for ( index_type iv = 0; iv < nv1; iv++ )
              {
                set( *(a.get_patch( iu, iv )), ju, iv );
              }
            }

            for ( index_type iu = 0; iu < nu2; iu++, ju++ )
            {
              for ( index_type iv = 0; iv < nv2; iv++ )
              {
                set( *(b.get_patch( iu, iv )), ju, iv );
              }
            }
          }

          void join_v( const piecewise<surface__, data_type, dim__, tol__> &a, const piecewise<surface__, data_type, dim__, tol__> &b )
          {
            tolerance_type tol;

            std::vector<data_type> upmap1, upmap2, vpmap1, vpmap2, upmap, vpmap;
            a.get_pmap_uv( upmap1, vpmap1 );
            b.get_pmap_uv( upmap2, vpmap2 );

            index_type nu1, nv1, nu2, nv2;
            nu1 = a.number_u_patches();
            nv1 = a.number_v_patches();
            nu2 = b.number_u_patches();
            nv2 = b.number_v_patches();

            vpmap.resize( nv1 + nv2 + 1 );

            // Append U parameter values
            index_type j = 0;
            for ( index_type i = 0; i < nv1 + 1; i++, j++ )
            {
              vpmap[j] = vpmap1[i];
            }
            for ( index_type i = 0; i < nv2; i++, j++ )
            {
              data_type dv = vpmap2[i+1] - vpmap2[i];
              vpmap[j] = vpmap[j-1] + dv;
            }

            // V should be identical.
            upmap = upmap1;

            init_uv( upmap, vpmap );


            for ( index_type iu = 0; iu < nu1; iu++ )
            {
              index_type jv = 0;
              for ( index_type iv = 0; iv < nv1; iv++, jv++ )
              {
                set( *(a.get_patch( iu, iv )), iu, jv );
              }

              for ( index_type iv = 0; iv < nv2; iv++, jv++ )
              {
                set( *(b.get_patch( iu, iv )), iu, jv );
              }
            }

          }

          void to_cubic_u(const data_type &ttol)
          {
            typename keymap_type::iterator uit, vit;

            // First pass to split patches until cubic approximation is within tolerance.
            for(uit = ukey.key.begin(); uit != ukey.key.end(); ++uit)
            {
              for(vit = vkey.key.begin(); vit != vkey.key.end(); ++vit)
              {
                index_type uk = uit->second;
                index_type vk = vit->second;

                surface_type sc = patches[uk][vk];

                sc.to_cubic_u();

                data_type d = patches[uk][vk].eqp_distance_bound(sc);

                while(d > ttol)
                {
                  data_type delta_u = ukey.get_delta_parm(uit);
                  data_type u_in = uit->first + static_cast<data_type>(0.5) * delta_u;

                  split_u(uk, uit, u_in, 0.5);

                  sc = patches[uk][vk];

                  sc.to_cubic_u();

                  d = patches[uk][vk].eqp_distance_bound(sc);
                }
              }
            }

            // Second pass to convert all patches to cubic.
            for (index_type uk=0; uk<nu; ++uk)
            {
              for (index_type vk=0; vk<nv; ++vk)
              {
                patches[uk][vk].to_cubic_u();
              }
            }
          }

          void to_cubic_v(const data_type &ttol)
          {
            typename keymap_type::iterator uit, vit;

            // First pass to split patches until cubic approximation is within tolerance.
            for(uit = ukey.key.begin(); uit != ukey.key.end(); ++uit)
            {
              for(vit = vkey.key.begin(); vit != vkey.key.end(); ++vit)
              {
                index_type uk = uit->second;
                index_type vk = vit->second;

                surface_type sc = patches[uk][vk];

                sc.to_cubic_v();

                data_type d = patches[uk][vk].eqp_distance_bound(sc);

                while(d > ttol)
                {
                  data_type delta_v = vkey.get_delta_parm(vit);
                  data_type v_in = vit->first + static_cast<data_type>(0.5) * delta_v;

                  split_v(vk, vit, v_in, 0.5);

                  sc = patches[uk][vk];

                  sc.to_cubic_v();

                  d = patches[uk][vk].eqp_distance_bound(sc);
                }
              }
            }

            // Second pass to convert all patches to cubic.
            for (index_type uk=0; uk<nu; ++uk)
            {
              for (index_type vk=0; vk<nv; ++vk)
              {
                patches[uk][vk].to_cubic_v();
              }
            }
          }

          void to_cubic(const data_type &ttol)
          {
            to_cubic_u(ttol);
            to_cubic_v(ttol);
          }

          void promote_u_to( const std::vector< index_type > ord )
          {
            assert ( ord.size() == nu );

            index_type uk;
            typename keymap_type::const_iterator uit;

            index_type i;

            for ( i = 0, uit = ukey.key.begin(); uit != ukey.key.end(); ++uit, ++i )
            {
              uk = uit->second;

              for (index_type j=0; j<nv; ++j)
              {
                patches[uk][j].promote_u_to( ord[i] );
              }
            }
          }

          void promote_v_to( const std::vector< index_type > ord )
          {
            assert ( ord.size() == nv );

            index_type vk;
            typename keymap_type::const_iterator vit;

            index_type j;

            for ( j = 0, vit = vkey.key.begin(); vit != vkey.key.end(); ++vit, ++j )
            {
              vk = vit->second;

              for (index_type i=0; i<nu; ++i)
              {
                patches[i][vk].promote_v_to( ord[j] );
              }
            }
          }

          void get_uconst_curve(piecewise_curve_type &pwc, const data_type &u) const
          {
            index_type uk, vk;
            typename keymap_type::const_iterator uit, vit;
            data_type uu(0), vv(0);
            data_type vmin = vkey.get_pmin();

            find_patch(uk, vk, uit, vit, uu, vv, u, vmin);

            assert ((uk != -1) && (vk != -1));

            pwc.clear();
            pwc.set_t0(vmin);

            for ( vit = vkey.key.begin(); vit != vkey.key.end(); ++vit )
            {
              vk = vit->second;

              data_type dv=vkey.get_delta_parm(vit);

              curve_type c;

              patches[uk][vk].get_uconst_curve(c, uu);

              pwc.push_back(c,dv);
            }
          }

          void get_umin_bndy_curve( piecewise_curve_type &pwc ) const
          {
            index_type uk, vk;
            typename keymap_type::const_iterator vit;
            data_type vmin = vkey.get_pmin();

            uk = ukey.key.begin()->second;

            pwc.clear();
            pwc.set_t0(vmin);

            for ( vit = vkey.key.begin(); vit != vkey.key.end(); ++vit )
            {
              vk = vit->second;

              data_type dv=vkey.get_delta_parm(vit);

              curve_type c;

              patches[uk][vk].get_umin_bndy_curve(c);

              pwc.push_back(c,dv);
            }
          }

          void get_umax_bndy_curve( piecewise_curve_type &pwc ) const
          {
            index_type uk, vk;
            typename keymap_type::const_iterator vit;
            data_type vmin = vkey.get_pmin();

            uk = ukey.key.rbegin()->second;

            pwc.clear();
            pwc.set_t0(vmin);

            for ( vit = vkey.key.begin(); vit != vkey.key.end(); ++vit )
            {
              vk = vit->second;

              data_type dv=vkey.get_delta_parm(vit);

              curve_type c;

              patches[uk][vk].get_umax_bndy_curve(c);

              pwc.push_back(c,dv);
            }
          }

          void get_vmin_bndy_curve( piecewise_curve_type &pwc ) const
          {
            index_type uk, vk;
            typename keymap_type::const_iterator uit;
            data_type umin = ukey.get_pmin();

            vk = vkey.key.begin()->second;

            pwc.clear();
            pwc.set_t0(umin);

            for ( uit = ukey.key.begin(); uit != ukey.key.end(); ++uit )
            {
              uk = uit->second;

              data_type du=ukey.get_delta_parm(uit);

              curve_type c;

              patches[uk][vk].get_vmin_bndy_curve(c);

              pwc.push_back(c,du);
            }
          }

          void get_vmax_bndy_curve( piecewise_curve_type &pwc ) const
          {
            index_type uk, vk;
            typename keymap_type::const_iterator uit;
            data_type umin = ukey.get_pmin();

            vk = vkey.key.rbegin()->second;

            pwc.clear();
            pwc.set_t0(umin);

            for ( uit = ukey.key.begin(); uit != ukey.key.end(); ++uit )
            {
              uk = uit->second;

              data_type du=ukey.get_delta_parm(uit);

              curve_type c;

              patches[uk][vk].get_vmax_bndy_curve(c);

              pwc.push_back(c,du);
            }
          }

          void get_umin_ndelta_pcurve( piecewise_curve_type &pwc ) const
          {
            index_type uk, vk;
            typename keymap_type::const_iterator vit;
            data_type vmin = vkey.get_pmin();

            uk = ukey.key.begin()->second;

            pwc.clear();
            pwc.set_t0(vmin);

            for ( vit = vkey.key.begin(); vit != vkey.key.end(); ++vit )
            {
              vk = vit->second;

              data_type dv=vkey.get_delta_parm(vit);

              curve_type c;

              patches[uk][vk].get_umin_ndelta_pcurve(c);

              pwc.push_back(c,dv);
            }
          }

          void get_umax_ndelta_pcurve( piecewise_curve_type &pwc ) const
          {
            index_type uk, vk;
            typename keymap_type::const_iterator vit;
            data_type vmin = vkey.get_pmin();

            uk = ukey.key.rbegin()->second;

            pwc.clear();
            pwc.set_t0(vmin);

            for ( vit = vkey.key.begin(); vit != vkey.key.end(); ++vit )
            {
              vk = vit->second;

              data_type dv=vkey.get_delta_parm(vit);

              curve_type c;

              patches[uk][vk].get_umax_ndelta_pcurve(c);

              pwc.push_back(c,dv);
            }
          }

          void get_vmin_ndelta_pcurve( piecewise_curve_type &pwc ) const
          {
            index_type uk, vk;
            typename keymap_type::const_iterator uit;
            data_type umin = ukey.get_pmin();

            vk = vkey.key.begin()->second;

            pwc.clear();
            pwc.set_t0(umin);

            for ( uit = ukey.key.begin(); uit != ukey.key.end(); ++uit )
            {
              uk = uit->second;

              data_type du=ukey.get_delta_parm(uit);

              curve_type c;

              patches[uk][vk].get_vmin_ndelta_pcurve(c);

              pwc.push_back(c,du);
            }
          }

          void get_vmax_ndelta_pcurve( piecewise_curve_type &pwc ) const
          {
            index_type uk, vk;
            typename keymap_type::const_iterator uit;
            data_type umin = ukey.get_pmin();

            vk = vkey.key.rbegin()->second;

            pwc.clear();
            pwc.set_t0(umin);

            for ( uit = ukey.key.begin(); uit != ukey.key.end(); ++uit )
            {
              uk = uit->second;

              data_type du=ukey.get_delta_parm(uit);

              curve_type c;

              patches[uk][vk].get_vmax_ndelta_pcurve(c);

              pwc.push_back(c,du);
            }
          }

          void get_uconst_f_u_curve(piecewise_curve_type &pwc, const data_type &u) const
          {
            index_type uk, vk;
            typename keymap_type::const_iterator uit, vit;
            data_type uu(0), vv(0);
            data_type vmin = vkey.get_pmin();

            find_patch(uk, vk, uit, vit, uu, vv, u, vmin);

            assert ((uk != -1) && (vk != -1));

            pwc.clear();
            pwc.set_t0(vmin);

            for ( vit = vkey.key.begin(); vit != vkey.key.end(); ++vit )
            {
              vk = vit->second;

              data_type dv=vkey.get_delta_parm(vit);

              curve_type c;

              patches[uk][vk].get_uconst_f_u_curve(c, uu);

              pwc.push_back(c,dv);
            }
          }

          void get_uconst_f_v_curve(piecewise_curve_type &pwc, const data_type &u) const
          {
            index_type uk, vk;
            typename keymap_type::const_iterator uit, vit;
            data_type uu(0), vv(0);
            data_type vmin = vkey.get_pmin();

            find_patch(uk, vk, uit, vit, uu, vv, u, vmin);

            assert ((uk != -1) && (vk != -1));

            pwc.clear();
            pwc.set_t0(vmin);

            for ( vit = vkey.key.begin(); vit != vkey.key.end(); ++vit )
            {
              vk = vit->second;

              data_type dv=vkey.get_delta_parm(vit);

              curve_type c;

              patches[uk][vk].get_uconst_f_v_curve(c, uu);

              pwc.push_back(c,dv);
            }
          }

          void get_vconst_curve(piecewise_curve_type &pwc, const data_type &v) const
          {
            index_type uk, vk;
            typename keymap_type::const_iterator uit, vit;
            data_type uu(0), vv(0);
            data_type umin = ukey.get_pmin();

            find_patch(uk, vk, uit, vit, uu, vv, umin, v);

            assert ((uk != -1) && (vk != -1));

            pwc.clear();
            pwc.set_t0(umin);

            for ( uit = ukey.key.begin(); uit != ukey.key.end(); ++uit )
            {
              uk = uit->second;

              data_type du=ukey.get_delta_parm(uit);

              curve_type c;

              patches[uk][vk].get_vconst_curve(c, vv);

              pwc.push_back(c,du);
            }
          }

          void get_vconst_f_u_curve(piecewise_curve_type &pwc, const data_type &v) const
          {
            index_type uk, vk;
            typename keymap_type::const_iterator uit, vit;
            data_type uu(0), vv(0);
            data_type umin = ukey.get_pmin();

            find_patch(uk, vk, uit, vit, uu, vv, umin, v);

            assert ((uk != -1) && (vk != -1));

            pwc.clear();
            pwc.set_t0(umin);

            for ( uit = ukey.key.begin(); uit != ukey.key.end(); ++uit )
            {
              uk = uit->second;

              data_type du=ukey.get_delta_parm(uit);

              curve_type c;

              patches[uk][vk].get_vconst_f_u_curve(c, vv);

              pwc.push_back(c,du);
            }
          }

          void get_vconst_f_v_curve(piecewise_curve_type &pwc, const data_type &v) const
          {
            index_type uk, vk;
            typename keymap_type::const_iterator uit, vit;
            data_type uu(0), vv(0);
            data_type umin = ukey.get_pmin();

            find_patch(uk, vk, uit, vit, uu, vv, umin, v);

            assert ((uk != -1) && (vk != -1));

            pwc.clear();
            pwc.set_t0(umin);

            for ( uit = ukey.key.begin(); uit != ukey.key.end(); ++uit )
            {
              uk = uit->second;

              data_type du=ukey.get_delta_parm(uit);

              curve_type c;

              patches[uk][vk].get_vconst_f_v_curve(c, vv);

              pwc.push_back(c,du);
            }
          }
          void find_interior_feature_edges(std::vector<data_type> &uconst, std::vector<data_type> &vconst, const data_type &angle_tol) const
          {
            index_type nu, nv, iu, iv;

            piecewise_curve_type c;
            std::vector<data_type> pmap, ldis, ldis_out;
            tolerance_type tol;

            // initialize the output
            uconst.clear();
            vconst.clear();

            // unnamed function to test if two parameters are close enough
            auto comp = [&tol](const data_type &x1, const data_type &x2)->bool
            {
              return tol.approximately_less_than(x1, x2);
            };

            // extract each v-const curve that is an edge of one of the patches to find u-parameters
            // that contain C0 only edges
            get_pmap_v(pmap);
            nv=pmap.size();
            assert(nv-1==number_v_patches());
            for (iv=0; iv<nv; ++iv)
            {
              get_vconst_curve(c, pmap[iv]);
              c.find_discontinuities(angle_tol, ldis);

              // merge these parameters with current list
              ldis_out.clear();
              std::set_union(uconst.begin(), uconst.end(), ldis.begin(), ldis.end(), std::back_inserter(ldis_out), comp);
              std::swap(uconst, ldis_out);
            }

            // extract each u-const curve that is an edge of one of the patches to find v-parameters
            // that contain C0 only edges
            pmap.clear();
            get_pmap_u(pmap);
            nu=pmap.size();
            assert(nu-1==number_u_patches());
            for (iu=0; iu<nu; ++iu)
            {
              get_uconst_curve(c, pmap[iu]);
              c.find_discontinuities(angle_tol, ldis);

              // merge these parameters with current list
              ldis_out.clear();
              std::set_union(vconst.begin(), vconst.end(), ldis.begin(), ldis.end(), std::back_inserter(ldis_out), comp);
              std::swap(vconst, ldis_out);
            }

            // TODO: Need to compare actual control points next to edges to catch cases where the
            //       patch corners satisfy the constraints but internally the constraints are not.

          }

          void find_interior_C0_edges(std::vector<data_type> &uconst, std::vector<data_type> &vconst) const
          {
            index_type nu, nv, iu, iv;

            piecewise_curve_type c;
            std::vector<data_type> pmap, ldis, ldis_out;
            tolerance_type tol;

            // initialize the output
            uconst.clear();
            vconst.clear();

            // unnamed function to test if two parameters are close enough
            auto comp = [&tol](const data_type &x1, const data_type &x2)->bool
            {
              return tol.approximately_less_than(x1, x2);
            };

            // extract each v-const curve that is an edge of one of the patches to find u-parameters
            // that contain C0 only edges
            get_pmap_v(pmap);
            nv=pmap.size();
            assert(nv-1==number_v_patches());
            for (iv=0; iv<nv; ++iv)
            {
              get_vconst_curve(c, pmap[iv]);
              c.find_discontinuities(eli::geom::general::G1, ldis);

              // merge these parameters with current list
              ldis_out.clear();
              std::set_union(uconst.begin(), uconst.end(), ldis.begin(), ldis.end(), std::back_inserter(ldis_out), comp);
              std::swap(uconst, ldis_out);
            }

            // extract each u-const curve that is an edge of one of the patches to find v-parameters
            // that contain C0 only edges
            pmap.clear();
            get_pmap_u(pmap);
            nu=pmap.size();
            assert(nu-1==number_u_patches());
            for (iu=0; iu<nu; ++iu)
            {
              get_uconst_curve(c, pmap[iu]);
              c.find_discontinuities(eli::geom::general::G1, ldis);

              // merge these parameters with current list
              ldis_out.clear();
              std::set_union(vconst.begin(), vconst.end(), ldis.begin(), ldis.end(), std::back_inserter(ldis_out), comp);
              std::swap(vconst, ldis_out);
            }

            // TODO: Need to compare actual control points next to edges to catch cases where the
            //       patch corners are continuous but internally it is not.
          }

          point_type f(const data_type &u, const data_type &v) const
          {
            // find patch that corresponds to given u & v
            index_type uk, vk;
            data_type uu(0), vv(0);

            find_patch(uk, vk, uu, vv, u, v);

            assert((uk != -1) && (vk != -1));

            return patches[uk][vk].f(uu, vv);
          }

          point_type f_u(const data_type &u, const data_type &v) const
          {
            // find patch that corresponds to given u & v
            index_type uk, vk;
            typename keymap_type::const_iterator uit, vit;
            data_type uu(0), vv(0);

            find_patch(uk, vk, uit, vit, uu, vv, u, v);

            assert((uk != -1) && (vk != -1));

            data_type delta_u = ukey.get_delta_parm(uit);

            return patches[uk][vk].f_u(uu, vv)/delta_u;
          }

          point_type f_v(const data_type &u, const data_type &v) const
          {
            // find patch that corresponds to given u & v
            index_type uk, vk;
            typename keymap_type::const_iterator uit, vit;
            data_type uu(0), vv(0);

            find_patch(uk, vk, uit, vit, uu, vv, u, v);

            assert((uk != -1) && (vk != -1));

            data_type delta_v = vkey.get_delta_parm(vit);

            return patches[uk][vk].f_v(uu, vv)/delta_v;
          }

          point_type f_uu(const data_type &u, const data_type &v) const
          {
            // find patch that corresponds to given u & v
            index_type uk, vk;
            typename keymap_type::const_iterator uit, vit;
            data_type uu(0), vv(0);

            find_patch(uk, vk, uit, vit, uu, vv, u, v);

            assert((uk != -1) && (vk != -1));

            data_type delta_u = ukey.get_delta_parm(uit);

            return patches[uk][vk].f_uu(uu, vv)/(delta_u*delta_u);
          }

          point_type f_uv(const data_type &u, const data_type &v) const
          {
            // find patch that corresponds to given u & v
            index_type uk, vk;
            typename keymap_type::const_iterator uit, vit;
            data_type uu(0), vv(0);

            find_patch(uk, vk, uit, vit, uu, vv, u, v);

            assert((uk != -1) && (vk != -1));

            data_type delta_u = ukey.get_delta_parm(uit);
            data_type delta_v = vkey.get_delta_parm(vit);

            return patches[uk][vk].f_uv(uu, vv)/(delta_u*delta_v);
          }

          point_type f_vv(const data_type &u, const data_type &v) const
          {
            // find patch that corresponds to given u & v
            index_type uk, vk;
            typename keymap_type::const_iterator uit, vit;
            data_type uu(0), vv(0);

            find_patch(uk, vk, uit, vit, uu, vv, u, v);

            assert((uk != -1) && (vk != -1));

            data_type delta_v = vkey.get_delta_parm(vit);

            return patches[uk][vk].f_vv(uu, vv)/(delta_v*delta_v);
          }

          point_type f_uuu(const data_type &u, const data_type &v) const
          {
            // find patch that corresponds to given u & v
            index_type uk, vk;
            typename keymap_type::const_iterator uit, vit;
            data_type uu(0), vv(0);

            find_patch(uk, vk, uit, vit, uu, vv, u, v);

            assert((uk != -1) && (vk != -1));

            data_type delta_u = ukey.get_delta_parm(uit);

            return patches[uk][vk].f_uuu(uu, vv)/(delta_u*delta_u*delta_u);
          }

          point_type f_uuv(const data_type &u, const data_type &v) const
          {
            // find patch that corresponds to given u & v
            index_type uk, vk;
            typename keymap_type::const_iterator uit, vit;
            data_type uu(0), vv(0);

            find_patch(uk, vk, uit, vit, uu, vv, u, v);

            assert((uk != -1) && (vk != -1));

            data_type delta_u = ukey.get_delta_parm(uit);
            data_type delta_v = vkey.get_delta_parm(vit);

            return patches[uk][vk].f_uuv(uu, vv)/(delta_u*delta_u*delta_v);
          }

          point_type f_uvv(const data_type &u, const data_type &v) const
          {
            // find patch that corresponds to given u & v
            index_type uk, vk;
            typename keymap_type::const_iterator uit, vit;
            data_type uu(0), vv(0);

            find_patch(uk, vk, uit, vit, uu, vv, u, v);

            assert((uk != -1) && (vk != -1));

            data_type delta_u = ukey.get_delta_parm(uit);
            data_type delta_v = vkey.get_delta_parm(vit);

            return patches[uk][vk].f_uvv(uu, vv)/(delta_u*delta_v*delta_v);
          }

          point_type f_vvv(const data_type &u, const data_type &v) const
          {
            // find patch that corresponds to given u & v
            index_type uk, vk;
            typename keymap_type::const_iterator uit, vit;
            data_type uu(0), vv(0);

            find_patch(uk, vk, uit, vit, uu, vv, u, v);

            assert((uk != -1) && (vk != -1));

            data_type delta_u = ukey.get_delta_parm(uit);
            data_type delta_v = vkey.get_delta_parm(vit);

            return patches[uk][vk].f_vvv(uu, vv)/(delta_v*delta_v*delta_v);
          }

          point_type normal(const data_type &u, const data_type &v) const
          {
            // find patch that corresponds to given u & v
            index_type uk, vk;
            data_type uu(0), vv(0);

            find_patch(uk, vk, uu, vv, u, v);

            assert((uk != -1) && (vk != -1));

            return patches[uk][vk].normal(uu, vv);
          }

          void f_pt_normal(const data_type &u, const data_type &v, point_type &pt, point_type &norm, const index_type &utie = 0, const index_type &vtie = 0 ) const
          {
            // find patch that corresponds to given u & v
            index_type uk, vk;
            data_type uu(0), vv(0);
            typename keymap_type::const_iterator uit, vit;

            patch_boundary_code_type bcode = find_patch(uk, vk, uit, vit, uu, vv, u, v);

            assert((uk != -1) && (vk != -1));

            // Handle u boundary code
            if ( bcode.first == -1 && utie == 1 ) // u is at start of interval, uu = 0.
            {
              if ( uit != ukey.key.begin() )
              {
                uit--;
                uk = uit->second;
                uu = 1.0;
              }
            }
            else if ( bcode.first == 1 && utie == -1 ) // u is at end of interval uu = 1;
            {
              typename keymap_type::const_iterator uitnext = uit;
              uitnext++;

              if ( uitnext != ukey.key.end() )
              {
                uit = uitnext;
                uk = uit->second;
                uu = 0.0;
              }
            }

            // Handle v boundary code
            if ( bcode.second == -1 && vtie == 1 ) // v is at start of interval, vv = 0.
            {
              if ( vit != vkey.key.begin() )
              {
                vit--;
                vk = vit->second;
                vv = 1.0;
              }
            }
            else if ( bcode.second == 1 && vtie == -1) // v is at end of interval vv = 1;
            {
              typename keymap_type::const_iterator vitnext = vit;
              vitnext++;

              if ( vitnext != vkey.key.end() )
              {
                vit = vitnext;
                vk = vit->second;
                vv = 0.0;
              }
            }

            pt = patches[uk][vk].f(uu, vv);
            norm = patches[uk][vk].normal(uu, vv);
          }

          void f_pt_normal_grid(const std::vector < data_type > &uvec, const std::vector < data_type > &vvec, std::vector < std::vector < point_type > > &ptmat, std::vector < std::vector < point_type > > &normmat ) const
          {
            typedef std::vector < data_type > pvec_type;

            index_type nuvec( uvec.size() );
            index_type nvvec( vvec.size() );

            pvec_type uuvec( nuvec );
            pvec_type vvvec( nvvec );
            std::vector < index_type > ukvec( nuvec );
            std::vector < index_type > vkvec( nvvec );

            int i, j;

            std::vector < index_type > ukbatch( nu, 0 );
            std::vector < index_type > vkbatch( nv, 0 );

            for ( i = 0; i < nuvec; i++ )
            {
              typename keymap_type::const_iterator uit;
              index_type uk;
              data_type uu;
              index_type code = ukey.find_segment( uk, uit, uu, uvec[i] );

              if ( code == -1 && i == ( uvec.size() - 1 ) ) // u is at start of interval, uu = 0.
              {
                if ( uit != ukey.key.begin() )
                {
                  uit--;
                  uk = uit->second;
                  uu = 1.0;
                }
              }
              else if ( code == 1 && i == 0 ) // u is at end of interval uu = 1;
              {
                typename keymap_type::const_iterator uitnext = uit;
                uitnext++;

                if ( uitnext != ukey.key.end() )
                {
                  uit = uitnext;
                  uk = uit->second;
                  uu = 0.0;
                }
              }

              uuvec[i] = uu;
              ukvec[i] = uk;
              ukbatch[uk] += 1;
            }

            for ( i = 0; i < nvvec; i++ )
            {
              typename keymap_type::const_iterator vit;
              index_type vk;
              data_type vv;
              index_type code = vkey.find_segment( vk, vit, vv, vvec[i] );

              if ( code == -1 && i == ( vvec.size() - 1 ) ) // v is at start of interval, vv = 0.
              {
                if ( vit != vkey.key.begin() )
                {
                  vit--;
                  vk = vit->second;
                  vv = 1.0;
                }
              }
              else if ( code == 1 && i == 0 ) // v is at end of interval vv = 1;
              {
                typename keymap_type::const_iterator vitnext = vit;
                vitnext++;

                if ( vitnext != vkey.key.end() )
                {
                  vit = vitnext;
                  vk = vit->second;
                  vv = 0.0;
                }
              }

              vvvec[i] = vv;
              vkvec[i] = vk;
              vkbatch[vk] += 1;
            }

            std::vector < std::vector < point_type > > S_u_mat;
            std::vector < std::vector < point_type > > S_v_mat;

            ptmat.resize( nuvec );
            normmat.resize( nuvec );
            S_u_mat.resize( nuvec );
            S_v_mat.resize( nuvec );
            for ( i = 0; i < nuvec; i++ )
            {
              ptmat[i].resize( nvvec );
              normmat[i].resize( nvvec );
              S_u_mat[i].resize( nvvec );
              S_v_mat[i].resize( nvvec );
            }

            for ( i = 0; i < nuvec; )
            {
              index_type uk = ukvec[i];
              index_type nuk = ukbatch[uk];

              for ( j = 0; j < nvvec; )
              {
                index_type vk = vkvec[j];
                index_type nvk = vkbatch[vk];

                patches[uk][vk].fbatch( i, j, nuk, nvk, uuvec, vvvec, ptmat );
                patches[uk][vk].normalbatch( i, j, nuk, nvk, uuvec, vvvec, S_u_mat, S_v_mat, normmat );

                j += nvk;
              }
              i += nuk;
            }
          }

          void f_pt_derivs(const data_type &u, const data_type &v, point_type &pt, point_type &pt_u, point_type &pt_v ) const
          {
            // find patch that corresponds to given u & v
            index_type uk, vk;
            typename keymap_type::const_iterator uit, vit;
            data_type uu(0), vv(0);

            find_patch(uk, vk, uit, vit, uu, vv, u, v);

            assert((uk != -1) && (vk != -1));

            pt = patches[uk][vk].f(uu, vv);

            data_type delta_u = ukey.get_delta_parm(uit);

            pt_u = patches[uk][vk].f_u(uu, vv)/delta_u;

            data_type delta_v = vkey.get_delta_parm(vit);

            pt_v = patches[uk][vk].f_v(uu, vv)/delta_v;
          }

          static void order_match_u( piecewise<surface__, data_type, dim__, tol__> &s1, piecewise<surface__, data_type, dim__, tol__> &s2 )
          {
            index_type nu1 = s1.number_u_patches();
            index_type nu2 = s2.number_u_patches();

            assert ( nu1 == nu2 );

            index_type nv1 = s1.number_v_patches();
            index_type nv2 = s2.number_v_patches();

            std::vector < index_type > deg1, deg2;
            s1.degree_u( deg1 );
            s2.degree_u( deg2 );

            for ( index_type i = 0; i < nu1; i++ )
            {
              if ( deg2[i] > deg1[i] )
              {
                for ( index_type j = 0; j < nv1; j++ )
                {
                  s1.get_patch( i, j )->promote_u_to( deg2[i] );
                }
              }
              else if ( deg1[i] > deg2[i] )
              {
                for ( index_type j = 0; j < nv2; j++ )
                {
                  s2.get_patch( i, j )->promote_u_to( deg1[i] );
                }
              }
            }
          }

          static void order_match_v( piecewise<surface__, data_type, dim__, tol__> &s1, piecewise<surface__, data_type, dim__, tol__> &s2 )
          {
            index_type nu1 = s1.number_u_patches();
            index_type nu2 = s2.number_u_patches();

            index_type nv1 = s1.number_v_patches();
            index_type nv2 = s2.number_v_patches();

            assert ( nv1 == nv2 );

            std::vector < index_type > deg1, deg2;
            s1.degree_v( deg1 );
            s2.degree_v( deg2 );

            for ( index_type j = 0; j < nv1; j++ )
            {
              if ( deg2[j] > deg1[j] )
              {
                for ( index_type i = 0; i < nu1; i++ )
                {
                  s1.get_patch( i, j )->promote_v_to( deg2[j] );
                }
              }
              else if ( deg1[j] > deg2[j] )
              {
                for ( index_type i = 0; i < nu2; i++ )
                {
                  s2.get_patch( i, j )->promote_v_to( deg1[j] );
                }
              }
            }
          }

          static void parm_match_u( piecewise<surface__, data_type, dim__, tol__> &s1, piecewise<surface__, data_type, dim__, tol__> &s2 )
          {
            std::vector<data_type> upmap1, upmap2, vpmap1, vpmap2, upmap, vpmap;
            s1.get_pmap_uv( upmap1, vpmap1 );
            s2.get_pmap_uv( upmap2, vpmap2 );

            tolerance_type ttol;
            // Comparison function for set_union.
            auto comp = [&ttol](const data_type &x1, const data_type &x2)->bool
            {
              return ttol.approximately_less_than(x1, x2);
            };

            // Place union of 1 and 2 into pmaps
            std::set_union( upmap1.begin(), upmap1.end(), upmap2.begin(), upmap2.end(), std::back_inserter(upmap), comp );

            for ( index_type i = 0; i < upmap.size(); i++ )
            {
              s1.split_u( upmap[i] );
              s2.split_u( upmap[i] );
            }
          }

          static void parm_match_v( piecewise<surface__, data_type, dim__, tol__> &s1, piecewise<surface__, data_type, dim__, tol__> &s2 )
          {
            std::vector<data_type> upmap1, upmap2, vpmap1, vpmap2, upmap, vpmap;
            s1.get_pmap_uv( upmap1, vpmap1 );
            s2.get_pmap_uv( upmap2, vpmap2 );

            tolerance_type ttol;
            // Comparison function for set_union.
            auto comp = [&ttol](const data_type &x1, const data_type &x2)->bool
            {
              return ttol.approximately_less_than(x1, x2);
            };

            // Place union of 1 and 2 into pmaps
            std::set_union( vpmap1.begin(), vpmap1.end(), vpmap2.begin(), vpmap2.end(), std::back_inserter(vpmap), comp );

            for ( index_type i = 0; i < vpmap.size(); i++ )
            {
              s1.split_v( vpmap[i] );
              s2.split_v( vpmap[i] );
            }
          }

          void sum( const piecewise<surface__, data_type, dim__, tol__> &a, const piecewise<surface__, data_type, dim__, tol__> &b )
          {
            typedef piecewise<surface__, data_type, dim__, tol__> piecewise_surf_type;
            tolerance_type tol;

            piecewise_surf_type s1(a);
            piecewise_surf_type s2(b);

            std::vector<data_type> upmap, vpmap;

            parm_match_u( s1, s2 );
            parm_match_v( s1, s2 );

            s1.get_pmap_uv( upmap, vpmap );

            init_uv( upmap, vpmap );

            for ( index_type iu = 0; iu < nu; iu++ )
            {
              for ( index_type iv = 0; iv < nv; iv++ )
              {
                surface_type *p1 = s1.get_patch( iu, iv );
                surface_type *p2 = s2.get_patch( iu, iv );
                surface_type *p = get_patch( iu, iv );

                p->sum( *p1, *p2 );
              }
            }
          }

          // RST is an alternate parameterization that assumes a u, v surface encloses a volume.
          // R [0, 1] runs the length of the volume (U-direction)
          // S [0, 1] runs the width of the volume (V-direction)
          // T [0, 1] runs the thickness of the volume
          point_type fRST(const data_type &r, const data_type &s, const data_type &t) const
          {
            data_type umax, umin, du, vmax, vmin, dv;
            get_parameter_min( umin, vmin );
            get_parameter_max( umax, vmax );
            du = umax - umin;
            dv = vmax - vmin;

            data_type u = umin + r * du;
            data_type vlow = vmin + 0.5 * s * dv;
            data_type vup = vmax - 0.5 * s * dv;

            point_type xup, xlow;

            xup = f( u, vup );
            xlow = f( u, vlow );

            return ( 1.0 - t ) * xlow + t * xup;
          }

          point_type f_R(const data_type &r, const data_type &s, const data_type &t) const
          {
            data_type umax, umin, du, vmax, vmin, dv;
            get_parameter_min( umin, vmin );
            get_parameter_max( umax, vmax );
            du = umax - umin;
            dv = vmax - vmin;

            data_type u = umin + r * du;
            data_type vlow = vmin + 0.5 * s * dv;
            data_type vup = vmax - 0.5 * s * dv;

            point_type dxup_du = f_u( u, vup );
            point_type dxlow_du = f_u( u, vlow );

            point_type dx_dR = ( 1.0 - t ) * dxlow_du * du + t * dxup_du * du;
            return dx_dR;
          }

          point_type f_S(const data_type &r, const data_type &s, const data_type &t) const
          {
            data_type umax, umin, du, vmax, vmin, dv;
            get_parameter_min( umin, vmin );
            get_parameter_max( umax, vmax );
            du = umax - umin;
            dv = vmax - vmin;

            data_type u = umin + r * du;
            data_type vlow = vmin + 0.5 * s * dv;
            data_type vup = vmax - 0.5 * s * dv;

            point_type dxup_dv = f_v( u, vup );
            point_type dxlow_dv = f_v( u, vlow );

            point_type dx_dS = 0.5 * ( ( 1.0 - t ) * dxlow_dv * dv - t * dxup_dv * dv );
            return dx_dS;
          }

          point_type f_T(const data_type &r, const data_type &s, const data_type &t) const
          {
            (void) t; // Silence un-used parameter warning.
            data_type umax, umin, du, vmax, vmin, dv;
            get_parameter_min( umin, vmin );
            get_parameter_max( umax, vmax );
            du = umax - umin;
            dv = vmax - vmin;

            data_type u = umin + r * du;
            data_type vlow = vmin + 0.5 * s * dv;
            data_type vup = vmax - 0.5 * s * dv;

            point_type xup, xlow;

            xup = f( u, vup );
            xlow = f( u, vlow );

            point_type dx_dT = xup - xlow;
            return dx_dT;
          }

          // TODO: NEED TO IMPLEMENT
          //       * fit
          //       * interpolate

        private:
//           template<template<typename, unsigned short, typename> class surf1__,
//                    typename data1__, unsigned short dim1__, typename tol1__>
//           friend void area(typename piecewise<surf1__, data1__, dim1__, tol1__>::data_type &len,
//                            const piecewise<surf1__, data1__, dim1__, tol1__> &pc,
//                            const typename piecewise<surf1__, data1__, dim1__, tol1__>::data_type &tol);
//           template<template<typename, unsigned short, typename> class surf1__,
//                             typename data1__, unsigned short dim1__, typename tol1__>
//           friend void area(typename piecewise<surf1__, data1__, dim1__, tol1__>::data_type &len,
//                            const piecewise<surf1__, data1__, dim1__, tol1__> &pc,
//                            const typename piecewise<surf1__, data1__, dim1__, tol1__>::data_type &t0,
//                            const typename piecewise<surf1__, data1__, dim1__, tol1__>::data_type &t1,
//                            const typename piecewise<surf1__, data1__, dim1__, tol1__>::data_type &tol);

          template<template<typename, unsigned short, typename> class surface1__, typename data1__, unsigned short dim1__, typename tol1__ >
          friend typename piecewise<surface1__, data1__, dim1__, tol1__>::data_type
            eli::geom::intersect::find_rst(
              typename piecewise<surface1__, data1__, dim1__, tol1__>::data_type &r,
              typename piecewise<surface1__, data1__, dim1__, tol1__>::data_type &s,
              typename piecewise<surface1__, data1__, dim1__, tol1__>::data_type &t,
              const piecewise<surface1__, data1__, dim1__, tol1__> &ps,
              const typename piecewise<surface1__, data1__, dim1__, tol1__>::point_type &pt,
              typename piecewise<surface1__, data1__, dim1__, tol1__>::index_type &ret );

          template<template<typename, unsigned short, typename> class surface1__, typename data1__, unsigned short dim1__, typename tol1__ >
          friend typename piecewise<surface1__, data1__, dim1__, tol1__>::data_type
            eli::geom::intersect::minimum_distance(
              typename piecewise<surface1__, data1__, dim1__, tol1__>::data_type &u,
              typename piecewise<surface1__, data1__, dim1__, tol1__>::data_type &v,
              const piecewise<surface1__, data1__, dim1__, tol1__> &ps,
              const typename piecewise<surface1__, data1__, dim1__, tol1__>::point_type &pt);

          template<template<typename, unsigned short, typename> class surface1__, typename data1__, unsigned short dim1__, typename tol1__ >
          friend typename piecewise<surface1__, data1__, dim1__, tol1__>::data_type
            eli::geom::intersect::intersect(
              typename piecewise<surface1__, data1__, dim1__, tol1__>::data_type &u,
              typename piecewise<surface1__, data1__, dim1__, tol1__>::data_type &v,
              typename piecewise<surface1__, data1__, dim1__, tol1__>::point_type &p,
              const piecewise<surface1__, data1__, dim1__, tol1__> &ps,
              const typename piecewise<surface1__, data1__, dim1__, tol1__>::point_type &p0,
              const typename piecewise<surface1__, data1__, dim1__, tol1__>::index_type &iproj);

          template<template<typename, unsigned short, typename> class surface1__, typename data1__, unsigned short dim1__, typename tol1__ >
          friend void
            eli::geom::intersect::intersect_segment( std::vector<data1__> &tvec,
              const piecewise<surface1__, data1__, dim1__, tol1__> &ps,
              const typename piecewise<surface1__, data1__, dim1__, tol1__>::point_type &pt,
              const typename piecewise<surface1__, data1__, dim1__, tol1__>::point_type &vec,
              const typename piecewise<surface1__, data1__, dim1__, tol1__>::bounding_box_type &bbox);

          typedef std::map< data_type, index_type > keymap_type;

          struct parameter_key
          {
            keymap_type key;
            data_type pmax;

            parameter_key() : pmax(0) {}
            parameter_key(const parameter_key &pk) : key(pk.key), pmax(pk.pmax) {}
            ~parameter_key() {}

            bool operator==(const parameter_key &pk) const
            {
              if (this==&pk)
                return true;
              if (pmax!=pk.pmax)
                return false;
              if (key!=pk.key)
                return false;

              return true;
            }

            bool operator!=(const parameter_key &pk) const
            {
              return !operator==(pk);
            }

            void clear()
            {
              pmax=0;
              key.clear();
            }

            data_type get_pmax() const
            {
              return pmax;
            }

            data_type get_pmin() const
            {
              if(!key.empty())
                return key.begin()->first;
              else
                return pmax;
            }

            void set_pmax(const data_type &pmax_in)
            {
              pmax = pmax_in;
            }

            void set_pmin(const data_type &pmin_in)
            {
              if(!key.empty())
              {
                if(pmin_in != key.begin()->first)
                {
                  data_type p = pmin_in;
                  keymap_type shiftkey;
                  for (typename keymap_type::iterator it=key.begin(); it!=key.end(); ++it)
                  {
                    data_type delta_p = get_delta_parm(it);

                    shiftkey.insert(shiftkey.end(), std::make_pair(p, it->second));

                    p+=delta_p;
                  }
                  key.swap(shiftkey);
                  pmax = p;
                }
              }
              else
              {
                pmax=pmin_in;
              }
            }

            void init(const index_type &nseg, const data_type &dp = 1, const data_type &p0 = 0)
            {
              key.clear();
              pmax = p0;
              append(nseg, dp);
            }

            void append(const data_type &dp = 1)
            {
              typename keymap_type::iterator itguess = key.end();
              itguess = key.insert(itguess, std::make_pair(pmax, key.size()));
              pmax += dp;
            }

            void append(const index_type &nseg, const data_type &dp = 1)
            {
              typename keymap_type::iterator itguess = key.end();
              index_type j = key.size();
              data_type p = pmax;
              for(index_type i = 0; i < nseg; ++i)
              {
                itguess = key.insert(itguess, std::make_pair(p, j));
                p += dp;
                ++j;
              }
              pmax = p;
            }

            template<typename it__>
            void init(const it__ &dps, const it__ &dpe, const data_type &p0 = 0)
            {
              key.clear();
              pmax = p0;
              append(dps, dpe);
            }

            template<typename it__>
            void append(const it__ &dps, const it__ &dpe)
            {
              typename keymap_type::iterator itguess = key.end();
              index_type j = key.size();
              data_type p = pmax;
              for(it__ dp = dps; dp != dpe; ++dp)
              {
                itguess = key.insert(itguess, std::make_pair(p, j));
                p += (*dp);
                ++j;
              }
              pmax = p;
            }

            void init( const std::vector<data_type> &pmap )
            {
              key.clear();
              typename keymap_type::iterator itguess = key.end();

              for ( index_type j = 0; j < pmap.size() - 1; j++ )
              {
                itguess = key.insert( itguess, std::make_pair( pmap[j], j) );
              }

              pmax = pmap.back();
            }

            void parameter_report() const
            {
              printf("Parameter report:\n");
              typename keymap_type::const_iterator it;

              int i = 0;
              // cycle through all segments to get each bounding box to add
              for (it=key.begin(); it!=key.end(); ++it)
              {
                printf(" seg: %d \t p: %f \t pk %d\n", i, it->first, it->second);
                ++i;
              }
              printf(" pmax: %f\n", pmax);
              printf("End report\n");
            }

            void get_pmap(std::vector<data_type> &pmap) const
            {
              pmap.clear();

              typename keymap_type::const_iterator it;
              for (it=key.cbegin(); it!=key.cend(); ++it)
              {
                pmap.push_back( it->first );
              }
              pmap.push_back( pmax );
            }

            void reverse_keymap()
            {
              keymap_type rkey;
              typename keymap_type::iterator itr;
              typename keymap_type::iterator itrguess = rkey.begin();

              data_type p = get_pmin();

              for (typename keymap_type::reverse_iterator it=key.rbegin(); it!=key.rend(); ++it)
              {
                itr = rkey.insert(itrguess, std::make_pair(p, it->second));

                data_type delta_p = get_delta_parm(it);
                p += delta_p;

                itrguess = itr;
              }
              key.swap(rkey);

              // Parametric length should stay the same.
              tolerance_type tol;
              assert(tol.approximately_equal(p, pmax));
            }

            void roll_keymap( const typename keymap_type::const_iterator &itstart )
            {
              keymap_type rollkey;
              typename keymap_type::const_iterator itr;
              typename keymap_type::const_iterator itguess = rollkey.begin();

              data_type p = get_pmin();

              for ( typename keymap_type::const_iterator it = itstart; it != key.end(); ++it )
              {
                itr = rollkey.insert( itguess, std::make_pair( p, it->second ) );

                data_type delta_p = get_delta_parm( it );
                p += delta_p;

                itguess = itr;
              }
              for ( typename keymap_type::const_iterator it = key.begin(); it != itstart; ++it )
              {
                itr = rollkey.insert( itguess, std::make_pair( p, it->second ) );

                data_type delta_p = get_delta_parm( it );
                p += delta_p;

                itguess = itr;
              }

              key.swap( rollkey );
            }

            void roll_keymap( const index_type &index )
            {
              index_type ikey;
              typename keymap_type::const_iterator itstart;

              find_segment( ikey, itstart, index );

              roll_keymap( itstart );
            }

            data_type get_delta_parm(const typename keymap_type::iterator &it) const
            {
              assert (it != key.end());

              typename keymap_type::iterator itnext = it;
              itnext++;

              data_type delta_p;

              if(itnext != key.end())
                delta_p = itnext->first - it->first;
              else
                delta_p = pmax - it->first;

              return delta_p;
            }

            data_type get_delta_parm(const typename keymap_type::const_iterator &it) const
            {
              assert (it != key.end());

              typename keymap_type::const_iterator itnext = it;
              itnext++;

              data_type delta_p;

              if(itnext != key.end())
                delta_p = itnext->first - it->first;
              else
                delta_p = pmax - it->first;

              return delta_p;
            }

            data_type get_delta_parm(const typename keymap_type::reverse_iterator &it) const
            {
              assert (it != key.rend());

              data_type delta_p;

              if(it != key.rbegin())
              {
                typename keymap_type::reverse_iterator itprev = it;
                itprev--;
                delta_p = itprev->first - it->first;
              }
              else
              {
                delta_p = pmax - it->first;
              }

              return delta_p;
            }

            data_type get_delta_parm(const typename keymap_type::const_reverse_iterator &it) const
            {
              assert (it != key.rend());

              data_type delta_p;

              if(it != key.rbegin())
              {
                typename keymap_type::const_reverse_iterator itprev = it;
                itprev--;
                delta_p = itprev->first - it->first;
              }
              else
              {
                delta_p = pmax - it->first;
              }

              return delta_p;
            }

            void find_segment(index_type &ikey, typename keymap_type::const_iterator &it, const index_type &index) const
            {
              if(index >= (int) key.size() || index < 0)
              {
                it=key.end();
                ikey=-1;
                return;
              }

              // advance to desired index
              index_type i;
              for (i=0, it=key.begin(); i<index; ++i, ++it) {}

              ikey=it->second;
            }

            void find_segment(index_type &ikey, typename keymap_type::iterator &it, const index_type &index) const
            {
              if(index >= (int) key.size() || index < 0)
              {
                it=key.end();
                ikey=-1;
                return;
              }

              // advance to desired index
              index_type i;
              for (i=0, it=key.begin(); i<index; ++i, ++it) {}

              ikey=it->second;
            }

            index_type find_segment(index_type &ikey, typename keymap_type::iterator &it, data_type &pp, const data_type &p_in)
            {
              tol__ tol;

              if(p_in>pmax)
              {
                it=key.end();
                ikey = -1;
                return 0;
              }

              data_type pmin = get_pmin();

              if(p_in<pmin)
              {
                it=key.end();
                ikey = -1;
                return 0;
              }

              // Use map::upper_bound for fast lookup of segment after p_in
              it=key.upper_bound(p_in);

              // Decrement to segment containing p_in
              if(it != key.begin())
                it--;

              ikey = it->second;

              // At start of segment
              if(p_in==it->first)
              {
                pp=static_cast<data_type>(0);
                return -1;
              }

              data_type delta_p = get_delta_parm(it);

              // At end of segment
              if(p_in == (it->first + delta_p))
              {
                pp=static_cast<data_type>(1);
                return 1;
              }

              // Typical case
              pp=(p_in-it->first)/delta_p;

              // Super careful checks
              if (pp>static_cast<data_type>(1))
                pp=static_cast<data_type>(1);
              if (pp<static_cast<data_type>(0))
                pp=static_cast<data_type>(0);

              return 0;
            }

            index_type find_segment(index_type &ikey, typename keymap_type::const_iterator &it, data_type &pp, const data_type &p_in) const
            {
              tol__ tol;

              if(p_in>pmax)
              {
                it=key.end();
                ikey = -1;
                return 0;
              }

              data_type pmin = get_pmin();

              if(p_in<pmin)
              {
                it=key.end();
                ikey = -1;
                return 0;
              }

              // Use map::upper_bound for fast lookup of segment after p_in
              it=key.upper_bound(p_in);

              // Decrement to segment containing p_in
              if(it != key.begin())
                it--;

              ikey = it->second;

              // At start of segment
              if(p_in == it->first)
              {
                pp=static_cast<data_type>(0);
                return -1;
              }

              data_type delta_p = get_delta_parm(it);

              // At end of segment
              if(p_in == (it->first + delta_p))
              {
                pp=static_cast<data_type>(1);
                return 1;
              }

              // Typical case
              pp=(p_in-it->first)/delta_p;

              // Super careful checks
              if (pp>static_cast<data_type>(1))
                pp=static_cast<data_type>(1);
              if (pp<static_cast<data_type>(0))
                pp=static_cast<data_type>(0);

              return 0;
            }

            index_type find_index(const typename keymap_type::const_iterator &it) const
            {
              return distance( key.begin(), it );
            }
          };


          typedef std::vector< surface_type > patch_strip_type;
          typedef std::vector< patch_strip_type > patch_collection_type;

          patch_collection_type patches;
          // By convention, patches[uk][vk]

          parameter_key ukey, vkey;
          index_type nu, nv;

          enum close_cache
          {
            UNKNOWN = 0,
            CLOSED = 1,
            OPEN = 2
          };

          mutable index_type uclosecache, vclosecache;

        protected:
          bool check_continuity(const eli::geom::general::continuity &/*cont*/) const
          {
            // TODO: Need to implement this
            return true;
          }

          bool check_u_continuity(const surface_type &/*s1*/, const surface_type &/*s2*/, const eli::geom::general::continuity &/*cont*/) const
          {
            // TODO: Need to implement this
            return true;
          }

          bool check_v_continuity(const surface_type &/*s1*/, const surface_type &/*s2*/, const eli::geom::general::continuity &/*cont*/) const
          {
            // TODO: Need to implement this
            return true;
          }

        private:

          void resize_store(const index_type &nu_in, const index_type &nv_in)
          {
            if ((nu_in<=0) || (nv_in<=0))
              return;

            patches.resize(nu_in);
            nu = nu_in;

            // Unconditionally do this to make sure newly added rows are properly sized.
            for(index_type i = 0; i < nu_in; i++)
              patches[i].resize(nv_in);

            nv = nv_in;
          }

          error_code subsurf(piecewise<surface__, data_type, dim__, tol__> &surf, const typename keymap_type::const_iterator &ustart, const typename keymap_type::const_iterator &uend, const typename keymap_type::const_iterator &vstart, const typename keymap_type::const_iterator &vend ) const
          {
            surf.clear();

            surf.set_u0( ustart->first );
            surf.set_v0( vstart->first );

            typename keymap_type::const_iterator uit, vit;

            index_type nusub = 0;
            for ( uit = ustart; uit != uend; uit++ )
            {
              data_type du = ukey.get_delta_parm( uit );
              surf.ukey.append( du );
              nusub++;
            }

            index_type nvsub = 0;
            for ( vit = vstart; vit != vend; vit++ )
            {
              data_type dv = vkey.get_delta_parm( vit );
              surf.vkey.append( dv );
                nvsub++;
            }

            surf.resize_store( nusub, nvsub );

            index_type ikstore, jkstore;
            ikstore = 0;
            for ( uit = ustart; uit != uend; uit++ )
            {
              jkstore = 0;
              for ( vit = vstart; vit != vend; vit++ )
              {
                surf.patches[ikstore][jkstore] = patches[(*uit).second][(*vit).second];
                jkstore++;
              }
              ikstore++;
            }

            return NO_ERRORS;
          }

          error_code split_u(const index_type &uk, const typename keymap_type::iterator &uit, const data_type &u_in, const data_type &uu)
          {
            tolerance_type tol;
            assert(!tol.approximately_equal(uu, 0));
            assert(!tol.approximately_equal(uu, 1));

            index_type ukr, vk;
            // Right half will be added at end of patch matrix.
            ukr=nu;
            ukey.key.insert(uit, std::make_pair(u_in, ukr));

            // Increase matrix size.
            resize_store(nu+1, nv);

            for (vk=0; vk<nv; ++vk)
            {
              surface_type s = patches[uk][vk];
              s.split_u(patches[uk][vk], patches[ukr][vk], uu);
            }

            return NO_ERRORS;
          }

          error_code split_v(const index_type &vk, const typename keymap_type::iterator &vit, const data_type &v_in, const data_type &vv)
          {
            tolerance_type tol;
            assert(!tol.approximately_equal(vv, 0));
            assert(!tol.approximately_equal(vv, 1));

            index_type uk, vkr;
            // Right half will be added at end of patch matrix.
            vkr=nv;
            vkey.key.insert(vit, std::make_pair(v_in, vkr));

            // Increase matrix size.
            resize_store(nu, nv+1);

            for (uk=0; uk<nu; ++uk)
            {
              surface_type s = patches[uk][vk];
              s.split_v(patches[uk][vk], patches[uk][vkr], vv);
            }

            return NO_ERRORS;
          }

          // Lookup based on i,j
          void find_patch(index_type &uk, index_type &vk,
                          typename keymap_type::iterator &uit, typename keymap_type::iterator &vit,
                          const index_type & ui, const index_type &vi)
          {
            ukey.find_segment(uk, uit, ui);
            vkey.find_segment(vk, vit, vi);
          }

          void find_patch(typename keymap_type::iterator &uit, typename keymap_type::iterator &vit,
                          const index_type & ui, const index_type &vi)
          {
            index_type uk, vk;
            find_patch(uk, vk, uit, vit, ui, vi);
          }

          void find_patch(index_type &uk, index_type &vk,
                          typename keymap_type::const_iterator &uit, typename keymap_type::const_iterator &vit,
                          const index_type & ui, const index_type &vi) const
          {
            ukey.find_segment(uk, uit, ui);
            vkey.find_segment(vk, vit, vi);
          }

          void find_patch(typename keymap_type::const_iterator &uit, typename keymap_type::const_iterator &vit,
                          const index_type & ui, const index_type &vi) const
          {
            index_type uk, vk;
            find_patch(uk, vk, uit, vit, ui, vi);
          }

          void find_patch(index_type &uk, index_type &vk,
                          const index_type & ui, const index_type &vi) const
          {
            typename keymap_type::const_iterator uit, vit;
            find_patch(uk, vk, uit, vit, ui, vi);
          }

          // Lookup based on u_in, v_in.
          patch_boundary_code_type find_patch(index_type &uk, index_type &vk,
                          typename keymap_type::iterator &uit, typename keymap_type::iterator &vit,
                          data_type &uu, data_type &vv,
                          const data_type &u_in, const data_type &v_in)
          {
            index_type ucode, vcode;
            ucode = ukey.find_segment(uk, uit, uu, u_in);
            vcode = vkey.find_segment(vk, vit, vv, v_in);
            return std::make_pair( ucode, vcode );
          }

          patch_boundary_code_type find_patch(typename keymap_type::iterator &uit, typename keymap_type::iterator &vit,
                          data_type &uu, data_type &vv,
                          const data_type &u_in, const data_type &v_in)
          {
            index_type uk, vk;
            return find_patch(uk, vk, uit, vit, uu, vv, u_in, v_in);
          }

          patch_boundary_code_type find_patch(index_type &uk, index_type &vk,
                          typename keymap_type::const_iterator &uit, typename keymap_type::const_iterator &vit,
                          data_type &uu, data_type &vv,
                          const data_type &u_in, const data_type &v_in) const
          {
            index_type ucode, vcode;
            ucode = ukey.find_segment(uk, uit, uu, u_in);
            vcode = vkey.find_segment(vk, vit, vv, v_in);
            return std::make_pair( ucode, vcode );
          }

          patch_boundary_code_type find_patch(typename keymap_type::const_iterator &uit, typename keymap_type::const_iterator &vit,
                          data_type &uu, data_type &vv,
                          const data_type &u_in, const data_type &v_in) const
          {
            index_type uk, vk;
            return find_patch(uk, vk, uit, vit, uu, vv, u_in, v_in);
          }

          patch_boundary_code_type find_patch(index_type &uk, index_type &vk,
                          data_type &uu, data_type &vv,
                          const data_type &u_in, const data_type &v_in) const
          {
            typename keymap_type::const_iterator uit, vit;
            return find_patch(uk, vk, uit, vit, uu, vv, u_in, v_in);
          }

      };
    }
  }
}
#endif
