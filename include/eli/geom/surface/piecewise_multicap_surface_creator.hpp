/*********************************************************************************
* Copyright (c) 2014 David D. Marshall <ddmarsha@calpoly.edu>
*
* All rights reserved. This program and the accompanying materials
* are made available under the terms of the Eclipse Public License v1.0
* which accompanies this distribution, and is available at
* http://www.eclipse.org/legal/epl-v10.html
*
* Contributors:
*    Rob McDonald -- based on piecewise_capped_surface_creator.hpp
********************************************************************************/

#ifndef eli_geom_surface_piecewise_multicap_surface_creator_hpp
#define eli_geom_surface_piecewise_multicap_surface_creator_hpp

#include <list>
#include <iterator>

#include "eli/code_eli.hpp"

#include "eli/constants/math.hpp"

#include "eli/util/tolerance.hpp"

#include "eli/geom/curve/bezier.hpp"
#include "eli/geom/curve/piecewise.hpp"
#include "eli/geom/curve/piecewise_creator.hpp"

#include "eli/geom/surface/bezier.hpp"
#include "eli/geom/surface/piecewise.hpp"
#include "eli/geom/surface/piecewise_creator_base.hpp"
#include "eli/geom/surface/piecewise_connection_data.hpp"
#include "eli/mutil/dm/binomial_coefficient.hpp"
#include "eli/mutil/opt/least_squares.hpp"

namespace eli
{
  namespace geom
  {
    namespace surface
    {
      template<typename data__, unsigned short dim__, typename tol__>
      class piecewise_multicap_surface_creator : public piecewise_creator_base<data__, dim__, tol__>
      {
        public:
          enum edge_cap_identifier
          {
            CAP_NONE,
            CAP_UMIN,
            CAP_UMAX,
            CAP_VMIN,
            CAP_VMAX
          };

          enum cap_type
          {
            FLAT,
            ROUND,
            EDGE,
            SHARP,
            POINT,
            ROUND_EXT,
            CONE,
            OGIVE,
            PROFILE
          };

          // The profiles a PROFILE cap's sides follow from the end to its point.
          enum profile_shape
          {
            POWER_LAW,
            VON_KARMAN,
            SEARS_HAACK,
          };

          typedef piecewise_creator_base<data__, dim__, tol__> base_class_type;
          typedef typename base_class_type::data_type data_type;
          typedef typename base_class_type::point_type point_type;
          typedef typename base_class_type::index_type index_type;
          typedef typename base_class_type::tolerance_type tolerance_type;
          typedef typename base_class_type::piecewise_surface_type piecewise_surface_type;
          typedef typename piecewise_surface_type::piecewise_curve_type piecewise_curve_type;
          typedef typename piecewise_curve_type::onedpiecewisecurve oned_piecewise_curve_type;

          typedef eli::geom::curve::piecewise_cubic_spline_creator<double, 1, tolerance_type> piecewise_oned_cubic_spline_creator_type;


          piecewise_multicap_surface_creator()
            : piecewise_creator_base<data__, dim__, tol__>(0, 0), curves_given(false), prep_done(false), surface_prep_done(false), delta_param(1), edge_to_cap(CAP_NONE)
          {
          }
          piecewise_multicap_surface_creator(const piecewise_surface_type &os, const data_type &dp, edge_cap_identifier ec)
            : piecewise_creator_base<data__, dim__, tol__>(0, 0), orig_surface(os), curves_given(false), prep_done(false), surface_prep_done(false), delta_param(dp), edge_to_cap(ec)
          {
          }
          piecewise_multicap_surface_creator(const piecewise_multicap_surface_creator<data_type, dim__, tolerance_type> & gs)
            : piecewise_creator_base<data_type, dim__, tolerance_type>(gs), orig_surface(gs.orig_surface),
              orig_edge(gs.orig_edge), orig_edge_next(gs.orig_edge_next), curves_given(gs.curves_given), prep_done(false), surface_prep_done(false), delta_param(gs.delta_param), edge_to_cap(gs.edge_to_cap)
          {
          }
          virtual ~piecewise_multicap_surface_creator()
          {
          }

          bool set_conditions(const piecewise_surface_type &os, const index_type & typ, const data_type &dp, edge_cap_identifier ec, const data_type &lf, const data_type &of, const data_type &sf, const point_type & pnt, const bool &sc )
          {
            // all deltas are positive, even on the min edges
            if (dp<=0)
              return false;

            cap_type = typ;
            orig_surface = os;
            curves_given = false;
            edge_to_cap = ec;
            delta_param = dp;
            len_factor = lf;
            off_factor = of;
            pnt_offset = pnt;
            str_factor = sf;
            sweep_correct = sc;
            rad_factor = 0;
            ogive_angle = 0;
            ogive_surfref = false;
            profile = POWER_LAW;
            profile_exp = 1;
            profile_haack = 0;
            truncate = false;
            frame_given = false;
            shear_offset = false;

            if ( orig_surface.number_u_patches() == 0 || orig_surface.number_v_patches() == 0 )
            {
              return false;
            }

            return true;
          }

          // As above, given the row of control points along the edge to cap and the row next to it, as
          // get_edge_curves finds them, in place of the surface.  The tip curve and the curves built on it
          // can be created; create has no surface to join a cap to, and fails.
          bool set_conditions(const piecewise_curve_type &e, const piecewise_curve_type &n, const index_type & typ, const data_type &dp, const data_type &lf, const data_type &of, const data_type &sf, const point_type & pnt, const bool &sc )
          {
            // all deltas are positive, even on the min edges
            if (dp<=0)
              return false;

            cap_type = typ;
            orig_surface.clear();
            orig_edge = e;
            orig_edge_next = n;
            curves_given = true;
            edge_to_cap = CAP_NONE;
            delta_param = dp;
            len_factor = lf;
            off_factor = of;
            pnt_offset = pnt;
            str_factor = sf;
            sweep_correct = sc;
            rad_factor = 0;
            ogive_angle = 0;
            ogive_surfref = false;
            profile = POWER_LAW;
            profile_exp = 1;
            profile_haack = 0;
            truncate = false;
            frame_given = false;
            shear_offset = false;

            if ( orig_edge.number_segments() == 0 || orig_edge.number_segments() != orig_edge_next.number_segments() )
            {
              return false;
            }

            return true;
          }

          void set_h_vals( const data_type & vm, const data_type & hm, const data_type &hs, const data_type &he )
          {
            v_hmax = vm;
            hmax = hm;
            h_start = hs;
            h_end = he;
          }

          void set_ext_conditions( const bool & le, const bool & te )
          {
            extendle = le;
            extendte = te;
          }

          // Radius of a cone cap's spherical tip, in units of the cap's length scale.
          void set_cone_conditions( const data_type & rf )
          {
            rad_factor = rf;
          }

          // Angle, in radians, at which an ogive cap's sides leave the end, turned about the end's tangent,
          // inwards for positive angles, from the surface's direction across the end where surfref is set, or
          // else from the end's normal.  The sphere's radius is set as a cone cap's.
          void set_ogive_conditions( const data_type & ang, const bool & surfref )
          {
            ogive_angle = ang;
            ogive_surfref = surfref;
          }

          // The profile a PROFILE cap's sides follow, for a power law its exponent, and for Von Karman the
          // Haack series' C.
          void set_profile_conditions( const index_type & shape, const data_type & exponent, const data_type & haack )
          {
            profile = shape;
            profile_exp = exponent;
            profile_haack = haack;
          }

          // Whether a cone, ogive or profile cap ends at a flat face in place of its point: the end scaled
          // about its middle by the radius factor and moved to the point.
          void set_truncate_conditions( const bool & tr )
          {
            truncate = tr;
          }

          // The frame the surface is built in.  A cap that closes on a point runs its length along axis, away
          // from the body, and its offset along up, made square to axis, each in units of the end's size
          // across axis.  axis points away from the body where the surface's direction across the end does
          // not say which way that is.  Without a frame the cap finds its direction from the end.
          void set_frame_conditions( const point_type &axis, const point_type &up )
          {
            tolerance_type tol;

            data_type alen = axis.norm();
            if ( tol.approximately_equal( alen, 0.0 ) )
            {
              frame_given = false;
              return;
            }
            frame_axis = axis / alen;

            frame_up = up - frame_axis * frame_axis.dot( up );
            data_type ulen = frame_up.norm();
            if ( tol.approximately_equal( ulen, 0.0 ) )
            {
              frame_given = false;
              return;
            }
            frame_up = frame_up / ulen;
            frame_given = true;
          }

          // Whether a cap that closes on a point in a frame takes its offset by shearing: built with none,
          // then each point moved along the frame's up direction by the offset in proportion to its distance
          // from the plane of the end, so the end stays put and the point reaches the offset.  This is how
          // DATCOM cambers a body.  A point or profile cap's offset is that shear already.
          void set_shear_conditions( const bool & sh )
          {
            shear_offset = sh;
          }

          // Whether this cap ends at a flat face: a truncated cone, ogive or profile cap with a face.
          bool truncated() const
          {
            return truncate && rad_factor > 0 && ( cap_type == CONE || cap_type == OGIVE || cap_type == PROFILE );
          }

          // Whether a cap of type typ closes on a single point.
          static bool closes_on_point( const index_type &typ )
          {
            switch ( typ )
            {
              case POINT:
              case CONE:
              case OGIVE:
              case PROFILE:
                return true;
              default:
                return false;
            }
          }

          // Whether this cap's sides run from the end onto a spherical tip that touches its point, as a cone
          // cap's do.
          // A truncated cap's sides run to its face instead.
          bool sphere_tip() const
          {
            switch ( cap_type )
            {
              case CONE:
              case OGIVE:
                return true;
              default:
                return false;
            }
          }

          // Whether this cap's rows are the end scaled along a profile, as a profile cap's are.
          bool profile_rows() const
          {
            return cap_type == PROFILE;
          }

          virtual bool create(piecewise_surface_type &ps) const
          {
            typedef typename eli::geom::curve::piecewise<eli::geom::curve::bezier, data_type, dim__, tolerance_type> piecewise_curve_type;
            typedef typename piecewise_curve_type::curve_type curve_type;
            typedef typename piecewise_surface_type::surface_type surface_type;

            tolerance_type tol;

            if ( curves_given )
            {
              return false;
            }

            if ( !prep_done )
            {
              bool rval = prep();
              if ( !rval )
              {
                return rval;
              }
            }

            if ( !surface_prep_done )
            {
              bool rval = prep_surface();
              if ( !rval )
              {
                return rval;
              }
            }

            data_type tsplit( ( edge.get_t0() + edge.get_tmax() ) * 0.5 );

            // split the curve at the mid-parameter location and reverse the second piece's parameterization
            piecewise_curve_type first_half, second_half;
            edge.split(first_half, second_half, tsplit);
            second_half.reverse();
            second_half.set_t0(first_half.get_t0());

            piecewise_curve_type seam;
            seam.sum( first_half, second_half );
            seam.scale( 0.5 );

            data_type tmidseam = ( seam.get_t0() + seam.get_tmax() ) * 0.5;
            point_type porigin = seam.f( tmidseam );

            // A cap that closes on a point closes on the middle of the edge cap's tip curve.
            point_type ppoint;
            if ( !cap_point( cap_type, porigin, tmidseam, len_factor, off_factor, ppoint ) )
            {
              return false;
            }

            // A sheared cap is built with no offset, and sheared once built.
            bool shear = shear_offset && frame_given && closes_on_point( cap_type );
            if ( shear && !cap_point( cap_type, porigin, tmidseam, len_factor, 0, ppoint ) )
            {
              return false;
            }

            // A truncated cap ends at a flat face centred on the point: the end scaled about its middle by
            // tface and moved there, so the face's radius is rad_factor times the length scale, as a cone
            // cap's sphere's is, where the end's is its root mean square distance from its middle.
            bool trunc = truncated();
            point_type pface = ppoint + pnt_offset;
            data_type tface( 0 );
            if ( trunc )
            {
              data_type rend = end_radius( porigin );
              if ( rend > 0 )
              {
                tface = rad_factor * length_scale( porigin, tmidseam ) / rend;
              }
            }

            // A spherical tip touches the point, and the sides gather on the sphere at its pole.
            point_type pcen, axis, pdir;
            data_type rad( 0 );
            if ( sphere_tip() )
            {
              if ( !cone_sphere( porigin, tmidseam, ppoint + pnt_offset, pcen, rad, axis, pdir ) )
              {
                return false;
              }
            }

            // Split the ndelta pseudocurve.  Since this isn't really a curve, it only works because the split
            // point already exists and this becomes reordering and manipulating of existing values, and should
            // not change any values of control points.
            piecewise_curve_type nd_first_half, nd_second_half;
            ndelta.split(nd_first_half, nd_second_half, tsplit);
            nd_second_half.reverse();
            nd_second_half.set_t0(nd_first_half.get_t0());

            // An ogive cap's sides leave the end at an angle from its normal, or from the surface's
            // direction across it.
            point_type pnorm, dmean;
            pnorm << 0, 0, 0;
            dmean = pnorm;
            piecewise_curve_type ed_first_half, ed_second_half;
            std::vector < point_type > otop, obot;
            if ( cap_type == OGIVE && sphere_tip() )
            {
              pnorm = end_normal( axis );
              edelta.split( ed_first_half, ed_second_half, tsplit );
              ed_second_half.reverse();
              ed_second_half.set_t0( ed_first_half.get_t0() );

              // The pole is ahead of the center along the mean line: the circular arc that leaves the middle
              // of the end along its normal, or the mean of the surface's directions across the end at the
              // middles of its halves, and runs through the center.
              dmean = pnorm;
              if ( ogive_surfref )
              {
                point_type d1 = ed_first_half.f( tmidseam );
                point_type d2 = ed_second_half.f( tmidseam );
                if ( !tol.approximately_equal( d1.norm(), 0.0 ) && !tol.approximately_equal( d2.norm(), 0.0 ) )
                {
                  point_type dm = d1 / d1.norm() + d2 / d2.norm();
                  if ( !tol.approximately_equal( dm.norm(), 0.0 ) )
                  {
                    dmean = dm / dm.norm();
                  }
                }
              }
              ogive_pole_and_sides( first_half, second_half, ed_first_half, ed_second_half, tmidseam, pnorm, porigin,
                                    dmean, pcen, rad, pdir, otop, obot );
            }

            // A profile cap's profile from its point to the end, in ( s, f ): the distance from the point as a
            // fraction of the cap's length, and the radius as a fraction of the end's.
            // A truncated cap's profile starts where its radius is tface, s0 of the way from a point beyond
            // the face, pproftip, so the face is at the cap's point.
            piecewise_curve_type pprof;
            std::vector<data_type> tprof;
            point_type pproftip = ppoint + pnt_offset;
            if ( profile_rows() )
            {
              profile_function fun;
              profile_setup( fun );

              data_type tstart( 0 );
              if ( trunc )
              {
                tstart = profile_start( tface, fun );
              }

              profile_fit( pprof, tstart, fun );
              pprof.get_pmap( tprof );

              data_type s0 = pprof.f( pprof.get_t0() )( 0 );
              if ( trunc && s0 < 1 )
              {
                pproftip = ( pface - porigin * s0 ) / ( 1 - s0 );
              }
            }

            std::vector<data_type> pmap, dvcap, ducap;
            // Get edge parameter map
            first_half.get_pmap( pmap );
            dvcap.resize( pmap.size() - 1 );

            for ( index_type i = 0; i < pmap.size() - 1; i++ )
            {
              dvcap[i] = pmap[ i + 1 ] - pmap[ i ];
            }

            index_type i_hmax = 0;
            if ( cap_type == ROUND_EXT )
            {
              ducap.resize( 4, delta_param / 2.0 );
              i_hmax = seam.find( v_hmax );
            }
            else if ( sphere_tip() )
            {
              ducap.resize( 4, delta_param / 2.0 );
            }
            else if ( profile_rows() )
            {
              // Each half of the cap has the profile's cubics, from the end to the point and back again.
              // A truncated cap's halves give the profile half their parameter and the face the other half.
              index_type nprof = pprof.number_segments();
              data_type dprof = delta_param / ( tprof[nprof] - tprof[0] );
              if ( trunc )
              {
                dprof = 0.5 * dprof;
              }
              ducap.resize( 2 * nprof );
              for ( index_type ip = 0; ip < nprof; ip++ )
              {
                ducap[ip] = dprof * ( tprof[nprof - ip] - tprof[nprof - ip - 1] );
                ducap[nprof + ip] = dprof * ( tprof[ip + 1] - tprof[ip] );
              }
              if ( trunc )
              {
                ducap.insert( ducap.begin() + nprof, 2, 0.5 * delta_param );
              }
            }
            else
            {
              ducap.resize( 2, delta_param );
            }

            // create surface connecting two edges with param spacing of 2*delta_param and points for other two edges
            // cap u-direction goes from second_half curve to first_half curve
            // cap v-direction follows first_half direction
            piecewise_surface_type cap;
            {
              cap.init_uv( ducap.begin(), ducap.end(), dvcap.begin(), dvcap.end(), 0.0, pmap[0] );

              index_type nseg = first_half.number_segments();

              // Vector along upper curve in average 'forward' direction
              point_type fwd;
              fwd = ( first_half.f( tsplit ) - first_half.f( first_half.get_t0() ) );
              data_type fwdlen = fwd.norm();
              if ( tol.approximately_equal( fwdlen, 0.0 ) )
              {
                fwd << 0, 0, 0;
              }
              else
              {
                fwd = fwd / fwdlen;
              }

              for ( index_type i = 0; i < nseg; i++ )
              {
                curve_type cfirst, clast;
                curve_type ndfirst, ndlast;

                first_half.get( cfirst, i );
                second_half.get( clast, i );
                nd_first_half.get( ndfirst, i );
                nd_second_half.get( ndlast, i );

                index_type deg = cfirst.degree();

                surface_type s1, s2;
                surface_type s3, s4;

                switch( cap_type ){
                case FLAT:
                  {
                    s1.resize( 1, deg );
                    s2.resize( 1, deg );

                    for ( index_type j = 0; j <= deg; j++ )
                    {

                      point_type pup, psplit, pdn;

                      pup = clast.get_control_point( j );
                      pdn = cfirst.get_control_point( j );

                      psplit = ( pup + pdn ) * 0.5;

                      // Build linear upper side curve
                      s1.set_control_point( pup, 0, j );
                      s1.set_control_point( psplit, 1, j );

                      // Build linear lower side curve
                      s2.set_control_point( psplit, 0, j );
                      s2.set_control_point( pdn, 1, j );

                    }
                  }
                  break;
                case POINT:
                  {
                    s1.resize( 1, deg );
                    s2.resize( 1, deg );

                    for ( index_type j = 0; j <= deg; j++ )
                    {
                      point_type pup, pdn;

                      pup = clast.get_control_point( j );
                      pdn = cfirst.get_control_point( j );

                      // Build linear upper side curve
                      s1.set_control_point( pup, 0, j );
                      s1.set_control_point( ppoint + pnt_offset, 1, j );

                      // Build linear lower side curve
                      s2.set_control_point( ppoint + pnt_offset, 0, j );
                      s2.set_control_point( pdn, 1, j );
                    }
                  }
                  break;
                case CONE:
                  if ( trunc )
                  {
                    // Straight sides from the end to the face, then the face
                    s3.resize( 1, deg );
                    s1.resize( 1, deg );
                    s2.resize( 1, deg );
                    s4.resize( 1, deg );

                    for ( index_type j = 0; j <= deg; j++ )
                    {
                      point_type pup = clast.get_control_point( j );
                      point_type pdn = cfirst.get_control_point( j );
                      point_type fup = pface + ( pup - porigin ) * tface;
                      point_type fdn = pface + ( pdn - porigin ) * tface;

                      s3.set_control_point( pup, 0, j );
                      s3.set_control_point( fup, 1, j );
                      s1.set_control_point( fup, 0, j );
                      s1.set_control_point( pface, 1, j );
                      s2.set_control_point( pface, 0, j );
                      s2.set_control_point( fdn, 1, j );
                      s4.set_control_point( fdn, 0, j );
                      s4.set_control_point( pdn, 1, j );
                    }
                  }
                  else
                  {
                    s3.resize( 1, deg );
                    s1.resize( 3, deg );
                    s2.resize( 3, deg );
                    s4.resize( 1, deg );

                    point_type ptip = pcen + pdir * rad;

                    std::vector < point_type > tupv, aup1v, aup2v, tdnv, adn1v, adn2v;
                    cone_sides( clast, porigin, pcen, rad, pdir, tupv, aup1v, aup2v );
                    cone_sides( cfirst, porigin, pcen, rad, pdir, tdnv, adn1v, adn2v );

                    for ( index_type j = 0; j <= deg; j++ )
                    {
                      point_type pup, pdn;
                      point_type tup, aup1, aup2;
                      point_type tdn, adn1, adn2;

                      pup = clast.get_control_point( j );
                      pdn = cfirst.get_control_point( j );

                      tup = tupv[j];
                      aup1 = aup1v[j];
                      aup2 = aup2v[j];
                      tdn = tdnv[j];
                      adn1 = adn1v[j];
                      adn2 = adn2v[j];

                      // Upper cone side, then the sphere from its tangent point to the tip
                      s3.set_control_point( pup, 0, j );
                      s3.set_control_point( tup, 1, j );

                      s1.set_control_point( tup, 0, j );
                      s1.set_control_point( aup1, 1, j );
                      s1.set_control_point( aup2, 2, j );
                      s1.set_control_point( ptip, 3, j );

                      // The sphere from the tip to the lower tangent point, then the lower cone side
                      s2.set_control_point( ptip, 0, j );
                      s2.set_control_point( adn2, 1, j );
                      s2.set_control_point( adn1, 2, j );
                      s2.set_control_point( tdn, 3, j );

                      s4.set_control_point( tdn, 0, j );
                      s4.set_control_point( pdn, 1, j );
                    }
                  }
                  break;
                case OGIVE:
                  if ( trunc )
                  {
                    // Circular arcs from the end to the face, then the face
                    s3.resize( 3, deg );
                    s1.resize( 1, deg );
                    s2.resize( 1, deg );
                    s4.resize( 3, deg );

                    curve_type cdfirst, cdlast;
                    ed_first_half.get( cdfirst, i );
                    ed_second_half.get( cdlast, i );

                    point_type tup0, tup1, tdn0, tdn1;
                    end_tangents( second_half, first_half, i, tup0, tup1 );
                    end_tangents( first_half, second_half, i, tdn0, tdn1 );

                    std::vector < std::vector < point_type > > up, dn;
                    ogive_truncated_sides( clast, cdlast, tup0, tup1, pnorm, porigin, pface, tface, up );
                    ogive_truncated_sides( cfirst, cdfirst, tdn0, tdn1, pnorm, porigin, pface, tface, dn );

                    for ( index_type j = 0; j <= deg; j++ )
                    {
                      point_type pup = clast.get_control_point( j );
                      point_type pdn = cfirst.get_control_point( j );
                      point_type fup = pface + ( pup - porigin ) * tface;
                      point_type fdn = pface + ( pdn - porigin ) * tface;

                      s3.set_control_point( pup, 0, j );
                      s3.set_control_point( up[0][j], 1, j );
                      s3.set_control_point( up[1][j], 2, j );
                      s3.set_control_point( fup, 3, j );

                      s1.set_control_point( fup, 0, j );
                      s1.set_control_point( pface, 1, j );
                      s2.set_control_point( pface, 0, j );
                      s2.set_control_point( fdn, 1, j );

                      s4.set_control_point( fdn, 0, j );
                      s4.set_control_point( dn[1][j], 1, j );
                      s4.set_control_point( dn[0][j], 2, j );
                      s4.set_control_point( pdn, 3, j );
                    }
                  }
                  else
                  {
                    s3.resize( 3, deg );
                    s1.resize( 3, deg );
                    s2.resize( 3, deg );
                    s4.resize( 3, deg );

                    point_type ptip = pcen + pdir * rad;

                    curve_type cdfirst, cdlast;
                    ed_first_half.get( cdfirst, i );
                    ed_second_half.get( cdlast, i );

                    point_type tup0, tup1, tdn0, tdn1;
                    end_tangents( second_half, first_half, i, tup0, tup1 );
                    end_tangents( first_half, second_half, i, tdn0, tdn1 );

                    std::vector < std::vector < point_type > > up, dn;
                    ogive_sides( clast, cdlast, tup0, tup1, pnorm, porigin, pcen, rad, pdir, otop, obot, up );
                    ogive_sides( cfirst, cdfirst, tdn0, tdn1, pnorm, porigin, pcen, rad, pdir, otop, obot, dn );

                    for ( index_type j = 0; j <= deg; j++ )
                    {
                      point_type pup, pdn;

                      pup = clast.get_control_point( j );
                      pdn = cfirst.get_control_point( j );

                      // Upper ogive arc, then the sphere from its tangent point to the tip
                      s3.set_control_point( pup, 0, j );
                      s3.set_control_point( up[0][j], 1, j );
                      s3.set_control_point( up[1][j], 2, j );
                      s3.set_control_point( up[2][j], 3, j );

                      s1.set_control_point( up[2][j], 0, j );
                      s1.set_control_point( up[3][j], 1, j );
                      s1.set_control_point( up[4][j], 2, j );
                      s1.set_control_point( ptip, 3, j );

                      // The sphere from the tip to the lower tangent point, then the lower ogive arc
                      s2.set_control_point( ptip, 0, j );
                      s2.set_control_point( dn[4][j], 1, j );
                      s2.set_control_point( dn[3][j], 2, j );
                      s2.set_control_point( dn[2][j], 3, j );

                      s4.set_control_point( dn[2][j], 0, j );
                      s4.set_control_point( dn[1][j], 1, j );
                      s4.set_control_point( dn[0][j], 2, j );
                      s4.set_control_point( pdn, 3, j );
                    }
                  }
                  break;
                case PROFILE:
                  {
                    index_type nprof = pprof.number_segments();

                    for ( index_type ip = 0; ip < nprof; ip++ )
                    {
                      // The upper side's cubic from the end towards the tip, and the lower side's from the tip
                      // towards the end
                      curve_type cpup, cpdn;
                      pprof.get( cpup, nprof - 1 - ip );
                      pprof.get( cpdn, ip );

                      surface_type sup, sdn;
                      sup.resize( 3, deg );
                      sdn.resize( 3, deg );

                      for ( index_type j = 0; j <= deg; j++ )
                      {
                        point_type pup, pdn;

                        pup = clast.get_control_point( j );
                        pdn = cfirst.get_control_point( j );

                        for ( index_type k = 0; k <= 3; k++ )
                        {
                          sup.set_control_point( profile_row( cpup.get_control_point( 3 - k ), pup, porigin,
                                                              pproftip ), k, j );
                          sdn.set_control_point( profile_row( cpdn.get_control_point( k ), pdn, porigin,
                                                              pproftip ), k, j );
                        }
                      }

                      cap.set( sup, ip, i );
                      cap.set( sdn, cap.number_u_patches() - nprof + ip, i );
                    }

                    // A truncated cap's face, from the profile's first rows to the point
                    if ( trunc )
                    {
                      curve_type cpfirst;
                      pprof.get( cpfirst, 0 );
                      point_type cface = cpfirst.get_control_point( 0 );

                      surface_type sup, sdn;
                      sup.resize( 1, deg );
                      sdn.resize( 1, deg );
                      for ( index_type j = 0; j <= deg; j++ )
                      {
                        point_type fup = profile_row( cface, clast.get_control_point( j ), porigin, pproftip );
                        point_type fdn = profile_row( cface, cfirst.get_control_point( j ), porigin, pproftip );

                        sup.set_control_point( fup, 0, j );
                        sup.set_control_point( pface, 1, j );
                        sdn.set_control_point( pface, 0, j );
                        sdn.set_control_point( fdn, 1, j );
                      }
                      cap.set( sup, nprof, i );
                      cap.set( sdn, nprof + 1, i );
                    }
                  }
                  break;
                case ROUND:
                  {
                    s1.resize( 3, deg );
                    s2.resize( 3, deg );

                    for ( index_type j = 0; j <= deg; j++ )
                    {
                      point_type pup, psplit, pdn;
                      point_type tanup, tandn;

                      pup = clast.get_control_point( j );
                      pdn = cfirst.get_control_point( j );
                      tanup = ndlast.get_control_point( j );
                      tandn = ndfirst.get_control_point( j );

                      data_type ksweep = 1.0;
                      if ( sweep_correct )
                      {
                        point_type avetan = ( tanup + tandn ) * 0.5;
                        avetan = avetan / avetan.norm();
                        ksweep = 1.0 / sin( acos( avetan.dot( fwd ) ) );

                        if ( tol.approximately_equal( ksweep, 0.0) )
                        {
                          ksweep = 1.0;
                        }
                      }

                      // Displacement vector
                      point_type disp = pup - pdn;

                      // Chord of circle to construct
                      data_type circ_chord = disp.norm();

                      if ( tol.approximately_equal( circ_chord, 0.0 ) )
                      {
                        disp << 0, 0, 0;
                      }
                      else
                      {
                        disp = disp/circ_chord;
                      }

                      // Find angles between tangents and displacement
                      data_type dottanup = disp.dot( tanup );
                      data_type thetaup = acos( dottanup );

                      data_type dottandn = disp.dot( tandn );
                      data_type thetadn = acos( dottandn );

                      // Find circular arc to include
                      data_type theta = eli::constants::math<data__>::pi() + thetadn - thetaup;

                      if ( tol.approximately_equal( theta, 0) ) // Make flat panel avoid division by zero.
                      {
                        psplit = ( pup + pdn ) * 0.5;

                        // Build cubic upper side curve
                        s1.set_control_point( pup, 0, j );
                        s1.set_control_point( ( 2 * pup + psplit ) / 3, 1, j );
                        s1.set_control_point( ( pup + 2 * psplit ) / 3, 2, j );
                        s1.set_control_point( psplit, 3, j );

                        // Build cubic lower side curve
                        s2.set_control_point( psplit, 0, j );
                        s2.set_control_point( ( pdn + 2 * psplit ) / 3, 1, j );
                        s2.set_control_point( ( 2 * pdn + psplit ) / 3, 2, j );
                        s2.set_control_point( pdn, 3, j );
                      }
                      else
                      {
                        // Find radius from chord length.
                        data_type radius = circ_chord / ( 2.0 * sin( theta / 2.0 ) );

                        // Distance to quad Bezier construction point.
                        data_type b = radius * tan( theta / 4.0 );

                        // Extrapolated tip up/down points.
                        point_type ptup, ptdn;
                        ptup = pup + tanup * b * len_factor * ksweep;
                        ptdn = pdn + tandn * b * len_factor * ksweep;

                        // Displacement vector at extrapolated tip
                        point_type dispt = ptup - ptdn;

                        // Point on center of circle.
                        psplit = ( ptup + ptdn ) * 0.5 + dispt * off_factor;

                        // Fraction of radius to place cubic control point
                        data_type f = eli::constants::math<data_type>::cubic_bezier_circle_const() * tan(theta / 8.0);

                        // Build cubic upper side curve
                        s1.set_control_point( pup, 0, j );
                        s1.set_control_point( pup +   tanup * radius * f * len_factor * ksweep, 1, j );
                        s1.set_control_point( psplit + disp * radius * f, 2, j );
                        s1.set_control_point( psplit, 3, j );

                        // Build cubic lower side curve
                        s2.set_control_point( psplit, 0, j );
                        s2.set_control_point( psplit - disp * radius * f, 1, j );
                        s2.set_control_point( pdn +   tandn * radius * f * len_factor * ksweep, 2, j );
                        s2.set_control_point( pdn, 3, j );
                      }
                    }

                  break;
                case ROUND_EXT:
                  {
                    s1.resize( 3, deg );
                    s2.resize( 3, deg );
                    s3.resize( 1, deg );
                    s4.resize( 1, deg );

                    for ( index_type j = 0; j <= deg; j++ )
                    {
                      point_type pup, psplit_base, pdn;
                      point_type tanup, tandn;

                      pup = clast.get_control_point( j );
                      pdn = cfirst.get_control_point( j );

                      psplit_base = ( pup + pdn ) * 0.5;

                      tanup = ndlast.get_control_point( j );
                      tandn = ndfirst.get_control_point( j );

                      point_type avetan = ( tanup + tandn ) * 0.5;
                      avetan = avetan / avetan.norm();


                      data_type sweepfactor = sin( acos( avetan.dot( fwd ) ) );

                      if ( tol.approximately_equal( sweepfactor, 0.0) )
                      {
                        sweepfactor = 1.0;
                      }

                      data_type ksweep = 1.0;
                      if ( sweep_correct )
                      {
                        ksweep = sweepfactor;
                      }

                      // Displacement vector
                      point_type disp = pup - pdn;

                      // Chord of circle to construct
                      data_type base_chord = disp.norm();

                      if ( tol.approximately_equal( base_chord, 0.0 ) )
                      {
                        disp << 0, 0, 0;
                      }
                      else
                      {
                        disp = disp/base_chord;
                      }

                      // Find angles between tangents and displacement
                      data_type dottanup = disp.dot( tanup );
                      data_type thetaup = acos( dottanup );

                      data_type dottandn = disp.dot( tandn );
                      data_type thetadn = acos( dottandn );

                      // Find circular arc to include
                      data_type theta = eli::constants::math<data__>::pi() + thetadn - thetaup;

                      if ( tol.approximately_equal( theta, 0) ) // Make flat panel avoid division by zero.
                      {
                        // Build cubic upper side curve
                        s1.set_control_point( pup, 0, j );
                        s1.set_control_point( ( 2 * pup + psplit_base ) / 3, 1, j );
                        s1.set_control_point( ( pup + 2 * psplit_base ) / 3, 2, j );
                        s1.set_control_point( psplit_base, 3, j );

                        s3.set_control_point( pup, 0, j );
                        s3.set_control_point( pup, 1, j );

                        // Build cubic lower side curve
                        s2.set_control_point( psplit_base, 0, j );
                        s2.set_control_point( ( pdn + 2 * psplit_base ) / 3, 1, j );
                        s2.set_control_point( ( 2 * pdn + psplit_base ) / 3, 2, j );
                        s2.set_control_point( pdn, 3, j );

                        s4.set_control_point( pdn, 0, j );
                        s4.set_control_point( pdn, 1, j );
                      }
                      else
                      {
                        point_type pfinal;

                        data_type h_target = 0;
                        { // Re-compute local h as target.
                          // Find radius from chord length.
                          data_type radius = base_chord / ( 2.0 * sin( theta / 2.0 ));

                          // Distance to quad Bezier construction point.
                          data_type b = radius * tan( theta / 4.0 );

                          // Extrapolated tip up/down points.
                          point_type ptup, ptdn;
                          ptup = pup + tanup * b * len_factor / ksweep;
                          ptdn = pdn + tandn * b * len_factor / ksweep;

                          // Point on center of circle.
                          point_type pcir = ( ptup + ptdn ) * 0.5;
                          pfinal = pcir;

                          point_type pmid = ( pup + pdn ) * 0.5;
                          h_target = ( pcir - pmid ).norm() * sweepfactor;
                        }

                        data_type hcorr = 0.0;

                        bool collapse = true;
                        if ( extendle )
                        {
                          // Better handle finite LE
                          if ( i == ( nseg - 1 ) )
                          {
                            hcorr = h_end;
                            base_chord = 0;
                          }

                          if ( i >= i_hmax || ( i == ( i_hmax - 1 ) && j > ( deg - 2 ) ) )
                          {
                            h_target = hmax;
                            collapse = false;
                          }
                        }
                        else
                        {
                          if ( i >= i_hmax ) // Collapse to LE
                          {
                            hcorr = hmax;
                            base_chord = 0;
                          }
                        }

                        if ( extendte )
                        {
                          // Better handle finite TE's.
                          if ( i == 0 )
                          {
                            hcorr = h_start;
                            base_chord = 0;
                          }

                          if ( i < i_hmax || ( i == i_hmax && j < 2 ) )
                          {
                            h_target = hmax;
                            collapse = false;
                          }
                        }
                        else
                        {
                          if ( i < i_hmax ) // Collapse to TE
                          {
                            hcorr = hmax;
                            base_chord = 0;
                          }
                        }

                        data_type h = ( hmax - hcorr ) / sweepfactor;

                        data_type w = len_factor * tan( theta / 4.0 ) / ksweep;
                        data_type a = base_chord / 2.0;
                        data_type x = ( h - w * a ) / (1.0 - w / tan( theta / 2.0 ) );

                        data_type e = x / sin( theta / 2.0 );

                        if ( e < 0 )
                        {
                          e = 0.0;
                        }

                        if ( collapse )
                        {
                          e = 0.0;
                        }

                        // Extrapolated extension up/down points.
                        point_type peup, pedn, pesplit;
                        peup = pup + tanup * e;
                        pedn = pdn + tandn * e;
                        pesplit = ( peup + pedn ) * 0.5;

                        s3.set_control_point( pup, 0, j );
                        s3.set_control_point( peup, 1, j );

                        s4.set_control_point( pedn, 0, j );
                        s4.set_control_point( pdn, 1, j );


                        data_type l = ( pesplit - psplit_base ).norm();

                        data_type m = h_target / sweepfactor - l;



                        data_type b = ( m / sin( theta / 2.0 ) ) / ( len_factor / ksweep );


                          // Displacement vector
                        disp = peup - pedn;

                        // Chord of circle to construct
                        data_type circ_chord = disp.norm();

                        if ( tol.approximately_equal( circ_chord, 0.0 ) )
                        {
                          disp << 0, 0, 0;
                          circ_chord = 0.0;
                        }
                        else
                        {
                          disp = disp/circ_chord;
                        }

                        // Re-calculate angles between tangents and displacement
                        dottanup = disp.dot( tanup );
                        thetaup = acos( dottanup );

                        dottandn = disp.dot( tandn );
                        thetadn = acos( dottandn );

                        // Find circular arc to include
                        theta = eli::constants::math<data__>::pi() + thetadn - thetaup;

                        // Find radius from chord length.
                        data_type radius = circ_chord / ( 2.0 * sin( theta / 2.0 ) );

                        // Distance to quad Bezier construction point.
                        // data_type b = radius * tan( theta / 4.0 );

                        radius = b / tan( theta / 4.0 );

                        if ( circ_chord == 0.0 )
                        {
                            b = 0.0;
                            radius = 0.0;
                        }

                        // Extrapolated tip up/down points.
                        point_type ptup, ptdn;
                        ptup = peup + tanup * b * len_factor / ksweep;
                        ptdn = pedn + tandn * b * len_factor / ksweep;

                        // Displacement vector at extrapolated tip
                        point_type dispt = ptup - ptdn;

                        // Point on center of new cap.
                        point_type psplit_cap = ( ptup + ptdn ) * 0.5 + dispt * off_factor;

                        // Fraction of radius to place cubic control point
                        data_type f = eli::constants::math<data_type>::cubic_bezier_circle_const() * tan(theta / 8.0);

                        // Build cubic upper side curve
                        s1.set_control_point( peup, 0, j );
                        s1.set_control_point( peup +   tanup * radius * f * len_factor / ksweep, 1, j );
                        s1.set_control_point( psplit_cap + disp * radius * f, 2, j );
                        s1.set_control_point( psplit_cap, 3, j );

                        // Build cubic lower side curve
                        s2.set_control_point( psplit_cap, 0, j );
                        s2.set_control_point( psplit_cap - disp * radius * f, 1, j );
                        s2.set_control_point( pedn +   tandn * radius * f * len_factor / ksweep, 2, j );
                        s2.set_control_point( pedn, 3, j );
                      }
                    }
                  }
                  break;
                case EDGE:
                  {
                      s1.resize( 1, deg );
                      s2.resize( 1, deg );

                      for ( index_type j = 0; j <= deg; j++ )
                      {
                        point_type pup, psplit, pdn;
                        point_type tanup, tandn;

                        pup = clast.get_control_point( j );
                        pdn = cfirst.get_control_point( j );
                        tanup = ndlast.get_control_point( j );
                        tandn = ndfirst.get_control_point( j );

                        data_type ksweep = 1.0;
                        if ( sweep_correct )
                        {
                          point_type avetan = ( tanup + tandn ) * 0.5;
                          avetan = avetan / avetan.norm();
                          ksweep = 1.0 / sin( acos( avetan.dot( fwd ) ) );

                          if ( tol.approximately_equal( ksweep, 0.0) )
                          {
                            ksweep = 1.0;
                          }
                        }

                        // Displacement vector
                        point_type disp = pup - pdn;

                        // Chord of circle to construct
                        data_type circ_chord = disp.norm();

                        if ( tol.approximately_equal( circ_chord, 0.0 ) )
                        {
                          disp << 0, 0, 0;
                        }
                        else
                        {
                          disp = disp/circ_chord;
                        }

                        // Extrapolated tip up/down points.
                        point_type ptup, ptdn;
                        ptup = pup + tanup * circ_chord * len_factor * ksweep * 0.5;
                        ptdn = pdn + tandn * circ_chord * len_factor * ksweep * 0.5;

                        // Displacement vector at extrapolated tip
                        point_type dispt = ptup - ptdn;

                        psplit = ( ptup + ptdn ) * 0.5 + dispt * off_factor;

                        // Build linear upper side curve
                        s1.set_control_point( pup, 0, j );
                        s1.set_control_point( psplit, 1, j );

                        // Build linear lower side curve
                        s2.set_control_point( psplit, 0, j );
                        s2.set_control_point( pdn, 1, j );

                      }
                    }
                    break;
                  case SHARP:
                    {
                      s1.resize( 2, deg );
                      s2.resize( 2, deg );

                      for ( index_type j = 0; j <= deg; j++ )
                      {
                        point_type pup, psplit, pdn;
                        point_type tanup, tandn;

                        pup = clast.get_control_point( j );
                        pdn = cfirst.get_control_point( j );
                        tanup = ndlast.get_control_point( j );
                        tandn = ndfirst.get_control_point( j );

                        data_type ksweep = 1.0;
                        if ( sweep_correct )
                        {
                          point_type avetan = ( tanup + tandn ) * 0.5;
                          avetan = avetan / avetan.norm();
                          ksweep = 1.0 / sin( acos( avetan.dot( fwd ) ) );

                          if ( tol.approximately_equal( ksweep, 0.0) )
                          {
                            ksweep = 1.0;
                          }
                        }

                        // Displacement vector
                        point_type disp = pup - pdn;

                        // Chord of circle to construct
                        data_type circ_chord = disp.norm();

                        if ( tol.approximately_equal( circ_chord, 0.0 ) )
                        {
                          disp << 0, 0, 0;
                        }
                        else
                        {
                          disp = disp/circ_chord;
                        }

                        // Extrapolated tip up/down points.
                        point_type ptup, ptdn;
                        ptup = pup + tanup * circ_chord * len_factor * ksweep * 0.5;
                        ptdn = pdn + tandn * circ_chord * len_factor * ksweep * 0.5;

                        // Displacement vector at extrapolated tip
                        point_type dispt = ptup - ptdn;

                        psplit = ( ptup + ptdn ) * 0.5 + dispt * off_factor;

                        // Upper/lower control point.
                        point_type pupcp, pdncp;

                        pupcp = pup + tanup * circ_chord * len_factor * ksweep * 0.5 * str_factor;
                        pdncp = pdn + tandn * circ_chord * len_factor * ksweep * 0.5 * str_factor;

                        // Build linear upper side curve
                        s1.set_control_point( pup, 0, j );
                        s1.set_control_point( pupcp, 1, j );
                        s1.set_control_point( psplit, 2, j );

                        // Build linear lower side curve
                        s2.set_control_point( psplit, 0, j );
                        s2.set_control_point( pdncp, 1, j );
                        s2.set_control_point( pdn, 2, j );

                      }
                    }
                    break;
                  }

                }

                if ( cap_type == ROUND_EXT || sphere_tip() )
                {
                  cap.set( s3, 0, i );
                  cap.set( s1, 1, i );
                  cap.set( s2, 2, i );
                  cap.set( s4, 3, i );
                }
                else if ( !profile_rows() )
                {
                  cap.set( s1, 0, i );
                  cap.set( s2, 1, i );
                }

              }
            }

            if ( shear )
            {
              shear_cap( cap, porigin );
            }

            // split the cap surface at mid point
            cap.split_u(surface_copy.get_u0()-delta_param);

            // resize output surface to needed size
            data_type u0, v0;
            std::vector<data_type> ucap_param, vcap_param, uparam, vparam, du, dv;
            index_type i, j, icap_mid, umin_offset(0), vmin_offset(0);

            surface_copy.get_pmap_u(uparam);
            cap.get_pmap_uv(ucap_param, vcap_param);
            // NOTE: This assumes that there are same number of patches on both sides of split
            icap_mid=ucap_param.size()/2;
            switch (edge_to_cap)
            {
              case (CAP_UMIN):
              {
                // establish the new u and v parameterization of surface
                const index_type nucap_patch(icap_mid), nvcap_patch(cap.number_v_patches());
                u0=ucap_param[icap_mid];
                du.resize(uparam.size()-1+nucap_patch);
                for (i=icap_mid; i<static_cast<index_type>(ucap_param.size())-1; ++i)
                {
                  du[i-icap_mid]=ucap_param[i+1]-ucap_param[i];
                }
                for (i=0; i<static_cast<index_type>(uparam.size())-1; ++i)
                {
                  du[i+nucap_patch]=uparam[i+1]-uparam[i];
                }

                // set the new v-parameterization
                vparam.resize(2*vcap_param.size()-1);
                data_type vmid(vcap_param[vcap_param.size()-1]);
                for (j=0; j<static_cast<index_type>(vcap_param.size()); ++j)
                {
                  vparam[j]=vcap_param[j];
                  vparam[vparam.size()-1-j]=2*vmid-vparam[j];
                }

                // set the v-parameters
                v0=vparam[0];
                dv.resize(vparam.size()-1);
                for (j=0; j<static_cast<index_type>(dv.size()); ++j)
                {
                  dv[j]=vparam[j+1]-vparam[j];
                }

                // split a copy original surface at any new v-coordinate locations
                piecewise_surface_type orig_copy(surface_copy);
                typename piecewise_surface_type::surface_type patch;

                for (j=0; j<static_cast<index_type>(vparam.size()); ++j)
                {
                  orig_copy.split_v(vparam[j]);
                }
                assert(orig_copy.number_v_patches()==(static_cast<index_type>(vparam.size())-1));

                // resize the output surface
                ps.init_uv(du.begin(), du.end(), dv.begin(), dv.end(), u0, v0);
                umin_offset=nucap_patch;

                // add first half of cap surfaces
                for (i=0; i<nucap_patch; ++i)
                {
                  for (j=0; j<nvcap_patch; ++j)
                  {
                    cap.get(patch, i + nucap_patch, j);
                    ps.set(patch, i, j);
                  }
                }

                // reverse the cap so that the second half
                cap.reverse_u();
                cap.reverse_v();

                // add second half of cap surfaces
                for (i=0; i<nucap_patch; ++i)
                {
                  for (j=nvcap_patch; j<static_cast<index_type>(dv.size()); ++j)
                  {
                    cap.get(patch, i + nucap_patch, j-nvcap_patch);
                    ps.set(patch, i, j);
                  }
                }

                // add non-cap surfaces
                for (i=0; i<orig_copy.number_u_patches(); ++i)
                {
                  for (j=0; j<orig_copy.number_v_patches(); ++j)
                  {
                    orig_copy.get(patch, i, j);
                    ps.set(patch, i+umin_offset, j+vmin_offset);
                  }
                }

                break;
              }
              case (CAP_UMAX):
              {
                // need to reverse u-direction
                cap.reverse_u();

                // establish the new u and v parameterization of surface
                const index_type nucap_patch(icap_mid), nu_patch(uparam.size()-1), nvcap_patch(cap.number_v_patches());
                u0=uparam[0];
                du.resize(uparam.size()-1+nucap_patch);
                for (i=0; i<static_cast<index_type>(uparam.size())-1; ++i)
                {
                  du[i]=uparam[i+1]-uparam[i];
                }
                for (i=nu_patch; i<static_cast<index_type>(du.size()); ++i)
                {
                  du[i]=ucap_param[(i-nu_patch)+1]-ucap_param[i-nu_patch];
                }

                // set the new v-parameterization
                vparam.resize(2*vcap_param.size()-1);
                data_type vmid(vcap_param[vcap_param.size()-1]);
                for (j=0; j<static_cast<index_type>(vcap_param.size()); ++j)
                {
                  vparam[j]=vcap_param[j];
                  vparam[vparam.size()-1-j]=2*vmid-vparam[j];
                }

                // set the v-parameters
                v0=vparam[0];
                dv.resize(vparam.size()-1);
                for (j=0; j<static_cast<index_type>(dv.size()); ++j)
                {
                  dv[j]=vparam[j+1]-vparam[j];
                }

                // split a copy original surface at any new v-coordinate locations
                piecewise_surface_type orig_copy(surface_copy);
                typename piecewise_surface_type::surface_type patch;

                for (j=0; j<static_cast<index_type>(vparam.size()); ++j)
                {
                  orig_copy.split_v(vparam[j]);
                }
                assert(orig_copy.number_v_patches()==(static_cast<index_type>(vparam.size())-1));

                // resize the output surface
                ps.init_uv(du.begin(), du.end(), dv.begin(), dv.end(), u0, v0);

                // add first half of cap surfaces
                for (i=0; i<nucap_patch; ++i)
                {
                  for (j=0; j<nvcap_patch; ++j)
                  {
                    cap.get(patch, i, j);
                    ps.set(patch, i + nu_patch, j);
                  }
                }

                // reverse the cap so that the second half
                cap.reverse_u();
                cap.reverse_v();

                // add second half of cap surfaces
                for (i=0; i<nucap_patch; ++i)
                {
                  for (j=nvcap_patch; j<static_cast<index_type>(dv.size()); ++j)
                  {
                    cap.get(patch, i, j-nvcap_patch);
                    ps.set(patch, i + nu_patch, j);
                  }
                }

                // add non-cap surfaces
                for (i=0; i<orig_copy.number_u_patches(); ++i)
                {
                  for (j=0; j<orig_copy.number_v_patches(); ++j)
                  {
                    orig_copy.get(patch, i, j);
                    ps.set(patch, i+umin_offset, j+vmin_offset);
                  }
                }

                break;
              }
              case (CAP_VMIN):
              {
                assert(false);
                return false;
                break;
              }
              case (CAP_VMAX):
              {
                assert(false);
                return false;
                break;
              }
              default:
              {
                // must not want any edge capped
                assert(edge_to_cap==CAP_NONE);
                return false;
              }
            }

            return true;
          }

          virtual bool create_curve( piecewise_curve_type &c ) const
          {
            return create_curve( c, cap_type, len_factor, off_factor );
          }

          // Curve where the cap seam comes together, for cap type typ, length factor lf and offset factor of.
          bool create_curve( piecewise_curve_type &c, const index_type &typ, const data_type &lf, const data_type &of ) const
          {
            typedef typename piecewise_curve_type::curve_type curve_type;
            typedef typename piecewise_surface_type::surface_type surface_type;

            c.clear();

            tolerance_type tol;

            if ( !prep_done )
            {
              bool rval = prep();
              if ( !rval )
              {
                return rval;
              }
            }

            data_type tsplit( ( edge.get_t0() + edge.get_tmax() ) * 0.5 );

            // split the curve at the mid-parameter location and reverse the second piece's parameterization
            piecewise_curve_type first_half, second_half;
            edge.split(first_half, second_half, tsplit);
            second_half.reverse();
            second_half.set_t0(first_half.get_t0());

            piecewise_curve_type seam;
            seam.sum( first_half, second_half );
            seam.scale( 0.5 );

            data_type tmidseam = ( seam.get_t0() + seam.get_tmax() ) * 0.5;
            point_type porigin = seam.f( tmidseam );

            // A cap that closes on a point closes on the middle of the edge cap's tip curve.
            point_type ppoint;
            if ( !cap_point( typ, porigin, tmidseam, lf, of, ppoint ) )
            {
              return false;
            }

            // Split the ndelta pseudocurve.  Since this isn't really a curve, it only works because the split
            // point already exists and this becomes reordering and manipulating of existing values, and should
            // not change any values of control points.
            piecewise_curve_type nd_first_half, nd_second_half;
            ndelta.split(nd_first_half, nd_second_half, tsplit);
            nd_second_half.reverse();
            nd_second_half.set_t0(nd_first_half.get_t0());

            std::vector<data_type> pmap, dvcap;
            // Get edge parameter map
            first_half.get_pmap( pmap );
            dvcap.resize( pmap.size() - 1 );

            for ( index_type i = 0; i < pmap.size() - 1; i++ )
            {
              dvcap[i] = pmap[ i + 1 ] - pmap[ i ];
            }

            index_type i_hmax = 0;
            if ( typ == ROUND_EXT )
            {
              i_hmax = seam.find( v_hmax );
            }

            c.set_t0( pmap[0] );

            // create curve where cap seam comes together.
            {
              index_type nseg = first_half.number_segments();

              // Vector along upper curve in average 'forward' direction
              point_type fwd;
              fwd = ( first_half.f( tsplit ) - first_half.f( first_half.get_t0() ) );
              data_type fwdlen = fwd.norm();
              if ( tol.approximately_equal( fwdlen, 0.0 ) )
              {
                fwd << 0, 0, 0;
              }
              else
              {
                fwd = fwd / fwdlen;
              }

              for ( index_type i = 0; i < nseg; i++ )
              {
                curve_type cfirst, clast;
                curve_type ndfirst, ndlast;

                first_half.get( cfirst, i );
                second_half.get( clast, i );
                nd_first_half.get( ndfirst, i );
                nd_second_half.get( ndlast, i );

                index_type deg = cfirst.degree();

                curve_type c1, c2;

                switch( typ ){
                case FLAT:
                  {
                    c1.resize( deg );

                    for ( index_type j = 0; j <= deg; j++ )
                    {
                      point_type pup, psplit, pdn;

                      pup = clast.get_control_point( j );
                      pdn = cfirst.get_control_point( j );

                      psplit = ( pup + pdn ) * 0.5;

                      // Build edge curve
                      c1.set_control_point( psplit, j );
                    }
                  }
                  break;
                case POINT:
                case CONE:
                case OGIVE:
                case PROFILE:
                  {
                    c1.resize( deg );

                    for ( index_type j = 0; j <= deg; j++ )
                    {
                      point_type pup, pdn;

                      pup = clast.get_control_point( j );
                      pdn = cfirst.get_control_point( j );

                      // Build linear upper side curve
                      c1.set_control_point( ppoint + pnt_offset, j );
                    }
                  }
                  break;
                case ROUND:
                  {
                    c1.resize( deg );

                    for ( index_type j = 0; j <= deg; j++ )
                    {
                      point_type pup, psplit, pdn;
                      point_type tanup, tandn;

                      pup = clast.get_control_point( j );
                      pdn = cfirst.get_control_point( j );
                      tanup = ndlast.get_control_point( j );
                      tandn = ndfirst.get_control_point( j );

                      data_type ksweep = 1.0;
                      if ( sweep_correct )
                      {
                        point_type avetan = ( tanup + tandn ) * 0.5;
                        avetan = avetan / avetan.norm();
                        ksweep = 1.0 / sin( acos( avetan.dot( fwd ) ) );

                        if ( tol.approximately_equal( ksweep, 0.0) )
                        {
                          ksweep = 1.0;
                        }
                      }

                      // Displacement vector
                      point_type disp = pup - pdn;

                      // Chord of circle to construct
                      data_type circ_chord = disp.norm();

                      if ( tol.approximately_equal( circ_chord, 0.0 ) )
                      {
                        disp << 0, 0, 0;
                      }
                      else
                      {
                        disp = disp/circ_chord;
                      }

                      // Find angles between tangents and displacement
                      data_type dottanup = disp.dot( tanup );
                      data_type thetaup = acos( dottanup );

                      data_type dottandn = disp.dot( tandn );
                      data_type thetadn = acos( dottandn );

                      // Find circular arc to include
                      data_type theta = eli::constants::math<data__>::pi() + thetadn - thetaup;

                      if ( tol.approximately_equal( theta, 0) ) // Make flat panel avoid division by zero.
                      {
                        psplit = ( pup + pdn ) * 0.5;

                        c1.set_control_point( psplit, j );
                      }
                      else
                      {
                        // Find radius from chord length.
                        data_type radius = circ_chord / ( 2.0 * sin( theta / 2.0 ) );

                        // Distance to quad Bezier construction point.
                        data_type b = radius * tan( theta / 4.0 );

                        // Extrapolated tip up/down points.
                        point_type ptup, ptdn;
                        ptup = pup + tanup * b * lf * ksweep;
                        ptdn = pdn + tandn * b * lf * ksweep;

                        // Displacement vector at extrapolated tip
                        point_type dispt = ptup - ptdn;

                        // Point on center of circle.
                        psplit = ( ptup + ptdn ) * 0.5 + dispt * of;

                        c1.set_control_point( psplit, j );
                      }
                    }

                  break;
                case ROUND_EXT:
                  {
                    c1.resize( deg );

                    for ( index_type j = 0; j <= deg; j++ )
                    {
                      point_type pup, psplit_base, pdn;
                      point_type tanup, tandn;

                      pup = clast.get_control_point( j );
                      pdn = cfirst.get_control_point( j );

                      psplit_base = ( pup + pdn ) * 0.5;

                      tanup = ndlast.get_control_point( j );
                      tandn = ndfirst.get_control_point( j );

                      point_type avetan = ( tanup + tandn ) * 0.5;
                      avetan = avetan / avetan.norm();


                      data_type sweepfactor = sin( acos( avetan.dot( fwd ) ) );

                      if ( tol.approximately_equal( sweepfactor, 0.0) )
                      {
                        sweepfactor = 1.0;
                      }

                      data_type ksweep = 1.0;
                      if ( sweep_correct )
                      {
                        ksweep = sweepfactor;
                      }

                      // Displacement vector
                      point_type disp = pup - pdn;

                      // Chord of circle to construct
                      data_type base_chord = disp.norm();

                      if ( tol.approximately_equal( base_chord, 0.0 ) )
                      {
                        disp << 0, 0, 0;
                      }
                      else
                      {
                        disp = disp/base_chord;
                      }

                      // Find angles between tangents and displacement
                      data_type dottanup = disp.dot( tanup );
                      data_type thetaup = acos( dottanup );

                      data_type dottandn = disp.dot( tandn );
                      data_type thetadn = acos( dottandn );

                      // Find circular arc to include
                      data_type theta = eli::constants::math<data__>::pi() + thetadn - thetaup;

                      if ( tol.approximately_equal( theta, 0) ) // Make flat panel avoid division by zero.
                      {
                        c1.set_control_point( psplit_base, j );
                      }
                      else
                      {
                        point_type pfinal;

                        data_type h_target = 0;
                        { // Re-compute local h as target.
                          // Find radius from chord length.
                          data_type radius = base_chord / ( 2.0 * sin( theta / 2.0 ));

                          // Distance to quad Bezier construction point.
                          data_type b = radius * tan( theta / 4.0 );

                          // Extrapolated tip up/down points.
                          point_type ptup, ptdn;
                          ptup = pup + tanup * b * lf / ksweep;
                          ptdn = pdn + tandn * b * lf / ksweep;

                          // Point on center of circle.
                          point_type pcir = ( ptup + ptdn ) * 0.5;
                          pfinal = pcir;

                          point_type pmid = ( pup + pdn ) * 0.5;
                          h_target = ( pcir - pmid ).norm() * sweepfactor;
                        }

                        data_type hcorr = 0.0;

                        bool collapse = true;
                        if ( extendle )
                        {
                          // Better handle finite LE
                          if ( i == ( nseg - 1 ) )
                          {
                            hcorr = h_end;
                            base_chord = 0;
                          }

                          if ( i >= i_hmax || ( i == ( i_hmax - 1 ) && j > ( deg - 2 ) ) )
                          {
                            h_target = hmax;
                            collapse = false;
                          }
                        }
                        else
                        {
                          if ( i >= i_hmax ) // Collapse to LE
                          {
                            hcorr = hmax;
                            base_chord = 0;
                          }
                        }

                        if ( extendte )
                        {
                          // Better handle finite TE's.
                          if ( i == 0 )
                          {
                            hcorr = h_start;
                            base_chord = 0;
                          }

                          if ( i < i_hmax || ( i == i_hmax && j < 2 ) )
                          {
                            h_target = hmax;
                            collapse = false;
                          }
                        }
                        else
                        {
                          if ( i < i_hmax ) // Collapse to TE
                          {
                            hcorr = hmax;
                            base_chord = 0;
                          }
                        }

                        data_type h = ( hmax - hcorr ) / sweepfactor;

                        data_type w = lf * tan( theta / 4.0 ) / ksweep;
                        data_type a = base_chord / 2.0;
                        data_type x = ( h - w * a ) / (1.0 - w / tan( theta / 2.0 ) );

                        data_type e = x / sin( theta / 2.0 );

                        if ( e < 0 )
                        {
                          e = 0.0;
                        }

                        if ( collapse )
                        {
                          e = 0.0;
                        }

                        // Extrapolated extension up/down points.
                        point_type peup, pedn, pesplit;
                        peup = pup + tanup * e;
                        pedn = pdn + tandn * e;
                        pesplit = ( peup + pedn ) * 0.5;


                        data_type l = ( pesplit - psplit_base ).norm();

                        data_type m = h_target / sweepfactor - l;

                        data_type b = ( m / sin( theta / 2.0 ) ) / ( lf / ksweep );


                          // Displacement vector
                        disp = peup - pedn;

                        // Chord of circle to construct
                        data_type circ_chord = disp.norm();

                        if ( tol.approximately_equal( circ_chord, 0.0 ) )
                        {
                          disp << 0, 0, 0;
                          circ_chord = 0.0;
                        }
                        else
                        {
                          disp = disp/circ_chord;
                        }

                        // Re-calculate angles between tangents and displacement
                        dottanup = disp.dot( tanup );
                        thetaup = acos( dottanup );

                        dottandn = disp.dot( tandn );
                        thetadn = acos( dottandn );

                        // Find circular arc to include
                        theta = eli::constants::math<data__>::pi() + thetadn - thetaup;

                        // Find radius from chord length.
                        data_type radius = circ_chord / ( 2.0 * sin( theta / 2.0 ) );

                        // Distance to quad Bezier construction point.
                        // data_type b = radius * tan( theta / 4.0 );

                        radius = b / tan( theta / 4.0 );

                        if ( circ_chord == 0.0 )
                        {
                            b = 0.0;
                            radius = 0.0;
                        }

                        // Extrapolated tip up/down points.
                        point_type ptup, ptdn;
                        ptup = peup + tanup * b * lf / ksweep;
                        ptdn = pedn + tandn * b * lf / ksweep;

                        // Displacement vector at extrapolated tip
                        point_type dispt = ptup - ptdn;

                        // Point on center of new cap.
                        point_type psplit_cap = ( ptup + ptdn ) * 0.5 + dispt * of;

                        c1.set_control_point( psplit_cap, j );
                      }
                    }
                  }
                  break;
                case EDGE:
                  {
                      c1.resize( deg );

                      for ( index_type j = 0; j <= deg; j++ )
                      {
                        point_type pup, psplit, pdn;
                        point_type tanup, tandn;

                        pup = clast.get_control_point( j );
                        pdn = cfirst.get_control_point( j );
                        tanup = ndlast.get_control_point( j );
                        tandn = ndfirst.get_control_point( j );

                        data_type ksweep = 1.0;
                        if ( sweep_correct )
                        {
                          point_type avetan = ( tanup + tandn ) * 0.5;
                          avetan = avetan / avetan.norm();
                          ksweep = 1.0 / sin( acos( avetan.dot( fwd ) ) );

                          if ( tol.approximately_equal( ksweep, 0.0) )
                          {
                            ksweep = 1.0;
                          }
                        }

                        // Displacement vector
                        point_type disp = pup - pdn;

                        // Chord of circle to construct
                        data_type circ_chord = disp.norm();

                        if ( tol.approximately_equal( circ_chord, 0.0 ) )
                        {
                          disp << 0, 0, 0;
                        }
                        else
                        {
                          disp = disp/circ_chord;
                        }

                        // Extrapolated tip up/down points.
                        point_type ptup, ptdn;
                        ptup = pup + tanup * circ_chord * lf * ksweep * 0.5;
                        ptdn = pdn + tandn * circ_chord * lf * ksweep * 0.5;

                        // Displacement vector at extrapolated tip
                        point_type dispt = ptup - ptdn;

                        psplit = ( ptup + ptdn ) * 0.5 + dispt * of;

                        // Build linear upper side curve
                        c1.set_control_point( psplit, j );
                      }
                    }
                    break;
                  case SHARP:
                    {
                      c1.resize( deg );

                      for ( index_type j = 0; j <= deg; j++ )
                      {
                        point_type pup, psplit, pdn;
                        point_type tanup, tandn;

                        pup = clast.get_control_point( j );
                        pdn = cfirst.get_control_point( j );
                        tanup = ndlast.get_control_point( j );
                        tandn = ndfirst.get_control_point( j );

                        data_type ksweep = 1.0;
                        if ( sweep_correct )
                        {
                          point_type avetan = ( tanup + tandn ) * 0.5;
                          avetan = avetan / avetan.norm();
                          ksweep = 1.0 / sin( acos( avetan.dot( fwd ) ) );

                          if ( tol.approximately_equal( ksweep, 0.0) )
                          {
                            ksweep = 1.0;
                          }
                        }

                        // Displacement vector
                        point_type disp = pup - pdn;

                        // Chord of circle to construct
                        data_type circ_chord = disp.norm();

                        if ( tol.approximately_equal( circ_chord, 0.0 ) )
                        {
                          disp << 0, 0, 0;
                        }
                        else
                        {
                          disp = disp/circ_chord;
                        }

                        // Extrapolated tip up/down points.
                        point_type ptup, ptdn;
                        ptup = pup + tanup * circ_chord * lf * ksweep * 0.5;
                        ptdn = pdn + tandn * circ_chord * lf * ksweep * 0.5;

                        // Displacement vector at extrapolated tip
                        point_type dispt = ptup - ptdn;

                        psplit = ( ptup + ptdn ) * 0.5 + dispt * of;

                        c1.set_control_point( psplit, j );
                      }
                    }
                    break;
                  }

                }

                c.push_back( c1, dvcap[i] );

              }
            }

            return true;
          }

          // Squared length and offset scales along the curve where the cap seam comes together.  The tip
          // is len_factor times a length curve from the tip with no length, plus off_factor times an
          // offset curve.  lensq is the square of the length curve's component across the forward direction,
          // and offsq the square of the offset curve's magnitude.  Both curves can depend on len_factor.
          virtual bool create_scale_sq_curves( oned_piecewise_curve_type &lensq, oned_piecewise_curve_type &offsq ) const
          {
            lensq.clear();
            offsq.clear();

            // The tip with no length: the middle of the uncapped end, which is the flat cap's tip, or for a
            // cap that closes on a point, the point at its middle.
            index_type midtyp = FLAT;
            if ( closes_on_point( cap_type ) )
            {
              midtyp = POINT;
            }
            piecewise_curve_type cmid, ctip, ctipoff;
            if ( !create_curve( cmid, midtyp, 0.0, 0.0 ) ||
                 !create_curve( ctip, cap_type, len_factor, 0.0 ) ||
                 !create_curve( ctipoff, cap_type, len_factor, 1.0 ) )
            {
              return false;
            }
            cmid.scale( -1.0 );

            // The tip moves linearly with off_factor.
            piecewise_curve_type offc, negtip( ctip );
            negtip.scale( -1.0 );
            offc.sum( ctipoff, negtip );

            piecewise_curve_type lenc;
            if ( len_factor > 0 )
            {
              lenc.sum( ctip, cmid );
              lenc.scale( 1.0 / len_factor );
            }
            else if ( cap_type != ROUND_EXT )
            {
              // These tips move linearly with len_factor, so a unit length gives the scale.
              if ( !create_curve( ctip, cap_type, 1.0, 0.0 ) )
              {
                return false;
              }
              lenc.sum( ctip, cmid );
            }
            else
            {
              // The extended round cap's height is set for its length.
              return false;
            }

            if ( length_across_forward() )
            {
              remove_forward_component( lenc );
            }

            piecewise_curve_type sq;
            sq.square( lenc );
            lensq = sq.sumcompcurve();
            sq.square( offc );
            offsq = sq.sumcompcurve();

            return true;
          }

          virtual bool create_sweepfactor_sq_curve( oned_piecewise_curve_type &c ) const
          {
            typedef typename piecewise_curve_type::curve_type curve_type;
            typedef typename piecewise_surface_type::surface_type surface_type;
            typedef typename oned_piecewise_curve_type::curve_type oned_curve_type;
            typedef typename oned_piecewise_curve_type::point_type oned_point_type;

            c.clear();

            tolerance_type tol;

            if ( !prep_done )
            {
              bool rval = prep();
              if ( !rval )
              {
                return rval;
              }
            }

            data_type tsplit( ( edge.get_t0() + edge.get_tmax() ) * 0.5 );

            // split the curve at the mid-parameter location and reverse the second piece's parameterization
            piecewise_curve_type first_half, second_half;
            edge.split(first_half, second_half, tsplit);
            second_half.reverse();
            second_half.set_t0(first_half.get_t0());

            piecewise_curve_type seam;
            seam.sum( first_half, second_half );
            seam.scale( 0.5 );

            data_type tmidseam = ( seam.get_t0() + seam.get_tmax() ) * 0.5;
            point_type porigin = seam.f( tmidseam );

            // Split the ndelta pseudocurve.  Since this isn't really a curve, it only works because the split
            // point already exists and this becomes reordering and manipulating of existing values, and should
            // not change any values of control points.
            piecewise_curve_type nd_first_half, nd_second_half;
            ndelta.split(nd_first_half, nd_second_half, tsplit);
            nd_second_half.reverse();
            nd_second_half.set_t0(nd_first_half.get_t0());

            std::vector<data_type> pmap, dvcap;
            // Get edge parameter map
            first_half.get_pmap( pmap );
            dvcap.resize( pmap.size() - 1 );

            index_type nseg = first_half.number_segments();

            piecewise_oned_cubic_spline_creator_type pcsc( nseg );

            pcsc.set_t0( pmap[0] );

            for ( index_type i = 0; i < pmap.size() - 1; i++ )
            {
              pcsc.set_segment_dt( pmap[ i + 1 ] - pmap[ i ], i );
            }

            // copy points over to new type
            vector< oned_point_type > pts( nseg + 1 );


            // Vector along upper curve in average 'forward' direction
            point_type fwd;
            fwd = ( first_half.f( tsplit ) - first_half.f( first_half.get_t0() ) );
            data_type fwdlen = fwd.norm();
            if ( tol.approximately_equal( fwdlen, 0.0 ) )
            {
              fwd << 0, 0, 0;
            }
            else
            {
              fwd = fwd / fwdlen;
            }

            for ( index_type i = 0; i <= nseg; i++ )
            {
              index_type iseg = i;

              if ( i == nseg )
              {
                iseg--;
              }

              curve_type cfirst, clast;
              curve_type ndfirst, ndlast;

              nd_first_half.get( ndfirst, iseg );
              nd_second_half.get( ndlast, iseg );

              index_type deg = ndfirst.degree();

              index_type j = 0;

              if ( i == nseg )
              {
                j = deg;
              }

              point_type tanup, tandn;

              tanup = ndlast.get_control_point( j );
              tandn = ndfirst.get_control_point( j );

              point_type avetan = ( tanup + tandn ) * 0.5;
              avetan = avetan / avetan.norm();
              data_type sweepfactor = sin( acos( avetan.dot( fwd ) ) );

              if ( tol.approximately_equal( sweepfactor, 0.0) )
              {
                sweepfactor = 1.0;
              }

              pts[i] << sweepfactor;
            }

            pcsc.set_chip( pts.begin(), eli::geom::general::NOT_CONNECTED );

            oned_piecewise_curve_type s;
            pcsc.create( s );

            c.square( s );

            return true;
          }

          // Copy the surface, and split and promote it so the two halves of the edge to cap match segment
          // for segment.
          virtual bool prep_surface() const
          {
            typedef typename piecewise_curve_type::curve_type curve_type;
            typedef typename piecewise_surface_type::surface_type surface_type;

            tolerance_type tol;

            surface_copy = orig_surface;

            // Ensure opposing surface faces have mid split and matching parameter splits
            switch (edge_to_cap)
            {
              case (CAP_UMIN):
              case (CAP_UMAX):
              {
                data_type vmin(surface_copy.get_v0()), vmax(surface_copy.get_vmax()), vsplit((vmin+vmax)/2);

                // Make sure split exists.
                surface_copy.split_v( vsplit );

                piecewise_curve_type edg;

                if ( edge_to_cap == CAP_UMIN )
                {
                  surface_copy.get_umin_bndy_curve( edg );
                }
                else
                {
                  surface_copy.get_umax_bndy_curve( edg );
                }

                // make sure that the split location is different than the start/end location otherwise
                // the edge is a point and no need to cap
                if (tol.approximately_equal(edg.f(vmin), edg.f(vsplit)))
                {
                  return false;
                }

                std::vector<data_type> pmap, rpmap, cmap;
                // Get original parameter map
                edg.get_pmap( pmap );

                // Get reversed parameter map
                edg.reverse();
                edg.get_pmap( rpmap );

                tolerance_type ttol;

                // Comparison function for set_union.
                auto comp = [&ttol](const data_type &x1, const data_type &x2)->bool
                {
                  return ttol.approximately_less_than(x1, x2);
                };

                // Place union of pmap and omap into cmap.
                std::set_union( pmap.begin(), pmap.end(), rpmap.begin(), rpmap.end(), std::back_inserter(cmap), comp );

                for ( index_type i = 0; i < cmap.size(); i++ )
                {
                  surface_copy.split_v( cmap[i] );
                }
                break;
              }
              case (CAP_VMIN):
              case (CAP_VMAX):
              {
                data_type umin(surface_copy.get_u0()), umax(surface_copy.get_umax()), usplit((umin+umax)/2);

                // Make sure split exists.
                surface_copy.split_u( usplit );

                piecewise_curve_type edg;

                if ( edge_to_cap == CAP_VMIN )
                {
                  surface_copy.get_vmin_bndy_curve( edg );
                }
                else
                {
                  surface_copy.get_vmax_bndy_curve( edg );
                }

                // make sure that the split location is different than the start/end location otherwise
                // the edge is a point and no need to cap
                if (tol.approximately_equal(edg.f(umin), edg.f(usplit)))
                {
                  return false;
                }

                std::vector<data_type> pmap, rpmap, cmap;
                // Get original parameter map
                edg.get_pmap( pmap );

                // Get reversed parameter map
                edg.reverse();
                edg.get_pmap( rpmap );

                tolerance_type ttol;

                // Comparison function for set_union.
                auto comp = [&ttol](const data_type &x1, const data_type &x2)->bool
                {
                  return ttol.approximately_less_than(x1, x2);
                };

                // Place union of pmap and omap into cmap.
                std::set_union( pmap.begin(), pmap.end(), rpmap.begin(), rpmap.end(), std::back_inserter(cmap), comp );

                for ( index_type i = 0; i < cmap.size(); i++ )
                {
                  surface_copy.split_u( cmap[i] );
                }
                break;
              }
              default:
              {
                // must not want any edge capped
                assert(edge_to_cap==CAP_NONE);
                return false;
              }
            }

            // Ensure opposing surface faces have matching order
            switch (edge_to_cap)
            {
              case (CAP_UMIN):
              case (CAP_UMAX):
              {
                piecewise_curve_type edg;
                surface_copy.get_umin_bndy_curve( edg );

                std::vector<index_type> ord;
                edg.degrees( std::back_inserter( ord ) );

                index_type n = ord.size();
                for ( index_type i = 0; i < n/2; i++ )
                {
                  index_type j = n-i-1;
                  ord[i] = std::max( ord[i], ord[j] );
                  ord[j] = ord[i];
                }

                surface_copy.promote_v_to( ord );

                break;
              }
              case (CAP_VMIN):
              case (CAP_VMAX):
              {
                piecewise_curve_type edg;
                surface_copy.get_vmin_bndy_curve( edg );

                std::vector<index_type> ord;
                edg.degrees( std::back_inserter( ord ) );

                index_type n = ord.size();
                for ( index_type i = 0; i < n/2; i++ )
                {
                  index_type j = n-i-1;
                  ord[i] = std::max( ord[i], ord[j] );
                  ord[j] = ord[i];
                }

                surface_copy.promote_u_to( ord );
                break;
              }
              default:
              {
                // must not want any edge capped
                assert(edge_to_cap==CAP_NONE);
                return false;
              }
            }

            surface_prep_done = true;

            return true;
          }

          // Find the edge to cap and the unit vectors across it, from the edge's two rows of control points
          // alone.  The rows are split and promoted as prep_surface splits and promotes the surface.
          virtual bool prep() const
          {
            typedef typename piecewise_curve_type::curve_type curve_type;

            tolerance_type tol;

            piecewise_curve_type edgnext;
            if ( curves_given )
            {
              edge = orig_edge;
              edgnext = orig_edge_next;
            }
            else if ( !get_edge_curves( orig_surface, edge_to_cap, edge, edgnext ) )
            {
              return false;
            }

            // Ensure opposing halves of the edge have mid split and matching parameter splits
            data_type tmin( edge.get_t0() ), tmax( edge.get_tmax() ), tsplit( ( tmin + tmax ) / 2 );

            // Make sure split exists.
            edge.split( tsplit );
            edgnext.split( tsplit );

            // make sure that the split location is different than the start/end location otherwise
            // the edge is a point and no need to cap
            if (tol.approximately_equal(edge.f(tmin), edge.f(tsplit)))
            {
              return false;
            }

            std::vector<data_type> pmap, rpmap, cmap;
            // Get original parameter map
            edge.get_pmap( pmap );

            // Get reversed parameter map
            piecewise_curve_type redg( edge );
            redg.reverse();
            redg.get_pmap( rpmap );

            tolerance_type ttol;

            // Comparison function for set_union.
            auto comp = [&ttol](const data_type &x1, const data_type &x2)->bool
            {
              return ttol.approximately_less_than(x1, x2);
            };

            // Place union of pmap and omap into cmap.
            std::set_union( pmap.begin(), pmap.end(), rpmap.begin(), rpmap.end(), std::back_inserter(cmap), comp );

            for ( index_type i = 0; i < cmap.size(); i++ )
            {
              edge.split( cmap[i] );
              edgnext.split( cmap[i] );
            }

            // Ensure opposing halves of the edge have matching order
            std::vector<index_type> ord;
            edge.degrees( std::back_inserter( ord ) );

            index_type n = ord.size();
            for ( index_type i = 0; i < n/2; i++ )
            {
              index_type j = n-i-1;
              ord[i] = std::max( ord[i], ord[j] );
              ord[j] = ord[i];
            }

            for ( index_type i = 0; i < n; i++ )
            {
              edge.degree_promote_to( i, ord[i] );
              edgnext.degree_promote_to( i, ord[i] );
            }

            // Unit vectors from the next row out to the edge, control point by control point.
            ndelta = edge;
            for ( index_type i = 0; i < n; i++ )
            {
              curve_type ce, cn;
              edge.get( ce, i );
              edgnext.get( cn, i );
              for ( index_type j = 0; j <= ce.degree(); j++ )
              {
                point_type p = ce.get_control_point( j ) - cn.get_control_point( j );

                data_type len = p.norm();
                if ( tol.approximately_equal(len, 0) )
                {
                  p << 0, 0, 0;
                  // A surface of no length across the end leaves it along the frame's axis, outward as given.
                  if ( frame_given )
                  {
                    p = frame_axis;
                  }
                }
                else
                {
                  p = p / len;
                }
                ce.set_control_point( p, j );
              }
              ndelta.replace( ce, i );
            }

            // The edge less the next row, control point by control point.  Along the end it is the
            // surface's derivative across the end, outward, to a positive factor.
            piecewise_curve_type negnext( edgnext );
            negnext.scale( -1.0 );
            edelta.sum( edge, negnext );

            prep_done = true;

            return true;
          }

          // The control points of the row of patches along an edge, and of the row next to it, as curves.
          // A patch of degree zero across the edge gives its edge row for both.
          static bool get_edge_curves( const piecewise_surface_type &s, edge_cap_identifier ec, piecewise_curve_type &e,
                                       piecewise_curve_type &n )
          {
            typedef typename piecewise_curve_type::curve_type curve_type;
            typedef typename piecewise_surface_type::surface_type surface_type;

            e.clear();
            n.clear();

            index_type nu( s.number_u_patches() ), nv( s.number_v_patches() );
            if ( nu == 0 || nv == 0 )
            {
              return false;
            }

            switch ( ec )
            {
              case ( CAP_UMIN ):
              case ( CAP_UMAX ):
              {
                index_type iu = 0;
                if ( ec == CAP_UMAX )
                {
                  iu = nu - 1;
                }

                e.set_t0( s.get_v0() );
                n.set_t0( s.get_v0() );
                for ( index_type iv = 0; iv < nv; iv++ )
                {
                  surface_type patch;
                  data_type du, dv;
                  s.get( patch, du, dv, iu, iv );

                  index_type deg( patch.degree_v() ), ie( 0 ), in( 1 );
                  if ( ec == CAP_UMAX )
                  {
                    ie = patch.degree_u();
                    in = ie - 1;
                  }
                  if ( patch.degree_u() <= 0 )
                  {
                    in = ie;
                  }

                  curve_type ce( deg ), cn( deg );
                  for ( index_type j = 0; j <= deg; j++ )
                  {
                    ce.set_control_point( patch.get_control_point( ie, j ), j );
                    cn.set_control_point( patch.get_control_point( in, j ), j );
                  }
                  e.push_back( ce, dv );
                  n.push_back( cn, dv );
                }
                break;
              }
              case ( CAP_VMIN ):
              case ( CAP_VMAX ):
              {
                index_type iv = 0;
                if ( ec == CAP_VMAX )
                {
                  iv = nv - 1;
                }

                e.set_t0( s.get_u0() );
                n.set_t0( s.get_u0() );
                for ( index_type iu = 0; iu < nu; iu++ )
                {
                  surface_type patch;
                  data_type du, dv;
                  s.get( patch, du, dv, iu, iv );

                  index_type deg( patch.degree_u() ), je( 0 ), jn( 1 );
                  if ( ec == CAP_VMAX )
                  {
                    je = patch.degree_v();
                    jn = je - 1;
                  }
                  if ( patch.degree_v() <= 0 )
                  {
                    jn = je;
                  }

                  curve_type ce( deg ), cn( deg );
                  for ( index_type i = 0; i <= deg; i++ )
                  {
                    ce.set_control_point( patch.get_control_point( i, je ), i );
                    cn.set_control_point( patch.get_control_point( i, jn ), i );
                  }
                  e.push_back( ce, du );
                  n.push_back( cn, du );
                }
                break;
              }
              default:
              {
                // must not want any edge capped
                assert(ec==CAP_NONE);
                return false;
              }
            }

            return true;
          }

          // Unit vector along the upper half of the end edge, as create_curve finds it.  Zero for a
          // collapsed edge.
          point_type forward_direction() const
          {
            tolerance_type tol;

            data_type tsplit( ( edge.get_t0() + edge.get_tmax() ) * 0.5 );

            piecewise_curve_type first_half, second_half;
            edge.split(first_half, second_half, tsplit);

            // Vector along upper curve in average 'forward' direction
            point_type fwd;
            fwd = ( first_half.f( tsplit ) - first_half.f( first_half.get_t0() ) );
            data_type fwdlen = fwd.norm();
            if ( tol.approximately_equal( fwdlen, 0.0 ) )
            {
              fwd << 0, 0, 0;
            }
            else
            {
              fwd = fwd / fwdlen;
            }
            return fwd;
          }

          // Whether a cap's length is the tip's height across the forward direction, which the extended round
          // cap holds.
          bool length_across_forward() const
          {
            // In a frame a cap that closes on a point moves along the frame's axis, and its move is the length.
            if ( frame_given && closes_on_point( cap_type ) )
            {
              return false;
            }
            return true;
          }

          // Remove the forward direction's component from a curve of displacements, leaving their height
          // across it.  The forward direction is constant, so this is exact along the whole curve.
          void remove_forward_component( piecewise_curve_type &c ) const
          {
            typedef typename piecewise_curve_type::curve_type curve_type;

            if ( !prep_done )
            {
              prep();
            }

            point_type fwd = forward_direction();
            for ( index_type i = 0; i < c.number_segments(); i++ )
            {
              curve_type cseg;
              c.get( cseg, i );
              for ( index_type j = 0; j <= cseg.degree(); j++ )
              {
                point_type cp = cseg.get_control_point( j );
                cseg.set_control_point( cp - fwd * fwd.dot( cp ), j );
              }
              c.replace( cseg, i );
            }
          }

          // The sphere at a cone cap's tip, which touches the point ptip.  Its center pcen is rad back from
          // ptip along axis, the unit vector the point moves along with length.  rad is rad_factor times the
          // length scale.
          // porigin is the middle of the end's seam, at tmidseam.  pdir points from the center to the pole,
          // where the cap's sides gather: along the line through porigin and the center, so that line runs
          // through the middle of the end, on the side ahead along the axis.  With no offset the pole is
          // ptip.
          bool cone_sphere( const point_type &porigin, const data_type &tmidseam, const point_type &ptip,
                            point_type &pcen, data_type &rad, point_type &axis, point_type &pdir ) const
          {
            tolerance_type tol;

            pcen = ptip;
            rad = 0;
            axis << 0, 0, 0;
            pdir = axis;

            if ( frame_given )
            {
              axis = frame_outward();
            }
            else
            {
              piecewise_curve_type cunit;
              if ( !create_curve( cunit, EDGE, 1.0, 0.0 ) )
              {
                return false;
              }

              point_type v = cunit.f( tmidseam ) - porigin;
              data_type vlen = v.norm();
              if ( tol.approximately_equal( vlen, 0.0 ) )
              {
                return true;
              }
              axis = v / vlen;
            }

            rad = rad_factor * length_scale( porigin, tmidseam );

            pcen = ptip - axis * rad;

            pdir = axis;
            point_type away = pcen - porigin;
            data_type alen = away.norm();
            if ( !tol.approximately_equal( alen, 0.0 ) )
            {
              pdir = away / alen;
              if ( pdir.dot( axis ) < 0 )
              {
                pdir = -pdir;
              }
            }
            return true;
          }

          // The point a cap of type typ closes on, for length and offset factors lf and of, from porigin, the
          // middle of the end's seam at tmidseam: the middle of the edge cap's tip curve.  A cap that does not
          // close on a point has porigin.
          // In a frame the point is along the frame's outward axis and its up direction, each times the end's
          // size across the axis.
          bool cap_point( const index_type &typ, const point_type &porigin, const data_type &tmidseam,
                          const data_type &lf, const data_type &of, point_type &ppoint ) const
          {
            ppoint = porigin;
            if ( !closes_on_point( typ ) )
            {
              return true;
            }

            if ( frame_given )
            {
              data_type sc = frame_scale( porigin );
              ppoint = porigin + frame_outward() * ( lf * sc ) + frame_up * ( of * sc );
              return true;
            }

            piecewise_curve_type cedge;
            if ( !create_curve( cedge, EDGE, lf, of ) )
            {
              return false;
            }
            ppoint = cedge.f( tmidseam );
            return true;
          }

          // Shear the cap built with no offset: each control point moves along the frame's up direction by the
          // offset times its distance from the plane of the end over the point's.  The end stays where it is,
          // even where it is not square to the axis, and the point moves by the offset.  The shear is linear,
          // so the surface moves with its control points.
          void shear_cap( piecewise_surface_type &cap, const point_type &porigin ) const
          {
            typedef typename piecewise_surface_type::surface_type surface_type;

            tolerance_type tol;

            data_type sc = frame_scale( porigin );
            data_type len = len_factor * sc;
            if ( !( len > 0 ) )
            {
              return;
            }
            point_type out = frame_outward();
            point_type nend = end_normal( out );

            // The point, built with no offset, is len along out from porigin.
            data_type ptlen = len * nend.dot( out );
            if ( tol.approximately_equal( ptlen, 0.0 ) )
            {
              return;
            }
            point_type step = frame_up * ( off_factor * sc / ptlen );

            for ( index_type iu = 0; iu < cap.number_u_patches(); iu++ )
            {
              for ( index_type iv = 0; iv < cap.number_v_patches(); iv++ )
              {
                surface_type patch;
                cap.get( patch, iu, iv );
                for ( index_type i = 0; i <= patch.degree_u(); i++ )
                {
                  for ( index_type j = 0; j <= patch.degree_v(); j++ )
                  {
                    point_type cp = patch.get_control_point( i, j );
                    patch.set_control_point( cp + step * nend.dot( cp - porigin ), i, j );
                  }
                }
                cap.set( patch, iu, iv );
              }
            }
          }

          // The frame's axis, pointing away from the body: the side of the end the surface's direction across
          // it leans to, its derivative across the end summed along it, or the axis as given where that is
          // square to the axis.
          point_type frame_outward() const
          {
            index_type n = 64;
            data_type t0 = edelta.get_t0();
            data_type dt = ( edelta.get_tmax() - t0 ) / n;

            point_type sum;
            sum << 0, 0, 0;
            for ( index_type k = 0; k < n; k++ )
            {
              point_type d = edelta.f( t0 + ( k + 0.5 ) * dt );
              data_type dlen = d.norm();
              if ( dlen > 0 )
              {
                sum = sum + d / dlen;
              }
            }

            if ( frame_axis.dot( sum ) < -1e-3 * n )
            {
              return -frame_axis;
            }
            return frame_axis;
          }

          // The end's size across the frame's axis: the root mean square distance of points along it from pend,
          // its middle, square to the axis.
          data_type frame_scale( const point_type &pend ) const
          {
            index_type n = 256;
            data_type t0 = edge.get_t0();
            data_type dt = ( edge.get_tmax() - t0 ) / n;

            data_type sum = 0;
            for ( index_type k = 0; k < n; k++ )
            {
              point_type d = edge.f( t0 + ( k + 0.5 ) * dt ) - pend;
              d = d - frame_axis * frame_axis.dot( d );
              sum += d.squaredNorm();
            }
            return std::sqrt( sum / n );
          }

          // The cap's length scale: in a frame the end's size across its axis, or else the height across the
          // forward direction of the point's move per unit length.  porigin is the middle of the end's seam, at
          // tmidseam.
          data_type length_scale( const point_type &porigin, const data_type &tmidseam ) const
          {
            if ( frame_given )
            {
              return frame_scale( porigin );
            }

            piecewise_curve_type cunit;
            if ( !create_curve( cunit, EDGE, 1.0, 0.0 ) )
            {
              return 0;
            }

            point_type v = cunit.f( tmidseam ) - porigin;
            point_type fwd = forward_direction();
            return ( v - fwd * fwd.dot( v ) ).norm();
          }

          // The end's root mean square distance from pend, its middle, over points along it.
          data_type end_radius( const point_type &pend ) const
          {
            index_type n = 256;
            data_type t0 = edge.get_t0();
            data_type dt = ( edge.get_tmax() - t0 ) / n;

            data_type sum = 0;
            for ( index_type k = 0; k < n; k++ )
            {
              sum += ( edge.f( t0 + ( k + 0.5 ) * dt ) - pend ).squaredNorm();
            }
            return std::sqrt( sum / n );
          }

          // The direction pdir of the circular arc from porigin, leaving it along unit vector d, where it
          // reaches pcen.  Where pcen lies on the line through porigin along d, pdir is d.
          void ogive_pole_direction( const point_type &porigin, const point_type &d, const point_type &pcen,
                                     point_type &pdir ) const
          {
            tolerance_type tol;

            point_type c = pcen - porigin;
            data_type clen = c.norm();
            if ( tol.approximately_equal( clen, 0.0 ) )
            {
              return;
            }
            c = c / clen;
            pdir = c * ( 2.0 * c.dot( d ) ) - d;
          }

          // An ogive cap's pole direction pdir, from the arc that leaves porigin along dmean and reaches the
          // sphere's center pcen, and its top and bottom sides otop and obot, from the middles of the end's
          // halves, which the other sides are blended from.
          // A truncated cap's sides run to its face, and it has neither.
          template < typename curve__ >
          void ogive_pole_and_sides( const curve__ &first_half, const curve__ &second_half,
                                     const curve__ &ed_first_half, const curve__ &ed_second_half,
                                     const data_type &tmidseam, const point_type &pnorm, const point_type &porigin,
                                     const point_type &dmean, const point_type &pcen, const data_type &rad,
                                     point_type &pdir, std::vector < point_type > &otop,
                                     std::vector < point_type > &obot ) const
          {
            if ( truncated() )
            {
              return;
            }

            ogive_pole_direction( porigin, dmean, pcen, pdir );
            ogive_reference_side( second_half, ed_second_half, tmidseam, pnorm, porigin, pcen, rad, pdir, otop );
            ogive_reference_side( first_half, ed_first_half, tmidseam, pnorm, porigin, pcen, rad, pdir, obot );
          }

          // The side of a cone cap from point q of the end: the point ptan where the line from q touches the
          // sphere of center pcen and radius rad, and the inner control points parc1 and parc2 of the cubic
          // arc of the sphere from ptan to the pole, pcen + pdir * rad.  The side lies in the plane of q and
          // the line through the pole and pcen, which runs through pend, the middle of the end, so the sides
          // from different points of the end never cross.  In that plane ( ptan - pcen ) . ( ptan - q ) = 0
          // has two roots, and the side takes one that leads on to the pole: on q's side of the line
          // through the pole and pcen, reached running along the meridian towards the pole, the one on
          // the outward side of the line from pcen to q, away from pend, where both do.  Where neither
          // does, the side runs straight to the pole.  The arc is the sphere's meridian from ptan to the
          // pole.  From a point inside the sphere the side runs straight out to it, and with no radius it runs
          // straight to the tip, ptan halfway along.
          // Where the meridian would turn through little, sphere_split parts the arc from the side further
          // back, and sphere_arc runs it on to the pole from there.
          void cone_side( const point_type &q, const point_type &pend, const point_type &pcen, const data_type &rad,
                          const point_type &pdir, point_type &ptan, point_type &parc1, point_type &parc2 ) const
          {
            tolerance_type tol;

            point_type ptip = pcen + pdir * rad;

            point_type dq = q - pcen;
            data_type dlen = dq.norm();
            if ( tol.approximately_equal( rad, 0.0 ) || tol.approximately_equal( dlen, 0.0 ) )
            {
              ptan = ( q + ptip ) * 0.5;
              parc1 = ( ptan * 2.0 + ptip ) / 3.0;
              parc2 = ( ptan + ptip * 2.0 ) / 3.0;
              return;
            }
            dq = dq / dlen;

            // Unit vector from the center to the tangent point.
            point_type u = dq;
            if ( dlen > rad )
            {
              point_type e = ( q - pend ) - dq * dq.dot( q - pend );
              data_type elen = e.norm();
              if ( tol.approximately_equal( elen, 0.0 ) )
              {
                e = pdir - dq * dq.dot( pdir );
                elen = e.norm();
              }
              if ( !tol.approximately_equal( elen, 0.0 ) )
              {
                e = e / elen;

                // The radius to the tangent point makes angle alpha with the line to q, cos( alpha ) = rad / |q - pcen|.
                data_type ca = rad / dlen;
                data_type sa = std::sqrt( 1.0 - ca * ca );
                // Each root and the direction from q to it, which stays well set as the roots close on q.
                u = dq * ca + e * sa;
                if ( !cone_root_leads( q, u, dq * ( -dlen * sa ) + e * rad, pcen, pdir ) )
                {
                  u = dq * ca - e * sa;
                  if ( !cone_root_leads( q, u, dq * ( -dlen * sa ) - e * rad, pcen, pdir ) )
                  {
                    u = pdir;
                  }
                }
              }
            }
            ptan = pcen + u * rad;

            if ( dlen <= rad )
            {
              sphere_arc( u, pcen, rad, pdir, parc1, parc2 );
              return;
            }

            point_type o1 = ( q * 2.0 + ptan ) / 3.0;
            point_type o2 = ( q + ptan * 2.0 ) / 3.0;
            point_type dt;
            sphere_split( q, o1, o2, ptan, dt, pcen, rad, pdir );
            sphere_arc( ptan, dt, pcen, rad, pdir, parc1, parc2 );
          }

          // Whether the point of the sphere of center pcen in unit direction u from it, where the line from q
          // along tq touches it, leads on to the pole, in direction pdir from the center: on q's side of the
          // line through the pole and pcen, and reached running along the meridian towards the pole.
          bool cone_root_leads( const point_type &q, const point_type &u, const point_type &tq, const point_type &pcen,
                                const point_type &pdir ) const
          {
            point_type eo = ( q - pcen ) - pdir * pdir.dot( q - pcen );
            if ( u.dot( eo ) < 0 )
            {
              return false;
            }
            point_type m = pdir - u * u.dot( pdir );
            return tq.dot( m ) >= 0;
          }

          // Where the arc of the sphere of center pcen and radius rad from a side's tangent point to the pole,
          // pcen + pdir * rad, would turn through less than a thirty-second of a half turn, the side's split
          // moves back along its first part, the cubic q, o1, o2, ptan, by the sphere's radius times the
          // shortfall, or at most halfway, so the arc's patch never collapses.  o1, o2 and ptan become the
          // first part up to the split, and dt is the unit direction it reaches the split along.
          void sphere_split( const point_type &q, point_type &o1, point_type &o2, point_type &ptan, point_type &dt,
                             const point_type &pcen, const data_type &rad, const point_type &pdir ) const
          {
            tolerance_type tol;

            data_type pi = eli::constants::math<data_type>::pi();
            data_type tmin = pi / 32.0;

            dt = ptan - o2;
            if ( tol.approximately_equal( dt.norm(), 0.0 ) )
            {
              dt = ptan - q;
            }
            if ( !tol.approximately_equal( dt.norm(), 0.0 ) )
            {
              dt = dt / dt.norm();
            }

            point_type v = ptan - pcen;
            data_type theta = 0;
            if ( !tol.approximately_equal( v.norm(), 0.0 ) )
            {
              theta = unit_angle( v / v.norm(), pdir );
            }
            data_type plen = ( o1 - q ).norm() + ( o2 - o1 ).norm() + ( ptan - o2 ).norm();
            if ( theta >= tmin || tol.approximately_equal( plen, 0.0 ) )
            {
              return;
            }

            data_type s = 1.0 - std::min( rad * ( tmin - theta ) / plen, static_cast<data_type>( 0.5 ) );
            point_type a = q + ( o1 - q ) * s;
            point_type b = o1 + ( o2 - o1 ) * s;
            point_type c = o2 + ( ptan - o2 ) * s;
            point_type ab = a + ( b - a ) * s;
            point_type bc = b + ( c - b ) * s;
            point_type split = ab + ( bc - ab ) * s;

            o1 = a;
            o2 = ab;
            ptan = split;
            dt = split - ab;
            if ( tol.approximately_equal( dt.norm(), 0.0 ) )
            {
              dt = bc - split;
            }
            if ( !tol.approximately_equal( dt.norm(), 0.0 ) )
            {
              dt = dt / dt.norm();
            }
          }

          // The inner control points parc1 and parc2 of the cubic arc of the sphere of center pcen and radius
          // rad along its meridian from the point in unit direction u from the center to the pole,
          // pcen + pdir * rad.
          void sphere_arc( const point_type &u, const point_type &pcen, const data_type &rad, const point_type &pdir,
                           point_type &parc1, point_type &parc2 ) const
          {
            tolerance_type tol;

            point_type ptan = pcen + u * rad;
            point_type ptip = pcen + pdir * rad;

            // The meridian's tangent at the point, towards the pole.
            data_type cth = std::max( std::min( u.dot( pdir ), static_cast<data_type>( 1 ) ), static_cast<data_type>( -1 ) );
            point_type ttan = pdir - u * cth;
            data_type alen = ttan.norm();
            if ( tol.approximately_equal( alen, 0.0 ) )
            {
              parc1 = ( ptan * 2.0 + ptip ) / 3.0;
              parc2 = ( ptan + ptip * 2.0 ) / 3.0;
              return;
            }

            sphere_arc( ptan, ttan / alen, pcen, rad, pdir, parc1, parc2 );
          }

          // The inner control points parc1 and parc2 of the cubic arc of the circle from ptan, which it leaves
          // along the unit direction dt, to the pole, pcen + pdir * rad.  With ptan on the sphere of center
          // pcen and radius rad and dt along it, the circle is where the plane of ptan, dt and the pole cuts
          // the sphere, and with dt along ptan's meridian it is the meridian.
          void sphere_arc( const point_type &ptan, const point_type &dt, const point_type &pcen, const data_type &rad,
                           const point_type &pdir, point_type &parc1, point_type &parc2 ) const
          {
            tolerance_type tol;

            point_type ptip = pcen + pdir * rad;
            point_type c = ptip - ptan;
            data_type clen = c.norm();
            if ( tol.approximately_equal( clen, 0.0 ) )
            {
              parc1 = ptan;
              parc2 = ptip;
              return;
            }
            point_type cu = c / clen;

            // The arc turns through psi, and reaches the pole along dt reflected in the chord.  Its handles
            // are those of no more than three quarters of a turn.
            data_type psi = 2.0 * unit_angle( dt, cu );
            point_type dtip = cu * ( 2.0 * cu.dot( dt ) ) - dt;
            data_type h = arc_handle( std::min( psi, 1.5 * eli::constants::math<data_type>::pi() ), clen );
            parc1 = ptan + dt * h;
            parc2 = ptip - dtip * h;
          }

          // The angle between unit vectors a and b.
          data_type unit_angle( const point_type &a, const point_type &b ) const
          {
            data_type c = std::max( std::min( a.dot( b ), static_cast<data_type>( 1 ) ), static_cast<data_type>( -1 ) );
            return std::acos( c );
          }

          // The length of the inner control polygon legs of the cubic arc of a circle through angle psi with
          // a chord of length clen.
          data_type arc_handle( const data_type &psi, const data_type &clen ) const
          {
            data_type s = std::sin( 0.5 * psi );
            if ( std::abs( s ) < 1e-8 )
            {
              return clen / 3.0;
            }
            data_type cbc = eli::constants::math<data_type>::cubic_bezier_circle_const();
            return cbc * std::tan( 0.25 * psi ) * clen / ( 2.0 * s );
          }

          // The sides of a cone cap from segment c of the end, whose middle is pend.  cone_side finds the
          // tangent point and arc control points exactly from points along the segment, which lie on the
          // end, and fit_side_rows fits the rows to them.
          template < typename curve__ >
          void cone_sides( const curve__ &c, const point_type &pend, const point_type &pcen, const data_type &rad,
                           const point_type &pdir,
                           std::vector < point_type > &ptan, std::vector < point_type > &parc1,
                           std::vector < point_type > &parc2 ) const
          {
            index_type nsamp = 8 * c.degree() + 1;

            std::vector < data_type > ts( nsamp );
            std::vector < std::vector < point_type > > vals( 3, std::vector < point_type >( nsamp ) );
            for ( index_type i = 0; i < nsamp; i++ )
            {
              ts[i] = static_cast<data_type>( i ) / ( nsamp - 1 );
              point_type q = c.f( ts[i] );
              cone_side( q, pend, pcen, rad, pdir, vals[0][i], vals[1][i], vals[2][i] );
            }

            std::vector < std::vector < point_type > > rows;
            fit_side_rows( c, pend, pdir, ts, vals, false, rows );
            ptan = rows[0];
            parc1 = rows[1];
            parc2 = rows[2];
          }

          // The direction d an ogive cap's side leaves point q of the end along: dref, the end's reference
          // direction at q, made square to qt, the end's tangent there, and turned about qt by angle ang,
          // towards pend for ang positive.  The surface's tangent plane at q holds d and qt.  Returns false
          // where dref runs along qt.
          bool ogive_direction( const point_type &q, const point_type &qt, const point_type &dref, const data_type &ang,
                                const point_type &pend, point_type &d ) const
          {
            tolerance_type tol;

            data_type qtlen = qt.norm();
            if ( tol.approximately_equal( qtlen, 0.0 ) )
            {
              return false;
            }
            point_type qtu = qt / qtlen;
            point_type dp = dref - qtu * qtu.dot( dref );
            data_type dplen = dp.norm();
            if ( tol.approximately_equal( dplen, 0.0 ) )
            {
              return false;
            }
            dp = dp / dplen;
            point_type nout = dp.cross( qtu );
            if ( nout.dot( q - pend ) < 0 )
            {
              nout = -nout;
            }
            d = dp * std::cos( ang ) - nout * std::sin( ang );
            return true;
          }

          // The side of an ogive cap from point q of the end: a circular arc from q to the sphere of center
          // pcen and radius rad, then an arc of the sphere on to the pole, pcen + pdir * rad.  The circular
          // arc leaves q along d, from ogive_direction, in the plane of q, d and pcen.  The circles through q
          // along d that touch the sphere's great circle in that plane come in two families, each holding one
          // of the two straight lines from q that touch it.  The arc is the one of the family holding the
          // line that has the sphere on its side towards pend, the middle of the end, as a cone cap's outward
          // side does: its curvature towards that side is 2 ( w . n + rad ) / ( rad^2 - |w|^2 ), with
          // w = q - pcen and n square to d in the plane towards pend, and it is that straight line where d
          // runs along it.  It touches the sphere at ptan, and sphere_arc runs on from there to the pole,
          // leaving along the arc.  With ang zero and dref square to the end it is a tangent ogive.  oarc1
          // and oarc2 are the circular arc's inner control points and parc1 and parc2 the sphere arc's, which
          // sphere_split parts from the circular arc.  With no radius the arc runs through the tip and is
          // split in two halfway along, at ptan.  From a point inside the sphere the side runs straight out to
          // it and on as from the point where it reaches it, which a side from outside the sphere becomes as
          // the point reaches the sphere.
          void ogive_side( const point_type &q, const point_type &qt, const point_type &dref, const data_type &ang,
                           const point_type &pend, const point_type &pcen, const data_type &rad,
                           const point_type &pdir, point_type &oarc1, point_type &oarc2, point_type &ptan,
                           point_type &parc1, point_type &parc2 ) const
          {
            tolerance_type tol;

            point_type ptip = pcen + pdir * rad;
            point_type w = q - pcen;

            point_type ctan, carc1, carc2;
            cone_side( q, pend, pcen, rad, pdir, ctan, carc1, carc2 );

            // The side leaves q along d, its direction from ogive_direction.
            point_type d;
            data_type wlen = w.norm();
            if ( tol.approximately_equal( wlen, 0.0 ) || !ogive_direction( q, qt, dref, ang, pend, d ) )
            {
              cone_ogive_side( q, ctan, carc1, carc2, oarc1, oarc2, ptan, parc1, parc2 );
              return;
            }

            // From a point inside the sphere the side runs straight out to it, at qs, and on from there as it
            // would from a point of the sphere.
            bool inside = wlen <= rad;
            point_type ws = w;
            if ( inside )
            {
              ws = w * ( rad / wlen );
            }
            point_type qs = pcen + ws;

            // The normal of the arc's plane, through qs along d and through pcen, or through pend where d runs
            // through pcen.
            point_type mnorm = d.cross( ws );
            if ( tol.approximately_equal( mnorm.norm(), 0.0 ) )
            {
              mnorm = d.cross( pend - q );
            }
            if ( tol.approximately_equal( mnorm.norm(), 0.0 ) )
            {
              cone_ogive_side( q, ctan, carc1, carc2, oarc1, oarc2, ptan, parc1, parc2 );
              return;
            }
            mnorm = mnorm / mnorm.norm();

            // Square to d in the arc's plane, towards pend.
            point_type n = mnorm.cross( d );
            if ( n.dot( pend - q ) < 0 )
            {
              n = -n;
            }

            // The arc's chord from qs to where it touches the sphere, - ( ws + n rad ) / ( 1 - kappa rad ) for
            // an arc of curvature kappa towards n.  On the sphere the arc has shrunk to qs, and the chord's
            // direction is its limit there.
            point_type cw = -( ws + n * rad );
            data_type cwlen = cw.norm();
            if ( tol.approximately_equal( cwlen, 0.0 ) )
            {
              cone_ogive_side( q, ctan, carc1, carc2, oarc1, oarc2, ptan, parc1, parc2 );
              return;
            }
            point_type cu = cw / cwlen;
            data_type clen = 0;
            if ( !inside )
            {
              data_type kappa = 2.0 * ( w.dot( n ) + rad ) / ( rad * rad - w.squaredNorm() );
              data_type den = 1.0 - kappa * rad;
              clen = cwlen / std::abs( den );
              if ( den < 0 )
              {
                cu = -cu;
              }
            }
            ptan = qs + cu * clen;
            point_type dt = cu * ( 2.0 * cu.dot( d ) ) - d;
            data_type psi = 2.0 * unit_angle( d, cu );

            // The arc's last leg is never shorter than a thousandth of the radius, so it always gives the way
            // the side reaches the sphere.
            data_type lmin = 1e-3 * rad;

            if ( inside )
            {
              // Out to the sphere, reaching it along dt.
              data_type slen = ( qs - q ).norm();
              oarc1 = ( q * 2.0 + qs ) / 3.0;
              oarc2 = qs - dt * std::max( slen / 3.0, lmin );
              sphere_arc( ptan, dt, pcen, rad, pdir, parc1, parc2 );
              return;
            }

            if ( tol.approximately_equal( rad, 0.0 ) )
            {
              // The arc through the tip, in two halves meeting halfway along it.
              data_type turn = 1.0;
              if ( mnorm.dot( d.cross( cu ) ) < 0 )
              {
                turn = -1.0;
              }
              point_type dm = rotate_about( d, mnorm, turn * 0.5 * psi );
              data_type hlen = half_chord( psi, clen );
              point_type pm = q + rotate_about( d, mnorm, turn * 0.25 * psi ) * hlen;
              data_type k = arc_handle( 0.5 * psi, hlen );
              oarc1 = q + d * k;
              oarc2 = pm - dm * k;
              ptan = pm;
              parc1 = pm + dm * k;
              parc2 = ptip - dt * k;
              return;
            }

            data_type k = arc_handle( psi, clen );
            oarc1 = q + d * k;
            oarc2 = ptan - dt * std::max( k, lmin );

            sphere_split( q, oarc1, oarc2, ptan, dt, pcen, rad, pdir );
            sphere_arc( ptan, dt, pcen, rad, pdir, parc1, parc2 );
          }

          // The cone cap's side from q, with tangent point ctan and sphere arc control points carc1 and
          // carc2, as an ogive cap's side whose arc is the straight line from q to ctan.
          void cone_ogive_side( const point_type &q, const point_type &ctan, const point_type &carc1,
                                const point_type &carc2, point_type &oarc1, point_type &oarc2, point_type &ptan,
                                point_type &parc1, point_type &parc2 ) const
          {
            ptan = ctan;
            parc1 = carc1;
            parc2 = carc2;
            oarc1 = ( q * 2.0 + ptan ) / 3.0;
            oarc2 = ( q + ptan * 2.0 ) / 3.0;
          }

          // The chord of half of an arc of a circle through angle psi with a chord of length clen.
          data_type half_chord( const data_type &psi, const data_type &clen ) const
          {
            data_type s = std::sin( 0.5 * psi );
            if ( std::abs( s ) < 1e-8 )
            {
              return 0.5 * clen;
            }
            return clen * std::sin( 0.25 * psi ) / s;
          }

          // v turned by angle a about the unit vector k.
          point_type rotate_about( const point_type &v, const point_type &k, const data_type &a ) const
          {
            return v * std::cos( a ) + k.cross( v ) * std::sin( a ) + k * ( k.dot( v ) * ( 1.0 - std::cos( a ) ) );
          }

          // The control points of the side of an ogive cap from the middle of half, an end's half with t0 at
          // its start, whose surface derivative across the end is dhalf: the point of the end, the ogive arc's
          // inner control points, the tangent point, the sphere arc's inner control points and the pole.
          template < typename curve__ >
          void ogive_reference_side( const curve__ &half, const curve__ &dhalf, const data_type &tmid,
                                     const point_type &pnorm, const point_type &pend, const point_type &pcen,
                                     const data_type &rad, const point_type &pdir,
                                     std::vector < point_type > &cp ) const
          {
            tolerance_type tol;

            cp.resize( 7 );
            cp[0] = half.f( tmid );

            point_type dref = pnorm;
            if ( ogive_surfref )
            {
              point_type d = dhalf.f( tmid );
              data_type dlen = d.norm();
              if ( !tol.approximately_equal( dlen, 0.0 ) )
              {
                dref = d / dlen;
              }
            }

            ogive_side( cp[0], half.fp( tmid ), dref, ogive_angle, pend, pcen, rad, pdir, cp[1], cp[2], cp[3],
                        cp[4], cp[5] );
            cp[6] = pcen + pdir * rad;
          }

          // The side of an ogive cap from point q of the end, blended from the top and bottom sides top and
          // bottom, each the seven control points from ogive_reference_side.  Their control points are framed
          // in pairs: m their middle, a half their difference, w the half length of a, and across the unit
          // vector square to the plane of the pair and pcen.  q's place in the end's frame is rho times the
          // unit vector at angle phi from a, and each of the side's control points is m plus its frame's unit
          // vector at phi times rho, eased from q's to one at the tangent point, or with no radius kept at
          // q's.  The tangent point's frame is a circle on the sphere, so any end's shape becomes the sphere's
          // circle there.  On the sphere the side is its arc from the tangent point to the pole, leaving along
          // the ogive arc's direction there.  The first ogive control point keeps q's rho and moves along q's
          // unit vector so the side leaves q in the plane of qt, the end's tangent there, and dq, the
          // direction the side in its own plane leaves q along.  Returns false where the frame is degenerate.
          bool ogive_blend_side( const point_type &q, const point_type &qt, const point_type &dq,
                                 const point_type &pcen, const data_type &rad, const point_type &pdir,
                                 const std::vector < point_type > &top,
                                 const std::vector < point_type > &bottom, point_type &oarc1, point_type &oarc2,
                                 point_type &ptan, point_type &parc1, point_type &parc2 ) const
          {
            tolerance_type tol;

            index_type n = top.size();
            std::vector < point_type > m( n ), a( n ), across( n );
            std::vector < data_type > alen( n ), w( n );
            for ( index_type k = 0; k < n; k++ )
            {
              m[k] = ( top[k] + bottom[k] ) * 0.5;
              a[k] = ( top[k] - bottom[k] ) * 0.5;
              alen[k] = a[k].norm();
            }
            if ( tol.approximately_equal( alen[0], 0.0 ) )
            {
              return false;
            }

            for ( index_type k = 0; k < n; k++ )
            {
              point_type c = ( top[k] - pcen ).cross( bottom[k] - pcen );
              if ( tol.approximately_equal( c.norm(), 0.0 ) )
              {
                c = a[k].cross( pdir );
              }
              data_type clen = c.norm();
              if ( tol.approximately_equal( clen, 0.0 ) )
              {
                if ( k == 0 )
                {
                  return false;
                }
                across[k] = across[k - 1];
              }
              else
              {
                across[k] = c / clen;
                if ( k > 0 && across[k].dot( across[k - 1] ) < 0 )
                {
                  across[k] = -across[k];
                }
              }
            }

            // Half lengths at the control points: the end values, and the slopes there for the inner ones.
            index_type ends[3] = { 0, 3, 6 };
            for ( index_type e = 0; e < 3; e++ )
            {
              w[ ends[e] ] = alen[ ends[e] ];
            }
            index_type inner[4][2] = { { 1, 0 }, { 2, 3 }, { 4, 3 }, { 5, 6 } };
            for ( index_type e = 0; e < 4; e++ )
            {
              index_type k = inner[e][0];
              index_type j = inner[e][1];
              if ( tol.approximately_equal( alen[j], 0.0 ) )
              {
                w[k] = ( a[k] - a[j] ).norm();
              }
              else
              {
                w[k] = alen[j] + a[j].dot( a[k] - a[j] ) / alen[j];
              }
            }

            // q's place in the end's frame.
            point_type r = q - m[0];
            data_type x = r.dot( a[0] ) / ( alen[0] * alen[0] );
            data_type y = r.dot( across[0] ) / alen[0];
            data_type rho = std::sqrt( x * x + y * y );
            data_type cphi = 1.0;
            data_type sphi = 0.0;
            if ( !tol.approximately_equal( rho, 0.0 ) )
            {
              cphi = x / rho;
              sphi = y / rho;
            }

            std::vector < point_type > u( n );
            for ( index_type k = 0; k < n; k++ )
            {
              u[k] = a[k] * cphi + across[k] * ( w[k] * sphi );
            }

            data_type rend = 1.0;
            if ( tol.approximately_equal( rad, 0.0 ) )
            {
              rend = rho;
            }

            oarc1 = q + ( m[1] - m[0] ) + ( u[1] - u[0] ) * rho;
            oarc2 = m[2] + u[2] * rend;
            ptan = m[3] + u[3] * rend;

            // Where one of the legs from the top and bottom sides' last ogive control points to their tangent
            // points is far shorter than the other, framing them as they are lets the longer one set the
            // side's direction on to the sphere.  There the side's leg is framed from the two at a common
            // length instead, and runs along the sphere's tangent plane at its length eased between theirs.
            data_type ltop = ( top[3] - top[2] ).norm();
            data_type lbot = ( bottom[3] - bottom[2] ).norm();
            data_type lmax = std::max( ltop, lbot );
            data_type keep = 1.0;
            if ( !tol.approximately_equal( lmax, 0.0 ) )
            {
              keep = ( std::min( ltop, lbot ) / lmax - 0.05 ) / 0.2;
            }
            if ( keep < 1.0 )
            {
              keep = std::max( keep, static_cast<data_type>( 0 ) );
              data_type lmid = 0.5 * ( ltop + lbot );
              point_type t2 = top[3];
              point_type b2 = bottom[3];
              if ( !tol.approximately_equal( ltop, 0.0 ) )
              {
                t2 = top[3] - ( top[3] - top[2] ) * ( lmid / ltop );
              }
              if ( !tol.approximately_equal( lbot, 0.0 ) )
              {
                b2 = bottom[3] - ( bottom[3] - bottom[2] ) * ( lmid / lbot );
              }
              point_type a2 = ( t2 - b2 ) * 0.5;
              point_type u2 = a2 * cphi + across[2] * ( a2.norm() * sphi );
              point_type leg = ptan - ( ( t2 + b2 ) * 0.5 + u2 * rend );
              if ( !tol.approximately_equal( rad, 0.0 ) && on_sphere( ptan, pcen, rad ) )
              {
                point_type nsph = ( ptan - pcen ) / rad;
                leg = leg - nsph * nsph.dot( leg );
              }
              data_type leglen = leg.norm();
              if ( !tol.approximately_equal( leglen, 0.0 ) )
              {
                point_type o2 = ptan - leg * ( ( lmid + 0.5 * ( ltop - lbot ) * cphi ) / leglen );
                oarc2 = oarc2 * keep + o2 * ( 1.0 - keep );
              }
            }
            parc1 = m[4] + u[4] * rend;
            parc2 = m[5] + u[5] * rend;

            // On the sphere the side is its arc from the tangent point to the pole, leaving along the ogive arc's
            // direction there.
            if ( !tol.approximately_equal( rad, 0.0 ) )
            {
              point_type lin = ptan - oarc2;
              data_type linlen = lin.norm();
              if ( !tol.approximately_equal( linlen, 0.0 ) )
              {
                sphere_arc( ptan, lin / linlen, pcen, rad, pdir, parc1, parc2 );
              }
            }

            // The side leaves q in the plane of qt and dq.
            if ( tol.approximately_equal( qt.norm(), 0.0 ) )
            {
              return true;
            }
            point_type nu = dq.cross( qt / qt.norm() );

            data_type den = nu.dot( u[0] );
            if ( std::abs( den ) > 1e-3 * u[0].norm() )
            {
              oarc1 = oarc1 - u[0] * ( nu.dot( oarc1 - q ) / den );
            }
            return true;
          }

          // Whether point p lies on the sphere of center pcen and radius rad.
          bool on_sphere( const point_type &p, const point_type &pcen, const data_type &rad ) const
          {
            return std::abs( ( p - pcen ).norm() - rad ) <= 1e-6 * rad;
          }

          // The unit tangents t0 and t1 of the end at the start and end of segment i of half, one of the end's
          // halves and other the other, each the mean of the tangents of the two segments that meet there,
          // so the segments on either side of a corner of the end turn their sides about the same line.
          template < typename pcurve__ >
          void end_tangents( const pcurve__ &half, const pcurve__ &other, const index_type &i, point_type &t0,
                             point_type &t1 ) const
          {
            index_type n = half.number_segments();
            typename pcurve__::curve_type seg, prev, next;
            half.get( seg, i );

            point_type a = seg.fp( 0 );
            point_type b;
            if ( i > 0 )
            {
              half.get( prev, i - 1 );
              b = prev.fp( 1 );
            }
            else
            {
              other.get( prev, 0 );
              b = prev.fp( 0 );
            }
            t0 = mean_line( a, b );

            a = seg.fp( 1 );
            if ( i < n - 1 )
            {
              half.get( next, i + 1 );
              b = next.fp( 0 );
            }
            else
            {
              other.get( next, n - 1 );
              b = next.fp( 1 );
            }
            t1 = mean_line( a, b );
          }

          // The unit vector along the mean of the lines along a and b, whichever way each points.
          point_type mean_line( const point_type &a, const point_type &b ) const
          {
            tolerance_type tol;

            point_type au = a;
            point_type bu = b;
            if ( !tol.approximately_equal( au.norm(), 0.0 ) )
            {
              au = au / au.norm();
            }
            if ( !tol.approximately_equal( bu.norm(), 0.0 ) )
            {
              bu = bu / bu.norm();
            }
            if ( au.dot( bu ) < 0 )
            {
              bu = -bu;
            }
            point_type m = au + bu;
            if ( tol.approximately_equal( m.norm(), 0.0 ) )
            {
              return au;
            }
            return m / m.norm();
          }

          // The sides of an ogive cap from segment c of the end, whose middle is pend, with cd the surface's
          // derivative across the end along the segment, and qt0 and qt1 the end's tangents to take at the
          // segment's ends.  The sides gather at the pole pcen + pdir * rad.
          // Their angle is measured from the surface's direction across the end where ogive_surfref is set,
          // or else from pnorm, the end's normal.  At points along the segment, which lie on the end,
          // ogive_side finds each side on its own, and ogive_blend_side blends it from the top and bottom
          // sides otop and obot, leaving along the side's own direction; the blended side is taken where the
          // blend's frame is sound, and the side on its own elsewhere.  fit_side_rows fits the rows to them:
          // the ogive arc's two inner control points, the tangent point and the sphere arc's two inner control
          // points.
          template < typename curve__ >
          void ogive_sides( const curve__ &c, const curve__ &cd, const point_type &qt0, const point_type &qt1,
                            const point_type &pnorm, const point_type &pend,
                            const point_type &pcen, const data_type &rad, const point_type &pdir,
                            const std::vector < point_type > &otop, const std::vector < point_type > &obot,
                            std::vector < std::vector < point_type > > &rows ) const
          {
            tolerance_type tol;

            index_type nsamp = 8 * c.degree() + 1;

            std::vector < data_type > ts( nsamp );
            std::vector < std::vector < point_type > > vals( 5, std::vector < point_type >( nsamp ) );
            for ( index_type i = 0; i < nsamp; i++ )
            {
              ts[i] = static_cast<data_type>( i ) / ( nsamp - 1 );
              point_type q = c.f( ts[i] );
              point_type qt = c.fp( ts[i] );
              if ( i == 0 )
              {
                qt = qt0;
              }
              if ( i == nsamp - 1 )
              {
                qt = qt1;
              }

              point_type dref = pnorm;
              if ( ogive_surfref )
              {
                point_type d = cd.f( ts[i] );
                data_type dlen = d.norm();
                if ( !tol.approximately_equal( dlen, 0.0 ) )
                {
                  dref = d / dlen;
                }
              }

              ogive_side( q, qt, dref, ogive_angle, pend, pcen, rad, pdir, vals[0][i], vals[1][i],
                          vals[2][i], vals[3][i], vals[4][i] );

              point_type b[5];
              point_type dq = vals[0][i] - q;
              if ( !tol.approximately_equal( dq.norm(), 0.0 ) )
              {
                dq = dq / dq.norm();
              }
              if ( ogive_blend_side( q, qt, dq, pcen, rad, pdir, otop, obot, b[0], b[1], b[2], b[3], b[4] ) )
              {
                for ( index_type k = 0; k < 5; k++ )
                {
                  vals[k][i] = b[k];
                }
              }
            }

            fit_side_rows( c, pend, pdir, ts, vals, true, rows );
          }

          // The side of a truncated ogive cap from point q of the end: the circular arc that leaves q along d,
          // from ogive_direction, and runs through q's point of the face, qf.  oarc1 and oarc2 are its inner
          // control points.  Where d runs along the chord the arc is the straight line to qf.
          void ogive_truncated_side( const point_type &q, const point_type &qt, const point_type &dref,
                                     const point_type &pend, const point_type &qf, point_type &oarc1,
                                     point_type &oarc2 ) const
          {
            tolerance_type tol;

            point_type c = qf - q;
            data_type clen = c.norm();
            oarc1 = ( q * 2.0 + qf ) / 3.0;
            oarc2 = ( q + qf * 2.0 ) / 3.0;
            if ( tol.approximately_equal( clen, 0.0 ) )
            {
              return;
            }
            point_type cu = c / clen;

            point_type d;
            if ( !ogive_direction( q, qt, dref, ogive_angle, pend, d ) )
            {
              return;
            }

            // The arc turns through psi, and reaches qf along d reflected in the chord.  Its handles are
            // those of no more than three quarters of a turn.
            data_type psi = 2.0 * unit_angle( d, cu );
            point_type dt = cu * ( 2.0 * cu.dot( d ) ) - d;
            data_type h = arc_handle( std::min( psi, 1.5 * eli::constants::math<data_type>::pi() ), clen );
            oarc1 = q + d * h;
            oarc2 = qf - dt * h;
          }

          // The sides of a truncated ogive cap from segment c of the end, whose middle is pend, to the face,
          // the end scaled about pend by tface and moved to pface.  cd, qt0, qt1 and pnorm are as for
          // ogive_sides.  At points along the segment, which lie on the end, ogive_truncated_side finds each
          // side, and fit_side_rows fits the rows to them: the arc's two inner control points.
          template < typename curve__ >
          void ogive_truncated_sides( const curve__ &c, const curve__ &cd, const point_type &qt0,
                                      const point_type &qt1, const point_type &pnorm, const point_type &pend,
                                      const point_type &pface, const data_type &tface,
                                      std::vector < std::vector < point_type > > &rows ) const
          {
            tolerance_type tol;

            // The sides run about the line from pend to pface.
            point_type pdir = pface - pend;
            if ( tol.approximately_equal( pdir.norm(), 0.0 ) )
            {
              pdir = pnorm;
            }
            else
            {
              pdir = pdir / pdir.norm();
            }

            index_type nsamp = 8 * c.degree() + 1;

            std::vector < data_type > ts( nsamp );
            std::vector < std::vector < point_type > > vals( 2, std::vector < point_type >( nsamp ) );
            for ( index_type i = 0; i < nsamp; i++ )
            {
              ts[i] = static_cast<data_type>( i ) / ( nsamp - 1 );
              point_type q = c.f( ts[i] );
              point_type qt = c.fp( ts[i] );
              if ( i == 0 )
              {
                qt = qt0;
              }
              if ( i == nsamp - 1 )
              {
                qt = qt1;
              }

              point_type dref = pnorm;
              if ( ogive_surfref )
              {
                point_type d = cd.f( ts[i] );
                data_type dlen = d.norm();
                if ( !tol.approximately_equal( dlen, 0.0 ) )
                {
                  dref = d / dlen;
                }
              }

              point_type qf = pface + ( q - pend ) * tface;
              ogive_truncated_side( q, qt, dref, pend, qf, vals[0][i], vals[1][i] );
            }

            fit_side_rows( c, pend, pdir, ts, vals, true, rows );
          }

          // Unit normal of the end's best plane, by Newell's method over points along it, on axis's side.
          point_type end_normal( const point_type &axis ) const
          {
            tolerance_type tol;

            index_type n = 256;
            data_type t0 = edge.get_t0();
            data_type dt = ( edge.get_tmax() - t0 ) / n;

            point_type nrm;
            nrm << 0, 0, 0;
            for ( index_type k = 0; k < n; k++ )
            {
              point_type a = edge.f( t0 + k * dt );
              point_type b = edge.f( t0 + ( k + 1 ) * dt );
              nrm = nrm + a.cross( b );
            }

            data_type len = nrm.norm();
            if ( tol.approximately_equal( len, 0.0 ) )
            {
              return axis;
            }
            nrm = nrm / len;
            if ( nrm.dot( axis ) < 0 )
            {
              nrm = -nrm;
            }
            return nrm;
          }

          // The control point of a profile cap's row at profile control point c, ( s, f ), for the point p of
          // the end: the end's point p scaled about its middle pend by f and moved towards the tip ptip, s of
          // the way from it.
          point_type profile_row( const point_type &c, const point_type &p, const point_type &pend,
                                  const point_type &ptip ) const
          {
            return ptip + ( pend - ptip ) * c( 0 ) + ( p - pend ) * c( 1 );
          }

          // A profile cap's profile from its point, t zero, to the end, t one: the point ( s, f ) and its
          // derivative in t, s the distance from the point as a fraction of the cap's length and f the radius as
          // a fraction of the end's.  A power law's radius is s raised to exponent, Von Karman's the Haack series
          // with C haack, and Sears-Haack's the front half of the Sears-Haack body.  Near the point each radius
          // rises as s to a power p, the exponent or three quarters.  Below one, s is t raised to one over p, so
          // f rises with t: the profile leaves the point square to the axis with a finite derivative.  From one
          // up, t is s, and the profile leaves the point along the axis.  Either way t runs along the flatter of
          // s and f, and a power law of exponent p and one of 1 / p, the same curve with s and f swapped, are
          // fitted alike.
          class profile_function
          {
            public:
              index_type shape;
              data_type exponent;
              data_type haack;

              void operator()( const data_type &t, point_type &p, point_type &m ) const
              {
                data_type x = std::min( std::max( t, static_cast<data_type>( 0 ) ), static_cast<data_type>( 1 ) );

                if ( shape == POWER_LAW )
                {
                  if ( exponent >= 1 )
                  {
                    data_type df = 0;
                    if ( x > 0 )
                    {
                      df = exponent * std::pow( x, exponent - 1.0 );
                    }
                    else if ( exponent == 1 )
                    {
                      df = 1;
                    }
                    p << x, std::pow( x, exponent ), 0;
                    m << 1, df, 0;
                    return;
                  }

                  data_type ds = 0;
                  if ( x > 0 )
                  {
                    ds = std::pow( x, 1.0 / exponent - 1.0 ) / exponent;
                  }
                  p << std::pow( x, 1.0 / exponent ), x, 0;
                  m << ds, 1, 0;
                  return;
                }

                data_type s = std::pow( x, static_cast<data_type>( 4.0 / 3.0 ) );
                data_type ds = 4.0 / 3.0 * std::cbrt( x );

                if ( shape == SEARS_HAACK )
                {
                  // ( s ( 2 - s ) ) to the three quarters is t ( 2 - s ) to the three quarters
                  data_type g = 2.0 - s;
                  data_type gq = std::pow( g, static_cast<data_type>( 0.75 ) );
                  p << s, x * gq, 0;
                  m << ds, gq - 0.75 * x * ds * gq / g, 0;
                  return;
                }

                // The Haack series at th, where s is the half of one less the cosine of th
                data_type pi = eli::constants::math<data_type>::pi();
                data_type th = 2.0 * std::asin( std::sqrt( s ) );
                data_type a;
                if ( th < 1e-2 )
                {
                  data_type th2 = th * th;
                  a = th * th2 * ( 2.0 / 3.0 - th2 * ( 2.0 / 15.0 - th2 * 4.0 / 315.0 ) );
                }
                else
                {
                  a = th - 0.5 * std::sin( 2.0 * th );
                }
                data_type sth = std::sin( th );
                data_type f = std::sqrt( std::max( a + haack * sth * sth * sth, static_cast<data_type>( 0 ) ) / pi );

                data_type df = std::sqrt( 8.0 * ( 2.0 / 3.0 + haack ) / pi );
                if ( x > 0 && f > 0 )
                {
                  df = sth * ( 2.0 + 3.0 * haack * std::cos( th ) ) / ( pi * f ) * ds;
                }
                p << s, f, 0;
                m << ds, df, 0;
              }
          };

          // The profile fun a profile cap's rows follow, as set_profile_conditions sets it.
          void profile_setup( profile_function &fun ) const
          {
            fun.shape = profile;
            fun.exponent = profile_exp;
            fun.haack = profile_haack;
          }

          // A profile cap's profile pc, along fun from its parameter tstart, zero at its point, to the end:
          // cubics matching the profile's points and derivatives at their ends, halved until each misses the
          // profile at its middle by less than a ten-thousandth, at most eight times.
          void profile_fit( piecewise_curve_type &pc, const data_type &tstart, const profile_function &fun ) const
          {
            eli::geom::curve::piecewise_binary_cubic_creator_abstract< profile_function, data_type, dim__,
                                                                       tolerance_type > pbcc;
            pbcc.setup( fun, tstart, 1, 1e-4, 1, 8 );
            pbcc.create( pc );
          }

          // The parameter along fun where its radius is the fraction f of the end's, by bisection, the radius
          // rising with the parameter.  Kept short of the end, where the profile would have no length left.
          data_type profile_start( const data_type &f, const profile_function &fun ) const
          {
            data_type tmax = 0.999;
            data_type lo = 0;
            data_type hi = tmax;
            point_type p, m;
            fun( hi, p, m );
            if ( p( 1 ) <= f )
            {
              return tmax;
            }
            for ( index_type k = 0; k < 60; k++ )
            {
              data_type mid = 0.5 * ( lo + hi );
              fun( mid, p, m );
              if ( p( 1 ) < f )
              {
                lo = mid;
              }
              else
              {
                hi = mid;
              }
            }
            return 0.5 * ( lo + hi );
          }

          // Rows of control points across segment c of the end, whose middle is pend, fitted to their exact
          // values vals at the segment's parameters ts.  Each exact value lies in the plane of its point of
          // the end and the line through pend along pdir, so it is pend plus a distance along pdir plus a
          // multiple of that point's offset across the line.  Each row is fitted with those two numbers per
          // control point of the segment, applied to the control point's own offset across the line: exact
          // at the segment's ends, least squares between.  The rows then go round the line as the end does,
          // and a circular end centred on it is reproduced exactly.  Where across is set the values may leave
          // those planes, and each control point takes a third number, a multiple of pdir crossed with its
          // offset.
          template < typename curve__ >
          void fit_side_rows( const curve__ &c, const point_type &pend, const point_type &pdir,
                              const std::vector < data_type > &ts, const std::vector < std::vector < point_type > > &vals,
                              const bool &across, std::vector < std::vector < point_type > > &rows ) const
          {
            typedef Eigen::Matrix< data_type, Eigen::Dynamic, Eigen::Dynamic > matrix_type;

            index_type deg = c.degree();
            index_type nsamp = ts.size();

            // Offsets across the line of the segment's control points, and those turned about the line.
            std::vector < point_type > cpoff( deg + 1 ), cpturn( deg + 1 );
            for ( index_type i = 0; i <= deg; i++ )
            {
              point_type o = c.get_control_point( i ) - pend;
              cpoff[i] = o - pdir * pdir.dot( o );
              cpturn[i] = pdir.cross( cpoff[i] );
            }

            index_type nper = 2;
            if ( across )
            {
              nper = 3;
            }

            rows.resize( vals.size() );
            for ( index_type k = 0; k < vals.size(); k++ )
            {
              std::vector < data_type > a( deg + 1, 0 ), b( deg + 1, 0 ), g( deg + 1, 0 );

              // The ends are on the end: their numbers come straight from the exact values.
              index_type iend[2] = { 0, deg };
              index_type isamp[2] = { 0, nsamp - 1 };
              for ( index_type e = 0; e < 2; e++ )
              {
                point_type x = vals[k][ isamp[e] ] - pend;
                a[ iend[e] ] = pdir.dot( x );
                data_type olen2 = cpoff[ iend[e] ].squaredNorm();
                if ( olen2 > 0 )
                {
                  b[ iend[e] ] = ( x - pdir * a[ iend[e] ] ).dot( cpoff[ iend[e] ] ) / olen2;
                  if ( across )
                  {
                    g[ iend[e] ] = x.dot( cpturn[ iend[e] ] ) / olen2;
                  }
                }
              }

              if ( deg >= 2 )
              {
                index_type ncol = nper * ( deg - 1 );
                matrix_type A( 3 * nsamp, ncol ), rhs( 3 * nsamp, 1 ), x( ncol, 1 );
                A.setZero();
                for ( index_type m = 0; m < nsamp; m++ )
                {
                  std::vector < data_type > bern( deg + 1 );
                  for ( index_type i = 0; i <= deg; i++ )
                  {
                    data_type bc;
                    eli::mutil::dm::binomial_coefficient( bc, deg, i );
                    bern[i] = bc * std::pow( ts[m], static_cast<data_type>( i ) ) * std::pow( 1 - ts[m], static_cast<data_type>( deg - i ) );
                  }

                  point_type r0 = vals[k][m] - pend;
                  for ( index_type e = 0; e < 2; e++ )
                  {
                    index_type i = iend[e];
                    r0 = r0 - ( pdir * a[i] + cpoff[i] * b[i] + cpturn[i] * g[i] ) * bern[i];
                  }
                  for ( index_type i = 1; i < deg; i++ )
                  {
                    for ( index_type d = 0; d < 3; d++ )
                    {
                      A( 3 * m + d, nper * ( i - 1 ) ) = bern[i] * pdir( d );
                      A( 3 * m + d, nper * ( i - 1 ) + 1 ) = bern[i] * cpoff[i]( d );
                      if ( across )
                      {
                        A( 3 * m + d, nper * ( i - 1 ) + 2 ) = bern[i] * cpturn[i]( d );
                      }
                    }
                  }
                  for ( index_type d = 0; d < 3; d++ )
                  {
                    rhs( 3 * m + d, 0 ) = r0( d );
                  }
                }

                eli::mutil::opt::least_squares_uncon( x, A, rhs );
                for ( index_type i = 1; i < deg; i++ )
                {
                  a[i] = x( nper * ( i - 1 ), 0 );
                  b[i] = x( nper * ( i - 1 ) + 1, 0 );
                  if ( across )
                  {
                    g[i] = x( nper * ( i - 1 ) + 2, 0 );
                  }
                }
              }

              rows[k].resize( deg + 1 );
              for ( index_type i = 0; i <= deg; i++ )
              {
                rows[k][i] = pend + pdir * a[i] + cpoff[i] * b[i] + cpturn[i] * g[i];
              }
            }
          }

          void dirty_prep()
          {
            prep_done = false;
            surface_prep_done = false;
            surface_copy.clear();
            edge.clear();
            ndelta.clear();
            edelta.clear();
          }

        private:
          piecewise_surface_type orig_surface;
          piecewise_curve_type orig_edge;
          piecewise_curve_type orig_edge_next;
          bool curves_given;
          mutable piecewise_surface_type surface_copy;
          mutable piecewise_curve_type edge;
          mutable piecewise_curve_type ndelta;
          mutable piecewise_curve_type edelta;
          mutable bool prep_done;
          mutable bool surface_prep_done;

          data_type delta_param;
          edge_cap_identifier edge_to_cap;
          data_type len_factor;
          data_type off_factor;
          data_type str_factor;
          bool sweep_correct;
          point_type pnt_offset;
          index_type cap_type;

          data_type v_hmax;
          data_type hmax;
          data_type h_start;
          data_type h_end;

          bool extendte;
          bool extendle;

          data_type rad_factor;
          data_type ogive_angle;
          bool ogive_surfref;
          index_type profile;
          data_type profile_exp;
          data_type profile_haack;
          bool truncate;
          bool frame_given;
          bool shear_offset;
          point_type frame_axis;
          point_type frame_up;
      };
    }
  }
}
#endif
