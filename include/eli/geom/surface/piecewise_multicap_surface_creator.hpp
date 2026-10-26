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
            frame_given = false;

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
            frame_given = false;

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

          // Whether a cap of type typ closes on a single point.
          static bool closes_on_point( const index_type &typ )
          {
            switch ( typ )
            {
              case POINT:
              case CONE:
                return true;
              default:
                return false;
            }
          }

          // Whether this cap's sides run from the end onto a spherical tip that touches its point, as a cone
          // cap's do.
          bool sphere_tip() const
          {
            switch ( cap_type )
            {
              case CONE:
                return true;
              default:
                return false;
            }
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
                else
                {
                  cap.set( s1, 0, i );
                  cap.set( s2, 1, i );
                }

              }
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
          bool frame_given;
          point_type frame_axis;
          point_type frame_up;
      };
    }
  }
}
#endif
