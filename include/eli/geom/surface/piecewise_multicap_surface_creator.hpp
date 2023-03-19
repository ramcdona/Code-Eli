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
            ROUND_EXT
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
            : piecewise_creator_base<data__, dim__, tol__>(0, 0), delta_param(1), edge_to_cap(CAP_NONE), prep_done(false)
          {
          }
          piecewise_multicap_surface_creator(const piecewise_surface_type &os, const data_type &dp, edge_cap_identifier ec)
            : piecewise_creator_base<data__, dim__, tol__>(0, 0), orig_surface(os), delta_param(dp), edge_to_cap(ec), prep_done(false)
          {
          }
          piecewise_multicap_surface_creator(const piecewise_multicap_surface_creator<data_type, dim__, tolerance_type> & gs)
            : piecewise_creator_base<data_type, dim__, tolerance_type>(gs), orig_surface(gs.orig_surface),
              delta_param(gs.delta_param), edge_to_cap(gs.edge_to_cap), prep_done(false)
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
            edge_to_cap = ec;
            delta_param = dp;
            len_factor = lf;
            off_factor = of;
            pnt_offset = pnt;
            str_factor = sf;
            sweep_correct = sc;

            if ( orig_surface.number_u_patches() == 0 || orig_surface.number_v_patches() == 0 )
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

          virtual bool create(piecewise_surface_type &ps) const
          {
            typedef typename eli::geom::curve::piecewise<eli::geom::curve::bezier, data_type, dim__, tolerance_type> piecewise_curve_type;
            typedef typename piecewise_curve_type::curve_type curve_type;
            typedef typename piecewise_surface_type::surface_type surface_type;

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
                      s1.set_control_point( porigin + pnt_offset, 1, j );

                      // Build linear lower side curve
                      s2.set_control_point( porigin + pnt_offset, 0, j );
                      s2.set_control_point( pdn, 1, j );
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

                if ( cap_type == ROUND_EXT )
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
            if ( cap_type == ROUND_EXT )
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

                switch( cap_type ){
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
                  {
                    c1.resize( deg );

                    for ( index_type j = 0; j <= deg; j++ )
                    {
                      point_type pup, pdn;

                      pup = clast.get_control_point( j );
                      pdn = cfirst.get_control_point( j );

                      // Build linear upper side curve
                      c1.set_control_point( porigin + pnt_offset, j );
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
                        ptup = pup + tanup * b * len_factor * ksweep;
                        ptdn = pdn + tandn * b * len_factor * ksweep;

                        // Displacement vector at extrapolated tip
                        point_type dispt = ptup - ptdn;

                        // Point on center of circle.
                        psplit = ( ptup + ptdn ) * 0.5 + dispt * off_factor;

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
                        ptup = pup + tanup * circ_chord * len_factor * ksweep * 0.5;
                        ptdn = pdn + tandn * circ_chord * len_factor * ksweep * 0.5;

                        // Displacement vector at extrapolated tip
                        point_type dispt = ptup - ptdn;

                        psplit = ( ptup + ptdn ) * 0.5 + dispt * off_factor;

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
                        ptup = pup + tanup * circ_chord * len_factor * ksweep * 0.5;
                        ptdn = pdn + tandn * circ_chord * len_factor * ksweep * 0.5;

                        // Displacement vector at extrapolated tip
                        point_type dispt = ptup - ptdn;

                        psplit = ( ptup + ptdn ) * 0.5 + dispt * off_factor;

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

          virtual bool prep() const
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

            // extract the curve corresponding to the edge of surface wanting to cap
            switch (edge_to_cap)
            {
              case (CAP_UMIN):
              {
                surface_copy.get_umin_bndy_curve( edge );
                surface_copy.get_umin_ndelta_pcurve( ndelta );
                break;
              }
              case (CAP_UMAX):
              {
                surface_copy.get_umax_bndy_curve( edge );
                surface_copy.get_umax_ndelta_pcurve( ndelta );
                break;
              }
              case (CAP_VMIN):
              {
                surface_copy.get_vmin_bndy_curve( edge );
                surface_copy.get_vmin_ndelta_pcurve( ndelta );
                break;
              }
              case (CAP_VMAX):
              {
                surface_copy.get_vmax_bndy_curve( edge );
                surface_copy.get_vmax_ndelta_pcurve( ndelta );
                break;
              }
              default:
              {
                // must not want any edge capped
                assert(edge_to_cap==CAP_NONE);
                return false;
              }
            }

            prep_done = true;

            return true;
          }

          void dirty_prep()
          {
            prep_done = false;
            surface_copy.clear();
            edge.clear();
            ndelta.clear();
          }

        private:
          piecewise_surface_type orig_surface;
          mutable piecewise_surface_type surface_copy;
          mutable piecewise_curve_type edge;
          mutable piecewise_curve_type ndelta;
          mutable bool prep_done;

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
      };
    }
  }
}
#endif
