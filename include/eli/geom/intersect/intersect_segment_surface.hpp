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

#ifndef ELI_INTERSECT_SEGMENT_SURFACE_HPP
#define ELI_INTERSECT_SEGMENT_SURFACE_HPP

class abouteq;

#include <cmath>
#include <limits>

#include "eli/code_eli.hpp"
#include "eli/geom/intersect/intersect_segment_triangle.hpp"

namespace eli
{
  namespace geom
  {
    namespace intersect
    {
      // ==== Triangle - Line Segment Intersection ====//
      // ==== A - Base Point on Triangle
      // ==== B - Vector for one   Side of Tri
      // ==== C - Vector for other Side of Tri
      // ==== D - Base Point for Line Seg
      // ==== E - Vector for Line Seg
      // ==============================================//

      template<typename surface__>
      void intersect_segment( const surface__ &s,
                              const typename surface__::point_type &p0,
                              const typename surface__::point_type &vec,
                              const typename surface__::bounding_box_type &line_box,
                              std::vector<typename surface__::data_type> &tvec )
      {
        typedef surface__ surface_type;
        typedef typename surface_type::index_type index_type;
        typedef typename surface_type::data_type data_type;
        typedef typename surface_type::point_type point_type;
        typedef typename surface_type::bounding_box_type bounding_box_type;

        bounding_box_type bbox;

        s.get_bounding_box( bbox );

        if ( !bbox.intersect( line_box ) )
        {
          return;
        }

        index_type n = s.degree_u();
        index_type m = s.degree_v();

        //==== Do Tri Seg intersection ====//
        if ( s.test_planar( 1.0e-5 ) )  // Uses a dimensional tolerance in test.
        {
          data_type u, w, t;
          point_type OA1 = s.get_control_point( 0, 0 );
          point_type A1 = s.get_control_point( n, m ) - OA1;
          point_type B1 = s.get_control_point( n, 0 ) - OA1;
          point_type C1 = s.get_control_point( 0, m ) - OA1;

          if ( seg_tri_intersect( OA1, A1, B1, p0, vec, u, w, t ) )
          {
            tvec.push_back( t );
          }
          if ( seg_tri_intersect( OA1, C1, A1, p0, vec, u, w, t ) )
          {
            tvec.push_back( t );
          }
          return;
        }

        surface_type bps0( n, m );
        surface_type bps1( n, m );
        surface_type bps2( n, m );
        surface_type bps3( n, m );

        s.simple_split_uv_quarter( bps0, bps1, bps2, bps3 );

        intersect_segment( bps0, p0, vec, line_box, tvec );
        intersect_segment( bps1, p0, vec, line_box, tvec );
        intersect_segment( bps2, p0, vec, line_box, tvec );
        intersect_segment( bps3, p0, vec, line_box, tvec );
      }

      template<template<typename, unsigned short, typename> class surface__, typename data__, unsigned short dim__, typename tol__ >
      void intersect_segment( std::vector <data__> &tvec,
          const surface::piecewise<surface__, data__, dim__, tol__> &ps,
          const typename surface::piecewise<surface__, data__, dim__, tol__>::point_type &p0,
          const typename surface::piecewise<surface__, data__, dim__, tol__>::point_type &vec,
          const typename surface::piecewise<surface__, data__, dim__, tol__>::bounding_box_type &bbox)
      {
        typedef surface::piecewise<surface__, data__, dim__, tol__> piecewise_type;
        typedef typename piecewise_type::index_type index_type;
        typedef typename piecewise_type::data_type data_type;
        typedef typename piecewise_type::bounding_box_type bounding_box_type;
        typedef typename piecewise_type::tolerance_type tolerance_type;

        typedef typename piecewise_type::keymap_type keymap_type;
        typedef typename keymap_type::const_iterator keyit;

        tvec.clear();

        bounding_box_type line_box;
        line_box.add( p0 );
        line_box.add( p0 + vec );

        if ( !bbox.intersect( line_box ) )
        {
          return;
        }
        index_type count = 0;

        for( keyit uit = ps.ukey.key.begin(); uit != ps.ukey.key.end(); ++uit )
        {
          for( keyit vit = ps.vkey.key.begin(); vit != ps.vkey.key.end(); ++vit )
          {
            index_type uk = uit->second;
            index_type vk = vit->second;

            intersect_segment( ps.patches[uk][vk], p0, vec, line_box, tvec );
            count++;
          }
        }

        // Sort first.
        std::sort( tvec.begin(), tvec.end() );

        // Then check for approximate uniqueness.
        typename std::vector< data_type >::iterator it;
        it = std::unique( tvec.begin(), tvec.end(),
                        [](const data_type& v1, const data_type& v2)
                        {
                          tolerance_type tol;
                          return tol.approximately_equal( v1, v2 );
                        });
        tvec.resize( std::distance( tvec.begin(), it ) );
      }

      template<template<typename, unsigned short, typename> class surface__, typename data__, unsigned short dim__, typename tol__ >
      void intersect_segment( std::vector <data__> &tvec,
          const surface::piecewise<surface__, data__, dim__, tol__> &ps,
          const typename surface::piecewise<surface__, data__, dim__, tol__>::point_type &p0,
          const typename surface::piecewise<surface__, data__, dim__, tol__>::point_type &vec )
      {
        typedef surface::piecewise<surface__, data__, dim__, tol__> piecewise_type;
        typedef typename piecewise_type::bounding_box_type bounding_box_type;

        bounding_box_type bbox;
        ps.get_bounding_box( bbox );
        intersect_segment( tvec, ps, p0, vec, bbox );
      }

      template<template<typename, unsigned short, typename> class surface__, typename data__, unsigned short dim__, typename tol__ >
      bool inside( const surface::piecewise<surface__, data__, dim__, tol__> &ps,
          const typename surface::piecewise<surface__, data__, dim__, tol__>::point_type &p0 )
      {
        typedef surface::piecewise<surface__, data__, dim__, tol__> piecewise_type;
        typedef typename piecewise_type::data_type data_type;
        typedef typename piecewise_type::bounding_box_type bounding_box_type;
        typedef typename piecewise_type::point_type point_type;

        bounding_box_type bbox;
        ps.get_bounding_box( bbox );

        // Build vector longer than longest dimension, pointing mostly in X, but with slight components
        // in Y, Z such that it is not perfectly axis aligned.
        point_type pmin = bbox.get_min();
        point_type pmax = bbox.get_max();
        data_type len = 1.1 * ( pmax - pmin ).norm();
        point_type vec;
        vec << 1.0, 1.2345e-4, 2.3456e-4;
        vec = len * vec;

        std::vector <data_type> tvec;

        intersect_segment( tvec, ps, p0, vec, bbox );

        return tvec.size() % 2 == 1;
      }


    }
  }
}



#endif //ELI_INTERSECT_SEGMENT_TRIANGLE_HPP
