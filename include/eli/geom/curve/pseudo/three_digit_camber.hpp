/*********************************************************************************
* Copyright (c) 2014 David D. Marshall <ddmarsha@calpoly.edu>
*
* All rights reserved. This program and the accompanying materials
* are made available under the terms of the Eclipse Public License v1.0
* which accompanies this distribution, and is available at
* http://www.eclipse.org/legal/epl-v10.html
*
* Contributors:
*    David D. Marshall - initial code and implementation
********************************************************************************/

#ifndef eli_geom_curve_pseudo_three_digit_camber_hpp
#define eli_geom_curve_pseudo_three_digit_camber_hpp

#include <string>    // std::string
#include <sstream>   // std::ostringstream, std::istringstream
#include <iomanip>   // std::setw
#include <algorithm> // std::transform
#include <complex>

#include "eli/code_eli.hpp"
#include "eli/geom/curve/pseudo/naca_af.hpp"

namespace eli
{
  namespace geom
  {
    namespace curve
    {
      namespace pseudo
      {
        template<typename data__>
        class three_digit_camber : virtual public naca_af<data__> {
        public:
            typedef naca_af<data__> base_class_type;
            typedef typename base_class_type::data_type data_type;
            typedef typename base_class_type::point_type point_type;
            typedef typename base_class_type::coefficient_type coefficient_type;
            typedef typename base_class_type::index_type index_type;
            typedef typename base_class_type::complex_data_type complex_data_type;

        public:

            bool set_cli(const data_type &c) {
                cli = c;
                return true;
            }

            data_type get_cli() const { return cli; }

            bool set_camber_loc(const data_type &xc) {
                camber_loc = xc;
                recalc_camber_coefficients();
                return true;
            }

            data_type get_camber_loc() { return camber_loc; }


          protected:

            void recalc_camber_coefficients()
            {
              complex_data_type xc = camber_loc; // Maximum camber location

              complex_data_type xt;  // Transition point of front to rear equation
              complex_data_type i(0,1); // imaginary number.

              data_type pi = eli::constants::math<data_type>::pi();
              complex_data_type xcsq = xc * xc;
              complex_data_type xccu = xc * xcsq;
              complex_data_type xc4 = xcsq * xcsq;

              xt = ((1.0 - i * sqrt(3.0)) * (18.0 * xc - 9.0))/(9.0 * pow( 2.0, 2.0/3.0 ) * pow( 3.0 * xcsq + sqrt( 9.0 * xc4 - 4.0 * xccu ) - 6.0 * xc + 2.0, 1.0/3.0)) -
                   ((1.0 + i * sqrt(3.0)) * pow(3.0 * xcsq + sqrt(9.0 * xc4 - 4.0 * xccu) - 6.0 * xc + 2.0, 1.0/3.0) )/(2.0 * pow(2.0,1.0/3.0))
                   + 1.0;

              xtrans = xt.real();

              if ( xtrans >= 1.0 )
              {
                xtrans = 1.0 - 1e-6;
              }

              data_type xt2 = xtrans * xtrans;
              data_type xt3 = xtrans * xt2;
              data_type xt4 = xt2 * xt2;

              data_type Q = ( 3.0 * xtrans - 7.0 * xt2 + 8.0 * xt3 - 4.0 * xt4 ) / sqrt( xtrans * ( 1.0-xtrans ) ) -
                            1.5 * ( 1.0 - 2.0 * xtrans ) * ( pi / 2.0 - asin( 1.0 - 2.0 * xtrans ) );

              k = 0.3 * 6.0 / Q;
            }

            void calc_camber(data_type &y, data_type &yp, data_type &ypp, data_type &yppp, const data_type &xi) const
            {
              // check to make sure given valid parametric value
              assert((xi>=0) && (xi<=1));

              data_type zero(0), one(1), two(2);

              // short circuit if no camber
              if ( cli==0 || camber_loc == 0 )
              {
                y=zero;
                yp=zero;
                ypp=zero;
                yppp=zero;

                return;
              }

              data_type r = cli / 0.3;

              if (xi<=xtrans)
              {
                data_type xi2 = xi * xi;
                data_type xi3 = xi2 * xi;

                data_type xt2 = xtrans * xtrans;
                data_type xt3 = xtrans * xt2;

                y = r * ( k / 6.0) * ( xi3 - 3.0 * xtrans * xi2 + xt2 * (3.0-xtrans) * xi );
                yp = - r * ( k / 6.0 ) * ( xt3 - 3.0 * xt2 + 6.0 * xtrans * xi - 3.0 * xi2 );
                ypp = r * k * ( xi - xtrans);
                yppp = zero;
              }
              else
              {
                yp = - r * ( k / 6.0 ) * xtrans * xtrans * xtrans;
                y = yp * ( xi - 1.0 );
                ypp = zero;
                yppp = zero;
              }
            }

          protected:
            data_type cli;          // Ideal lift coefficient.
            data_type camber_loc;   // chord-wise location of maximum camber index

            data_type xtrans;
            data_type k;
        };
      }
    }
  }
}
#endif
