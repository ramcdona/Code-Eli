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

#ifndef eli_geom_curve_pseudo_one_six_series_hpp
#define eli_geom_curve_pseudo_one_six_series_hpp

#include <string>    // std::string
#include <sstream>   // std::ostringstream, std::istringstream
#include <iomanip>   // std::setw
#include <algorithm> // std::transform
#include <complex>

#include "eli/code_eli.hpp"
#include "eli/geom/curve/pseudo/naca_af.hpp"
#include "eli/geom/curve/pseudo/four_digit_mod_thickness.hpp"

namespace eli
{
  namespace geom
  {
    namespace curve
    {
      namespace pseudo
      {
        template<typename data__>
        class one_six_series : virtual public naca_af<data__>, virtual public four_digit_mod_thickness<data__> {
          public:
            typedef naca_af<data__> base_class_type;
            typedef typename base_class_type::data_type data_type;
            typedef typename base_class_type::point_type point_type;
            typedef typename base_class_type::coefficient_type coefficient_type;
            typedef typename base_class_type::index_type index_type;
            typedef typename base_class_type::complex_data_type complex_data_type;

          public:
            one_six_series() : naca_af<data__>()
            {
              this->thickness = 0.12;
              cli = 0.3;
              this->sharp_te = false;
              this->thickness_loc = 0.5;
              this->le_radius_indx = 4.0;

              this->recalc_thickness_coefficients();
            }

            one_six_series(const one_six_series <data_type> &fs)
            {
              this->thickness = fs.thickness;
              cli = fs.cli;
              this->sharp_te = fs.sharp_te;
              this->thickness_loc = 0.5;
              this->le_radius_indx = 4.0;

              this->recalc_thickness_coefficients();
            }

            one_six_series( const data_type &t, const data_type &c, bool fl )
            {
              set( t, c, fl );
              this->recalc_thickness_coefficients();
            }

            bool set_cli(const data_type &c)
            {
              cli = c;
              return true;
            }

            data_type get_cli() const { return cli; }

            bool set( const data_type &t, const data_type &c, bool fl )
            {
                this->thickness_loc = 0.5;
                this->le_radius_indx = 4.0;

                if ( (t < 0) || (t > 1) ) {
                    return false;
                }
                this->thickness = t;

                cli = c;

                this->sharp_te = fl;

                return true;
            }

          protected:

            void calc_camber(data_type &y, data_type &yp, data_type &ypp, data_type &yppp, const data_type &xi) const
            {
              // check to make sure given valid parametric value
              assert((xi>=0) && (xi<=1));

              data_type zero(0), one(1), two(2);

              // short circuit if no camber
              if (cli==0)
              {
                y=zero;
                yp=zero;
                ypp=zero;
                yppp=zero;

                return;
              }

              data_type k = cli / ( 4.0 * eli::constants::math<data_type>::pi() );

              data_type x, omx;
              if ( xi == 0 )
              {
                x = xi + 1e-6;
                omx = 1.0 - x;
                y = zero;
              }
              else if ( xi == 1 )
              {
                x = xi - 1e-6;
                omx = 1.0 - x;
                y = zero;
              }
              else
              {
                x = xi;
                omx = 1.0 - x;
                y = -k * ( omx * log( omx ) + x * log( x ) );
              }

              yp = k * ( log( omx ) - log( x ) );
              ypp = -k / ( omx * x );
              yppp = ( k - 2.0 * k * x ) / ( omx * omx * x * x );
            }

          private:
            data_type cli;          // Ideal lift coefficient.

        };
      }
    }
  }
}
#endif
