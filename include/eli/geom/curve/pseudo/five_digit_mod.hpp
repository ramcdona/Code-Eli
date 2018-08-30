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

#ifndef eli_geom_curve_pseudo_five_digit_mod_hpp
#define eli_geom_curve_pseudo_five_digit_mod_hpp

#include <string>    // std::string
#include <sstream>   // std::ostringstream, std::istringstream
#include <iomanip>   // std::setw
#include <algorithm> // std::transform
#include <complex>

#include "eli/code_eli.hpp"
#include "eli/geom/curve/pseudo/three_digit_camber.hpp"
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
        class five_digit_mod : virtual public four_digit_mod_thickness<data__>, virtual public three_digit_camber<data__> {
        public:
            typedef naca_af<data__> base_class_type;
            typedef typename base_class_type::data_type data_type;
            typedef typename base_class_type::point_type point_type;
            typedef typename base_class_type::coefficient_type coefficient_type;
            typedef typename base_class_type::index_type index_type;
            typedef typename base_class_type::complex_data_type complex_data_type;

          public:
            five_digit_mod()
            {
              this->thickness = 0.12;
              this->sharp_te = false;
              this->camber_loc = 0.15;
              this->thickness_loc = 0.3;
              this->le_radius_indx = 6.0;

              this->recalc_thickness_coefficients();
              this->recalc_camber_coefficients();
            }

            five_digit_mod(const five_digit_mod <data_type> &fs)
            {
              this->thickness = fs.thickness;
              this->sharp_te = fs.sharp_te;
              this->camber_loc = fs.camber_loc;
              this->thickness_loc = fs.thickness_loc;
              this->le_radius_indx = fs.le_radius_indx;

              this->recalc_thickness_coefficients();
              this->recalc_camber_coefficients();
            }

            five_digit_mod( const data_type &t, const data_type &c, const data_type &xc, const data_type &lei, const data_type &t_loc, bool fl )
            {
              set( t, c, xc, lei, t_loc, fl );
            }

            bool set( const data_type &t, const data_type &c, const data_type &xc, const data_type &lei, const data_type &t_loc, bool fl )
            {
                if ( (t < 0) || (t > 1) ) {
                    return false;
                }
                this->thickness = t;

                this->cli = c;

                this->camber_loc = xc;

                this->sharp_te = fl;

                if ((lei<0) || (lei>9))
                {
                  return false;
                }
                this->le_radius_indx=lei;

                if ((t_loc<.2) || (t_loc>.6))
                {
                  return false;
                }
                this->thickness_loc=t_loc;


                this->recalc_thickness_coefficients();
                this->recalc_camber_coefficients();

                return true;
            }

          protected:

          private:
        };
      }
    }
  }
}
#endif
