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

#ifndef eli_geom_curve_pseudo_four_digit_hpp
#define eli_geom_curve_pseudo_four_digit_hpp

#include <string>    // std::string
#include <sstream>   // std::ostringstream, std::istringstream
#include <iomanip>   // std::setw
#include <algorithm> // std::transform

#include "eli/code_eli.hpp"
#include "eli/geom/curve/pseudo/two_digit_camber.hpp"
#include "eli/geom/curve/pseudo/four_digit_thickness.hpp"

namespace eli
{
  namespace geom
  {
    namespace curve
    {
      namespace pseudo
      {
        template<typename data__>
        class four_digit : virtual public four_digit_thickness<data__>, public two_digit_camber<data__> {
        public:
            typedef naca_af<data__> base_class_type;
            typedef typename base_class_type::data_type data_type;
            typedef typename base_class_type::point_type point_type;
            typedef typename base_class_type::coefficient_type coefficient_type;
            typedef typename base_class_type::index_type index_type;
            typedef typename base_class_type::complex_data_type complex_data_type;

          public:
            four_digit()
            {
              this->thickness = 0.10;
              this->sharp_te = false;
              this->camber = 0.0;
              this->camber_loc = 0.0;

              this->recalc_thickness_coefficients();
            }

            four_digit(const four_digit<data_type> &fs)
            {
              this->thickness = fs.thickness;
              this->sharp_te = fs.sharp_te;
              this->camber = fs.camber;
              this->camber_loc = fs.camber_loc;

              this->recalc_thickness_coefficients();
            }

          protected:

          private:

        };
      }
    }
  }
}
#endif
