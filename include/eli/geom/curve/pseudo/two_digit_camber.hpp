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

#ifndef eli_geom_curve_pseudo_two_digit_camber_hpp
#define eli_geom_curve_pseudo_two_digit_camber_hpp

#include <string>    // std::string
#include <sstream>   // std::ostringstream, std::istringstream
#include <iomanip>   // std::setw
#include <algorithm> // std::transform

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
        class two_digit_camber : virtual public naca_af<data__> {
        public:
            typedef naca_af<data__> base_class_type;
            typedef typename base_class_type::data_type data_type;
            typedef typename base_class_type::point_type point_type;
            typedef typename base_class_type::coefficient_type coefficient_type;
            typedef typename base_class_type::index_type index_type;
            typedef typename base_class_type::complex_data_type complex_data_type;

          public:

            // valid camber values are greater than or equal to zero and less or equal to 9
            // valid camber location values are zero or greater than or equal to 1 and less than or equal to 9
            //
            // note that if either is zero then they both should be zero
            bool set_camber(const data_type &cam, const data_type &cam_loc)
            {
              if ((cam == 0) || (cam_loc == 0))
              {
                camber = 0;
                camber_loc = 0;

                return true;
              }

              if ((cam<=0) || (cam>=.09))
              {
                return false;
              }
              if ((cam_loc<.1) || (cam_loc>.9))
              {
                return false;
              }

              camber=cam;
              camber_loc=cam_loc;

              return true;
            }
            data_type get_maximum_camber() const {return camber;}
            data_type get_maximum_camber_location() const {return camber_loc;}

          protected:

            void calc_camber(data_type &y, data_type &yp, data_type &ypp, data_type &yppp, const data_type &xi) const
            {
              // check to make sure given valid parametric value
              assert((xi>=0) && (xi<=1));

              data_type zero(0), one(1), two(2);

              // short circuit if no camber
              if (camber==0)
              {
                y=zero;
                yp=zero;
                ypp=zero;
                yppp=zero;

                return;
              }

              data_type p = camber;
              data_type m = camber_loc;

              if (xi<=m)
              {
                data_type pm2(p/(m*m));

                y=pm2*(xi*(two*m-xi));
                yp=two*pm2*(m-xi);
                ypp=-two*pm2;
                yppp=zero;
              }
              else
              {
                data_type p1m2(p/(one+m*(m-two)));

                y=p1m2*(one-two*m+xi*(two*m-xi));
                yp=two*p1m2*(m-xi);
                ypp=-two*p1m2;
                yppp=zero;
              }
            }


          protected:
            data_type camber;       // maximum camber [0,.09]
            data_type camber_loc;   // chord-wise location of maximum camber [0,.9]
        };
      }
    }
  }
}
#endif
