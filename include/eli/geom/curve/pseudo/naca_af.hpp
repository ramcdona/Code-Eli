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

#ifndef eli_geom_curve_pseudo_naca_af_hpp
#define eli_geom_curve_pseudo_naca_af_hpp

#include <string>    // std::string
#include <sstream>   // std::ostringstream, std::istringstream
#include <iomanip>   // std::setw
#include <algorithm> // std::transform
#include <complex>

#include "eli/code_eli.hpp"

namespace eli
{
  namespace geom
  {
    namespace curve
    {
      namespace pseudo
      {
        template<typename data__>
        class naca_af {
          public:
            typedef data__ data_type;
            typedef Eigen::Matrix<data_type, 1, 2> point_type;
            typedef Eigen::Matrix<data_type, 5, 1> coefficient_type;
            typedef typename point_type::Index index_type;
            typedef std::complex <data_type> complex_data_type;

          public:

            data_type get_t0() const { return static_cast<data_type>(-1); }

            data_type get_tmax() const { return static_cast<data_type>(1); }

            void set_sharp_trailing_edge(bool fl)
            {
              sharp_te = fl;
              recalc_thickness_coefficients();
            }

            bool sharp_trailing_edge() const { return sharp_te; }

            // Valid values of thickness are greater than 0 and less than 100
            bool set_thickness(const data_type &t)
            {
              if ((t >= 0) && (t <= 1))
              {
                thickness = t;
                return true;
              }
              return false;
            }

            data_type get_thickness() const { return thickness; }

            point_type f(const data_type &xi) const
            {
              point_type x, xp, xpp;

              evaluate(x, xp, xpp, xi);

              return x;
            }

            point_type fp(const data_type &xi) const
            {
              point_type x, xp, xpp;

              evaluate(x, xp, xpp, xi);

              return xp;
            }

            point_type fpp(const data_type &xi) const
            {
              point_type x, xp, xpp;

              evaluate(x, xp, xpp, xi);

              return xpp;
            }

            void evaluate(point_type &x, point_type &xp, point_type &xpp, const data_type &xi) const
            {
              // check to make sure given valid parametric value
              assert((xi>=-1) && (xi<=1));

              data_type xc, yc, ycp, ycpp, ycppp, delta, deltap, deltapp;
              const data_type one(1), two(2), three(3);
              index_type surf_sign;

              // calculate the lower surface
              if (xi<0)
              {
                xc=-xi;
                surf_sign=-1;
              }
              // calculate the upper surface
              else
              {
                xc=xi;
                surf_sign=1;
              }

              // calculate the supporting quantities needed
              calc_camber(yc, ycp, ycpp, ycppp, xc);
              calc_thickness(delta, deltap, deltapp, xc);

              data_type tmp1, tmp13, cos_theta, sin_theta, curv, curvp;

              tmp1=std::sqrt(one+ycp*ycp);
              tmp13=tmp1*tmp1*tmp1;
              cos_theta=one/tmp1;
              sin_theta=ycp/tmp1;
              curv=ycpp/tmp13;
              curvp=ycppp/tmp13-three*ycp*tmp1*curv*curv;

              // calculate the info
              x(0)=surf_sign*(xi-delta*sin_theta);
              x(1)=yc+surf_sign*delta*cos_theta;
              xp(0)=surf_sign-deltap*sin_theta-delta*curv;
              xp(1)=ycp*(surf_sign-delta*curv)+deltap*cos_theta;
              xpp(0)=-surf_sign*(deltapp*sin_theta+two*deltap*curv+delta*curvp);
              xpp(1)=ycpp*(one-surf_sign*delta*curv)+surf_sign*(deltapp*cos_theta-ycp*(two*deltap*curv+delta*curvp));
            }

            point_type tangent(const data_type &xi) const
            {
              point_type tgt(fp(xi));

              tgt.normalize();
              return tgt;
            }

          protected:

            virtual void recalc_thickness_coefficients(){};
            virtual void recalc_camber_coefficients(){};

            virtual void calc_camber(data_type &y, data_type &yp, data_type &ypp, data_type &yppp, const data_type &xi) const = 0;
            virtual void calc_thickness(data_type &y, data_type &yp, data_type &ypp, const data_type &xi) const = 0;

        protected:
            data_type thickness;    // thickness index (integer [00,99] with practical limit of [00, 30].
                                    // Index is interpreted as 100 times the maximum thickness.
            bool sharp_te; // flag to indicate if the trailing edge should be sharp
        };
      }
    }
  }
}
#endif
