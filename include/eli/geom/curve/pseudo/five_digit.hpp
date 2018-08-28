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

#ifndef eli_geom_curve_pseudo_five_digit_hpp
#define eli_geom_curve_pseudo_five_digit_hpp

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
        class five_digit {
        public:
            typedef data__ data_type;
            typedef Eigen::Matrix<data_type, 1, 2> point_type;
            typedef Eigen::Matrix<data_type, 5, 1> coefficient_type;
            typedef typename point_type::Index index_type;
            typedef std::complex <data_type> complex_data_type;

        public:
            five_digit() : thickness(.12), cli(0.3), camber_loc(0.15), sharp_te(false) {
                recalc_thickness_coefficients();
                recalc_camber_coefficients();
            }

            five_digit(const five_digit <data_type> &fs)
                    : thickness(fs.thickness), cli(fs.cli), camber_loc(fs.camber_loc),
                      sharp_te(fs.sharp_te) {
                recalc_thickness_coefficients();
                recalc_camber_coefficients();
            }

            five_digit( const data_type &t, const data_type &c, const data_type &xc, bool fl )
            {
                set( t, c, xc, fl );
            }

            coefficient_type get_thickness_coefficients() const { return a; }

            data_type get_t0() const { return static_cast<data_type>(-1); }

            data_type get_tmax() const { return static_cast<data_type>(1); }

            void set_sharp_trailing_edge(bool fl) {
                sharp_te = fl;
                recalc_thickness_coefficients();
            }

            bool sharp_trailing_edge() const { return sharp_te; }

            // Valid values of thickness are greater than 0 and less than 100
            bool set_thickness(const data_type &t) {
                if ((t >= 0) && (t <= 1)) {
                    thickness = t;
                    return true;
                }
                return false;
            }

            data_type get_thickness() const { return thickness; }

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

            bool set( const data_type &t, const data_type &c, const data_type &xc, bool fl )
            {
                if ( (t < 0) || (t > 1) ) {
                    return false;
                }
                thickness = t;

                cli = c;

                camber_loc = xc;

                sharp_te = fl;

                recalc_thickness_coefficients();
                recalc_camber_coefficients();

                return true;
            }

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

            void recalc_thickness_coefficients()
            {
              typedef Eigen::Matrix<data_type, 5, 5> coefficient_matrix_type;

              coefficient_matrix_type coef_mat;
              coefficient_type rhs;
              coefficient_type orig_a;

              // set the specified coefficients
              orig_a << static_cast<data_type>(0.2969),
                        static_cast<data_type>(-0.1260),
                        static_cast<data_type>(-0.3516),
                        static_cast<data_type>(0.2843),
                        static_cast<data_type>(-0.1015);

              // if blunt trailing edge, then use specified coefficients
              if (!sharp_trailing_edge())
              {
                a=orig_a;
                return;
              }

              // if want sharp trailing edge then find "actual" constraints that
              // are placed on thickness distribution

              // calculate the constraint coefficients
              calc_four_digit_args(coef_mat, 0, static_cast<data_type>(1));        // (1) trailing edge location
              calc_four_digit_der_args(coef_mat, 1, static_cast<data_type>(1));    // (2) trailing edge slope
              calc_four_digit_args(coef_mat, 2, static_cast<data_type>(0.1));      // (3) leading edge shape
              calc_four_digit_args(coef_mat, 3, static_cast<data_type>(0.3));      // (4) thickness at x/c=0.3 (should be max thickness location, but isn't)
              calc_four_digit_der_args(coef_mat, 4, static_cast<data_type>(0.3));  // (5) slope at x/c=0.3 (should be zero slope at max thickness, but isn't)

              // calculate the corresponding constraints for the blunt trailing edge
              rhs=coef_mat*orig_a;

              // correct the trailing edge thickness constraint to zero while leaving the rest un changed
              rhs(0)=static_cast<data_type>(0);
              a=coef_mat.lu().solve(rhs);
            }

            template<typename Derived1>
            static void calc_four_digit_args(Eigen::MatrixBase<Derived1> &A, const typename Derived1::Index &i, const data_type &xi)
            {
              data_type xi2(xi*xi), xi3(xi2*xi), xi4(xi3*xi);

              A(i,0)=std::sqrt(xi);
              A(i,1)=xi;
              A(i,2)=xi2;
              A(i,3)=xi3;
              A(i,4)=xi4;
            }

            template<typename Derived1>
            static void calc_four_digit_der_args(Eigen::MatrixBase<Derived1> &A, const typename Derived1::Index &i, const data_type &xi)
            {
              data_type xi2(xi*xi), xi3(xi2*xi);

              A(i,0)=static_cast<data_type>(0.5)/std::sqrt(xi);
              A(i,1)=static_cast<data_type>(1);
              A(i,2)=static_cast<data_type>(2)*xi;
              A(i,3)=static_cast<data_type>(3)*xi2;
              A(i,4)=static_cast<data_type>(4)*xi3;
            }

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
              if (cli==0)
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

            void calc_thickness(data_type &y, data_type &yp, data_type &ypp, const data_type &xi) const
            {
              // check to make sure given valid parametric value
              assert((xi>=0) && (xi<=1));

              const data_type zero(0), one(1), two(2), three(3), four(4), six(6), twelve(12), half(one/two), quarter(one/four);
              const data_type xi2(xi*xi), xi3(xi*xi2), xi4(xi2*xi2), sqrtxi(std::sqrt(xi));
              const data_type trat(thickness/static_cast<data_type>(0.20));

              // short circuit for no thickness
              if (thickness==0)
              {
                y=zero;
                yp=zero;
                ypp=zero;
                return;
              }

              if (xi==0)
              {
                y=zero;
                yp=one/std::numeric_limits<data_type>::epsilon();
                ypp=yp;
                return;
              }
              else if ((xi==1) && sharp_trailing_edge())
              {
                y=zero;
                yp=trat*(a.sum()-half*a(0));
                ypp=trat*(-quarter*a(0)+two*a(2)+six*a(3)+twelve*a(4));
                return;
              }

              y=trat*(a(0)*sqrtxi+a(1)*xi+a(2)*xi2+a(3)*xi3+a(4)*xi4);
              yp=trat*(half*a(0)/sqrtxi+a(1)+two*a(2)*xi+three*a(3)*xi2+four*a(4)*xi3);
              ypp=trat*(-quarter*a(0)/sqrtxi/xi+two*a(2)+six*a(3)*xi+twelve*a(4)*xi2);
            }

          private:
            data_type thickness;    // thickness index (integer [00,99] with practical limit of [00, 30].
                                    // Index is interpreted as 100 times the maximum thickness.
            data_type cli;          // Ideal lift coefficient.
            data_type camber_loc;   // chord-wise location of maximum camber index

            data_type xtrans;
            data_type k;
            bool sharp_te; // flag to indicate if the trailing edge should be sharp
            coefficient_type a; // coefficients for thickness distribution
        };
      }
    }
  }
}
#endif
