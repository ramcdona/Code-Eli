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

#ifndef eli_geom_curve_pseudo_four_digit_thickness_hpp
#define eli_geom_curve_pseudo_four_digit_thickness_hpp

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
        class four_digit_thickness : virtual public naca_af<data__> {
        public:
            typedef naca_af<data__> base_class_type;
            typedef typename base_class_type::data_type data_type;
            typedef typename base_class_type::point_type point_type;
            typedef typename base_class_type::coefficient_type coefficient_type;
            typedef typename base_class_type::index_type index_type;
            typedef typename base_class_type::complex_data_type complex_data_type;

          public:

            coefficient_type get_thickness_coefficients() const {return a;}

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
              if (!this->sharp_trailing_edge())
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

            void calc_thickness(data_type &y, data_type &yp, data_type &ypp, const data_type &xi) const
            {
              // check to make sure given valid parametric value
              assert((xi>=0) && (xi<=1));

              const data_type zero(0), one(1), two(2), three(3), four(4), six(6), twelve(12), half(one/two), quarter(one/four);
              const data_type xi2(xi*xi), xi3(xi*xi2), xi4(xi2*xi2), sqrtxi(std::sqrt(xi));
              const data_type trat(this->thickness/static_cast<data_type>(0.20));

              // short circuit for no thickness
              if (this->thickness==0)
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
              else if ((xi==1) && this->sharp_trailing_edge())
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

            coefficient_type a; // coefficients for thickness distribution
        };
      }
    }
  }
}
#endif
