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

namespace eli
{
  namespace geom
  {
    namespace curve
    {
      namespace pseudo
      {
        template<typename data__>
        class five_digit_mod {
          public:
            typedef data__ data_type;
            typedef Eigen::Matrix<data_type, 1, 2> point_type;
            typedef Eigen::Matrix<data_type, 5, 1> coefficient_type;
            typedef typename point_type::Index index_type;
            typedef std::complex <data_type> complex_data_type;

          public:
            five_digit_mod() : thickness(.12), cli(0.3), camber_loc(0.15), thickness_loc(.3), le_radius_indx(6), sharp_te(false)
            {
              recalc_thickness_coefficients();
              recalc_camber_coefficients();
            }

            five_digit_mod(const five_digit_mod <data_type> &fs)
                    : thickness(fs.thickness), cli(fs.cli), camber_loc(fs.camber_loc), thickness_loc(fs.thickness_loc),
                      le_radius_indx(fs.le_radius_index), sharp_te(fs.sharp_te)
            {
              recalc_thickness_coefficients();
              recalc_camber_coefficients();
            }

            five_digit_mod( const data_type &t, const data_type &c, const data_type &xc, const data_type &lei, const data_type &t_loc, bool fl )
            {
              set( t, c, xc, lei, t_loc, fl );
            }

            coefficient_type get_thickness_coefficients() const { return a; }

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

            bool set_cli(const data_type &c)
            {
              cli = c;
              return true;
            }

            data_type get_cli() const { return cli; }

            bool set_camber_loc(const data_type &xc)
            {
              camber_loc = xc;
              recalc_camber_coefficients();
              return true;
            }

            data_type get_camber_loc() { return camber_loc; }

            // Valid values of thickness location greater than 0 and less than 100
            bool set_thickness_loc(const data_type &t_loc)
            {
                if ((t_loc>=.2) && (t_loc<=.6))
                {
                    thickness_loc=t_loc;
                    recalc_thickness_coefficients();
                    return true;
                }
                return false;
            }
            data_type get_maximum_thickness_location() const {return thickness_loc;}

            // Valid values of le radius index are greater than 0 and less than 8
            bool set_leradindex(const data_type &lei)
            {
                if ((lei>=0) && (lei<=9))
                {
                    le_radius_indx=lei;
                    recalc_thickness_coefficients();
                    return true;
                }
                return false;
            }
            data_type get_leradindex() const {return le_radius_indx;}


            bool set( const data_type &t, const data_type &c, const data_type &xc, const data_type &lei, const data_type &t_loc, bool fl )
            {
                if ( (t < 0) || (t > 1) ) {
                    return false;
                }
                thickness = t;

                cli = c;

                camber_loc = xc;

                sharp_te = fl;

                if ((lei<0) || (lei>9))
                {
                  return false;
                }
                le_radius_indx=lei;

                if ((t_loc<.2) || (t_loc>.6))
                {
                  return false;
                }
                thickness_loc=t_loc;


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
              // Calculate 4-digit modified curve coefficients for 20% thick airfoil
              // With given maximum thickness position and leading edge radius index.
              // Flag for sharp or 'thick' trailing edge.

              data_type a0, a1, a2, a3;
              data_type d0, d1, d2, d3;

              data_type tocref = 0.2;

              data_type mt = thickness_loc;

              data_type xmtsq = mt * mt;
              data_type xmtcu = xmtsq * mt;
              data_type srxmt = sqrt( mt );

              data_type rle = ( 1.1019 / 36.0 ) * tocref * tocref * le_radius_indx * le_radius_indx;

              if ( le_radius_indx > 8.0 )
              {
                data_type rle8 = ( 1.1019 / 36.0 ) * tocref * tocref * 64.0; // Quadratic up to index of 8
                data_type rle9 = ( 1.1019 ) * tocref * tocref * 3.0;         // Three times normal radius.

                // Linear interpolation of radius between index 8 and 9.
                data_type di = le_radius_indx - 8.0;
                data_type drle = rle9 - rle8;
                rle = rle8 + di * drle;
              }

              if ( !sharp_trailing_edge() ) // 'Official' airfoil has TE semi-bluntness of 1% of thickness.
              {
                d0 = 0.01 * tocref;
              }
              else
              {
                d0 = 0.0;
              }

              // The 'official' definition of modified airfoils only exists at five locations of maximum thickness.
              // The d1 coefficient is tabulated at these five points.
              //
              //   m     d1   Riegels  Current
              //  0.2  0.200  0.19990  0.20025
              //  0.3  0.234  0.23364  0.23312
              //  0.4  0.315  0.31442  0.31609
              //  0.5  0.465  0.46434  0.46439
              //  0.6  0.700  0.72189  0.70012
              //
              // Riegels developed an approximate equation that allows continuous variation of the maximum thickness
              // location over the same range.  Riegels' equation has been used (uncited) in NASA published computer
              // programs and reports (TM X 3284 and TM 4741) for generating the NACA modified thickness forms.  It is
              // identified by name in Bill Mason's course notes on airfoil geometry.  Ralph Carmichael distributes
              // an airfoil generating program via PDAS that uses a polynomial fit as an improvement over Riegels.
              //
              // Considerable effort was extended on developing an improved approximation with the hopes of enabling
              // airfoils with maximum thickness location beyond the [0.2, 0.6] range.  Those efforts ultimately proved
              // unsuccessful at generating reasonable airfoils over the extended range -- it appears that the d1
              // coefficient was not the problem.
              //
              // A Matlab program was written to find the nonlinear least squares best fit to Riegels equation.  This
              // optimal sef of coefficients is used here.
              //
              // Riegels equation.
              // d1 = ( 0.1 * ( 2.24 - 5.42 * mt + 12.3 * xmtsq ) / ( 1.0 - 0.878 * mt ) );
              //
              // Riegels, F. W., "Aerofoil Sections, Results from Wind-Tunnel Investigations, Theoretical Foundations",
              // Translated from the German by D. G. Randall, Butterworths, London, 1961.

              // Best fit to Riegels form
              data_type c1, c2, c3, c4, c5;
              c1 =  0.244364095382286;
              c2 = -0.677199764201376;
              c3 =  1.609809719636767;
              c4 = -0.672612098884539;
              d1 = ( c1 + c2 * mt + c3 * xmtsq ) / ( 1.0 + c4 * mt );

              // With d1 (which determines trailing edge angle) determined by the location of the maximum thickness, and
              // a0 determined by the leading edge radius, the geometric constraints of C2 matching at the max thickness
              // point is enough to set all the remaining coefficients.  These equations had to be re-derived to work
              // for non-default values of the TE thickness.  This also allows for more precision.

              /*
              Maxima computer algebra commands to derive airfoil coefficients.

              set_display('none)$
              yr(x):=d0+d1*(1-x)+d2*(1-x)^2+d3*(1-x)^3;
              yrp(x):=diff(yr(x),x,1);
              yrpp(x):=diff(yr(x),x,2);

              sol:solve([yr(m)=t/2, yrp(m)=0],[d2,d3]);
              s2:rhs(subst(sol[1][1],d2,yrpp(m)));
              s3:ratsimp(rhs(subst(sol[1][2],d3,s2)));

              yf(x):= a0 * sqrt(x) + a1 * x + a2 * x^2 + a3 * x^3;
              yfp(x):=diff(yf(x),x,1);
              yfpp(x):=diff(yf(x),x,2);
              solve([yf(m)=t/2, yfp(m)=0, yfpp(m)=K],[a3,a2,a1]);
              solve([yf(m)=t/2, yfp(m)=0],[a2,a1]);
              solve([yf(m)=t/2],[a1]);
              */

              d3 = ( tocref + d1 * mt - d1 - 2.0 * d0 ) / ( xmtcu - 3.0 * xmtsq + 3.0 * mt - 1.0 );
              d2 = ( 3.0 * tocref + 4.0 * d1 * mt - 4.0 * d1 - 6.0 * d0 ) / ( 2.0 * xmtsq - 4.0 * mt + 2.0 );

              a0 = sqrt( 2.0 * rle );

              data_type ypp;
              ypp = -( 3.0 * tocref + 2.0 * d1 * mt - 2.0 * d1 - 6.0 * d0 ) / ( xmtsq - 2.0 * mt + 1.0 );

              a3 = ( 4.0 * tocref + 4.0 * ypp * xmtsq - 3.0 * a0 * srxmt ) / ( 8.0 * xmtcu );
              a2 = -( srxmt * tocref + 4.0 * a3 * srxmt * xmtcu - a0 * mt ) / ( 2.0 * srxmt * xmtsq );
              a1 = ( tocref - 2.0 * a3 * xmtcu - 2.0 * a2 * xmtsq - 2.0 * a0 * srxmt ) / ( 2.0 * mt );

              a(0) = a0;
              a(1) = a1;
              a(2) = a2;
              a(3) = a3;

              d(0) = d0;
              d(1) = d1;
              d(2) = d2;
              d(3) = d3;
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

              const data_type zero(0), one(1), two(2), three(3), four(4), six(6), half(one/two), quarter(one/four);
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
              else if ((xi==1) && sharp_trailing_edge())  // This is clearly wrong.
              {
                y=zero;
                yp=trat*(a.sum()-half*a(0));
                ypp=trat*(-quarter*a(0)+two*a(2)+six*a(3));
                return;
              }

              if (xi<thickness_loc)
              {
                const data_type xi2(xi*xi), xi3(xi*xi2), sqrtxi(std::sqrt(xi));

                y=trat*(a(0)*sqrtxi+a(1)*xi+a(2)*xi2+a(3)*xi3);
                yp=trat*(half*a(0)/sqrtxi+a(1)+two*a(2)*xi+three*a(3)*xi2);
                ypp=trat*(-quarter*a(0)/sqrtxi/xi+two*a(2)+six*a(3)*xi);
              }
              else
              {
                const data_type omx( one - xi );
                const data_type omx2(omx*omx), omx3(omx*omx2);

                y=trat*( d(0)+d(1)*omx+d(2)*omx2+d(3)*omx3 );
                yp=trat*( -d(1) - d(2)*two*omx - d(3)*three*omx2 );
                ypp=trat*( d(2)*two + d(3)*six*omx );
              }
            }

          private:
            data_type thickness;    // thickness index (integer [00,99] with practical limit of [00, 30].
                                    // Index is interpreted as 100 times the maximum thickness.
            data_type cli;          // Ideal lift coefficient.
            data_type camber_loc;   // chord-wise location of maximum camber index
            data_type thickness_loc; // chord-wise location of maximum thickness index
            data_type le_radius_indx; // Leading edge radius index

            data_type xtrans;
            data_type k;
            bool sharp_te; // flag to indicate if the trailing edge should be sharp
            coefficient_type a; // coefficients for thickness distribution
            coefficient_type d; // coefficients for thickness distribution
        };
      }
    }
  }
}
#endif
