/*********************************************************************************
* Copyright (c) 2013 David D. Marshall <ddmarsha@calpoly.edu>
*
* All rights reserved. This program and the accompanying materials
* are made available under the terms of the Eclipse Public License v1.0
* which accompanies this distribution, and is available at
* http://www.eclipse.org/legal/epl-v10.html
*
* Contributors:
*    David D. Marshall - initial code and implementation
********************************************************************************/

#ifndef eli_constants_math_hpp
#define eli_constants_math_hpp

#include <cmath>

#include "eli/code_eli.hpp"

namespace eli
{
  namespace constants
  {
    template <typename T__>
    class math
    {
    };

    template<>
    class math <float>
    {
      public:
        static float exp()     {return 2.7182818f;}
        static float ln_two()  {return 0.6931472f;}
        static float log_exp() {return 0.43429447f;}

        static float pi()             {return 3.1415927f;}
        static float two_pi()         {return 6.2831853f;}
        static float pi_by_two()      {return 1.57079633f;}
        static float pi_by_four()     {return 0.78539816f;}
        static float pi_squared()     {return 9.8696046f;}
        static float pi_cubed()       {return 31.006279f;}
        static float sqrt_pi()        {return 1.7724539f;}
        static float cbrt_pi()        {return 1.4645919f;}
        static float one_by_pi()      {return 0.31830988f;}
        static float two_by_pi()      {return 0.63661977f;}
        static float one_by_sqrt_pi() {return 0.56418958f;}
        static float two_by_sqrt_pi() {return 1.12837916f;}

        static float sqrt_two()        {return 1.41421356f;}
        static float sqrt_two_by_two() {return 0.70710678f;}

        // Classical constant value
        // k = 4./3.;
        // Fraction of radius to place cubic control point
        // f = k * tan( theta * 0.25 );
        // Where theta is the angle of the arc.
        //
        // Improved constant due to analysis similar to: http://spencermortensen.com/articles/bezier-circle/
        // conducted symbolically in Matlab, leading to expression
        // (3*2^(1/2)*c)/8 + 2^(1/2)/2 + (abs(3*c - 1)*(3*c^4 + 8*c^3 + 12*c^2 - 24*c + 8)^(1/2))/(3*c - 2)^2 - 2 == 0
        // Which was solved numerically using Wolfram Alpha to 100 digits:
        // FindRoot[-2 + 1/Sqrt[2] + (3 c)/(4 Sqrt[2]) + (Sqrt[8 - 24 c + 12 c^2 + 8 c^3 + 3 c^4] Abs[-1 + 3 c])/(-2 + 3 c)^2 == 0, {c, 0.55166, 0.552473}, WorkingPrecision -> 100]
        // Improved constant value to 100 digits (for 90 deg arc).
        // c = 0.5519150244935105707435627227925666423361803947243088973369805374675870988527781759268533834535800161;
        // k = c / ( sqrt_two()-1 );
        static double cubic_bezier_circle_const()  {return 1.332440737409712163877261103501904586378054505346634846763928262809148692267066998533212237794404131f;}
    };

    template <>
    class math <double>
    {
      public:
        static const double exp1;

      public:
        static double exp()     {return 2.718281828459045;}
        static double ln_two()  {return 0.6931471805599453;}
        static double log_exp() {return 0.4342944819032518;}

        static double pi()             {return 3.141592653589793;}
        static double two_pi()         {return 6.283185307179586;}
        static double pi_by_two()      {return 1.5707963267948966;}
        static double pi_by_four()     {return 0.7853981633974483;}
        static double pi_squared()     {return 9.8696044010893586;}
        static double pi_cubed()       {return 31.006276680299817;}
        static double sqrt_pi()        {return 1.7724538509055159;}
#if defined(_MSC_VER)
# if (_MSC_VER==1800)
		static double cbrt_pi()        {return 1.4645918875615233;}
# else
		static double cbrt_pi()        { return 1.4645918875615231; }
# endif
#elif defined(__INTEL_COMPILER)
		static double cbrt_pi()        {return 1.4645918875615231;}
#elif defined(__clang__)
        static double cbrt_pi()        {return 1.4645918875615233;}
#elif defined(__GNUC__)
# if (__GNUC__==4) && (__GNUC_MINOR__>=7) && (defined NDEBUG)
        static double cbrt_pi()        {return 1.4645918875615231;}
# else
        static double cbrt_pi()        {return 1.4645918875615233;}
# endif
#else
        static double cbrt_pi()        {return 1.4645918875615231;}
#endif
        static double one_by_pi()      {return 0.31830988618379067;}
        static double two_by_pi()      {return 0.63661977236758134;}
        static double one_by_sqrt_pi() {return 0.56418958354775629;}
        static double two_by_sqrt_pi() {return 1.1283791670955126;}

        static double sqrt_two()        {return 1.41421356237309504;}
        static double sqrt_two_by_two() {return 0.70710678118654752;}

        // Classical constant value
        // k = 4./3.;
        // Fraction of radius to place cubic control point
        // f = k * tan( theta * 0.25 );
        // Where theta is the angle of the arc.
        //
        // Improved constant due to analysis similar to: http://spencermortensen.com/articles/bezier-circle/
        // conducted symbolically in Matlab, leading to expression
        // (3*2^(1/2)*c)/8 + 2^(1/2)/2 + (abs(3*c - 1)*(3*c^4 + 8*c^3 + 12*c^2 - 24*c + 8)^(1/2))/(3*c - 2)^2 - 2 == 0
        // Which was solved numerically using Wolfram Alpha to 100 digits:
        // FindRoot[-2 + 1/Sqrt[2] + (3 c)/(4 Sqrt[2]) + (Sqrt[8 - 24 c + 12 c^2 + 8 c^3 + 3 c^4] Abs[-1 + 3 c])/(-2 + 3 c)^2 == 0, {c, 0.55166, 0.552473}, WorkingPrecision -> 100]
        // Improved constant value to 100 digits (for 90 deg arc).
        // c = 0.5519150244935105707435627227925666423361803947243088973369805374675870988527781759268533834535800161;
        // k = c / ( sqrt_two()-1 );
        static double cubic_bezier_circle_const()  {return 1.332440737409712163877261103501904586378054505346634846763928262809148692267066998533212237794404131;}
    };

    template<>
    class math <long double>
    {
      public:
        static long double exp()     {return 2.7182818284590452354L;}
        static long double ln_two()  {return 0.69314718055994530942L;}
        static long double log_exp() {return 0.43429448190325182766L;}

        static long double pi()             {return 3.1415926535897932385L;}
        static long double two_pi()         {return 6.2831853071795864770L;}
        static long double pi_by_two()      {return 1.57079632679489661923L;}
        static long double pi_by_four()     {return 0.78539816339744830962L;}
        static long double pi_squared()     {return 9.869604401089358619L;}
#ifdef _MSC_VER
        static long double pi_cubed()       {return 31.006276680299817L;}
        static long double sqrt_pi()        {return 1.7724538509055159L;}
# if (_MSC_VER==1800)
		static long double cbrt_pi()        {return 1.4645918875615233L;}
# else
		static long double cbrt_pi()        {return 1.4645918875615231L;}
# endif
#else
#ifdef __aarch64__
        static long double pi_cubed()       {return 31.0062766802998162063L;}
        static long double sqrt_pi()        {return 1.7724538509055158819L;}
#else
        static long double pi_cubed()       {return 31.0062766802998201763L;}
        static long double sqrt_pi()        {return 1.7724538509055160273L;}
#endif
        static long double cbrt_pi()        {return 1.464591887561523263L;}
#endif
        static long double one_by_pi()      {return 0.31830988618379067154L;}
        static long double two_by_pi()      {return 0.6366197723675813431L;}
        static long double one_by_sqrt_pi() {return 0.5641895835477562869L;}
        static long double two_by_sqrt_pi() {return 1.1283791670955125738L;}

        static long double sqrt_two()        {return 1.4142135623730950488L;}
        static long double sqrt_two_by_two() {return 0.7071067811865475244L;}

        // Classical constant value
        // k = 4./3.;
        // Fraction of radius to place cubic control point
        // f = k * tan( theta * 0.25 );
        // Where theta is the angle of the arc.
        //
        // Improved constant due to analysis similar to: http://spencermortensen.com/articles/bezier-circle/
        // conducted symbolically in Matlab, leading to expression
        // (3*2^(1/2)*c)/8 + 2^(1/2)/2 + (abs(3*c - 1)*(3*c^4 + 8*c^3 + 12*c^2 - 24*c + 8)^(1/2))/(3*c - 2)^2 - 2 == 0
        // Which was solved numerically using Wolfram Alpha to 100 digits:
        // FindRoot[-2 + 1/Sqrt[2] + (3 c)/(4 Sqrt[2]) + (Sqrt[8 - 24 c + 12 c^2 + 8 c^3 + 3 c^4] Abs[-1 + 3 c])/(-2 + 3 c)^2 == 0, {c, 0.55166, 0.552473}, WorkingPrecision -> 100]
        // Improved constant value to 100 digits (for 90 deg arc).
        // c = 0.5519150244935105707435627227925666423361803947243088973369805374675870988527781759268533834535800161;
        // k = c / ( sqrt_two()-1 );
        static double cubic_bezier_circle_const()  {return 1.332440737409712163877261103501904586378054505346634846763928262809148692267066998533212237794404131L;}
    };
  }
}

#endif
