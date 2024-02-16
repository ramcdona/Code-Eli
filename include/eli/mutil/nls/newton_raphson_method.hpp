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

#ifndef eli_mutil_nls_newton_raphson_method_hpp
#define eli_mutil_nls_newton_raphson_method_hpp

#include "eli/code_eli.hpp"

#include "eli/mutil/nls/iterative_root_base_constrained.hpp"

namespace eli
{
  namespace mutil
  {
    namespace nls
    {
      template<typename data__>
      class newton_raphson_method : public iterative_root_base_constrained<data__>
      {
        public:
          typedef data__ data_type;

          static const int hit_constraint = 101;

        private:
          template<typename f__, typename g__>
          struct wrapper_functor
          {
            f__ fun;
            g__ fprime;

            void operator()( data_type &fx, data_type &fpx, const data_type &x ) const
            {
              fx = fun( x );
              fpx = fprime( x );
            }
          };

        public:
          newton_raphson_method() : iterative_root_base_constrained<data_type>(), x0(0)
          {
          }

          newton_raphson_method(const newton_raphson_method<data_type> &nrm)
          : iterative_root_base_constrained<data_type>(nrm), x0(nrm.x0)
          {
          }

          ~newton_raphson_method()
          {
          }

          void set_initial_guess(const data_type &xg)
          {
            x0=xg;
          }

          const data_type & get_initial_guess() const
          {
            return x0;
          }

          template<typename f__, typename g__>
          int find_root( data_type &root, const f__ &fun, const g__ &fprime, const data_type &f0 ) const
          {
            wrapper_functor< f__, g__ > ffprime;
            ffprime.fun = fun;
            ffprime.fprime = fprime;
            return find_root( root, ffprime, f0 );
          }

          template<typename fg__>
          int find_root(data_type &root, const fg__ &ffprime, const data_type &f0) const
          {
            data_type fx, fpx;
            data_type x(x0), eval, eval_abs, eval_abs2, dx(1);
            ffprime( fx, fpx, x0 );
            typename iterative_root_base<data__>::iteration_type count;

            // calculate the function evaluated at the initial location
            eval=fx-f0;
            eval_abs=std::abs(eval);
            if (this->test_converged(0, eval_abs/f0, eval_abs, 0, 0))
            {
              root=x;
              return iterative_root_base<data__>::converged;
            }

            eval_abs2=0;
            count=0;
            while (!this->test_converged(count, eval_abs/f0, eval_abs, eval_abs2/x0, eval_abs2) && (std::abs(dx)>0))
            {
              if (fpx==0)
                return iterative_root_base<data__>::no_root_found;

              dx=this->calculate_delta_factor(x, -eval/fpx);
              x+=dx;
              ffprime( fx, fpx, x );
              eval=fx-f0;
              eval_abs=std::abs(eval);
              eval_abs2=std::abs(dx);

              ++count;
            }

            root=x;
            if (this->max_iteration_reached(count))
              return this->max_iteration; // could not converge
            if (dx==0)
              return this->hit_constraint; // constraints limited convergence

            return this->converged;
          }

        private:
          data_type x0;
      };
    }
  }
}
#endif
