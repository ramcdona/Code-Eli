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

#ifndef eli_geom_curve_piecewise_linear_creator_hpp
#define eli_geom_curve_piecewise_linear_creator_hpp

#include <vector>

#include "eli/code_eli.hpp"

#include "eli/geom/curve/piecewise_creator_base.hpp"
#include "eli/geom/curve/piecewise.hpp"
#include "eli/geom/curve/bezier.hpp"

namespace eli
{
  namespace geom
  {
    namespace curve
    {
      template<typename data__, unsigned short dim__, typename tol__>
      class piecewise_linear_creator : public piecewise_creator_base<data__, dim__, tol__>
      {
        public:
          typedef piecewise_creator_base<data__, dim__, tol__> base_class_type;
          typedef typename base_class_type::data_type data_type;
          typedef typename base_class_type::point_type point_type;
          typedef typename base_class_type::index_type index_type;
          typedef typename base_class_type::tolerance_type tolerance_type;

          piecewise_linear_creator() : piecewise_creator_base<data_type, dim__, tolerance_type>(2, 0), corner(3) {}
          piecewise_linear_creator(const index_type &ns) : piecewise_creator_base<data_type, dim__, tolerance_type>(ns, 0), corner(ns+1) {}
          piecewise_linear_creator(const piecewise_linear_creator<data_type, dim__, tolerance_type> &ppc)
            : piecewise_creator_base<data_type, dim__, tolerance_type>(ppc), corner(ppc.corner) {}
          ~piecewise_linear_creator() {}

          void set_corner(const point_type &c, const index_type &i)
          {
            if ((i>=0) && (i<static_cast<index_type>(corner.size())))
              corner[i]=c;
            else
              assert(false);
          }
          point_type get_corner(const index_type &i) const
          {
            if ((i<0) || (i>=static_cast<index_type>(corner.size())))
            {
              return corner[0];
              assert(false);
            }
            return corner[i];
          }

          void set_number_segments(const index_type &ns)
          {
            this->set_num_segs(ns);
            corner.resize(ns+1);
          }

          virtual bool create(piecewise<bezier, data_type, dim__, tolerance_type> &pc) const
          {
            typedef piecewise<bezier, data_type, dim__, tolerance_type> piecewise_curve_type;
            typedef typename piecewise_curve_type::curve_type curve_type;
            typedef typename piecewise_curve_type::error_code error_code;

            pc.clear();

            curve_type c(1);
            error_code err;
            index_type nsegs(this->get_number_segments());

            // do sanity check
            if (corner.size()!=static_cast<size_t>(nsegs+1))
            {
              assert(false);
              return false;
            }

            // set the start parameter
            pc.set_t0(this->get_t0());

            // set the first n edges
            for (index_type i=0; i<nsegs; ++i)
            {
              c.set_control_point(corner[i], 0);
              c.set_control_point(corner[i+1], 1);
              err=pc.push_back(c, this->get_segment_dt(i));
              if (err!=piecewise_curve_type::NO_ERRORS)
              {
                pc.clear();
                pc.set_t0(0);
                return false;
              }
            }

            return true;
          }

        private:
          std::vector<point_type, Eigen::aligned_allocator<point_type> > corner;
      };
    }
  }
}
#endif
