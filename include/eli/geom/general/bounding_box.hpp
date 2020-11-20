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

#ifndef eli_geom_general_bounding_box_hpp
#define eli_geom_general_bounding_box_hpp

#include "eli/code_eli.hpp"

#include "eli/util/tolerance.hpp"

namespace eli
{
  namespace geom
  {
    namespace general
    {
      template<typename data__, unsigned short dim__, typename tol__=eli::util::tolerance<data__> >
      class bounding_box
      {
        public:
          typedef data__ data_type;
          typedef unsigned short dimension_type;
          typedef Eigen::Matrix<data_type, 1, dim__> point_type;
          typedef typename point_type::Index index_type;
          typedef tol__ tolerance_type;

        public:
          bounding_box() : empty(true) {}
          bounding_box(const point_type &p) : empty(false), pmin(p), pmax(p) {}
          bounding_box(const bounding_box<data_type, dim__> &bb) : empty(bb.empty), pmin(bb.pmin), pmax(bb.pmax) {}
          ~bounding_box() {}

          bounding_box<data_type, dim__> & operator=(const bounding_box<data_type, dim__> &bb)
          {
            if (this!=&bb)
            {
              empty=bb.empty;
              pmin=bb.pmin;
              pmax=bb.pmax;
            }

            return (*this);
          }

          bool operator==(const bounding_box<data_type, dim__> &bb)
          {
            // if compared against self then true
            if (this==&bb)
              return true;

            // if both are empty sets then true
            if (empty && bb.empty)
              return true;

            // start comparing member data
            if (empty!=bb.empty)
              return false;
            if (pmin!=bb.pmin)
              return false;
            if (pmax!=bb.pmax)
              return false;

            return true;
          }

          bool operator!=(const bounding_box<data_type, dim__> &bb)
          {
            return !operator==(bb);
          }

          static dimension_type dimension() {return dim__;}

          bool empty_set() const {return empty;}

          void set_min(const point_type &pm)
          {
            if (empty_set())
            {
              pmax=pm;
              empty=false;
            }
            pmin=pm;
          }
          point_type get_min() const {return pmin;}

          void set_max(const point_type &pm)
          {
            if (empty_set())
            {
              pmin=pm;
              empty=false;
            }
            pmax=pm;
          }
          point_type get_max() const {return pmax;}

          // Get negative and positive points in relation to a specified normal vector direction
          void get_nppts( const point_type &nvec, point_type &npt, point_type &ppt ) const
          {
            for (index_type i=0; i<dim__; ++i)
            {
              npt(i) = pmax(i);
              ppt(i) = pmin(i);
              if ( nvec(i) >= 0.0 )
              {
                npt(i) = pmin(i);
                ppt(i) = pmax(i);
              }
            }
          }

          // Test whether bounding box intersects a plane given as a point/vector.
          bool intersect_plane( const point_type &pt, const point_type &nvec ) const
          {
            point_type npt, ppt;
            get_nppts( nvec, npt, ppt );

            data_type d = - nvec.dot( pt );

            if ( d + nvec.dot( ppt ) < 0.0 )
            {
              return false;
            }

            if ( d + nvec.dot( npt ) > 0.0 )
            {
              return false;
            }

            return true;
          }

          bool intersect_axis( const point_type &p, const index_type iax ) const
          {
            for ( index_type i = 0; i < dim__; ++i )
            {
              if ( i != iax )  // Ignore iax direction.  Otherwise a point intersection test.
              {
                if ( ( p(i) < pmin(i) ) || ( p(i) > pmax(i) ) )
                  return false;
              }
            }

            return true;
          }

          void clear()
          {
            empty=true;
            pmin.setZero();
            pmax.setZero();
          }

          bool add(const point_type &p)
          {
            index_type i;
            bool changed(false);

            // handle special case of no bounding box
            if (empty_set())
            {
              set_min(p);
              set_max(p);
              return true;
            }

            for (i=0; i<dim__; ++i)
            {
              if (p(i)<pmin(i))
              {
                changed=true;
                pmin(i)=p(i);
              }
              if (p(i)>pmax(i))
              {
                changed=true;
                pmax(i)=p(i);
              }
            }

            return changed;
          }

          bool add(const bounding_box<data_type, dim__> &bb)
          {
            bool rtn1=add(bb.pmin), rtn2=add(bb.pmax);
            return (rtn1 || rtn2);
          }

          bool inside(const point_type &p) const
          {
            for (index_type i=0; i<dim__; ++i)
            {
              if ( (p(i)<pmin(i)) || (p(i)>pmax(i)) )
                return false;
            }

            return true;
          }

          bool intersect(const bounding_box<data_type, dim__> &bb) const
          {
            const index_type nc(1<<dim__);
            point_type c[nc], bb_min(bb.get_min()), bb_max(bb.get_max());

            // set the first and last corner
            c[0]=bb_min;
            c[nc-1]=bb_max;

            // set the corners for 2 or 3 dimensional cases
            if (dim__>1)
            {
              c[1](0)=bb_max(0);
              c[1](1)=bb_min(1);
              c[2](0)=bb_min(0);
              c[2](1)=bb_max(1);
              if (dim__==3)
              {
                c[1](2)=bb_min(2);
                c[2](2)=bb_min(2);
              }
            }

            // set the remaining corners for 3 dimensional cases
            if (dim__>2)
            {
              c[3] << bb_max(0), bb_max(1), bb_min(2);
              c[4] << bb_min(0), bb_min(1), bb_max(2);
              c[5] << bb_max(0), bb_min(1), bb_max(2);
              c[6] << bb_min(0), bb_max(1), bb_max(2);
            }

            // if any corner is inside this bbox, then they intersect
            for (index_type i=0; i<nc; ++i)
            {
              if (inside(c[i]))
                return true;
            }

            // at this point know that no bb edge crosses this bbox, but need
            // to check if this bbox in entirely inside bb
            return bb.inside(pmin);
          }

        private:
          bool empty;
          point_type pmin, pmax;
      };
    }
  }
}
#endif
