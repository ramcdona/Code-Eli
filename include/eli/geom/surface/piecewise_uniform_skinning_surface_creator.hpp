/*********************************************************************************
* Copyright (c) 2026 Rob McDonald <rob.a.mcdonald@gmail.com>
*
* All rights reserved. This program and the accompanying materials
* are made available under the terms of the Eclipse Public License v1.0
* which accompanies this distribution, and is available at
* http://www.eclipse.org/legal/epl-v10.html
********************************************************************************/

#ifndef eli_geom_surface_piecewise_uniform_skinning_surface_creator_hpp
#define eli_geom_surface_piecewise_uniform_skinning_surface_creator_hpp

#include <cassert>
#include <vector>

#include "eli/code_eli.hpp"

#include <Eigen/SparseCore>
#include <Eigen/SparseLU>

#include "eli/geom/surface/piecewise_general_skinning_surface_creator.hpp"

namespace eli
{
  namespace geom
  {
    namespace surface
    {
      // Skinning surface creator specialized for the case where every control point strip
      // shares one constraint structure.
      //
      // The general skinning creator solves an independent linear system per control point
      // strip.  With the current rib API the constraint structure -- which conditions
      // exist, the segment degrees, and the coefficient values -- is determined entirely
      // by rib level metadata (continuity, which derivative curves are specified) and the
      // u parameterization, none of which vary between strips.  Only the right hand side
      // values (rib control points) differ.  This creator exploits that: it compiles the
      // condition list once, assembles one sparse N x N coefficient matrix directly (the
      // conditions are identical for every dimension, so the general creator's
      // (N*dim) x (N*dim) block diagonal system is redundant), factorizes it once, and
      // back-substitutes a fresh right hand side per strip.
      //
      // If ribs ever gain per-sector structure (different conditions for different parts
      // of the cross section), that case cannot be represented here -- use
      // piecewise_general_skinning_surface_creator, which remains the reference
      // implementation.  Enable set_validate(true) to check this creator's results
      // against the general implementation at runtime.
      template<typename data__, unsigned short dim__, typename tol__>
      class piecewise_uniform_skinning_surface_creator : public piecewise_general_skinning_surface_creator<data__, dim__, tol__>
      {
        public:
          typedef piecewise_general_skinning_surface_creator<data__, dim__, tol__> base_class_type;
          typedef typename base_class_type::data_type data_type;
          typedef typename base_class_type::point_type point_type;
          typedef typename base_class_type::index_type index_type;
          typedef typename base_class_type::tolerance_type tolerance_type;
          typedef typename base_class_type::piecewise_surface_type piecewise_surface_type;
          typedef typename base_class_type::rib_data_type rib_data_type;

          piecewise_uniform_skinning_surface_creator()
            : base_class_type(), validate(false)
          {
          }
          piecewise_uniform_skinning_surface_creator(const data_type &uu0, const data_type &vv0)
            : base_class_type(uu0, vv0), validate(false)
          {
          }

          void set_validate(bool v) {validate=v;}
          bool get_validate() const {return validate;}

          virtual bool create(piecewise_surface_type &ps) const
          {
            typedef typename piecewise_surface_type::surface_type surface_type;
            typedef eli::geom::curve::piecewise<eli::geom::curve::bezier, data_type, dim__, tolerance_type> piecewise_curve_type;
            typedef typename piecewise_curve_type::curve_type curve_type;

            index_type nribs(this->get_number_u_segments()+1), i, j;
            std::vector<rib_data_type> rib_states(this->ribs);

            // The general creator does not handle closed surfaces either.
            assert(!this->closed);

            ps.clear();

            // split ribs so have same number of curves (with same joint parameters) for
            // all ribs and get degree -- identical to the general creator prologue.
            index_type njoints(this->get_number_v_segments()+1);
            std::vector<data_type> joint_params(njoints);
            std::vector<index_type> max_jdegs(njoints-1, 0);

            joint_params[0]=this->get_v0();
            for (j=0; j<(njoints-1); ++j)
            {
              joint_params[j+1]=joint_params[j]+this->get_segment_dv(j);
            }

            for (i=0; i<nribs; ++i)
            {
              std::vector<index_type> jdegs;
              rib_states[i].split(joint_params.begin(), joint_params.end(), std::back_inserter(jdegs));
              for (j=0; j<(njoints-1); ++j)
              {
                if (jdegs[j]>max_jdegs[j])
                {
                  max_jdegs[j]=jdegs[j];
                }
              }
            }

            for (i=0; i<nribs; ++i)
            {
              rib_states[i].promote(max_jdegs.begin(), max_jdegs.end());
            }

            index_type u, v, nu(nribs-1), nv(njoints-1);

            ps.init_uv(this->du_begin(), this->du_end(), this->dv_begin(), this->dv_end(), this->get_u0(), this->get_v0());

            // ---- structure pass: compile the condition program once ----
            structure_type st;
            if (!compile_structure(rib_states, st))
            {
              return false;
            }

            // assemble and factorize the coefficient matrix once
            Eigen::SparseMatrix<data_type> A(st.nrow, st.nrow);
            A.setFromTriplets(st.triplets.begin(), st.triplets.end());
            Eigen::SparseLU< Eigen::SparseMatrix<data_type> > solver;
            solver.compute(A);
            if (solver.info()!=Eigen::Success)
            {
              assert(false);
              return false;
            }

            // ---- per v segment: gather every strip's right hand side into one wide
            //      matrix and solve them all with a single back substitution, so the
            //      factored solve runs blocked (level 3) across all strips at once ----
            for (v=0; v<nv; ++v)
            {
              std::vector<surface_type> surfs(nu);
              const index_type nstrip(max_jdegs[v]+1);

              // hoist the per-segment curves for each rib out of the strip loop
              std::vector<curve_type> fcrv(nribs), lfpcrv(nribs), rfpcrv(nribs), lfppcrv(nribs), rfppcrv(nribs);
              for (u=0; u<nribs; ++u)
              {
                rib_states[u].get_f().get(fcrv[u], v);
                if (rib_states[u].use_left_fp())
                {
                  rib_states[u].get_left_fp().get(lfpcrv[u], v);
                }
                if (rib_states[u].use_right_fp())
                {
                  rib_states[u].get_right_fp().get(rfpcrv[u], v);
                }
                if (rib_states[u].use_left_fpp())
                {
                  rib_states[u].get_left_fpp().get(lfppcrv[u], v);
                }
                if (rib_states[u].use_right_fpp())
                {
                  rib_states[u].get_right_fpp().get(rfppcrv[u], v);
                }
              }

              // strip j occupies right hand side columns [j*dim, j*dim+dim)
              Eigen::Matrix<data_type, Eigen::Dynamic, Eigen::Dynamic> B(st.nrow, dim__*nstrip);
              for (j=0; j<nstrip; ++j)
              {
                for (index_type c=0; c<st.nrow; ++c)
                {
                  const condition_type &cnd(st.conditions[c]);
                  point_type val;
                  switch (cnd.src)
                  {
                    case SRC_ZERO:
                      val.setZero();
                      break;
                    case SRC_F:
                      val=fcrv[cnd.joint].get_control_point(j);
                      break;
                    case SRC_FP_LEFT:
                      val=lfpcrv[cnd.joint].get_control_point(j);
                      break;
                    case SRC_FP_RIGHT:
                      val=rfpcrv[cnd.joint].get_control_point(j);
                      break;
                    case SRC_FPP_LEFT:
                      val=lfppcrv[cnd.joint].get_control_point(j);
                      break;
                    case SRC_FPP_RIGHT:
                      val=rfppcrv[cnd.joint].get_control_point(j);
                      break;
                    default:
                      assert(false);
                      val.setZero();
                      break;
                  }
                  B.block(c, j*dim__, 1, dim__)=val;
                }
              }

              Eigen::Matrix<data_type, Eigen::Dynamic, Eigen::Dynamic> X(solver.solve(B));
              if (solver.info()!=Eigen::Success)
              {
                assert(false);
                return false;
              }

              // place control points directly into the temporary surfaces
              for (u=0; u<nu; ++u)
              {
                surfs[u].resize(st.seg_degree[u], max_jdegs[v]);
                for (j=0; j<nstrip; ++j)
                {
                  for (i=0; i<=st.seg_degree[u]; ++i)
                  {
                    point_type cp(X.block(st.seg_ind[u]+i, j*dim__, 1, dim__));
                    surfs[u].set_control_point(cp, i, j);
                  }
                }
              }

              // put these surfaces into piecewise surface
              typename piecewise_surface_type::error_code ec;
              for (u=0; u<nu; ++u)
              {
                ec=ps.set(surfs[u], u, v);
                if (ec!=piecewise_surface_type::NO_ERRORS)
                {
                  assert(false);
                  return false;
                }
              }
            }

            if (validate)
            {
              if (!validate_against_general(ps))
              {
                return false;
              }
            }

            return true;
          }

        private:
          bool validate;

          enum condition_src
          {
            SRC_ZERO,
            SRC_F,
            SRC_FP_LEFT,
            SRC_FP_RIGHT,
            SRC_FPP_LEFT,
            SRC_FPP_RIGHT
          };

          struct condition_type
          {
            index_type joint;
            condition_src src;
          };

          struct structure_type
          {
            index_type nrow;
            std::vector<index_type> seg_degree;
            std::vector<index_type> seg_ind;
            std::vector<condition_type> conditions;
            std::vector< Eigen::Triplet<data_type> > triplets;
          };

          // Compile the strip-invariant condition program.  This reproduces, joint for
          // joint, the degree determination and C1/C2 derivative mirroring logic of
          // piecewise_general_creator::create() -- including its exact tie-breaking
          // branches -- but records coefficient scalars and right hand side sources
          // instead of assembling a dense per-strip system.
          bool compile_structure(const std::vector<rib_data_type> &rib_states, structure_type &st) const
          {
            typedef eli::geom::general::continuity continuity_type;

            index_type nsegs(this->get_number_u_segments()), i;
            index_type njnt(nsegs+1);

            // effective flags and value sources per joint
            std::vector<bool> u_lfp(njnt), u_rfp(njnt), u_lfpp(njnt), u_rfpp(njnt);
            std::vector<condition_src> s_lfp(njnt, SRC_FP_LEFT), s_rfp(njnt, SRC_FP_RIGHT);
            std::vector<condition_src> s_lfpp(njnt, SRC_FPP_LEFT), s_rfpp(njnt, SRC_FPP_RIGHT);
            std::vector<continuity_type> cont(njnt);

            for (i=0; i<njnt; ++i)
            {
              u_lfp[i]=rib_states[i].use_left_fp();
              u_rfp[i]=rib_states[i].use_right_fp();
              u_lfpp[i]=rib_states[i].use_left_fpp();
              u_rfpp[i]=rib_states[i].use_right_fpp();
              cont[i]=static_cast<continuity_type>(rib_states[i].get_continuity());
            }

            // open-curve end conditions, as enforced by the general creator
            assert(!u_lfp[0] && !u_lfpp[0] && cont[0]==eli::geom::general::C0);
            assert(!u_rfp[nsegs] && !u_rfpp[nsegs] && cont[nsegs]==eli::geom::general::C0);

            // minimum degree per segment from original flags
            st.seg_degree.assign(nsegs, 1);
            for (i=0; i<nsegs; ++i)
            {
              if (u_rfp[i])   st.seg_degree[i]+=1;
              if (u_lfp[i+1]) st.seg_degree[i]+=1;
              if (u_rfpp[i])  st.seg_degree[i]+=1;
              if (u_lfpp[i+1])st.seg_degree[i]+=1;

              if (!valid_degree(st.seg_degree[i], this->max_degree[i]))
              {
                return false;
              }
            }

            // joint continuity handling: mirror derivative specifications and bump degrees.
            // This replicates the general creator's branches exactly, including the
            // asymmetry between the C1 and C2 "neither side specified" cases.
            for (i=0; i<njnt; ++i)
            {
              if (cont[i]==eli::geom::general::C0)
              {
                continue;
              }

              // C1 handling
              if (u_lfp[i])
              {
                if (!u_rfp[i])
                {
                  u_rfp[i]=true;
                  s_rfp[i]=s_lfp[i];
                  st.seg_degree[i]+=1;
                }
              }
              else
              {
                if (u_rfp[i])
                {
                  u_lfp[i]=true;
                  s_lfp[i]=s_rfp[i];
                  st.seg_degree[i-1]+=1;
                }
                else
                {
                  bool hit_max_im1(this->max_degree[i-1]>0 && st.seg_degree[i-1]>=this->max_degree[i-1]);
                  bool hit_max_i(this->max_degree[i]>0 && st.seg_degree[i]>=this->max_degree[i]);

                  if ((i==0) & (!this->closed))
                  {
                    st.seg_degree[i]+=1;
                  }
                  else if ((i==nsegs) && (!this->closed))
                  {
                    st.seg_degree[i-1]+=1;
                  }
                  else if ( hit_max_i || (!hit_max_im1 && (st.seg_degree[i-1]<=st.seg_degree[i])) )
                  {
                    st.seg_degree[i-1]+=1;
                  }
                  else
                  {
                    st.seg_degree[i]+=1;
                  }
                }
              }

              if (cont[i]==eli::geom::general::C1)
              {
                continue;
              }

              // C2 handling
              if (u_lfpp[i])
              {
                if (!u_rfpp[i])
                {
                  u_rfpp[i]=true;
                  s_rfpp[i]=s_lfpp[i];
                  st.seg_degree[i]+=1;
                }
              }
              else
              {
                if (u_rfpp[i])
                {
                  u_lfpp[i]=true;
                  s_lfpp[i]=s_rfpp[i];
                  st.seg_degree[i-1]+=1;
                }
                else
                {
                  bool hit_max_im1(this->max_degree[i-1]>0 && st.seg_degree[i-1]>=this->max_degree[i-1]);
                  bool hit_max_i(this->max_degree[i]>0 && st.seg_degree[i]>=this->max_degree[i]);

                  if ((i==0) && (!this->closed))
                  {
                    st.seg_degree[i]+=1;
                  }
                  else if ((i==nsegs) && (!this->closed))
                  {
                    st.seg_degree[i-1]+=1;
                  }
                  // note: intentionally a separate if, matching the general creator
                  if ( hit_max_i || (!hit_max_im1 && (st.seg_degree[i-1]<=st.seg_degree[i])) )
                  {
                    st.seg_degree[i-1]+=1;
                  }
                  else
                  {
                    st.seg_degree[i]+=1;
                  }
                }
              }
            }

            // final degree check and control point indexing (no fit points here)
            st.seg_ind.assign(nsegs+1, 0);
            for (i=0; i<nsegs; ++i)
            {
              if (!valid_degree(st.seg_degree[i], this->max_degree[i]))
              {
                return false;
              }
              st.seg_ind[i+1]=st.seg_ind[i]+st.seg_degree[i]+1;
            }
            st.nrow=st.seg_ind[nsegs];

            // ---- emit conditions in the general creator's order ----
            st.conditions.clear();
            st.conditions.reserve(st.nrow);
            st.triplets.clear();
            index_type row(0);

            for (i=0; i<nsegs; ++i)
            {
              const index_type i0(st.seg_ind[i]);
              const index_type deg(st.seg_degree[i]);
              const data_type dt(this->get_segment_du(i));

              // end point conditions
              add_condition(st, row, i, SRC_F);
              st.triplets.push_back(Eigen::Triplet<data_type>(row-1, i0, 1));

              add_condition(st, row, i+1, SRC_F);
              st.triplets.push_back(Eigen::Triplet<data_type>(row-1, i0+deg, 1));

              // first derivative conditions
              if (u_rfp[i])
              {
                const data_type c(deg/dt);
                add_condition(st, row, i, s_rfp[i]);
                st.triplets.push_back(Eigen::Triplet<data_type>(row-1, i0, -c));
                st.triplets.push_back(Eigen::Triplet<data_type>(row-1, i0+1, c));
              }
              if (u_lfp[i+1])
              {
                const data_type c(deg/dt);
                add_condition(st, row, i+1, s_lfp[i+1]);
                st.triplets.push_back(Eigen::Triplet<data_type>(row-1, i0+deg-1, -c));
                st.triplets.push_back(Eigen::Triplet<data_type>(row-1, i0+deg, c));
              }

              // second derivative conditions
              if (u_rfpp[i])
              {
                const data_type c(deg*(deg-1)/dt/dt);
                add_condition(st, row, i, s_rfpp[i]);
                st.triplets.push_back(Eigen::Triplet<data_type>(row-1, i0, c));
                st.triplets.push_back(Eigen::Triplet<data_type>(row-1, i0+1, -2*c));
                st.triplets.push_back(Eigen::Triplet<data_type>(row-1, i0+2, c));
              }
              if (u_lfpp[i+1])
              {
                const data_type c(deg*(deg-1)/dt/dt);
                add_condition(st, row, i+1, s_lfpp[i+1]);
                st.triplets.push_back(Eigen::Triplet<data_type>(row-1, i0+deg-2, c));
                st.triplets.push_back(Eigen::Triplet<data_type>(row-1, i0+deg-1, -2*c));
                st.triplets.push_back(Eigen::Triplet<data_type>(row-1, i0+deg, c));
              }
            }

            // interior joint continuity conditions without specified values
            for (i=0; i<nsegs; ++i)
            {
              if ( (cont[i]>eli::geom::general::C0) && !u_lfp[i] && !u_rfp[i] )
              {
                const index_type lI(st.seg_ind[i-1]);
                const index_type ldeg(st.seg_degree[i-1]);
                const index_type rI(st.seg_ind[i]);
                const index_type rdeg(st.seg_degree[i]);
                const data_type cl(ldeg/this->get_segment_du(i-1));
                const data_type cr(rdeg/this->get_segment_du(i));

                add_condition(st, row, i, SRC_ZERO);
                st.triplets.push_back(Eigen::Triplet<data_type>(row-1, lI+ldeg-1, -cl));
                st.triplets.push_back(Eigen::Triplet<data_type>(row-1, lI+ldeg, cl));
                st.triplets.push_back(Eigen::Triplet<data_type>(row-1, rI, cr));
                st.triplets.push_back(Eigen::Triplet<data_type>(row-1, rI+1, -cr));
              }

              if ( (cont[i]>eli::geom::general::C1) && !u_lfpp[i] && !u_rfpp[i] )
              {
                const index_type lI(st.seg_ind[i-1]);
                const index_type ldeg(st.seg_degree[i-1]);
                const index_type rI(st.seg_ind[i]);
                const index_type rdeg(st.seg_degree[i]);

                add_condition(st, row, i, SRC_ZERO);
                if (ldeg>1)
                {
                  const data_type cl(ldeg*(ldeg-1)/this->get_segment_du(i-1)/this->get_segment_du(i-1));
                  st.triplets.push_back(Eigen::Triplet<data_type>(row-1, lI+ldeg-2, cl));
                  st.triplets.push_back(Eigen::Triplet<data_type>(row-1, lI+ldeg-1, -2*cl));
                  st.triplets.push_back(Eigen::Triplet<data_type>(row-1, lI+ldeg, cl));
                }
                if (rdeg>1)
                {
                  const data_type cr(rdeg*(rdeg-1)/this->get_segment_du(i)/this->get_segment_du(i));
                  st.triplets.push_back(Eigen::Triplet<data_type>(row-1, rI, -cr));
                  st.triplets.push_back(Eigen::Triplet<data_type>(row-1, rI+1, 2*cr));
                  st.triplets.push_back(Eigen::Triplet<data_type>(row-1, rI+2, -cr));
                }
              }
            }

            // the system must be exactly determined, as in the general direct solve
            if (row!=st.nrow)
            {
              assert(false);
              return false;
            }

            return true;
          }

          static void add_condition(structure_type &st, index_type &row, const index_type &joint, condition_src src)
          {
            condition_type c;
            c.joint=joint;
            c.src=src;
            st.conditions.push_back(c);
            ++row;
          }

          static bool valid_degree(const index_type &deg, const index_type &max_deg)
          {
            if (max_deg<=0)
              return true;
            if (deg<=max_deg)
              return true;

            return false;
          }

          // Compare this creator's surface against the reference general implementation.
          bool validate_against_general(const piecewise_surface_type &ps) const
          {
            typedef typename piecewise_surface_type::surface_type surface_type;

            base_class_type gc(*this);
            piecewise_surface_type ps_ref;

            if (!gc.create(ps_ref))
            {
              assert(false);
              return false;
            }

            if ( (ps.number_u_patches()!=ps_ref.number_u_patches()) ||
                 (ps.number_v_patches()!=ps_ref.number_v_patches()) )
            {
              assert(false);
              return false;
            }

            data_type scale(0), err(0);
            for (index_type iu=0; iu<ps.number_u_patches(); ++iu)
            {
              for (index_type iv=0; iv<ps.number_v_patches(); ++iv)
              {
                const surface_type *s(ps.get_patch(iu, iv)), *sr(ps_ref.get_patch(iu, iv));
                if ( (s->degree_u()!=sr->degree_u()) || (s->degree_v()!=sr->degree_v()) )
                {
                  assert(false);
                  return false;
                }
                for (index_type ii=0; ii<=s->degree_u(); ++ii)
                {
                  for (index_type jj=0; jj<=s->degree_v(); ++jj)
                  {
                    point_type d(s->get_control_point(ii, jj)-sr->get_control_point(ii, jj));
                    err=std::max(err, d.norm());
                    scale=std::max(scale, sr->get_control_point(ii, jj).norm());
                  }
                }
              }
            }

            if (scale==0)
            {
              scale=1;
            }
            if (err>1e-8*scale)
            {
              assert(false);
              return false;
            }

            return true;
          }
      };
    }
  }
}
#endif
