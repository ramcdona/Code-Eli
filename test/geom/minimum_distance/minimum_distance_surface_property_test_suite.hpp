/*********************************************************************************
* Property based tests for the seeded surface minimum distance search.
*
* The existing minimum distance tests compare against coordinates written out to
* many digits, and assert that every starting guess returns the same u, v and
* distance.  Both make the tests brittle in ways that have nothing to do with
* correctness.  A hard coded coordinate breaks whenever the iteration path
* changes, even when the answer is just as good; and requiring every seed to
* agree is only true while every seed succeeds -- improve the search for some
* seeds and not others, and previously agreeing failures start to disagree.
*
* These tests assert properties instead.  Every one of them is a statement that
* has to hold of any correct answer, whichever of several equally good answers
* is returned, and none of them contains a number that came out of a previous
* run:
*
*   - a seeded search never returns something worse than the seed it was given;
*   - seeded from the best cell of a coarse sweep, it does at least as well as
*     that sweep;
*   - the reported distance really is the distance to the reported point;
*   - at an interior answer the residual is perpendicular to the surface;
*   - a seed that is not already at a stationary point gets moved closer.
*
* The last of those is the one with teeth: it is the property that was being
* violated, and it fails on the search as it stood before this series.  "No
* worse than the seed" is what makes a seeded search worth calling at all, but
* it is already guaranteed by a guard at the end of the iteration, so it holds
* on the unfixed code too and is here to keep that guard honest rather than to
* catch the original defect.
*
* Every one of them calls the seeded form.  Reaching for the unseeded overload
* here is an easy mistake and a silent one -- it tests a different solver, and
* would pass with this whole series reverted.
*********************************************************************************/

#ifndef minimum_distance_surface_property_test_suite_hpp
#define minimum_distance_surface_property_test_suite_hpp

#include <cmath>
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <limits>

#include "eli/constants/math.hpp"
#include "eli/geom/surface/bezier.hpp"
#include "eli/geom/intersect/minimum_distance_surface.hpp"

template<typename data__>
class minimum_distance_surface_property_test_suite : public Test::Suite
{
  public:
    typedef data__ data_type;
    typedef eli::geom::surface::bezier<data_type, 3> surface_type;
    typedef typename surface_type::point_type point_type;
    typedef typename surface_type::index_type index_type;
    typedef typename surface_type::tolerance_type tolerance_type;

    tolerance_type tol;

  protected:
    void AddTests(const float &)
    {
      TEST_ADD(minimum_distance_surface_property_test_suite<float>::seeded_never_worse_than_seed_test);
      TEST_ADD(minimum_distance_surface_property_test_suite<float>::seeded_beats_coarse_sweep_test);
      TEST_ADD(minimum_distance_surface_property_test_suite<float>::answer_is_self_consistent_test);
      TEST_ADD(minimum_distance_surface_property_test_suite<float>::interior_answer_is_orthogonal_test);
      TEST_ADD(minimum_distance_surface_property_test_suite<float>::seeded_makes_progress_test);
    }
    void AddTests(const double &)
    {
      TEST_ADD(minimum_distance_surface_property_test_suite<double>::seeded_never_worse_than_seed_test);
      TEST_ADD(minimum_distance_surface_property_test_suite<double>::seeded_beats_coarse_sweep_test);
      TEST_ADD(minimum_distance_surface_property_test_suite<double>::answer_is_self_consistent_test);
      TEST_ADD(minimum_distance_surface_property_test_suite<double>::interior_answer_is_orthogonal_test);
      TEST_ADD(minimum_distance_surface_property_test_suite<double>::seeded_makes_progress_test);
    }
    void AddTests(const long double &)
    {
      TEST_ADD(minimum_distance_surface_property_test_suite<long double>::seeded_never_worse_than_seed_test);
      TEST_ADD(minimum_distance_surface_property_test_suite<long double>::seeded_beats_coarse_sweep_test);
      TEST_ADD(minimum_distance_surface_property_test_suite<long double>::answer_is_self_consistent_test);
      TEST_ADD(minimum_distance_surface_property_test_suite<long double>::interior_answer_is_orthogonal_test);
      TEST_ADD(minimum_distance_surface_property_test_suite<long double>::seeded_makes_progress_test);
    }

  public:
    minimum_distance_surface_property_test_suite()
    {
      AddTests(data__());
    }
    ~minimum_distance_surface_property_test_suite()
    {
    }

  private:
    // A bicubic patch with enough shape that the tangent plane is a poor guide far from the
    // answer, which is the situation the seeded search has to cope with.
    void build_curved_patch(surface_type &s) const
    {
      index_type n(3), m(3), i, j;
      point_type cp[4][4];

      cp[0][0] << -15,   0,  15;
      cp[1][0] <<  -5,   5,  15;
      cp[2][0] <<   5,  -5,  15;
      cp[3][0] <<  15,   0,  15;
      cp[0][1] << -15,   5,   5;
      cp[1][1] <<  -5,  10,   5;
      cp[2][1] <<   5,  -8,   5;
      cp[3][1] <<  15,   5,   5;
      cp[0][2] << -15,  -5,  -5;
      cp[1][2] <<  -5,   9,  -5;
      cp[2][2] <<   5,  -9,  -5;
      cp[3][2] <<  15,  -5,  -5;
      cp[0][3] << -15,   0, -15;
      cp[1][3] <<  -5,   6, -15;
      cp[2][3] <<   5,  -6, -15;
      cp[3][3] <<  15,   0, -15;

      s.resize(n, m);
      for (i=0; i<=n; ++i)
      {
        for (j=0; j<=m; ++j)
        {
          s.set_control_point(cp[i][j], i, j);
        }
      }
    }

    // Several points, off the surface in different directions and at different distances.
    void build_probe_points(std::vector<point_type> &pts) const
    {
      point_type p;

      p <<   0,  30,   0;  pts.push_back(p);
      p <<   0, -30,   0;  pts.push_back(p);
      p <<  20,  10,  10;  pts.push_back(p);
      p << -20,  10, -10;  pts.push_back(p);
      p <<   0,   2,   0;  pts.push_back(p);
      p <<  40,   0,   0;  pts.push_back(p);
      p <<   3,  -4,   7;  pts.push_back(p);
    }

    data_type distance_at(const surface_type &s, const data_type &u, const data_type &v, const point_type &pt) const
    {
      return eli::geom::point::distance(s.f(u, v), pt);
    }

  public:
    // The property that makes a seeded search worth calling: it must never hand back something
    // worse than the guess it was given.  A search that cannot improve on its guess should return
    // the guess, not wander off.
    void seeded_never_worse_than_seed_test()
    {
      surface_type s;
      std::vector<point_type> pts;
      data_type umin, umax, vmin, vmax;

      build_curved_patch(s);
      build_probe_points(pts);

      s.get_parameter_min(umin, vmin);
      s.get_parameter_max(umax, vmax);

      const index_type n(7);

      for (size_t ipt=0; ipt<pts.size(); ++ipt)
      {
        for (index_type i=0; i<n; ++i)
        {
          for (index_type j=0; j<n; ++j)
          {
            data_type u0(umin + (umax-umin)*(i+static_cast<data_type>(0.5))/n);
            data_type v0(vmin + (vmax-vmin)*(j+static_cast<data_type>(0.5))/n);
            data_type u, v;

            data_type dseed = distance_at(s, u0, v0, pts[ipt]);
            data_type dist = eli::geom::intersect::minimum_distance(u, v, s, pts[ipt], u0, v0);

            // Allow for round off, but nothing more than that.
            data_type slack = 100*std::numeric_limits<data_type>::epsilon()*(1+std::abs(dseed));

            TEST_ASSERT(dist <= dseed + slack);
            TEST_ASSERT(u >= umin);
            TEST_ASSERT(u <= umax);
            TEST_ASSERT(v >= vmin);
            TEST_ASSERT(v <= vmax);
          }
        }
      }
    }

    // A seeded search should at least find what a coarse sweep of the surface would have found
    // from somewhere nearby.  Stated against the best of the sweep rather than against the true
    // minimum, so a genuine second local minimum does not make the test lie.
    void seeded_beats_coarse_sweep_test()
    {
      surface_type s;
      std::vector<point_type> pts;
      data_type umin, umax, vmin, vmax;

      build_curved_patch(s);
      build_probe_points(pts);

      s.get_parameter_min(umin, vmin);
      s.get_parameter_max(umax, vmax);

      const index_type nsweep(40);

      for (size_t ipt=0; ipt<pts.size(); ++ipt)
      {
        // Best the surface has to offer, to sweep resolution, and where it was found.
        data_type dsweep(std::numeric_limits<data_type>::max());
        data_type usweep(umin), vsweep(vmin);
        for (index_type i=0; i<=nsweep; ++i)
        {
          for (index_type j=0; j<=nsweep; ++j)
          {
            data_type uu(umin + (umax-umin)*i/nsweep);
            data_type vv(vmin + (vmax-vmin)*j/nsweep);
            data_type dd = distance_at(s, uu, vv, pts[ipt]);
            if (dd < dsweep)
            {
              dsweep = dd;
              usweep = uu;
              vsweep = vv;
            }
          }
        }

        // Seeded from the best cell of that sweep, which is how a seeded search is actually
        // used -- sweep coarsely, then refine.  Calling the unseeded form here would exercise
        // the unseeded solver instead, and would pass with the tangent plane work reverted.
        data_type u, v;
        data_type dist = eli::geom::intersect::minimum_distance(u, v, s, pts[ipt], usweep, vsweep);

        TEST_ASSERT(dist <= dsweep + std::sqrt(std::numeric_limits<data_type>::epsilon())*(1+dsweep));
      }
    }

    // Whatever it returns, the distance it reports has to be the distance to the point it
    // reports.  Cheap, and it catches a whole class of bookkeeping mistake -- returning the
    // distance from one iterate and the coordinates of another.
    void answer_is_self_consistent_test()
    {
      surface_type s;
      std::vector<point_type> pts;
      data_type umin, umax, vmin, vmax;

      build_curved_patch(s);
      build_probe_points(pts);

      s.get_parameter_min(umin, vmin);
      s.get_parameter_max(umax, vmax);

      const index_type n(5);

      for (size_t ipt=0; ipt<pts.size(); ++ipt)
      {
        for (index_type i=0; i<n; ++i)
        {
          for (index_type j=0; j<n; ++j)
          {
            data_type u0(umin + (umax-umin)*(i+static_cast<data_type>(0.5))/n);
            data_type v0(vmin + (vmax-vmin)*(j+static_cast<data_type>(0.5))/n);
            data_type u, v;

            data_type dist = eli::geom::intersect::minimum_distance(u, v, s, pts[ipt], u0, v0);

            TEST_ASSERT(tol.approximately_equal(dist, distance_at(s, u, v, pts[ipt])));
          }
        }
      }
    }

    // The one the others miss.  "Never worse than the seed" is satisfied by a search that gives up
    // and hands the seed straight back, which is exactly the failure this suite was written to
    // chase: on a body of revolution there were starting guesses from which the seeded search
    // returned the guess untouched while a sweep found a point less than half as far away.
    //
    // So state the stronger thing: if the seed is not already a stationary point -- if the residual
    // there is not perpendicular to the surface -- then the search has to make some progress.  A
    // search that cannot improve on a guess that is plainly not the answer has not done its job.
    void seeded_makes_progress_test()
    {
      surface_type s;
      std::vector<point_type> pts;
      data_type umin, umax, vmin, vmax;

      build_curved_patch(s);
      build_probe_points(pts);

      s.get_parameter_min(umin, vmin);
      s.get_parameter_max(umax, vmax);

      const index_type n(7);
      data_type edge = std::sqrt(std::numeric_limits<data_type>::epsilon());

      for (size_t ipt=0; ipt<pts.size(); ++ipt)
      {
        for (index_type i=0; i<n; ++i)
        {
          for (index_type j=0; j<n; ++j)
          {
            data_type u0(umin + (umax-umin)*(i+static_cast<data_type>(0.5))/n);
            data_type v0(vmin + (vmax-vmin)*(j+static_cast<data_type>(0.5))/n);

            point_type r  = s.f(u0, v0) - pts[ipt];
            point_type su = s.f_u(u0, v0);
            point_type sv = s.f_v(u0, v0);

            data_type rmag(r.norm()), sumag(su.norm()), svmag(sv.norm());

            if ((rmag <= edge) || (sumag <= 0) || (svmag <= 0))
            {
              continue;
            }

            data_type cosu = std::abs(su.dot(r))/(sumag*rmag);
            data_type cosv = std::abs(sv.dot(r))/(svmag*rmag);

            // Only demand progress where there is obviously progress to be made.  The threshold is
            // deliberately generous: this is not measuring how well the search converges, it is
            // catching a search that does not start.
            data_type notstationary = static_cast<data_type>(0.01);

            if ((cosu > notstationary) || (cosv > notstationary))
            {
              data_type u, v;
              data_type dist = eli::geom::intersect::minimum_distance(u, v, s, pts[ipt], u0, v0);

              TEST_ASSERT(dist < rmag);
            }
          }
        }
      }
    }

    // At an answer that is not against a parameter boundary, the residual must be perpendicular to
    // the surface -- that is the condition being solved for.  Expressed as a cosine so it is
    // dimensionless and one tolerance covers any model scale.
    void interior_answer_is_orthogonal_test()
    {
      surface_type s;
      std::vector<point_type> pts;
      data_type umin, umax, vmin, vmax;

      build_curved_patch(s);
      build_probe_points(pts);

      s.get_parameter_min(umin, vmin);
      s.get_parameter_max(umax, vmax);

      data_type edge = std::sqrt(std::numeric_limits<data_type>::epsilon());

      const index_type nseed(4);

      for (size_t ipt=0; ipt<pts.size(); ++ipt)
      {
        // Seeded, and from a spread of seeds rather than one.  The property has to hold for
        // whatever the seeded search returns, not only for the unseeded solver's answer -- and
        // calling the unseeded form here would not touch the tangent plane iteration at all.
        for (index_type iseed=0; iseed<=nseed; ++iseed)
        {
          data_type u0(umin + (umax-umin)*iseed/nseed);
          data_type v0(vmin + (vmax-vmin)*iseed/nseed);

          data_type u, v;
          data_type dist = eli::geom::intersect::minimum_distance(u, v, s, pts[ipt], u0, v0);

          bool interior = (u > umin+edge) && (u < umax-edge) && (v > vmin+edge) && (v < vmax-edge);

          if (interior && dist > edge)
          {
            point_type r  = s.f(u, v) - pts[ipt];
            point_type su = s.f_u(u, v);
            point_type sv = s.f_v(u, v);

            data_type rmag(r.norm()), sumag(su.norm()), svmag(sv.norm());

            if ((sumag > 0) && (svmag > 0) && (rmag > 0))
            {
              data_type cosu = std::abs(su.dot(r))/(sumag*rmag);
              data_type cosv = std::abs(sv.dot(r))/(svmag*rmag);

              // Loose on purpose.  This is checking that the answer is a stationary point at
              // all, not how many digits the solver polished it to.
              data_type costol = 100*std::sqrt(std::numeric_limits<data_type>::epsilon());

              TEST_ASSERT(cosu < costol);
              TEST_ASSERT(cosv < costol);
            }
          }
        }
      }
    }
};

#endif
