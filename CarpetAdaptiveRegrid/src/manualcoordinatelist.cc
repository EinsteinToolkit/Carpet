#include <cassert>
#include <cmath>
#include <cstring>
#include <sstream>
#include <vector>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "gh.hh"

#include "CAR.hh"
#include "carpet.hh"

namespace CarpetAdaptiveRegrid {

using namespace std;
using namespace Carpet;

int ManualCoordinateList(cGH const *const cctkGH, gh const &hh,
                         gh::mexts &bbsss, gh::rbnds &obss, gh::rbnds &rbss,
                         gh::rprocs &pss, gh::mexts &local_bbsss,
                         gh::rbnds &local_obss, gh::rbnds &local_rbss) {
  DECLARE_CCTK_PARAMETERS;
  int ierr;

  assert(refinement_levels >= 1);

  // do nothing if the levels already exist
  if (reflevel == refinement_levels)
    return 0;

  assert(bbsss.size() >= 1);
  vector<vector<ibbox> > bbss = bbsss.at(0);
  vector<vector<ibbox> > local_bbss = local_bbsss.at(0);

  jjvect nboundaryzones, is_internal, is_staggered, shiftout;
  ierr = GetBoundarySpecification(2 * dim, &nboundaryzones[0][0],
                                  &is_internal[0][0], &is_staggered[0][0],
                                  &shiftout[0][0]);
  assert(!ierr);
  rvect physical_min, physical_max;
  rvect interior_min, interior_max;
  rvect exterior_min, exterior_max;
  rvect base_spacing;
  ierr = GetDomainSpecification(
      dim, &physical_min[0], &physical_max[0], &interior_min[0],
      &interior_max[0], &exterior_min[0], &exterior_max[0], &base_spacing[0]);
  assert(!ierr);

  bbss.resize(refinement_levels);
  local_bbss.resize(refinement_levels);
  obss.resize(refinement_levels);
  local_obss.resize(refinement_levels);
  rbss.resize(refinement_levels);
  local_rbss.resize(refinement_levels);
  pss.resize(refinement_levels);

  vector<vector<rbbox> > newbbss;
  if (strcmp(coordinates, "") != 0) {
    istringstream gp_str(coordinates);
    try {
      gp_str >> newbbss;
    } catch (input_error) {
      CCTK_WARN(0, "Could not parse parameter \"coordinates\"");
    }
  }

  vector<vector<bbvect> > newobss;
  vector<vector<bbvect> > newrbss;

  newobss.resize(newbbss.size());
  newrbss.resize(newbbss.size());
  for (size_t rl = 0; rl < newobss.size(); ++rl) {
    newobss.at(rl).resize(newbbss.at(rl).size());
    newrbss.at(rl).resize(newbbss.at(rl).size());
    for (size_t c = 0; c < newobss.at(rl).size(); ++c) {
      for (int d = 0; d < dim; ++d) {
        assert(mglevel == 0);
        rvect const spacing = base_spacing *
                              ipow((CCTK_REAL)mgfact, basemglevel) /
                              rvect(spacereffacts.at(rl + 1));
        ierr = ConvertFromPhysicalBoundary(
            dim, &physical_min[0], &physical_max[0], &interior_min[0],
            &interior_max[0], &exterior_min[0], &exterior_max[0], &spacing[0]);
        assert(!ierr);
        newobss.at(rl).at(c)[d][0] = abs(newbbss.at(rl).at(c).lower()[d] -
                                         physical_min[d]) < 1.0e-6 * spacing[d];
        newrbss.at(rl).at(c)[d][0] = !newobss.at(rl).at(c)[d][0];
        if (newobss.at(rl).at(c)[d][0]) {
          rvect lo = newbbss.at(rl).at(c).lower();
          rvect up = newbbss.at(rl).at(c).upper();
          rvect str = newbbss.at(rl).at(c).stride();
          lo[d] = exterior_min[d];
          newbbss.at(rl).at(c) = rbbox(lo, up, str);
        }
        newobss.at(rl).at(c)[d][1] =
            abs(newbbss.at(rl).at(c).upper()[d] - physical_max[d]) <
            1.0e-6 * base_spacing[d] / spacereffacts.at(rl)[d];
        newrbss.at(rl).at(c)[d][1] = !newobss.at(rl).at(c)[d][1];
        if (newobss.at(rl).at(c)[d][1]) {
          rvect lo = newbbss.at(rl).at(c).lower();
          rvect up = newbbss.at(rl).at(c).upper();
          rvect str = newbbss.at(rl).at(c).stride();
          up[d] = exterior_max[d];
          newbbss.at(rl).at(c) = rbbox(lo, up, str);
        }
      }
    }
  }

  if (newbbss.size() < refinement_levels - 1) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "The parameter \"coordinates\" must contain at least "
               "\"refinement_levels-1\" (here: %d) levels",
               int(refinement_levels - 1));
  }

  for (size_t rl = 1; rl < refinement_levels; ++rl) {

    vector<ibbox> bbs;
    gh::cbnds obs;
    gh::cbnds rbs;

    bbs.reserve(newbbss.at(rl - 1).size());
    obs.reserve(newbbss.at(rl - 1).size());
    rbs.reserve(newbbss.at(rl - 1).size());

    for (size_t c = 0; c < newbbss.at(rl - 1).size(); ++c) {
      rbbox const &ext = newbbss.at(rl - 1).at(c);
      bbvect const &ob = newobss.at(rl - 1).at(c);
      bbvect const &rb = newrbss.at(rl - 1).at(c);
      // TODO: why can basemglevel not be used here?
      // rvect const spacing = base_spacing * ipow(CCTK_REAL(mgfact),
      // basemglevel) / ipow(reffact, rl);
      rvect const spacing = base_spacing / rvect(spacereffacts.at(rl));
      if (!all(abs(ext.stride() - spacing) < spacing * (CCTK_REAL)1.0e-10)) {
        assert(dim == 3);
        CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "The grid spacing on refinement level %d is incorrect.  I "
                   "expected [%g,%g,%g], but I found [%g,%g,%g].",
                   int(rl), double(spacing[0]), double(spacing[1]),
                   double(spacing[2]), double(ext.stride()[0]),
                   double(ext.stride()[1]), double(ext.stride()[2]));
      }
      assert(all(abs(ext.stride() - spacing) < spacing * (CCTK_REAL)1.0e-10));

      ManualCoordinates_OneLevel(cctkGH, hh, rl, refinement_levels, ext.lower(),
                                 ext.upper(), ob, rb, bbs, obs, rbs);
    }

    local_bbss.at(rl) = bbs;
    local_obss.at(rl) = obs;
    local_rbss.at(rl) = rbs;

    if (verbose) {
      ostringstream buf;
      buf << "About to make multigrid aware for first time whilst doing "
          << "manual coordinate list, level " << rl
          << ". Total list is:" << endl
          << bbss << endl
          << "with obss being " << local_obss;
      CCTK_INFO(buf.str().c_str());
    }

    // make multigrid aware
    MakeMultigridBoxes(cctkGH, local_bbss, local_obss, local_bbsss);

    if (verbose) {
      ostringstream buf;
      buf << "Doing manual coordinate list, level " << rl
          << ". Total list is:" << endl
          << local_bbsss << endl
          << "with obss being " << local_obss;
      CCTK_INFO(buf.str().c_str());
    }

    // make multiprocessor aware
    gh::cprocs ps;
    SplitRegions(cctkGH, bbs, obs, rbs, ps);

    bbss.at(rl) = bbs;
    obss.at(rl) = obs;
    rbss.at(rl) = rbs;
    pss.at(rl) = ps;

    if (verbose) {
      ostringstream buf;
      buf << "Doing manual coordinate list, multiprocs, level " << rl << "."
          << endl
          << "Total list is:" << endl
          << bbss << endl
          << "with obss being " << obss << endl
          << "Sizes are (bb) " << bbss.at(rl).size() << " and (ob) "
          << obss.at(rl).size();
      CCTK_INFO(buf.str().c_str());
    }

    // make multigrid aware
    MakeMultigridBoxes(cctkGH, bbss, obss, bbsss);

    if (verbose) {
      ostringstream buf;
      buf << "Did manual coordinate list, multiprocs, level " << rl << "."
          << endl
          << "Total list is:" << endl
          << bbsss;
      CCTK_INFO(buf.str().c_str());
    }

  } // for rl

  return 1;
}

void ManualCoordinates_OneLevel(const cGH *const cctkGH, const gh &hh,
                                const int rl, const int numrl,
                                const rvect lower, const rvect upper,
                                const bbvect obound, const bbvect rbound,
                                vector<ibbox> &bbs, vector<bbvect> &obs,
                                vector<bbvect> &rbs) {
  if (rl >= numrl)
    return;

  jvect const ilower = pos2int(cctkGH, hh, lower, rl);
  jvect const iupper = pos2int(cctkGH, hh, upper, rl);

  ManualGridpoints_OneLevel(cctkGH, hh, rl, numrl, ilower, iupper, obound,
                            rbound, bbs, obs, rbs);
}

void ManualGridpoints_OneLevel(const cGH *const cctkGH, const gh &hh,
                               const int rl, const int numrl,
                               const ivect ilower, const ivect iupper,
                               const bbvect obound, const bbvect rbound,
                               vector<ibbox> &bbs, vector<bbvect> &obs,
                               vector<bbvect> &rbs) {
  const ivect rstr = hh.baseextent.stride();
  const ivect rlb = hh.baseextent.lower();
  const ivect rub = hh.baseextent.upper();

  const ivect levfac = hh.reffacts.at(rl);
  assert(all(rstr % levfac == 0));
  const ivect str(rstr / levfac);
  const ivect lb(ilower);
  const ivect ub(iupper);
  if (!all(lb >= rlb && ub <= rub)) {
    ostringstream buf;
    buf << "The refinement region boundaries for refinement level #" << rl
        << " are not within the main grid.  Allowed are the grid point "
           "boundaries "
        << rlb << " - " << rub << "; specified were " << lb << " - " << ub
        << ends;
    CCTK_WARN(0, buf.str().c_str());
  }
  if (!all(lb <= ub)) {
    ostringstream buf;
    buf << "The refinement region boundaries for refinement level #" << rl
        << " have the upper boundary (" << ub
        << ") less than the lower boundary (" << lb << ")" << ends;
    CCTK_WARN(0, buf.str().c_str());
  }
  if (!all(lb % str == 0 && ub % str == 0)) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "The refinement region boundaries for refinement level #%d are "
               "not a multiple of the stride for that level",
               rl);
  }
  assert(all(lb >= rlb && ub <= rub));
  assert(all(lb <= ub));
  assert(all(lb % str == 0 && ub % str == 0));

  bbs.push_back(ibbox(lb, ub, str));
  obs.push_back(obound);
  rbs.push_back(rbound);
}

} // namespace CarpetRegrid
