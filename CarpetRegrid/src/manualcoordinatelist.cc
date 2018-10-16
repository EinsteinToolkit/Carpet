#include <cassert>
#include <cmath>
#include <cstring>
#include <sstream>
#include <vector>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "gh.hh"

#include "carpet.hh"
#include "regrid.hh"

namespace CarpetRegrid {

using namespace std;
using namespace Carpet;

int ManualCoordinateList(cGH const *const cctkGH, gh const &hh,
                         gh::rregs &regss) {
  DECLARE_CCTK_PARAMETERS;
  int ierr;

  assert(refinement_levels >= 1);

  // do nothing if the levels already exist
  if (reflevel == refinement_levels && !tracking)
    return 0;

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

  regss.resize(refinement_levels);

  vector<vector<rbbox> > newbbss;
  if (strcmp(coordinates, "") != 0) {
    istringstream gp_str(coordinates);
    try {
      gp_str >> newbbss;
    } catch (input_error) {
      CCTK_WARN(0, "Could not parse parameter \"coordinates\"");
    }
    if (newbbss.size() >= spacereffacts.size()) {
      CCTK_WARN(0, "Parameter \"coordinates\" defines too many refinement "
                   "levels; at most Carpet::max_refinement_levels - 1 may be "
                   "defined");
    }
  }

  if ((int)newbbss.size() < refinement_levels - 1) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "The parameter \"coordinates\" must contain at least "
               "\"refinement_levels-1\" (here: %d) levels",
               int(refinement_levels - 1));
  }

  // Remove superfluous boxes
  newbbss.resize(refinement_levels - 1);

  vector<vector<b2vect> > newobss;
  if (smart_outer_boundaries) {
    // TODO:
    // assert (domain_from_coordbase);

    newobss.resize(newbbss.size());
    for (int rl = 0; rl < (int)newobss.size(); ++rl) {
      ivect const spacereffact = spacereffacts.at(rl + 1);
      assert(mglevel == 0);
      rvect const spacing = base_spacing *
                            ipow((CCTK_REAL)mgfact, basemglevel) /
                            rvect(spacereffact);
      ierr = ConvertFromPhysicalBoundary(
          dim, &physical_min[0], &physical_max[0], &interior_min[0],
          &interior_max[0], &exterior_min[0], &exterior_max[0], &spacing[0]);
      assert(!ierr);

      newobss.at(rl).resize(newbbss.at(rl).size());
      for (int c = 0; c < (int)newobss.at(rl).size(); ++c) {
        rvect lo = newbbss.at(rl).at(c).lower();
        rvect up = newbbss.at(rl).at(c).upper();
        rvect str = newbbss.at(rl).at(c).stride();
        b2vect ob;
        for (int d = 0; d < dim; ++d) {
          ob[0][d] = (abs(newbbss.at(rl).at(c).lower()[d] - physical_min[d]) <
                      1.0e-6 * spacing[d]);
          if (ob[0][d]) {
            lo[d] = exterior_min[d];
          }
          ob[1][d] = (abs(newbbss.at(rl).at(c).upper()[d] - physical_max[d]) <
                      1.0e-6 * base_spacing[d] / spacereffact[d]);
          if (ob[1][d]) {
            up[d] = exterior_max[d];
          }
          str[d] *= ipow((CCTK_REAL)mgfact, basemglevel);
        }
        newbbss.at(rl).at(c) = rbbox(lo, up, str);
        newobss.at(rl).at(c) = ob;
      } // for c
    }   // for rl

  } else { // if not smart_outer_boundaries

    vector<vector<bbvect> > newobss1;
    if (strcmp(outerbounds, "") != 0) {
      istringstream ob_str(outerbounds);
      try {
        ob_str >> newobss1;
      } catch (input_error) {
        CCTK_WARN(0, "Could not parse parameter \"outerbounds\"");
      }
      if (newobss1.size() >= spacereffacts.size()) {
        CCTK_WARN(0, "Parameter \"outerbounds\" defines too many refinement "
                     "levels; at most Carpet::max_refinement_levels - 1 may be "
                     "defined");
      }
      bool good = newobss1.size() == newbbss.size();
      if (good) {
        for (int rl = 0; rl < (int)newobss1.size(); ++rl) {
          good = good and newobss1.at(rl).size() == newbbss.at(rl).size();
        }
      }
      if (!good) {
        cout << "coordinates: " << newbbss << endl;
        cout << "outerbounds: " << newobss1 << endl;
        CCTK_WARN(0, "The parameters \"outerbounds\" and \"coordinates\" must "
                     "have the same structure");
      }

      newobss.resize(newobss1.size());
      for (int rl = 0; rl < (int)newobss.size(); ++rl) {
        newobss.at(rl).resize(newobss1.at(rl).size());
        for (int c = 0; c < (int)newobss.at(rl).size(); ++c) {
          newobss.at(rl).at(c) = xpose(newobss1.at(rl).at(c));
        }
      }

    } else {
      newobss.resize(newbbss.size());
      for (int rl = 0; rl < (int)newobss.size(); ++rl) {
        newobss.at(rl).resize(newbbss.at(rl).size());
        for (int c = 0; c < (int)newobss.at(rl).size(); ++c) {
          newobss.at(rl).at(c) = b2vect(bvect(false));
        }
      }
    }

  } // if not smart_outer_boundaries

  for (int rl = 1; rl < refinement_levels; ++rl) {

    vector<region_t> regs;
    regs.reserve(newbbss.at(rl - 1).size());

    for (int c = 0; c < (int)newbbss.at(rl - 1).size(); ++c) {
      rbbox const &ext = newbbss.at(rl - 1).at(c);
      b2vect const &ob = newobss.at(rl - 1).at(c);
      // TODO:
      // assert (domain_from_coordbase);
      ivect const spacereffact = spacereffacts.at(rl);
      rvect const spacing = base_spacing *
                            ipow(CCTK_REAL(mgfact), basemglevel) /
                            rvect(spacereffact);
      if (!all(abs(ext.stride() - spacing) < spacing * CCTK_REAL(1.0e-10))) {
        assert(dim == 3);
        CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "The grid spacing on refinement level %d is incorrect.  I "
                   "expected [%g,%g,%g], but I found [%g,%g,%g].",
                   int(rl), double(spacing[0]), double(spacing[1]),
                   double(spacing[2]), double(ext.stride()[0]),
                   double(ext.stride()[1]), double(ext.stride()[2]));
      }
      assert(all(abs(ext.stride() - spacing) < spacing * CCTK_REAL(1.0e-10)));

      rvect offset = rvect(0);
      if (c < num_offsets) {
        if (rl >= offset_firstlevel) {
          assert(dim == 3);
          offset = rvect(offsetx[c], offsety[c], offsetz[c]);
        }
      }

      region_t reg;
      reg.map = Carpet::map;
      reg.outer_boundaries = ob;

      ManualCoordinates_OneLevel(cctkGH, hh, rl, refinement_levels,
                                 ext.lower() + offset, ext.upper() + offset,
                                 reg, regs);
    }

    if (merge_overlapping_components) {

    // Check if one or more of our components touch
    again:
      // Loop over all pairs of components
      for (int c = 0; c < (int)regs.size(); ++c) {
        for (int cc = 0; cc < c; ++cc) {
          ibbox const overlap = regs.at(c).extent & regs.at(cc).extent;
          if (not overlap.empty()) {
            ibbox const combined =
                regs.at(c).extent.expanded_containing(regs.at(cc).extent);
            regs.at(c).extent = combined;
            assert(all(all(regs.at(c).outer_boundaries ==
                           regs.at(cc).outer_boundaries)));
            regs.erase(regs.begin() + cc);
            goto again;
          }
        }
      }

    } // if merge_overlapping_components

    regss.at(rl) = regs;

  } // for rl

  return 1;
}

} // namespace CarpetRegrid
