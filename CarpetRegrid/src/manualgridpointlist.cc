#include <assert.h>
#include <string.h>

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

int ManualGridpointList(cGH const *const cctkGH, gh const &hh,
                        gh::rregs &regss) {
  DECLARE_CCTK_PARAMETERS;

  assert(refinement_levels >= 1);

  // do nothing if the levels already exist
  if (reflevel == refinement_levels)
    return 0;

  regss.resize(refinement_levels);

  vector<vector<ibbox> > newbbss;
  if (strcmp(gridpoints, "") != 0) {
    istringstream gp_str(gridpoints);
    try {
      gp_str >> newbbss;
    } catch (input_error) {
      CCTK_WARN(0, "Could not parse parameter \"gridpoints\"");
    }
  }

  vector<vector<bbvect> > newobss;
  if (strcmp(outerbounds, "") != 0) {
    istringstream ob_str(outerbounds);
    try {
      ob_str >> newobss;
    } catch (input_error) {
      CCTK_WARN(0, "Could not parse parameter \"outerbounds\"");
    }
    bool good = newobss.size() == newbbss.size();
    if (good) {
      for (size_t rl = 0; rl < newobss.size(); ++rl) {
        good = good and newobss.at(rl).size() == newbbss.at(rl).size();
      }
    }
    if (!good) {
      cout << "gridpoints: " << newbbss << endl;
      cout << "outerbounds: " << newobss << endl;
      CCTK_WARN(0, "The parameters \"outerbounds\" and \"gridpoints\" must "
                   "have the same structure");
    }
  } else {
    newobss.resize(newbbss.size());
    for (size_t rl = 0; rl < newobss.size(); ++rl) {
      newobss.at(rl).resize(newbbss.at(rl).size());

      for (size_t c = 0; c < newobss.at(rl).size(); ++c) {
        newobss.at(rl).at(c) = bbvect(vect<bool, 2>(false));
      }
    }
  }

  vector<vector<bbvect> > newrbss;
  newrbss.resize(newobss.size());
  for (int rl = 0; rl < (int)newobss.size(); ++rl) {
    newrbss.at(rl).resize(newbbss.at(rl).size());
    for (int c = 0; c < (int)newobss.at(rl).size(); ++c) {
      bbvect const &ob = newobss.at(rl).at(c);
      bbvect &rb = newrbss.at(rl).at(c);
      for (int d = 0; d < dim; ++d) {
        for (int f = 0; f < 2; ++f) {
          rb[d][f] = !ob[d][f];
        }
      }
    }
  }

  if ((int)newbbss.size() < refinement_levels - 1) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "The parameter \"gridpoints\" must contain at least "
               "\"refinement_levels-1\" (here: %d) levels",
               (int)refinement_levels - 1);
  }

  vector<vector<region_t> > newregs(newbbss.size());
  for (int rl = 1; rl < refinement_levels; ++rl) {

    vector<region_t> regs;

    regs.reserve(newbbss.at(rl - 1).size());

    for (size_t c = 0; c < newbbss.at(rl - 1).size(); ++c) {
      ibbox const ext = newbbss.at(rl - 1).at(c);
      b2vect const ob = xpose(newobss.at(rl - 1).at(c));
      region_t reg;
      reg.extent = ext;
      reg.map = Carpet::map;
      reg.outer_boundaries = ob;
      ManualGridpoints_OneLevel(cctkGH, hh, rl, refinement_levels, ext.lower(),
                                ext.upper(), reg, regs);
    }

    regss.at(rl) = regs;

  } // for rl

  return 1;
}

} // namespace CarpetRegrid
