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
  
  
  
  int ManualGridpointList (cGH const * const cctkGH,
                           gh const & hh,
                           gh::mexts  & bbsss,
                           gh::rbnds  & obss,
                           gh::rprocs & pss)
  {
    DECLARE_CCTK_PARAMETERS;
    
    assert (refinement_levels >= 1);
    
    // do nothing if the levels already exist
    if (reflevel == refinement_levels) return 0;
    
    assert (bbsss.size() >= 1);
    vector<vector<ibbox> > bbss = bbsss.at(0);
    
    bbss.resize (refinement_levels);
    obss.resize (refinement_levels);
    pss.resize (refinement_levels);
    
    vector<vector<ibbox> > newbbss;
    if (strcmp(gridpoints, "") != 0) {
      istringstream gp_str(gridpoints);
      try {
        gp_str >> newbbss;
      } catch (input_error) {
        CCTK_WARN (0, "Could not parse parameter \"gridpoints\"");
      }
    }
    
    vector<vector<bbvect> > newobss;
    if (strcmp(outerbounds, "") !=0 ) {
      istringstream ob_str (outerbounds);
      try {
        ob_str >> newobss;
      } catch (input_error) {
        CCTK_WARN (0, "Could not parse parameter \"outerbounds\"");
      }
      bool good = newobss.size() == newbbss.size();
      if (good) {
        for (size_t rl=0; rl<newobss.size(); ++rl) {
          good = good && newobss.at(rl).size() == newbbss.at(rl).size();
        }
      }
      if (! good) {
        cout << "gridpoints: " << newbbss << endl;
        cout << "outerbounds: " << newobss << endl;
        CCTK_WARN (0, "The parameters \"outerbounds\" and \"gridpoints\" must have the same structure");
      }
    } else {
      newobss.resize(newbbss.size());
      for (size_t rl=0; rl<newobss.size(); ++rl) {
        newobss.at(rl).resize(newbbss.at(rl).size());
        for (size_t c=0; c<newobss.at(rl).size(); ++c) {
          newobss.at(rl).at(c) = bbvect(false);
        }
      }
    }
    
    if (newbbss.size() < refinement_levels-1) {
      CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "The parameter \"gridpoints\" must contain at least \"refinement_levels-1\" (here: %d) levels", (int)refinement_levels-1);
    }
    
    for (size_t rl=1; rl<refinement_levels; ++rl) {
      
      vector<ibbox> bbs;
      gh::cbnds obs;
      
      bbs.reserve (newbbss.at(rl-1).size());
      obs.reserve (newbbss.at(rl-1).size());
      
      for (size_t c=0; c<newbbss.at(rl-1).size(); ++c) {
        ibbox const & ext = newbbss.at(rl-1).at(c);
        bbvect const & ob = newobss.at(rl-1).at(c);
        ManualGridpoints_OneLevel
          (cctkGH, hh, rl, refinement_levels,
           ext.lower(), ext.upper(), ob, bbs, obs);
      }
      
      // make multiprocessor aware
      gh::cprocs ps;
      SplitRegions (cctkGH, bbs, obs, ps);
      
      bbss.at(rl) = bbs;
      obss.at(rl) = obs;
      pss.at(rl) = ps;
      
      // make multigrid aware
      MakeMultigridBoxes (cctkGH, bbss, obss, bbsss);
      
    } // for rl
    
    return 1;
  }
  
} // namespace CarpetRegrid
