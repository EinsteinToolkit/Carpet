#include <cassert>
#include <cmath>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "gh.hh"

#include "carpet.hh"
#include "regrid.hh"



namespace CarpetRegrid {
  
  using namespace std;
  using namespace Carpet;
  
  
  
  int Moving (cGH const * const cctkGH,
              gh const & hh,
              gh::mexts  & bbsss,
              gh::rbnds  & obss,
              gh::rprocs & pss)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    assert (refinement_levels >= 1);
    
    assert (bbsss.size() >= 1);
    vector<vector<ibbox> > bbss = bbsss.at(0);
    
    bbss.resize (refinement_levels);
    obss.resize (refinement_levels);
    pss.resize (refinement_levels);
    
    bvect const symmetric (symmetry_x, symmetry_y, symmetry_z);
    bbvect const ob (false);
    
    assert (! smart_outer_boundaries);
    
    for (size_t rl=1; rl<bbss.size(); ++rl) {
      
      // calculate new extent
      CCTK_REAL const argument = 2*M_PI * moving_circle_frequency * cctk_time;
      rvect const pos
        (moving_centre_x + moving_circle_radius * cos(argument),
         moving_centre_y + moving_circle_radius * sin(argument),
         moving_centre_z);
      rvect const radius
        (rvect(moving_region_radius) / spacereffacts.at(rl-1));
      
      rvect const rlb (symmetric.ifthen (rvect(0), pos - radius));
      rvect const rub (symmetric.ifthen (radius  , pos + radius));
      
      vector<ibbox> bbs;
      gh::cbnds obs;
      
      ManualCoordinates_OneLevel
        (cctkGH, hh, rl, refinement_levels, rlb, rub, ob, bbs, obs);
      
      // make multiprocessor aware
      gh::cprocs ps;
      SplitRegions (cctkGH, bbs, obs, ps);
      
      bbss.at(rl) = bbs;
      obss.at(rl) = obs;
      pss.at(rl) = ps;
      
    } // for rl
    
    // make multigrid aware
    MakeMultigridBoxes (cctkGH, bbss, obss, bbsss);
    
    return 1;
  }
  
} // namespace CarpetRegrid
