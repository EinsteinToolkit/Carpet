#include <cassert>
#include <cmath>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "gh.hh"

#include "carpet.hh"
#include "regrid.hh"

extern "C" {
  static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetRegrid/src/moving.cc,v 1.5 2004/08/02 11:42:20 schnetter Exp $";
  CCTK_FILEVERSION(Carpet_CarpetRegrid_moving_cc);
}



namespace CarpetRegrid {
  
  using namespace std;
  using namespace Carpet;
  
  
  
  int Moving (cGH const * const cctkGH,
              gh<dim> const & hh,
              gh<dim>::rexts  & bbsss,
              gh<dim>::rbnds  & obss,
              gh<dim>::rprocs & pss)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    assert (refinement_levels >= 1);
    
    assert (bbsss.size() >= 1);
    
    bbsss.resize (refinement_levels);
    obss.resize (refinement_levels);
    pss.resize (refinement_levels);
    
    bvect const symmetric (symmetry_x, symmetry_y, symmetry_z);
    bbvect const ob (false);
    
    assert (! smart_outer_boundaries);
    
    for (size_t rl=1; rl<bbsss.size(); ++rl) {
      
      // calculate new extent
      CCTK_REAL const argument = 2*M_PI * moving_circle_frequency * cctk_time;
      rvect const pos
        (moving_centre_x + moving_circle_radius * cos(argument),
         moving_centre_y + moving_circle_radius * sin(argument),
         moving_centre_z);
      CCTK_REAL const radius = moving_region_radius / ipow(reffact, rl-1);
      
      rvect const rlb (symmetric.ifthen (rvect(0),      pos - rvect(radius)));
      rvect const rub (symmetric.ifthen (rvect(radius), pos + rvect(radius)));
      
      vector<ibbox> bbs;
      gh<dim>::cbnds obs;
      
      ManualCoordinates_OneLevel
        (cctkGH, hh, rl, refinement_levels, rlb, rub, ob, bbs, obs);
      
      // make multiprocessor aware
      gh<dim>::cprocs ps;
      SplitRegions (cctkGH, bbs, obs, ps);
      
      // make multigrid aware
      vector<vector<ibbox> > bbss;
      MakeMultigridBoxes (cctkGH, bbs, obs, bbss);
      
      bbsss.at(rl) = bbss;
      obss.at(rl) = obs;
      pss.at(rl) = ps;
      
    } // for rl
    
    return 1;
  }
  
} // namespace CarpetRegrid
