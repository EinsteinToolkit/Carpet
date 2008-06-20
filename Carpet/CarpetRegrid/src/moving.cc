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
              gh::rregs & regss)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    assert (refinement_levels >= 1);
    
    regss.resize (refinement_levels);
    
    bvect const symmetric (symmetry_x, symmetry_y, symmetry_z);
    b2vect const ob (false);
    b2vect const rb (true);
    
    assert (! smart_outer_boundaries);
    
    for (size_t rl=1; rl<regss.size(); ++rl) {
      
      // calculate new extent
      CCTK_REAL const argument = 2*M_PI * moving_circle_frequency * cctk_time;
      rvect const pos
        (moving_centre_x + moving_circle_radius * cos(argument),
         moving_centre_y + moving_circle_radius * sin(argument),
         moving_centre_z);
      rvect const radius
        (rvect(moving_region_radius) / rvect(spacereffacts.at(rl-1)));
      
      rvect const rlb (either (symmetric, rvect(0), pos - radius));
      rvect const rub (either (symmetric, radius  , pos + radius));
      
      region_t reg;
      reg.map = Carpet::map;
      reg.outer_boundaries = b2vect (false);
      
      vector<region_t> regs;
      ManualCoordinates_OneLevel
        (cctkGH, hh, rl, refinement_levels, rlb, rub, reg, regs);
      
      regss.at(rl) = regs;
      
    } // for rl
    
    return 1;
  }
  
} // namespace CarpetRegrid
