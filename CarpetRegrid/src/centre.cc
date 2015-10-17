#include <assert.h>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "gh.hh"

#include "carpet.hh"
#include "regrid.hh"



namespace CarpetRegrid {
  
  using namespace std;
  using namespace Carpet;
  
  
  
  int Centre (cGH const * const cctkGH,
              gh const & hh,
              gh::rregs & regss)
  {
    DECLARE_CCTK_PARAMETERS;
    
    assert (refinement_levels >= 1);
    
    // do nothing if the levels already exist
    if (reflevel == refinement_levels) return 0;
    
    regss.resize (refinement_levels);
    
    bvect const symmetric (symmetry_x, symmetry_y, symmetry_z);
    ivect const zero(0), one(1), two(2);
    
    ivect rstr = hh.baseextents.at(0).at(0).stride();
    ivect rlb  = hh.baseextents.at(0).at(0).lower();
    ivect rub  = hh.baseextents.at(0).at(0).upper();
    
    assert (not smart_outer_boundaries);
    
    for (size_t rl=1; rl<regss.size(); ++rl) {
      
      // save old values
      ivect const oldrlb = rlb;
      ivect const oldrub = rub;
      
      // refined boxes have smaller stride
      assert (all(rstr%(hh.reffacts.at(rl)/hh.reffacts.at(rl-1)) == 0));
      rstr /= hh.reffacts.at(rl)/hh.reffacts.at(rl-1);
      
      // calculate new extent
      ivect const quarter = (rub - rlb) / 4 / rstr * rstr;
      ivect const half    = (rub - rlb) / 2 / rstr * rstr;
      rlb = oldrlb + either (symmetric, zero, quarter);
      rub = oldrub - either (symmetric, half, quarter);
      assert (all(rlb >= oldrlb and rub <= oldrub));
      
      vector<region_t> regs (1);
      
      ibbox const ext (rlb, rub, rstr);
      regs.at(0).extent = ext;
      
      b2vect const ob (false);
      regs.at(0).outer_boundaries = ob;
      
      regs.at(0).map = Carpet::map;
      
      regss.at(rl) = regs;
      
    } // for rl
    
    return 1;
  }
  
} // namespace CarpetRegrid
