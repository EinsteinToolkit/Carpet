#include <assert.h>

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
  
  
  
  int ManualGridpoints (cGH const * const cctkGH,
                        gh const & hh,
                        gh::rregs & regss)
  {
    DECLARE_CCTK_PARAMETERS;
    
    if (refinement_levels > 4) {
      CCTK_WARN (0, "Cannot currently specify manual refinement regions for more than 4 refinement levels");
    }
    assert (refinement_levels >= 1 and refinement_levels <= 4);
    
    // do nothing if the levels already exist
    if (reflevel == refinement_levels) return 0;
    
    regss.resize (refinement_levels);
    
    vector<ivect> ilower(3), iupper(3);
    ilower.at(0) = ivect (l1ixmin, l1iymin, l1izmin);
    iupper.at(0) = ivect (l1ixmax, l1iymax, l1izmax);
    ilower.at(1) = ivect (l2ixmin, l2iymin, l2izmin);
    iupper.at(1) = ivect (l2ixmax, l2iymax, l2izmax);
    ilower.at(2) = ivect (l3ixmin, l3iymin, l3izmin);
    iupper.at(2) = ivect (l3ixmax, l3iymax, l3izmax);
    
    assert (not smart_outer_boundaries);
    
    for (size_t rl=1; rl<regss.size(); ++rl) {
      
      b2vect const ob (false);
      region_t reg;
      reg.map = Carpet::map;
      reg.outer_boundaries = ob;
      
      vector<region_t> regs;
      ManualGridpoints_OneLevel
        (cctkGH, hh, rl,refinement_levels,
         ilower.at(rl-1), iupper.at(rl-1), reg, regs);
      
      regss.at(rl) = regs;
      
    } // for rl
    
    return 1;
  }
  
  
  
  void ManualGridpoints_OneLevel (const cGH * const cctkGH,
                                  const gh & hh,
                                  const int rl,
                                  const int numrl,
                                  const ivect ilower,
                                  const ivect iupper,
                                  const region_t & reg,
                                  vector<region_t> & regs)
  {
    const ivect rstr = hh.baseextents.at(0).at(0).stride();
    const ivect rlb  = hh.baseextents.at(0).at(0).lower();
    const ivect rub  = hh.baseextents.at(0).at(0).upper();
    
    const ivect levfac = hh.reffacts.at(rl);
    assert (all (rstr % levfac == 0));
    const ivect str (rstr / levfac);
    const ivect lb  (ilower);
    const ivect ub  (iupper);
    if (! all(lb>=rlb and ub<=rub)) {
      ostringstream buf;
      buf << "The refinement region boundaries for refinement level #" << rl << " are not within the main grid.  Allowed are the grid point boundaries " << rlb << " - " << rub << "; specified were " << lb << " - " << ub << ends;
      CCTK_WARN (0, buf.str().c_str());
    }
    if (! all(lb<=ub)) {
      ostringstream buf;
      buf << "The refinement region boundaries for refinement level #" << rl << " have the upper boundary (" << ub << ") less than the lower boundary (" << lb << ")" << ends;
      CCTK_WARN (0, buf.str().c_str());
    }
    if (! all(lb%str==0 and ub%str==0)) {
      CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
		  "The refinement region boundaries for refinement level #%d are not a multiple of the stride for that level", rl);
    }
    assert (all(lb>=rlb and ub<=rub));
    assert (all(lb<=ub));
    assert (all(lb%str==0 and ub%str==0));
    
    region_t newreg (reg);
    newreg.extent = ibbox(lb, ub, str);
    regs.push_back (newreg);
  }
  
} // namespace CarpetRegrid
