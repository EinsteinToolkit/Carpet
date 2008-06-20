#include <cassert>
#include <vector>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "gh.hh"

#include "carpet.hh"
#include "regrid.hh"



namespace CarpetRegrid {
  
  using namespace std;
  using namespace Carpet;
  
  
  
  int ManualCoordinates (cGH const * const cctkGH,
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
    
    vector<rvect> lower(3), upper(3);
    lower.at(0) = rvect (l1xmin, l1ymin, l1zmin);
    upper.at(0) = rvect (l1xmax, l1ymax, l1zmax);
    lower.at(1) = rvect (l2xmin, l2ymin, l2zmin);
    upper.at(1) = rvect (l2xmax, l2ymax, l2zmax);
    lower.at(2) = rvect (l3xmin, l3ymin, l3zmin);
    upper.at(2) = rvect (l3xmax, l3ymax, l3zmax);
    
    assert (! smart_outer_boundaries);
    
    for (size_t rl=1; rl<regss.size(); ++rl) {
      
      b2vect const ob (false);
      region_t reg;
      reg.map = Carpet::map;
      reg.outer_boundaries = ob;
      
      vector<region_t> regs;
      ManualCoordinates_OneLevel
        (cctkGH, hh, rl, refinement_levels,
         lower.at(rl-1), upper.at(rl-1), reg, regs);
      
      regss.at(rl) = regs;
      
    } // for rl
    
    return 1;
  }
  
  
  
  void ManualCoordinates_OneLevel (const cGH * const cctkGH,
                                   const gh & hh,
                                   const int rl,
                                   const int numrl,
                                   const rvect lower,
                                   const rvect upper,
                                   const region_t & reg,
                                   vector<region_t> & regs)
  {
    if (rl >= numrl) return;
    
    jvect const ilower = pos2int (cctkGH, hh, lower, rl);
    jvect const iupper = pos2int (cctkGH, hh, upper, rl);
    
    ManualGridpoints_OneLevel
      (cctkGH, hh, rl, numrl, ilower, iupper, reg, regs);
  }
  
  
  
  ivect delta2int (const cGH * const cctkGH, const gh& hh,
                   const rvect & rpos, const int rl)
  {
    rvect global_lower, global_upper;
    for (int d=0; d<dim; ++d) {
      const int ierr = CCTK_CoordRange
	(cctkGH, &global_lower[d], &global_upper[d], d+1, 0, "cart3d");
      if (ierr<0) {
	global_lower[d] = 0;
	global_upper[d] = 1;
      }
    }
    const ivect global_extent (hh.baseextents.at(0).at(0).upper() -
                               hh.baseextents.at(0).at(0).lower());
    
    const rvect scale  = rvect(global_extent) / (global_upper - global_lower);
    const ivect levfac = hh.reffacts.at(rl);
    assert (all (hh.baseextents.at(0).at(0).stride() % levfac == 0));
    const ivect istride = hh.baseextents.at(0).at(0).stride() / levfac;
    
    const ivect ipos
      = ivect(floor(rpos * scale / rvect(istride) + (CCTK_REAL) 0.5)) * istride;
    
    const rvect apos = rpos * scale;
    assert (all(abs(apos - rvect(ipos)) < rvect(istride) * (CCTK_REAL) 0.01));
    
    return ipos;
  }
  
  
  
  ivect pos2int (const cGH* const cctkGH, const gh& hh,
                 const rvect & rpos, const int rl)
  {
    rvect global_lower, global_upper;
#if 0
    for (int d=0; d<dim; ++d) {
      const int ierr = CCTK_CoordRange
	(cctkGH, &global_lower[d], &global_upper[d], d+1, 0, "cart3d");
      if (ierr<0) {
	global_lower[d] = 0;
	global_upper[d] = 1;
      }
    }
#endif
    assert (Carpet::map >= 0);
    global_lower = domainspecs.at(Carpet::map).exterior_min;
    global_upper = domainspecs.at(Carpet::map).exterior_max;
    const ivect global_extent (hh.baseextents.at(0).at(0).upper() -
                               hh.baseextents.at(0).at(0).lower());
    
    const rvect scale  = rvect(global_extent) / (global_upper - global_lower);
    const ivect levfac = hh.reffacts.at(rl);
    assert (all (hh.baseextents.at(0).at(0).stride() % levfac == 0));
    const ivect istride = hh.baseextents.at(0).at(0).stride() / levfac;
    
    const ivect ipos
      = (ivect(floor((rpos - global_lower) * scale / rvect(istride)
                     + (CCTK_REAL) 0.5))
         * istride);
    
    return ipos;
  }
  
} // namespace CarpetRegrid
