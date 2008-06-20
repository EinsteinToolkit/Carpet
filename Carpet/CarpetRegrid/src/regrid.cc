#include <cassert>
#include <sstream>
#include <string>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "gh.hh"
#include "vect.hh"

#include "carpet.hh"
#include "regrid.hh"



namespace CarpetRegrid {
  
  using namespace std;
  using namespace Carpet;
  
  CCTK_INT CarpetRegrid_Regrid (CCTK_POINTER_TO_CONST const cctkGH_,
                                CCTK_POINTER const superregss_,
                                CCTK_POINTER const regsss_,
				CCTK_INT force)
  {
    DECLARE_CCTK_PARAMETERS;
    
    const cGH * const cctkGH = (const cGH *) cctkGH_;
    
    gh::rregs & superregss = * (gh::rregs *) superregss_;
    gh::mregs & regsss = * (gh::mregs *) regsss_;
    
    gh const & hh = *vhh.at(Carpet::map);
    
    assert (is_singlemap_mode());
    
    // In force mode (force == true) we do not check the CarpetRegrid
    // parameters

    if (!force) {

      assert (regrid_every == -1 or regrid_every == 0
	      or regrid_every % maxmglevelfact == 0);
    
      // Return if no regridding is desired
      if (regrid_every == -1) return 0;
      
      // Return if we want to regrid during initial data only, and this
      // is not the time for initial data
      if (regrid_every == 0 and cctkGH->cctk_iteration != 0) return 0;

      // Return if we want to regrid regularly, but not at this time
      if ((regrid_every > 0 and
           cctkGH->cctk_iteration != 0 and
           (cctkGH->cctk_iteration-1) % regrid_every != 0)
          or
	  ((cctkGH->cctk_iteration-1) == 0 and
           regrid_every > 1))
      {
	return 0;
      }
          
      // Steer parameters
      const int oldnumlevels = refinement_levels;
      if (CCTK_EQUALS(activate_levels_on_regrid, "none")) {
      
	// do nothing
      
      } else if (CCTK_EQUALS(activate_levels_on_regrid, "fixed")) {
      
	if (cctkGH->cctk_iteration-1 >= activate_next) {
	  const int newnumlevels
	    = min((int)(refinement_levels + num_new_levels), maxreflevels);
	  assert (newnumlevels>0 and newnumlevels<=maxreflevels);
        
	  *const_cast<CCTK_INT*>(&activate_next) = cctkGH->cctk_iteration;
	  ostringstream next;
	  next << activate_next;
	  CCTK_ParameterSet
	    ("activate_next", "CarpetRegrid", next.str().c_str());
        
	  *const_cast<CCTK_INT*>(&refinement_levels) = newnumlevels;
	  ostringstream param;
	  param << refinement_levels;
	  CCTK_ParameterSet
	    ("refinement_levels", "CarpetRegrid", param.str().c_str());
        
	  if (verbose) {
	    ostringstream buf1, buf2;
	    buf1 << "Activating " << newnumlevels - oldnumlevels << " new refinement levels";
	    buf2 << "There are now " << newnumlevels << " refinement levels";
	    CCTK_INFO (buf1.str().c_str());
	    CCTK_INFO (buf2.str().c_str());
	  }
	}
      
      } else if (CCTK_EQUALS(activate_levels_on_regrid, "function")) {
      
	if (! CCTK_IsFunctionAliased("RegridLevel")) {
	  CCTK_WARN (0, "No thorn has provided the function \"RegridLevel\"");
	}
	const int newnumlevels
	  = RegridLevel (cctkGH, refinement_levels, maxreflevels);
	if (newnumlevels>0 and newnumlevels<=maxreflevels) {
        
	  *const_cast<CCTK_INT*>(&refinement_levels) = newnumlevels;
	  ostringstream param;
	  param << refinement_levels;
	  CCTK_ParameterSet
	    ("refinement_levels", "CarpetRegrid", param.str().c_str());
        
	  if (verbose) {
	    ostringstream buf1, buf2;
	    buf1 << "Activating " << newnumlevels - oldnumlevels << " new refinement levels";
	    buf2 << "There are now " << newnumlevels << " refinement levels";
	    CCTK_INFO (buf1.str().c_str());
	    CCTK_INFO (buf2.str().c_str());
	  }
        
	} else {
	  CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
		      "The aliased function \"RegridLevel\" returned an illegal number of refinement levels (%d).  No levels will be activated or deactivated.", newnumlevels);
	}
      
      } else {
      
	assert (0);
      
      }
    
    
    
      // Return if this is not during initial data generation, and if no
      // change in the grid structure is desired
      if (cctkGH->cctk_iteration != 0) {
	if (keep_same_grid_structure and refinement_levels == reflevels && !tracking) {
          return 0;
        }
      }
    
    } else {

      // If force is active, steer activate_next to current iteration

      ostringstream next;
      next << cctkGH->cctk_iteration;
      CCTK_ParameterSet
	("activate_next", "CarpetRegrid", next.str().c_str());      

    } // if (!force)    
    
    int do_recompose;
    
    if (CCTK_EQUALS(refined_regions, "none")) {
      
      do_recompose = BaseLevel (cctkGH, hh, superregss);
                 
    } else if (CCTK_EQUALS(refined_regions, "centre")) {
      
      do_recompose = Centre (cctkGH, hh, superregss);
                 
    } else if (CCTK_EQUALS(refined_regions, "manual-gridpoints")) {
      
      do_recompose
        = ManualGridpoints (cctkGH, hh, superregss);
                 
    } else if (CCTK_EQUALS(refined_regions, "manual-coordinates")) {
      
      do_recompose
        = ManualCoordinates (cctkGH, hh, superregss);
                 
    } else if (CCTK_EQUALS(refined_regions, "manual-gridpoint-list")) {
      
      do_recompose
        = ManualGridpointList (cctkGH, hh, superregss);
                 
    } else if (CCTK_EQUALS(refined_regions, "manual-coordinate-list")) {
      
      do_recompose
        = ManualCoordinateList (cctkGH, hh, superregss);
                 
    } else if (CCTK_EQUALS(refined_regions, "moving")) {
      
      do_recompose = Moving (cctkGH, hh, superregss);
                 
    } else if (CCTK_EQUALS(refined_regions, "automatic")) {
      
      do_recompose = Automatic (cctkGH, hh, superregss);
                 
    } else {
      assert (0);
    }
    
    // make multiprocessor aware
    vector<vector<region_t> > regss(superregss.size());
    for (size_t rl=0; rl<superregss.size(); ++rl) {
#warning "TODO: delete .processors"
      SplitRegions (cctkGH, superregss.at(rl), regss.at(rl));
    }
    
    // make multigrid aware
    MakeMultigridBoxes (cctkGH, Carpet::map, regss, regsss);
    
    assert (regsss.size() > 0);
    
    return do_recompose;
  }
  
} // namespace CarpetRegrid
