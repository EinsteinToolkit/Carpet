#include <assert.h>

#include <sstream>
#include <string>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "gh.hh"
#include "vect.hh"

#include "carpet.hh"
#include "regrid.hh"

extern "C" {
  static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetRegrid/src/regrid.cc,v 1.40 2004/04/19 18:48:07 schnetter Exp $";
  CCTK_FILEVERSION(Carpet_CarpetRegrid_regrid_cc);
}



namespace CarpetRegrid {
  
  using namespace std;
  using namespace Carpet;
  
  CCTK_INT CarpetRegrid_Regrid (CCTK_POINTER_TO_CONST const cctkGH_,
                                CCTK_POINTER const bbsss_,
                                CCTK_POINTER const obss_,
                                CCTK_POINTER const pss_)
  {
    DECLARE_CCTK_PARAMETERS;
    
    const cGH * const cctkGH = (const cGH *) cctkGH_;
    
    gh<dim>::rexts  & bbsss = * (gh<dim>::rexts  *) bbsss_;
    gh<dim>::rbnds  & obss  = * (gh<dim>::rbnds  *) obss_;
    gh<dim>::rprocs & pss   = * (gh<dim>::rprocs *) pss_;
    
    gh<dim> const & hh = *vhh.at(Carpet::map);
    
    assert (is_singlemap_mode());
    
    
    
    assert (regrid_every == -1 || regrid_every == 0
	    || regrid_every % maxmglevelfact == 0);
    
    // Return if no regridding is desired
    if (regrid_every == -1) return 0;
    
    // Return if we want to regrid during initial data only, and this
    // is not the time for initial data
    if (regrid_every == 0 && cctkGH->cctk_iteration != 0) return 0;

    // Return if we want to regrid regularly, but not at this time
    if (regrid_every > 0 && cctkGH->cctk_iteration != 0
	&& cctkGH->cctk_iteration % regrid_every != 0) {
      return 0;
    }
    
    
    
    // Steer parameters
    if (CCTK_EQUALS(activate_levels_on_regrid, "none")) {
      
      // do nothing
      
    } else if (CCTK_EQUALS(activate_levels_on_regrid, "fixed")) {
      
      if (cctkGH->cctk_iteration >= activate_next) {
        const int newnumlevels
          = min(refinement_levels + num_new_levels, maxreflevels);
        assert (newnumlevels>0 && newnumlevels<=maxreflevels);
        
        *const_cast<CCTK_INT*>(&activate_next) = cctkGH->cctk_iteration + 1;
        ostringstream next;
        next << activate_next;
        CCTK_ParameterSet
          ("activate_next", "CarpetRegrid", next.str().c_str());
        
        *const_cast<CCTK_INT*>(&refinement_levels) = newnumlevels;
        ostringstream param;
        param << refinement_levels;
        CCTK_ParameterSet
          ("refinement_levels", "CarpetRegrid", param.str().c_str());
        
      }
      
    } else if (CCTK_EQUALS(activate_levels_on_regrid, "function")) {
      
      if (! CCTK_IsFunctionAliased("RegridLevel")) {
        CCTK_WARN (0, "No thorn has provided the function \"RegridLevel\"");
      }
      const int newnumlevels
        = RegridLevel (cctkGH, refinement_levels, maxreflevels);
      if (newnumlevels>0 && newnumlevels<=maxreflevels) {
        
        *const_cast<CCTK_INT*>(&refinement_levels) = newnumlevels;
        ostringstream param;
        param << refinement_levels;
        CCTK_ParameterSet
          ("refinement_levels", "CarpetRegrid", param.str().c_str());
        
      } else {
        CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "The aliased function \"RegridLevel\" returned an illegal number of refinement levels (%d).  No levels will be activated or deactivated.", newnumlevels);
      }
      
    } else {
      
      assert (0);
      
    }
    
    
    
    int do_recompose;
    
    if (CCTK_EQUALS(refined_regions, "none")) {
      
      do_recompose = BaseLevel (cctkGH, hh, bbsss, obss, pss);
                 
    } else if (CCTK_EQUALS(refined_regions, "centre")) {
      
      do_recompose = Centre (cctkGH, hh, bbsss, obss, pss);
                 
    } else if (CCTK_EQUALS(refined_regions, "manual-gridpoints")) {
      
      do_recompose
        = ManualGridpoints (cctkGH, hh, bbsss, obss, pss);
                 
    } else if (CCTK_EQUALS(refined_regions, "manual-coordinates")) {
      
      do_recompose
        = ManualCoordinates (cctkGH, hh, bbsss, obss, pss);
                 
    } else if (CCTK_EQUALS(refined_regions, "manual-gridpoint-list")) {
      
      do_recompose
        = ManualGridpointList (cctkGH, hh, bbsss, obss, pss);
                 
    } else if (CCTK_EQUALS(refined_regions, "manual-coordinate-list")) {
      
      do_recompose
        = ManualCoordinateList (cctkGH, hh, bbsss, obss, pss);
                 
    } else if (CCTK_EQUALS(refined_regions, "moving")) {
      
      do_recompose = Moving (cctkGH, hh, bbsss, obss, pss);
                 
    } else if (CCTK_EQUALS(refined_regions, "automatic")) {
      
      do_recompose = Automatic (cctkGH, hh, bbsss, obss, pss);
                 
    } else {
      assert (0);
    }
    
    return do_recompose;
  }
  
} // namespace CarpetRegrid
