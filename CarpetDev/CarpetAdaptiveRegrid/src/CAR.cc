#include <assert.h>

#include <sstream>
#include <string>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "gh.hh"
#include "vect.hh"

#include "carpet.hh"
#include "CAR.hh"

extern "C" {
  static const char* rcsid = "$Header:$";
  CCTK_FILEVERSION(Carpet_CarpetAdaptiveregrid_regrid_cc);
}



namespace CarpetAdaptiveRegrid {
  
  using namespace std;
  using namespace Carpet;

  extern "C" {
    void CCTK_FCALL CCTK_FNAME(copy_mask)
      (const CCTK_INT& snx, const CCTK_INT& sny, const CCTK_INT& snz,
       const CCTK_INT* smask, const CCTK_INT sbbox[3][3],
       const CCTK_INT& dnx, const CCTK_INT& dny, const CCTK_INT& dnz,
       CCTK_INT* dmask, const CCTK_INT dbbox[3][3]);
    void CCTK_FCALL CCTK_FNAME(check_box)
      (const CCTK_INT& nx, const CCTK_INT& ny, const CCTK_INT& nz,
       const CCTK_INT* mask,
       CCTK_INT* sum_x, CCTK_INT* sum_y, CCTK_INT* sum_z,
       CCTK_INT* sig_x, CCTK_INT* sig_y, CCTK_INT* sig_z,
       const CCTK_INT bbox[3][3],
       CCTK_INT newbbox1[3][3], CCTK_INT newbbox2[3][3],
       const CCTK_INT& min_width, const CCTK_REAL& min_density,
       CCTK_INT& didit);
  }
  
  CCTK_INT CarpetAdaptiveRegrid_Regrid (CCTK_POINTER_TO_CONST const cctkGH_,
                                CCTK_POINTER const bbsss_,
                                CCTK_POINTER const obss_,
                                CCTK_POINTER const pss_,
				CCTK_INT force)
  {
    DECLARE_CCTK_PARAMETERS;
    
    const cGH * const cctkGH = (const cGH *) cctkGH_;
    
    gh<dim>::rexts  & bbsss = * (gh<dim>::rexts  *) bbsss_;
    gh<dim>::rbnds  & obss  = * (gh<dim>::rbnds  *) obss_;
    gh<dim>::rprocs & pss   = * (gh<dim>::rprocs *) pss_;
    
    gh<dim> const & hh = *vhh.at(Carpet::map);
    
    assert (is_singlemap_mode());
    
    // In force mode (force == true) we do not check the
    // CarpetAdaptiveregrid parameters

    if (!force) {

      assert (regrid_every == -1 || regrid_every == 0
	      || regrid_every % maxmglevelfact == 0);
    
      // Return if no regridding is desired
      if (regrid_every == -1) return 0;
      
      // Return if we want to regrid during initial data only, and this
      // is not the time for initial data
      if (regrid_every == 0 && cctkGH->cctk_iteration != 0) return 0;

      // Return if we want to regrid regularly, but not at this time
      if (regrid_every > 0 && cctkGH->cctk_iteration != 0
	  && (cctkGH->cctk_iteration-1) % regrid_every != 0)
      {
	return 0;
      }

    }    

    int do_recompose;

    // Make a cboxlist (or set) from the boxes on the level.

    // cboxlist is list of cbox's.

    // cbox is a bbox 
    // plus a dim-D array of booleans mask(i,j,k)
    // plus 3 sum arrays (sigma_i = sum_{j,k} mask(i,:,:) etc)
    // plus 3 signature arrays (1D laplacian of sigma arrays)

    // First thing to do; make a cbox containing all active boxes
    // on this level.
    // cbox list contains this cbox.
    // Create the mask for this cbox by looking at the error GF.

    // Then recursively for every cbox in the list until no changes:

    //   Prune the cbox (remove all zeros from the ends)

    //   Split the cbox (find zero in sigma_:, split along that line)

    //   if density too low then find zero crossings of signature and
    //   split there

    // if none of these steps are taken the cbox is finished.

    // Most of the actual work is done by the utility functions in
    // CAR_utils.F77.
    // What we have to do is

    // loop over every box in the list
    //   for each box call check_box
    //   if it returns zero the box is done
    //   if it returns one then the box needs shrinking
    //   if it returns two then the box needs splitting

    // so, the method will be:

    // Create a cboxlist with one member

    // loop over the list:
    //   call check_box
    //   if the return is two then create a new box from newbbox2 and add it
    //     to the end of the list
    //   if the return value is one or two shrink the current box to the
    //     values given in newbbox1
    //   if the return value is one or two go redo this loop again.
    
    return do_recompose;
  }
  
} // namespace CarpetAdaptiveRegrid
