// $Header:$

#ifndef CARPETADAPTIVEREGRID_HH
#define CARPETADAPTIVEREGRID_HH

#include <list>

#include "cctk.h"
#include "cctk_Arguments.h"

#include "bbox.hh"
#include "gf.hh"
#include "gh.hh"
#include "vect.hh"

#include "carpet.hh"



namespace CarpetAdaptiveRegrid {
  
  using namespace std;
  using namespace Carpet;
  
  
  
  extern "C" {
    
    /* Scheduled functions */
    int CarpetAdaptiveRegridParamcheck (CCTK_ARGUMENTS);
    
    /* Aliased functions */
//     CCTK_INT CarpetAdaptiveRegrid_Regrid (const cGH * const cctkGH,
//                                   gh<dim>::rexts  * bbsss,
//                                   gh<dim>::rbnds  * obss,
//                                   gh<dim>::rprocs * pss);
    CCTK_INT CarpetAdaptiveRegrid_Regrid (CCTK_POINTER_TO_CONST const cctkGH_,
                                  CCTK_POINTER const bbsss_,
                                  CCTK_POINTER const obss_,
                                  CCTK_POINTER const pss_,
				  CCTK_INT force);
  }

} // namespace CarpetAdaptiveregrid

#endif // !defined(CARPETADAPTIVEREGRID_HH)
