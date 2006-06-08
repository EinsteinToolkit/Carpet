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
    void CarpetAdaptiveRegridParamcheck (CCTK_ARGUMENTS);
    
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

  int ManualCoordinateList (cGH const * const cctkGH,
                            gh const & hh,
                            gh::mexts  & bbsss,
                            gh::rbnds  & obss,
                            gh::rprocs & pss,
                            gh::mexts  & local_bbsss,
                            gh::rbnds  & local_obss);

    
  void ManualCoordinates_OneLevel (const cGH * const cctkGH,
                                   const gh & hh,
                                   const int rl,
                                   const int numrl,
                                   const rvect lower,
                                   const rvect upper,
                                   const bbvect obound,
                                   vector<ibbox> & bbs,
                                   vector<bbvect> & obs);
  
  void ManualGridpoints_OneLevel (const cGH * const cctkGH,
                                  const gh & hh,
                                  const int rl,
                                  const int numrl,
                                  const ivect ilower,
                                  const ivect iupper,
                                  const bbvect obound,
                                  vector<ibbox> & bbs,
                                  vector<bbvect> & obs);

  rvect int2pos (const cGH* const cctkGH, const gh& hh,
                 const ivect & ipos, const int rl);

  ivect pos2int (const cGH* const cctkGH, const gh& hh,
                 const rvect & rpos, const int rl);

} // namespace CarpetAdaptiveregrid

#endif // !defined(CARPETADAPTIVEREGRID_HH)
