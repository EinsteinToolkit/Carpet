#ifndef CARPETREGRID_HH
#define CARPETREGRID_HH

#include <list>
#include <vector>

#include "cctk.h"
#include "cctk_Arguments.h"

#include "bbox.hh"
#include "gf.hh"
#include "gh.hh"
#include "vect.hh"

#include "carpet.hh"



namespace CarpetRegrid {
  
  using namespace std;
  using namespace Carpet;
  
  
  
  extern "C" {
    
    /* Scheduled functions */
    int CarpetRegridParamcheck (CCTK_ARGUMENTS);
    
    /* Aliased functions */
//     CCTK_INT CarpetRegrid_Regrid (const cGH * const cctkGH,
//                                   gh<dim>::mexts  * bbsss,
//                                   gh<dim>::rbnds  * obss,
//                                   gh<dim>::rprocs * pss);
    CCTK_INT CarpetRegrid_Regrid (CCTK_POINTER_TO_CONST const cctkGH_,
                                  CCTK_POINTER const bbsss_,
                                  CCTK_POINTER const obss_,
                                  CCTK_POINTER const pss_,
				  CCTK_INT force);
  }
  
  
  
  int BaseLevel (cGH const * const cctkGH,
                 gh const & hh,
                 gh::mexts  & bbsss,
                 gh::rbnds  & obss,
                 gh::rprocs & pss);

  int Centre (cGH const * const cctkGH,
              gh const & hh,
              gh::mexts  & bbsss,
              gh::rbnds  & obss,
              gh::rprocs & pss);

  int ManualGridpoints (cGH const * const cctkGH,
                        gh const & hh,
                        gh::mexts  & bbsss,
                        gh::rbnds  & obss,
                        gh::rprocs & pss);
  
  void ManualGridpoints_OneLevel (const cGH * const cctkGH,
                                  const gh & hh,
                                  const int rl,
                                  const int numrl,
                                  const ivect ilower,
                                  const ivect iupper,
                                  const bbvect obound,
                                  vector<ibbox> & bbs,
                                  vector<bbvect> & obs);

  int ManualCoordinates (cGH const * const cctkGH,
                         gh const & hh,
                         gh::mexts  & bbsss,
                         gh::rbnds  & obss,
                         gh::rprocs & pss);
  
  void ManualCoordinates_OneLevel (const cGH * const cctkGH,
                                   const gh & hh,
                                   const int rl,
                                   const int numrl,
                                   const rvect lower,
                                   const rvect upper,
                                   const bbvect obound,
                                   vector<ibbox> & bbs,
                                   vector<bbvect> & obs);
  
  ivect delta2int (const cGH * const cctkGH,
                   const gh& hh,
                   const rvect & rpos,
                   const int rl);
  ivect pos2int (const cGH* const cctkGH,
                 const gh& hh,
                 const rvect & rpos,
                 const int rl);

  int ManualGridpointList (cGH const * const cctkGH,
                           gh const & hh,
                           gh::mexts  & bbsss,
                           gh::rbnds  & obss,
                           gh::rprocs & pss);
  
  int ManualCoordinateList (cGH const * const cctkGH,
                            gh const & hh,
                            gh::mexts  & bbsss,
                            gh::rbnds  & obss,
                            gh::rprocs & pss);

  int Moving (cGH const * const cctkGH,
              gh const & hh,
              gh::mexts  & bbsss,
              gh::rbnds  & obss,
              gh::rprocs & pss);
  
  int Automatic (cGH const * const cctkGH,
                 gh const & hh,
                 gh::mexts  & bbsss,
                 gh::rbnds  & obss,
                 gh::rprocs & pss);
  
  void Automatic_OneLevel (const cGH * const cctkGH,
                           const gh & hh,
                           const int rl,
                           const int numrl,
                           const gf<CCTK_REAL> & errorvar,
                           vector<ibbox> & bbs,
                           vector<bbvect> & obs);
  
  void Automatic_Recursive (const cGH * const cctkGH,
                            const gh & hh,
                            const data<CCTK_REAL> & errorvar,
                            list<ibbox> & bbl,
                            const ibbox & region,
                            const ivect & reffact);
  
  void Automatic_Recombine (list<ibbox> & bbl1,
                            list<ibbox> & bbl2,
                            list<ibbox> & bbl,
                            const ibbox & iface,
                            const int dir);
  
} // namespace CarpetRegrid

#endif // !defined(CARPETREGRID_HH)
