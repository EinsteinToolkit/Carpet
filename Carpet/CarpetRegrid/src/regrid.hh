// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetRegrid/src/regrid.hh,v 1.14 2004/06/02 07:08:52 bzink Exp $

#ifndef CARPETREGRID_HH
#define CARPETREGRID_HH

#include <list>

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
//                                   gh<dim>::rexts  * bbsss,
//                                   gh<dim>::rbnds  * obss,
//                                   gh<dim>::rprocs * pss);
    CCTK_INT CarpetRegrid_Regrid (CCTK_POINTER_TO_CONST const cctkGH_,
                                  CCTK_POINTER const bbsss_,
                                  CCTK_POINTER const obss_,
                                  CCTK_POINTER const pss_,
				  CCTK_INT force);
  }
  
  
  
  int BaseLevel (cGH const * const cctkGH,
                 gh<dim> const & hh,
                 gh<dim>::rexts  & bbsss,
                 gh<dim>::rbnds  & obss,
                 gh<dim>::rprocs & pss);

  int Centre (cGH const * const cctkGH,
              gh<dim> const & hh,
              gh<dim>::rexts  & bbsss,
              gh<dim>::rbnds  & obss,
              gh<dim>::rprocs & pss);

  int ManualGridpoints (cGH const * const cctkGH,
                        gh<dim> const & hh,
                        gh<dim>::rexts  & bbsss,
                        gh<dim>::rbnds  & obss,
                        gh<dim>::rprocs & pss);
  
  void ManualGridpoints_OneLevel (const cGH * const cctkGH,
                                  const gh<dim> & hh,
                                  const int rl,
                                  const int numrl,
                                  const ivect ilower,
                                  const ivect iupper,
                                  const bbvect obound,
                                  vector<ibbox> & bbs,
                                  vector<bbvect> & obs);

  int ManualCoordinates (cGH const * const cctkGH,
                         gh<dim> const & hh,
                         gh<dim>::rexts  & bbsss,
                         gh<dim>::rbnds  & obss,
                         gh<dim>::rprocs & pss);
  
  void ManualCoordinates_OneLevel (const cGH * const cctkGH,
                                   const gh<dim> & hh,
                                   const int rl,
                                   const int numrl,
                                   const rvect lower,
                                   const rvect upper,
                                   const bbvect obound,
                                   vector<ibbox> & bbs,
                                   vector<bbvect> & obs);
  
  ivect delta2int (const cGH * const cctkGH,
                   const gh<dim>& hh,
                   const rvect & rpos,
                   const int rl);
  ivect pos2int (const cGH* const cctkGH,
                 const gh<dim>& hh,
                 const rvect & rpos,
                 const int rl);

  int ManualGridpointList (cGH const * const cctkGH,
                           gh<dim> const & hh,
                           gh<dim>::rexts  & bbsss,
                           gh<dim>::rbnds  & obss,
                           gh<dim>::rprocs & pss);
  
  int ManualCoordinateList (cGH const * const cctkGH,
                            gh<dim> const & hh,
                            gh<dim>::rexts  & bbsss,
                            gh<dim>::rbnds  & obss,
                            gh<dim>::rprocs & pss);

  int Moving (cGH const * const cctkGH,
              gh<dim> const & hh,
              gh<dim>::rexts  & bbsss,
              gh<dim>::rbnds  & obss,
              gh<dim>::rprocs & pss);
  
  int Automatic (cGH const * const cctkGH,
                 gh<dim> const & hh,
                 gh<dim>::rexts  & bbsss,
                 gh<dim>::rbnds  & obss,
                 gh<dim>::rprocs & pss);
  
  void Automatic_OneLevel (const cGH * const cctkGH,
                           const gh<dim> & hh,
                           const int rl,
                           const int numrl,
                           const int minwidth,
                           const CCTK_REAL minfraction,
                           const CCTK_REAL maxerror,
                           const gf<CCTK_REAL,dim> & errorvar,
                           vector<ibbox> & bbs,
                           vector<bbvect> & obs);
  
  void Automatic_Recursive (const cGH * const cctkGH,
                            const gh<dim> & hh,
                            const int minwidth,
                            const CCTK_REAL minfraction,
                            const CCTK_REAL maxerror,
                            const data<CCTK_REAL,dim> & errorvar,
                            list<ibbox> & bbl,
                            const ibbox & region);
  
  void Automatic_Recombine (list<ibbox> & bbl1,
                            list<ibbox> & bbl2,
                            list<ibbox> & bbl,
                            const ibbox & iface,
                            const int dir);
  
} // namespace CarpetRegrid

#endif // !defined(CARPETREGRID_HH)
