// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetRegrid/src/regrid.hh,v 1.12 2004/04/14 22:19:44 schnetter Exp $

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
                                  CCTK_INT const reflevel,
                                  CCTK_INT const map,
                                  CCTK_INT const size,
                                  CCTK_INT const * const nboundaryzones,
                                  CCTK_INT const * const is_internal,
                                  CCTK_INT const * const is_staggered,
                                  CCTK_INT const * const shiftout,
                                  CCTK_POINTER const bbsss_,
                                  CCTK_POINTER const obss_,
                                  CCTK_POINTER const pss_);
  }
  
  
  
  int BaseLevel (cGH const * const cctkGH,
                 gh<dim> const & hh,
                 int const reflevel,
                 int const map,
                 int const size,
                 jjvect const & nboundaryzones,
                 jjvect const & is_internal,
                 jjvect const & is_staggered,
                 jjvect const & shiftout,
                 gh<dim>::rexts  & bbsss,
                 gh<dim>::rbnds  & obss,
                 gh<dim>::rprocs & pss);

  int Centre (cGH const * const cctkGH,
              gh<dim> const & hh,
              int const reflevel,
              int const map,
              int const size,
              jjvect const & nboundaryzones,
              jjvect const & is_internal,
              jjvect const & is_staggered,
              jjvect const & shiftout,
              gh<dim>::rexts  & bbsss,
              gh<dim>::rbnds  & obss,
              gh<dim>::rprocs & pss);

  int ManualGridpoints (cGH const * const cctkGH,
                        gh<dim> const & hh,
                        int const reflevel,
                        int const map,
                        int const size,
                        jjvect const & nboundaryzones,
                        jjvect const & is_internal,
                        jjvect const & is_staggered,
                        jjvect const & shiftout,
                        gh<dim>::rexts  & bbsss,
                        gh<dim>::rbnds  & obss,
                        gh<dim>::rprocs & pss);
  
  void ManualGridpoints_OneLevel (const cGH * const cctkGH,
                                  const gh<dim> & hh,
                                  const int reflevel,
                                  const int reflevels,
                                  const ivect ilower,
                                  const ivect iupper,
                                  const bbvect obound,
                                  vector<ibbox> & bbs,
                                  vector<bbvect> & obs);

  int ManualCoordinates (cGH const * const cctkGH,
                         gh<dim> const & hh,
                         int const reflevel,
                         int const map,
                         int const size,
                         jjvect const & nboundaryzones,
                         jjvect const & is_internal,
                         jjvect const & is_staggered,
                         jjvect const & shiftout,
                         gh<dim>::rexts  & bbsss,
                         gh<dim>::rbnds  & obss,
                         gh<dim>::rprocs & pss);
  
  void ManualCoordinates_OneLevel (const cGH * const cctkGH,
                                   const gh<dim> & hh,
                                   const int reflevel,
                                   const int reflevels,
                                   const rvect ilower,
                                   const rvect iupper,
                                   const bbvect obound,
                                   vector<ibbox> & bbs,
                                   vector<bbvect> & obs);
  
  ivect delta2int (const cGH * const cctkGH, const gh<dim>& hh,
                   const rvect & rpos, const int reflevel);
  ivect pos2int (const cGH* const cctkGH, const gh<dim>& hh,
                 const rvect & rpos, const int reflevel);

  int ManualGridpointList (cGH const * const cctkGH,
                           gh<dim> const & hh,
                           int const reflevel,
                           int const map,
                           int const size,
                           jjvect const & nboundaryzones,
                           jjvect const & is_internal,
                           jjvect const & is_staggered,
                           jjvect const & shiftout,
                           gh<dim>::rexts  & bbsss,
                           gh<dim>::rbnds  & obss,
                           gh<dim>::rprocs & pss);
  
  int ManualCoordinateList (cGH const * const cctkGH,
                            gh<dim> const & hh,
                            int const reflevel,
                            int const map,
                            int const size,
                            jjvect const & nboundaryzones,
                            jjvect const & is_internal,
                            jjvect const & is_staggered,
                            jjvect const & shiftout,
                            gh<dim>::rexts  & bbsss,
                            gh<dim>::rbnds  & obss,
                            gh<dim>::rprocs & pss);

  int Moving (cGH const * const cctkGH,
              gh<dim> const & hh,
              int const reflevel,
              int const map,
              int const size,
              jjvect const & nboundaryzones,
              jjvect const & is_internal,
              jjvect const & is_staggered,
              jjvect const & shiftout,
              gh<dim>::rexts  & bbsss,
              gh<dim>::rbnds  & obss,
              gh<dim>::rprocs & pss);
  
  int Automatic (cGH const * const cctkGH,
                 gh<dim> const & hh,
                 int const reflevel,
                 int const map,
                 int const size,
                 jjvect const & nboundaryzones,
                 jjvect const & is_internal,
                 jjvect const & is_staggered,
                 jjvect const & shiftout,
                 gh<dim>::rexts  & bbsss,
                 gh<dim>::rbnds  & obss,
                 gh<dim>::rprocs & pss);
  
  void Automatic_OneLevel (const cGH * const cctkGH,
                           const gh<dim> & hh,
                           const int reflevel,
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
