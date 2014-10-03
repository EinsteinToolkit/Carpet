#ifndef CARPETREGRID_HH
#define CARPETREGRID_HH

#include <list>
#include <vector>

#include "cctk.h"
#include "cctk_Arguments.h"

#include "bbox.hh"
#include "gf.hh"
#include "gh.hh"
#include "region.hh"
#include "vect.hh"

#include "carpet.hh"



namespace CarpetRegrid {
  
  using namespace std;
  using namespace Carpet;
  
  
  
  extern "C" {
    
    // Scheduled functions
    void CarpetRegridParamcheck (CCTK_ARGUMENTS);
    
    // Aliased functions
//     CCTK_INT CarpetRegrid_Regrid (const cGH * const cctkGH,
//                                   gh<dim>::rregs * superregss,
//                                   gh<dim>::mregs * regsss);
    CCTK_INT CarpetRegrid_Regrid (CCTK_POINTER_TO_CONST const cctkGH_,
                                  CCTK_POINTER const superregss_,
                                  CCTK_POINTER const regsss_,
				  CCTK_INT force);
  }
  
  
  
  int BaseLevel (cGH const * const cctkGH,
                 gh const & hh,
                 gh::rregs & regss);

  int Centre (cGH const * const cctkGH,
              gh const & hh,
              gh::rregs & regss);

  int ManualGridpoints (cGH const * const cctkGH,
                        gh const & hh,
                        gh::rregs & regss);
  
  void ManualGridpoints_OneLevel (const cGH * const cctkGH,
                                  const gh & hh,
                                  const int rl,
                                  const int numrl,
                                  const ivect ilower,
                                  const ivect iupper,
                                  const region_t & reg,
                                  vector<region_t> & regs);

  int ManualCoordinates (cGH const * const cctkGH,
                         gh const & hh,
                         gh::rregs & regss);
  
  void ManualCoordinates_OneLevel (const cGH * const cctkGH,
                                   const gh & hh,
                                   const int rl,
                                   const int numrl,
                                   const rvect lower,
                                   const rvect upper,
                                   const region_t & reg,
                                   vector<region_t> & regs);
  
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
                           gh::rregs & regss);
  
  int ManualCoordinateList (cGH const * const cctkGH,
                            gh const & hh,
                            gh::rregs & regss);

  int Moving (cGH const * const cctkGH,
              gh const & hh,
              gh::rregs & regss);
  
  int Automatic (cGH const * const cctkGH,
                 gh const & hh,
                 gh::rregs & regss);
  
  void Automatic_OneLevel (const cGH * const cctkGH,
                           const gh & hh,
                           const int rl,
                           const int numrl,
                           const gf<CCTK_REAL> & errorvar,
                           vector<region_t> & regs);
  
  void Automatic_Recursive (const cGH * const cctkGH,
                            const gh & hh,
                            const data<CCTK_REAL> & errorvar,
                            list<region_t> & regl,
                            const region_t & region,
                            const ivect & reffact);
  
  void Automatic_Recombine (list<region_t> & regl1,
                            list<region_t> & regl2,
                            list<region_t> & regl,
                            const ibbox & iface,
                            const int dir);
  
} // namespace CarpetRegrid

#endif // #ifndef CARPETREGRID_HH
