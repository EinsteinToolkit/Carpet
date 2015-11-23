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
void CarpetAdaptiveRegridParamcheck(CCTK_ARGUMENTS);

/* Aliased functions */
//     CCTK_INT CarpetAdaptiveRegrid_Regrid (const cGH * const cctkGH,
//                                   gh<dim>::rregs * regsss);
CCTK_INT CarpetAdaptiveRegrid_Regrid(CCTK_POINTER_TO_CONST const cctkGH_,
                                     CCTK_POINTER const regsss_,
                                     CCTK_INT force);
}

int ManualCoordinateList(cGH const *const cctkGH, gh const &hh,
                         gh::mregs &regsss, gh::mregs &local_regsss);

void ManualCoordinates_OneLevel(const cGH *const cctkGH, const gh &hh,
                                const int rl, const int numrl,
                                const rvect lower, const rvect upper,
                                const bbvect obound, const bbvect rbound,
                                vector<ibbox> &bbs, vector<bbvect> &obs,
                                vector<bbvect> &rbs);

void ManualGridpoints_OneLevel(const cGH *const cctkGH, const gh &hh,
                               const int rl, const int numrl,
                               const ivect ilower, const ivect iupper,
                               const bbvect obound, const bbvect rbound,
                               vector<ibbox> &bbs, vector<bbvect> &obs,
                               vector<bbvect> &rbs);

rvect int2pos(const cGH *const cctkGH, const gh &hh, const ivect &ipos,
              const int rl);

ivect pos2int(const cGH *const cctkGH, const gh &hh, const rvect &rpos,
              const int rl);

} // namespace CarpetAdaptiveregrid

#endif // !defined(CARPETADAPTIVEREGRID_HH)
