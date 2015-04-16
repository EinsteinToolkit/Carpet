#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

void CustomGhosts_Initialize(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  const int rl = GetRefinementLevel(cctkGH);
  const int r = CCTK_MyProc(cctkGH);

#pragma omp parallel
  CCTK_LOOP3_ALL(CustomGhost_Initialize, cctkGH, i,j,k) {
    const CCTK_INT idx = CCTK_GFINDEX3D(cctkGH, i,j,k);
    custom[idx] = rl * 100. + r;
    regular[idx] = rl * 100. + r;
  } CCTK_ENDLOOP3_ALL(CustomGhost_Initialize);
}

void CustomGhosts_Sync(CCTK_ARGUMENTS)
{
  // nop, only exists so that I can tie the SYNC to it
}
