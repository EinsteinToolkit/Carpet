#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <loopcontrol.h>

#include "bits.h"



void
MaskBase_SetMask (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose) {
    int const reflevel = GetRefinementLevel(cctkGH);
    CCTK_VInfo (CCTK_THORNSTRING,
                "Finalise the weight on level %d", reflevel);
  }
  
  CCTK_REAL const factor = 1.0 / BMSK(cctk_dim);
#pragma omp parallel
  LC_LOOP3(MaskBase_InitMask_interior,
           i,j,k,
           0,0,0, cctk_lsh[0],cctk_lsh[1],cctk_lsh[2],
           cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
    weight[ind] = factor * BCNT(iweight[ind]);
    one[ind] = 1.0;
  } LC_ENDLOOP3(MaskBase_InitMask_interior);
}
