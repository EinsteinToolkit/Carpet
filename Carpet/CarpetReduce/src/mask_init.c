#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <loopcontrol.h>



void
MaskBase_InitMask (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  /* Initialise the weight to 1 everywhere */
  if (verbose) {
    CCTK_INFO ("Initialising weight to 1");
  }
  
#pragma omp parallel
  LC_LOOP3(MaskBase_InitMask_interior,
           i,j,k,
           0,0,0, cctk_lsh[0],cctk_lsh[1],cctk_lsh[2],
           cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
    weight[ind] = 1.0;
  } LC_ENDLOOP3(MaskBase_InitMask_interior);
}
