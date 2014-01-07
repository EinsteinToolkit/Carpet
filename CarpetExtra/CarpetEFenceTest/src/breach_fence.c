#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void CarpetEFenceTest_BreachFence(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  int idx;
  if(CCTK_EQUALS(breach_at,"lower")) {
    idx = -1;
  } else if(CCTK_EQUALS(breach_at,"upper")) {
    idx = CCTK_GFINDEX3D(cctkGH, cctk_lsh[0], cctk_lsh[1]-1, cctk_lsh[2]-1);
  } else {
    /* this should never happen */
    CCTK_ERROR("Uknown breach_at value");
  }

  sheep[idx] = 42;
}
