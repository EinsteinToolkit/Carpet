#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>



void
MaskBase_AllocateMask (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  /* Allocate helpers for the weight function */
  if (verbose) {
    CCTK_INFO ("Allocating weight function helpers");
  }
  
  CCTK_EnableGroupStorage (cctkGH, "CarpetReduce::iweight");
  CCTK_EnableGroupStorage (cctkGH, "CarpetReduce::one");
}
