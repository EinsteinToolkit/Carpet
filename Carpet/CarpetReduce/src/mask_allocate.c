#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <assert.h>



void
MaskBase_AllocateMask (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  /* Allocate helpers for the weight function */
  if (verbose) {
    CCTK_INFO ("Allocating weight function helpers");
  }
  
  int const ierr1 = CCTK_EnableGroupStorage (cctkGH, "CarpetReduce::iweight");
  int const ierr2 = CCTK_EnableGroupStorage (cctkGH, "CarpetReduce::one");
  assert (!ierr1);
  assert (!ierr2);
}
