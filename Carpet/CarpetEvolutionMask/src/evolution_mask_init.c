#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "util_Table.h"



void
EvolutionMaskBase_InitEvolutionMask (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  int i, j, k;
  
  
  
  /* Initialise the mask to 1 everywhere */
  if (verbose) {
    CCTK_INFO ("Initialising the mask to 1");
  }
  
  for (k=0; k<cctk_lsh[2]; ++k) {
    for (j=0; j<cctk_lsh[1]; ++j) {
      for (i=0; i<cctk_lsh[0]; ++i) {
        
        int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
        evolution_mask[ind] = 1.0;
        
      }
    }
  }

}
