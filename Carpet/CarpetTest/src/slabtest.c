#include <assert.h>
#include <stdlib.h>

#include "cctk.h"
#include "cctk_Arguments.h"

void
carpettest_slabtest (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  
  int vi;
  
  CCTK_INT directions[3][3];
  CCTK_INT origin[3];
  CCTK_INT extent[3];
  CCTK_INT hsize[3];
  
  int mapping;
  
  int d, dd;
  
  int ierr;
  
  
  
  assert (CCTK_nProcs(cctkGH) == 1);
  
  vi = CCTK_VarIndex ("grid::y");
  assert (vi >= 0);
  
  for (d=0; d<3; ++d) {
    for (dd=0; dd<3; ++dd) {
      directions[d][dd] = d == dd;
    }
    origin[d] = cctk_lbnd[d];
    extent[d] = cctk_lsh[d];
  }
  
  mapping = Hyperslab_GlobalMappingByIndex
    (cctkGH, vi, 3, &directions[0][0], origin, extent, NULL, -1, NULL, hsize);
  assert (mapping >= 0);
  
  for (d=0; d<3; ++d) {
    assert (hsize[d] == cctk_lsh[d]);
  }
  
  ierr = Hyperslab_Get (cctkGH, mapping, -1, vi, 0, CCTK_VARIABLE_REAL, yy);
  assert (! ierr);
  
  ierr = Hyperslab_FreeMapping (mapping);
  assert (! ierr);
}
