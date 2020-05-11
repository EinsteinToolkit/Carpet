#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "util_Table.h"
#include <assert.h>
#include <math.h>

void CoM2_Local(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_CoM2_Local;
  DECLARE_CCTK_PARAMETERS;

  int i, j, k, index;

  int nx = cctk_lsh[0];
  int ny = cctk_lsh[1];
  int nz = cctk_lsh[2];

  for (k = 0; k < nz; k++)
    for (j = 0; j < ny; j++)
      for (i = 0; i < nx; i++) {

        index = CCTK_GFINDEX3D(cctkGH, i, j, k);

        redvar[index] = 0.45e0;
      }
}
