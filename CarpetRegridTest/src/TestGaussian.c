#include <math.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void CarpetRegrid_TestGaussian(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS_CarpetRegrid_TestGaussian;

  int i, j, k;

  int index;
  CCTK_REAL X, Y, Z, R;

  for (k = 0; k < cctk_lsh[2]; k++) {
    for (j = 0; j < cctk_lsh[1]; j++) {
      for (i = 0; i < cctk_lsh[0]; i++) {
        index = CCTK_GFINDEX3D(cctkGH, i, j, k);

        X = x[index];
        Y = y[index];
        Z = z[index];

        R = sqrt(X * X + Y * Y + Z * Z);

        phi_error[index] =
            phi[index] - amplitude * exp(-pow((R - radius) / sigma, 2.0));
        phi_relerror[index] =
            phi_error[index] /
            (amplitude * exp(-pow((R - radius) / sigma, 2.0)));
      }
    }
  }

  /*  CCTK_VInfo(CCTK_THORNSTRING,"Performed CarpetRegridTest"); */
}
