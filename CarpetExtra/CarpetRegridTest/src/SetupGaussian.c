#include <math.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"


static const char *rcsid = "$Header: /home/eschnett/C/carpet/Carpet/CarpetExtra/CarpetRegridTest/src/SetupGaussian.c,v 1.2 2004/05/27 13:22:49 schnetter Exp $";

CCTK_FILEVERSION(CarpetExtra_CarpetRegridTest_SetupGaussian_c);

void CarpetRegrid_SetupGaussian(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  int i,j,k;

  CCTK_REAL omega;
  int index;
  CCTK_REAL X, Y, Z, R;
  CCTK_REAL pi;

  for(k=0; k<cctk_lsh[2]; k++)
    {
      for(j=0; j<cctk_lsh[1]; j++)
	{
	  for(i=0; i<cctk_lsh[0]; i++)
	    {
	      index =  CCTK_GFINDEX3D(cctkGH,i,j,k);

	      X = x[index];
	      Y = y[index];
	      Z = z[index];

	      R = sqrt(X*X + Y*Y + Z*Z);

	      phi[index] = amplitude*exp( - pow( (R - radius) / sigma, 2.0 ) );

	      phi_p[index] = phi[index];
	      phi_p_p[index] = phi[index];
	    }
	}
    }

  CCTK_VInfo(CCTK_THORNSTRING,"Gaussian initial data have been set.");

}
