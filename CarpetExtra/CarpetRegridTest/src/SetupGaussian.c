#include <math.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"


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

	      R = r[index];

	      phi[index] = amplitude*exp( - pow( (R - radius) / sigma, 2.0 ) );
	    }
	}
    }

  CCTK_VInfo(CCTK_THORNSTRING,"Gaussian initial data have been set.");

}
