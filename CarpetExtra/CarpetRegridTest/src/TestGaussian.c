#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"


static const char *rcsid = "$Header: /home/eschnett/C/carpet/Carpet/CarpetExtra/CarpetRegridTest/src/TestGaussian.c,v 1.1 2004/05/23 15:56:14 cott Exp $";

CCTK_FILEVERSION(CarpetExtra_CarpetRegridTest_TestGaussian_c);

void CarpetRegrid_TestGaussian(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  int i,j,k;

  int index;
  CCTK_REAL X, Y, Z, R;


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

	      phi_error[index] = phi[index] - amplitude*exp( - pow( (R - radius) / sigma, 2.0 ) );

	    }
	}
    }

  //  CCTK_VInfo(CCTK_THORNSTRING,"Performed CarpetRegridTest");

}
