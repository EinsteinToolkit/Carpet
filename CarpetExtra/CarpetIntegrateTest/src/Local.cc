// $Header:$

#include <math.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include "carpet.hh"
#include "defs.hh"
#include "vect.hh"

#include "util_Table.h"



extern "C" { CCTK_INT CarpetIntegrate_Local(CCTK_ARGUMENTS);
}

CCTK_INT CarpetIntegrate_Local(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  using namespace std;
  using namespace Carpet;

  int nx = cctk_lsh[0];
  int ny = cctk_lsh[1];
  int nz = cctk_lsh[2];

  for (int k=0;k<nz;k++) {
    for (int j=0;j<ny;j++) {
      for (int i=0;i<nx;i++) {
	int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
        
	integrand[index] = 1;
      }
    }
  }
  
  return 0;
}
