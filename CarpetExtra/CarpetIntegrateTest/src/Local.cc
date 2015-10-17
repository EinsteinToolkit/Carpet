#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"



extern "C" void CarpetIntegrate_Local(CCTK_ARGUMENTS);

void CarpetIntegrate_Local(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int nx = cctk_lsh[0];
  int ny = cctk_lsh[1];
  int nz = cctk_lsh[2];

  for (int k=0;k<nz;k++) {
    for (int j=0;j<ny;j++) {
      for (int i=0;i<nx;i++) {
	int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
        
	integrand[index] = constant + timefact * cctk_time;
      }
    }
  }
}
