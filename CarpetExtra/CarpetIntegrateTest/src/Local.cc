// $Header: /home/eschnett/C/carpet/Carpet/CarpetExtra/CarpetIntegrateTest/src/Local.cc,v 1.1 2004/09/01 17:47:53 schnetter Exp $

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include "carpet.hh"
#include "defs.hh"
#include "vect.hh"


#include "util_Table.h"
#include "math.h"

extern "C" { CCTK_INT CarpetIntegrate_Local(CCTK_ARGUMENTS);
}

static void SpatialDeterminant(CCTK_REAL gxx,
                               CCTK_REAL gxy,
                               CCTK_REAL gxz,
                               CCTK_REAL gyy,
                               CCTK_REAL gyz,
                               CCTK_REAL gzz,
                               CCTK_REAL *detg);

CCTK_INT CarpetIntegrate_Local(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  using namespace std;
  using namespace Carpet;

  int nx = cctk_lsh[0];
  int ny = cctk_lsh[1];
  int nz = cctk_lsh[2];

  for (int k=0;k<nz;k++)
    for (int j=0;j<ny;j++)
      for (int i=0;i<nx;i++) {
	int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

	/*
	SpatialDeterminant(gxx[index],gxy[index],gxz[index],
			   gyy[index],gyz[index],gzz[index],
			   &detg);
	*/

	integrand[index] = 1;

      }

  return 0;
}


void SpatialDeterminant(CCTK_REAL gxx,
			CCTK_REAL gxy,
			CCTK_REAL gxz,
			CCTK_REAL gyy,
			CCTK_REAL gyz,
			CCTK_REAL gzz,
			CCTK_REAL *detg)
{

  *detg = -gxz*gxz*gyy + 2.0*gxy*gxz*gyz 
    - gxx*gyz*gyz - gxy*gxy*gzz + gxx*gyy*gzz;

  return;
}
