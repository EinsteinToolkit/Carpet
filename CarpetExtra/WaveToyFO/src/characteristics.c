/* $Header: /home/eschnett/C/carpet/Carpet/CarpetExtra/WaveToyFO/src/Attic/characteristics.c,v 1.1 2004/05/07 20:55:47 schnetter Exp $ */

#include <assert.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"



void CCTK_FCALL
CCTK_FNAME (WaveToyFO_calc_inv4) (CCTK_REAL xform[4][4],
				  CCTK_REAL xform1[4][4]);



CCTK_INT
WaveToyFO_MultiPatch_Prim2Char (CCTK_POINTER_TO_CONST const cctkGH_,
				CCTK_INT const normal[],
				CCTK_INT const lbnd[],
				CCTK_INT const lsh[],
				CCTK_INT const rhs_flag,
				CCTK_INT const num_modes,
				CCTK_POINTER const modes[],
				CCTK_REAL speeds[])
{
  cGH const * restrict const cctkGH = cctkGH_;
  DECLARE_CCTK_ARGUMENTS;
  
  CCTK_INT tangent[2][3];
  
  CCTK_REAL const * restrict prims[4];
  CCTK_REAL       * restrict chars[4];
  
  CCTK_REAL xform[4][4];
  
  int n;			/* mode */
  int i, j, k;
  int a, b;
  int d;
  
  assert (cctk_dim == 3);
  assert (num_modes == 4);
  for (d=0; d<3; ++d) {
    assert (lsh[d] >= 0);
  }
  for (n=0; n<num_modes; ++n) {
    assert (modes[n]);
  }
  assert (speeds);
  
  for (d=0; d<3; ++d) {
    tangent[0][d] = normal[(d+1)%3];
    tangent[1][d] = normal[(d+2)%3];
  }
  
  xform[0][0] = +1;
  xform[1][0] = +1;
  xform[2][0] =  0;
  xform[3][0] =  0;
  for (d=0; d<3; ++d) {
    xform[0][1+d] = +normal[d];
    xform[1][1+d] = -normal[d];
    xform[2][1+d] = tangent[0][d];
    xform[3][1+d] = tangent[1][d];
  }
  
  if (rhs_flag) {
    prims[0] = phidot;
    prims[1] = psixdot;
    prims[2] = psiydot;
    prims[3] = psizdot;
  } else {
    prims[0] = phi;
    prims[1] = psix;
    prims[2] = psiy;
    prims[3] = psiz;
  }
  for (n=0; n<4; ++n) {
    chars[n] = modes[n];
  }
  
  for (k=0; k<lsh[2]; ++k) {
    for (j=0; j<lsh[1]; ++j) {
      for (i=0; i<lsh[0]; ++i) {
	int const ind3d
	  = CCTK_GFINDEX3D (cctkGH, lbnd[0]+i, lbnd[1]+j, lbnd[2]+k);
	int const ind2d = i + lsh[0] * (j + lsh[1] * k);
	
	for (a=0; a<4; ++a) {
	  chars[a][ind2d] = 0;
	  for (b=0; b<4; ++b) {
	    chars[a][ind2d] += xform[a][b] * prims[b][ind3d];
	  }
	}
	
      }
    }
  }
  
  speeds[0] = +1.0; 		/* incoming */
  speeds[1] = -1.0;		/* outgoing */
  speeds[2] =  0.0;		/* zero speed */
  speeds[3] =  0.0;
  
  return 0;
}



CCTK_INT
WaveToyFO_MultiPatch_Char2Prim (CCTK_POINTER_TO_CONST const cctkGH_,
				CCTK_INT const normal[],
				CCTK_INT const lbnd[],
				CCTK_INT const lsh[],
				CCTK_INT const rhs_flag,
				CCTK_INT const num_modes,
				CCTK_POINTER_TO_CONST const modes[])
{
  cGH const * restrict const cctkGH = cctkGH_;
  DECLARE_CCTK_ARGUMENTS;
  
  CCTK_INT tangent[2][3];
  
  CCTK_REAL const * restrict chars[4];
  CCTK_REAL       * restrict prims[4];
  
  CCTK_REAL xform[4][4];
  CCTK_REAL xform1[4][4];
  
  int n;			/* mode */
  int i, j, k;
  int a, b;
  int d;
  
  assert (cctk_dim == 3);
  assert (num_modes == 4);
  for (d=0; d<3; ++d) {
    assert (lsh[d] >= 0);
  }
  for (n=0; n<num_modes; ++n) {
    assert (modes[n]);
  }
  
  for (d=0; d<3; ++d) {
    tangent[0][d] = normal[(d+1)%3];
    tangent[1][d] = normal[(d+2)%3];
  }
  
  xform[0][0] = +1;
  xform[1][0] = +1;
  xform[2][0] =  0;
  xform[3][0] =  0;
  for (d=0; d<3; ++d) {
    xform[0][1+d] = +normal[d];
    xform[1][1+d] = -normal[d];
    xform[2][1+d] = tangent[0][d];
    xform[3][1+d] = tangent[1][d];
  }
  
  CCTK_FNAME (WaveToyFO_calc_inv4) (xform, xform1);
  
  for (n=0; n<4; ++n) {
    chars[n] = modes[n];
  }
  if (rhs_flag) {
    prims[0] = phidot;
    prims[1] = psixdot;
    prims[2] = psiydot;
    prims[3] = psizdot;
  } else {
    prims[0] = phi;
    prims[1] = psix;
    prims[2] = psiy;
    prims[3] = psiz;
  }
  
  for (k=0; k<lsh[2]; ++k) {
    for (j=0; j<lsh[1]; ++j) {
      for (i=0; i<lsh[0]; ++i) {
	int const ind3d
	  = CCTK_GFINDEX3D (cctkGH, lbnd[0]+i, lbnd[1]+j, lbnd[2]+k);
	int const ind2d = i + lsh[0] * (j + lsh[1] * k);
	
	for (a=0; a<4; ++a) {
	  prims[a][ind3d] = 0;
	  for (b=0; b<4; ++b) {
	    prims[a][ind3d] += xform1[a][b] * chars[b][ind2d];
	  }
	}
	
      }
    }
  }
  
  return 0;
}
