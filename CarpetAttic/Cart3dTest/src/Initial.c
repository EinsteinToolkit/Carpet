#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "AlphaThorns/Cart3d/src/Cart3d.h"



int Cart3dTest_Initial (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  int r;
  int i,j,k;
  
  int vi[9];
  
#if 0
  r = Cart3dSetTensorTypeVN (cctkGH, "Cart3dTest::s", "SCALAR");
  assert (r>=0);
  
  r = Cart3dSetTensorTypeVN
    (cctkGH, "Cart3dTest::vx Cart3dTest::vy Cart3dTest::vz", "VECTOR");
  assert (r>=0);
  
  r = Cart3dSetTensorTypeVN
    (cctkGH, "Cart3dTest::txx Cart3dTest::txy Cart3dTest::txz"
     " Cart3dTest::tyy Cart3dTest::tyz Cart3dTest::tzz", "SYMMTENSOR");
  assert (r>=0);
  
  r = Cart3dSetTensorTypeVN
    (cctkGH, "Cart3dTest::az Cart3dTest::ay Cart3dTest::ax", "ANTISYMMTENSOR");
  assert (r>=0);
#endif
  
  vi[0] = CCTK_VarIndex("Cart3dTest::s");
  r = Cart3dSetTensorTypeVI (cctkGH, 1, vi, "SCALAR");
  assert (r>=0);
  
  vi[0] = CCTK_VarIndex("Cart3dTest::vx");
  vi[1] = CCTK_VarIndex("Cart3dTest::vy");
  vi[2] = CCTK_VarIndex("Cart3dTest::vz");
  r = Cart3dSetTensorTypeVI (cctkGH, 3, vi, "VECTOR");
  assert (r>=0);
  
  vi[0] = CCTK_VarIndex("Cart3dTest::txx");
  vi[1] = CCTK_VarIndex("Cart3dTest::txy");
  vi[2] = CCTK_VarIndex("Cart3dTest::txz");
  vi[3] = CCTK_VarIndex("Cart3dTest::tyy");
  vi[4] = CCTK_VarIndex("Cart3dTest::tyz");
  vi[5] = CCTK_VarIndex("Cart3dTest::tzz");
  r = Cart3dSetTensorTypeVI (cctkGH, 6, vi, "SYMMTENSOR");
  assert (r>=0);
  
  vi[0] = CCTK_VarIndex("Cart3dTest::fxx");
  vi[1] = CCTK_VarIndex("Cart3dTest::fxy");
  vi[2] = CCTK_VarIndex("Cart3dTest::fxz");
  vi[3] = CCTK_VarIndex("Cart3dTest::fyx");
  vi[4] = CCTK_VarIndex("Cart3dTest::fyy");
  vi[5] = CCTK_VarIndex("Cart3dTest::fyz");
  vi[6] = CCTK_VarIndex("Cart3dTest::fzx");
  vi[7] = CCTK_VarIndex("Cart3dTest::fzy");
  vi[8] = CCTK_VarIndex("Cart3dTest::fzz");
  r = Cart3dSetTensorTypeVI (cctkGH, 9, vi, "TENSOR");
  assert (r>=0);
  
  vi[0] = CCTK_VarIndex("Cart3dTest::az"); /* watch the order! */
  vi[1] = CCTK_VarIndex("Cart3dTest::ay");
  vi[2] = CCTK_VarIndex("Cart3dTest::ax");
  r = Cart3dSetTensorTypeVI (cctkGH, 3, vi, "ANTISYMMTENSOR");
  assert (r>=0);
  
  for (k=0; k<cctk_lsh[2]; ++k) {
    for (j=0; j<cctk_lsh[1]; ++j) {
      for (i=0; i<cctk_lsh[0]; ++i) {
	
	const int ii = cctk_lbnd[0] + i;
	const int jj = cctk_lbnd[1] + j;
	const int kk = cctk_lbnd[2] + k;
	
	const int ind = CCTK_GFINDEX3D(cctkGH, i,j,k);
	
	
	
	/* store the position in the components */
	
	s[ind]   = (( 1 * 100 + ii) * 100 + jj) * 100 + kk;
	
	vx[ind]  = ((11 * 100 + ii) * 100 + jj) * 100 + kk;
	vy[ind]  = ((12 * 100 + ii) * 100 + jj) * 100 + kk;
	vz[ind]  = ((13 * 100 + ii) * 100 + jj) * 100 + kk;
	
	txx[ind] = ((21 * 100 + ii) * 100 + jj) * 100 + kk;
	txy[ind] = ((22 * 100 + ii) * 100 + jj) * 100 + kk;
	txz[ind] = ((23 * 100 + ii) * 100 + jj) * 100 + kk;
	tyy[ind] = ((24 * 100 + ii) * 100 + jj) * 100 + kk;
	tyz[ind] = ((25 * 100 + ii) * 100 + jj) * 100 + kk;
	tzz[ind] = ((26 * 100 + ii) * 100 + jj) * 100 + kk;
	
	ax[ind]  = ((31 * 100 + ii) * 100 + jj) * 100 + kk;
	ay[ind]  = ((32 * 100 + ii) * 100 + jj) * 100 + kk;
	az[ind]  = ((33 * 100 + ii) * 100 + jj) * 100 + kk;
	
	fxx[ind] = ((41 * 100 + ii) * 100 + jj) * 100 + kk;
	fxy[ind] = ((42 * 100 + ii) * 100 + jj) * 100 + kk;
	fxz[ind] = ((43 * 100 + ii) * 100 + jj) * 100 + kk;
	fyx[ind] = ((44 * 100 + ii) * 100 + jj) * 100 + kk;
	fyy[ind] = ((45 * 100 + ii) * 100 + jj) * 100 + kk;
	fyz[ind] = ((46 * 100 + ii) * 100 + jj) * 100 + kk;
	fzx[ind] = ((47 * 100 + ii) * 100 + jj) * 100 + kk;
	fzy[ind] = ((48 * 100 + ii) * 100 + jj) * 100 + kk;
	fzz[ind] = ((49 * 100 + ii) * 100 + jj) * 100 + kk;
	
      }
    }
  }
  
  return 0;
}
