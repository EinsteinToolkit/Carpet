#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>
using std::sqrt;



void TestAlignedAlloc0(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  const int di = 1;
  const int dj = di * cctk_ash[0];
  const int dk = dj * cctk_ash[1];
  
#pragma omp parallel
  CCTK_LOOP3_ALL(all3, cctkGH, i,j,k) {
    const ptrdiff_t ind3d = CCTK_GFINDEX3D(cctkGH, i,j,k);
    const CCTK_REAL xl = x[ind3d];
    const CCTK_REAL yl = y[ind3d];
    const CCTK_REAL zl = z[ind3d];
    const CCTK_REAL ul = sqrt(xl*xl + yl*yl + zl*zl);
    u0[ind3d] = ul;
  } CCTK_ENDLOOP3_ALL(all3);
  
#pragma omp parallel
  CCTK_LOOP3_INT(int3, cctkGH, i,j,k) {
    const ptrdiff_t ind3d = CCTK_GFINDEX3D(cctkGH, i,j,k);
    const CCTK_REAL rl = r[ind3d];
    const CCTK_REAL rm0l = r[ind3d-di];
    const CCTK_REAL rp0l = r[ind3d+di];
    const CCTK_REAL rm1l = r[ind3d-dj];
    const CCTK_REAL rp1l = r[ind3d+dj];
    const CCTK_REAL rm2l = r[ind3d-dk];
    const CCTK_REAL rp2l = r[ind3d+dk];
    const CCTK_REAL dd0r = rm0l - 2*rl + rp0l;
    const CCTK_REAL dd1r = rm1l - 2*rl + rp1l;
    const CCTK_REAL dd2r = rm2l - 2*rl + rp2l;
    const CCTK_REAL ddul = dd0r + dd1r + dd2r;
    ddu0[ind3d] = ddul;
  } CCTK_ENDLOOP3_INT(int3);
  
#pragma omp parallel
  CCTK_LOOP3_BND(bnd3, cctkGH, i,j,k, ni,nj,nk) {
    const ptrdiff_t ind3d = CCTK_GFINDEX3D(cctkGH, i,j,k);
    const CCTK_REAL ddul = 0;
    ddu0[ind3d] = ddul;
  } CCTK_ENDLOOP3_BND(bnd3);
}
