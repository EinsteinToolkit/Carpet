#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <loopcontrol.h>
#include <vectors.h>

#include <cassert>
#include <cmath>
#include <cstdlib>
using std::ptrdiff_t;
using std::size_t;
using std::sqrt;



typedef vectype<CCTK_REAL> CCTK_VREAL;



void TestAlignedAlloc(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  const ptrdiff_t stride = CCTK_REAL_VEC_SIZE;
  assert(size_t(u) % sizeof(CCTK_REAL) == 0);
  const ptrdiff_t offset = size_t(u) / sizeof(CCTK_REAL) % stride;
  
  const ptrdiff_t di = 1;
  const ptrdiff_t dj = di * cctk_ash[0];
  const ptrdiff_t dk = dj * cctk_ash[1];
  
#pragma omp parallel
  CCTK_LOOP3STROFF_ALL(all3, cctkGH, i,j,k, imin,imax, stride,offset) {
    // const ptrdiff_t ind3d = CCTK_GFINDEX3D(cctkGH, i,j,k);
    const ptrdiff_t ind3d = i*di + j*dj + k*dk;
    const bool have_error = size_t(&u[ind3d]) % sizeof(CCTK_REAL_VEC) != 0;
    if (have_error) {
      CCTK_ERROR("Alignment error");
    }
    vec_store_partial_prepare(i, imin, imax);
    const CCTK_VREAL xl = CCTK_VREAL::load(x[ind3d]);
    const CCTK_VREAL yl = CCTK_VREAL::load(y[ind3d]);
    const CCTK_VREAL zl = CCTK_VREAL::load(z[ind3d]);
    const CCTK_VREAL ul = sqrt(xl*xl + yl*yl + zl*zl);
    vec_store_nta_partial(u[ind3d], ul);
  } CCTK_ENDLOOP3STROFF_ALL(all3);
  
#pragma omp parallel
  CCTK_LOOP3STROFF_INT(int3, cctkGH, i,j,k, imin,imax, stride,offset) {
    // const ptrdiff_t ind3d = CCTK_GFINDEX3D(cctkGH, i,j,k);
    const ptrdiff_t ind3d = i*di + j*dj + k*dk;
    const bool have_error = size_t(&u[ind3d]) % sizeof(CCTK_REAL_VEC) != 0;
    if (have_error) {
      CCTK_ERROR("Alignment error");
    }
    vec_store_partial_prepare(i, imin, imax);
    const CCTK_VREAL rl = vec_load(r[ind3d]);
    const CCTK_VREAL rm0l = vec_loadu_maybe3(-1,0,0,r[ind3d-di]);
    const CCTK_VREAL rp0l = vec_loadu_maybe3(+1,0,0,r[ind3d+di]);
    const CCTK_VREAL rm1l = vec_loadu_maybe3(0,-1,0,r[ind3d-dj]);
    const CCTK_VREAL rp1l = vec_loadu_maybe3(0,+1,0,r[ind3d+dj]);
    const CCTK_VREAL rm2l = vec_loadu_maybe3(0,0,-1,r[ind3d-dk]);
    const CCTK_VREAL rp2l = vec_loadu_maybe3(0,0,+1,r[ind3d+dk]);
    const CCTK_VREAL dd0r = rm0l - CCTK_VREAL(2.0)*rl + rp0l;
    const CCTK_VREAL dd1r = rm1l - CCTK_VREAL(2.0)*rl + rp1l;
    const CCTK_VREAL dd2r = rm2l - CCTK_VREAL(2.0)*rl + rp2l;
    const CCTK_VREAL ddul = dd0r + dd1r + dd2r;
    vec_store_nta_partial(ddu[ind3d], ddul);
  } CCTK_ENDLOOP3STROFF_INT(int3);
  
#pragma omp parallel
  CCTK_LOOP3STROFF_BND(bnd3, cctkGH, i,j,k, ni,nj,nk, imin,imax, stride,offset)
  {
    // const ptrdiff_t ind3d = CCTK_GFINDEX3D(cctkGH, i,j,k);
    const ptrdiff_t ind3d = i*di + j*dj + k*dk;
    const bool have_error = size_t(&u[ind3d]) % sizeof(CCTK_REAL_VEC) != 0;
    if (have_error) {
      CCTK_ERROR("Alignment error");
    }
    vec_store_partial_prepare(i, imin, imax);
    const CCTK_VREAL ddul = CCTK_VREAL(0.0);
    vec_store_nta_partial(ddu[ind3d], ddul);
  } CCTK_ENDLOOP3STROFF_BND(bnd3);
}
