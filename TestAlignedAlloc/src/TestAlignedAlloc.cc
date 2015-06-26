#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <loopcontrol.h>
#include <vectors.h>

#include <cassert>
#include <cstdlib>
using std::ptrdiff_t;
using std::size_t;



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
  
  bool have_error;
  
  have_error = false;
#pragma omp parallel reduction(||: have_error)
  CCTK_LOOP3STRMOD_ALL(all3, cctkGH, i,j,k, imin,imax, stride,offset) {
    // const ptrdiff_t ind3d = CCTK_GFINDEX3D(cctkGH, i,j,k);
    const ptrdiff_t ind3d = i*di + j*dj + k*dk;
    have_error = have_error ||
      ptrdiff_t(&u[ind3d]) % sizeof(CCTK_REAL_VEC) != 0;
    vec_store_partial_prepare(i, imin, imax);
    const CCTK_REAL_VEC xl = vec_load(x[ind3d]);
    const CCTK_REAL_VEC yl = vec_load(y[ind3d]);
    const CCTK_REAL_VEC zl = vec_load(z[ind3d]);
    const CCTK_REAL_VEC ul = ksqrt(xl*xl + yl*yl + zl*zl);
    vec_store_nta_partial(u[ind3d], ul);
  } CCTK_ENDLOOP3STRMOD_ALL(all3);
  if (have_error) {
    CCTK_ERROR("Alignment error");
  }
  
  have_error = false;
#pragma omp parallel reduction(||: have_error)
  CCTK_LOOP3STRMOD_INT(int3, cctkGH, i,j,k, imin,imax, stride,offset) {
    // const ptrdiff_t ind3d = CCTK_GFINDEX3D(cctkGH, i,j,k);
    const ptrdiff_t ind3d = i*di + j*dj + k*dk;
    have_error = have_error ||
      ptrdiff_t(&u[ind3d]) % sizeof(CCTK_REAL_VEC) != 0;
    assert(!have_error);
    vec_store_partial_prepare(i, imin, imax);
    const CCTK_REAL_VEC rl = vec_load(r[ind3d]);
    const CCTK_REAL_VEC rm0l = vec_loadu_maybe3(-1,0,0,r[ind3d-di]);
    const CCTK_REAL_VEC rp0l = vec_loadu_maybe3(+1,0,0,r[ind3d+di]);
    const CCTK_REAL_VEC rm1l = vec_loadu_maybe3(0,-1,0,r[ind3d-dj]);
    const CCTK_REAL_VEC rp1l = vec_loadu_maybe3(0,+1,0,r[ind3d+dj]);
    const CCTK_REAL_VEC rm2l = vec_loadu_maybe3(0,0,-1,r[ind3d-dk]);
    const CCTK_REAL_VEC rp2l = vec_loadu_maybe3(0,0,+1,r[ind3d+dk]);
    const CCTK_REAL_VEC dd0r = rm0l - 2*rl + rp0l;
    const CCTK_REAL_VEC dd1r = rm1l - 2*rl + rp1l;
    const CCTK_REAL_VEC dd2r = rm2l - 2*rl + rp2l;
    const CCTK_REAL_VEC ddul = dd0r + dd1r + dd2r;
    vec_store_nta_partial(ddu[ind3d], ddul);
  } CCTK_ENDLOOP3STRMOD_INT(int3);
  if (have_error) {
    CCTK_ERROR("Alignment error");
  }
  
  have_error = false;
#pragma omp parallel reduction(||: have_error)
  CCTK_LOOP3STRMOD_BND(bnd3, cctkGH, i,j,k, ni,nj,nk, imin,imax, stride,offset)
  {
    // const ptrdiff_t ind3d = CCTK_GFINDEX3D(cctkGH, i,j,k);
    const ptrdiff_t ind3d = i*di + j*dj + k*dk;
    have_error = have_error ||
      ptrdiff_t(&u[ind3d]) % sizeof(CCTK_REAL_VEC) != 0;
    assert(!have_error);
    vec_store_partial_prepare(i, imin, imax);
    const CCTK_REAL_VEC ddul = vec_set1(0);
    vec_store_nta_partial(ddu[ind3d], ddul);
  } CCTK_ENDLOOP3STRMOD_BND(bnd3);
  if (have_error) {
    CCTK_ERROR("Alignment error");
  }
}
