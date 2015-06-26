#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <loopcontrol.h>

#include <cstdlib>
#include <vector>
using std::vector;
using std::size_t;



void TestLoopControlCXXVectorized(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  // We choose a stride of 4, so that with a grid size of 8 that is
  // defined in the test case parameter file, all j and k loops are
  // aligned as well.
  int const stride = 4;
  
  // Alignment in bytes
  int const alignment = stride * sizeof(CCTK_REAL);
  
  // We allocate a new grid function to ensure it is aligned. We
  // allocate "stride" more elements so that we can safely align
  // upwards.
  vector<CCTK_REAL> r2m(cctk_ash[0] * cctk_ash[1] * cctk_ash[2] + stride);
  // Align grid function
  CCTK_REAL *restrict const r2 =
    (CCTK_REAL*)(size_t(&r2m[stride]) / alignment * alignment);
  
#pragma omp parallel
  CCTK_LOOP3_ALL(copy, cctkGH, i,j,k) {
    int const ind3d = CCTK_GFINDEX3D(cctkGH, i,j,k);
    r2[ind3d] = r[ind3d];
  } CCTK_ENDLOOP3_ALL(copy);
  
  bool alignment_error = false;
  
  *cxxsum_all = 0.0;
#pragma omp parallel reduction(||: alignment_error)
  CCTK_LOOP3STR_ALL(loop3str_all, cctkGH, i,j,k, imin,imax, stride) {
    int const ind3d = CCTK_GFINDEX3D(cctkGH, imin,j,k) + i-imin;
    alignment_error = alignment_error || size_t(&r2[ind3d]) % alignment != 0;
    for (int v=0; v<stride; ++v) {
      bool const m = i+v>=imin && i+v<imax;
      if (m) {
        *cxxsum_all += r2[ind3d+v];
      }
    }
  } CCTK_ENDLOOP3STR_ALL(loop3str_all);
  if (alignment_error) {
    CCTK_ERROR("Alignment error");
  }
  
  *cxxsum_int = 0.0;
#pragma omp parallel
  CCTK_LOOP3STR_INT(loop3str_int, cctkGH, i,j,k, imin,imax, stride) {
    int const ind3d = CCTK_GFINDEX3D(cctkGH, imin,j,k) + i-imin;
    alignment_error = alignment_error || size_t(&r2[ind3d]) % alignment != 0;
    for (int v=0; v<stride; ++v) {
      bool const m = i+v>=imin && i+v<imax;
      if (m) {
        *cxxsum_int += r2[ind3d+v];
      }
    }
  } CCTK_ENDLOOP3STR_INT(loop3str_int);
  if (alignment_error) {
    CCTK_ERROR("Alignment error");
  }
  
  *cxxsum_bnd = 0.0;
#pragma omp parallel
  CCTK_LOOP3STR_BND(loop3str_bnd, cctkGH, i,j,k, ni,nj,nk, imin,imax, stride) {
    int const ind3d = CCTK_GFINDEX3D(cctkGH, imin,j,k) + i-imin;
    alignment_error = alignment_error || size_t(&r2[ind3d]) % alignment != 0;
    for (int v=0; v<stride; ++v) {
      bool const m = i+v>=imin && i+v<imax;
      if (m) {
        *cxxsum_bnd += r2[ind3d+v];
      }
    }
  } CCTK_ENDLOOP3STR_BND(loop3str_bnd);
  if (alignment_error) {
    CCTK_ERROR("Alignment error");
  }
  
  *cxxsum_intbnd = 0.0;
#pragma omp parallel
  CCTK_LOOP3STR_INTBND(loop3str_intbnd,
                       cctkGH, i,j,k, ni,nj,nk, imin,imax, stride)
  {
    int const ind3d = CCTK_GFINDEX3D(cctkGH, imin,j,k) + i-imin;
    alignment_error = alignment_error || size_t(&r2[ind3d]) % alignment != 0;
    for (int v=0; v<stride; ++v) {
      bool const m = i+v>=imin && i+v<imax;
      if (m) {
        *cxxsum_intbnd += r2[ind3d+v];
      }
    }
  } CCTK_ENDLOOP3STR_INTBND(loop3str_intbnd);
  if (alignment_error) {
    CCTK_ERROR("Alignment error");
  }
}



void TestLoopControlCXXVectorized2(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  // We choose a stride of 4, so that with a grid size of 8 that is
  // defined in the test case parameter file, all j and k loops are
  // aligned as well.
  int const stride = 4;
  
  // Alignment in bytes
  int const alignment = stride * sizeof(CCTK_REAL);
  
  // We allocate a new grid function to ensure it is aligned. We
  // allocate "stride" more elements so that we can safely align
  // upwards, and we align another "stride" more elements to ensure we
  // can access dummy elements before the first.
  vector<CCTK_REAL> r2m(cctk_ash[0] * cctk_ash[1] * cctk_ash[2] + 2*stride);
  // Offset of the first interior point
  int const int_off = CCTK_GFINDEX3D(cctkGH,
                                     cctk_nghostzones[0],
                                     cctk_nghostzones[1],
                                     cctk_nghostzones[2]);
  // Align grid function
  CCTK_REAL *restrict const r2 =
    (CCTK_REAL*)(size_t(&r2m[int_off + 2*stride]) / alignment * alignment) -
    int_off;
  
#pragma omp parallel
  CCTK_LOOP3_ALL(copy, cctkGH, i,j,k) {
    int const ind3d = CCTK_GFINDEX3D(cctkGH, i,j,k);
    r2[ind3d] = r[ind3d];
  } CCTK_ENDLOOP3_ALL(copy);
  
  int const offset = size_t(r2) / sizeof(CCTK_REAL) % stride;
  
  bool alignment_error = false;
  
  *cxx2sum_all = 0.0;
#pragma omp parallel reduction(||: alignment_error)
  CCTK_LOOP3STRMOD_ALL(loop3strmod_all,
                       cctkGH, i,j,k, imin,imax, stride,offset)
  {
    int const ind3d = CCTK_GFINDEX3D(cctkGH, 0,j,k) + i;
    alignment_error = alignment_error || size_t(&r2[ind3d]) % alignment != 0;
    for (int v=0; v<stride; ++v) {
      bool const m = i+v>=imin && i+v<imax;
      if (m) {
        *cxx2sum_all += r2[ind3d+v];
      }
    }
  } CCTK_ENDLOOP3STRMOD_ALL(loop3strmod_all);
  if (alignment_error) {
    CCTK_ERROR("Alignment error");
  }
  
  *cxx2sum_int = 0.0;
#pragma omp parallel reduction(||: alignment_error)
  CCTK_LOOP3STRMOD_INT(loop3strmod_int,
                       cctkGH, i,j,k, imin,imax, stride,offset)
  {
    int const ind3d = CCTK_GFINDEX3D(cctkGH, 0,j,k) + i;
    alignment_error = alignment_error || size_t(&r2[ind3d]) % alignment != 0;
    for (int v=0; v<stride; ++v) {
      bool const m = i+v>=imin && i+v<imax;
      if (m) {
        *cxx2sum_int += r2[ind3d+v];
      }
    }
  } CCTK_ENDLOOP3STRMOD_INT(loop3strmod_int);
  if (alignment_error) {
    CCTK_ERROR("Alignment error");
  }
  
  *cxx2sum_bnd = 0.0;
#pragma omp parallel reduction(||: alignment_error)
  CCTK_LOOP3STRMOD_BND(loop3strmod_bnd,
                       cctkGH, i,j,k, ni,nj,nk, imin,imax, stride,offset)
  {
    int const ind3d = CCTK_GFINDEX3D(cctkGH, 0,j,k) + i;
    alignment_error = alignment_error || size_t(&r2[ind3d]) % alignment != 0;
    for (int v=0; v<stride; ++v) {
      bool const m = i+v>=imin && i+v<imax;
      if (m) {
        *cxx2sum_bnd += r2[ind3d+v];
      }
    }
  } CCTK_ENDLOOP3STRMOD_BND(loop3strmod_bnd);
  if (alignment_error) {
    CCTK_ERROR("Alignment error");
  }
  
  *cxx2sum_intbnd = 0.0;
#pragma omp parallel reduction(||: alignment_error)
  CCTK_LOOP3STRMOD_INTBND(loop3strmod_intbnd,
                          cctkGH, i,j,k, ni,nj,nk, imin,imax, stride,offset)
  {
    int const ind3d = CCTK_GFINDEX3D(cctkGH, 0,j,k) + i;
    alignment_error = alignment_error || size_t(&r2[ind3d]) % alignment != 0;
    for (int v=0; v<stride; ++v) {
      bool const m = i+v>=imin && i+v<imax;
      if (m) {
        *cxx2sum_intbnd += r2[ind3d+v];
      }
    }
  } CCTK_ENDLOOP3STRMOD_INTBND(loop3strmod_intbnd);
  if (alignment_error) {
    CCTK_ERROR("Alignment error");
  }
}
