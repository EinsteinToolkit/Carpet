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
  // allocate more elements so that we can safely align upwards, and
  // can access additional elements past the end.
  size_t const npoints = cctk_ash[0] * cctk_ash[1] * cctk_ash[2];
  size_t const last_index =
    CCTK_GFINDEX3D(cctkGH, cctk_lsh[0]-1, cctk_lsh[1]-1, cctk_lsh[2]-1);
  vector<unsigned char>
    r2m(sizeof(CCTK_REAL) * (npoints + stride-1) + alignment-1);
  size_t const r2ptr = size_t(&r2m[0]);
  // Align grid function
  size_t const r2ptr_aligned = (r2ptr + alignment-1) / alignment * alignment;
  assert(r2ptr_aligned >= r2ptr);
  assert(r2ptr_aligned < r2ptr + alignment);
  assert(r2ptr_aligned % alignment == 0);
  // Convert pointer
  CCTK_REAL *restrict const r2 = (CCTK_REAL*)r2ptr_aligned;
  assert((unsigned char*)&r2[0] >= &r2m.begin()[0]);
  assert((unsigned char*)&r2[last_index + stride] <= &r2m.end()[0]);
  assert(size_t(&r2[0]) % alignment == 0);
  
#pragma omp parallel
  CCTK_LOOP3_ALL(copy, cctkGH, i,j,k) {
    int const ind3d = CCTK_GFINDEX3D(cctkGH, i,j,k);
    r2[ind3d] = r[ind3d];
  } CCTK_ENDLOOP3_ALL(copy);
  
  CCTK_REAL sum_all = 0.0;
#pragma omp parallel reduction(+: sum_all)
  CCTK_LOOP3STR_ALL(loop3str_all, cctkGH, i,j,k, imin,imax, stride) {
    int const ind3d = CCTK_GFINDEX3D(cctkGH, imin,j,k) + i-imin;
    bool const alignment_error = size_t(&r2[ind3d]) % alignment != 0;
    if (alignment_error) {
      CCTK_ERROR("Alignment error");
    }
    bool any_m = false;
    for (int v=0; v<stride; ++v) {
      bool const m = i+v>=imin && i+v<imax;
      if (m) {
        any_m = true;
        sum_all += r2[ind3d+v];
      }
    }
    assert(any_m);
  } CCTK_ENDLOOP3STR_ALL(loop3str_all);
  *cxxsum_all = sum_all;
  
  CCTK_REAL sum_int = 0.0;
#pragma omp parallel reduction(+: sum_int)
  CCTK_LOOP3STR_INT(loop3str_int, cctkGH, i,j,k, imin,imax, stride) {
    int const ind3d = CCTK_GFINDEX3D(cctkGH, imin,j,k) + i-imin;
    bool const alignment_error = size_t(&r2[ind3d]) % alignment != 0;
    if (alignment_error) {
      CCTK_ERROR("Alignment error");
    }
    bool any_m = false;
    for (int v=0; v<stride; ++v) {
      bool const m = i+v>=imin && i+v<imax;
      if (m) {
        any_m = true;
        sum_int += r2[ind3d+v];
      }
    }
    assert(any_m);
  } CCTK_ENDLOOP3STR_INT(loop3str_int);
  *cxxsum_int = sum_int;
  
  CCTK_REAL sum_bnd = 0.0;
#pragma omp parallel reduction(+: sum_bnd)
  CCTK_LOOP3STR_BND(loop3str_bnd, cctkGH, i,j,k, ni,nj,nk, imin,imax, stride) {
    int const ind3d = CCTK_GFINDEX3D(cctkGH, imin,j,k) + i-imin;
    bool const alignment_error = size_t(&r2[ind3d]) % alignment != 0;
    if (alignment_error) {
      CCTK_ERROR("Alignment error");
    }
    bool any_m = false;
    for (int v=0; v<stride; ++v) {
      bool const m = i+v>=imin && i+v<imax;
      if (m) {
        any_m = true;
        sum_bnd += r2[ind3d+v];
      }
    }
    assert(any_m);
  } CCTK_ENDLOOP3STR_BND(loop3str_bnd);
  *cxxsum_bnd = sum_bnd;
  
  CCTK_REAL sum_intbnd = 0.0;
#pragma omp parallel reduction(+: sum_intbnd)
  CCTK_LOOP3STR_INTBND(loop3str_intbnd,
                       cctkGH, i,j,k, ni,nj,nk, imin,imax, stride)
  {
    int const ind3d = CCTK_GFINDEX3D(cctkGH, imin,j,k) + i-imin;
    bool const alignment_error = size_t(&r2[ind3d]) % alignment != 0;
    if (alignment_error) {
      CCTK_ERROR("Alignment error");
    }
    bool any_m = false;
    for (int v=0; v<stride; ++v) {
      bool const m = i+v>=imin && i+v<imax;
      if (m) {
        any_m = true;
        sum_intbnd += r2[ind3d+v];
      }
    }
    assert(any_m);
  } CCTK_ENDLOOP3STR_INTBND(loop3str_intbnd);
  *cxxsum_intbnd = sum_intbnd;
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
  // allocate more elements so that we can safely align upwards, and
  // can access additional elements before the beginning and past the
  // end.
  size_t const npoints = cctk_ash[0] * cctk_ash[1] * cctk_ash[2];
  size_t const interior_index = CCTK_GFINDEX3D(cctkGH,
                                               cctk_nghostzones[0],
                                               cctk_nghostzones[1],
                                               cctk_nghostzones[2]);
  size_t const last_index =
    CCTK_GFINDEX3D(cctkGH, cctk_lsh[0]-1, cctk_lsh[1]-1, cctk_lsh[2]-1);
  vector<unsigned char>
    r2m(sizeof(CCTK_REAL) * (npoints + 2*(stride-1)) + alignment-1);
  size_t const r2ptr = size_t(&r2m[0]);
  // Offset of the first interior point
  size_t const int_off = sizeof(CCTK_REAL) * interior_index;
  // Add space before first grid point
  size_t const r2ptr2 = r2ptr + sizeof(CCTK_REAL) * (stride-1);
  // Align grid function
  size_t const r2ptr_aligned =
    (r2ptr2 + alignment-1 + int_off) / alignment * alignment - int_off;
  assert(r2ptr_aligned >= r2ptr2);
  assert(r2ptr_aligned < r2ptr2 + alignment);
  assert((r2ptr_aligned + int_off) % alignment == 0);
  // Convert pointer
  CCTK_REAL *restrict const r2 = (CCTK_REAL*)r2ptr_aligned;
  assert((unsigned char*)&r2[-(stride-1)] >= &r2m.begin()[0]);
  assert((unsigned char*)&r2[last_index + stride] <= &r2m.end()[0]);
  assert(size_t(&r2[interior_index]) % alignment == 0);
  
#pragma omp parallel
  CCTK_LOOP3_ALL(copy, cctkGH, i,j,k) {
    int const ind3d = CCTK_GFINDEX3D(cctkGH, i,j,k);
    r2[ind3d] = r[ind3d];
  } CCTK_ENDLOOP3_ALL(copy);
  
  int const offset = size_t(r2) / sizeof(CCTK_REAL) % stride;
  
  CCTK_REAL sum_all = 0.0;
#pragma omp parallel reduction(+: sum_all)
  CCTK_LOOP3STROFF_ALL(loop3stroff_all,
                       cctkGH, i,j,k, imin,imax, stride,offset)
  {
    int const ind3d = CCTK_GFINDEX3D(cctkGH, 0,j,k) + i;
    bool const alignment_error = size_t(&r2[ind3d]) % alignment != 0;
    if (alignment_error) {
      CCTK_ERROR("Alignment error");
    }
    bool any_m = false;
    for (int v=0; v<stride; ++v) {
      bool const m = i+v>=imin && i+v<imax;
      if (m) {
        any_m = true;
        sum_all += r2[ind3d+v];
      }
    }
    assert(any_m);
  } CCTK_ENDLOOP3STROFF_ALL(loop3stroff_all);
  *cxx2sum_all = sum_all;
  
  CCTK_REAL sum_int = 0.0;
#pragma omp parallel reduction(+: sum_int)
  CCTK_LOOP3STROFF_INT(loop3stroff_int,
                       cctkGH, i,j,k, imin,imax, stride,offset)
  {
    int const ind3d = CCTK_GFINDEX3D(cctkGH, 0,j,k) + i;
    bool const alignment_error = size_t(&r2[ind3d]) % alignment != 0;
    if (alignment_error) {
      CCTK_ERROR("Alignment error");
    }
    bool any_m = false;
    for (int v=0; v<stride; ++v) {
      bool const m = i+v>=imin && i+v<imax;
      if (m) {
        any_m = true;
        sum_int += r2[ind3d+v];
      }
    }
    assert(any_m);
  } CCTK_ENDLOOP3STROFF_INT(loop3stroff_int);
  *cxx2sum_int = sum_int;
  
  CCTK_REAL sum_bnd = 0.0;
#pragma omp parallel reduction(+: sum_bnd)
  CCTK_LOOP3STROFF_BND(loop3stroff_bnd,
                       cctkGH, i,j,k, ni,nj,nk, imin,imax, stride,offset)
  {
    int const ind3d = CCTK_GFINDEX3D(cctkGH, 0,j,k) + i;
    bool const alignment_error = size_t(&r2[ind3d]) % alignment != 0;
    if (alignment_error) {
      CCTK_ERROR("Alignment error");
    }
    bool any_m = false;
    for (int v=0; v<stride; ++v) {
      bool const m = i+v>=imin && i+v<imax;
      if (m) {
        any_m = true;
        sum_bnd += r2[ind3d+v];
      }
    }
    assert(any_m);
  } CCTK_ENDLOOP3STROFF_BND(loop3stroff_bnd);
  *cxx2sum_bnd = sum_bnd;
  
  CCTK_REAL sum_intbnd = 0.0;
#pragma omp parallel reduction(+: sum_intbnd)
  CCTK_LOOP3STROFF_INTBND(loop3stroff_intbnd,
                          cctkGH, i,j,k, ni,nj,nk, imin,imax, stride,offset)
  {
    int const ind3d = CCTK_GFINDEX3D(cctkGH, 0,j,k) + i;
    bool const alignment_error = size_t(&r2[ind3d]) % alignment != 0;
    if (alignment_error) {
      CCTK_ERROR("Alignment error");
    }
    bool any_m = false;
    for (int v=0; v<stride; ++v) {
      bool const m = i+v>=imin && i+v<imax;
      if (m) {
        any_m = true;
        sum_intbnd += r2[ind3d+v];
      }
    }
    assert(any_m);
  } CCTK_ENDLOOP3STROFF_INTBND(loop3stroff_intbnd);
  *cxx2sum_intbnd = sum_intbnd;
}
