#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <loopcontrol.h>

void TestLoopControlC(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL csum;

  switch (cctk_dim) {

  case 1: {
    /* 1D */

    csum = 0.;
#pragma omp parallel reduction(+: csum)
    CCTK_LOOP1_ALL(loop1_all, cctkGH, i) {
      int const ind1d = CCTK_GFINDEX1D(cctkGH, i);
      csum += r[ind1d];
    } CCTK_ENDLOOP1_ALL(loop1_all);
    *csum_all = csum;

    csum = 0.;
#pragma omp parallel reduction(+: csum)
    CCTK_LOOP1_INT(loop1_int, cctkGH, i) {
      int const ind1d = CCTK_GFINDEX1D(cctkGH, i);
      csum += r[ind1d];
    } CCTK_ENDLOOP1_INT(loop1_int);
    *csum_int = csum;

    csum = 0.;
#pragma omp parallel reduction(+: csum)
    CCTK_LOOP1_BND(loop1_bnd, cctkGH, i, ni) {
      int const ind1d = CCTK_GFINDEX1D(cctkGH, i);
      csum += r[ind1d];
    } CCTK_ENDLOOP1_BND(loop1_bnd);
    *csum_bnd = csum;

    csum = 0.;
#pragma omp parallel reduction(+: csum)
    CCTK_LOOP1_INTBND(loop1_intbnd, cctkGH, i, ni) {
      int const ind1d = CCTK_GFINDEX1D(cctkGH, i);
      csum += r[ind1d];
    } CCTK_ENDLOOP1_INTBND(loop1_intbnd);
    *csum_intbnd = csum;

    break;
  }



  case 2:{
    /* 2D */

    csum = 0.;
#pragma omp parallel reduction(+: csum)
    CCTK_LOOP2_ALL(loop2_all, cctkGH, i,j) {
      int const ind2d = CCTK_GFINDEX2D(cctkGH, i,j);
      csum += r[ind2d];
    } CCTK_ENDLOOP2_ALL(loop2_all);
    *csum_all = csum;

    csum = 0.;
#pragma omp parallel reduction(+: csum)
    CCTK_LOOP2_INT(loop2_int, cctkGH, i,j) {
      int const ind2d = CCTK_GFINDEX2D(cctkGH, i,j);
      csum += r[ind2d];
    } CCTK_ENDLOOP2_INT(loop2_int);
    *csum_int = csum;

    csum = 0.;
#pragma omp parallel reduction(+: csum)
    CCTK_LOOP2_BND(loop2_bnd, cctkGH, i,j, ni,nj) {
      int const ind2d = CCTK_GFINDEX2D(cctkGH, i,j);
      csum += r[ind2d];
    } CCTK_ENDLOOP2_BND(loop2_bnd);
    *csum_bnd = csum;

    csum = 0.;
#pragma omp parallel reduction(+: csum)
    CCTK_LOOP2_INTBND(loop2_intbnd, cctkGH, i,j, ni,nj) {
      int const ind2d = CCTK_GFINDEX2D(cctkGH, i,j);
      csum += r[ind2d];
    } CCTK_ENDLOOP2_INTBND(loop2_intbnd);
    *csum_intbnd = csum;

    break;
  }



  case 3: {
    /* 3D */

    csum = 0.;
#pragma omp parallel reduction(+: csum)
    CCTK_LOOP3_ALL(loop3_all, cctkGH, i,j,k) {
      int const ind3d = CCTK_GFINDEX3D(cctkGH, i,j,k);
      csum += r[ind3d];
    } CCTK_ENDLOOP3_ALL(loop3_all);
    *csum_all = csum;

    csum = 0.;
#pragma omp parallel reduction(+: csum)
    CCTK_LOOP3_INT(loop3_int, cctkGH, i,j,k) {
      int const ind3d = CCTK_GFINDEX3D(cctkGH, i,j,k);
      csum += r[ind3d];
    } CCTK_ENDLOOP3_INT(loop3_int);
    *csum_int = csum;

    csum = 0.;
#pragma omp parallel reduction(+: csum)
    CCTK_LOOP3_BND(loop3_bnd, cctkGH, i,j,k, ni,nj,nk) {
      int const ind3d = CCTK_GFINDEX3D(cctkGH, i,j,k);
      csum += r[ind3d];
    } CCTK_ENDLOOP3_BND(loop3_bnd);
    *csum_bnd = csum;

    csum = 0.;
#pragma omp parallel reduction(+: csum)
    CCTK_LOOP3_INTBND(loop3_intbnd, cctkGH, i,j,k, ni,nj,nk) {
      int const ind3d = CCTK_GFINDEX3D(cctkGH, i,j,k);
      csum += r[ind3d];
    } CCTK_ENDLOOP3_INTBND(loop3_intbnd);
    *csum_intbnd = csum;

    break;
  }



  case 4: {
    /* 4D */

    csum = 0.;
#pragma omp parallel reduction(+: csum)
    CCTK_LOOP4_ALL(loop4_all, cctkGH, i,j,k,l) {
      int const ind4d = CCTK_GFINDEX4D(cctkGH, i,j,k,l);
      csum += r[ind4d];
    } CCTK_ENDLOOP4_ALL(loop4_all);
    *csum_all = csum;

    csum = 0.;
#pragma omp parallel reduction(+: csum)
    CCTK_LOOP4_INT(loop4_int, cctkGH, i,j,k,l) {
      int const ind4d = CCTK_GFINDEX4D(cctkGH, i,j,k,l);
      csum += r[ind4d];
    } CCTK_ENDLOOP4_INT(loop4_int);
    *csum_int = csum;

    csum = 0.;
#pragma omp parallel reduction(+: csum)
    CCTK_LOOP4_BND(loop4_bnd, cctkGH, i,j,k,l, ni,nj,nk,nl) {
      int const ind4d = CCTK_GFINDEX4D(cctkGH, i,j,k,l);
      csum += r[ind4d];
    } CCTK_ENDLOOP4_BND(loop4_bnd);
    *csum_bnd = csum;

    csum = 0.;
#pragma omp parallel reduction(+: csum)
    CCTK_LOOP4_INTBND(loop4_intbnd, cctkGH, i,j,k,l, ni,nj,nk,nl) {
      int const ind4d = CCTK_GFINDEX4D(cctkGH, i,j,k,l);
      csum += r[ind4d];
    } CCTK_ENDLOOP4_INTBND(loop4_intbnd);
    *csum_intbnd = csum;

    break;
  }



  default:
    CCTK_ERROR("cctk_dim out of range");

  }
}
