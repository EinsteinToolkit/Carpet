#include <cctk.h>
#include <cctk_Arguments.h>

#include <vectors.h>

#include <loopcontrol.h>

void lc_selftest (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  
#ifdef CCTK_REAL_VEC_SIZE
  /* Vectorisation is enabled */
  int const vector_size = CCTK_REAL_VEC_SIZE;
#else
  /* Vectorisation is disabled */
  int const vector_size = 1;
#endif
  
  /* Initialise (without using LoopControl) */
#pragma omp parallel for
  for (int k=0; k<cctk_lsh[2]; ++k) {
    for (int j=0; j<cctk_lsh[1]; ++j) {
      for (int i=0; i<cctk_lsh[0]; ++i) {
        int const ind3d = CCTK_GFINDEX3D(cctkGH, i,j,k);
        var1[ind3d] = 0;
        var2[ind3d] = 0;
      }
    }
  }
  
  /* Test 1: Loop over all points */
  {
    int imin[3], imax[3];
    for (int d=0; d<3; ++d) {
      imin[d] = 0;
      imax[d] = cctk_lsh[d];
    }
#pragma omp parallel
    LC_LOOP3VEC(lc_selftest1,
                i,j,k,
                imin[0],imin[1],imin[2],
                imax[0],imax[1],imax[2],
                cctk_lsh[0],cctk_lsh[1],cctk_lsh[2],
                vector_size)
    {
      int const ind3d = CCTK_GFINDEX3D(cctkGH, i,j,k);
#pragma omp atomic
      ++ var1[ind3d];
    } LC_ENDLOOP3VEC(lc_selftest1);
  }
  
  /* Test 2: Loop over interior points, then loop over boundary
     points */
  {
    int imin[3], imax[3];
    for (int d=0; d<3; ++d) {
      imin[d] = cctk_nghostzones[d];
      imax[d] = cctk_lsh[d] - cctk_nghostzones[d];
    }
#pragma omp parallel
    LC_LOOP3VEC(lc_selftest2a,
                i,j,k,
                imin[0],imin[1],imin[2],
                imax[0],imax[1],imax[2],
                cctk_lsh[0],cctk_lsh[1],cctk_lsh[2],
                vector_size)
    {
      int const ind3d = CCTK_GFINDEX3D(cctkGH, i,j,k);
#pragma omp atomic
      ++ var2[ind3d];
    } LC_ENDLOOP3VEC(lc_selftest2a);
    
    for (int dir=0; dir<3; ++dir) {
      for (int face=0; face<2; ++face) {
        for (int d=0; d<dir; ++d) {
          imin[d] = 0;
          imax[d] = cctk_lsh[d];
        }
        if (face==0) {
          imin[dir] = 0;
          imax[dir] = cctk_nghostzones[dir];
        } else {
          imin[dir] = cctk_lsh[dir] - cctk_nghostzones[dir];
          imax[dir] = cctk_lsh[dir];
        }
        for (int d=dir+1; d<3; ++d) {
          imin[d] = cctk_nghostzones[d];
          imax[d] = cctk_lsh[d] - cctk_nghostzones[d];
        }
#pragma omp parallel
        LC_LOOP3VEC(lc_selftest2b,
                    i,j,k,
                    imin[0],imin[1],imin[2],
                    imax[0],imax[1],imax[2],
                    cctk_lsh[0],cctk_lsh[1],cctk_lsh[2],
                    vector_size)
        {
          int const ind3d = CCTK_GFINDEX3D(cctkGH, i,j,k);
#pragma omp atomic
          ++ var2[ind3d];
        } LC_ENDLOOP3VEC(lc_selftest2b);
      }
    }
  }
  
  /* Evaluate tests (without using LoopControl) */
  int failure = 0;
#pragma omp parallel for reduction(+: failure)
  for (int k=0; k<cctk_lsh[2]; ++k) {
    for (int j=0; j<cctk_lsh[1]; ++j) {
      for (int i=0; i<cctk_lsh[0]; ++i) {
        int const ind3d = CCTK_GFINDEX3D(cctkGH, i,j,k);
        failure += var1[ind3d] != 1;
        failure += var2[ind3d] != 1;
      }
    }
  }
  
  if (failure) {
    CCTK_WARN (CCTK_WARN_ABORT, "LoopControl self-test failed");
  }
}
