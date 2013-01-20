#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <stdio.h>

#include <loopcontrol.h>



static void TestLoopControlPointwise_All(CCTK_ARGUMENTS);
static void TestLoopControlPointwise_Int(CCTK_ARGUMENTS);
static void TestLoopControlPointwise_Bnd(CCTK_ARGUMENTS);
static void TestLoopControlPointwise_IntBnd(CCTK_ARGUMENTS);



void TestLoopControlPointwise_All(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  for (int k=0; k<cctk_lsh[2]; ++k) {
    for (int j=0; j<cctk_lsh[1]; ++j) {
      for (int i=0; i<cctk_lsh[0]; ++i) {
        int const ind3d = CCTK_GFINDEX3D(cctkGH, i,j,k);
        pointtypes[ind3d] = 1;
      }
    }
  }
  
#pragma omp parallel
  CCTK_LOOP3_ALL(loop3_all, cctkGH, i,j,k) {
    int const ind3d = CCTK_GFINDEX3D(cctkGH, i,j,k);
    pointtypes[ind3d] -= 1;
  } CCTK_ENDLOOP3_ALL(loop3_all);
  
  int num_errors = 0;
  for (int k=0; k<cctk_lsh[2]; ++k) {
    for (int j=0; j<cctk_lsh[1]; ++j) {
      for (int i=0; i<cctk_lsh[0]; ++i) {
        int const ind3d = CCTK_GFINDEX3D(cctkGH, i,j,k);
        num_errors += pointtypes[ind3d] != 0;
      }
    }
  }
  if (num_errors > 0) {
    CCTK_WARN(CCTK_WARN_ABORT, "TestLoopControlPointwise_All failed");
  }
}



void TestLoopControlPointwise_Int(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT bndsize    [6];
  CCTK_INT is_ghostbnd[6];
  CCTK_INT is_symbnd  [6];
  CCTK_INT is_physbnd [6];
  GetBoundarySizesAndTypes
    (cctkGH, 6, bndsize, is_ghostbnd, is_symbnd, is_physbnd);
  
  int imin[3], imax[3];
  for (int d=0; d<3; ++d) {
    imin[d] = cctk_bbox[2*d  ] ? bndsize[2*d  ] : cctk_nghostzones[d];
    imax[d] = (cctk_lsh[d] -
               (cctk_bbox[2*d+1] ? bndsize[2*d+1] : cctk_nghostzones[d]));
  }
  
  for (int k=0; k<cctk_lsh[2]; ++k) {
    for (int j=0; j<cctk_lsh[1]; ++j) {
      for (int i=0; i<cctk_lsh[0]; ++i) {
        int const ind3d = CCTK_GFINDEX3D(cctkGH, i,j,k);
        pointtypes[ind3d] = 0;
      }
    }
  }
  
  for (int k=imin[2]; k<imax[2]; ++k) {
    for (int j=imin[1]; j<imax[1]; ++j) {
      for (int i=imin[0]; i<imax[0]; ++i) {
        int const ind3d = CCTK_GFINDEX3D(cctkGH, i,j,k);
        pointtypes[ind3d] = 1;
      }
    }
  }
  
#pragma omp parallel
  CCTK_LOOP3_INT(loop3_int, cctkGH, i,j,k) {
    int const ind3d = CCTK_GFINDEX3D(cctkGH, i,j,k);
    pointtypes[ind3d] -= 1;
  } CCTK_ENDLOOP3_INT(loop3_int);
  
  int num_errors = 0;
  for (int k=0; k<cctk_lsh[2]; ++k) {
    for (int j=0; j<cctk_lsh[1]; ++j) {
      for (int i=0; i<cctk_lsh[0]; ++i) {
        int const ind3d = CCTK_GFINDEX3D(cctkGH, i,j,k);
        num_errors += pointtypes[ind3d] != 0;
      }
    }
  }
  if (num_errors > 0) {
    CCTK_WARN(CCTK_WARN_ABORT, "TestLoopControlPointwise_Int failed");
  }
}



void TestLoopControlPointwise_Bnd(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT bndsize    [6];
  CCTK_INT is_ghostbnd[6];
  CCTK_INT is_symbnd  [6];
  CCTK_INT is_physbnd [6];
  GetBoundarySizesAndTypes
    (cctkGH, 6, bndsize, is_ghostbnd, is_symbnd, is_physbnd);
  
  for (int k=0; k<cctk_lsh[2]; ++k) {
    for (int j=0; j<cctk_lsh[1]; ++j) {
      for (int i=0; i<cctk_lsh[0]; ++i) {
        int const ind3d = CCTK_GFINDEX3D(cctkGH, i,j,k);
        pointtypes[ind3d] = 0;
      }
    }
  }
  
  for (int dir=0; dir<3; ++dir) {
    for (int face=0; face<2; ++face) {
      if (is_physbnd[2*dir+face]) {
        int imin[3], imax[3];
        for (int d=0; d<3; ++d) {
          imin[d] = 0;
          imax[d] = cctk_lsh[d];
        }
        if (face==0) {
          imax[dir] = bndsize[2*dir];
        } else {
          imin[dir] = cctk_lsh[dir] - bndsize[2*dir+1];
        }
        for (int k=imin[2]; k<imax[2]; ++k) {
          for (int j=imin[1]; j<imax[1]; ++j) {
            for (int i=imin[0]; i<imax[0]; ++i) {
              int const ind3d = CCTK_GFINDEX3D(cctkGH, i,j,k);
              pointtypes[ind3d] = 1;
            }
          }
        }
      }
    }
  }
  
#pragma omp parallel
  CCTK_LOOP3_BND(loop3_bnd, cctkGH, i,j,k, ni,nj,nk) {
    int const ind3d = CCTK_GFINDEX3D(cctkGH, i,j,k);
    pointtypes[ind3d] -= 1;
  } CCTK_ENDLOOP3_BND(loop3_bnd);
  
  int num_errors = 0;
  for (int k=0; k<cctk_lsh[2]; ++k) {
    for (int j=0; j<cctk_lsh[1]; ++j) {
      for (int i=0; i<cctk_lsh[0]; ++i) {
        int const ind3d = CCTK_GFINDEX3D(cctkGH, i,j,k);
        num_errors += pointtypes[ind3d] != 0;
        if (pointtypes[ind3d] != 0) {
          printf("[%d %d %d] %d\n", i,j,k, (int)pointtypes[ind3d]);
        }
      }
    }
  }
  if (num_errors > 0) {
    CCTK_WARN(CCTK_WARN_ABORT, "TestLoopControlPointwise_Bnd failed");
  }
}



void TestLoopControlPointwise_IntBnd(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT bndsize    [6];
  CCTK_INT is_ghostbnd[6];
  CCTK_INT is_symbnd  [6];
  CCTK_INT is_physbnd [6];
  GetBoundarySizesAndTypes
    (cctkGH, 6, bndsize, is_ghostbnd, is_symbnd, is_physbnd);
  
  for (int k=0; k<cctk_lsh[2]; ++k) {
    for (int j=0; j<cctk_lsh[1]; ++j) {
      for (int i=0; i<cctk_lsh[0]; ++i) {
        int const ind3d = CCTK_GFINDEX3D(cctkGH, i,j,k);
        pointtypes[ind3d] = 0;
      }
    }
  }
  
  for (int dir=0; dir<3; ++dir) {
    for (int face=0; face<2; ++face) {
      if (is_physbnd[2*dir+face] &&
          !is_ghostbnd[2*dir+face] && !is_symbnd[2*dir+face])
      {
        int imin[3], imax[3];
        for (int d=0; d<3; ++d) {
          imin[d] = 0;
          imax[d] = cctk_lsh[d];
        }
        if (face==0) {
          imax[dir] = bndsize[2*dir];
        } else {
          imin[dir] = cctk_lsh[dir] - bndsize[2*dir+1];
        }
        for (int k=imin[2]; k<imax[2]; ++k) {
          for (int j=imin[1]; j<imax[1]; ++j) {
            for (int i=imin[0]; i<imax[0]; ++i) {
              int const ind3d = CCTK_GFINDEX3D(cctkGH, i,j,k);
              pointtypes[ind3d] = 1;
            }
          }
        }
      }
    }
  }
  
#pragma omp parallel
  CCTK_LOOP3_INTBND(loop3_intbnd, cctkGH, i,j,k, ni,nj,nk) {
    int const ind3d = CCTK_GFINDEX3D(cctkGH, i,j,k);
    pointtypes[ind3d] += 1;
  } CCTK_ENDLOOP3_INTBND(loop3_intbnd);
  
  int num_errors = 0;
  for (int k=0; k<cctk_lsh[2]; ++k) {
    for (int j=0; j<cctk_lsh[1]; ++j) {
      for (int i=0; i<cctk_lsh[0]; ++i) {
        int const ind3d = CCTK_GFINDEX3D(cctkGH, i,j,k);
        num_errors += pointtypes[ind3d] != 0 && pointtypes[ind3d] != 2;
      }
    }
  }
  if (num_errors > 0) {
    CCTK_WARN(CCTK_WARN_ABORT, "TestLoopControlPointwise_IntBnd failed");
  }
}



void TestLoopControlPointwise(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_dim != 3) {
    CCTK_WARN(CCTK_WARN_ABORT, "cctk_dim out of range");
  }
  
  TestLoopControlPointwise_All(CCTK_PASS_CTOC);
  TestLoopControlPointwise_Int(CCTK_PASS_CTOC);
  TestLoopControlPointwise_Bnd(CCTK_PASS_CTOC);
  TestLoopControlPointwise_IntBnd(CCTK_PASS_CTOC);
}
