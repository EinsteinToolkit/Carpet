#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "util_Table.h"
#include <assert.h>
#include <math.h>

void CoM_Local(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int i, j, k, index;

  int nx = cctk_lsh[0];
  int ny = cctk_lsh[1];
  int nz = cctk_lsh[2];

  for (k = 0; k < nz; k++)
    for (j = 0; j < ny; j++)
      for (i = 0; i < nx; i++) {

        index = CCTK_GFINDEX3D(cctkGH, i, j, k);

        dMx[index] = dens[index] * x[index];
        dMy[index] = dens[index] * y[index];
        dMz[index] = dens[index] * z[index];
      }
}

void CoM(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int varindex = -1;
  int ierr = 0;
  int reduction_handle;

  CCTK_REAL denstotal;

  CCTK_REAL sym_factor1, sym_factor2, sym_factor3;

  if (CCTK_EQUALS(domain, "bitant")) {
    sym_factor1 = 2.0e0;
    sym_factor2 = 2.0e0;
    sym_factor3 = 0.0e0;
  } else if (CCTK_EQUALS(domain, "octant")) {
    sym_factor1 = 8.0e0;
    sym_factor2 = 0.0e0;
    sym_factor3 = 0.0e0;
  } else {
    sym_factor1 = 1.0e0;
    sym_factor2 = 1.0e0;
    sym_factor3 = 1.0e0;
  }

  reduction_handle = CCTK_ReductionHandle("sum");

  varindex = CCTK_VarIndex("GRHydro::dens");
  assert(varindex >= 0);
  ierr = CCTK_Reduce(cctkGH, -1, reduction_handle, 1, CCTK_VARIABLE_REAL,
                     (void *)&denstotal, 1, varindex);
  assert(!ierr);

  varindex = CCTK_VarIndex("ReductionTest::dMx");
  assert(varindex >= 0);
  ierr = CCTK_Reduce(cctkGH, -1, reduction_handle, 1, CCTK_VARIABLE_REAL,
                     (void *)Mx, 1, varindex);
  assert(!ierr);

  varindex = CCTK_VarIndex("ReductionTest::dMy");
  assert(varindex >= 0);
  ierr = CCTK_Reduce(cctkGH, -1, reduction_handle, 1, CCTK_VARIABLE_REAL,
                     (void *)My, 1, varindex);
  assert(!ierr);

  varindex = CCTK_VarIndex("ReductionTest::dMz");
  assert(varindex >= 0);
  ierr = CCTK_Reduce(cctkGH, -1, reduction_handle, 1, CCTK_VARIABLE_REAL,
                     (void *)Mz, 1, varindex);
  assert(!ierr);

  denstotal = sym_factor1 * denstotal;
  *Mx = sym_factor2 * (*Mx) / (denstotal);
  *My = sym_factor2 * (*My) / (denstotal);
  *Mz = sym_factor3 * (*Mz) / (denstotal);
  *Mr = sqrt((*Mx) * (*Mx) + (*My) * (*My) + (*Mz) * (*Mz));

  CCTK_VInfo(CCTK_THORNSTRING, "Mr: %15.6E Mx: %15.6E My: %15.6E Mz %15.6E",
             *Mr, *Mx, *My, *Mz);
}
