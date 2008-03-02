#include <assert.h>
#include <stdio.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include "carpet.hh"
#include "defs.hh"
#include "vect.hh"

#include "util_Table.h"



extern "C" void CarpetIntegrate_Global(CCTK_ARGUMENTS);

void CarpetIntegrate_Global(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  using namespace std;
  using namespace Carpet;

  int varindex = -1;
  int ierr = 0;

  CCTK_REAL integral;

  int reduction_handle = CCTK_ReductionHandle("sum");

  varindex = CCTK_VarIndex("CarpetIntegrateTest::integrand");
  assert(varindex>=0);
  ierr = CCTK_Reduce(cctkGH, -1, reduction_handle, 
		     1, CCTK_VARIABLE_REAL, &integral, 1, varindex);
  assert(!ierr);

  CCTK_REAL d3x = cctk_delta_space[0]*cctk_delta_space[1]*cctk_delta_space[2];
  integral *= d3x;

  printf("Integral: %g\n", (double)integral);
}
