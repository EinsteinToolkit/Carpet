#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "AlphaThorns/Cart3d/src/Cart3d.h"

int Cart3dTest_Evolve (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  return Cart3dSymGN (cctkGH, "Cart3dTest::quantities");
}
