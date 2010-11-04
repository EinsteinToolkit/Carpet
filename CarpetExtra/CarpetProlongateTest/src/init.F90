#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

subroutine CarpetProlongateTest_Init (CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  u = x**power_x + y**power_y + z**power_z
end subroutine CarpetProlongateTest_Init

subroutine CarpetProlongateTest_Diff (CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  du = u - (x**power_x + y**power_y + z**power_z)
end subroutine CarpetProlongateTest_Diff
