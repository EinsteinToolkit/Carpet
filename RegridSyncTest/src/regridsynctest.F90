#include <cctk.h>
#include <cctk_Parameters.h>
#include <cctk_Arguments.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Functions.h>



subroutine RegridSyncTest_setup(CCTK_ARGUMENTS)

  implicit none
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS_RegridSyncTest_setup

  myregridtestgf = 0.9d0

end subroutine RegridSyncTest_setup

subroutine RegridSyncTest_sync(CCTK_ARGUMENTS)

  implicit none
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS_RegridSyncTest_sync

end subroutine RegridSyncTest_sync



subroutine RegridSyncTest_do_something(CCTK_ARGUMENTS)

  implicit none
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS_RegridSyncTest_do_something

  myregridtestgf = 1.0d0

end subroutine RegridSyncTest_do_something



subroutine RegridSyncTest_evolve(CCTK_ARGUMENTS)

  implicit none
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS_RegridSyncTest_evolve

  myregridtestgf = 0.8d0

end subroutine RegridSyncTest_evolve
