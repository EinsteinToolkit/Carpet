#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

      subroutine WaveToyMoL_RegisterVars (CCTK_ARGUMENTS)
      implicit none
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_FUNCTIONS
      DECLARE_CCTK_PARAMETERS
      integer group, groupdot
      integer ierr
      call CCTK_GroupIndex(group, "wavetoymol::scalarevolve")
      if (group.lt.0) call CCTK_WARN (0, "internal error")
      call CCTK_GroupIndex(groupdot, "wavetoymol::scalarevolvedot")
      if (groupdot.lt.0) call CCTK_WARN (0, "internal error")
      ierr = MoLRegisterEvolvedGroup (group, groupdot)
      if (ierr.ne.0) call CCTK_WARN (0, "internal error")
      end
