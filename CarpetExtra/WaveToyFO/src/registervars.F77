c     $Header: /home/eschnett/C/carpet/Carpet/CarpetExtra/WaveToyFO/src/registervars.F77,v 1.1 2003/06/18 18:24:30 schnetter Exp $

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

      subroutine WaveToyFO_RegisterVars (CCTK_ARGUMENTS)
      implicit none
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_FUNCTIONS
      DECLARE_CCTK_PARAMETERS
      external MolRegisterEvolvedGroup
      CCTK_INT MolRegisterEvolvedGroup
      integer group, groupdot
      integer ierr
      call CCTK_GroupIndex(group, "WaveToyFO::scalarevolve")
      if (group.lt.0) call CCTK_WARN (0, "internal error")
      call CCTK_GroupIndex(groupdot, "WaveToyFO::scalarevolvedot")
      if (groupdot.lt.0) call CCTK_WARN (0, "internal error")
      ierr = MolRegisterEvolvedGroup (group, groupdot)
      if (ierr.ne.0) call CCTK_WARN (0, "internal error")
      end
