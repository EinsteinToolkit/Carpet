c     $Header:$

#include "cctk.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

      subroutine WaveToyFO_Startup
      implicit none
      DECLARE_CCTK_FUNCTIONS
      DECLARE_CCTK_PARAMETERS
      integer ierr
      call CCTK_RegisterBanner (ierr, "WaveToyFO: Evolving a scalar field as a first order system")
      if (ierr.ne.0) call CCTK_WARN (0, "internal error")
      end
