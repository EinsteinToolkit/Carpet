c     $Header: /home/eschnett/C/carpet/Carpet/CarpetExtra/WaveToyMoL/src/initsymbound.F77,v 1.3 2004/04/07 17:08:21 schnetter Exp $

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

      subroutine WaveToyMoL_InitSymBound (CCTK_ARGUMENTS)
      implicit none
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_FUNCTIONS
      DECLARE_CCTK_PARAMETERS
      integer sym(3)
      integer d
      integer ierr
      do d=1,3
         sym(d) = +1
      end do
      call SetCartSymVN (ierr, cctkGH, sym, "WaveToyMoL::phi")
      if (ierr.ne.0) call CCTK_WARN (0, "internal error")
      call SetCartSymVN (ierr, cctkGH, sym, "WaveToyMoL::phidot")
      if (ierr.ne.0) call CCTK_WARN (0, "internal error")
      call SetCartSymVN (ierr, cctkGH, sym, "WaveToyMoL::error")
      if (ierr.ne.0) call CCTK_WARN (0, "internal error")
      end
