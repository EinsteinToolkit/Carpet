c     $Header: /home/eschnett/C/carpet/Carpet/CarpetExtra/WaveToyExpl/src/startup.F77,v 1.1 2003/06/18 18:24:30 schnetter Exp $

#include "cctk.h"
#include "cctk_Parameters.h"

      subroutine WaveToyExpl_Startup
      implicit none
      DECLARE_CCTK_FUNCTIONS
      DECLARE_CCTK_PARAMETERS
      integer ierr
      call CCTK_RegisterBanner (ierr, "WaveToyExpl: Evolving a scalar field with explicit intermediate integration steps")
      if (ierr.ne.0) call CCTK_WARN (0, "internal error")
      end
