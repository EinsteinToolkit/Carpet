c     $Header: /home/eschnett/C/carpet/Carpet/CarpetExtra/WaveToyExpl/src/startup.F77,v 1.2 2003/10/27 15:31:40 schnetter Exp $

#include "cctk.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

      subroutine WaveToyExpl_Startup
      implicit none
      DECLARE_CCTK_FUNCTIONS
      DECLARE_CCTK_PARAMETERS
      integer ierr
      call CCTK_RegisterBanner (ierr, "WaveToyExpl: Evolving a scalar field with explicit intermediate integration steps")
      if (ierr.ne.0) call CCTK_WARN (0, "internal error")
      end
