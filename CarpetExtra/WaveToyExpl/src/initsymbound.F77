c     $Header: /home/eschnett/C/carpet/Carpet/CarpetExtra/WaveToyExpl/src/initsymbound.F77,v 1.1 2003/06/18 18:24:30 schnetter Exp $

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

      subroutine WaveToyExpl_InitSymBound (CCTK_ARGUMENTS)
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
      call SetCartSymGN (ierr, cctkGH, sym, "WaveToyExpl::scalarevolve")
      if (ierr.ne.0) call CCTK_WARN (0, "internal error")
      call SetCartSymGN (ierr, cctkGH, sym, "WaveToyExpl::scalarevolveint")
      if (ierr.ne.0) call CCTK_WARN (0, "internal error")
      end
