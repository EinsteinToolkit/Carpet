c     $Header:$

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
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
