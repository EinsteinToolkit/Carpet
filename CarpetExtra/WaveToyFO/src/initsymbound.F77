c     $Header: /home/eschnett/C/carpet/Carpet/CarpetExtra/WaveToyFO/src/initsymbound.F77,v 1.2 2003/10/27 15:31:41 schnetter Exp $

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

      subroutine WaveToyFO_InitSymBound (CCTK_ARGUMENTS)
      implicit none
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_FUNCTIONS
      DECLARE_CCTK_PARAMETERS
      integer sym(3), symx(3), symy(3), symz(3)
      data sym  /+1, +1, +1/
      data symx /-1, +1, +1/
      data symy /+1, -1, +1/
      data symz /+1, +1, -1/
      integer ierr
      call SetCartSymVN (ierr, cctkGH, sym, "WaveToyFO::phi")
      if (ierr.ne.0) call CCTK_WARN (0, "internal error")
      call SetCartSymVN (ierr, cctkGH, symx, "WaveToyFO::psix")
      if (ierr.ne.0) call CCTK_WARN (0, "internal error")
      call SetCartSymVN (ierr, cctkGH, symy, "WaveToyFO::psiy")
      if (ierr.ne.0) call CCTK_WARN (0, "internal error")
      call SetCartSymVN (ierr, cctkGH, symz, "WaveToyFO::psiz")
      if (ierr.ne.0) call CCTK_WARN (0, "internal error")
      end
