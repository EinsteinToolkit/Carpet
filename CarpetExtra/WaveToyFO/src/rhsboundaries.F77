c     $Header:$

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

      subroutine WaveToyFO_RHSBoundaries (CCTK_ARGUMENTS)
      implicit none
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_FUNCTIONS
      DECLARE_CCTK_PARAMETERS
      CCTK_INT bndwidth
      parameter (bndwidth = 1)
      integer ierr
      ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, bndwidth, -1,
     $     "WaveToyFO::scalarevolvedot", "none");
      if (ierr.ne.0) call CCTK_WARN (0, "internal error")
      end
