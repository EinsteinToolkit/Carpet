c     $Header: /home/eschnett/C/carpet/Carpet/CarpetExtra/WaveToyMoL/src/boundaries.F77,v 1.1 2003/06/18 18:24:31 schnetter Exp $

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Faces.h"
#include "cctk_Parameters.h"

      subroutine WaveToyMoL_Boundaries (CCTK_ARGUMENTS)
      implicit none
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_FUNCTIONS
      DECLARE_CCTK_PARAMETERS
      external Boundary_SelectGroupForBC
      CCTK_INT Boundary_SelectGroupForBC
      CCTK_INT bndwidth
      parameter (bndwidth = 1)
      character fbound*100
      CCTK_INT  fboundlen
      integer ierr
      call CCTK_FortranString (fboundlen, bound, fbound)
      if (fboundlen.lt.0) call CCTK_WARN (0, "internal error")
      ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, bndwidth, -1,
     $     "WaveToyMoL::scalarevolve", fbound);
      if (ierr.ne.0) call CCTK_WARN (0, "internal error")
      call CartSymGN (ierr, cctkGH, "WaveToyMoL::scalarevolve")
      if (ierr.ne.0) call CCTK_WARN (0, "internal error")
      end
