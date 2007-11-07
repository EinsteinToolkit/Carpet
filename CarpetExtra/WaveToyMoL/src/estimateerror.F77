#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

      subroutine WaveToyMoL_EstimateError (CCTK_ARGUMENTS)
      implicit none
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_FUNCTIONS
      DECLARE_CCTK_PARAMETERS
      CCTK_INT bndwidth
      parameter (bndwidth = 1)
      integer imin(3), imax(3)
      CCTK_REAL dx(3)
      integer i, j, k
      integer d
      do d=1,3
         imin(d) = 1+cctk_nghostzones(d)
         imax(d) = cctk_lsh(d)-cctk_nghostzones(d)
         if (cctk_bbox(2*d-1).ne.0) imin(d) = 1+bndwidth
         if (cctk_bbox(2*d  ).ne.0) imax(d) = cctk_lsh(d)-bndwidth
      end do
      do d=1,3
         dx(d) = CCTK_DELTA_SPACE(d)
      end do
      do k=imin(3),imax(3)
         do j=imin(2),imax(2)
            do i=imin(1),imax(1)
               error(i,j,k) = max (
     $              abs(phi(i+1,j,k) - phi(i-1,j,k)),
     $              abs(phi(i,j+1,k) - phi(i,j-1,k)),
     $              abs(phi(i,j,k+1) - phi(i,j,k-1)))
            end do
         end do
      end do
      end

      subroutine WaveToyMoL_EstimateErrorBoundaries (CCTK_ARGUMENTS)
      implicit none
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_FUNCTIONS
      DECLARE_CCTK_PARAMETERS
      CCTK_INT bndwidth
      parameter (bndwidth = 1)
      integer ierr
      ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, bndwidth, -1,
     $     "WaveToyMoL::scalarevolveerrorestimate", "scalar");
      if (ierr.ne.0) call CCTK_WARN (0, "internal error")
      end
