#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

      subroutine WaveToyMoL_CalcRHS (CCTK_ARGUMENTS)
      implicit none
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_FUNCTIONS
      DECLARE_CCTK_PARAMETERS
      CCTK_INT bndwidth
      parameter (bndwidth = 1)
      integer imin(3), imax(3)
      CCTK_REAL dx(3)
      integer i, j, k
      logical isnearboundary
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
               phidot(i,j,k) = psi(i,j,k)
               isnearboundary = .false.
               isnearboundary = isnearboundary .or. i.eq.imin(1) .or. i.eq.imax(1)
               isnearboundary = isnearboundary .or. j.eq.imin(2) .or. j.eq.imax(2)
               isnearboundary = isnearboundary .or. k.eq.imin(3) .or. k.eq.imax(3)
               if (order.eq.2 .or. isnearboundary) then
                  psidot(i,j,k) = (
     $                 + (phi(i-1,j,k) - 2*phi(i,j,k) + phi(i+1,j,k)) / dx(1)**2
     $                 + (phi(i,j-1,k) - 2*phi(i,j,k) + phi(i,j+1,k)) / dx(2)**2
     $                 + (phi(i,j,k-1) - 2*phi(i,j,k) + phi(i,j,k+1)) / dx(3)**2)
               else if (order.eq.4) then
                  psidot(i,j,k) = (
     $                 + (- phi(i-2,j,k) + 16*phi(i-1,j,k) - 30*phi(i,j,k) + 16*phi(i+1,j,k) - phi(i+2,j,k)) / (12 * dx(1)**2)
     $                 + (- phi(i,j-2,k) + 16*phi(i,j-1,k) - 30*phi(i,j,k) + 16*phi(i,j+1,k) - phi(i,j+2,k)) / (12 * dx(2)**2)
     $                 + (- phi(i,j,k-2) + 16*phi(i,j,k-1) - 30*phi(i,j,k) + 16*phi(i,j,k+1) - phi(i,j,k+2)) / (12 * dx(3)**2))
               else
                  call CCTK_WARN (0, "internal error")
               end if
            end do
         end do
      end do
      end
