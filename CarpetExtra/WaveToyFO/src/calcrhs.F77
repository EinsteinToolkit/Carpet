c     $Header: /home/eschnett/C/carpet/Carpet/CarpetExtra/WaveToyFO/src/calcrhs.F77,v 1.2 2003/06/18 20:08:29 schnetter Exp $

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

      subroutine WaveToyFO_CalcRHS (CCTK_ARGUMENTS)
      implicit none
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_FUNCTIONS
      DECLARE_CCTK_PARAMETERS
      CCTK_REAL dx(3)
      integer bndwidth(3)
      integer imin(3), imax(3)
      integer i, j, k
      integer d
      integer ierr
      do d=1,3
         bndwidth(d) = 1
      end do
      do d=1,3
         if (cctk_bbox(2*d).ne.0) then
            imin(d) = 1+bndwidth(d)
         else
            imin(d) = 1+cctk_nghostzones(d)
         end if
         if (cctk_bbox(2*d+1).ne.0) then
            imax(d) = cctk_lsh(d)-bndwidth(d)
         else
            imax(d) = cctk_lsh(d)-cctk_nghostzones(d)
         end if
      end do
      do d=1,3
         dx(d) = CCTK_DELTA_SPACE(3)
      end do
      do k=imin(3),imax(3)
         do j=imin(2),imax(2)
            do i=imin(1),imax(1)
               phidot(i,j,k) =
     $              ( (psix(i+1,j,k) - psix(i-1,j,k)) / (2*dx(1))
     $              + (psiy(i,j+1,k) - psiy(i,j-1,k)) / (2*dx(2))
     $              + (psiz(i,j,k+1) - psiz(i,j,k-1)) / (2*dx(3)))
               psixdot(i,j,k) = (phi(i+1,j,k) - phi(i-1,j,k)) / (2*dx(1))
               psiydot(i,j,k) = (phi(i,j+1,k) - phi(i,j-1,k)) / (2*dx(2))
               psizdot(i,j,k) = (phi(i,j,k+1) - phi(i,j,k-1)) / (2*dx(3))
            end do
         end do
      end do
      end
