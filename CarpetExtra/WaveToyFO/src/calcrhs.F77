c     $Header: /home/eschnett/C/carpet/Carpet/CarpetExtra/WaveToyFO/src/calcrhs.F77,v 1.4 2004/05/06 14:21:18 schnetter Exp $

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
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
      integer im, jm, km
      integer ip, jp, kp
      integer d
      do d=1,3
         bndwidth(d) = 0
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
               im = max(i-1,1)
               jm = max(j-1,1)
               km = max(k-1,1)
               ip = min(i+1,cctk_lsh(1))
               jp = min(j+1,cctk_lsh(2))
               kp = min(k+1,cctk_lsh(3))
               phidot(i,j,k) =
     $              ( (psix(ip,j,k) - psix(im,j,k)) / ((ip-im)*dx(1))
     $              + (psiy(i,jp,k) - psiy(i,jm,k)) / ((jp-jm)*dx(2))
     $              + (psiz(i,j,kp) - psiz(i,j,km)) / ((kp-km)*dx(3)))
               psixdot(i,j,k) = (phi(ip,j,k) - phi(im,j,k)) / ((ip-im)*dx(1))
               psiydot(i,j,k) = (phi(i,jp,k) - phi(i,jm,k)) / ((jp-jm)*dx(2))
               psizdot(i,j,k) = (phi(i,j,kp) - phi(i,j,km)) / ((kp-km)*dx(3))
            end do
         end do
      end do
      end
