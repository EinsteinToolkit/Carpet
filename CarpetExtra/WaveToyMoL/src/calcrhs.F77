c     $Header:$

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
      integer d
      do d=1,3
         imin(d) = 1+bndwidth
         imax(d) = cctk_lsh(d)-bndwidth
      end do
      do d=1,3
         dx(d) = CCTK_DELTA_SPACE(3)
      end do
      do k=imin(3),imax(3)
         do j=imin(2),imax(2)
            do i=imin(1),imax(1)
               phidot(i,j,k) = psi(i,j,k)
               psidot(i,j,k) =
     $              ( (phi(i-1,j,k) - 2*phi(i,j,k) + phi(i+1,j,k)) / dx(1)**2
     $              + (phi(i,j-1,k) - 2*phi(i,j,k) + phi(i,j+1,k)) / dx(2)**2
     $              + (phi(i,j,k-1) - 2*phi(i,j,k) + phi(i,j,k+1)) / dx(3)**2)
            end do
         end do
      end do
      end
