c     $Header:$

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

      subroutine WaveToyExpl_FinalStep (CCTK_ARGUMENTS)
      implicit none
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_FUNCTIONS
      DECLARE_CCTK_PARAMETERS
      CCTK_REAL dx(3), dt
      integer i, j, k
      integer d
      do d=1,3
         dx(d) = CCTK_DELTA_SPACE(d)
      end do
      dt = CCTK_DELTA_TIME
      do k=1+cctk_nghostzones(3),cctk_lsh(3)-cctk_nghostzones(3)
         do j=1+cctk_nghostzones(2),cctk_lsh(2)-cctk_nghostzones(2)
            do i=1+cctk_nghostzones(1),cctk_lsh(1)-cctk_nghostzones(1)
c     d/dt phi = psi
               phi(i,j,k) = phi_p(i,j,k) + dt/2 * (psi_p(i,j,k) + psiint(i,j,k))
c     d/dt psi = d^2/dx^2 phi
               psi(i,j,k) = psi_p(i,j,k) + dt/2 *
     $              ( (phi_p(i-1,j,k) - 2*phi_p(i,j,k) + phi_p(i+1,j,k)) / dx(1)**2
     $              + (phi_p(i,j-1,k) - 2*phi_p(i,j,k) + phi_p(i,j+1,k)) / dx(2)**2
     $              + (phi_p(i,j,k-1) - 2*phi_p(i,j,k) + phi_p(i,j,k+1)) / dx(3)**2
     $              + (phiint(i-1,j,k) - 2*phiint(i,j,k) + phiint(i+1,j,k)) / dx(1)**2
     $              + (phiint(i,j-1,k) - 2*phiint(i,j,k) + phiint(i,j+1,k)) / dx(2)**2
     $              + (phiint(i,j,k-1) - 2*phiint(i,j,k) + phiint(i,j,k+1)) / dx(3)**2)
            end do
         end do
      end do
      end
