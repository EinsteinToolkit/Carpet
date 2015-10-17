c     $Header:$

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

      subroutine WaveToyExpl_IntInit (CCTK_ARGUMENTS)
      implicit none
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_FUNCTIONS
      DECLARE_CCTK_PARAMETERS
      integer i, j, k
      do k=1,cctk_lsh(3)
         do j=1,cctk_lsh(2)
            do i=1,cctk_lsh(1)
               phiint(i,j,k) = phi(i,j,k)
               phiint_p(i,j,k) = phi_p(i,j,k)
               phiint_p_p(i,j,k) = phi_p_p(i,j,k)
               psiint(i,j,k) = psi(i,j,k)
               psiint_p(i,j,k) = psi_p(i,j,k)
               psiint_p_p(i,j,k) = psi_p_p(i,j,k)
            end do
         end do
      end do
      end
