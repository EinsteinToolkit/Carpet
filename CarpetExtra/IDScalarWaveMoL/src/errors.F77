c     $Header: /home/eschnett/C/carpet/Carpet/CarpetExtra/IDScalarWaveMoL/src/errors.F77,v 1.1 2003/06/18 18:24:29 schnetter Exp $

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
      
      subroutine IDScalarWaveMoL_Errors (CCTK_ARGUMENTS)
      implicit none
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_FUNCTIONS
      DECLARE_CCTK_PARAMETERS
      CCTK_REAL tmp
      integer i, j, k
      do k=1,cctk_lsh(3)
         do j=1,cctk_lsh(2)
            do i=1,cctk_lsh(1)
               phierror(i,j,k) = phi(i,j,k)
               psierror(i,j,k) = psi(i,j,k)
            end do
         end do
      end do
      call IDScalarWaveMol_InitialData (CCTK_PASS_FTOF)
      do k=1,cctk_lsh(3)
         do j=1,cctk_lsh(2)
            do i=1,cctk_lsh(1)
               tmp = phi(i,j,k)
               phi(i,j,k) = phierror(i,j,k)
               phierror(i,j,k) = phi(i,j,k) - tmp
               tmp = psi(i,j,k)
               psi(i,j,k) = psierror(i,j,k)
               psierror(i,j,k) = psi(i,j,k) - tmp
            end do
         end do
      end do
      end
