c     $Header: /home/eschnett/C/carpet/Carpet/CarpetExtra/IDScalarWaveExpl/src/initialdata.F77,v 1.2 2003/11/05 16:18:40 schnetter Exp $

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
      
      subroutine IDScalarWaveExpl_InitialData (CCTK_ARGUMENTS)
      implicit none
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_FUNCTIONS
      DECLARE_CCTK_PARAMETERS
      CCTK_REAL pi
      parameter (pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068d0)
      CCTK_REAL omega
      integer i, j, k
      omega = sqrt(wave_number(1)**2 + wave_number(2)**2 + wave_number(3)**2)
      do k=1,cctk_lsh(3)
         do j=1,cctk_lsh(2)
            do i=1,cctk_lsh(1)
               phi(i,j,k) = amplitude * cos (2*pi *
     $              ( wave_number(1)*(x(i,j,k)-phase_offset(1))
     $              + wave_number(2)*(y(i,j,k)-phase_offset(2))
     $              + wave_number(3)*(z(i,j,k)-phase_offset(3))
     $              + omega*(cctk_time-time_offset)))
               psi(i,j,k) = - amplitude * 2*pi * omega * sin (2*pi *
     $              ( wave_number(1)*(x(i,j,k)-phase_offset(1))
     $              + wave_number(2)*(y(i,j,k)-phase_offset(2))
     $              + wave_number(3)*(z(i,j,k)-phase_offset(3))
     $              + omega*(cctk_time-time_offset)))
            end do
         end do
      end do
      end
