c     $Header: /home/eschnett/C/carpet/Carpet/CarpetExtra/IDScalarWaveFO/src/initialdata.F77,v 1.5 2004/05/07 22:51:08 schnetter Exp $

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
      
      subroutine IDScalarWaveFO_InitialData (CCTK_ARGUMENTS)
      implicit none
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_FUNCTIONS
      DECLARE_CCTK_PARAMETERS
      CCTK_REAL pi
      parameter (pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068d0)
      CCTK_REAL omega
      integer i, j, k
      
      if (CCTK_EQUALS(initial_data, "plane")) then
         omega = sqrt(wave_number(1)**2 + wave_number(2)**2 + wave_number(3)**2)
         do k=1,cctk_lsh(3)
            do j=1,cctk_lsh(2)
               do i=1,cctk_lsh(1)
                  phi(i,j,k) = amplitude * cos (2*pi *
     $                 ( wave_number(1)*(x(i,j,k)-phase_offset(1))
     $                 + wave_number(2)*(y(i,j,k)-phase_offset(2))
     $                 + wave_number(3)*(z(i,j,k)-phase_offset(3))
     $                 + omega*(cctk_time-time_offset)))
                  psix(i,j,k) = amplitude * wave_number(1) / omega * cos (2*pi *
     $                 ( wave_number(1)*(x(i,j,k)-phase_offset(1))
     $                 + wave_number(2)*(y(i,j,k)-phase_offset(2))
     $                 + wave_number(3)*(z(i,j,k)-phase_offset(3))
     $                 + omega*(cctk_time-time_offset)))
                  psiy(i,j,k) = amplitude * wave_number(2) / omega * cos (2*pi *
     $                 ( wave_number(1)*(x(i,j,k)-phase_offset(1))
     $                 + wave_number(2)*(y(i,j,k)-phase_offset(2))
     $                 + wave_number(3)*(z(i,j,k)-phase_offset(3))
     $                 + omega*(cctk_time-time_offset)))
                  psiz(i,j,k) = amplitude * wave_number(3) / omega * cos (2*pi *
     $                 ( wave_number(1)*(x(i,j,k)-phase_offset(1))
     $                 + wave_number(2)*(y(i,j,k)-phase_offset(2))
     $                 + wave_number(3)*(z(i,j,k)-phase_offset(3))
     $                 + omega*(cctk_time-time_offset)))
               end do
            end do
         end do
         
      else if (CCTK_EQUALS(initial_data, "noise")) then
         
         do k=1,cctk_lsh(3)
            do j=1,cctk_lsh(2)
               do i=1,cctk_lsh(1)
                  call random_number (phi(i,j,k))
                  call random_number (psix(i,j,k))
                  call random_number (psiy(i,j,k))
                  call random_number (psiz(i,j,k))
               end do
            end do
         end do
         
      end if
      end
