#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
      
      subroutine IDScalarWaveMoL_InitialData (CCTK_ARGUMENTS)
      implicit none
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_FUNCTIONS
      DECLARE_CCTK_PARAMETERS
      CCTK_REAL one, two
      parameter (one = 1)
      parameter (two = 2)
      CCTK_REAL pi
      parameter (pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068d0)
      CCTK_REAL omega
      CCTK_REAL tmp
      integer i, j, k
      if (CCTK_EQUALS(initial_data, "plane wave")) then
         omega = sqrt(wave_number(1)**2 + wave_number(2)**2
     $        + wave_number(3)**2)
         do k=1,cctk_lsh(3)
            do j=1,cctk_lsh(2)
               do i=1,cctk_lsh(1)
                  tmp =  wave_number(1)*(x(i,j,k)-phase_offset(1))
     $                 + wave_number(2)*(y(i,j,k)-phase_offset(2))
     $                 + wave_number(3)*(z(i,j,k)-phase_offset(3))
     $                 + omega*(cctk_time-time_offset)
                  phi(i,j,k) = amplitude * cos (2*pi * tmp)
                  psi(i,j,k) = - amplitude * 2*pi * omega * sin (2*pi * tmp)
               end do
            end do
         end do
      else if (CCTK_EQUALS(initial_data, "Gaussian pulse")) then
         omega = sqrt(pulse_direction(1)**2 + pulse_direction(2)**2
     $        + pulse_direction(3)**2)
         do k=1,cctk_lsh(3)
            do j=1,cctk_lsh(2)
               do i=1,cctk_lsh(1)
                  tmp =  pulse_direction(1)*(x(i,j,k)-pulse_offset(1))
     $                 + pulse_direction(2)*(y(i,j,k)-pulse_offset(2))
     $                 + pulse_direction(3)*(z(i,j,k)-pulse_offset(3))
     $                 + omega*(cctk_time-time_offset)
                  phi(i,j,k) = amplitude * exp (-0.5d0 * tmp**2)
                  psi(i,j,k) = - tmp * omega * phi(i,j,k)
               end do
            end do
         end do
      else if (CCTK_EQUALS(initial_data, "Gaussian")) then
         do k=1,cctk_lsh(3)
            do j=1,cctk_lsh(2)
               do i=1,cctk_lsh(1)
                  tmp = sqrt(  (x(i,j,k)-origin(1))**2
     $                       + (y(i,j,k)-origin(2))**2
     $                       + (z(i,j,k)-origin(3))**2)
                  if (tmp .gt. 1.0d-6 * sigma) then
                     phi(i,j,k) = amplitude/2
     $                    * (  (cctk_time + tmp)/tmp * exp (-0.5d0 * ((cctk_time + tmp - radius) / sigma)**2)
     $                       - (cctk_time - tmp)/tmp * exp (-0.5d0 * ((cctk_time - tmp - radius) / sigma)**2))
                     psi(i,j,k) = amplitude/2
     $                    * (+ (1 + (cctk_time + tmp) * (cctk_time + tmp - radius) / sigma**2) / tmp * exp (-0.5d0 * ((cctk_time + tmp - radius) / sigma)**2)
     $                       - (1 + (cctk_time - tmp) * (cctk_time - tmp - radius) / sigma**2) / tmp * exp (-0.5d0 * ((cctk_time - tmp - radius) / sigma)**2))
                  else
                     phi(i,j,k) = amplitude
     $                    * (  (radius * cctk_time - cctk_time**2 + sigma**2) / sigma**2
     $                       + (radius**3 * cctk_time - cctk_time**4 + 6 * cctk_time**2 * sigma**2 - 3 * sigma**4 - 3 * radius**2 * (cctk_time - sigma)**2 + 3 * radius * (cctk_time**3 - 3 * cctk_time * sigma**2)) / (6 * sigma**6) * tmp**2)
     $                    * exp (-0.5d0 * (radius - cctk_time)**2 / sigma**2)
                     psi(i,j,k) = amplitude
     $                    * (  (radius**2 * cctk_time - 2 * radius * cctk_time**2 + cctk_time**3 + 2 * radius * sigma**2 - 3 * cctk_time * sigma**2) / sigma**4
     $                       + (radius**4 * cctk_time + cctk_time**5 - 10 * cctk_time**3 * sigma**2 + 15 * cctk_time * sigma**4 - 4 * radius**4 * (cctk_time**2 - sigma**2) + 6 * radius**2 * (cctk_time**3 - 3 * cctk_time * sigma**2) - 4 * radius * (cctk_time**4 - 6 * cctk_time**2 * sigma**2 + 3 * sigma**4)) / (6 * sigma**8) * tmp**2)
     $                    * exp (-0.5d0 * (radius - cctk_time)**2 / sigma**2)
                  end if
               end do
            end do
         end do
      else if (CCTK_EQUALS(initial_data, "level index")) then
         do k=1,cctk_lsh(3)
            do j=1,cctk_lsh(2)
               do i=1,cctk_lsh(1)
                  phi(i,j,k) = log (one * cctk_levfac(1)) / log (two)
                  psi(i,j,k) = log (one * cctk_levfac(1)) / log (two)
               end do
            end do
         end do
      end if
      end
