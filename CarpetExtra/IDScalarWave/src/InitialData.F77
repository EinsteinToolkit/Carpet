c     -*-Fortran-*-

 /*@@
   @file      InitialData.F77
   @date      
   @author    Tom Goodale
   @desc 
              Initial data for the 3D Wave Equation
   @enddesc 
 @@*/

#include "cctk.h" 
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"


 /*@@
   @routine    IDScalarWave_InitialData
   @date       
   @author     Tom Goodale
   @desc 
               Set up initial data for the wave equation
   @enddesc 
   @calls      
   @calledby   
   @history 
 
   @endhistory 

@@*/

      subroutine IDScalarWave_InitialData (CCTK_ARGUMENTS)

      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS

      INTEGER    i,j,k
      CCTK_REAL  dt,omega, pi
      
c     call CCTK_INFO ("IDScalarWave_InitialData")
      
      pi = 4.d0*atan(1.d0)

      dt = CCTK_DELTA_TIME

      omega = sqrt(kx**2+ky**2+kz**2)

      if (CCTK_EQUALS(initial_data,"plane")) then 
         
         do k=1,cctk_lsh(3)
            do j=1,cctk_lsh(2)
               do i=1,cctk_lsh(1)

                  phi(i,j,k) = amplitude
     $                 * cos((kx*cart3d_x(i,j,k) + ky*cart3d_y(i,j,k)
     $                 + kz*cart3d_z(i,j,k) + omega*cctk_time) * pi)
                  phi_p(i,j,k) = amplitude
     $                 * cos((kx*cart3d_x(i,j,k) + ky*cart3d_y(i,j,k)
     $                 + kz*cart3d_z(i,j,k) + omega*(cctk_time - dt)) * pi)
                  
               end do
            end do
         end do

      else if (CCTK_EQUALS(initial_data,"gaussian")) then 

         do k=1, cctk_lsh(3)
            do j=1, cctk_lsh(2)
               do i=1, cctk_lsh(1)

                  phi(i,j,k) = amplitude
     $                 * exp(- (spher3d_r(i,j,k) - radius)**2 / sigma**2)
                  phi_p(i,j,k) = amplitude
     $                 * exp(- (spher3d_r(i,j,k)**2 - radius - dt)**2 / sigma**2)

               end do
            end do
         end do

      else if (CCTK_EQUALS(initial_data, "box")) then

c        Use kx,ky,kz as number of modes in each direction.

         do k=1,cctk_lsh(3)
            do j=1,cctk_lsh(2)
               do i=1,cctk_lsh(1)

                  phi(i,j,k) = amplitude
     $                 * sin(kx * (cart3d_x(i,j,k) - 0.5d0) * pi)
     $                 * sin(ky * (cart3d_y(i,j,k) - 0.5d0) * pi)
     $                 * sin(kz * (cart3d_z(i,j,k) - 0.5d0) * pi)
     $                 * cos(omega * cctk_time * pi)
                  
                  phi_p(i,j,k) = amplitude
     $                 * sin(kx * (cart3d_x(i,j,k) - 0.5d0) * pi)
     $                 * sin(ky * (cart3d_y(i,j,k) - 0.5d0) * pi)
     $                 * sin(kz * (cart3d_z(i,j,k) - 0.5d0) * pi)
     $                 * cos(omega * (cctk_time-dt) * pi)

               end do
            end do
         end do

      else

         do k=1,cctk_lsh(3)
            do j=1,cctk_lsh(2)
               do i=1,cctk_lsh(1)
                  
                  phi(i,j,k)   = 0.0d0
                  phi_p(i,j,k) = 0.0d0

               end do
            end do
         end do 
         
      end if
      
      end 
