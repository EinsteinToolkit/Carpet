c     -*-Fortran-*-
c     $Header:$

 /*@@
   @file      InitialData.F77
   @date      
   @author    Tom Goodale, Erik Schnetter
   @desc 
              Initial data for the 3D Wave Equation
   @enddesc 
 @@*/

#include "cctk.h" 
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"


 /*@@
   @routine    IDScalarWave_InitialData
   @date       
   @author     Tom Goodale, Erik Schnetter
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
      CCTK_REAL  t
      CCTK_REAL  x,y,z, r
      
c     call CCTK_INFO ("IDScalarWave_InitialData")
      
      pi = 4.d0*atan(1.d0)

      t = cctk_time
      dt = CCTK_DELTA_TIME

      omega = sqrt(kx**2+ky**2+kz**2)

      if (CCTK_EQUALS(initial_data,"plane")) then 
         
         do k=1,cctk_lsh(3)
            do j=1,cctk_lsh(2)
               do i=1,cctk_lsh(1)

                  x = cart3d_x(i,j,k)
                  y = cart3d_y(i,j,k)
                  z = cart3d_z(i,j,k)
                  
                  phi(i,j,k) = amplitude
     $                 * cos((kx*x + ky*y + kz*z + omega*cctk_time) * pi)

                  phi_p(i,j,k) = amplitude
     $                 * cos((kx*x + ky*y + kz*z + omega*(cctk_time - dt)) * pi)

                  phi_p_p(i,j,k) = amplitude
     $                 * cos((kx*x + ky*y + kz*z + omega*(cctk_time - 2*dt)) * pi)
                  
               end do
            end do
         end do

      else if (CCTK_EQUALS(initial_data,"gaussian")) then 

         do k=1, cctk_lsh(3)
            do j=1, cctk_lsh(2)
               do i=1, cctk_lsh(1)
                  
                  r = spher3d_r(i,j,k)
                  
                  if (r .ne. 0) then
                     phi(i,j,k) = amplitude/2 * (r - radius - t) / r
     $                    * exp(- (r - radius - t)**2 / sigma**2)
     $                    + amplitude/2 * (r + radius + t) / r
     $                    * exp(- (r + radius + t)**2 / sigma**2)
                     
                     phi_p(i,j,k) = amplitude/2 * (r - radius - t + dt) / r
     $                    * exp(- (r - radius - t + dt)**2 / sigma**2)
     $                    + amplitude/2 * (r + radius + t - dt) / r
     $                    * exp(- (r + radius + t - dt)**2 / sigma**2)
                     
                     phi_p_p(i,j,k) = amplitude/2 * (r - radius - t + 2*dt) / r
     $                    * exp(- (r - radius - t + 2*dt)**2 / sigma**2)
     $                    + amplitude/2 * (r + radius + t - 2*dt) / r
     $                    * exp(- (r + radius + t - 2*dt)**2 / sigma**2)
                     
                  else
                     
                     phi(i,j,k) = amplitude
     $                    * (1 - 2 * (radius + t)**2 / sigma**2)
     $                    * exp(- (radius + t)**2 / sigma**2)
                     
                     phi_p(i,j,k) = amplitude
     $                    * (1 - 2 * (radius + t - dt)**2 / sigma**2)
     $                    * exp(- (radius + t - dt)**2 / sigma**2)
                     
                     phi_p_p(i,j,k) = amplitude
     $                    * (1 - 2 * (radius + t - 2*dt)**2 / sigma**2)
     $                    * exp(- (radius + t - 2*dt)**2 / sigma**2)
                     
                  end if
                  
               end do
            end do
         end do

      else if (CCTK_EQUALS(initial_data, "box")) then

c        Use kx,ky,kz as number of modes in each direction.

         do k=1,cctk_lsh(3)
            do j=1,cctk_lsh(2)
               do i=1,cctk_lsh(1)

                  x = cart3d_x(i,j,k)
                  y = cart3d_y(i,j,k)
                  z = cart3d_z(i,j,k)

                  phi(i,j,k) = amplitude
     $                 * sin(kx * (x - 0.5d0) * pi)
     $                 * sin(ky * (y - 0.5d0) * pi)
     $                 * sin(kz * (z - 0.5d0) * pi)
     $                 * cos(omega * cctk_time * pi)
                  
                  phi_p(i,j,k) = amplitude
     $                 * sin(kx * (x - 0.5d0) * pi)
     $                 * sin(ky * (y - 0.5d0) * pi)
     $                 * sin(kz * (z - 0.5d0) * pi)
     $                 * cos(omega * (cctk_time - dt) * pi)
                  
                  phi_p_p(i,j,k) = amplitude
     $                 * sin(kx * (x - 0.5d0) * pi)
     $                 * sin(ky * (y - 0.5d0) * pi)
     $                 * sin(kz * (z - 0.5d0) * pi)
     $                 * cos(omega * (cctk_time - 2*dt) * pi)

               end do
            end do
         end do

      else if (CCTK_EQUALS(initial_data, "1/r")) then
         
         do k=1,cctk_lsh(3)
            do j=1,cctk_lsh(2)
               do i=1,cctk_lsh(1)
                  
                  x = cart3d_x(i,j,k)
                  y = cart3d_y(i,j,k)
                  z = cart3d_z(i,j,k)

                  phi(i,j,k)     = 1 / sqrt((x-cx)**2 + (y-cy)**2 + (z-cz)**2)
                  phi_p(i,j,k)   = phi(i,j,k)
                  phi_p_p(i,j,k) = phi(i,j,k)
                  
               end do
            end do
         end do 
         
      else

         do k=1,cctk_lsh(3)
            do j=1,cctk_lsh(2)
               do i=1,cctk_lsh(1)
                  
                  phi(i,j,k)     = 0.0d0
                  phi_p(i,j,k)   = 0.0d0
                  phi_p_p(i,j,k) = 0.0d0

               end do
            end do
         end do 
         
      end if
      
      end 
