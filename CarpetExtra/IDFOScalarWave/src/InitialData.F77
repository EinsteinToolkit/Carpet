c     -*-Fortran-*-
c     $Header: /home/eschnett/C/carpet/Carpet/CarpetExtra/IDFOScalarWave/src/InitialData.F77,v 1.1 2002/02/18 11:26:34 shawley Exp $

 /*@@
   @file      InitialData.F77
   @date      
   @author    Scott Hawley, using code from Tom Goodale, Erik Schnetter
   @desc 
              Initial data for the 3D Wave Equation
   @enddesc 
 @@*/

#include "cctk.h" 
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"


 /*@@
   @routine    IDFOScalarWave_InitialData
   @date       
   @author     Scott Hawley, using code from Tom Goodale, Erik Schnetter
   @desc 
               Set up initial data for the wave equation
   @enddesc 
   @calls      
   @calledby   
   @history 
 
   @endhistory 

@@*/

      subroutine IDFOScalarWave_InitialData (CCTK_ARGUMENTS)

      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS

      INTEGER    i,j,k
      CCTK_REAL  dt,omega, cpi
      CCTK_REAL  x,y,z, r, ri3
      
c     call CCTK_INFO ("IDFOScalarWave_InitialData")
      
      cpi = 4.d0*atan(1.d0)

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
     $                 * cos((kx*x + ky*y + kz*z + omega*cctk_time) * cpi)

                  phi_p(i,j,k) = amplitude
     $              * cos((kx*x + ky*y + kz*z + omega*(cctk_time - dt)) * cpi)

                  phi_p_p(i,j,k) = amplitude
     $             * cos((kx*x + ky*y + kz*z + omega*(cctk_time - 2*dt)) * cpi)
                  
               end do
            end do
         end do

      else if (CCTK_EQUALS(initial_data,"gaussian")) then 

         do k=1, cctk_lsh(3)
            do j=1, cctk_lsh(2)
               do i=1, cctk_lsh(1)

                  r = spher3d_r(i,j,k)

                  phi(i,j,k) = amplitude
     $                 * exp(- (r - radius)**2 / sigma**2)

                  phi_p(i,j,k) = amplitude/2 * (r - dt) / r
     $                 * exp(- (r - radius - dt)**2 / sigma**2)
     $                 + amplitude/2 * (r + dt) / r
     $                 * exp(- (r - radius + dt)**2 / sigma**2)

                  phi_p_p(i,j,k) = amplitude/2 * (r - 2*dt) / r
     $                 * exp(- (r - radius - 2*dt)**2 / sigma**2)
     $                 + amplitude/2 * (r + 2*dt) / r
     $                 * exp(- (r - radius + 2*dt)**2 / sigma**2)

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
     $                 * sin(kx * (x - 0.5d0) * cpi)
     $                 * sin(ky * (y - 0.5d0) * cpi)
     $                 * sin(kz * (z - 0.5d0) * cpi)
     $                 * cos(omega * cctk_time * cpi)
                  
                  phi_p(i,j,k) = amplitude
     $                 * sin(kx * (x - 0.5d0) * cpi)
     $                 * sin(ky * (y - 0.5d0) * cpi)
     $                 * sin(kz * (z - 0.5d0) * cpi)
     $                 * cos(omega * (cctk_time - dt) * cpi)
                  
                  phi_p_p(i,j,k) = amplitude
     $                 * sin(kx * (x - 0.5d0) * cpi)
     $                 * sin(ky * (y - 0.5d0) * cpi)
     $                 * sin(kz * (z - 0.5d0) * cpi)
     $                 * cos(omega * (cctk_time - 2*dt) * cpi)

               end do
            end do
         end do

      else if (CCTK_EQUALS(initial_data, "1/r")) then
         
         do k=1,cctk_lsh(3)
            do j=1,cctk_lsh(2)
               do i=1,cctk_lsh(1)
                  
                  pi(i,j,k)  = 0.0
                  phi(i,j,k) = 1 / spher3d_r(i,j,k)
                  ri3 = phi(i,j,k)**3
                  phix(i,j,k) = - cart3d_x(i,j,k) * ri3
                  phiy(i,j,k) = - cart3d_y(i,j,k) * ri3
                  phiz(i,j,k) = - cart3d_z(i,j,k) * ri3

                  pi_p(i,j,k)     = pi(i,j,k)
                  pi_p_p(i,j,k)   = pi(i,j,k)
                  phix_p(i,j,k)   = phix(i,j,k)
                  phix_p_p(i,j,k) = phix(i,j,k)
                  phiy_p(i,j,k)   = phiy(i,j,k)
                  phiy_p_p(i,j,k) = phiy(i,j,k)
                  phiz_p(i,j,k)   = phiz(i,j,k)
                  phiz_p_p(i,j,k) = phiz(i,j,k)
                  phi_p(i,j,k)    = phi(i,j,k)
                  phi_p_p(i,j,k)  = phi(i,j,k)
                  
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
