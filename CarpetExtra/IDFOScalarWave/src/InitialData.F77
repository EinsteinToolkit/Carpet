c     -*-Fortran-*-
c     $Header: /home/eschnett/C/carpet/Carpet/CarpetExtra/IDFOScalarWave/src/InitialData.F77,v 1.4 2002/02/26 19:43:32 shawley Exp $

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
      CCTK_REAL  dxi,dyi,dzi
      CCTK_REAL  tmp
      
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

                  x = cart3d_x(i,j,k)
                  y = cart3d_y(i,j,k)
                  z = cart3d_z(i,j,k)

                  phi(i,j,k) = amplitude 
     $                 * exp(- (x - radius)**2 / sigma**2)
     $                 * exp(- (y - radius)**2 / sigma**2)
     $                 * exp(- (z - radius)**2 / sigma**2)

                  pi(i,j,k) = 0.0
                  phix(i,j,k) = phi(i,j,k)*  (-2) * (x - radius) / sigma**2
                  phiy(i,j,k) = phi(i,j,k)*  (-2) * (y - radius) / sigma**2
                  phiz(i,j,k) = phi(i,j,k)*  (-2) * (z - radius) / sigma**2
                  resid(i,j,k) = 0.0

                  pi_p(i,j,k) = amplitude *
     &                    exp( -(x - radius - dt)**2 / sigma**2 )* 
     &                    exp( -(y - radius - dt)**2 / sigma**2 )* 
     &                    exp( -(z - radius - dt)**2 / sigma**2 )* 
     &               ((x-radius-dt) + (y-radius-dt) + (z-radius-dt))/sigma**2
     &                        - amplitude * 
     &                    exp( -(x - radius + dt)**2 / sigma**2 )* 
     &                    exp( -(y - radius + dt)**2 / sigma**2 )* 
     &                    exp( -(z - radius + dt)**2 / sigma**2 )* 
     &               ((x-radius+dt) + (y-radius+dt) + (z-radius+dt))/sigma**2

                  pi_p_p(i,j,k) = amplitude * 
     &                    exp( -(x - radius - 2*dt)**2 / sigma**2 )* 
     &                    exp( -(y - radius - 2*dt)**2 / sigma**2 )* 
     &                    exp( -(z - radius - 2*dt)**2 / sigma**2 )* 
     &          ((x-radius-2*dt) + (y-radius-2*dt) + (z-radius-2*dt))/sigma**2
     &                        - amplitude * 
     &                    exp( -(x - radius + 2*dt)**2 / sigma**2 )* 
     &                    exp( -(y - radius + 2*dt)**2 / sigma**2 )* 
     &                    exp( -(z - radius + 2*dt)**2 / sigma**2 )* 
     &          ((x-radius+2*dt) + (y-radius+2*dt) + (z-radius+2*dt))/sigma**2

                  phix_p(i,j,k) = - amplitude * (x - radius - dt) / sigma**2 
     &                 *  exp( -(x - radius - dt)**2 / sigma**2 )* 
     &                    exp( -(y - radius - dt)**2 / sigma**2 )* 
     &                    exp( -(z - radius - dt)**2 / sigma**2 )
     &                            - amplitude * (x - radius + dt) / sigma**2
     &                 *  exp( -(x - radius + dt)**2 / sigma**2 )* 
     &                    exp( -(y - radius + dt)**2 / sigma**2 )*
     &                    exp( -(z - radius + dt)**2 / sigma**2 )
     
                  phix_p_p(i,j,k) = - amplitude * (x - radius - 2*dt) / sigma**2
     &                 *  exp( -(x - radius - 2*dt)**2 / sigma**2 )* 
     &                    exp( -(y - radius - 2*dt)**2 / sigma**2 )* 
     &                    exp( -(z - radius - 2*dt)**2 / sigma**2 )
     &                            - amplitude * (x - radius + 2*dt) / sigma**2
     &                 *  exp( -(x - radius + 2*dt)**2 / sigma**2 )* 
     &                    exp( -(y - radius + 2*dt)**2 / sigma**2 )*
     &                    exp( -(z - radius + 2*dt)**2 / sigma**2 )

                  phiy_p(i,j,k) = - amplitude * (y - radius - dt) / sigma**2
     &                 *  exp( -(x - radius - dt)**2 / sigma**2 )*
     &                    exp( -(y - radius - dt)**2 / sigma**2 )*
     &                    exp( -(z - radius - dt)**2 / sigma**2 ) 
     &                            - amplitude * (y - radius + dt) / sigma**2
     &                 *  exp( -(x - radius + dt)**2 / sigma**2 )*
     &                    exp( -(y - radius + dt)**2 / sigma**2 )*
     &                    exp( -(z - radius + dt)**2 / sigma**2 ) 
    
                  phiy_p_p(i,j,k) = - amplitude * (y - radius - 2*dt) / sigma**2
     &                 *  exp( -(x - radius - 2*dt)**2 / sigma**2 )*
     &                    exp( -(y - radius - 2*dt)**2 / sigma**2 )*
     &                    exp( -(z - radius - 2*dt)**2 / sigma**2 )
     &                            - amplitude * (y - radius + 2*dt) / sigma**2
     &                 *  exp( -(x - radius + 2*dt)**2 / sigma**2 )*
     &                    exp( -(y - radius + 2*dt)**2 / sigma**2 )*
     &                    exp( -(z - radius + 2*dt)**2 / sigma**2 )

                  phiz_p(i,j,k) = - amplitude * (z - radius - dt) / sigma**2
     &                 *  exp( -(x - radius - dt)**2 / sigma**2 )*
     &                    exp( -(y - radius - dt)**2 / sigma**2 )*
     &                    exp( -(z - radius - dt)**2 / sigma**2 ) 
     &                            - amplitude * (z - radius + dt) / sigma**2
     &                 *  exp( -(x - radius + dt)**2 / sigma**2 )*
     &                    exp( -(y - radius + dt)**2 / sigma**2 )*
     &                    exp( -(z - radius + dt)**2 / sigma**2 ) 
    
                  phiz_p_p(i,j,k) = - amplitude * (z - radius - 2*dt) / sigma**2
     &                 *  exp( -(x - radius - 2*dt)**2 / sigma**2 )*
     &                    exp( -(y - radius - 2*dt)**2 / sigma**2 )*
     &                    exp( -(z - radius - 2*dt)**2 / sigma**2 )
     &                            - amplitude * (z - radius + 2*dt) / sigma**2
     &                 *  exp( -(x - radius + 2*dt)**2 / sigma**2 )*
     &                    exp( -(y - radius + 2*dt)**2 / sigma**2 )*
     &                    exp( -(z - radius + 2*dt)**2 / sigma**2 )

                  phi_p(i,j,k) = amplitude / 2 *
     &                    exp( -(x - radius - dt)**2 / sigma**2 )*
     &                    exp( -(y - radius - dt)**2 / sigma**2 )*
     &                    exp( -(z - radius - dt)**2 / sigma**2 )
     &                          + amplitude / 2 *
     &                    exp( -(x - radius + dt)**2 / sigma**2 )*
     &                    exp( -(y - radius + dt)**2 / sigma**2 )*
     &                    exp( -(z - radius + dt)**2 / sigma**2 )

                  phi_p_p(i,j,k) = amplitude / 2 *
     &                    exp( -(x - radius - 2*dt)**2 / sigma**2 )*
     &                    exp( -(y - radius - 2*dt)**2 / sigma**2 )*
     &                    exp( -(z - radius - 2*dt)**2 / sigma**2 )
     &                          + amplitude / 2 *
     &                    exp( -(x - radius + 2*dt)**2 / sigma**2 )*
     &                    exp( -(y - radius + 2*dt)**2 / sigma**2 )*
     &                    exp( -(z - radius + 2*dt)**2 / sigma**2 )

                  resid_p(i,j,k) = 0.0
                  resid_p_p(i,j,k) = 0.0

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
                  
                  phi_p(i,j,k) = phi(i,j,k)
     $                 * cos(omega * (cctk_time - dt) * cpi)
                  
                  phi_p_p(i,j,k) = phi(i,j,k)
     $                 * cos(omega * (cctk_time - 2*dt) * cpi)

                  phi(i,j,k) = phi(i,j,k) * cos(omega * cctk_time * cpi)

               end do
            end do
         end do

      else if (CCTK_EQUALS(initial_data, "1/r")) then
  
         dxi = 1/CCTK_DELTA_SPACE(1)
         dyi = 1/CCTK_DELTA_SPACE(2)
         dzi = 1/CCTK_DELTA_SPACE(3)

         do k=1,cctk_lsh(3)
            do j=1,cctk_lsh(2)
               do i=1,cctk_lsh(1)
                  
                  pi(i,j,k)  = 0.0
                  phi(i,j,k) = 1 / spher3d_r(i,j,k)
                  ri3 = phi(i,j,k)**3
                  phix(i,j,k) = - cart3d_x(i,j,k) * ri3
                  phiy(i,j,k) = - cart3d_y(i,j,k) * ri3
                  phiz(i,j,k) = - cart3d_z(i,j,k) * ri3
                  resid(i,j,k) = 0.0

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
                  resid_p(i,j,k)  = resid(i,j,k)
                  resid_p_p(i,j,k)= resid(i,j,k)
                  
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
