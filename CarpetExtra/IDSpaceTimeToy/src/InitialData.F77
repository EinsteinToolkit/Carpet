c     -*-Fortran-*-
c     $Header: /home/eschnett/C/carpet/Carpet/CarpetExtra/IDSpaceTimeToy/src/InitialData.F77,v 1.4 2001/03/26 02:28:54 eschnett Exp $
      
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
      
      subroutine IDSpaceTimeToy_InitialData (CCTK_ARGUMENTS)
      
      implicit none
      
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS
      
      CCTK_REAL pi
      CCTK_REAL omega
      CCTK_REAL dt
      integer i,j,k
      
      pi = 4*atan(1.d0)
      
      omega = sqrt(kx**2+ky**2+kz**2)
      
      dt = CCTK_DELTA_TIME
      
      if (CCTK_EQUALS(initial_data,"plane")) then 
         
         do k=1,cctk_lsh(3)
            do j=1,cctk_lsh(2)
               do i=1,cctk_lsh(1)
                  
                  phi(i,j,k) = amplitude
     $                 * cos((kx*cart3d_x(i,j,k) + ky*cart3d_y(i,j,k)
     $                 + kz*cart3d_z(i,j,k) + omega*cctk_time) * pi)
                  psi(i,j,k) = - amplitude
     $                 * sin((kx*cart3d_x(i,j,k) + ky*cart3d_y(i,j,k)
     $                 + kz*cart3d_z(i,j,k) + omega*cctk_time) * pi)
     $                 * pi * omega
                  
                  phi_p(i,j,k) = amplitude
     $                 * cos((kx*cart3d_x(i,j,k) + ky*cart3d_y(i,j,k)
     $                 + kz*cart3d_z(i,j,k) + omega*(cctk_time-dt)) * pi)
                  psi_p(i,j,k) = - amplitude
     $                 * sin((kx*cart3d_x(i,j,k) + ky*cart3d_y(i,j,k)
     $                 + kz*cart3d_z(i,j,k) + omega*(cctk_time-dt)) * pi)
     $                 * pi * omega
                  
               end do
            end do
         end do
         
      else if (CCTK_EQUALS(initial_data,"gaussian")) then 
         
         do k=1,cctk_lsh(3)
            do j=1,cctk_lsh(2)
               do i=1,cctk_lsh(1)
                  
                  phi(i,j,k) = amplitude / spher3d_r(i,j,k)
     $                 * exp(- (spher3d_r(i,j,k) - radius - cctk_time)**2 / sigma**2)
                  psi(i,j,k) = phi(i,j,k)
     $                 * 2 * (spher3d_r(i,j,k) - radius - cctk_time) / sigma**2
                  
                  phi_p(i,j,k) = amplitude / spher3d_r(i,j,k)
     $                 * exp(- (spher3d_r(i,j,k) - radius - (cctk_time-dt))**2 / sigma**2)
                  psi_p(i,j,k) = phi(i,j,k)
     $                 * 2 * (spher3d_r(i,j,k) - radius - (cctk_time-dt)) / sigma**2
                  
               end do
            end do
         end do
         
      else if (CCTK_EQUALS(initial_data, "box")) then
         
c     Use kx,ky,kz as number of modes in each direction.
         
         do k=1,cctk_lsh(3)
            do j=1,cctk_lsh(2)
               do i=1,cctk_lsh(1)
                  
                  phi(i,j,k) = amplitude
     $                 * sin(kx * (cart3d_x(i,j,k) - 0.5d0) * pi)
     $                 * sin(ky * (cart3d_y(i,j,k) - 0.5d0) * pi)
     $                 * sin(kz * (cart3d_z(i,j,k) - 0.5d0) * pi)
     $                 * cos(omega * cctk_time * pi)
                  psi(i,j,k) = - amplitude
     $                 * sin(kx * (cart3d_x(i,j,k) - 0.5d0) * pi)
     $                 * sin(ky * (cart3d_y(i,j,k) - 0.5d0) * pi)
     $                 * sin(kz * (cart3d_z(i,j,k) - 0.5d0) * pi)
     $                 * sin(omega * cctk_time * pi)
     $                 * omega * pi
                  
                  phi_p(i,j,k) = amplitude
     $                 * sin(kx * (cart3d_x(i,j,k) - 0.5d0) * pi)
     $                 * sin(ky * (cart3d_y(i,j,k) - 0.5d0) * pi)
     $                 * sin(kz * (cart3d_z(i,j,k) - 0.5d0) * pi)
     $                 * cos(omega * (cctk_time-dt) * pi)
                  psi_p(i,j,k) = - amplitude
     $                 * sin(kx * (cart3d_x(i,j,k) - 0.5d0) * pi)
     $                 * sin(ky * (cart3d_y(i,j,k) - 0.5d0) * pi)
     $                 * sin(kz * (cart3d_z(i,j,k) - 0.5d0) * pi)
     $                 * sin(omega * (cctk_time-dt) * pi)
     $                 * omega * pi
                  
               end do
            end do
         end do
         
      else
         
         do k=1,cctk_lsh(3)
            do j=1,cctk_lsh(2)
               do i=1,cctk_lsh(1)
                  
                  phi(i,j,k) = 0
                  psi(i,j,k) = 0
                  
                  phi_p(i,j,k) = 0
                  psi_p(i,j,k) = 0
                  
               end do
            end do
         end do 
         
      end if
      
      if (cctk_isthornactive("hydrotoy").eq.1) then
         
         do k=1,cctk_lsh(3)
            do j=1,cctk_lsh(2)
               do i=1,cctk_lsh(1)
                  
                  psi(i,j,k) = psi(i,j,k) - u(i,j,k)
                  psi_p(i,j,k) = psi_p(i,j,k) - u_p(i,j,k)
                  
               end do
            end do
         end do
         
      end if
      
      end 
