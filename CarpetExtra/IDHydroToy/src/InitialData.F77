c     -*-Fortran-*-
c     $Header: /home/eschnett/C/carpet/Carpet/CarpetExtra/IDHydroToy/src/InitialData.F77,v 1.2 2001/03/19 21:30:37 eschnett Exp $
      
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
      
      subroutine IDHydroToy_InitialData (CCTK_ARGUMENTS)
      
      implicit none
      
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS
      
      CCTK_REAL pi
      CCTK_REAL omega
      integer i,j,k
      
      CCTK_REAL vr
      
      external erf
      real*8 erf
      
      pi = 4*atan(1.d0)
      
      omega = sqrt(kx**2+ky**2+kz**2)
      
      if (CCTK_EQUALS(initial_data,"plane")) then 
         
         do k=1,cctk_lsh(3)
            do j=1,cctk_lsh(2)
               do i=1,cctk_lsh(1)
                  
                  u(i,j,k) = amplitude
     $                 * cos((kx*cart3d_x(i,j,k) + ky*cart3d_y(i,j,k)
     $                 + kz*cart3d_z(i,j,k) + omega*cctk_time) * pi)
                  vx(i,j,k) = u(i,j,k) * kx / omega
                  vy(i,j,k) = u(i,j,k) * ky / omega
                  vz(i,j,k) = u(i,j,k) * kz / omega
                  
               end do
            end do
         end do
         
      else if (CCTK_EQUALS(initial_data,"gaussian")) then 
         
         do k=1,cctk_lsh(3)
            do j=1,cctk_lsh(2)
               do i=1,cctk_lsh(1)
                  
                  u(i,j,k) = amplitude / spher3d_r(i,j,k)
     $                 * exp(- (spher3d_r(i,j,k) - radius + cctk_time)**2 / sigma**2)

c     vr = -1/2*A*(sqrt(Pi)*sigma*erf((t-r+R)/sigma)+2*exp(-(t-r+R)^2/(sigma^2))*r)/(r^2)
                  vr = -0.5d0 * amplitude
     $                 * (sqrt(pi) * sigma * erf((cctk_time - spher3d_r(i,j,k) + radius) / sigma)
     $                    + 2*exp(-(cctk_time - spher3d_r(i,j,k) + radius)**2 / sigma**2) * spher3d_r(i,j,k))
     $                 / spher3d_r(i,j,k)**2
                  vx(i,j,k) = vr * cart3d_x(i,j,k)/spher3d_r(i,j,k)
                  vy(i,j,k) = vr * cart3d_y(i,j,k)/spher3d_r(i,j,k)
                  vz(i,j,k) = vr * cart3d_z(i,j,k)/spher3d_r(i,j,k)
                  
               end do
            end do
         end do
         
      else if (CCTK_EQUALS(initial_data, "box")) then
         
c     Use kx,ky,kz as number of modes in each direction.
         
         do k=1,cctk_lsh(3)
            do j=1,cctk_lsh(2)
               do i=1,cctk_lsh(1)
                  
                  u(i,j,k) = amplitude
     $                 * sin(kx * (cart3d_x(i,j,k) - 0.5d0) * pi)
     $                 * sin(ky * (cart3d_y(i,j,k) - 0.5d0) * pi)
     $                 * sin(kz * (cart3d_z(i,j,k) - 0.5d0) * pi)
     $                 * cos(omega * cctk_time * pi)
                  
                  vx(i,j,k) = amplitude
     $                 * cos(kx * (cart3d_x(i,j,k) - 0.5d0) * pi)
     $                 * sin(ky * (cart3d_y(i,j,k) - 0.5d0) * pi)
     $                 * sin(kz * (cart3d_z(i,j,k) - 0.5d0) * pi)
     $                 * sin(omega * cctk_time * pi)
     $                 * kx / omega
                  vy(i,j,k) = amplitude
     $                 * sin(kx * (cart3d_x(i,j,k) - 0.5d0) * pi)
     $                 * cos(ky * (cart3d_y(i,j,k) - 0.5d0) * pi)
     $                 * sin(kz * (cart3d_z(i,j,k) - 0.5d0) * pi)
     $                 * sin(omega * cctk_time * pi)
     $                 * ky / omega
                  vz(i,j,k) = amplitude
     $                 * sin(kx * (cart3d_x(i,j,k) - 0.5d0) * pi)
     $                 * sin(ky * (cart3d_y(i,j,k) - 0.5d0) * pi)
     $                 * cos(kz * (cart3d_z(i,j,k) - 0.5d0) * pi)
     $                 * sin(omega * cctk_time * pi)
     $                 * kz / omega
                  
               end do
            end do
         end do
         
      else
         
         do k=1,cctk_lsh(3)
            do j=1,cctk_lsh(2)
               do i=1,cctk_lsh(1)
                  
                  u(i,j,k)  = 0
                  vx(i,j,k) = 0
                  vy(i,j,k) = 0
                  vz(i,j,k) = 0
                  
               end do
            end do
         end do 
         
      end if
      
      end 
