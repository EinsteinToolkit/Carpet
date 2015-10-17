c     -*-Fortran-*-
c     $Header:$
      
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
      
      subroutine IDHydroToy_InitialData (CCTK_ARGUMENTS)
      
      implicit none
      
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS
      
      CCTK_REAL pi
      CCTK_REAL omega
      CCTK_REAL dt
      CCTK_REAL x,y,z, r
      integer i,j,k
      
      CCTK_REAL vr
      
      external erf
      real*8 erf
      
      pi = 4*atan(1.d0)
      
      omega = sqrt(kx**2+ky**2+kz**2)
      
      dt = CCTK_DELTA_TIME
      
      if (CCTK_EQUALS(initial_data,"plane")) then 
         
         do k=1,cctk_lsh(3)
            do j=1,cctk_lsh(2)
               do i=1,cctk_lsh(1)
                  
                  x = cart3d_x(i,j,k)
                  y = cart3d_y(i,j,k)
                  z = cart3d_z(i,j,k)
                  
                  u(i,j,k) = amplitude
     $                 * cos((kx*x + ky*y + kz*z + omega*cctk_time) * pi)
                  vx(i,j,k) = u(i,j,k) * kx / omega
                  vy(i,j,k) = u(i,j,k) * ky / omega
                  vz(i,j,k) = u(i,j,k) * kz / omega
                  
                  u_p(i,j,k) = amplitude
     $                 * cos((kx*x + ky*y + kz*z + omega*(cctk_time-dt)) * pi)
                  vx_p(i,j,k) = u_p(i,j,k) * kx / omega
                  vy_p(i,j,k) = u_p(i,j,k) * ky / omega
                  vz_p(i,j,k) = u_p(i,j,k) * kz / omega
                  
                  u_p_p(i,j,k) = amplitude
     $                 * cos((kx*x + ky*y + kz*z + omega*(cctk_time-2*dt)) * pi)
                  vx_p_p(i,j,k) = u_p_p(i,j,k) * kx / omega
                  vy_p_p(i,j,k) = u_p_p(i,j,k) * ky / omega
                  vz_p_p(i,j,k) = u_p_p(i,j,k) * kz / omega
                  
               end do
            end do
         end do
         
      else if (CCTK_EQUALS(initial_data,"gaussian")) then 
         
         do k=1,cctk_lsh(3)
            do j=1,cctk_lsh(2)
               do i=1,cctk_lsh(1)
                  
                  x = cart3d_x(i,j,k)
                  y = cart3d_y(i,j,k)
                  z = cart3d_z(i,j,k)
                  r = spher3d_r(i,j,k)
                  
                  u(i,j,k) = amplitude
     $                 * exp(- (r - radius)**2 / sigma**2)
                  
                  vr = - 2*amplitude * (r - radius) / sigma**2
     $                 * exp(- (r - radius)**2 / sigma**2)
                  vx(i,j,k) = vr * x/r
                  vy(i,j,k) = vr * y/r
                  vz(i,j,k) = vr * z/r
                  
                  u_p(i,j,k) = amplitude/2 * (r - dt) / r
     $                 * exp(- (r - radius - dt)**2 / sigma**2)
     $                 + amplitude/2 * (r + dt) / r
     $                 * exp(- (r - radius + dt)**2 / sigma**2)
                  
                  vr = - amplitude/2 * (-dt / r**2 + (r - dt) * (r - radius - dt) / (r * sigma**2))
     $                 * exp(- (r - radius - dt)**2 / sigma**2)
     $                 - amplitude/2 * ( dt / r**2 + (r + dt) * (r - radius + dt) / (r * sigma**2))
     $                 * exp(- (r - radius - dt)**2 / sigma**2)
                  vx_p(i,j,k) = vr * x/r
                  vy_p(i,j,k) = vr * y/r
                  vz_p(i,j,k) = vr * z/r
                  
                  u_p_p(i,j,k) = amplitude/2 * (r - 2*dt) / r
     $                 * exp(- (r - radius - 2*dt)**2 / sigma**2)
     $                 + amplitude/2 * (r + 2*dt) / r
     $                 * exp(- (r - radius + 2*dt)**2 / sigma**2)
                  
                  vr = - amplitude/2 * (-2*dt / r**2 + (r - 2*dt) * (r - radius - 2*dt) / (r * sigma**2))
     $                 * exp(- (r - radius - 2*dt)**2 / sigma**2)
     $                 - amplitude/2 * ( 2*dt / r**2 + (r + 2*dt) * (r - radius + 2*dt) / (r * sigma**2))
     $                 * exp(- (r - radius - 2*dt)**2 / sigma**2)
                  vx_p_p(i,j,k) = vr * x/r
                  vy_p_p(i,j,k) = vr * y/r
                  vz_p_p(i,j,k) = vr * z/r
                  
               end do
            end do
         end do
         
      else if (CCTK_EQUALS(initial_data, "box")) then
         
c     Use kx,ky,kz as number of modes in each direction.
         
         do k=1,cctk_lsh(3)
            do j=1,cctk_lsh(2)
               do i=1,cctk_lsh(1)
                  
                  x = cart3d_x(i,j,k)
                  y = cart3d_y(i,j,k)
                  z = cart3d_z(i,j,k)
                  
                  u(i,j,k) = amplitude
     $                 * sin(kx * (x - 0.5d0) * pi)
     $                 * sin(ky * (y - 0.5d0) * pi)
     $                 * sin(kz * (z - 0.5d0) * pi)
     $                 * cos(omega * cctk_time * pi)
                  vx(i,j,k) = amplitude
     $                 * cos(kx * (x - 0.5d0) * pi)
     $                 * sin(ky * (y - 0.5d0) * pi)
     $                 * sin(kz * (z - 0.5d0) * pi)
     $                 * sin(omega * cctk_time * pi)
     $                 * kx / omega
                  vy(i,j,k) = amplitude
     $                 * sin(kx * (x - 0.5d0) * pi)
     $                 * cos(ky * (y - 0.5d0) * pi)
     $                 * sin(kz * (z - 0.5d0) * pi)
     $                 * sin(omega * cctk_time * pi)
     $                 * ky / omega
                  vz(i,j,k) = amplitude
     $                 * sin(kx * (x - 0.5d0) * pi)
     $                 * sin(ky * (y - 0.5d0) * pi)
     $                 * cos(kz * (z - 0.5d0) * pi)
     $                 * sin(omega * cctk_time * pi)
     $                 * kz / omega
                  
                  u_p(i,j,k) = amplitude
     $                 * sin(kx * (x - 0.5d0) * pi)
     $                 * sin(ky * (y - 0.5d0) * pi)
     $                 * sin(kz * (z - 0.5d0) * pi)
     $                 * cos(omega * (cctk_time - dt) * pi)
                  vx_p(i,j,k) = amplitude
     $                 * cos(kx * (x - 0.5d0) * pi)
     $                 * sin(ky * (y - 0.5d0) * pi)
     $                 * sin(kz * (z - 0.5d0) * pi)
     $                 * sin(omega * (cctk_time - dt) * pi)
     $                 * kx / omega
                  vy_p(i,j,k) = amplitude
     $                 * sin(kx * (x - 0.5d0) * pi)
     $                 * cos(ky * (y - 0.5d0) * pi)
     $                 * sin(kz * (z - 0.5d0) * pi)
     $                 * sin(omega * (cctk_time - dt) * pi)
     $                 * ky / omega
                  vz_p(i,j,k) = amplitude
     $                 * sin(kx * (x - 0.5d0) * pi)
     $                 * sin(ky * (y - 0.5d0) * pi)
     $                 * cos(kz * (z - 0.5d0) * pi)
     $                 * sin(omega * (cctk_time - dt) * pi)
     $                 * kz / omega
                  
                  u_p_p(i,j,k) = amplitude
     $                 * sin(kx * (x - 0.5d0) * pi)
     $                 * sin(ky * (y - 0.5d0) * pi)
     $                 * sin(kz * (z - 0.5d0) * pi)
     $                 * cos(omega * (cctk_time - 2*dt) * pi)
                  vx_p_p(i,j,k) = amplitude
     $                 * cos(kx * (x - 0.5d0) * pi)
     $                 * sin(ky * (y - 0.5d0) * pi)
     $                 * sin(kz * (z - 0.5d0) * pi)
     $                 * sin(omega * (cctk_time - 2*dt) * pi)
     $                 * kx / omega
                  vy_p_p(i,j,k) = amplitude
     $                 * sin(kx * (x - 0.5d0) * pi)
     $                 * cos(ky * (y - 0.5d0) * pi)
     $                 * sin(kz * (z - 0.5d0) * pi)
     $                 * sin(omega * (cctk_time - 2*dt) * pi)
     $                 * ky / omega
                  vz_p_p(i,j,k) = amplitude
     $                 * sin(kx * (x - 0.5d0) * pi)
     $                 * sin(ky * (y - 0.5d0) * pi)
     $                 * cos(kz * (z - 0.5d0) * pi)
     $                 * sin(omega * (cctk_time - 2*dt) * pi)
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
                  
                  u_p(i,j,k)  = 0
                  vx_p(i,j,k) = 0
                  vy_p(i,j,k) = 0
                  vz_p(i,j,k) = 0
                  
                  u_p_p(i,j,k)  = 0
                  vx_p_p(i,j,k) = 0
                  vy_p_p(i,j,k) = 0
                  vz_p_p(i,j,k) = 0
                  
               end do
            end do
         end do 
         
      end if
      
      end 
