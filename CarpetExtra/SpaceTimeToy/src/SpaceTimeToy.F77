c     -*-Fortran-*-
c     $Header: /home/eschnett/C/carpet/Carpet/CarpetExtra/SpaceTimeToy/src/SpaceTimeToy.F77,v 1.11 2003/10/27 15:31:40 schnetter Exp $
      
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
      
      
      
      subroutine SpaceTimeToy_EulerStep (CCTK_ARGUMENTS)
      
      implicit none
      
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_FUNCTIONS
      DECLARE_CCTK_PARAMETERS
      
      integer i,j,k
      
c     Copy
      do k=1,cctk_lsh(3)
         do j=1,cctk_lsh(2)
            do i=1,cctk_lsh(1)
               
               phi_i(i,j,k) = phi_p(i,j,k)
               psi_i(i,j,k) = psi_p(i,j,k)
               
            end do
         end do
      end do
      
      if (hydrotoy_active.eq.1) then
         
         do k=1,cctk_lsh(3)
            do j=1,cctk_lsh(2)
               do i=1,cctk_lsh(1)
               
                  u_i(i,j,k)  = u_p(i,j,k)
                  vx_i(i,j,k) = vx_p(i,j,k)
                  vy_i(i,j,k) = vy_p(i,j,k)
                  vz_i(i,j,k) = vz_p(i,j,k)
                  
               end do
            end do
         end do
         
      else
         
         do k=1,cctk_lsh(3)
            do j=1,cctk_lsh(2)
               do i=1,cctk_lsh(1)
                  
                  u_i(i,j,k)  = 0
                  vx_i(i,j,k) = 0
                  vy_i(i,j,k) = 0
                  vz_i(i,j,k) = 0
                  
               end do
            end do
         end do
         
      end if
      
c     Evolve
      call SpaceTimeToy_Step (CCTK_PASS_FTOF)
      
c     Initialise ICN iterations
      icn_iteration = 0
      do_iterate = 1
      if (icn_iteration .eq. icn_iterations) then
         do_iterate = 0
      end if
      
      end
      
      
      
      subroutine SpaceTimeToy_ICNStep (CCTK_ARGUMENTS)
      
      implicit none
      
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_FUNCTIONS
      DECLARE_CCTK_PARAMETERS
      
      CCTK_REAL two, half
      parameter (two=2, half=1/two)
      
      integer i,j,k
      
c     Average
      do k=1,cctk_lsh(3)
         do j=1,cctk_lsh(2)
            do i=1,cctk_lsh(1)
               
               phi_i(i,j,k) = half * (phi_p(i,j,k) + phi(i,j,k))
               psi_i(i,j,k) = half * (psi_p(i,j,k) + psi(i,j,k))
               
            end do
         end do
      end do
      
      if (hydrotoy_active.eq.1) then
         
         do k=1,cctk_lsh(3)
            do j=1,cctk_lsh(2)
               do i=1,cctk_lsh(1)
                  
                  u_i(i,j,k)  = half * (u_p(i,j,k)  + u(i,j,k))
                  vx_i(i,j,k) = half * (vx_p(i,j,k) + vx(i,j,k))
                  vy_i(i,j,k) = half * (vy_p(i,j,k) + vy(i,j,k))
                  vz_i(i,j,k) = half * (vz_p(i,j,k) + vz(i,j,k))
                  
               end do
            end do
         end do
         
      else
         
         do k=1,cctk_lsh(3)
            do j=1,cctk_lsh(2)
               do i=1,cctk_lsh(1)
                  
                  u_i(i,j,k)  = 0
                  vx_i(i,j,k) = 0
                  vy_i(i,j,k) = 0
                  vz_i(i,j,k) = 0
                  
               end do
            end do
         end do
         
      end if
      
c     Evolve
      call SpaceTimeToy_Step (CCTK_PASS_FTOF)
      
c     Step ICN iterations
      icn_iteration = icn_iteration + 1
      if (icn_iteration .eq. icn_iterations) then
         do_iterate = 0
      end if
      
      end
      
      
      
      subroutine SpaceTimeToy_Step (CCTK_ARGUMENTS)
      
      implicit none
      
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_FUNCTIONS
      DECLARE_CCTK_PARAMETERS
      
      CCTK_REAL dx,dy,dz,dt
      integer i,j,k
      
      dx = CCTK_DELTA_SPACE(1)
      dy = CCTK_DELTA_SPACE(2)
      dz = CCTK_DELTA_SPACE(3)
      dt = CCTK_DELTA_TIME
      
c     Evolve
      do k=1+cctk_nghostzones(3),cctk_lsh(3)-cctk_nghostzones(3)
         do j=1+cctk_nghostzones(2),cctk_lsh(2)-cctk_nghostzones(2)
            do i=1+cctk_nghostzones(1),cctk_lsh(1)-cctk_nghostzones(1)
               
               phi(i,j,k) = phi_p(i,j,k)
     $              + dt * psi_i(i,j,k)
     $              + dt * u_i(i,j,k)
               
               psi(i,j,k) = psi_p(i,j,k)
     $              + dt * (phi_i(i-1,j,k) - 2*phi_i(i,j,k) + phi_i(i+1,j,k)) / dx**2
     $              + dt * (phi_i(i,j-1,k) - 2*phi_i(i,j,k) + phi_i(i,j+1,k)) / dy**2
     $              + dt * (phi_i(i,j,k-1) - 2*phi_i(i,j,k) + phi_i(i,j,k+1)) / dz**2
     $              - dt * (vx_i(i+1,j,k) - vx_i(i-1,j,k)) / (2*dx)
     $              - dt * (vy_i(i,j+1,k) - vy_i(i,j-1,k)) / (2*dy)
     $              - dt * (vz_i(i,j,k+1) - vz_i(i,j,k-1)) / (2*dz)
               
            end do
         end do
      end do
      
      end
      
      
      
      subroutine SpaceTimeToy_Boundaries (CCTK_ARGUMENTS)
      
      implicit none
      
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_FUNCTIONS
      DECLARE_CCTK_PARAMETERS
      
      character fbound*1000
      CCTK_INT fboundlen
      
      integer options
      
      CCTK_INT boundary_width
      CCTK_INT options1
      
      integer d
      integer ierr
      
      boundary_width = cctk_nghostzones(1)
      do d=1,3
         if (cctk_nghostzones(d) .ne. boundary_width) then
            call CCTK_WARN (0, "internal error")
         end if
      end do
      
      call Util_TableCreateFromString (options, "")
      if (options .lt. 0) call CCTK_WARN (0, "internal error")
      
      call CCTK_FortranString (fboundlen, bound, fbound)
      
      options1 = options
      ierr = Boundary_SelectGroupForBC (cctkGH, CCTK_ALL_FACES, boundary_width, options1, "spacetimetoy::spacetimeevolve", fbound)
      if (ierr .lt. 0) then
         call CCTK_WARN (0, "Error while selecting boundary condition")
      end if
      
      call Util_TableDestroy (ierr, options)
      if (ierr .lt. 0) call CCTK_WARN (0, "internal error")
      
      call Cart3dSymGN (ierr, cctkGH, "spacetimetoy::spacetimeevolve")
      if (ierr .lt. 0) then
         call CCTK_WARN (0, "Error while applying symmetry condition")
      end if
    
      end
