c     -*-Fortran-*-
c     $Header: /home/eschnett/C/carpet/Carpet/CarpetExtra/SpaceTimeToy/src/SpaceTimeToy.F77,v 1.9 2001/08/26 13:59:58 schnetter Exp $
      
#include "cctk.h" 
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
      
      
      
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
      
c     Apply boundaries
      call SpaceTimeToy_Boundaries (CCTK_PASS_FTOF)
      
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
      
c     Apply boundaries
      call SpaceTimeToy_Boundaries (CCTK_PASS_FTOF)
      
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
      
      CCTK_REAL zero, one
      parameter (zero=0, one=1)
      
      CCTK_REAL finf
      parameter (finf=1)
      integer npow
      parameter (npow=1)
      
      integer sw(3)
      integer ierr
      
      sw(1) = cctk_nghostzones(1)
      sw(2) = cctk_nghostzones(2)
      sw(3) = cctk_nghostzones(3)
      
c     Apply boundary condition
      if (CCTK_EQUALS(bound, "flat")) then
         call BndFlatGN (ierr, cctkGH, sw, "spacetimetoy::spacetimeevolve")
      else if (CCTK_EQUALS(bound, "zero")) then
         call BndScalarGN (ierr, cctkGH, sw, zero,
     $        "spacetimetoy::spacetimeevolve")
      else if (CCTK_EQUALS(bound, "radiation")) then
         call BndRadiativeGN (ierr, cctkGH, sw, zero, one,
     $        "spacetimetoy::spacetimeevolve", "spacetimetoy::spacetimeevolve")
      else if (CCTK_EQUALS(bound, "robin")) then
         call BndRobinGN (ierr, cctkGH, sw, finf, npow,
     $        "spacetimetoy::spacetimeevolve")
      else if (CCTK_EQUALS(bound, "none")) then
         ierr = 0
      else
         call CCTK_WARN (0, "internal error")
      end if
      if (ierr .lt. 0)  then
         call CCTK_WARN (0, "Error while applying boundary condition")
      end if
      
      call Cart3dSymGN (ierr, cctkGH, "spacetimetoy::spacetimeevolve")
      if (ierr .lt. 0)  then
         call CCTK_WARN (0, "Error while applying symmetry condition")
      end if
      
      end
