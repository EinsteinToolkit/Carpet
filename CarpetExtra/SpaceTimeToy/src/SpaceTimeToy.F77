c     -*-Fortran-*-
c     $Header: /home/eschnett/C/carpet/Carpet/CarpetExtra/SpaceTimeToy/src/SpaceTimeToy.F77,v 1.4 2001/03/20 16:29:55 eschnett Exp $
      
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
               
               phi_i(i,j,k) = phi(i,j,k)
               psi_i(i,j,k) = psi(i,j,k)
               
               u_i(i,j,k)  = u(i,j,k)
               vx_i(i,j,k) = vx(i,j,k)
               vy_i(i,j,k) = vy(i,j,k)
               vz_i(i,j,k) = vz(i,j,k)
               
            end do
         end do
      end do
      
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
               
               phi_i(i,j,k) = half * (phi(i,j,k) + phi_n(i,j,k))
               psi_i(i,j,k) = half * (psi(i,j,k) + psi_n(i,j,k))
               
               u_i(i,j,k)  = half * (u(i,j,k)  + u_n(i,j,k))
               vx_i(i,j,k) = half * (vx(i,j,k) + vx_n(i,j,k))
               vy_i(i,j,k) = half * (vy(i,j,k) + vy_n(i,j,k))
               vz_i(i,j,k) = half * (vz(i,j,k) + vz_n(i,j,k))
               
            end do
         end do
      end do
      
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
      do k=2,cctk_lsh(3)-1
         do j=2,cctk_lsh(2)-1
            do i=2,cctk_lsh(1)-1
               
               phi_n(i,j,k) = phi(i,j,k)
     $              + dt * psi_i(i,j,k)
     $              + dt * u_i(i,j,k)
               
               psi_n(i,j,k) = psi(i,j,k)
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
      
      CCTK_REAL finf, npow
      parameter (finf=1, npow=1)
      
      integer sw(3)
      integer ierr
      
      sw(1) = 1
      sw(2) = 1
      sw(3) = 1
      
c     Apply boundary condition
      if (CCTK_EQUALS(bound, "flat")) then
         call BndFlatGN (ierr, cctkGH, sw, "spacetimetoy::spacetimeevolve")
      else if (CCTK_EQUALS(bound, "zero")) then
         call BndScalarGN (ierr, cctkGH, zero, sw,
     $        "spacetimetoy::spacetimeevolve")
      else if (CCTK_EQUALS(bound, "radiation")) then
         call BndRadiativeGN (ierr, cctkGH, sw, zero, one,
     $        "spacetimetoy::spacetimeevolve", "spacetimetoy::spacetimeevolve")
      else if (CCTK_EQUALS(bound, "robin")) then
         call BndRobinGN (ierr, cctkGH, sw, finf, npow,
     $        "spacetimetoy::spacetimeevolve")
      else
         call CCTK_WARN (0, "Internal error")
      end if
      if (ierr .lt. 0)  then
         call CCTK_WARN (0, "Error while applying boundary condition")
      end if
      
      call Cart3dSymGN (ierr, cctkGH, "spacetimetoy::spacetimeevolve")
      if (ierr .lt. 0)  then
         call CCTK_WARN (0, "Error while applying boundary condition")
      end if
      
      end
