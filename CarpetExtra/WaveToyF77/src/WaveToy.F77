c     -*-Fortran-*-

 /*@@
   @file      WaveToy.F77
   @date      
   @author    Tom Goodale, Erik Schnetter
   @desc 
              Evolution routines for the wave equation solver
   @enddesc 
 @@*/
  
#include "cctk.h" 
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"



 /*@@
   @routine    WaveToyF77_Evolution
   @date       
   @author     Tom Goodale, Erik Schnetter
   @desc 
               Evolution for the wave equation
   @enddesc 
   @calls      CCTK_SyncGroup
   @calledby   
   @history 
 
   @endhistory 

@@*/

      subroutine WaveToyF77_Evolution (CCTK_ARGUMENTS)
      
      implicit none
      
c     Declare variables in argument list
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_FUNCTIONS
      DECLARE_CCTK_PARAMETERS
      
      INTEGER   i,j,k
      INTEGER   istart, jstart, kstart, iend, jend, kend
      CCTK_REAL dx,dy,dz,dt
      CCTK_REAL dx2,dy2,dz2,dt2
      CCTK_REAL dx2i,dy2i,dz2i
      
      CCTK_REAL factor
      
c     call CCTK_INFO ("WaveToyF77_Evolution")
      
c     Set up shorthands
c     -----------------
      dx = CCTK_DELTA_SPACE(1)
      dy = CCTK_DELTA_SPACE(2)
      dz = CCTK_DELTA_SPACE(3)
      dt = CCTK_DELTA_TIME
      
      dx2 = dx**2
      dy2 = dy**2
      dz2 = dz**2
      dt2 = dt**2
      
      dx2i = 1/dx2
      dy2i = 1/dy2
      dz2i = 1/dz2
      
      istart = 1+cctk_nghostzones(1)
      jstart = 1+cctk_nghostzones(2)
      kstart = 1+cctk_nghostzones(3)
      
      iend = cctk_lsh(1)-cctk_nghostzones(1)
      jend = cctk_lsh(2)-cctk_nghostzones(2)
      kend = cctk_lsh(3)-cctk_nghostzones(3)
      
      factor = 2 * (1 - dt2 * (dx2i + dy2i + dz2i))
      
c     Do the evolution
c     ----------------
      do k = kstart, kend
         do j = jstart, jend
            do i = istart, iend
               
               phi(i,j,k) = factor*phi_p(i,j,k) -
     $              phi_p_p(i,j,k) + (dt2) * 
     $              ((phi_p(i+1,j,k)+phi_p(i-1,j,k))*dx2i
     $              +(phi_p(i,j+1,k)+phi_p(i,j-1,k))*dy2i
     $              +(phi_p(i,j,k+1)+phi_p(i,j,k-1))*dz2i)
               
            end do
         end do
      end do
      
      end



 /*@@
   @routine    WaveToyF77_Boundaries
   @date       
   @author     Tom Goodale, Erik Schnetter
   @desc 
               Boundary conditions for the wave equation
   @enddesc 
   @history 
 
   @endhistory 

@@*/

      subroutine WaveToyF77_Boundaries (CCTK_ARGUMENTS)
      
      implicit none
      
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_FUNCTIONS
      DECLARE_CCTK_PARAMETERS
      
c     Local declarations
      CCTK_REAL zero, one
      parameter (zero=0, one=1)
      
      CCTK_REAL finf
      integer npow
      parameter (finf = 1)
      parameter (npow = 1)
      
      integer i,j,k
      
      integer ierr
      integer sw(3)
      
c     call CCTK_INFO ("WaveToyF77_Boundaries")
      
c     Set the stencil width
c     ---------------------
      sw(1) = cctk_nghostzones(1)
      sw(2) = cctk_nghostzones(2)
      sw(3) = cctk_nghostzones(3)
      
c     Apply the excision boundary condition
c     -------------------------------------
      if (CCTK_EQUALS(excision_bound, "none")) then
c     do nothing
      else if (CCTK_EQUALS(excision_bound, "1/r")) then
         do k=1,cctk_lsh(3)
            do j=1,cctk_lsh(2)
               do i=1,cctk_lsh(1)
                  if (spher3d_r(i,j,k) .le. excision_radius) then
                     phi(i,j,k) = 1 / spher3d_r(i,j,k)
                  end if
               end do
            end do
         end do 
      else
         call CCTK_WARN (0, "internal error")
      end if
      
c     Apply the outer boundary conditions
c     -----------------------------------
      if (CCTK_EQUALS(bound, "flat")) then
         call BndFlatVN (ierr, cctkGH, sw, "wavetoy::phi")
      else if (CCTK_EQUALS(bound, "zero")) then
         call BndScalarVN (ierr, cctkGH, sw, zero, "wavetoy::phi")
      else if (CCTK_EQUALS(bound, "static")) then
         call BndStaticVN (ierr, cctkGH, sw, 1, 0, "wavetoy::phi")
      else if (CCTK_EQUALS(bound, "radiation")) then
         call BndRadiativeVN (ierr, cctkGH, sw, zero, one,
     $        "wavetoy::phi", "wavetoy::phi")
      else if (CCTK_EQUALS(bound, "robin")) then
         call BndRobinVN (ierr, cctkGH, sw, finf, npow, "wavetoy::phi")

      else if (CCTK_EQUALS(bound, "none")) then
         ierr = 0
      else
         call CCTK_WARN (0, "internal error")
      end if
      if (ierr .lt. 0)  then
        call CCTK_WARN (0, "Boundary conditions not applied - giving up!")
      end if
      
c     Apply the symmetry boundary conditions on any coordinate axes
c     -------------------------------------------------------------
      call Cart3dSymGN (ierr, cctkGH, "wavetoy::scalarevolve")
      
      if (ierr .lt. 0)  then
        call CCTK_WARN (0, "Symmetry conditions not applied - giving up!")
      end if
      
      end
