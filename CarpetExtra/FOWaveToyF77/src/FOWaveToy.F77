c     -*-Fortran-*-

 /*@@
   @file      FOWaveToy.F77
   @date      
   @author    Scott Hawley, based on WaveToy by Tom Goodale, Erik Schnetter
   @desc 
              Evolution routines for the wave equation solver
              written as a first order system.
              The field phi itself does not enter into the evolution,
              system but we carry it along using pi = partial phi / partial t.
   @enddesc 
 @@*/
  
#include "cctk.h" 
#include "cctk_Faces.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"



 /*@@
   @routine    FOWaveToyF77_Evolution
   @date       
   @author     Scott Hawley, using code of Tom Goodale, Erik Schnetter
   @desc 
               Evolution for the wave equation as a 1st-order system
   @enddesc 
   @calls      CCTK_SyncGroup
   @calledby   
   @history 
 
   @endhistory 

@@*/

      subroutine FOWaveToyF77_Evolution (CCTK_ARGUMENTS)
      
      implicit none
      
c     Declare variables in argument list
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_FUNCTIONS
      DECLARE_CCTK_PARAMETERS
      
      INTEGER   i,j,k
      INTEGER   istart, jstart, kstart, iend, jend, kend
      CCTK_REAL dx,dy,dz,dt
      CCTK_REAL dxi,dyi,dzi
      
c     call CCTK_INFO ("WaveToyF77_Evolution")
      
c     Set up shorthands
c     -----------------
      dx = 2*CCTK_DELTA_SPACE(1)
      dy = 2*CCTK_DELTA_SPACE(2)
      dz = 2*CCTK_DELTA_SPACE(3)
      dt = 2*CCTK_DELTA_TIME
      
      dxi = 1/dx
      dyi = 1/dy
      dzi = 1/dz
      
      istart = 1+cctk_nghostzones(1)
      jstart = 1+cctk_nghostzones(2)
      kstart = 1+cctk_nghostzones(3)
      
      iend = cctk_lsh(1)-cctk_nghostzones(1)
      jend = cctk_lsh(2)-cctk_nghostzones(2)
      kend = cctk_lsh(3)-cctk_nghostzones(3)
      
c     Do the evolution
c     ----------------
      do k = kstart, kend
         do j = jstart, jend
            do i = istart, iend
               
               pi(i,j,k) = pi_p_p(i,j,k) + dt * ( 
     $                      (phix_p(i+1,j,k) - phix_p(i-1,j,k))*dxi
     $                    + (phiy_p(i,j+1,k) - phiy_p(i,j-1,k))*dyi
     $                    + (phiz_p(i,j,k+1) - phiz_p(i,j,k-1))*dzi
     $                                          )

               phix(i,j,k) = phix_p_p(i,j,k) + dt * dxi * ( 
     $                       pi_p(i+1,j,k) - pi_p(i-1,j,k)  )

               phiy(i,j,k) = phiy_p_p(i,j,k) + dt * dyi * ( 
     $                       pi_p(i,j+1,k) - pi_p(i,j-1,k)  )

               phiz(i,j,k) = phiz_p_p(i,j,k) + dt * dzi * ( 
     $                       pi_p(i,j,k+1) - pi_p(i,j,k-1)  )

               phi(i,j,k) = phi_p_p(i,j,k) + 2*dt * pi_p(i,j,i) 
               
            end do
         end do
      end do
      
      end



 /*@@
   @routine    FOWaveToyF77_Boundaries
   @date       
   @author     Scott Hawley stealing from Tom Goodale, Erik Schnetter
   @desc 
               Boundary conditions for the wave equation
   @enddesc 
   @history 
 
   @endhistory 

@@*/

      subroutine FOWaveToyF77_Boundaries (CCTK_ARGUMENTS)
      
      implicit none
      
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS
      CCTK_INT Boundary_SelectVarForBC 
      
c     Local declarations
      CCTK_REAL zero, one
      parameter (zero=0, one=1)
      
      CCTK_REAL finf
      integer npow
      parameter (finf = 1)
      parameter (npow = 1)
      
      integer i,j,k
      CCTK_REAL  spher3d_r
      
      integer ierr
      integer sw(3)
      CCTK_REAL ri3
      
c     call CCTK_INFO ("FOWaveToyF77_Boundaries")
      
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
                  spher3d_r = sqrt(x(i,j,k)**2 + y(i,j,k)**2 + z(i,j,k)**2)
                  if (spher3d_r .le. excision_radius) then
                     pi(i,j,k)  = 0.0
                     phi(i,j,k) = 1.0 / spher3d_r
                     ri3 = phi(i,j,k)**3
                     phix(i,j,k) = - x(i,j,k) * ri3
                     phiy(i,j,k) = - y(i,j,k) * ri3
                     phiz(i,j,k) = - z(i,j,k) * ri3
                  end if
               end do
            end do
         end do 
      else
         call CCTK_WARN (0, "internal error")
      end if

c     Apply the symmetry boundary conditions on any coordinate axes
c     -------------------------------------------------------------
      call CartSymGN(ierr,cctkGH,"wavetoy::scalarevolve")
      call CartSymGN(ierr,cctkGH,"wavetoy::scalarevolve_derivs")

      
c     Apply the outer boundary conditions
c     Only "flat" and "zero" and "none" are currently supported
c     -----------------------------------
      write(0,*) ' FOWaveToyF77: bound = ', bound

        call CCTK_PrintString(bound)
      if (CCTK_EQUALS(bound, "flat")) then
         write(0,*) ' FOWavetoy: bound is flat'
         call BndFlatVN (ierr, cctkGH, sw, "wavetoy::pi")
         call BndScalarVN (ierr, cctkGH, sw, zero, "wavetoy::phix")
         call BndScalarVN (ierr, cctkGH, sw, zero, "wavetoy::phiy")
         call BndScalarVN (ierr, cctkGH, sw, zero, "wavetoy::phiz")
         call BndFlatVN (ierr, cctkGH, sw, "wavetoy::phi")
      else if (CCTK_EQUALS(bound, "zero")) then
         write(0,*) ' FOWavetoy: bound is zero'
         call BndScalarVN (ierr, cctkGH, sw, zero, "wavetoy::pi")
         call BndScalarVN (ierr, cctkGH, sw, zero, "wavetoy::phix")
         call BndScalarVN (ierr, cctkGH, sw, zero, "wavetoy::phiy")
         call BndScalarVN (ierr, cctkGH, sw, zero, "wavetoy::phiz")
         call BndScalarVN (ierr, cctkGH, sw, zero, "wavetoy::phi")
      else if (CCTK_EQUALS(bound, "static")) then
         write(0,*) ' FOWavetoy: bound is static'
         ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
     $        "wavetoy::phi", "Static");
         ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
     $        "wavetoy::pi", "Static");
         ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
     $        "wavetoy::phix", "Static");
         ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
     $        "wavetoy::phiy", "Static");
         ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
     $        "wavetoy::phiz", "Static");
      else if (CCTK_EQUALS(bound, "radiation")) then
         write(0,*) ' FOWavetoy: bound is radiation'
         ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
     $        "wavetoy::phi", "Radiation");
      else if (CCTK_EQUALS(bound, "robin")) then
         write(0,*) ' FOWavetoy: bound is robin'
         ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
     $        "wavetoy::phi", "Robin");
      else if (CCTK_EQUALS(bound, "none")) then
         ierr = 0
      else
         call CCTK_WARN (0, "internal error")
      end if
      if (ierr .lt. 0)  then
        write(0,*) ' ierr = ',ierr
        call CCTK_WARN (0, "Boundary conditions not applied - giving up!")
      end if
      
c  SHH: Taking these out on June 26, 2003.  Presumably symmetry BCs are
c       now handled by some other thorn?
c     Apply the symmetry boundary conditions on any coordinate axes
c     -------------------------------------------------------------
c      call Cart3dSymGN (ierr, cctkGH, "wavetoy::scalarevolve")
c      call Cart3dSymGN (ierr, cctkGH, "wavetoy::scalarevolve_derivs")
c      
c      if (ierr .lt. 0)  then
c        call CCTK_WARN (0, "Symmetry conditions not applied - giving up!")
c      end if
      
      end
