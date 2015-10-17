c     $Header:$

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
      
      
      
      function WaveToyFO_Deriv (fm2, fm1, f0, fp1, fp2, dx, dir)
      implicit none
      CCTK_REAL WaveToyFO_Deriv
      CCTK_REAL fm2, fm1, f0, fp1, fp2
      CCTK_REAL dx
      integer dir
      
      if (dir .eq. 0) then
         WaveToyFO_Deriv = (fp1 - fm1) / (2 * dx)
      else if (dir .eq. -1) then
c     Use first order one-sided derivatives here.
c     These are less dissipative than second order one-sided derivatives.
c     (Sic.)
         WaveToyFO_Deriv = (f0 - fm1) / dx
c$$$         WaveToyFO_Deriv = (+ 3*f0 - 4*fm1 + fm2) / (2*dx)
      else if (dir .eq. +1) then
         WaveToyFO_Deriv = (fp1 - f0) / dx
c$$$         WaveToyFO_Deriv = (- 3*f0 + 4*fp1 - fp2) / (2*dx)
      else
         call CCTK_WARN (0, "internal error")
      end if
      
      end
      
      
      
#define CHOOSEDIR(i,imin,imax) ((min(i+1,imax) - (i+1)) + (max(i-1,imin) - (i-1)))
      
#define DIFFX(var,i,j,k) WaveToyFO_Deriv(var((i)-2,j,k), var((i)-1,j,k), var(i,j,k), var((i)+1,j,k), var((i)+2,j,k), dx(1), idir(1))
#define DIFFY(var,i,j,k) WaveToyFO_Deriv(var(i,(j)-2,k), var(i,(j)-1,k), var(i,j,k), var(i,(j)+1,k), var(i,(j)+2,k), dx(2), idir(2))
#define DIFFZ(var,i,j,k) WaveToyFO_Deriv(var(i,j,(k)-2), var(i,j,(k)-1), var(i,j,k), var(i,j,(k)+1), var(i,j,(k)+2), dx(3), idir(3))
      
      
      
      subroutine WaveToyFO_CalcRHS (CCTK_ARGUMENTS)
      implicit none
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_FUNCTIONS
      DECLARE_CCTK_PARAMETERS
      
      CCTK_REAL dx(3)
      integer bndwidth(3)
      integer imin(3), imax(3), idir(3)
      integer i, j, k
      integer d
      
      external WaveToyFO_Deriv
      CCTK_REAL WaveToyFO_Deriv
      
      do d=1,3
         bndwidth(d) = 0
      end do
      
      do d=1,3
         if (cctk_bbox(2*d).ne.0) then
            imin(d) = 1+bndwidth(d)
         else
            imin(d) = 1+cctk_nghostzones(d)
         end if
         if (cctk_bbox(2*d+1).ne.0) then
            imax(d) = cctk_lsh(d)-bndwidth(d)
         else
            imax(d) = cctk_lsh(d)-cctk_nghostzones(d)
         end if
      end do
      
      do d=1,3
         dx(d) = CCTK_DELTA_SPACE(d)
      end do
      
      do k=imin(3),imax(3)
         idir(3) = CHOOSEDIR(k,imin(3),imax(3))
         do j=imin(2),imax(2)
            idir(2) = CHOOSEDIR(j,imin(2),imax(2))
            do i=imin(1),imax(1)
               idir(1) = CHOOSEDIR(i,imin(1),imax(1))
               
               phidot(i,j,k) =   DIFFX (psix,i,j,k)
     $                         + DIFFY (psiy,i,j,k)
     $                         + DIFFZ (psiz,i,j,k)
               psixdot(i,j,k) = DIFFX (phi,i,j,k)
               psiydot(i,j,k) = DIFFY (phi,i,j,k)
               psizdot(i,j,k) = DIFFZ (phi,i,j,k)
               
            end do
         end do
      end do
      
      end
