c     -*-Fortran-*-
c     $Header:$
      
#include "cctk.h"


c$$$ This routine performs "TVD" prolongation. It is intended to be used 
c$$$ with GFs that are not expected to be smooth, particularly those
c$$$ that must also obey certain constraints. The obvious example is the 
c$$$ density in hydrodynamics, which may be discontinuous yet must be
c$$$ strictly positive.
c$$$
c$$$ To ensure that this prolongation method is used you should add the
c$$$ tag
c$$$
c$$$      tags='Prolongation="TVD"'
c$$$
c$$$ to the interface.ccl on the appropriate group.
c$$$
c$$$ This applies minmod type limiting to the slope, checking over the
c$$$ entire coarse grid cell for the minimum modulus in each direction.
c$$$
c$$$ The actual minmod function is defined in the routine
c$$$
c$$$      prolongate_3d_real8_minmod.F77


#define CHKIDX(i,j,k, imax,jmax,kmax, where)                                  \
      if ((i).lt.1 .or. (i).gt.(imax)                                         \
           .or. (j).lt.1 .or. (j).gt.(jmax)                                   \
           .or. (k).lt.1 .or. (k).gt.(kmax)) then                           &&\
         write (msg, '(a, " array index out of bounds: shape is (",i4,",",i4,",",i4,"), index is (",i4,",",i4,",",i4,")")') \
              (where), (imax), (jmax), (kmax), (i), (j), (k)                &&\
         call CCTK_WARN (0, msg(1:len_trim(msg)))                           &&\
      end if

      function minmod(a, b)

      implicit none

      CCTK_REAL8 minmod
      CCTK_REAL8 a, b
      CCTK_REAL8 zero
      parameter (zero = 0)

      if (a * b .lt. zero) then
        minmod = zero
      else if (abs(a) < abs(b)) then
        minmod = a
      else
        minmod = b
      end if

      end function
      
      subroutine prolongate_3d_real8_minmod (
     $     src, srciext, srcjext, srckext,
     $     dst, dstiext, dstjext, dstkext,
     $     srcbbox, dstbbox, regbbox)
      
      implicit none
      
      CCTK_REAL8 one
      parameter (one = 1)
      
      integer srciext, srcjext, srckext
      CCTK_REAL8 src(srciext,srcjext,srckext)
      integer dstiext, dstjext, dstkext
      CCTK_REAL8 dst(dstiext,dstjext,dstkext)
c     bbox(:,1) is lower boundary (inclusive)
c     bbox(:,2) is upper boundary (inclusive)
c     bbox(:,3) is stride
      integer srcbbox(3,3), dstbbox(3,3), regbbox(3,3)
      
      integer offsetlo, offsethi
      
      integer regiext, regjext, regkext
      
      integer dstifac, dstjfac, dstkfac
      
      integer srcioff, srcjoff, srckoff
      integer dstioff, dstjoff, dstkoff
      
      integer i, j, k
      integer i0, j0, k0
      integer fi, fj, fk
      integer ifac(4), jfac(4), kfac(4)
      integer ii, jj, kk
      integer fac
      CCTK_REAL8 res
      integer d
      
      character msg*1000
      
      external minmod
      CCTK_REAL8 minmod

      CCTK_REAL8 half, zero
      parameter (half = 0.5)
      parameter (zero = 0)
      CCTK_REAL8 dupw, dloc, slopex, slopey, slopez
      logical firstloop
      

      do d=1,3
         if (srcbbox(d,3).eq.0 .or. dstbbox(d,3).eq.0
     $        .or. regbbox(d,3).eq.0) then
            call CCTK_WARN (0, "Internal error: stride is zero")
         end if
         if (srcbbox(d,3).le.regbbox(d,3)
     $        .or. dstbbox(d,3).ne.regbbox(d,3)) then
            call CCTK_WARN (0, "Internal error: strides disagree")
         end if
         if (mod(srcbbox(d,3), dstbbox(d,3)).ne.0) then
            call CCTK_WARN (0, "Internal error: destination strides are not integer multiples of the source strides")
         end if
         if (mod(srcbbox(d,1), srcbbox(d,3)).ne.0
     $        .or. mod(dstbbox(d,1), dstbbox(d,3)).ne.0
     $        .or. mod(regbbox(d,1), regbbox(d,3)).ne.0) then
            call CCTK_WARN (0, "Internal error: array origins are not integer multiples of the strides")
         end if
         if (regbbox(d,1).gt.regbbox(d,2)) then
c     This could be handled, but is likely to point to an error elsewhere
            call CCTK_WARN (0, "Internal error: region extent is empty")
         end if
         if (mod(srcbbox(d,2) - srcbbox(d,1), srcbbox(d,3)).ne.0
     $        .or. mod(dstbbox(d,2) - dstbbox(d,1), dstbbox(d,3)).ne.0
     $        .or. mod(regbbox(d,2) - regbbox(d,1), regbbox(d,3)).ne.0) then
            call CCTK_WARN (0, "Internal error: array extents are not integer multiples of the strides")
         end if
         regkext = (regbbox(d,2) - regbbox(d,1)) / regbbox(d,3) + 1
         dstkfac = srcbbox(d,3) / dstbbox(d,3)
         srckoff = (regbbox(d,1) - srcbbox(d,1)) / dstbbox(d,3)
         offsetlo = regbbox(d,3)
         if (mod(srckoff + 0, dstkfac).eq.0) then
            offsetlo = 0
            if (regkext.gt.1) then
               offsetlo = regbbox(d,3)
            end if
         end if
         offsethi = regbbox(d,3)
         if (mod(srckoff + regkext-1, dstkfac).eq.0) then
            offsethi = 0
            if (regkext.gt.1) then
               offsethi = regbbox(d,3)
            end if
         end if
         if (regbbox(d,1)-offsetlo.lt.srcbbox(d,1)
     $        .or. regbbox(d,2)+offsethi.gt.srcbbox(d,2)
     $        .or. regbbox(d,1).lt.dstbbox(d,1)
     $        .or. regbbox(d,2).gt.dstbbox(d,2)) then
            call CCTK_WARN (0, "Internal error: region extent is not contained in array extent")
         end if
      end do
      
      if (srciext.ne.(srcbbox(1,2)-srcbbox(1,1))/srcbbox(1,3)+1
     $     .or. srcjext.ne.(srcbbox(2,2)-srcbbox(2,1))/srcbbox(2,3)+1
     $     .or. srckext.ne.(srcbbox(3,2)-srcbbox(3,1))/srcbbox(3,3)+1
     $     .or. dstiext.ne.(dstbbox(1,2)-dstbbox(1,1))/dstbbox(1,3)+1
     $     .or. dstjext.ne.(dstbbox(2,2)-dstbbox(2,1))/dstbbox(2,3)+1
     $     .or. dstkext.ne.(dstbbox(3,2)-dstbbox(3,1))/dstbbox(3,3)+1) then
         call CCTK_WARN (0, "Internal error: array sizes don't agree with bounding boxes")
      end if
      
      
      
      regiext = (regbbox(1,2) - regbbox(1,1)) / regbbox(1,3) + 1
      regjext = (regbbox(2,2) - regbbox(2,1)) / regbbox(2,3) + 1
      regkext = (regbbox(3,2) - regbbox(3,1)) / regbbox(3,3) + 1
      
      dstifac = srcbbox(1,3) / dstbbox(1,3)
      dstjfac = srcbbox(2,3) / dstbbox(2,3)
      dstkfac = srcbbox(3,3) / dstbbox(3,3)
      
      srcioff = (regbbox(1,1) - srcbbox(1,1)) / dstbbox(1,3)
      srcjoff = (regbbox(2,1) - srcbbox(2,1)) / dstbbox(2,3)
      srckoff = (regbbox(3,1) - srcbbox(3,1)) / dstbbox(3,3)
      
      dstioff = (regbbox(1,1) - dstbbox(1,1)) / dstbbox(1,3)
      dstjoff = (regbbox(2,1) - dstbbox(2,1)) / dstbbox(2,3)
      dstkoff = (regbbox(3,1) - dstbbox(3,1)) / dstbbox(3,3)
      
      
      
c     Loop over fine region

      do k = 0, regkext-1
         k0 = (srckoff + k) / dstkfac
         fk = mod(srckoff + k, dstkfac)
         
         do j = 0, regjext-1
            j0 = (srcjoff + j) / dstjfac
            fj = mod(srcjoff + j, dstjfac)
            
            do i = 0, regiext-1
               i0 = (srcioff + i) / dstifac
               fi = mod(srcioff + i, dstifac)

               slopex = zero
               slopey = zero
               slopez = zero

               firstloop = .true.

               do kk = 1, 2
                 do jj = 1, 2

                   dupw = src(i0+1 ,j0+jj,k0+kk) - src(i0+0 ,j0+jj,k0+kk)
                   dloc = src(i0+2 ,j0+jj,k0+kk) - src(i0+1 ,j0+kk,k0+kk)
                   if (firstloop) then
                     slopex = half * dble(fi) * minmod(dupw,dloc)
                     firstloop = .false.
                   else
                     slopex = 
     $                    minmod(slopex, half * dble(fi) * minmod(dupw,dloc))
                   end if
                 end do
               end do
               
               firstloop = .true.

               do kk = 1, 2
                 do ii = 1, 2

                   dupw = src(i0+ii,j0+1 ,k0+kk) - src(i0+ii,j0+0 ,k0+kk)
                   dloc = src(i0+ii,j0+2 ,k0+kk) - src(i0+ii,j0+1 ,k0+kk)
                   if (firstloop) then
                     slopey = half * dble(fj) * minmod(dupw,dloc)
                     firstloop = .false.
                   else
                     slopey = 
     $                    minmod(slopey, half * dble(fj) * minmod(dupw,dloc))
                   end if
                 end do
               end do

               firstloop = .true.

               do jj = 1, 2
                 do ii = 1, 2

                   dupw = src(i0+ii,j0+jj,k0+1 ) - src(i0+ii,j0+jj,k0+0 )
                   dloc = src(i0+ii,j0+jj,k0+2 ) - src(i0+ii,j0+jj,k0+1 )
                   if (firstloop) then
                     slopez = half * dble(fk) * minmod(dupw,dloc)
                     firstloop = .false.
                   else
                     slopez = 
     $                    minmod(slopez, half * dble(fk) * minmod(dupw,dloc))
                   end if
                 end do
               end do

               CHKIDX (dstioff+i+1, dstjoff+j+1, dstkoff+k+1, \
                       dstiext,dstjext,dstkext, "destination")
               dst (dstioff+i+1, dstjoff+j+1, dstkoff+k+1) = 
     .                      src(i0+1,j0+1,k0+1) + slopex + slopey + slopez
               

            end do
         end do
      end do
      
      end
