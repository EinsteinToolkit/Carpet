c     -*-Fortran-*-
c     $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/prolongate_3d_real8_3tl.F77,v 1.2 2001/12/09 16:43:11 schnetter Exp $
      
#include "cctk.h"
      
      subroutine prolongate_3d_real8_3tl (
     $     src1, t1, src2, t2, src3, t3, srciext, srcjext, srckext,
     $     dst, t, dstiext, dstjext, dstkext,
     $     srcbbox, dstbbox, regbbox)
      
      implicit none
      
      CCTK_REAL8 one
      parameter (one = 1)
      
      integer srciext, srcjext, srckext
      CCTK_REAL8 src1(srciext,srcjext,srckext)
      integer t1
      CCTK_REAL8 src2(srciext,srcjext,srckext)
      integer t2
      CCTK_REAL8 src3(srciext,srcjext,srckext)
      integer t3
      integer dstiext, dstjext, dstkext
      CCTK_REAL8 dst(dstiext,dstjext,dstkext)
      integer t
c     bbox(:,1) is lower boundary (inclusive)
c     bbox(:,2) is upper boundary (inclusive)
c     bbox(:,3) is stride
      integer srcbbox(3,3), dstbbox(3,3), regbbox(3,3)
      
      integer regiext, regjext, regkext
      
      integer dstifac, dstjfac, dstkfac
      
      integer srcioff, srcjoff, srckoff
      integer dstioff, dstjoff, dstkoff
      
      CCTK_REAL8 s1fac, s2fac, s3fac
      
      integer i, j, k
      integer i0, j0, k0
      integer fi, fj, fk
      integer ii, jj, kk
      integer fac
      CCTK_REAL8 res
      integer d
      
      
      
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
         if (regbbox(d,1).lt.srcbbox(d,1)
     $        .or. regbbox(d,1).lt.dstbbox(d,1)
     $        .or. regbbox(d,2).gt.srcbbox(d,2)
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
      
      
      
c     Quadratic (second order) interpolation
      if (t1.eq.t2 .or. t1.eq.t3 .or. t2.eq.t3) then
         call CCTK_WARN (0, "Internal error: arrays have same time")
      end if
      
      s1fac = (t - t2) * (t - t3) * one / ((t1 - t2) * (t1 - t3))
      s2fac = (t - t1) * (t - t3) * one / ((t2 - t1) * (t2 - t3))
      s3fac = (t - t1) * (t - t2) * one / ((t3 - t1) * (t3 - t2))
      
      
      
c     Loop over fine region
      do k = 0, regkext-1
         do j = 0, regjext-1
            do i = 0, regiext-1
               
               i0 = (srcioff + i) / dstifac + 1
               j0 = (srcjoff + j) / dstjfac + 1
               k0 = (srckoff + k) / dstkfac + 1
               
               fi = mod(srcioff + i, dstifac)
               fj = mod(srcjoff + j, dstjfac)
               fk = mod(srckoff + k, dstkfac)
               
               res = 0
               
               do kk=0,1
                  do jj=0,1
                     do ii=0,1
                        
                        fac = 1
                        
                        if (ii.eq.0) then
                           fac = fac * (dstifac - fi)
                        else
                           fac = fac *            fi
                        end if
                        if (jj.eq.0) then
                           fac = fac * (dstjfac - fj)
                        else
                           fac = fac *            fj
                        end if
                        if (kk.eq.0) then
                           fac = fac * (dstkfac - fk)
                        else
                           fac = fac *            fk
                        end if
                        
                        if (fac.ne.0) then
                           res = res
     $                          + fac * s1fac * src1(i0+ii,j0+jj,k0+kk)
     $                          + fac * s2fac * src2(i0+ii,j0+jj,k0+kk)
     $                          + fac * s3fac * src3(i0+ii,j0+jj,k0+kk)
                        end if
                        
                     end do
                  end do
               end do
               
               dst (dstioff+i+1, dstjoff+j+1, dstkoff+k+1)
     $              = res / (dstifac * dstjfac * dstkfac)
               
            end do
         end do
      end do
      
      end
