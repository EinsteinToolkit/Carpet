c     -*-Fortran-*-
c     $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/prolongate_3d_real8.F77,v 1.9 2003/02/24 17:43:10 schnetter Exp $
      
#include "cctk.h"
      
      subroutine prolongate_3d_real8 (
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
      
      integer regiext, regjext, regkext
      
      integer dstifac, dstjfac, dstkfac
      
      integer srcioff, srcjoff, srckoff
      integer dstioff, dstjoff, dstkoff
      
      CCTK_REAL8 dstdiv
      integer i, j, k
      integer i0, j0, k0
      integer fi, fj, fk
      integer ifac(2), jfac(2), kfac(2)
c$$$      integer ii, jj, kk
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
         if (regbbox(d,1).gt.regbbox(d,2)) then
c     This could be handled, but is likely to point to an error elsewhere
            call CCTK_WARN (0, "Internal error: region extent is empty")
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
      dstdiv = one / (dstifac * dstjfac * dstkfac)
      
      do k = 0, regkext-1
         k0 = (srckoff + k) / dstkfac
         fk = mod(srckoff + k, dstkfac)
         kfac(1) = (fk-dstkfac) * (-1)
         kfac(2) = (fk        ) * 1
         
         do j = 0, regjext-1
            j0 = (srcjoff + j) / dstjfac
            fj = mod(srcjoff + j, dstjfac)
            jfac(1) = (fj-dstjfac) * (-1)
            jfac(2) = (fj        ) * 1
            
            do i = 0, regiext-1
               i0 = (srcioff + i) / dstifac
               fi = mod(srcioff + i, dstifac)
               ifac(1) = (fi-dstifac) * (-1)
               ifac(2) = (fi        ) * 1
               
               res = 0
               
c$$$               do kk=1,2
c$$$                  do jj=1,2
c$$$                     do ii=1,2
c$$$                        
c$$$                        fac = ifac(ii) * jfac(jj) * kfac(kk)
c$$$                        
c$$$                        if (fac.ne.0) then
c$$$                           res = res + fac * src(i0+ii, j0+jj, k0+kk)
c$$$                        end if
c$$$                        
c$$$                     end do
c$$$                  end do
c$$$               end do
               
               fac = ifac(1) * jfac(1) * kfac(1)
               if (fac.ne.0) res = res + fac * src(i0+1, j0+1, k0+1)
               
               fac = ifac(2) * jfac(1) * kfac(1)
               if (fac.ne.0) res = res + fac * src(i0+2, j0+1, k0+1)
               
               fac = ifac(1) * jfac(2) * kfac(1)
               if (fac.ne.0) res = res + fac * src(i0+1, j0+2, k0+1)
               
               fac = ifac(2) * jfac(2) * kfac(1)
               if (fac.ne.0) res = res + fac * src(i0+2, j0+2, k0+1)
               
               fac = ifac(1) * jfac(1) * kfac(2)
               if (fac.ne.0) res = res + fac * src(i0+1, j0+1, k0+2)
               
               fac = ifac(2) * jfac(1) * kfac(2)
               if (fac.ne.0) res = res + fac * src(i0+2, j0+1, k0+2)
               
               fac = ifac(1) * jfac(2) * kfac(2)
               if (fac.ne.0) res = res + fac * src(i0+1, j0+2, k0+2)
               
               fac = ifac(2) * jfac(2) * kfac(2)
               if (fac.ne.0) res = res + fac * src(i0+2, j0+2, k0+2)
               
               dst (dstioff+i+1, dstjoff+j+1, dstkoff+k+1) = dstdiv * res
               
            end do
         end do
      end do
      
      end
