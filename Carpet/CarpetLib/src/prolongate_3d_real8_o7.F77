c     -*-Fortran-*-
      
#include "cctk.h"



#define CHKIDX(i,j,k, imax,jmax,kmax, where)                                  \
      if ((i).lt.1 .or. (i).gt.(imax)                                         \
           .or. (j).lt.1 .or. (j).gt.(jmax)                                   \
           .or. (k).lt.1 .or. (k).gt.(kmax)) then                           &&\
         write (msg, '(a, " array index out of bounds: shape is (",i4,",",i4,",",i4,"), index is (",i4,",",i4,",",i4,")")') \
              (where), (imax), (jmax), (kmax), (i), (j), (k)                &&\
         call CCTK_WARN (0, msg)                                            &&\
      end if


      
      subroutine prolongate_3d_real8_o7 (
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
      
      CCTK_REAL8 dstdiv
      integer i, j, k
      integer i0, j0, k0
      integer fi, fj, fk
      integer ifac(8), jfac(8), kfac(8)
      integer ii, jj, kk
      CCTK_REAL8 fac
      CCTK_REAL8 res
      integer d
      
      character msg*1000
      
      
      
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
c     (This expression cannot be evaluated as integer)
      dstdiv = one / (5040*dstifac**7) / (5040*dstjfac**7) / (5040*dstkfac**7)
      
      do k = 0, regkext-1
         k0 = (srckoff + k) / dstkfac
         fk = mod(srckoff + k, dstkfac)
         kfac(1) = (fk+2*dstkfac) * (fk+  dstkfac) * (fk        ) * (fk-dstkfac) * (fk-2*dstkfac) * (fk-3*dstkfac) * (fk-4*dstkfac) * (- 1)
         kfac(2) = (fk+3*dstkfac) * (fk+  dstkfac) * (fk        ) * (fk-dstkfac) * (fk-2*dstkfac) * (fk-3*dstkfac) * (fk-4*dstkfac) * (  7)
         kfac(3) = (fk+3*dstkfac) * (fk+2*dstkfac) * (fk        ) * (fk-dstkfac) * (fk-2*dstkfac) * (fk-3*dstkfac) * (fk-4*dstkfac) * (-21)
         kfac(4) = (fk+3*dstkfac) * (fk+2*dstkfac) * (fk+dstkfac) * (fk-dstkfac) * (fk-2*dstkfac) * (fk-3*dstkfac) * (fk-4*dstkfac) * ( 35)
         kfac(5) = (fk+3*dstkfac) * (fk+2*dstkfac) * (fk+dstkfac) * (fk        ) * (fk-2*dstkfac) * (fk-3*dstkfac) * (fk-4*dstkfac) * (-35)
         kfac(6) = (fk+3*dstkfac) * (fk+2*dstkfac) * (fk+dstkfac) * (fk        ) * (fk-  dstkfac) * (fk-3*dstkfac) * (fk-4*dstkfac) * ( 21)
         kfac(7) = (fk+3*dstkfac) * (fk+2*dstkfac) * (fk+dstkfac) * (fk        ) * (fk-  dstkfac) * (fk-2*dstkfac) * (fk-4*dstkfac) * (- 7)
         kfac(8) = (fk+3*dstkfac) * (fk+2*dstkfac) * (fk+dstkfac) * (fk        ) * (fk-  dstkfac) * (fk-2*dstkfac) * (fk-3*dstkfac) * (  1)
         
         do j = 0, regjext-1
            j0 = (srcjoff + j) / dstjfac
            fj = mod(srcjoff + j, dstjfac)
            jfac(1) = (fj+2*dstjfac) * (fj+  dstjfac) * (fj        ) * (fj-dstjfac) * (fj-2*dstjfac) * (fj-3*dstjfac) * (fj-4*dstjfac) * (- 1)
            jfac(2) = (fj+3*dstjfac) * (fj+  dstjfac) * (fj        ) * (fj-dstjfac) * (fj-2*dstjfac) * (fj-3*dstjfac) * (fj-4*dstjfac) * (  7)
            jfac(3) = (fj+3*dstjfac) * (fj+2*dstjfac) * (fj        ) * (fj-dstjfac) * (fj-2*dstjfac) * (fj-3*dstjfac) * (fj-4*dstjfac) * (-21)
            jfac(4) = (fj+3*dstjfac) * (fj+2*dstjfac) * (fj+dstjfac) * (fj-dstjfac) * (fj-2*dstjfac) * (fj-3*dstjfac) * (fj-4*dstjfac) * ( 35)
            jfac(5) = (fj+3*dstjfac) * (fj+2*dstjfac) * (fj+dstjfac) * (fj        ) * (fj-2*dstjfac) * (fj-3*dstjfac) * (fj-4*dstjfac) * (-35)
            jfac(6) = (fj+3*dstjfac) * (fj+2*dstjfac) * (fj+dstjfac) * (fj        ) * (fj-  dstjfac) * (fj-3*dstjfac) * (fj-4*dstjfac) * ( 21)
            jfac(7) = (fj+3*dstjfac) * (fj+2*dstjfac) * (fj+dstjfac) * (fj        ) * (fj-  dstjfac) * (fj-2*dstjfac) * (fj-4*dstjfac) * (- 7)
            jfac(8) = (fj+3*dstjfac) * (fj+2*dstjfac) * (fj+dstjfac) * (fj        ) * (fj-  dstjfac) * (fj-2*dstjfac) * (fj-3*dstjfac) * (  1)
            
            do i = 0, regiext-1
               i0 = (srcioff + i) / dstifac
               fi = mod(srcioff + i, dstifac)
               ifac(1) = (fi+2*dstifac) * (fi+  dstifac) * (fi        ) * (fi-dstifac) * (fi-2*dstifac) * (fi-3*dstifac) * (fi-4*dstifac) * (- 1)
               ifac(2) = (fi+3*dstifac) * (fi+  dstifac) * (fi        ) * (fi-dstifac) * (fi-2*dstifac) * (fi-3*dstifac) * (fi-4*dstifac) * (  7)
               ifac(3) = (fi+3*dstifac) * (fi+2*dstifac) * (fi        ) * (fi-dstifac) * (fi-2*dstifac) * (fi-3*dstifac) * (fi-4*dstifac) * (-21)
               ifac(4) = (fi+3*dstifac) * (fi+2*dstifac) * (fi+dstifac) * (fi-dstifac) * (fi-2*dstifac) * (fi-3*dstifac) * (fi-4*dstifac) * ( 35)
               ifac(5) = (fi+3*dstifac) * (fi+2*dstifac) * (fi+dstifac) * (fi        ) * (fi-2*dstifac) * (fi-3*dstifac) * (fi-4*dstifac) * (-35)
               ifac(6) = (fi+3*dstifac) * (fi+2*dstifac) * (fi+dstifac) * (fi        ) * (fi-  dstifac) * (fi-3*dstifac) * (fi-4*dstifac) * ( 21)
               ifac(7) = (fi+3*dstifac) * (fi+2*dstifac) * (fi+dstifac) * (fi        ) * (fi-  dstifac) * (fi-2*dstifac) * (fi-4*dstifac) * (- 7)
               ifac(8) = (fi+3*dstifac) * (fi+2*dstifac) * (fi+dstifac) * (fi        ) * (fi-  dstifac) * (fi-2*dstifac) * (fi-3*dstifac) * (  1)
               
               res = 0
               
               do kk=1,8
                  do jj=1,8
                     do ii=1,8
                        
                        if (ifac(ii).ne.0 .and. jfac(jj).ne.0 .and. kfac(kk).ne.0) then
c     (This expression cannot be evaluated as integer)
                           fac = one * ifac(ii) * jfac(jj) * kfac(kk)
                           
                           CHKIDX (i0+ii-3, j0+jj-3, k0+kk-3, \
                                   srciext,srcjext,srckext, "source")
                           res = res + fac * src(i0+ii-3, j0+jj-3, k0+kk-3)
                        end if
                        
                     end do
                  end do
               end do
               
               CHKIDX (dstioff+i+1, dstjoff+j+1, dstkoff+k+1, \
                       dstiext,dstjext,dstkext, "destination")
               dst (dstioff+i+1, dstjoff+j+1, dstkoff+k+1) = dstdiv * res
               
            end do
         end do
      end do
      
      end
