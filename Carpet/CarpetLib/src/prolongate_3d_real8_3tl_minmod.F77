c     -*-Fortran-*-
c     $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/prolongate_3d_real8_3tl_minmod.F77,v 1.2 2003/10/15 11:51:55 hawke Exp $
      
#include "cctk.h"



#define CHKIDX(i,j,k, imax,jmax,kmax, where)                                  \
      if ((i).lt.1 .or. (i).gt.(imax)                                         \
           .or. (j).lt.1 .or. (j).gt.(jmax)                                   \
           .or. (k).lt.1 .or. (k).gt.(kmax)) then                           &&\
         write (msg, '(a, " array index out of bounds: shape is (",i4,",",i4,",",i4,"), index is (",i4,",",i4,",",i4,")")') \
              (where), (imax), (jmax), (kmax), (i), (j), (k)                &&\
         call CCTK_WARN (0, msg(1:len_trim(msg)))                           &&\
      end if


      
      subroutine prolongate_3d_real8_3tl_minmod (
     $     src1, t1, src2, t2, src3, t3, srciext, srcjext, srckext,
     $     dst, t, dstiext, dstjext, dstkext,
     $     srcbbox, dstbbox, regbbox)
      
      implicit none
      
      CCTK_REAL8 one
      parameter (one = 1)
      
      CCTK_REAL8 eps
      parameter (eps = 1.0d-10)
      
      integer srciext, srcjext, srckext
      CCTK_REAL8 src1(srciext,srcjext,srckext)
      CCTK_REAL8 t1
      CCTK_REAL8 src2(srciext,srcjext,srckext)
      CCTK_REAL8 t2
      CCTK_REAL8 src3(srciext,srcjext,srckext)
      CCTK_REAL8 t3
      integer dstiext, dstjext, dstkext
      CCTK_REAL8 dst(dstiext,dstjext,dstkext)
      CCTK_REAL8 t
c     bbox(:,1) is lower boundary (inclusive)
c     bbox(:,2) is upper boundary (inclusive)
c     bbox(:,3) is stride
      integer srcbbox(3,3), dstbbox(3,3), regbbox(3,3)
      
      integer offsetlo, offsethi
      
      integer regiext, regjext, regkext
      
      integer dstifac, dstjfac, dstkfac
      
      integer srcioff, srcjoff, srckoff
      integer dstioff, dstjoff, dstkoff
      
      CCTK_REAL8 s1fac, s2fac, s3fac
      
      CCTK_REAL8 dstdiv
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
      CCTK_REAL8 dupw, dloc, slopex(3), slopey(3), slopez(3)
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
      
      
      
c     Quadratic (second order) interpolation
      if (t1.eq.t2 .or. t1.eq.t3 .or. t2.eq.t3) then
         call CCTK_WARN (0, "Internal error: arrays have same time")
      end if
      if (t.lt.min(t1,t2,t3)-eps .or. t.gt.max(t1,t2,t3)+eps) then
         call CCTK_WARN (0, "Internal error: extrapolation")
      end if
      
      s1fac = (t - t2) * (t - t3) / ((t1 - t2) * (t1 - t3))
      s2fac = (t - t1) * (t - t3) / ((t2 - t1) * (t2 - t3))
      s3fac = (t - t1) * (t - t2) / ((t3 - t1) * (t3 - t2))
      
      
      
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


               slopex(1) = zero
               slopey(1) = zero
               slopez(1) = zero

c$$$               firstloop = .true.
c$$$
c$$$               do kk = 1, 3, 2
c$$$                 do jj = 1, 3, 2
c$$$
c$$$                   dupw = src1(i0+1 ,j0+jj,k0+kk) - src1(i0+0 ,j0+jj,k0+kk)
c$$$                   dloc = src1(i0+2 ,j0+jj,k0+kk) - src1(i0+1 ,j0+kk,k0+kk)
c$$$                   if (firstloop) then
c$$$                     slopex(1) = half * dble(fi) * minmod(dupw,dloc)
c$$$                     firstloop = .false.
c$$$                   else
c$$$                     slopex(1) = 
c$$$     $                    minmod(slopex(1), half * dble(fi) * minmod(dupw,dloc))
c$$$                   end if
c$$$                 end do
c$$$               end do

               dupw = src1(i0+1,j0+1,k0+1) - src1(i0  ,j0+1,k0+1)
               dloc = src1(i0+2,j0+1,k0+1) - src1(i0+1,j0+1,k0+1)
               slopex(1) = half * dble(fi) * minmod(dupw,dloc)

c$$$               firstloop = .true.
c$$$
c$$$               do kk = 1, 3, 2
c$$$                 do ii = 1, 3, 2
c$$$
c$$$                   dupw = src1(i0+ii,j0+1 ,k0+kk) - src1(i0+ii,j0+0 ,k0+kk)
c$$$                   dloc = src1(i0+ii,j0+2 ,k0+kk) - src1(i0+ii,j0+1 ,k0+kk)
c$$$                   if (firstloop) then
c$$$                     slopey(1) = half * dble(fj) * minmod(dupw,dloc)
c$$$                     firstloop = .false.
c$$$                   else
c$$$                     slopey(1) = 
c$$$     $                    minmod(slopey(1), half * dble(fj) * minmod(dupw,dloc))
c$$$                   end if
c$$$                 end do
c$$$               end do

               dupw = src1(i0+1,j0+1,k0+1) - src1(i0+1,j0  ,k0+1)
               dloc = src1(i0+1,j0+2,k0+1) - src1(i0+1,j0+1,k0+1)
               slopey(1) = half * dble(fj) * minmod(dupw,dloc)
c$$$
c$$$               firstloop = .true.
c$$$
c$$$               do jj = 1, 3, 2
c$$$                 do ii = 1, 3, 2
c$$$
c$$$                   dupw = src1(i0+ii,j0+jj,k0+1 ) - src1(i0+ii,j0+jj,k0+0 )
c$$$                   dloc = src1(i0+ii,j0+jj,k0+2 ) - src1(i0+ii,j0+jj,k0+1 )
c$$$                   if (firstloop) then
c$$$                     slopez(1) = half * dble(fk) * minmod(dupw,dloc)
c$$$                     firstloop = .false.
c$$$                   else
c$$$                     slopez(1) = 
c$$$     $                    minmod(slopez(1), half * dble(fk) * minmod(dupw,dloc))
c$$$                   end if
c$$$
c$$$                 end do
c$$$               end do

               dupw = src1(i0+1,j0+1,k0+1) - src1(i0+1,j0+1,k0  )
               dloc = src1(i0+1,j0+1,k0+2) - src1(i0+1,j0+1,k0+1)
               slopez(1) = half * dble(fk) * minmod(dupw,dloc)

               slopex(2) = zero
               slopey(2) = zero
               slopez(2) = zero

c$$$               firstloop = .true.
c$$$
c$$$               do kk = 1, 3, 2
c$$$                 do jj = 1, 3, 2
c$$$
c$$$                   dupw = src2(i0+1 ,j0+jj,k0+kk) - src2(i0+0 ,j0+jj,k0+kk)
c$$$                   dloc = src2(i0+2 ,j0+jj,k0+kk) - src2(i0+1 ,j0+kk,k0+kk)
c$$$                   if (firstloop) then
c$$$                     slopex(2) = half * dble(fi) * minmod(dupw,dloc)
c$$$                     firstloop = .false.
c$$$                   else
c$$$                     slopex(2) = 
c$$$     $                    minmod(slopex(2), half * dble(fi) * minmod(dupw,dloc))
c$$$                   end if
c$$$                 end do
c$$$               end do

               dupw = src2(i0+1,j0+1,k0+1) - src2(i0  ,j0+1,k0+1)
               dloc = src2(i0+2,j0+1,k0+1) - src2(i0+1,j0+1,k0+1)
               slopex(2) = half * dble(fi) * minmod(dupw,dloc)
c$$$
c$$$               firstloop = .true.
c$$$
c$$$               do kk = 1, 3, 2
c$$$                 do ii = 1, 3, 2
c$$$
c$$$                   dupw = src2(i0+ii,j0+1 ,k0+kk) - src2(i0+ii,j0+0 ,k0+kk)
c$$$                   dloc = src2(i0+ii,j0+2 ,k0+kk) - src2(i0+ii,j0+1 ,k0+kk)
c$$$                   if (firstloop) then
c$$$                   slopey(2) = half * dble(fj) * minmod(dupw,dloc)
c$$$                     firstloop = .false.
c$$$                   else
c$$$                     slopey(2) = 
c$$$     $                    minmod(slopey(2), half * dble(fj) * minmod(dupw,dloc))
c$$$                   end if
c$$$                 end do
c$$$               end do

               dupw = src2(i0+1,j0+1,k0+1) - src2(i0+1,j0  ,k0+1)
               dloc = src2(i0+1,j0+2,k0+1) - src2(i0+1,j0+1,k0+1)
               slopey(2) = half * dble(fj) * minmod(dupw,dloc)

c$$$               firstloop = .true.
c$$$
c$$$               do jj = 1, 3, 2
c$$$                 do ii = 1, 3, 2
c$$$
c$$$                   dupw = src2(i0+ii,j0+jj,k0+1 ) - src2(i0+ii,j0+jj,k0+0 )
c$$$                   dloc = src2(i0+ii,j0+jj,k0+2 ) - src2(i0+ii,j0+jj,k0+1 )
c$$$                   if (firstloop) then
c$$$                     slopez(2) = half * dble(fk) * minmod(dupw,dloc)
c$$$                     firstloop = .false.
c$$$                   else
c$$$                     slopez(2) = 
c$$$     $                    minmod(slopez(2), half * dble(fk) * minmod(dupw,dloc))
c$$$                   end if
c$$$                 end do
c$$$               end do

               dupw = src2(i0+1,j0+1,k0+1) - src2(i0+1,j0+1,k0  )
               dloc = src2(i0+1,j0+1,k0+2) - src2(i0+1,j0+1,k0+1)
               slopez(2) = half * dble(fk) * minmod(dupw,dloc)

c$$$               firstloop = .true.

               slopex(3) = zero
               slopey(3) = zero
               slopez(3) = zero

c$$$               do kk = 1, 3, 2
c$$$                 do jj = 1, 3, 2
c$$$
c$$$                   dupw = src3(i0+1 ,j0+jj,k0+kk) - src3(i0+0 ,j0+jj,k0+kk)
c$$$                   dloc = src3(i0+2 ,j0+jj,k0+kk) - src3(i0+1 ,j0+kk,k0+kk)
c$$$                   if (firstloop) then
c$$$                     slopex(3) = half * dble(fi) * minmod(dupw,dloc)
c$$$                     firstloop = .false.
c$$$                   else
c$$$                     slopex(3) = 
c$$$     $                    minmod(slopex(3), half * dble(fi) * minmod(dupw,dloc))
c$$$                   end if
c$$$                 end do
c$$$               end do

               dupw = src3(i0+1,j0+1,k0+1) - src3(i0  ,j0+1,k0+1)
               dloc = src3(i0+2,j0+1,k0+1) - src3(i0+1,j0+1,k0+1)
               slopex(3) = half * dble(fi) * minmod(dupw,dloc)

c$$$               firstloop = .true.
c$$$
c$$$               do kk = 1, 3, 2
c$$$                 do ii = 1, 3, 2
c$$$
c$$$                   dupw = src3(i0+ii,j0+1 ,k0+kk) - src3(i0+ii,j0+0 ,k0+kk)
c$$$                   dloc = src3(i0+ii,j0+2 ,k0+kk) - src3(i0+ii,j0+1 ,k0+kk)
c$$$                   if (firstloop) then
c$$$                     slopey(3) = half * dble(fj) * minmod(dupw,dloc)
c$$$                     firstloop = .false.
c$$$                   else
c$$$                     slopey(3) = 
c$$$     $                    minmod(slopey(3), half * dble(fj) * minmod(dupw,dloc))
c$$$                   end if
c$$$                 end do
c$$$               end do

               dupw = src3(i0+1,j0+1,k0+1) - src3(i0+1,j0  ,k0+1)
               dloc = src3(i0+1,j0+2,k0+1) - src3(i0+1,j0+1,k0+1)
               slopey(3) = half * dble(fj) * minmod(dupw,dloc)
c$$$
c$$$               firstloop = .true.
c$$$
c$$$               do jj = 1, 3, 2
c$$$                 do ii = 1, 3, 2
c$$$
c$$$                   dupw = src3(i0+ii,j0+jj,k0+1 ) - src3(i0+ii,j0+jj,k0+0 )
c$$$                   dloc = src3(i0+ii,j0+jj,k0+2 ) - src3(i0+ii,j0+jj,k0+1 )
c$$$                   if (firstloop) then
c$$$                     slopez(3) = half * dble(fk) * minmod(dupw,dloc)
c$$$                     firstloop = .false.
c$$$                   else
c$$$                     slopez(3) = 
c$$$     $                    minmod(slopez(3), half * dble(fk) * minmod(dupw,dloc))
c$$$                   end if
c$$$                 end do
c$$$               end do

               dupw = src3(i0+1,j0+1,k0+1) - src3(i0+1,j0+1,k0  )
               dloc = src3(i0+1,j0+1,k0+2) - src3(i0+1,j0+1,k0+1)
               slopez(3) = half * dble(fk) * minmod(dupw,dloc)
                              
               CHKIDX (dstioff+i+1, dstjoff+j+1, dstkoff+k+1, \
                       dstiext,dstjext,dstkext, "destination")
               dst (dstioff+i+1, dstjoff+j+1, dstkoff+k+1) = 
     $                      s1fac * (src1(i0+1,j0+1,k0+1) + 
     $                      slopex(1) + slopey(1) + slopez(1)) +
     $                      s2fac * (src2(i0+1,j0+1,k0+1) + 
     $                      slopex(2) + slopey(2) + slopez(2)) +
     $                      s3fac * (src3(i0+1,j0+1,k0+1) + 
     $                      slopex(3) + slopey(3) + slopez(3))
               
            end do
         end do
      end do
      
      end
