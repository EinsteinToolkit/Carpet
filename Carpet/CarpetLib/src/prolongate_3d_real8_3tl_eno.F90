!!$     $Header:$

#include "cctk.h"


!!$ This routine performs "ENO" prolongation. It is intended to be used 
!!$ with GFs that are not expected to be smooth, particularly those
!!$ that must also obey certain constraints. The obvious example is the 
!!$ density in hydrodynamics, which may be discontinuous yet must be
!!$ strictly positive.
!!$
!!$ To ensure that this prolongation method is used you should add the
!!$ tag
!!$
!!$      tags='Prolongation="ENO"'
!!$
!!$ to the interface.ccl on the appropriate group.
!!$
!!$ This applies ENO2 type limiting to the slope, checking over the
!!$ entire coarse grid cell for the least oscillatory quadratic in each 
!!$ direction. If the slope changes sign over the extrema, linear
!!$ interpolation is used instead.
!!$
!!$ The actual eno1d function is defined in the routine
!!$
!!$      prolongate_3d_real8_eno.F77


#define CHKIDX(i,j,k, imax,jmax,kmax, where)                                  \
if ((i).lt.1 .or. (i).gt.(imax)                                         \
  .or. (j).lt.1 .or. (j).gt.(jmax)                                   \
  .or. (k).lt.1 .or. (k).gt.(kmax)) then                           &&\
  write (msg, '(a, " array index out of bounds: shape is (",i4,",",i4,",",i4,"), index is (",i4,",",i4,",",i4,")")') \
  (where), (imax), (jmax), (kmax), (i), (j), (k)                &&\
  call CCTK_WARN (0, msg(1:len_trim(msg)))                           &&\
end if

subroutine prolongate_3d_real8_3tl_eno (src1, t1, src2, t2, &
     src3, t3, srciext, srcjext, srckext, dst, t, dstiext, &
     dstjext, dstkext, srcbbox, dstbbox, regbbox)
      
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
!!$     bbox(:,1) is lower boundary (inclusive)
!!$     bbox(:,2) is upper boundary (inclusive)
!!$     bbox(:,3) is stride
  integer srcbbox(3,3), dstbbox(3,3), regbbox(3,3)

  integer offsetlo, offsethi

  integer regiext, regjext, regkext

  integer dstifac, dstjfac, dstkfac

  integer srcioff, srcjoff, srckoff
  integer dstioff, dstjoff, dstkoff
      
  CCTK_REAL8 s1fac, s2fac, s3fac, tmps1fac, tmps2fac, tmps3fac

  integer i, j, k
  integer i0, j0, k0
  integer fi, fj, fk
  integer ifac(4), jfac(4), kfac(4)
  integer ii, jj, kk
  integer fac
  CCTK_REAL8 res
  integer d

  character msg*1000

  CCTK_REAL8, dimension(0:3,0:3) :: tmp1
  CCTK_REAL8, dimension(0:3) :: tmp2
  CCTK_REAL8 :: dsttmp1, dsttmp2, dsttmp3

  external eno1d
  CCTK_REAL8 eno1d

  CCTK_REAL8 half, zero
  parameter (half = 0.5)
  parameter (zero = 0)

  do d=1,3
    if (srcbbox(d,3).eq.0 .or. dstbbox(d,3).eq.0 &
         .or. regbbox(d,3).eq.0) then
      call CCTK_WARN (0, "Internal error: stride is zero")
    end if
    if (srcbbox(d,3).le.regbbox(d,3) &
         .or. dstbbox(d,3).ne.regbbox(d,3)) then
      call CCTK_WARN (0, "Internal error: strides disagree")
    end if
    if (mod(srcbbox(d,3), dstbbox(d,3)).ne.0) then
      call CCTK_WARN (0, "Internal error: destination strides are not integer multiples of the source strides")
    end if
    if (mod(srcbbox(d,1), srcbbox(d,3)).ne.0 &
         .or. mod(dstbbox(d,1), dstbbox(d,3)).ne.0 &
         .or. mod(regbbox(d,1), regbbox(d,3)).ne.0) then
      call CCTK_WARN (0, "Internal error: array origins are not integer multiples of the strides")
    end if
    if (regbbox(d,1).gt.regbbox(d,2)) then
!!$     This could be handled, but is likely to point to an error elsewhere
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
    if (regbbox(d,1)-offsetlo.lt.srcbbox(d,1) &
         .or. regbbox(d,2)+offsethi.gt.srcbbox(d,2) &
         .or. regbbox(d,1).lt.dstbbox(d,1) &
         .or. regbbox(d,2).gt.dstbbox(d,2)) then
      call CCTK_WARN (0, "Internal error: region extent is not contained in array extent")
    end if
  end do

  if (srciext.ne.(srcbbox(1,2)-srcbbox(1,1))/srcbbox(1,3)+1 &
       .or. srcjext.ne.(srcbbox(2,2)-srcbbox(2,1))/srcbbox(2,3)+1 &
       .or. srckext.ne.(srcbbox(3,2)-srcbbox(3,1))/srcbbox(3,3)+1 &
       .or. dstiext.ne.(dstbbox(1,2)-dstbbox(1,1))/dstbbox(1,3)+1 &
       .or. dstjext.ne.(dstbbox(2,2)-dstbbox(2,1))/dstbbox(2,3)+1 &
       .or. dstkext.ne.(dstbbox(3,2)-dstbbox(3,1))/dstbbox(3,3)+1) then
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
      
!!$     Quadratic (second order) interpolation
  if (t1.eq.t2 .or. t1.eq.t3 .or. t2.eq.t3) then
    call CCTK_WARN (0, "Internal error: arrays have same time")
  end if
  if (t.lt.min(t1,t2,t3)-eps .or. t.gt.max(t1,t2,t3)+eps) then
    call CCTK_WARN (0, "Internal error: extrapolation in time")
  end if
  
  s1fac = (t - t2) * (t - t3) / ((t1 - t2) * (t1 - t3))
  s2fac = (t - t1) * (t - t3) / ((t2 - t1) * (t2 - t3))
  s3fac = (t - t1) * (t - t2) / ((t3 - t1) * (t3 - t2))
      
!!$     Loop over fine region

  do k = 0, regkext-1
    k0 = (srckoff + k) / dstkfac
    fk = mod(srckoff + k, dstkfac)

    do j = 0, regjext-1
      j0 = (srcjoff + j) / dstjfac
      fj = mod(srcjoff + j, dstjfac)

      do i = 0, regiext-1
        i0 = (srcioff + i) / dstifac
        fi = mod(srcioff + i, dstifac)

!!$        Where is the fine grid point w.r.t the coarse grid?

!!$        write(*,*) i,j,k,fi,fj,fk

        select case (fi + 10*fj + 100*fk)
        case (0)
!!$            On a coarse grid point exactly!

          dsttmp1 = src1(i0+1,j0+1,k0+1)
          dsttmp2 = src2(i0+1,j0+1,k0+1)
          dsttmp3 = src3(i0+1,j0+1,k0+1)

        case (1)
!!$          Interpolate only in x

          dsttmp1 = eno1d(src1(i0:i0+3,j0+1,k0+1))
          dsttmp2 = eno1d(src2(i0:i0+3,j0+1,k0+1))
          dsttmp3 = eno1d(src3(i0:i0+3,j0+1,k0+1))

        case (10)
!!$          Interpolate only in y

          dsttmp1 = eno1d(src1(i0+1,j0:j0+3,k0+1))
          dsttmp2 = eno1d(src2(i0+1,j0:j0+3,k0+1))
          dsttmp3 = eno1d(src3(i0+1,j0:j0+3,k0+1))

        case (11)
!!$          Interpolate only in x and y

          do jj = 0, 3
            tmp2(jj) = eno1d(src1(i0:i0+3,j0+jj,k0+1))
          end do

          dsttmp1 = eno1d(tmp2(0:3))

          do jj = 0, 3
            tmp2(jj) = eno1d(src2(i0:i0+3,j0+jj,k0+1))
          end do

          dsttmp2 = eno1d(tmp2(0:3))

          do jj = 0, 3
            tmp2(jj) = eno1d(src3(i0:i0+3,j0+jj,k0+1))
          end do

          dsttmp3 = eno1d(tmp2(0:3))

        case (100)
!!$          Interpolate only in z

          dsttmp1 = eno1d(src1(i0+1,j0+1,k0:k0+3))
          dsttmp2 = eno1d(src2(i0+1,j0+1,k0:k0+3))
          dsttmp3 = eno1d(src3(i0+1,j0+1,k0:k0+3))

        case (101)
!!$          Interpolate only in x and z

          do kk = 0, 3
            tmp2(kk) = eno1d(src1(i0:i0+3,j0+1,k0+kk))
          end do

          dsttmp1 = eno1d(tmp2(0:3))

          do kk = 0, 3
            tmp2(kk) = eno1d(src2(i0:i0+3,j0+1,k0+kk))
          end do

          dsttmp2 = eno1d(tmp2(0:3))

          do kk = 0, 3
            tmp2(kk) = eno1d(src3(i0:i0+3,j0+1,k0+kk))
          end do

          dsttmp3 = eno1d(tmp2(0:3))

        case (110)
!!$          Interpolate only in y and z

          do kk = 0, 3
            tmp2(kk) = eno1d(src1(i0+1,j0:j0+3,k0+kk))
          end do

          dsttmp1 = eno1d(tmp2(0:3))

          do kk = 0, 3
            tmp2(kk) = eno1d(src2(i0+1,j0:j0+3,k0+kk))
          end do

          dsttmp2 = eno1d(tmp2(0:3))

          do kk = 0, 3
            tmp2(kk) = eno1d(src3(i0+1,j0:j0+3,k0+kk))
          end do

          dsttmp3 = eno1d(tmp2(0:3))

        case (111)
!!$          Interpolate in all of x, y, and z

          do jj = 0, 3
            do kk = 0, 3
              tmp1(jj,kk) = eno1d(src1(i0:i0+3,j0+jj,k0+kk))
            end do
          end do
          do ii = 0, 3
            tmp2(ii) = eno1d(tmp1(0:3,ii))
          end do
          
          dsttmp1 = eno1d(tmp2(0:3))

          do jj = 0, 3
            do kk = 0, 3
              tmp1(jj,kk) = eno1d(src2(i0:i0+3,j0+jj,k0+kk))
            end do
          end do
          do ii = 0, 3
            tmp2(ii) = eno1d(tmp1(0:3,ii))
          end do
          
          dsttmp2 = eno1d(tmp2(0:3))

          do jj = 0, 3
            do kk = 0, 3
              tmp1(jj,kk) = eno1d(src3(i0:i0+3,j0+jj,k0+kk))
            end do
          end do
          do ii = 0, 3
            tmp2(ii) = eno1d(tmp1(0:3,ii))
          end do
          
          dsttmp3 = eno1d(tmp2(0:3))

        case default
          call CCTK_WARN(0, "Internal error in ENO prolongation. Should only be used with refinement factor 2!")
        end select

        dst (dstioff+i+1, dstjoff+j+1, dstkoff+k+1) = &
             s1fac * dsttmp1 + s2fac * dsttmp2 + s3fac * dsttmp3

!!$        write(*,*) i,j,k,dst (dstioff+i+1, dstjoff+j+1, dstkoff+k+1),&
!!$             s1fac,s2fac,s3fac,dsttmp1,dsttmp2,dsttmp3

        if ( (dst(dstioff+i+1, dstjoff+j+1, dstkoff+k+1) - &
             max(dsttmp1, dsttmp2, dsttmp3)) * &
             (dst(dstioff+i+1, dstjoff+j+1, dstkoff+k+1) - &
             min(dsttmp1, dsttmp2, dsttmp3)) .lt. 0 ) then

!!$          Do linear interpolation in time instead

!!$          write(*,*) t,t1,t2,t3

          if (t < t2) then

            tmps2fac = (t - t3) / (t2 - t3)
            tmps3fac = (t - t2) / (t3 - t2)

            dst (dstioff+i+1, dstjoff+j+1, dstkoff+k+1) = &
                 tmps2fac * dsttmp2 + tmps3fac * dsttmp3

          else
            
            tmps1fac = (t - t2) / (t1 - t2)
            tmps2fac = (t - t1) / (t2 - t1)
            
            dst (dstioff+i+1, dstjoff+j+1, dstkoff+k+1) = &
                 tmps1fac * dsttmp1 + tmps2fac * dsttmp2          

          end if
          
        end if
        
      end do
    end do
  end do

end subroutine prolongate_3d_real8_3tl_eno
