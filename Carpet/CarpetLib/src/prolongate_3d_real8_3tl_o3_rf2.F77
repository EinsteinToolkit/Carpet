c     -*-Fortran-*-
c     $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/prolongate_3d_real8_3tl_o3_rf2.F77,v 1.1 2003/11/05 16:18:39 schnetter Exp $
      
#include "cctk.h"
#include "cctk_Parameters.h"
      


      subroutine prolongate_3d_real8_3tl_o3_rf2 (
     $     src1, t1, src2, t2, src3, t3, srciext, srcjext, srckext,
     $     dst, t, dstiext, dstjext, dstkext,
     $     srcbbox, dstbbox, regbbox)
      
      implicit none
      
      DECLARE_CCTK_PARAMETERS
      
      CCTK_REAL8 eps
      parameter (eps = 1.0d-10)
      
      CCTK_REAL8 one, half, fourth, eighth, sixteenth
      parameter (one = 1)
      parameter (half = one/2)
      parameter (fourth = one/4)
      parameter (eighth = one/8)
      parameter (sixteenth = one/16)
      CCTK_REAL8 f1, f2, f3, f4
      parameter (f1 = - sixteenth)
      parameter (f2 = 9*sixteenth)
      parameter (f3 = 9*sixteenth)
      parameter (f4 = - sixteenth)
      
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
      
      integer regiext, regjext, regkext
      
      integer srcioff, srcjoff, srckoff
      integer dstioff, dstjoff, dstkoff
      
      integer offsetlo, offsethi
      
      CCTK_REAL8 s1fac, s2fac, s3fac
      
      integer i0, j0, k0
      integer fi, fj, fk
      integer is, js, ks
      integer id, jd, kd
      integer i, j, k
      
      CCTK_REAL8 res1, res2, res3
      
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
         if (srcbbox(d,3).ne.dstbbox(d,3)*2) then
            call CCTK_WARN (0, "Internal error: source strides are not twice the destination strides")
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
         srckoff = (regbbox(d,1) - srcbbox(d,1)) / dstbbox(d,3)
         offsetlo = regbbox(d,3)
         if (mod(srckoff, 2).eq.0) then
            offsetlo = 0
            if (regkext.gt.1) then
               offsetlo = regbbox(d,3)
            end if
         end if
         offsethi = regbbox(d,3)
         if (mod(srckoff + regkext-1, 2).eq.0) then
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
      
      srcioff = (regbbox(1,1) - srcbbox(1,1)) / dstbbox(1,3)
      srcjoff = (regbbox(2,1) - srcbbox(2,1)) / dstbbox(2,3)
      srckoff = (regbbox(3,1) - srcbbox(3,1)) / dstbbox(3,3)
      
      dstioff = (regbbox(1,1) - dstbbox(1,1)) / dstbbox(1,3)
      dstjoff = (regbbox(2,1) - dstbbox(2,1)) / dstbbox(2,3)
      dstkoff = (regbbox(3,1) - dstbbox(3,1)) / dstbbox(3,3)
      
      
      
c     Quadratic (second order) time interpolation
      if (t1.eq.t2 .or. t1.eq.t3 .or. t2.eq.t3) then
         call CCTK_WARN (0, "Internal error: arrays have same time")
      end if
      if (t.lt.min(t1,t2,t3)-eps .or. t.gt.max(t1,t2,t3)+eps) then
         call CCTK_WARN (0, "Internal error: extrapolation")
      end if
      
      s1fac = (t - t2) * (t - t3) / ((t1 - t2) * (t1 - t3))
      s2fac = (t - t1) * (t - t3) / ((t2 - t1) * (t2 - t3))
      s3fac = (t - t1) * (t - t2) / ((t3 - t1) * (t3 - t2))
      
      
      
      fi = mod(srcioff, 2)
      fj = mod(srcjoff, 2)
      fk = mod(srckoff, 2)
      
      i0 = srcioff / 2
      j0 = srcjoff / 2
      k0 = srckoff / 2
      
      
      
c     Loop over fine region
c     Label scheme: 8 fk fj fi
      
c     begin k loop
 8    continue
      k = 0
      ks = k0+1
      kd = dstkoff+1
      if (fk.eq.0) goto 80
      if (fk.eq.1) goto 81
      stop
      
c     begin j loop
 80   continue
      j = 0
      js = j0+1
      jd = dstjoff+1
      if (fj.eq.0) goto 800
      if (fj.eq.1) goto 801
      stop
      
c     begin i loop
 800  continue
      i = 0
      is = i0+1
      id = dstioff+1
      if (fi.eq.0) goto 8000
      if (fi.eq.1) goto 8001
      stop
      
c     kernel
 8000 continue
      if (check_array_accesses.ne.0) then
         call checkindex (is,js,ks, 1,1,1, srciext,srcjext,srckext, "source")
         call checkindex (id,jd,kd, 1,1,1, dstiext,dstjext,dstkext, "destination")
      end if
      dst(id,jd,kd) =
     $     + s1fac * src1(is,js,ks)
     $     + s2fac * src2(is,js,ks)
     $     + s3fac * src3(is,js,ks)
      i = i+1
      id = id+1
      if (i.lt.regiext) goto 8001
      goto 900
      
c     kernel
 8001 continue
      if (check_array_accesses.ne.0) then
         call checkindex (is-1,js,ks, 4,1,1, srciext,srcjext,srckext, "source")
         call checkindex (id,jd,kd, 1,1,1, dstiext,dstjext,dstkext, "destination")
      end if
      dst(id,jd,kd) =
     $     + f1 * s1fac * src1(is-1,js,ks) + f2 * s1fac * src1(is  ,js,ks)
     $     + f3 * s1fac * src1(is+1,js,ks) + f4 * s1fac * src1(is+2,js,ks)
     $     + f1 * s2fac * src2(is-1,js,ks) + f2 * s2fac * src2(is  ,js,ks)
     $     + f3 * s2fac * src2(is+1,js,ks) + f4 * s2fac * src2(is+2,js,ks)
     $     + f1 * s3fac * src3(is-1,js,ks) + f2 * s3fac * src3(is  ,js,ks)
     $     + f3 * s3fac * src3(is+1,js,ks) + f4 * s3fac * src3(is+2,js,ks)
      i = i+1
      id = id+1
      is = is+1
      if (i.lt.regiext) goto 8000
      goto 900
      
c     end i loop
 900  continue
      j = j+1
      jd = jd+1
      if (j.lt.regjext) goto 801
      goto 90
      
c     begin i loop
 801  continue
      i = 0
      is = i0+1
      id = dstioff+1
      if (fi.eq.0) goto 8010
      if (fi.eq.1) goto 8011
      stop
      
c     kernel
 8010 continue
      if (check_array_accesses.ne.0) then
         call checkindex (is,js-1,ks, 1,4,1, srciext,srcjext,srckext, "source")
         call checkindex (id,jd,kd, 1,1,1, dstiext,dstjext,dstkext, "destination")
      end if
      dst(id,jd,kd) =
     $     + f1 * s1fac * src1(is,js-1,ks) + f2 * s1fac * src1(is,js  ,ks)
     $     + f3 * s1fac * src1(is,js+1,ks) + f4 * s1fac * src1(is,js+2,ks)
     $     + f1 * s2fac * src2(is,js-1,ks) + f2 * s2fac * src2(is,js  ,ks)
     $     + f3 * s2fac * src2(is,js+1,ks) + f4 * s2fac * src2(is,js+2,ks)
     $     + f1 * s3fac * src3(is,js-1,ks) + f2 * s3fac * src3(is,js  ,ks)
     $     + f3 * s3fac * src3(is,js+1,ks) + f4 * s3fac * src3(is,js+2,ks)
      i = i+1
      id = id+1
      if (i.lt.regiext) goto 8011
      goto 901
      
c     kernel
 8011 continue
      if (check_array_accesses.ne.0) then
         call checkindex (is-1,js-1,ks, 4,4,1, srciext,srcjext,srckext, "source")
         call checkindex (id,jd,kd, 1,1,1, dstiext,dstjext,dstkext, "destination")
      end if
      dst(id,jd,kd) =
     $     + f1*f1 * s1fac * src1(is-1,js-1,ks)
     $     + f2*f1 * s1fac * src1(is  ,js-1,ks)
     $     + f3*f1 * s1fac * src1(is+1,js-1,ks)
     $     + f4*f1 * s1fac * src1(is+2,js-1,ks)
     $     + f1*f2 * s1fac * src1(is-1,js  ,ks)
     $     + f2*f2 * s1fac * src1(is  ,js  ,ks)
     $     + f3*f2 * s1fac * src1(is+1,js  ,ks)
     $     + f4*f2 * s1fac * src1(is+2,js  ,ks)
     $     + f1*f3 * s1fac * src1(is-1,js+1,ks)
     $     + f2*f3 * s1fac * src1(is  ,js+1,ks)
     $     + f3*f3 * s1fac * src1(is+1,js+1,ks)
     $     + f4*f3 * s1fac * src1(is+2,js+1,ks)
     $     + f1*f4 * s1fac * src1(is-1,js+2,ks)
     $     + f2*f4 * s1fac * src1(is  ,js+2,ks)
     $     + f3*f4 * s1fac * src1(is+1,js+2,ks)
     $     + f4*f4 * s1fac * src1(is+2,js+2,ks)
     $     
     $     + f1*f1 * s2fac * src2(is-1,js-1,ks)
     $     + f2*f1 * s2fac * src2(is  ,js-1,ks)
     $     + f3*f1 * s2fac * src2(is+1,js-1,ks)
     $     + f4*f1 * s2fac * src2(is+2,js-1,ks)
     $     + f1*f2 * s2fac * src2(is-1,js  ,ks)
     $     + f2*f2 * s2fac * src2(is  ,js  ,ks)
     $     + f3*f2 * s2fac * src2(is+1,js  ,ks)
     $     + f4*f2 * s2fac * src2(is+2,js  ,ks)
     $     + f1*f3 * s2fac * src2(is-1,js+1,ks)
     $     + f2*f3 * s2fac * src2(is  ,js+1,ks)
     $     + f3*f3 * s2fac * src2(is+1,js+1,ks)
     $     + f4*f3 * s2fac * src2(is+2,js+1,ks)
     $     + f1*f4 * s2fac * src2(is-1,js+2,ks)
     $     + f2*f4 * s2fac * src2(is  ,js+2,ks)
     $     + f3*f4 * s2fac * src2(is+1,js+2,ks)
     $     + f4*f4 * s2fac * src2(is+2,js+2,ks)
     $     
     $     + f1*f1 * s3fac * src3(is-1,js-1,ks)
     $     + f2*f1 * s3fac * src3(is  ,js-1,ks)
     $     + f3*f1 * s3fac * src3(is+1,js-1,ks)
     $     + f4*f1 * s3fac * src3(is+2,js-1,ks)
     $     + f1*f2 * s3fac * src3(is-1,js  ,ks)
     $     + f2*f2 * s3fac * src3(is  ,js  ,ks)
     $     + f3*f2 * s3fac * src3(is+1,js  ,ks)
     $     + f4*f2 * s3fac * src3(is+2,js  ,ks)
     $     + f1*f3 * s3fac * src3(is-1,js+1,ks)
     $     + f2*f3 * s3fac * src3(is  ,js+1,ks)
     $     + f3*f3 * s3fac * src3(is+1,js+1,ks)
     $     + f4*f3 * s3fac * src3(is+2,js+1,ks)
     $     + f1*f4 * s3fac * src3(is-1,js+2,ks)
     $     + f2*f4 * s3fac * src3(is  ,js+2,ks)
     $     + f3*f4 * s3fac * src3(is+1,js+2,ks)
     $     + f4*f4 * s3fac * src3(is+2,js+2,ks)
      i = i+1
      id = id+1
      is = is+1
      if (i.lt.regiext) goto 8010
      goto 901
      
c     end i loop
 901  continue
      j = j+1
      jd = jd+1
      js = js+1
      if (j.lt.regjext) goto 800
      goto 90
      
c     end j loop
 90   continue
      k = k+1
      kd = kd+1
      if (k.lt.regkext) goto 81
      goto 9
      
c     begin j loop
 81   continue
      j = 0
      js = j0+1
      jd = dstjoff+1
      if (fj.eq.0) goto 810
      if (fj.eq.1) goto 811
      stop
      
c     begin i loop
 810  continue
      i = 0
      is = i0+1
      id = dstioff+1
      if (fi.eq.0) goto 8100
      if (fi.eq.1) goto 8101
      stop
      
c     kernel
 8100 continue
      if (check_array_accesses.ne.0) then
         call checkindex (is,js,ks-1, 1,1,4, srciext,srcjext,srckext, "source")
         call checkindex (id,jd,kd, 1,1,1, dstiext,dstjext,dstkext, "destination")
      end if
      dst(id,jd,kd) =
     $     + f1 * s1fac * src1(is,js,ks-1) + f2 * s1fac * src1(is,js,ks  )
     $     + f3 * s1fac * src1(is,js,ks+1) + f4 * s1fac * src1(is,js,ks+2)
     $     + f1 * s2fac * src2(is,js,ks-1) + f2 * s2fac * src2(is,js,ks  )
     $     + f3 * s2fac * src2(is,js,ks+1) + f4 * s2fac * src2(is,js,ks+2)
     $     + f1 * s3fac * src3(is,js,ks-1) + f2 * s3fac * src3(is,js,ks  )
     $     + f3 * s3fac * src3(is,js,ks+1) + f4 * s3fac * src3(is,js,ks+2)
      i = i+1
      id = id+1
      if (i.lt.regiext) goto 8101
      goto 910
      
c     kernel
 8101 continue
      if (check_array_accesses.ne.0) then
         call checkindex (is-1,js,ks-1, 4,1,4, srciext,srcjext,srckext, "source")
         call checkindex (id,jd,kd, 1,1,1, dstiext,dstjext,dstkext, "destination")
      end if
      dst(id,jd,kd) =
     $     + f1*f1 * s1fac * src1(is-1,js,ks-1)
     $     + f2*f1 * s1fac * src1(is  ,js,ks-1)
     $     + f3*f1 * s1fac * src1(is+1,js,ks-1)
     $     + f4*f1 * s1fac * src1(is+2,js,ks-1)
     $     + f1*f2 * s1fac * src1(is-1,js,ks  )
     $     + f2*f2 * s1fac * src1(is  ,js,ks  )
     $     + f3*f2 * s1fac * src1(is+1,js,ks  )
     $     + f4*f2 * s1fac * src1(is+2,js,ks  )
     $     + f1*f3 * s1fac * src1(is-1,js,ks+1)
     $     + f2*f3 * s1fac * src1(is  ,js,ks+1)
     $     + f3*f3 * s1fac * src1(is+1,js,ks+1)
     $     + f4*f3 * s1fac * src1(is+2,js,ks+1)
     $     + f1*f4 * s1fac * src1(is-1,js,ks+2)
     $     + f2*f4 * s1fac * src1(is  ,js,ks+2)
     $     + f3*f4 * s1fac * src1(is+1,js,ks+2)
     $     + f4*f4 * s1fac * src1(is+2,js,ks+2)
     $     
     $     + f1*f1 * s2fac * src2(is-1,js,ks-1)
     $     + f2*f1 * s2fac * src2(is  ,js,ks-1)
     $     + f3*f1 * s2fac * src2(is+1,js,ks-1)
     $     + f4*f1 * s2fac * src2(is+2,js,ks-1)
     $     + f1*f2 * s2fac * src2(is-1,js,ks  )
     $     + f2*f2 * s2fac * src2(is  ,js,ks  )
     $     + f3*f2 * s2fac * src2(is+1,js,ks  )
     $     + f4*f2 * s2fac * src2(is+2,js,ks  )
     $     + f1*f3 * s2fac * src2(is-1,js,ks+1)
     $     + f2*f3 * s2fac * src2(is  ,js,ks+1)
     $     + f3*f3 * s2fac * src2(is+1,js,ks+1)
     $     + f4*f3 * s2fac * src2(is+2,js,ks+1)
     $     + f1*f4 * s2fac * src2(is-1,js,ks+2)
     $     + f2*f4 * s2fac * src2(is  ,js,ks+2)
     $     + f3*f4 * s2fac * src2(is+1,js,ks+2)
     $     + f4*f4 * s2fac * src2(is+2,js,ks+2)
     $     
     $     + f1*f1 * s3fac * src3(is-1,js,ks-1)
     $     + f2*f1 * s3fac * src3(is  ,js,ks-1)
     $     + f3*f1 * s3fac * src3(is+1,js,ks-1)
     $     + f4*f1 * s3fac * src3(is+2,js,ks-1)
     $     + f1*f2 * s3fac * src3(is-1,js,ks  )
     $     + f2*f2 * s3fac * src3(is  ,js,ks  )
     $     + f3*f2 * s3fac * src3(is+1,js,ks  )
     $     + f4*f2 * s3fac * src3(is+2,js,ks  )
     $     + f1*f3 * s3fac * src3(is-1,js,ks+1)
     $     + f2*f3 * s3fac * src3(is  ,js,ks+1)
     $     + f3*f3 * s3fac * src3(is+1,js,ks+1)
     $     + f4*f3 * s3fac * src3(is+2,js,ks+1)
     $     + f1*f4 * s3fac * src3(is-1,js,ks+2)
     $     + f2*f4 * s3fac * src3(is  ,js,ks+2)
     $     + f3*f4 * s3fac * src3(is+1,js,ks+2)
     $     + f4*f4 * s3fac * src3(is+2,js,ks+2)
      i = i+1
      id = id+1
      is = is+1
      if (i.lt.regiext) goto 8100
      goto 910
      
c     end i loop
 910  continue
      j = j+1
      jd = jd+1
      if (j.lt.regjext) goto 811
      goto 91
      
c     begin i loop
 811  continue
      i = 0
      is = i0+1
      id = dstioff+1
      if (fi.eq.0) goto 8110
      if (fi.eq.1) goto 8111
      stop
      
c     kernel
 8110 continue
      if (check_array_accesses.ne.0) then
         call checkindex (is,js-1,ks-1, 1,4,4, srciext,srcjext,srckext, "source")
         call checkindex (id,jd,kd, 1,1,1, dstiext,dstjext,dstkext, "destination")
      end if
      dst(id,jd,kd) =
     $     + f1*f1 * s1fac * src1(is,js-1,ks-1)
     $     + f2*f1 * s1fac * src1(is,js  ,ks-1)
     $     + f3*f1 * s1fac * src1(is,js+1,ks-1)
     $     + f4*f1 * s1fac * src1(is,js+2,ks-1)
     $     + f1*f2 * s1fac * src1(is,js-1,ks  )
     $     + f2*f2 * s1fac * src1(is,js  ,ks  )
     $     + f3*f2 * s1fac * src1(is,js+1,ks  )
     $     + f4*f2 * s1fac * src1(is,js+2,ks  )
     $     + f1*f3 * s1fac * src1(is,js-1,ks+1)
     $     + f2*f3 * s1fac * src1(is,js  ,ks+1)
     $     + f3*f3 * s1fac * src1(is,js+1,ks+1)
     $     + f4*f3 * s1fac * src1(is,js+2,ks+1)
     $     + f1*f4 * s1fac * src1(is,js-1,ks+2)
     $     + f2*f4 * s1fac * src1(is,js  ,ks+2)
     $     + f3*f4 * s1fac * src1(is,js+1,ks+2)
     $     + f4*f4 * s1fac * src1(is,js+2,ks+2)
     $     
     $     + f1*f1 * s2fac * src2(is,js-1,ks-1)
     $     + f2*f1 * s2fac * src2(is,js  ,ks-1)
     $     + f3*f1 * s2fac * src2(is,js+1,ks-1)
     $     + f4*f1 * s2fac * src2(is,js+2,ks-1)
     $     + f1*f2 * s2fac * src2(is,js-1,ks  )
     $     + f2*f2 * s2fac * src2(is,js  ,ks  )
     $     + f3*f2 * s2fac * src2(is,js+1,ks  )
     $     + f4*f2 * s2fac * src2(is,js+2,ks  )
     $     + f1*f3 * s2fac * src2(is,js-1,ks+1)
     $     + f2*f3 * s2fac * src2(is,js  ,ks+1)
     $     + f3*f3 * s2fac * src2(is,js+1,ks+1)
     $     + f4*f3 * s2fac * src2(is,js+2,ks+1)
     $     + f1*f4 * s2fac * src2(is,js-1,ks+2)
     $     + f2*f4 * s2fac * src2(is,js  ,ks+2)
     $     + f3*f4 * s2fac * src2(is,js+1,ks+2)
     $     + f4*f4 * s2fac * src2(is,js+2,ks+2)
     $     
     $     + f1*f1 * s3fac * src3(is,js-1,ks-1)
     $     + f2*f1 * s3fac * src3(is,js  ,ks-1)
     $     + f3*f1 * s3fac * src3(is,js+1,ks-1)
     $     + f4*f1 * s3fac * src3(is,js+2,ks-1)
     $     + f1*f2 * s3fac * src3(is,js-1,ks  )
     $     + f2*f2 * s3fac * src3(is,js  ,ks  )
     $     + f3*f2 * s3fac * src3(is,js+1,ks  )
     $     + f4*f2 * s3fac * src3(is,js+2,ks  )
     $     + f1*f3 * s3fac * src3(is,js-1,ks+1)
     $     + f2*f3 * s3fac * src3(is,js  ,ks+1)
     $     + f3*f3 * s3fac * src3(is,js+1,ks+1)
     $     + f4*f3 * s3fac * src3(is,js+2,ks+1)
     $     + f1*f4 * s3fac * src3(is,js-1,ks+2)
     $     + f2*f4 * s3fac * src3(is,js  ,ks+2)
     $     + f3*f4 * s3fac * src3(is,js+1,ks+2)
     $     + f4*f4 * s3fac * src3(is,js+2,ks+2)
      i = i+1
      id = id+1
      if (i.lt.regiext) goto 8111
      goto 911
      
c     kernel
 8111 continue
      if (check_array_accesses.ne.0) then
         call checkindex (is-1,js-1,ks-1, 4,4,4, srciext,srcjext,srckext, "source")
         call checkindex (id,jd,kd, 1,1,1, dstiext,dstjext,dstkext, "destination")
      end if
      res1 =
     $     + f1*f1*f1 * s1fac * src1(is-1,js-1,ks-1)
     $     + f2*f1*f1 * s1fac * src1(is  ,js-1,ks-1)
     $     + f3*f1*f1 * s1fac * src1(is+1,js-1,ks-1)
     $     + f4*f1*f1 * s1fac * src1(is+2,js-1,ks-1)
     $     + f1*f2*f1 * s1fac * src1(is-1,js  ,ks-1)
     $     + f2*f2*f1 * s1fac * src1(is  ,js  ,ks-1)
     $     + f3*f2*f1 * s1fac * src1(is+1,js  ,ks-1)
     $     + f4*f2*f1 * s1fac * src1(is+2,js  ,ks-1)
     $     + f1*f3*f1 * s1fac * src1(is-1,js+1,ks-1)
     $     + f2*f3*f1 * s1fac * src1(is  ,js+1,ks-1)
     $     + f3*f3*f1 * s1fac * src1(is+1,js+1,ks-1)
     $     + f4*f3*f1 * s1fac * src1(is+2,js+1,ks-1)
     $     + f1*f4*f1 * s1fac * src1(is-1,js+2,ks-1)
     $     + f2*f4*f1 * s1fac * src1(is  ,js+2,ks-1)
     $     + f3*f4*f1 * s1fac * src1(is+1,js+2,ks-1)
     $     + f4*f4*f1 * s1fac * src1(is+2,js+2,ks-1)
     $     
     $     + f1*f1*f2 * s1fac * src1(is-1,js-1,ks  )
     $     + f2*f1*f2 * s1fac * src1(is  ,js-1,ks  )
     $     + f3*f1*f2 * s1fac * src1(is+1,js-1,ks  )
     $     + f4*f1*f2 * s1fac * src1(is+2,js-1,ks  )
     $     + f1*f2*f2 * s1fac * src1(is-1,js  ,ks  )
     $     + f2*f2*f2 * s1fac * src1(is  ,js  ,ks  )
     $     + f3*f2*f2 * s1fac * src1(is+1,js  ,ks  )
     $     + f4*f2*f2 * s1fac * src1(is+2,js  ,ks  )
     $     + f1*f3*f2 * s1fac * src1(is-1,js+1,ks  )
     $     + f2*f3*f2 * s1fac * src1(is  ,js+1,ks  )
     $     + f3*f3*f2 * s1fac * src1(is+1,js+1,ks  )
     $     + f4*f3*f2 * s1fac * src1(is+2,js+1,ks  )
     $     + f1*f4*f2 * s1fac * src1(is-1,js+2,ks  )
     $     + f2*f4*f2 * s1fac * src1(is  ,js+2,ks  )
     $     + f3*f4*f2 * s1fac * src1(is+1,js+2,ks  )
     $     + f4*f4*f2 * s1fac * src1(is+2,js+2,ks  )
     $     
     $     + f1*f1*f3 * s1fac * src1(is-1,js-1,ks+1)
     $     + f2*f1*f3 * s1fac * src1(is  ,js-1,ks+1)
     $     + f3*f1*f3 * s1fac * src1(is+1,js-1,ks+1)
     $     + f4*f1*f3 * s1fac * src1(is+2,js-1,ks+1)
     $     + f1*f2*f3 * s1fac * src1(is-1,js  ,ks+1)
     $     + f2*f2*f3 * s1fac * src1(is  ,js  ,ks+1)
     $     + f3*f2*f3 * s1fac * src1(is+1,js  ,ks+1)
     $     + f4*f2*f3 * s1fac * src1(is+2,js  ,ks+1)
     $     + f1*f3*f3 * s1fac * src1(is-1,js+1,ks+1)
     $     + f2*f3*f3 * s1fac * src1(is  ,js+1,ks+1)
     $     + f3*f3*f3 * s1fac * src1(is+1,js+1,ks+1)
     $     + f4*f3*f3 * s1fac * src1(is+2,js+1,ks+1)
     $     + f1*f4*f3 * s1fac * src1(is-1,js+2,ks+1)
     $     + f2*f4*f3 * s1fac * src1(is  ,js+2,ks+1)
     $     + f3*f4*f3 * s1fac * src1(is+1,js+2,ks+1)
     $     + f4*f4*f3 * s1fac * src1(is+2,js+2,ks+1)
     $     
     $     + f1*f1*f4 * s1fac * src1(is-1,js-1,ks+2)
     $     + f2*f1*f4 * s1fac * src1(is  ,js-1,ks+2)
     $     + f3*f1*f4 * s1fac * src1(is+1,js-1,ks+2)
     $     + f4*f1*f4 * s1fac * src1(is+2,js-1,ks+2)
     $     + f1*f2*f4 * s1fac * src1(is-1,js  ,ks+2)
     $     + f2*f2*f4 * s1fac * src1(is  ,js  ,ks+2)
     $     + f3*f2*f4 * s1fac * src1(is+1,js  ,ks+2)
     $     + f4*f2*f4 * s1fac * src1(is+2,js  ,ks+2)
     $     + f1*f3*f4 * s1fac * src1(is-1,js+1,ks+2)
     $     + f2*f3*f4 * s1fac * src1(is  ,js+1,ks+2)
     $     + f3*f3*f4 * s1fac * src1(is+1,js+1,ks+2)
     $     + f4*f3*f4 * s1fac * src1(is+2,js+1,ks+2)
     $     + f1*f4*f4 * s1fac * src1(is-1,js+2,ks+2)
     $     + f2*f4*f4 * s1fac * src1(is  ,js+2,ks+2)
     $     + f3*f4*f4 * s1fac * src1(is+1,js+2,ks+2)
     $     + f4*f4*f4 * s1fac * src1(is+2,js+2,ks+2)
      res2 =
     $     + f1*f1*f1 * s2fac * src2(is-1,js-1,ks-1)
     $     + f2*f1*f1 * s2fac * src2(is  ,js-1,ks-1)
     $     + f3*f1*f1 * s2fac * src2(is+1,js-1,ks-1)
     $     + f4*f1*f1 * s2fac * src2(is+2,js-1,ks-1)
     $     + f1*f2*f1 * s2fac * src2(is-1,js  ,ks-1)
     $     + f2*f2*f1 * s2fac * src2(is  ,js  ,ks-1)
     $     + f3*f2*f1 * s2fac * src2(is+1,js  ,ks-1)
     $     + f4*f2*f1 * s2fac * src2(is+2,js  ,ks-1)
     $     + f1*f3*f1 * s2fac * src2(is-1,js+1,ks-1)
     $     + f2*f3*f1 * s2fac * src2(is  ,js+1,ks-1)
     $     + f3*f3*f1 * s2fac * src2(is+1,js+1,ks-1)
     $     + f4*f3*f1 * s2fac * src2(is+2,js+1,ks-1)
     $     + f1*f4*f1 * s2fac * src2(is-1,js+2,ks-1)
     $     + f2*f4*f1 * s2fac * src2(is  ,js+2,ks-1)
     $     + f3*f4*f1 * s2fac * src2(is+1,js+2,ks-1)
     $     + f4*f4*f1 * s2fac * src2(is+2,js+2,ks-1)
     $     
     $     + f1*f1*f2 * s2fac * src2(is-1,js-1,ks  )
     $     + f2*f1*f2 * s2fac * src2(is  ,js-1,ks  )
     $     + f3*f1*f2 * s2fac * src2(is+1,js-1,ks  )
     $     + f4*f1*f2 * s2fac * src2(is+2,js-1,ks  )
     $     + f1*f2*f2 * s2fac * src2(is-1,js  ,ks  )
     $     + f2*f2*f2 * s2fac * src2(is  ,js  ,ks  )
     $     + f3*f2*f2 * s2fac * src2(is+1,js  ,ks  )
     $     + f4*f2*f2 * s2fac * src2(is+2,js  ,ks  )
     $     + f1*f3*f2 * s2fac * src2(is-1,js+1,ks  )
     $     + f2*f3*f2 * s2fac * src2(is  ,js+1,ks  )
     $     + f3*f3*f2 * s2fac * src2(is+1,js+1,ks  )
     $     + f4*f3*f2 * s2fac * src2(is+2,js+1,ks  )
     $     + f1*f4*f2 * s2fac * src2(is-1,js+2,ks  )
     $     + f2*f4*f2 * s2fac * src2(is  ,js+2,ks  )
     $     + f3*f4*f2 * s2fac * src2(is+1,js+2,ks  )
     $     + f4*f4*f2 * s2fac * src2(is+2,js+2,ks  )
     $     
     $     + f1*f1*f3 * s2fac * src2(is-1,js-1,ks+1)
     $     + f2*f1*f3 * s2fac * src2(is  ,js-1,ks+1)
     $     + f3*f1*f3 * s2fac * src2(is+1,js-1,ks+1)
     $     + f4*f1*f3 * s2fac * src2(is+2,js-1,ks+1)
     $     + f1*f2*f3 * s2fac * src2(is-1,js  ,ks+1)
     $     + f2*f2*f3 * s2fac * src2(is  ,js  ,ks+1)
     $     + f3*f2*f3 * s2fac * src2(is+1,js  ,ks+1)
     $     + f4*f2*f3 * s2fac * src2(is+2,js  ,ks+1)
     $     + f1*f3*f3 * s2fac * src2(is-1,js+1,ks+1)
     $     + f2*f3*f3 * s2fac * src2(is  ,js+1,ks+1)
     $     + f3*f3*f3 * s2fac * src2(is+1,js+1,ks+1)
     $     + f4*f3*f3 * s2fac * src2(is+2,js+1,ks+1)
     $     + f1*f4*f3 * s2fac * src2(is-1,js+2,ks+1)
     $     + f2*f4*f3 * s2fac * src2(is  ,js+2,ks+1)
     $     + f3*f4*f3 * s2fac * src2(is+1,js+2,ks+1)
     $     + f4*f4*f3 * s2fac * src2(is+2,js+2,ks+1)
     $     
     $     + f1*f1*f4 * s2fac * src2(is-1,js-1,ks+2)
     $     + f2*f1*f4 * s2fac * src2(is  ,js-1,ks+2)
     $     + f3*f1*f4 * s2fac * src2(is+1,js-1,ks+2)
     $     + f4*f1*f4 * s2fac * src2(is+2,js-1,ks+2)
     $     + f1*f2*f4 * s2fac * src2(is-1,js  ,ks+2)
     $     + f2*f2*f4 * s2fac * src2(is  ,js  ,ks+2)
     $     + f3*f2*f4 * s2fac * src2(is+1,js  ,ks+2)
     $     + f4*f2*f4 * s2fac * src2(is+2,js  ,ks+2)
     $     + f1*f3*f4 * s2fac * src2(is-1,js+1,ks+2)
     $     + f2*f3*f4 * s2fac * src2(is  ,js+1,ks+2)
     $     + f3*f3*f4 * s2fac * src2(is+1,js+1,ks+2)
     $     + f4*f3*f4 * s2fac * src2(is+2,js+1,ks+2)
     $     + f1*f4*f4 * s2fac * src2(is-1,js+2,ks+2)
     $     + f2*f4*f4 * s2fac * src2(is  ,js+2,ks+2)
     $     + f3*f4*f4 * s2fac * src2(is+1,js+2,ks+2)
     $     + f4*f4*f4 * s2fac * src2(is+2,js+2,ks+2)
      res3 =
     $     + f1*f1*f1 * s3fac * src3(is-1,js-1,ks-1)
     $     + f2*f1*f1 * s3fac * src3(is  ,js-1,ks-1)
     $     + f3*f1*f1 * s3fac * src3(is+1,js-1,ks-1)
     $     + f4*f1*f1 * s3fac * src3(is+2,js-1,ks-1)
     $     + f1*f2*f1 * s3fac * src3(is-1,js  ,ks-1)
     $     + f2*f2*f1 * s3fac * src3(is  ,js  ,ks-1)
     $     + f3*f2*f1 * s3fac * src3(is+1,js  ,ks-1)
     $     + f4*f2*f1 * s3fac * src3(is+2,js  ,ks-1)
     $     + f1*f3*f1 * s3fac * src3(is-1,js+1,ks-1)
     $     + f2*f3*f1 * s3fac * src3(is  ,js+1,ks-1)
     $     + f3*f3*f1 * s3fac * src3(is+1,js+1,ks-1)
     $     + f4*f3*f1 * s3fac * src3(is+2,js+1,ks-1)
     $     + f1*f4*f1 * s3fac * src3(is-1,js+2,ks-1)
     $     + f2*f4*f1 * s3fac * src3(is  ,js+2,ks-1)
     $     + f3*f4*f1 * s3fac * src3(is+1,js+2,ks-1)
     $     + f4*f4*f1 * s3fac * src3(is+2,js+2,ks-1)
     $     
     $     + f1*f1*f2 * s3fac * src3(is-1,js-1,ks  )
     $     + f2*f1*f2 * s3fac * src3(is  ,js-1,ks  )
     $     + f3*f1*f2 * s3fac * src3(is+1,js-1,ks  )
     $     + f4*f1*f2 * s3fac * src3(is+2,js-1,ks  )
     $     + f1*f2*f2 * s3fac * src3(is-1,js  ,ks  )
     $     + f2*f2*f2 * s3fac * src3(is  ,js  ,ks  )
     $     + f3*f2*f2 * s3fac * src3(is+1,js  ,ks  )
     $     + f4*f2*f2 * s3fac * src3(is+2,js  ,ks  )
     $     + f1*f3*f2 * s3fac * src3(is-1,js+1,ks  )
     $     + f2*f3*f2 * s3fac * src3(is  ,js+1,ks  )
     $     + f3*f3*f2 * s3fac * src3(is+1,js+1,ks  )
     $     + f4*f3*f2 * s3fac * src3(is+2,js+1,ks  )
     $     + f1*f4*f2 * s3fac * src3(is-1,js+2,ks  )
     $     + f2*f4*f2 * s3fac * src3(is  ,js+2,ks  )
     $     + f3*f4*f2 * s3fac * src3(is+1,js+2,ks  )
     $     + f4*f4*f2 * s3fac * src3(is+2,js+2,ks  )
     $     
     $     + f1*f1*f3 * s3fac * src3(is-1,js-1,ks+1)
     $     + f2*f1*f3 * s3fac * src3(is  ,js-1,ks+1)
     $     + f3*f1*f3 * s3fac * src3(is+1,js-1,ks+1)
     $     + f4*f1*f3 * s3fac * src3(is+2,js-1,ks+1)
     $     + f1*f2*f3 * s3fac * src3(is-1,js  ,ks+1)
     $     + f2*f2*f3 * s3fac * src3(is  ,js  ,ks+1)
     $     + f3*f2*f3 * s3fac * src3(is+1,js  ,ks+1)
     $     + f4*f2*f3 * s3fac * src3(is+2,js  ,ks+1)
     $     + f1*f3*f3 * s3fac * src3(is-1,js+1,ks+1)
     $     + f2*f3*f3 * s3fac * src3(is  ,js+1,ks+1)
     $     + f3*f3*f3 * s3fac * src3(is+1,js+1,ks+1)
     $     + f4*f3*f3 * s3fac * src3(is+2,js+1,ks+1)
     $     + f1*f4*f3 * s3fac * src3(is-1,js+2,ks+1)
     $     + f2*f4*f3 * s3fac * src3(is  ,js+2,ks+1)
     $     + f3*f4*f3 * s3fac * src3(is+1,js+2,ks+1)
     $     + f4*f4*f3 * s3fac * src3(is+2,js+2,ks+1)
     $     
     $     + f1*f1*f4 * s3fac * src3(is-1,js-1,ks+2)
     $     + f2*f1*f4 * s3fac * src3(is  ,js-1,ks+2)
     $     + f3*f1*f4 * s3fac * src3(is+1,js-1,ks+2)
     $     + f4*f1*f4 * s3fac * src3(is+2,js-1,ks+2)
     $     + f1*f2*f4 * s3fac * src3(is-1,js  ,ks+2)
     $     + f2*f2*f4 * s3fac * src3(is  ,js  ,ks+2)
     $     + f3*f2*f4 * s3fac * src3(is+1,js  ,ks+2)
     $     + f4*f2*f4 * s3fac * src3(is+2,js  ,ks+2)
     $     + f1*f3*f4 * s3fac * src3(is-1,js+1,ks+2)
     $     + f2*f3*f4 * s3fac * src3(is  ,js+1,ks+2)
     $     + f3*f3*f4 * s3fac * src3(is+1,js+1,ks+2)
     $     + f4*f3*f4 * s3fac * src3(is+2,js+1,ks+2)
     $     + f1*f4*f4 * s3fac * src3(is-1,js+2,ks+2)
     $     + f2*f4*f4 * s3fac * src3(is  ,js+2,ks+2)
     $     + f3*f4*f4 * s3fac * src3(is+1,js+2,ks+2)
     $     + f4*f4*f4 * s3fac * src3(is+2,js+2,ks+2)
      dst(id,jd,kd) = res1 + res2 + res3
      i = i+1
      id = id+1
      is = is+1
      if (i.lt.regiext) goto 8110
      goto 911
      
c     end i loop
 911  continue
      j = j+1
      jd = jd+1
      js = js+1
      if (j.lt.regjext) goto 810
      goto 91
      
c     end j loop
 91   continue
      k = k+1
      kd = kd+1
      ks = ks+1
      if (k.lt.regkext) goto 80
      goto 9
      
c     end k loop
 9    continue
      
      end
