c     -*-Fortran-*-
      
#include "cctk.h"
#include "cctk_Parameters.h"


      
      subroutine copy_3d_real8 (
     $     src, srciext, srcjext, srckext,
     $     dst, dstiext, dstjext, dstkext,
     $     srcbbox, dstbbox, regbbox)
      
      implicit none
      
      DECLARE_CCTK_PARAMETERS
      
      integer srciext, srcjext, srckext
      CCTK_REAL8 src(srciext,srcjext,srckext)
      integer dstiext, dstjext, dstkext
      CCTK_REAL8 dst(dstiext,dstjext,dstkext)
c     bbox(:,1) is lower boundary (inclusive)
c     bbox(:,2) is upper boundary (inclusive)
c     bbox(:,3) is stride
      integer srcbbox(3,3), dstbbox(3,3), regbbox(3,3)
      
      integer regiext, regjext, regkext
      
      integer srcioff, srcjoff, srckoff
      integer dstioff, dstjoff, dstkoff
      
      integer i, j, k
      integer d
      
      character msg*1000
      
      
      
      do d=1,3
         if (srcbbox(d,3).eq.0 .or. dstbbox(d,3).eq.0
     $        .or. regbbox(d,3).eq.0) then
            call CCTK_WARN (0, "Internal error: stride is zero")
         end if
         if (srcbbox(d,3).ne.regbbox(d,3)
     $        .or. dstbbox(d,3).ne.regbbox(d,3)) then
            call CCTK_WARN (0, "Internal error: strides disagree")
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
      
      srcioff = (regbbox(1,1) - srcbbox(1,1)) / srcbbox(1,3)
      srcjoff = (regbbox(2,1) - srcbbox(2,1)) / srcbbox(2,3)
      srckoff = (regbbox(3,1) - srcbbox(3,1)) / srcbbox(3,3)
      
      dstioff = (regbbox(1,1) - dstbbox(1,1)) / dstbbox(1,3)
      dstjoff = (regbbox(2,1) - dstbbox(2,1)) / dstbbox(2,3)
      dstkoff = (regbbox(3,1) - dstbbox(3,1)) / dstbbox(3,3)
      
      
      
c     Loop over region
      do k = 1, regkext
         do j = 1, regjext
            do i = 1, regiext
               
               if (check_array_accesses.ne.0) then
                  call checkindex (srcioff+i, srcjoff+j+1, srckoff+k+1, 1,1,1,
     $                 "source")
                  call checkindex (dstioff+i, dstjoff+j+1, dstkoff+k+1, 1,1,1,
     $                 "destination")
               end if
               
               dst (dstioff+i, dstjoff+j, dstkoff+k)
     $              = src (srcioff+i, srcjoff+j, srckoff+k)
               
            end do
         end do
      end do
      
      end
