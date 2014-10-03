c     -*-Fortran-*-

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
      
      subroutine carpettest_check_arguments (CCTK_ARGUMENTS)
      implicit none
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_FUNCTIONS
      DECLARE_CCTK_PARAMETERS
      integer i,j,k
      print *, "Xgfg ", cctk_lsh(1), cctk_lsh(2), cctk_lsh(3)
      print *, "Xarrg3 ", i3,i4,i5
      print *, "Xarrg2 ", i6,i7
      print *, "Xarrg1 ", i8
      print *, "Xscg"
      print *
      do k=1,cctk_lsh(3)
         do j=1,cctk_lsh(2)
            do i=1,cctk_lsh(1)
               gf(i,j,k) = i*10000 + j*100 + k
            end do
         end do
      end do
      do k=1,SIZE(arr3,3)
         do j=1,SIZE(arr3,2)
            do i=1,SIZE(arr3,1)
               arr3(i,j,k) = i*10000 + j*100 + k
            end do
         end do
      end do
      do j=1,SIZE(arr2,2)
         do i=1,SIZE(arr2,1)
            arr2(i,j) = i*100 + j
         end do
      end do
      do i=1,SIZE(arr1,1)
         arr1(i) = i
      end do
      sc = 42
      end
