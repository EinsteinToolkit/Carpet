c     -*-Fortran-*-
c     $Header:$

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
      print *, "Xgfg ", Xgfg0, Xgfg1, Xgfg2
      print *, "Xarrg3 ", Xarrg30, Xarrg31, Xarrg32
      print *, "Xarrg2 ", Xarrg20, Xarrg21
      print *, "Xarrg1 ", Xarrg10
      print *, "Xscg"
      print *
      do k=1,Xgfg2
         do j=1,Xgfg1
            do i=1,Xgfg0
               gf(i,j,k) = i*10000 + j*100 + k
            end do
         end do
      end do
      do k=1,Xarrg32
         do j=1,Xarrg31
            do i=1,Xarrg30
               arr3(i,j,k) = i*10000 + j*100 + k
            end do
         end do
      end do
      do j=1,Xarrg21
         do i=1,Xarrg20
            arr2(i,j) = i*100 + j
         end do
      end do
      do i=1,Xarrg10
         arr1(i) = i
      end do
      sc = 42
      end
