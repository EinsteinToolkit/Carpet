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
      print *, "Xgfg ", cctk_lsh0, cctk_lsh1, cctk_lsh2
      print *, "Xarrg3 ", X0arrg3, X1arrg3, X2arrg3
      print *, "Xarrg2 ", X0arrg2, X1arrg2
      print *, "Xarrg1 ", X0arrg1
      print *, "Xscg"
      print *
      do k=1,cctk_lsh2
         do j=1,cctk_lsh1
            do i=1,cctk_lsh0
               gf(i,j,k) = i*10000 + j*100 + k
            end do
         end do
      end do
      do k=1,X2arrg3
         do j=1,X1arrg3
            do i=1,X0arrg3
               arr3(i,j,k) = i*10000 + j*100 + k
            end do
         end do
      end do
      do j=1,X1arrg2
         do i=1,X0arrg2
            arr2(i,j) = i*100 + j
         end do
      end do
      do i=1,X0arrg1
         arr1(i) = i
      end do
      sc = 42
      end
