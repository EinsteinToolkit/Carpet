c     -*-Fortran-*-
c     $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetTest/src/carpettest_check_arguments.F77,v 1.3 2001/07/09 09:00:27 schnetter Exp $

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
      
      subroutine carpettest_check_arguments (CCTK_ARGUMENTS)
      implicit none
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS
      integer i,j,k
      print *, "gfg ", gfg0, gfg1, gfg2
      print *, "arrg3 ", arrg30, arrg31, arrg32
      print *, "arrg2 ", arrg20, arrg21
      print *, "arrg1 ", arrg10
      print *, "scg"
      print *
      do k=1,gfg2
         do j=1,gfg1
            do i=1,gfg0
               gf(i,j,k) = i*10000 + j*100 + k
            end do
         end do
      end do
      do k=1,arrg32
         do j=1,arrg31
            do i=1,arrg30
               arr3(i,j,k) = i*10000 + j*100 + k
            end do
         end do
      end do
      do j=1,arrg21
         do i=1,arrg20
            arr2(i,j) = i*100 + j
         end do
      end do
      do i=1,arrg10
         arr1(i) = i
      end do
      sc = 42
      end
