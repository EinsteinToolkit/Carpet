c     $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetTest/src/carpettest_check_arguments.F77,v 1.1 2001/07/04 12:29:54 schnetter Exp $

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
      
      subroutine carpettest_check_arguments (CCTK_ARGUMENTS)
      implicit none
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS
      print *, "gfg ", gfg0, gfg1, gfg2
      print *, "arrg3 ", arrg30, arrg31, arrg32
      print *, "arrg2 ", arrg20, arrg21
      print *, "arrg1 ", arrg10
      print *
      end
