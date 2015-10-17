c     -*-Fortran-*-
c     $Header:$
      
#include "cctk.h"
      
      subroutine IDSpaceTimeToy_Startup
      
      implicit none
      
      integer ierr
      
      call CCTK_RegisterBanner
     $     (ierr, "IDSpaceTimeToy: Initial data for two Scalar Fields")
      
      end 
