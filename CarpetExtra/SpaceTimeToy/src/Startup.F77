c     -*-Fortran-*-
c     $Header:$
      
#include "cctk.h"
      
      subroutine SpaceTimeToy_Startup
      
      implicit none
      
      integer ierr
      
      call CCTK_RegisterBanner
     $     (ierr, "SpaceTimeToy: Evolutions of two Scalar Fields")
      
      end 
