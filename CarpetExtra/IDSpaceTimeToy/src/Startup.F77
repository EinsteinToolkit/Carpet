c     -*-Fortran-*-
c     $Header: /home/eschnett/C/carpet/Carpet/CarpetExtra/IDSpaceTimeToy/src/Startup.F77,v 1.1 2001/03/17 16:05:56 eschnett Exp $
      
#include "cctk.h"
      
      subroutine IDSpaceTimeToy_Startup
      
      implicit none
      
      integer ierr
      
      call CCTK_RegisterBanner
     $     (ierr, "IDSpaceTimeToy: Initial data for two Scalar Fields")
      
      end 
