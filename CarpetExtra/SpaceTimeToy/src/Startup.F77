c     -*-Fortran-*-
c     $Header: /home/eschnett/C/carpet/Carpet/CarpetExtra/SpaceTimeToy/src/Startup.F77,v 1.1 2001/03/17 16:05:56 eschnett Exp $
      
#include "cctk.h"
      
      subroutine SpaceTimeToy_Startup
      
      implicit none
      
      integer ierr
      
      call CCTK_RegisterBanner
     $     (ierr, "SpaceTimeToy: Evolutions of two Scalar Fields")
      
      end 
