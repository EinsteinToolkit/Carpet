c     -*-Fortran-*-
c     $Header: /home/eschnett/C/carpet/Carpet/CarpetExtra/HydroToy/src/Startup.F77,v 1.1 2001/03/18 22:37:09 eschnett Exp $
      
#include "cctk.h"
      
      subroutine HydroToy_Startup
      
      implicit none
      
      integer ierr
      
      call CCTK_RegisterBanner
     $     (ierr, "HydroToy: Evolutions of a Scalar and a Vector Field")
      
      end 
