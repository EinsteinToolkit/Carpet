c     -*-Fortran-*-
c     $Header: /home/eschnett/C/carpet/Carpet/CarpetExtra/IDHydroToy/src/Startup.F77,v 1.1 2001/03/18 22:37:10 eschnett Exp $
      
#include "cctk.h"
      
      subroutine IDHydroToy_Startup
      
      implicit none
      
      integer ierr
      
      call CCTK_RegisterBanner
     $     (ierr, "IDHydroToy: Initial data for a Scalar and a Vector Field")
      
      end 
