c     -*-Fortran-*-
c     $Header:$
      
#include "cctk.h"
      
      subroutine IDHydroToy_Startup
      
      implicit none
      
      integer ierr
      
      call CCTK_RegisterBanner
     $     (ierr, "IDHydroToy: Initial data for a Scalar and a Vector Field")
      
      end 
