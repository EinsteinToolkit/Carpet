c     -*-Fortran-*-
c     $Header:$
      
#include "cctk.h"
      
      subroutine HydroToy_Startup
      
      implicit none
      
      integer ierr
      
      call CCTK_RegisterBanner
     $     (ierr, "HydroToy: Evolutions of a Scalar and a Vector Field")
      
      end 
