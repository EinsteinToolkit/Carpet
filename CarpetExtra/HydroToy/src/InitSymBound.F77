c     -*-Fortran-*-
c     $Header: /home/eschnett/C/carpet/Carpet/CarpetExtra/HydroToy/src/InitSymBound.F77,v 1.1 2001/03/18 22:37:09 eschnett Exp $

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
      
      subroutine HydroToy_InitSymBound (CCTK_ARGUMENTS)
      
      implicit none
      
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_FUNCTIONS
      DECLARE_CCTK_PARAMETERS
      
      integer ierr
      
      call Cart3dSetTensorTypeVN (ierr, cctkGH, "HydroToy::u", "scalar")
      call Cart3dSetTensorTypeVN (ierr, cctkGH, "HydroToy::vx HydroToy::vy HydroToy::vz", "vector")
      
      end
