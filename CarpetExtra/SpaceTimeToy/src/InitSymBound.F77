c     -*-Fortran-*-
c     $Header:$

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
      
      subroutine SpaceTimeToy_InitSymBound (CCTK_ARGUMENTS)
      
      implicit none
      
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_FUNCTIONS
      DECLARE_CCTK_PARAMETERS
      
      integer ierr
      
      call Cart3dSetTensorTypeVN (ierr, cctkGH, "SpaceTimeToy::phi", "scalar")
      call Cart3dSetTensorTypeVN (ierr, cctkGH, "SpaceTimeToy::psi", "scalar")
      
      end
