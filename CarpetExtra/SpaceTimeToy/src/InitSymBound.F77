c     -*-Fortran-*-
c     $Header: /home/eschnett/C/carpet/Carpet/CarpetExtra/SpaceTimeToy/src/InitSymBound.F77,v 1.1 2001/03/17 16:05:56 eschnett Exp $

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
      
      subroutine SpaceTimeToy_InitSymBound (CCTK_ARGUMENTS)
      
      implicit none
      
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_FUNCTIONS
      DECLARE_CCTK_PARAMETERS
      
      integer ierr
      integer sym(3)
      
      sym(1) = 1
      sym(2) = 1
      sym(3) = 1
      
      call Cart3dSetTensorTypeVN (ierr, cctkGH, "SpaceTimeToy::phi", "scalar")
      call Cart3dSetTensorTypeVN (ierr, cctkGH, "SpaceTimeToy::psi", "scalar")
      
      end
