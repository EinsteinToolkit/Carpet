c	-*-Fortran-*-

 /*@@
   @file    InitSymBound.F77
   @date    
   @author  Gabrielle Allen, Erik Schnetter
   @desc
            Sets the symmetries for Wave Toy
   @enddesc
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

 /*@@
   @routine    FOWaveToyF77_InitSymBound
   @date      
   @author     Gabrielle Allen, Erik Schnetter
   @desc 
               Sets the symmetries for Wave Toy
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

      subroutine FOWaveToyF77_InitSymBound (CCTK_ARGUMENTS)
      
      implicit none
      
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_FUNCTIONS
      DECLARE_CCTK_PARAMETERS

      INTEGER one
      PARAMETER (one = 1)
      INTEGER ierr
      INTEGER sym(3)

      sym(1) = one
      sym(2) = one
      sym(3) = one

      call SetCartSymVN(ierr, cctkGH, sym,'wavetoy::pi')
      call SetCartSymVN(ierr, cctkGH, sym,'wavetoy::phi')
      call SetCartSymVN(ierr, cctkGH, sym,'wavetoy::phix')
      call SetCartSymVN(ierr, cctkGH, sym,'wavetoy::phiy')
      call SetCartSymVN(ierr, cctkGH, sym,'wavetoy::phiz')
      
      end
