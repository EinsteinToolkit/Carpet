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
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

 /*@@
   @routine    WaveToyF77_InitSymBound
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

      subroutine WaveToyF77_InitSymBound (CCTK_ARGUMENTS)
      
      implicit none
      
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_FUNCTIONS
      DECLARE_CCTK_PARAMETERS
      
      integer ierr
      
      call Cart3dSetTensorTypeVN (ierr, cctkGH, 'wavetoy::phi', 'SCALAR')
      
      end
