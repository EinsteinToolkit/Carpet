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
   @author     Erik Schnetter
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

      INTEGER sym(3), symx(3), symy(3), symz(3)
      DATA sym  /+1, +1, +1/
      DATA symx /-1, +1, +1/
      DATA symy /+1, -1, +1/
      DATA symz /+1, +1, -1/
      INTEGER ierr

      call SetCartSymVN(ierr, cctkGH, sym ,'fowavetoy::pi')
      call SetCartSymVN(ierr, cctkGH, sym ,'fowavetoy::phi')
      call SetCartSymVN(ierr, cctkGH, symx,'fowavetoy::phix')
      call SetCartSymVN(ierr, cctkGH, symy,'fowavetoy::phiy')
      call SetCartSymVN(ierr, cctkGH, symz,'fowavetoy::phiz')
      
      end
