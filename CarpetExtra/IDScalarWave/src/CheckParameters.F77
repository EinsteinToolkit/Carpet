c	-*-Fortran-*-

 /*@@
   @file      CheckParameters.F77
   @date      
   @author    Gabrielle Allen
   @desc 
              Check parameters for the wave equation initial data
   @enddesc 
 @@*/

#include "cctk.h" 
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"


 /*@@
   @routine    IDScalarWave_CheckParameters
   @date       
   @author     Gabrielle Allen
   @desc 
               Check parameters for the wave equation initial data
   @enddesc 
   @calls      
   @calledby   
   @history 
 
   @endhistory 

@@*/

      subroutine IDScalarWave_CheckParameters(CCTK_ARGUMENTS)

      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS

      INTEGER CCTK_Equals

      if (CCTK_Equals(initial_data,"box").eq.1) then

         if (CCTK_Equals(type, "box").eq.0) then
            call CCTK_PARAMWARN("Must have a box grid with box initial data")
         end if

         if (kx.eq.0 .or. ky.eq.0 .or. kz.eq.0) then
            call CCTK_PARAMWARN("Cannot have zero kx,ky,kz for box initial data")
         end if

      end if

      return
      end
