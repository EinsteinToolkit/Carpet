c     -*-Fortran-*-
c     $Header:$

 /*@@
   @file      CheckParameters.F77
   @date      
   @author    Gabrielle Allen
   @desc 
              Check parameters for the wave equation initial data
   @enddesc 
 @@*/

#include "cctk.h" 
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"


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
      DECLARE_CCTK_FUNCTIONS
      DECLARE_CCTK_PARAMETERS

      if (CCTK_EQUALS(initial_data,"box")) then

         if (kx.eq.0 .or. ky.eq.0 .or. kz.eq.0) then
            call CCTK_PARAMWARN("Cannot have zero kx,ky,kz for box initial data")
         end if

      end if

      return
      end
