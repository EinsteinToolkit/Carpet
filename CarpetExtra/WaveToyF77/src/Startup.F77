c	-*-Fortran-*-

 /*@@
   @file      Startup.F77
   @date      
   @author    Gabrielle Allen
   @desc 
              Register banner 
   @enddesc 
   @version $Header:$
 @@*/

#include "cctk.h"

      integer function WaveToyF77_Startup()

      implicit none
   
      integer ierr
      call CCTK_RegisterBanner(ierr, "WaveToyF77: Evolutions of a Scalar Field")

      WaveToyF77_Startup = 0

      end 
