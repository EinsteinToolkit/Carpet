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

      subroutine FOWaveToyF77_Startup

      implicit none
      integer ierr
   
      call CCTK_RegisterBanner(ierr, "FOWaveToyF77: Evolutions of a Scalar Field")

      end 
