c	-*-Fortran-*-

 /*@@
   @file      Startup.F77
   @date      
   @author    Gabrielle Allen
   @desc 
              Register banner 
   @enddesc 
   @version $Header: /home/eschnett/C/carpet/Carpet/CarpetExtra/WaveToyF77/src/Startup.F77,v 1.1 2001/03/01 13:40:11 eschnett Exp $
 @@*/

#include "cctk.h"

      integer function WaveToyF77_Startup()

      implicit none
   
      integer ierr
      call CCTK_RegisterBanner(ierr, "WaveToyF77: Evolutions of a Scalar Field")

      WaveToyF77_Startup = 0

      end 
