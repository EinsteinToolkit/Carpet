c	-*-Fortran-*-

 /*@@
   @file      Startup.F77
   @date      
   @author    Gabrielle Allen
   @desc 
              Register banner 
   @enddesc 
   @version $Header: /home/eschnett/C/carpet/Carpet/CarpetExtra/FOWaveToyF77/src/Startup.F77,v 1.1 2002/02/18 11:25:34 shawley Exp $
 @@*/

#include "cctk.h"

      integer function FOWaveToyF77_Startup()

      implicit none
   
      integer ierr
      call CCTK_RegisterBanner(ierr, "FOWaveToyF77: Evolutions of a Scalar Field")

      FOWaveToyF77_Startup = 0

      end 
