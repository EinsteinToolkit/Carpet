c	-*-Fortran-*-

 /*@@
   @file      Startup.F77
   @date      
   @author    Gabrielle Allen
   @desc 
              Register banner 
   @enddesc 
   @version $Header: /home/eschnett/C/carpet/Carpet/CarpetExtra/FOWaveToyF77/src/Startup.F77,v 1.2 2003/06/27 15:54:58 schnetter Exp $
 @@*/

#include "cctk.h"

      subroutine FOWaveToyF77_Startup

      implicit none
      integer ierr
   
      call CCTK_RegisterBanner(ierr, "FOWaveToyF77: Evolutions of a Scalar Field")

      end 
