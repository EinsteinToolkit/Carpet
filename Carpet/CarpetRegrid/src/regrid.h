/* $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetRegrid/src/regrid.h,v 1.1 2002/09/01 14:52:29 schnetter Exp $ */

#ifndef CARPETREGRID_H
#define CARPETREGRID_H

#include "cctk_Arguments.h"

#ifdef __cplusplus
namespace CarpetRegrid {
  extern "C" {
#endif
    
    /* Scheduled functions */
    int CarpetRegridStartup ();
    int CarpetRegridParamcheck (CCTK_ARGUMENTS);

#ifdef __cplusplus
  } /* extern "C" */
} /* namespace CarpetRegrid */
#endif

#endif /* !defined(CARPETREGRID_H) */
