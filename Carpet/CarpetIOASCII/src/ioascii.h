/* $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetIOASCII/src/ioascii.h,v 1.6 2004/01/25 14:57:29 schnetter Exp $ */

#ifndef CARPETIOASCII_H
#define CARPETIOASCII_H

#include "cctk_Arguments.h"

#ifdef __cplusplus
namespace CarpetIOASCII {
  extern "C" {
#endif
    
    /* Scheduled functions */
    void CarpetIOASCIIStartup (void);
    void CarpetIOASCIIInit (CCTK_ARGUMENTS);
    
#ifdef __cplusplus
  } /* extern "C" */
} /* namespace CarpetIOASCII */
#endif

#endif /* !defined(CARPETIOASCII_H) */
