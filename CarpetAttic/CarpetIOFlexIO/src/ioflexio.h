/* $Header:$ */

#ifndef CARPETIOFLEXIO_H
#define CARPETIOFLEXIO_H

#include "cctk_Arguments.h"

#ifdef __cplusplus
namespace CarpetIOFlexIO {
  extern "C" {
#endif
    
    /* Scheduled functions */
    int CarpetIOFlexIOStartup (void);
    int CarpetIOFlexIOReadData (CCTK_ARGUMENTS);
    
#ifdef __cplusplus
  } /* extern "C" */
} /* namespace CarpetIOFlexIO */
#endif

#endif /* !defined(CARPETIOFLEXIO_H) */
