#ifndef CARPETIOASCII_H
#define CARPETIOASCII_H

#include "cctk_Arguments.h"

#ifdef __cplusplus
namespace CarpetIOASCII {
  extern "C" {
#endif
    
    /* Scheduled functions */
    int CarpetIOASCIIStartup (void);
    void CarpetIOASCIIInit (CCTK_ARGUMENTS);
    
#ifdef __cplusplus
  } /* extern "C" */
} /* namespace CarpetIOASCII */
#endif

#endif /* !defined(CARPETIOASCII_H) */
