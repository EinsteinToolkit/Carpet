#ifndef CARPETREGRID_H
#define CARPETREGRID_H

#include "cctk_Arguments.h"

#ifdef __cplusplus
namespace CarpetRegrid {
  extern "C" {
#endif
    
    /* Scheduled functions */
    int CarpetRegridStartup ();
    void CarpetRegridParamcheck (CCTK_ARGUMENTS);

#ifdef __cplusplus
  } /* extern "C" */
} /* namespace CarpetRegrid */
#endif

#endif /* !defined(CARPETREGRID_H) */
