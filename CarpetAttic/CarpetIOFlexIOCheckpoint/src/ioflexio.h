/* $Header: /home/eschnett/C/carpet/Carpet/CarpetAttic/CarpetIOFlexIOCheckpoint/src/ioflexio.h,v 1.1 2003/05/16 14:02:18 hawke Exp $ */

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

#ifdef __cplusplus
namespace CarpetCheckpointRestart {
  extern "C" {
#endif
    
    /* Scheduled functions */
    void CarpetChReEvolutionCheckpoint (const cGH*);
#ifdef __cplusplus
  } /* extern "C" */
} /* namespace CarpetIOFlexIO */
#endif


#endif /* !defined(CARPETIOFLEXIO_H) */
